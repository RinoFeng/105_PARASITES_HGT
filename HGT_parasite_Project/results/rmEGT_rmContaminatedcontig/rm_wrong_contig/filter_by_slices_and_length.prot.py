import pandas as pd
import argparse
import sys
import numpy as np
from tqdm import tqdm

# 注册 tqdm 到 pandas
tqdm.pandas(desc="Processing")

def parse_slice_accession_vectorized(series):
    """ 向量化解析切片 Accession 获取 Contig ID """
    s = series.astype(str)
    res_dash = s.str.rsplit('-', n=2)
    mask_dash = res_dash.str.len() == 3
    res_under = s.str.rsplit('_', n=2)
    mask_under = res_under.str.len() == 3
    
    final_ids = s.copy()
    final_ids[mask_dash] = res_dash[mask_dash].str[0]
    mask_final_under = (~mask_dash) & (mask_under)
    final_ids[mask_final_under] = res_under[mask_final_under].str[0]
    return final_ids

def main():
    parser = argparse.ArgumentParser(description="[蛋白版V2] 基于Contig质量过滤HGT蛋白 (无映射即剔除)")
    
    # 核心输入
    parser.add_argument('--slice_info', required=True, help="基因组切片筛选文件")
    parser.add_argument('--hgt_list', required=True, help="HGT 蛋白候选列表")
    parser.add_argument('--blast_map', required=True, help="蛋白比对到基因组的 BLAST 结果")
    parser.add_argument('--lengths', required=True, help="基因组 Contig 长度文件")
    parser.add_argument('-o', '--output', required=True, help="输出文件前缀")
    
    # 过滤参数
    parser.add_argument('--slice_eff_len', type=int, default=1300, help="切片有效步长")
    parser.add_argument('--min_euk_anchor', type=int, default=1, help="最少真核锚点数")
    parser.add_argument('--min_len', type=int, default=0, help="Contig 最小长度阈值")
    
    args = parser.parse_args()

    # ==========================================
    # 1. 加载 Contig 长度
    # ==========================================
    print(f"[1/6] 加载 Contig 长度: {args.lengths} ...")
    len_map = {}
    try:
        with open(args.lengths, 'r') as f:
            header = next(f) 
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    len_map[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"[ERROR] 读取长度文件失败: {e}")
        sys.exit(1)

    # ==========================================
    # 2. 分析切片文件
    # ==========================================
    print(f"[2/6] 分析切片数据，鉴定脏 Contig...")
    try:
        df_slice = pd.read_csv(args.slice_info, sep='\t', engine='python', on_bad_lines='skip')
        if 'Accession' not in df_slice.columns and df_slice.index.name == 'Accession':
            df_slice = df_slice.reset_index()
        else:
            df_slice.columns = [c.strip() for c in df_slice.columns]
            
        df_slice['Contig_ID'] = parse_slice_accession_vectorized(df_slice['Accession'])
        
        for col in ['alien_index', 'Out_pct', 'Close_num', 'Self_num']:
            df_slice[col] = pd.to_numeric(df_slice[col], errors='coerce').fillna(0)
        df_slice['First_distant_taxinfo'] = df_slice['First_distant_taxinfo'].astype(str).fillna('')
        
        mask_bac = ((df_slice['alien_index'] > 0) & (df_slice['Out_pct'] > 0.8)) | \
                   ((df_slice['First_distant_taxinfo'].str.contains('Bacteria', case=False)) & (df_slice['Out_pct'] > 0.8))
        mask_euk = (df_slice['alien_index'] < -5) | (df_slice['Close_num'] > 0) | (df_slice['Self_num'] > 0)
        
        df_slice['is_bac'] = mask_bac.astype(int)
        df_slice['is_euk'] = mask_euk.astype(int)
        
        stats = df_slice.groupby('Contig_ID').agg({
            'is_bac': 'sum', 'is_euk': 'sum', 'Accession': 'size'
        }).rename(columns={'Accession': 'total_hits'})
        
        stats['real_len'] = stats.index.map(len_map).fillna(10000)
        stats['mapped_len'] = stats['total_hits'] * args.slice_eff_len
        stats['coverage'] = (stats['mapped_len'] / stats['real_len']).clip(upper=1.0)
        stats['bac_ratio'] = stats['is_bac'] / stats['total_hits']
        
        stats['is_dirty'] = False
        stats['reason'] = 'Pass'
        
        # 判定逻辑
        if args.min_len > 0:
            stats.loc[stats['real_len'] < args.min_len, ['is_dirty', 'reason']] = [True, 'Too_Short']
            
        mask_check = ~stats['is_dirty']
        mask_high_cov = mask_check & (stats['coverage'] > 0.3) & (stats['bac_ratio'] > 0.5)
        stats.loc[mask_high_cov, ['is_dirty', 'reason']] = [True, 'Contam_HighCov']
        
        mask_check = ~stats['is_dirty']
        mask_no_anchor = mask_check & (stats['coverage'] <= 0.3) & (stats['is_bac'] > 0) & (stats['is_euk'] < args.min_euk_anchor)
        stats.loc[mask_no_anchor, ['is_dirty', 'reason']] = [True, 'No_Anchor']
        
        dirty_contigs_set = set(stats[stats['is_dirty']].index)
        print(f"      Contig 判定完成。脏 Contig: {len(dirty_contigs_set)}")

    except Exception as e:
        print(f"[ERROR] 切片分析失败: {e}")
        sys.exit(1)

    # ==========================================
    # 3. 加载 BLAST 映射
    # ==========================================
    print(f"[3/6] 读取 BLAST 映射: {args.blast_map} ...")
    try:
        df_blast = pd.read_csv(args.blast_map, sep='\t', header=None, usecols=[0, 1], names=['Protein', 'Contig'])
        df_blast = df_blast.drop_duplicates(subset=['Protein'], keep='first')
        prot_to_contig = df_blast.set_index('Protein')['Contig'].to_dict()
        print(f"      已建立 {len(prot_to_contig)} 条映射。")
    except Exception as e:
        print(f"[ERROR] BLAST 读取失败: {e}")
        sys.exit(1)

    # ==========================================
    # 4. 过滤 HGT 蛋白列表 (核心修改处)
    # ==========================================
    print(f"[4/6] 清洗 HGT 蛋白列表...")
    try:
        df_hgt = pd.read_csv(args.hgt_list, sep='\t', on_bad_lines='skip')
    except:
        sys.exit(1)

    acc_col = 'Accession' if 'Accession' in df_hgt.columns else df_hgt.columns[0]
    initial_count = len(df_hgt)
    
    # 统计剔除原因
    drop_reasons = {'Dirty_Contig': 0, 'No_Mapping': 0, 'Keep': 0}
    clean_rows = []
    
    for _, row in tqdm(df_hgt.iterrows(), total=initial_count):
        prot_acc = str(row[acc_col])
        contig_id = prot_to_contig.get(prot_acc)
        
        # --- 修改开始: 无映射直接剔除 ---
        if not contig_id:
            drop_reasons['No_Mapping'] += 1
            continue  # 只有 continue，没有 append
        # --- 修改结束 ---
            
        # 检查是否脏
        if contig_id in dirty_contigs_set:
            drop_reasons['Dirty_Contig'] += 1
            continue
            
        # 检查长度 (即使不在dirty_set里)
        c_len = len_map.get(contig_id, 0)
        if args.min_len > 0 and c_len < args.min_len:
             drop_reasons['Dirty_Contig'] += 1
             continue
             
        # 保留
        drop_reasons['Keep'] += 1
        clean_rows.append(row)
        
    df_clean = pd.DataFrame(clean_rows)
    
    print(f"\n[5/6] 过滤统计:")
    print(f"      原始数量: {initial_count}")
    print(f"      剔除 (脏Contig/过短): {drop_reasons['Dirty_Contig']}")
    print(f"      剔除 (无基因组映射): {drop_reasons['No_Mapping']}  <-- 此处已剔除")
    print(f"      最终保留: {len(df_clean)}")

    # ==========================================
    # 5. 输出
    # ==========================================
    out_file = f"{args.output}.prot_validated.txt"
    df_clean.to_csv(out_file, sep='\t', index=False)
    print(f"[SUCCESS] 结果保存至: {out_file}")

if __name__ == "__main__":
    main()
