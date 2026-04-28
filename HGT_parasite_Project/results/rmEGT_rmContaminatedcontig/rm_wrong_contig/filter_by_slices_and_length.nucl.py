import pandas as pd
import argparse
import sys
import numpy as np
from tqdm import tqdm

# 注册 tqdm 到 pandas
tqdm.pandas(desc="Processing IDs")

def parse_slice_accession_vectorized(series):
    """
    针对 Series 的向量化 ID 解析
    """
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
    parser = argparse.ArgumentParser(description="[极速版V2] 结合长度、锚定原则与硬阈值剔除污染")
    parser.add_argument('--slice_info', required=True, help="基因组切片筛选文件")
    parser.add_argument('--hgt_list', required=True, help="HGT 候选列表")
    parser.add_argument('--lengths', required=True, help="Step1 生成的 all_contig_lengths.txt")
    parser.add_argument('-o', '--output', required=True, help="输出文件前缀")
    
    # 原有参数
    parser.add_argument('--slice_eff_len', type=int, default=1300, help="切片有效长度 (步长)")
    parser.add_argument('--min_euk_anchor', type=int, default=1, help="最少需要几个真核切片作为锚点")
    
    # === 新增参数 ===
    parser.add_argument('--min_len', type=int, default=0, help="Contig 最小长度阈值 (例如 10000)，低于此长度直接剔除")
    
    args = parser.parse_args()

    # ==========================================
    # 1. 加载长度信息
    # ==========================================
    print(f"[1/5] 加载 Contig 长度: {args.lengths} ...")
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
    print(f"      已加载 {len(len_map)} 条 Contig 长度。")

    # ==========================================
    # 2. 极速加载切片文件
    # ==========================================
    print(f"[2/5] 加载并分析切片文件...")
    try:
        df_slice = pd.read_csv(
            args.slice_info, 
            sep='\t', 
            engine='python', 
            on_bad_lines='skip'
        )
        if 'Accession' not in df_slice.columns and df_slice.index.name == 'Accession':
            df_slice = df_slice.reset_index()
        else:
            df_slice.columns = [c.strip() for c in df_slice.columns]
    except Exception as e:
        print(f"[ERROR] 读取切片文件失败: {e}")
        sys.exit(1)

    # ==========================================
    # 3. 向量化计算
    # ==========================================
    print("      -> 正在进行向量化计算...")
    df_slice['Contig_ID'] = parse_slice_accession_vectorized(df_slice['Accession'])

    cols_to_numeric = ['alien_index', 'Out_pct', 'Close_num', 'Self_num']
    for col in cols_to_numeric:
        df_slice[col] = pd.to_numeric(df_slice[col], errors='coerce').fillna(0)
    
    df_slice['First_distant_taxinfo'] = df_slice['First_distant_taxinfo'].astype(str).fillna('')

    # 判定逻辑
    cond_bac_1 = (df_slice['alien_index'] > 0) & (df_slice['Out_pct'] > 0.8)
    cond_bac_2 = (df_slice['First_distant_taxinfo'].str.contains('Bacteria', case=False)) & (df_slice['Out_pct'] > 0.8)
    mask_bac = cond_bac_1 | cond_bac_2
    mask_euk = (df_slice['alien_index'] < -5) | (df_slice['Close_num'] > 0) | (df_slice['Self_num'] > 0)

    # 聚合
    df_slice['is_bac'] = mask_bac.astype(int)
    df_slice['is_euk'] = mask_euk.astype(int)
    
    stats = df_slice.groupby('Contig_ID').agg({
        'is_bac': 'sum',
        'is_euk': 'sum',
        'Accession': 'size' 
    }).rename(columns={'Accession': 'total_hits'})

    # ==========================================
    # 4. 长度映射与最终判定 (核心修改区)
    # ==========================================
    print(f"[3/5] 执行多重规则判定...")
    
    stats['real_len'] = stats.index.map(len_map).fillna(10000) # 找不到长度的默认给个安全值
    stats['mapped_len'] = stats['total_hits'] * args.slice_eff_len
    stats['coverage'] = (stats['mapped_len'] / stats['real_len']).clip(upper=1.0)
    stats['bac_ratio'] = stats['is_bac'] / stats['total_hits']

    stats['is_dirty'] = False
    stats['reason'] = 'Pass'

    # --- 规则 1: 绝对长度硬阈值 (用户新增需求) ---
    if args.min_len > 0:
        mask_too_short = stats['real_len'] < args.min_len
        stats.loc[mask_too_short, 'is_dirty'] = True
        stats.loc[mask_too_short, 'reason'] = 'Too_Short'

    # --- 规则 2: 高覆盖度下的污染 ---
    # 只有当还没有被判定为 Too_Short 时才应用
    mask_check = ~stats['is_dirty']
    
    mask_high_cov_contam = mask_check & (stats['coverage'] > 0.3) & (stats['bac_ratio'] > 0.5)
    stats.loc[mask_high_cov_contam, 'is_dirty'] = True
    stats.loc[mask_high_cov_contam, 'reason'] = 'Contam_HighCov'

    # --- 规则 3: 低覆盖度下的锚定 ---
    mask_check = ~stats['is_dirty'] # 更新 mask
    mask_no_anchor = mask_check & (stats['coverage'] <= 0.3) & (stats['is_bac'] > 0) & (stats['is_euk'] < args.min_euk_anchor)
    stats.loc[mask_no_anchor, 'is_dirty'] = True
    stats.loc[mask_no_anchor, 'reason'] = 'No_Anchor'

    dirty_contigs = set(stats[stats['is_dirty']].index)
    
    print("      判定详情 (被剔除的原因分布):")
    print(stats['reason'].value_counts())

    # ==========================================
    # 5. 清洗 HGT 表
    # ==========================================
    print(f"[4/5] 清洗 HGT 列表...")
    try:
        df_hgt = pd.read_csv(args.hgt_list, sep='\t', on_bad_lines='skip')
    except:
        sys.exit(1)

    acc_col = 'Accession' if 'Accession' in df_hgt.columns else df_hgt.columns[0]
    initial_count = len(df_hgt)
    hgt_contig_ids = df_hgt[acc_col].astype(str).str.split(':').str[0]
    
    # 过滤
    mask_keep = ~hgt_contig_ids.isin(dirty_contigs)
    
    # 如果开启了 min_len，我们还需要检查那些没有出现在切片文件里的 HGT 所在的 Contig 长度
    if args.min_len > 0:
        # 获取每一行 HGT 对应 Contig 的长度
        hgt_lens = hgt_contig_ids.map(len_map).fillna(0) # 找不到长度视为0
        mask_len_ok = hgt_lens >= args.min_len
        # 合并条件: 不脏 AND 长度够
        mask_keep = mask_keep & mask_len_ok
    
    df_clean = df_hgt[mask_keep]
    
    print(f"[5/5] 输出结果: 原始 {initial_count} -> 清洗后 {len(df_clean)}")
    df_clean.to_csv(f"{args.output}.len_validated.txt", sep='\t', index=False)

if __name__ == "__main__":
    main()
