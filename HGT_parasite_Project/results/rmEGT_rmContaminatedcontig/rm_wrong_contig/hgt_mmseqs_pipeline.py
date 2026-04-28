import pandas as pd
import argparse
import sys
import os
import subprocess
import shutil
from Bio import SeqIO

def run_mmseqs(input_fasta, output_prefix, min_seq_id, cov, threads=8):
    """
    调用 mmseqs easy-cluster
    生成文件: {output_prefix}_cluster.tsv, {output_prefix}_rep_seq.fasta, etc.
    """
    tmp_dir = f"{output_prefix}_tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    
    cmd = [
        "mmseqs", "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov),
        "--cov-mode", "1", # 0: 双向覆盖 (更严格), 1: 目标覆盖, 2: 查询覆盖
        "--threads", str(threads),
        "-v", "1" # 减少日志输出
    ]
    
    try:
        subprocess.run(cmd, check=True)
        # 清理临时文件夹
        shutil.rmtree(tmp_dir)
        return f"{output_prefix}_cluster.tsv", f"{output_prefix}_rep_seq.fasta"
    except subprocess.CalledProcessError:
        print(f"Error: MMseqs2 failed for {input_fasta}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: mmseqs command not found. Please install via: conda install -c bioconda mmseqs2")
        sys.exit(1)

def parse_mmseqs_cluster(cluster_tsv):
    """
    解析 mmseqs 的 _cluster.tsv 文件 (格式: Rep \t Member)
    返回: {Rep_ID: [List of Members]}
    """
    cluster_map = {}
    with open(cluster_tsv, 'r') as f:
        for line in f:
            rep, member = line.strip().split('\t')
            if rep not in cluster_map:
                cluster_map[rep] = []
            cluster_map[rep].append(member)
    return cluster_map

def generate_output_table(original_df, rep_list, cluster_map, output_file, id_col='Accession'):
    """
    根据聚类结果生成新的表格
    只保留代表序列的行，并增加 Members 列
    """
    # 筛选出代表序列的行
    rep_df = original_df[original_df[id_col].isin(rep_list)].copy()
    
    # 建立 Rep -> Members String 的映射
    members_str_map = {}
    size_map = {}
    
    for rep, members in cluster_map.items():
        # 将 list 转为逗号分隔字符串
        members_str_map[rep] = ",".join(members)
        size_map[rep] = len(members)
        
    rep_df['Cluster_Size'] = rep_df[id_col].map(size_map)
    rep_df['Cluster_Members'] = rep_df[id_col].map(members_str_map)
    
    # 排序：按 Cluster Size 降序
    rep_df = rep_df.sort_values('Cluster_Size', ascending=False)
    
    rep_df.to_csv(output_file, sep='\t', index=False)
    print(f"表格生成完毕: {output_file} (Rows: {len(rep_df)})")
    return rep_df

def main():
    parser = argparse.ArgumentParser(description="HGT MMseqs2 两步聚类流程")
    parser.add_argument('--fasta', required=True, help="总蛋白 Fasta 文件")
    parser.add_argument('--table', required=True, help="总信息表格 (.txt)")
    parser.add_argument('-o', '--outdir', required=True, help="输出目录")
    
    # 可调参数
    parser.add_argument('--id1', type=float, default=0.8, help="Step 1 物种内聚类 Identity (默认 0.8)")
    parser.add_argument('--cov1', type=float, default=0.8, help="Step 1 物种内聚类 Coverage (默认 0.8)")
    parser.add_argument('--id2', type=float, default=0.5, help="Step 2 全局聚类 Identity (默认 0.5)")
    parser.add_argument('--cov2', type=float, default=0.8, help="Step 2 全局聚类 Coverage (默认 0.8)")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. 读取数据
    print(">>> [1/5] 读取输入文件...")
    df = pd.read_csv(args.table, sep='\t')
    # 确保Accession列名为字符串，防止纯数字ID报错
    df['Accession'] = df['Accession'].astype(str)
    
    # 读取序列到内存字典
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    
    # 2. Step 1: 物种内聚类
    print(f">>> [2/5] 开始 Step 1: 物种内聚类 (ID={args.id1}, Cov={args.cov1})...")
    step1_dir = os.path.join(args.outdir, "step1_intra_species")
    os.makedirs(step1_dir, exist_ok=True)
    
    # 获取物种列表
    species_col = 'Specie' # 根据你的表头调整
    if species_col not in df.columns:
        print(f"Error: 列名 '{species_col}' 不在表格中，请检查表头。")
        sys.exit(1)
        
    step1_reps_fasta_list = [] # 存储所有物种的代表序列文件路径
    
    for sp, group in df.groupby(species_col):
        safe_sp = sp.replace(" ", "_").replace("/", "_")
        sp_fasta_path = os.path.join(step1_dir, f"{safe_sp}.in.fasta")
        sp_out_prefix = os.path.join(step1_dir, f"{safe_sp}")
        
        # 写入该物种的临时 fasta
        count = 0
        with open(sp_fasta_path, 'w') as f:
            for acc in group['Accession']:
                if acc in seq_dict:
                    SeqIO.write(seq_dict[acc], f, "fasta")
                    count += 1
        
        if count == 0: continue
        
        # 运行 mmseqs
        _, rep_fasta = run_mmseqs(sp_fasta_path, sp_out_prefix, args.id1, args.cov1)
        step1_reps_fasta_list.append(rep_fasta)

    # 3. 合并 Step 1 的结果并生成表格
    print(">>> [3/5] 合并 Step 1 结果...")
    # 我们需要读取所有的 cluster.tsv 来生成 Step 1 的总表
    step1_cluster_map = {} # 全局映射
    step1_rep_ids = []
    
    # 合并所有 Step 1 的代表序列用于 Step 2
    step1_combined_fasta = os.path.join(args.outdir, "step1_combined_reps.fasta")
    
    with open(step1_combined_fasta, 'w') as outfile:
        for fname in os.listdir(step1_dir):
            if fname.endswith("_cluster.tsv"):
                # 解析聚类关系
                partial_map = parse_mmseqs_cluster(os.path.join(step1_dir, fname))
                step1_cluster_map.update(partial_map)
                step1_rep_ids.extend(partial_map.keys())
            
            if fname.endswith("_rep_seq.fasta"):
                # 合并序列
                with open(os.path.join(step1_dir, fname), 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
    
    # 生成 Step 1 表格
    step1_table_path = os.path.join(args.outdir, "Step1_IntraSpecies_Clustered.txt")
    generate_output_table(df, step1_rep_ids, step1_cluster_map, step1_table_path)
    
    # 4. Step 2: 全局聚类
    print(f">>> [4/5] 开始 Step 2: 全局跨物种聚类 (ID={args.id2}, Cov={args.cov2})...")
    step2_dir = os.path.join(args.outdir, "step2_global")
    os.makedirs(step2_dir, exist_ok=True)
    
    step2_prefix = os.path.join(step2_dir, "global")
    
    # 对 Step 1 的代表序列进行聚类
    cluster_tsv_2, _ = run_mmseqs(step1_combined_fasta, step2_prefix, args.id2, args.cov2)
    
    # 5. 生成最终表格
    print(">>> [5/5] 解析全局聚类并生成最终表格...")
    
    # 解析 Step 2 关系: Step2_Rep -> [Step1_Reps]
    step2_map_temp = parse_mmseqs_cluster(cluster_tsv_2)
    
    # 关键步骤：我们需要“展开”这个关系，变成 Step2_Rep -> [Original_Accessions]
    # 因为 Step 1 的 Rep 背后其实包含了好几个原始 ID
    step2_final_map = {}
    
    for global_rep, member_reps in step2_map_temp.items():
        all_original_accs = []
        for rep1 in member_reps:
            # 找到 rep1 在第一步中代表了哪些原始序列
            if rep1 in step1_cluster_map:
                all_original_accs.extend(step1_cluster_map[rep1])
            else:
                # 理论上不应该发生，除非序列ID匹配有问题
                all_original_accs.append(rep1)
        
        step2_final_map[global_rep] = all_original_accs
        
    step2_final_rep_ids = list(step2_final_map.keys())
    
    # 生成 Step 2 表格
    step2_table_path = os.path.join(args.outdir, "Step2_Global_Clustered.txt")
    generate_output_table(df, step2_final_rep_ids, step2_final_map, step2_table_path)
    
    print("\n==========================================")
    print("流程完成！")
    print(f"1. Step 1 (物种内去冗余) 表格: {step1_table_path}")
    print(f"2. Step 2 (HGT家族/跨物种) 表格: {step2_table_path}")
    print("==========================================")

if __name__ == "__main__":
    main()
