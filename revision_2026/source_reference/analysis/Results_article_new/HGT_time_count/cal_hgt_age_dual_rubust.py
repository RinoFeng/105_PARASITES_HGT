import os
import sys
import csv
from Bio import Phylo
from collections import defaultdict

# ================= 核心配置区域 =================

# 1. 树文件路径
TREE_FILE = "phylogenetic_tree/105_parasites_timetree.timetree.nex"

# 2. 物种列表文件
NAME_MAP_FILE = "105_name2accfile.txt"

# 3. 任务列表
TASKS = [
    {
        "name": "TBLASTN_Protein",
        "input_dir": "out_tblastn",
        "suffix": ".tblastn.txt",
        "output_file": "Final_HGT_Robust_TBLASTN.csv"
    },
    {
        "name": "TBLASTX_Nucleotide",
        "input_dir": "out_tblastx",
        "suffix": ".tblastx.txt",
        "output_file": "Final_HGT_Robust_TBLASTX.csv"
    }
]

# 4. 过滤阈值
MIN_IDENTITY = 50.0   # 最小一致性 (%)
MAX_EVALUE = 1e-5     # 最大 E-value
MIN_COVERAGE = 80.0    # [要求] 最小覆盖度 (qcovs), 若有需要可设为 30.0 等

# 5. 鲁棒性分析参数
#    如果基因在超过 85% 的物种中存在，视为垂直遗传，标记为 Excluded
VERTICAL_THRESHOLD_RATIO = 1.0 

#    [离群值检测]
#    如果去掉某个物种后，时间缩短了 50% 以上，则认为该物种是离群值
AGE_DROP_THRESHOLD = 0.5 
MIN_AGE_TO_CHECK = 200.0

# ===============================================

def load_species_mapping(map_file):
    """建立文件名 -> 树上物种名 的映射"""
    file_to_name = {}
    print(f"[Init] Loading species mapping from {map_file}...")
    try:
        with open(map_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split() 
                if len(parts) < 2: continue
                
                filename = parts[-1] 
                tree_name = "_".join(parts[:-1])
                file_to_name[filename] = tree_name
    except Exception as e:
        print(f"[Error] Failed to read mapping file: {e}")
        sys.exit(1)
    return file_to_name

def parse_blast_results(blast_dir, mapping, suffix):
    """解析 BLAST 结果，包含 qcovs 过滤"""
    hgt_distribution = defaultdict(set)
    processed_count = 0
    
    print(f"   -> Scanning directory: {blast_dir}")
    
    if not os.path.exists(blast_dir):
        print(f"   [Warning] Directory {blast_dir} does not exist! Skipping.")
        return {}

    for fname in os.listdir(blast_dir):
        if not fname.endswith(suffix): continue
            
        core_filename = fname.replace(suffix, "")
        if core_filename not in mapping: continue
            
        species_name = mapping[core_filename]
        file_path = os.path.join(blast_dir, fname)
        processed_count += 1
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    # 确保有足够的列 (qcovs 通常是第 13 列，索引 12)
                    if len(cols) < 13: continue 
                    
                    qseqid = cols[0]
                    pident = float(cols[2])
                    evalue = float(cols[10])
                    qcovs = float(cols[12]) # 获取覆盖度
                    
                    if evalue <= MAX_EVALUE and pident >= MIN_IDENTITY and qcovs >= MIN_COVERAGE:
                        hgt_distribution[qseqid].add(species_name)
        except Exception as e:
            continue

    print(f"   -> Parsed {processed_count} files. Found {len(hgt_distribution)} unique candidates.")
    return hgt_distribution

def get_node_age(tree, node):
    """计算节点到末端的距离 (时间)"""
    if node.is_terminal(): return 0.0
    terminals = node.get_terminals()
    if not terminals: return 0.0
    return tree.distance(node, terminals[0])

def calculate_robust_age(tree, species_list, gene_id):
    """
    计算时间并检测离群值 (Leave-One-Out)
    """
    count = len(species_list)
    
    # 1. 基础情况
    if count == 1:
        return 0.0, "Tip (Single Species)", []
        
    mrca = tree.common_ancestor(species_list)
    initial_age = get_node_age(tree, mrca)
    
    if initial_age < MIN_AGE_TO_CHECK or count <= 2:
        return initial_age, "Clade MRCA", []

    # 2. 离群值检测
    best_age = initial_age
    outlier_removed = None
    
    for sp_to_remove in species_list:
        sub_list = [s for s in species_list if s != sp_to_remove]
        sub_mrca = tree.common_ancestor(sub_list)
        sub_age = get_node_age(tree, sub_mrca)
        
        # 判定：时间是否大幅缩短
        if sub_age < initial_age * (1 - AGE_DROP_THRESHOLD):
            if sub_age < best_age:
                best_age = sub_age
                outlier_removed = sp_to_remove

    if outlier_removed:
        return best_age, "Robust (Outlier Removed)", [outlier_removed]
    else:
        return initial_age, "Clade MRCA", []

def process_task(task, mapping, tree, valid_species, total_species):
    """执行单个任务"""
    print(f"\n=== Starting Task: {task['name']} ===")
    
    hgt_data = parse_blast_results(task['input_dir'], mapping, task['suffix'])
    
    if not hgt_data:
        print("   -> No data found.")
        return

    print(f"   -> Analyzing ages and filtering outliers...")
    threshold_count = total_species * VERTICAL_THRESHOLD_RATIO
    
    with open(task['output_file'], 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        # [修改] 增加了最后一列 All_Species_List
        header = ["Gene_ID", "Age_MYA", "Sp_Count", "Ratio", "Status", "Method", "Removed_Outliers", "All_Species_List"]
        writer.writerow(header)
        
        count_kept = 0
        count_excluded = 0
        
        for gene, sp_set in hgt_data.items():
            # 过滤掉不在树上的无效物种
            present_sp = [s for s in sp_set if s in valid_species]
            if not present_sp: continue
            
            count = len(present_sp)
            ratio = count / total_species
            
            # 生成全物种字符串 (分号分隔)
            species_str = ";".join(present_sp)
            
            # === 过滤器 1: 过于普遍 (可能是垂直遗传) ===
            if count > threshold_count:
                writer.writerow([gene, "NA", count, f"{ratio:.2f}", "Excluded", "Too Frequent", "", species_str])
                count_excluded += 1
                continue
            
            # === 过滤器 2: 计算时间并剔除离群值 ===
            try:
                age, method, outliers = calculate_robust_age(tree, present_sp, gene)
                
                outlier_str = ";".join(outliers) if outliers else ""
                
                writer.writerow([gene, f"{age:.4f}", count, f"{ratio:.2f}", "Kept", method, outlier_str, species_str])
                count_kept += 1
                
            except Exception as e:
                pass

    print(f"   -> Task Finished.")
    print(f"      Kept (HGT): {count_kept}")
    print(f"      Excluded (Vertical): {count_excluded}")
    print(f"      Result saved to: {task['output_file']}")

def main():
    file2name = load_species_mapping(NAME_MAP_FILE)
    
    print(f"[Init] Loading phylogenetic tree from {TREE_FILE}...")
    try:
        trees = list(Phylo.parse(TREE_FILE, "nexus"))
        tree = trees[0]
        valid_tree_species = set(term.name for term in tree.get_terminals())
        total_species = len(valid_tree_species)
        print(f"       Tree contains {total_species} species.")
    except Exception as e:
        print(f"[Error] Could not load tree: {e}")
        sys.exit(1)

    for task in TASKS:
        process_task(task, file2name, tree, valid_tree_species, total_species)

    print("\n[Success] All tasks completed!")

if __name__ == "__main__":
    main()
