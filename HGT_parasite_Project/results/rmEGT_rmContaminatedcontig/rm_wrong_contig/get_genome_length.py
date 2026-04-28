import sys
import os
import glob
import gzip

def parse_fasta_length(file_path):
    """
    读取 fasta 文件并返回 {seq_id: length} 字典
    支持 .fa, .fna, .fasta 以及 .gz 压缩格式
    """
    lengths = {}
    current_id = None
    current_len = 0
    
    # 自动处理 gzip
    open_func = gzip.open if file_path.endswith('.gz') else open
    
    try:
        with open_func(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                if line.startswith('>'):
                    # 保存上一个序列
                    if current_id:
                        lengths[current_id] = current_len
                    
                    # 开始新序列
                    # 解析 ID: >UYYG01000001.1 Plasmodium... -> UYYG01000001.1
                    full_header = line[1:]
                    current_id = full_header.split()[0]  # 取空格前第一部分
                    current_len = 0
                else:
                    current_len += len(line)
            
            # 保存最后一个序列
            if current_id:
                lengths[current_id] = current_len
                
    except Exception as e:
        print(f"[ERROR] 读取文件出错 {file_path}: {e}", file=sys.stderr)
        
    return lengths

def main():
    if len(sys.argv) < 3:
        print("Usage: python get_genome_lengths.py <genome_dir_or_file> <output_file>")
        print("Example: python get_genome_lengths.py ./genomes/ all_contig_lengths.txt")
        sys.exit(1)
        
    input_path = sys.argv[1]
    output_file = sys.argv[2]
    
    all_lengths = {}
    files_to_process = []
    
    # 判断输入是文件夹还是文件
    if os.path.isdir(input_path):
        # 递归查找常见后缀
        for ext in ['*.fna']:
            files_to_process.extend(glob.glob(os.path.join(input_path, '**', ext), recursive=True))
    else:
        files_to_process = [input_path]
        
    print(f"[INFO] 找到 {len(files_to_process)} 个基因组文件，开始扫描...")
    
    count = 0
    for fp in files_to_process:
        count += 1
        if count % 10 == 0:
            print(f"Processing ({count}/{len(files_to_process)}): {os.path.basename(fp)} ...")
        
        file_lens = parse_fasta_length(fp)
        # 合并字典 (注意：如果有重名ID，后面会覆盖前面，但在不同物种间应该没问题)
        all_lengths.update(file_lens)
        
    print(f"[INFO] 扫描完成。共统计 {len(all_lengths)} 条 Contig 的长度。")
    print(f"[INFO] 正在写入结果: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("Contig_ID\tLength\n")
        for acc, length in all_lengths.items():
            f.write(f"{acc}\t{length}\n")

if __name__ == "__main__":
    main()
