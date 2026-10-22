import sys

# 读取文件A，构建accession对到比对起止位置的映射
def load_accession_to_alignment(file_a):
    accession_to_alignment = {}
    try:
        with open(file_a, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                query_accession = parts[0]
                target_accession = parts[1]
                alignment_start = int(parts[6])  # 修正为第7列
                alignment_end = int(parts[7])    # 修正为第8列

                # 判断起始和终止位置
                if alignment_start > alignment_end:
                    alignment_start, alignment_end = alignment_end, alignment_start

                accession_to_alignment[(query_accession, target_accession)] = (alignment_start, alignment_end)
    except FileNotFoundError:
        print(f"错误: 文件 {file_a} 未找到！")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件 {file_a} 时发生错误: {e}")
        sys.exit(1)
    return accession_to_alignment

# 解析文件D，根据文件D中的物种信息去文件C中查找对应的taxid
def process_file_d(file_d, accession_to_alignment):
    try:
        with open(file_d, 'r', encoding='utf-8') as f:
            for line in f:
                #print(line)
                parts = line.strip().split('\t')
                query_accession = parts[0]
                target_accession = parts[5]
                first_distant_taxid = parts[6]
                first_distant_taxinfo = parts[7]



                if first_distant_taxid:
                    if query_accession and target_accession:
                        # 获取比对的起始和结束位置
                        alignment_start, alignment_end = accession_to_alignment.get((query_accession, target_accession), (None, None))
                        if alignment_start is not None and alignment_end is not None:
                            # 输出最终的比对位置
                            accession = query_accession.split('-')[0]
                            start = int(query_accession.split('-')[1]) + alignment_start - 1
                            end = int(query_accession.split('-')[1]) + alignment_end - 1 
                            print(f'{accession}\t{start}\t{end}\t{first_distant_taxid}\t{first_distant_taxinfo}')
                        else:
                            print(f'未找到比对信息：{query_accession} 与 {target_accession}')
                    else:
                        print(f'未在文件B中找到对应的accession：{first_distant_taxid}')

    except FileNotFoundError:
        print(f"错误: 文件 {file_d} 未找到！")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件 {file_d} 时发生错误: {e}")
        sys.exit(1)

# 主程序执行部分
def main():
    if len(sys.argv) != 3:
        print("错误: 请提供两个文件路径：origin blast和HGT info文件！")
        sys.exit(1)

    file_a = sys.argv[1]
    file_d = sys.argv[2]

    # 加载文件数据
    accession_to_alignment = load_accession_to_alignment(file_a)

    # 处理文件D，输出最终结果
    process_file_d(file_d, accession_to_alignment)

if __name__ == "__main__":
    main()
