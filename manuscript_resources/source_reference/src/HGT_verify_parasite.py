import csv
import sys
import math

def read_taxonomy_file(file_path):
    taxonomy_info = {}
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                taxid = row[0]
                taxonomy = row[1].split(' ', 1)[1] if ' ' in row[1] else ''
                taxonomy_short = row[2]
                taxonomy_info[taxid] = {"detail":taxonomy,"short":taxonomy_short}
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)
    return taxonomy_info


def calculate_entropy(data_dict):
    """
    计算字典的熵值 (单位为 bits)

    :param data_dict: 格式为 {key: value}，value应为非负数
    :return: 熵值 (float)
    """

     # 检查输入是否有效
    if not isinstance(data_dict, dict):
        raise TypeError("输入必须为字典")

    total = sum(data_dict.values())
    #print(total)
    # 处理特殊情况
    if total <= 0:
        total = 1
        #raise ValueError("所有值的总和必须大于零")
    if any(v < 0 for v in data_dict.values()):
        raise ValueError("字典值不能包含负数")

    entropy = 0.0
    for value in data_dict.values():
        probability = value / total
        if probability > 0:  # 忽略零概率项
            entropy -= probability * math.log2(probability)

    return entropy

def analyze_target_file(taxonomy_info, file_path, self_taxid):
    result = {}
    #first we set the distant group as Microbiome and green plants.
    distant_levels = ["Bacteria;","Viruses;",";Fungi;",";Archaea;",";Viridiplantae;"]
    #set their own phylums as self_group.
    if ";Metazoa;" in taxonomy_info.get(self_taxid,'')['detail']:
        self_level = taxonomy_info.get(self_taxid,'')['short'].split(";")[1]
    else:
        self_level = taxonomy_info.get(self_taxid,'')['detail'].split(";")[2]
    #according Metazoa divide into protozoa and metazoa
    if ";Lepeophtheirus salmonis" in taxonomy_info.get(self_taxid,'')['detail']:
        #for Metazoa we use ;Protostomia; as close group, phylum as self_group
        distant_host_levels = [";Actinopterygii"]
    elif ";Varroa destructor" in taxonomy_info.get(self_taxid,'')['detail']:
        distant_host_levels = [";Insecta;"]
    elif ";Aphanomyces astaci" in taxonomy_info.get(self_taxid,'')['detail']:
        distant_host_levels = [";Malacostraca;"]
    elif ";Perkinsus marinus" in taxonomy_info.get(self_taxid,'')['detail']:
        distant_host_levels = [";Bivalvia;"]
    elif ";Leishmania" in taxonomy_info.get(self_taxid,'')['detail']:
        distant_host_levels = [";Diptera;",";Mammalia;"]
    elif any(name in taxonomy_info.get(self_taxid,'')['detail'] for name in [";Angiostrongylus",";Fasciola",";Schistosoma"]):
        distant_host_levels = [";Gastropoda;",";Mammalia;"]
    elif any(name in taxonomy_info.get(self_taxid,'')['detail'] for name in [";Clonorchis",";Opisthorchis",";Anisakis"]):
        distant_host_levels = [";Gastropoda;",";Mammalia;",";Actinopterygii;"]
    elif ";Babesia" in taxonomy_info.get(self_taxid,'')['detail']:
        #for babesia the tick and mammalia are its host.
        distant_host_levels = [";Acari;",";Mammalia;"]
    elif any(name in taxonomy_info.get(self_taxid,'')['detail'] for name in [";Plasmodium",";Brugia",";Wuchereria"]):
        #for Plasmodium, the mosquito and mammals are its host.
        distant_host_levels = [";Culicidae;",";Mammalia;"]
    elif ";Cryptosporidium" in taxonomy_info.get(self_taxid,'')['detail']:
        #for Cryptosporidium, birds and mammals are mian hosts.
        distant_host_levels = [";Aves;",";Mammalia;"]
    else:
        distant_host_levels = [";Mammalia;"]
    #print(close_level)
    print(f"{taxonomy_info.get(self_taxid,'')['detail']}\t{distant_host_levels}")
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                query_accession = row[0]
                target_accession = row[2]
                target_taxid = row[3]
                iden = float(row[4])/100
                length = int(row[5])
                score = float(row[6])
                
                if query_accession not in result:
                    result[query_accession] = {
                        'self_num': 0, 
                        'close_num': 0, 
                        'distant_num': 0, 
                        'processed_taxids': set(),
                        'self_length': -1.0,
                        'first_distant_taxid': '',
                        'first_distant_target_accession': '',
                        'first_distant_taxid_info': '',
                        'first_distant_score': 0.0,
                        'first_distant_length': 0,
                        'first_distant_coverage': 0.0,
                        'first_distant_identity': 0.0,
                        'first_close_info': '',
                        'first_close_score': 0.0,
                        'first_close_length': 0,
                        'first_close_coverage': 0.0,
                        'first_close_identity': 0.0,
                        'alien_index' : 0.0,
                        'distant_host_num' : 0,
                        'distant_other_num' : 0,
                        'class_level_entropy' : 0.0,
                        'class_level_num' : {}
                    }

                if target_taxid in result[query_accession]['processed_taxids']:
                    continue

                result[query_accession]['processed_taxids'].add(target_taxid)

                taxonomy = taxonomy_info.get(target_taxid, '')["detail"]
                taxonomy_short = taxonomy_info.get(target_taxid, '')["short"]
                if taxonomy_short == "":
                    continue
                #print(taxonomy_short)
                #add the class level and num into result[class_level_num] then to calculate the entropy
                if taxonomy_short.split(";")[2] not in result[query_accession]['class_level_num']:
                    result[query_accession]['class_level_num'][taxonomy_short.split(";")[2]] = 1
                else:
                    result[query_accession]['class_level_num'][taxonomy_short.split(";")[2]] += 1

                if f';{self_level};' in taxonomy:
                    result[query_accession]['self_num'] += 1
                    if result[query_accession]['self_num'] == 1:
                        len1 = int(query_accession.split(":")[1].split("-")[0])
                        len2 = int(query_accession.split(":")[1].split("-")[1])
                        result[query_accession]['self_length'] = (len2 - len1)/3
                elif 'organisms;' not in taxonomy:
                    continue
                elif any(distant_host_name in taxonomy for distant_host_name in distant_host_levels):
                    result[query_accession]['distant_host_num'] += 1
                    if result[query_accession]['first_distant_taxid'] == '':
                        result[query_accession]['first_distant_target_accession'] = target_accession
                        result[query_accession]['first_distant_taxid'] = target_taxid
                        result[query_accession]['first_distant_taxid_info'] = taxonomy
                        result[query_accession]['first_distant_score'] = score
                        result[query_accession]['first_distant_length'] = length
                        result[query_accession]['first_distant_identity'] = iden
                elif any(distant_other_name in taxonomy for distant_other_name in distant_levels):
                    result[query_accession]['distant_other_num'] += 1
                    if result[query_accession]['first_distant_taxid'] == '':
                        result[query_accession]['first_distant_target_accession'] = target_accession
                        result[query_accession]['first_distant_taxid'] = target_taxid
                        result[query_accession]['first_distant_taxid_info'] = taxonomy
                        result[query_accession]['first_distant_score'] = score
                        result[query_accession]['first_distant_length'] = length
                        result[query_accession]['first_distant_identity'] = iden
                else:
                    result[query_accession]['close_num'] += 1
                    if result[query_accession]['first_close_info'] == '':
                        #print(score)
                        result[query_accession]['first_close_score'] = score
                        #print(result[query_accession]['first_close_socre'])
                        result[query_accession]['first_close_info'] = taxonomy
                        result[query_accession]['first_close_length'] = length
                        result[query_accession]['first_close_identity'] = iden  

    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
#    except Exception as e:
#        print(f"Error reading file {file_path}: {e}")
#        sys.exit(1)
    
    for query_accession, counts in result.items():
        distant_num = counts['distant_host_num'] + counts['distant_other_num']
        result[query_accession]['distant_num'] = distant_num
        close_num = counts['close_num']
        ratio = distant_num / (distant_num + close_num) if (distant_num + close_num) > 0 else 0
        result[query_accession]['ratio'] = ratio
#        print(counts['first_distant_score'])
#        print(counts['first_close_score'])
        ai = counts['first_distant_score'] - counts['first_close_score']
        result[query_accession]['alien_index'] = ai
        result[query_accession]['first_distant_coverage'] = result[query_accession]['first_distant_length']/result[query_accession]['self_length']
        result[query_accession]['first_close_coverage'] = result[query_accession]['first_close_length']/result[query_accession]['self_length']
        # calculate the entrofy
        entropy = calculate_entropy(counts['class_level_num'])
        result[query_accession]['class_level_entropy'] = entropy
    
    return result

# 调用示例
if len(sys.argv) != 4:
    print("Usage: python script.py <taxonomy_file> <target_file> <self_taxid>")
    sys.exit(1)

taxonomy_file = sys.argv[1]
target_file = sys.argv[2]
self_taxid = sys.argv[3]


taxonomy_info = read_taxonomy_file(taxonomy_file)
result = analyze_target_file(taxonomy_info, target_file, self_taxid)

# 输出结果
print(f"Accession\tOut_pct\tSelf_num\tClose_num\tDistant_num\tFisrt_distant_target_accession\tFirst_distant_taxid\tFirst_distant_taxinfo\tFirst_distant_coverage\tFirst_distant_identity\talien_index\tdistant_host_num\tdistant_other_num\tfirst_close_taxonomy_info\tfirst_close_coverage\tfirst_close_identity\tentropy")
for query_accession, counts in result.items():
    print(f"{query_accession}\t{counts['ratio']:.3f}\t{counts['self_num']}\t{counts['close_num']}\t{counts['distant_num']}\t{counts['first_distant_target_accession']}\t{counts['first_distant_taxid']}\t{counts['first_distant_taxid_info']}\t{counts['first_distant_coverage']}\t{counts['first_distant_identity']}\t{counts['alien_index']}\t{counts['distant_host_num']}\t{counts['distant_other_num']}\t{counts['first_close_info']}\t{counts['first_close_coverage']}\t{counts['first_close_identity']}\t{counts['class_level_entropy']}")