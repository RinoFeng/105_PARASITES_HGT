import sys
import argparse
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import os

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='HGT候选基因筛选工具',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', 
                      help='输入文件路径，格式要求：第二行为列头，空格分隔')
    parser.add_argument('-o', '--output', 
                      default='hgt_output.csv',
                      help='输出文件路径')
    return parser.parse_args()

def read_data(file_path):
    """读取数据并进行基本校验"""
    try:
        # 检查文件是否存在
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"输入文件不存在：{file_path}")
        
        # 读取数据
        df = pd.read_csv(
            file_path,
            sep='\t',
            skiprows=1,
            header=0,
            engine='python',
            na_values=['NA', '-', '']
        )
        
        # 检查必要列是否存在
        required_columns = ['entropy', 'Out_pct','Self_num', 'Distant_num', 
                          'First_distant_taxinfo','First_distant_coverage',
                            'First_distant_identity','alien_index','first_close_coverage']
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            raise ValueError(f"缺少必要列：{missing_cols}")

        # 转换数值类型
        numeric_cols = ['Out_pct', 'Distant_num', 'entropy','Self_num','First_distant_coverage','first_close_coverage','First_distant_identity','alien_index']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            if df[col].isnull().all():
                raise ValueError(f"列 {col} 包含无法转换的数值")

        return df.dropna(subset=numeric_cols)
    
    except pd.errors.ParserError as e:
        sys.exit(f"文件解析失败：{str(e)}\n请检查文件分隔符是否正确")

def find_valley_threshold(series, width_percent=0.05):
    """带错误处理的密度谷计算"""
    try:
        # 过滤无效值
        valid_series = series.dropna()
        if len(valid_series) < 10:
            raise ValueError("有效熵值数据不足(<10条)，无法计算密度估计")
        
        kde = gaussian_kde(valid_series)
        x = np.linspace(valid_series.min(), valid_series.max(), 1024)
        y = kde(x)
        
        # 寻找密度谷
        valley_idx = None
        for i in range(1, len(y)-1):
            if y[i-1] > y[i] < y[i+1]:
                valley_idx = i
                break
        
        if valley_idx is None:
            print("未检测到密度谷，可能为单峰分布,使用极端entropy的0.2宽度作为阈值。")
            return 0.2 * (x[-1] - x[0])
            
        base = x[valley_idx]
        threshold = base + width_percent * (x[-1] - x[0])
        return threshold
    
    except np.linalg.LinAlgError:
        sys.exit("密度估计失败：数据可能包含重复值或分布异常")

def main():
    args = parse_arguments()
    
    try:
        # 读取数据
        df = read_data(args.input_file)
        if df.empty:
            sys.exit("错误：数据读取后为空，请检查输入文件格式")
        
        # 新增步骤：过滤Distant_num=0的数据用于计算阈值
        filtered_df = df[df['Distant_num'] > 0]
        if len(filtered_df) < 10:
            sys.exit("错误：Distant_num>0的有效数据不足10条，无法计算密度谷")

        # 使用过滤后的数据计算阈值
        entropy_series = filtered_df['entropy']
        threshold = find_valley_threshold(entropy_series)
        print(f"[计算完成] 基于Distant_num>0数据的熵阈值：{threshold:.4f}")
    
        # 构建筛选条件
        cond_normal = (
            (df['entropy'] < threshold) &
            (df['Out_pct'] > 0.9) &
            (df['Distant_num'] > 1) &
            (df['Distant_num'] >= df['Self_num']) &
            (df['Self_num'] > 0) &
            (df['First_distant_coverage'] >= 0.7) &
            (df['First_distant_identity'] >= 0.3) &
            (df['alien_index'] > 0)
        )
        
        cond_special = (
            df['First_distant_taxinfo'].str.contains(
                'Bacteria', case=False, na=False) &
            (df['Out_pct'] > 0.9) &
            (df['Distant_num'] > 1) &
            (df['Distant_num'] >= df['Self_num']) &
            (df['Self_num'] > 0) &
            (df['alien_index'] > 0) &
            (df['First_distant_coverage'] >= 0.5) &
            (df['First_distant_identity'] >= 0.3) 
        )
        
        result = df[cond_normal | cond_special]
        # 保存结果
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        result.to_csv(args.output, index=False,sep = "\t")
        print(f"结果已保存至：{args.output}")
        print(f"匹配记录数：{len(result)}/{len(df)}")
        
    except Exception as e:
        sys.exit(f"程序运行出错：{str(e)}")

if __name__ == "__main__":
    main()
