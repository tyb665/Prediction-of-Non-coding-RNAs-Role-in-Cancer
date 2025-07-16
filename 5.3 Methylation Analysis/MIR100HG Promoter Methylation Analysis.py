#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os

# ======== 配置路径（当前目录下）========
base_dir = "./"  # 所有文件都在当前目录下
probe_file = os.path.join(base_dir, "Probe_Map.csv")
promoter_file = os.path.join(base_dir, "Promoter_Region.csv")
output_file = os.path.join(base_dir, "MIR100HG_Promoter_Probes.csv")

# ======== 读取数据 ========
probe_df = pd.read_csv(probe_file)
promoter_df = pd.read_csv(promoter_file)

# ======== 提取 MIR100HG 启动子区域 ========
mir100hg_promoter = promoter_df[promoter_df['HGNC_Symbol'] == 'MIR100HG']

if mir100hg_promoter.empty:
    print("⚠️ 未找到 MIR100HG 的启动子区域，请检查 Promoter_Region.csv 是否包含该基因。")
else:
    matched_probes = []

    for _, prom_row in mir100hg_promoter.iterrows():
        chr_match = probe_df['Probe_Chr'] == prom_row['Promoter_Chr']
        start_in_range = probe_df['Probe_Chrom_Start'] >= prom_row['Promoter_Start']
        end_in_range = probe_df['Probe_Chrom_End'] <= prom_row['Promoter_End']

        probes_in_promoter = probe_df[chr_match & start_in_range & end_in_range].copy()
        probes_in_promoter["HGNC_Symbol"] = prom_row['HGNC_Symbol']
        matched_probes.append(probes_in_promoter)

    # 合并结果
    if matched_probes:
        result_df = pd.concat(matched_probes, ignore_index=True)
        result_df = result_df[['Probe_ID', 'HGNC_Symbol', 'Probe_Chr', 'Probe_Chrom_Start', 'Probe_Chrom_End']]
        result_df.to_csv(output_file, index=False)
        print(f"✅ 已成功筛选 MIR100HG 启动子区域探针，共 {len(result_df)} 条")
        print(f"📁 结果保存至：{output_file}")
        display(result_df.head())  # Jupyter 环境下直接展示前几行
    else:
        print("⚠️ 未找到任何探针落入 MIR100HG 启动子区域。")


# In[10]:


import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import re

# 定义要分析的癌症类型
cancer_types = ["STAD", "LUAD", "PRAD", "PAAD", "SKCM"]

# 配置路径
base_dir = "./"
# 为不同癌症类型定义特定的文件名模式
methylation_file_patterns = {
    "PRAD": ["{}_Model_Methylation_Features_Test.csv"],
    "SKCM": ["{}_Model_Methylation_Features_Test.csv"],
    "default": [
        "{}_Model_Methylation_Features.csv",
        "{}_Methylation_Probes_in_Promoters.csv",
        "{}_Methylation_Features.csv",
        "{}_DMA_Significant_Results.csv"
    ]
}
expression_file_pattern = "{}_Model_MIR100HG_Expression_Levels.csv"
output_file_pattern = "{}_MIR100HG_Methylation_Correlation.csv"
promoter_probes_file = os.path.join(base_dir, "MIR100HG_Promoter_Probes.csv")

def benjamini_hochberg_correction(p_values):
    """
    实现Benjamini-Hochberg FDR校正方法，不依赖statsmodels
    """
    n = len(p_values)
    if n == 0:
        return np.array([])
        
    # 将p值与其索引一起排序
    indices = np.arange(n)
    sorted_data = sorted(zip(p_values, indices))
    sorted_p, sorted_indices = zip(*sorted_data)
    
    # 计算校正后的p值
    sorted_adjusted_p = np.zeros(n)
    prev_adjusted_p = 1.0
    
    for i in range(n-1, -1, -1):
        adjusted_p = min(prev_adjusted_p, sorted_p[i] * n / (i + 1))
        sorted_adjusted_p[i] = adjusted_p
        prev_adjusted_p = adjusted_p
    
    # 将校正后的p值重新排序回原始顺序
    adjusted_p = np.zeros(n)
    for i, idx in enumerate(sorted_indices):
        adjusted_p[idx] = sorted_adjusted_p[i]
    
    return adjusted_p

def load_expression_data(cancer_type):
    """
    加载MIR100HG表达数据
    """
    file_path = os.path.join(base_dir, expression_file_pattern.format(cancer_type))
    
    try:
        # 读取CSV文件
        df = pd.read_csv(file_path)
        
        # 检查必要的列是否存在
        if 'Sample_ID' not in df.columns:
            # 尝试查找类似的列名
            sample_cols = [col for col in df.columns if 'sample' in col.lower() or 'id' in col.lower()]
            if sample_cols:
                print(f"使用'{sample_cols[0]}'作为样本ID列")
                df.rename(columns={sample_cols[0]: 'Sample_ID'}, inplace=True)
            else:
                print(f"错误：{cancer_type}表达数据中找不到样本ID列")
                print(f"可用列：{df.columns.tolist()}")
                return None
                
        if 'MIR100HG_Expression' not in df.columns:
            # 尝试查找类似的列名
            expr_cols = [col for col in df.columns if 'mir100hg' in col.lower() or 'expression' in col.lower()]
            if expr_cols:
                print(f"使用'{expr_cols[0]}'作为表达值列")
                df.rename(columns={expr_cols[0]: 'MIR100HG_Expression'}, inplace=True)
            else:
                print(f"错误：{cancer_type}表达数据中找不到MIR100HG表达列")
                print(f"可用列：{df.columns.tolist()}")
                return None
        
        # 创建以Sample_ID为索引的DataFrame
        expr_df = df[['Sample_ID', 'MIR100HG_Expression']].set_index('Sample_ID')
        
        # 确保表达值为数值类型
        expr_df['MIR100HG_Expression'] = pd.to_numeric(expr_df['MIR100HG_Expression'], errors='coerce')
        
        # 删除缺失值
        expr_df.dropna(inplace=True)
        
        print(f"加载了{cancer_type}的{len(expr_df)}个样本的表达数据")
        return expr_df
        
    except Exception as e:
        print(f"加载{cancer_type}的表达数据时出错: {str(e)}")
        return None

def find_methylation_file(cancer_type):
    """
    查找当前目录中所有可能的甲基化数据文件
    """
    # 获取当前癌症类型的文件名模式
    if cancer_type in methylation_file_patterns:
        patterns = methylation_file_patterns[cancer_type]
    else:
        patterns = methylation_file_patterns["default"]
    
    # 尝试所有可能的文件名模式
    for pattern in patterns:
        file_path = os.path.join(base_dir, pattern.format(cancer_type))
        print(f"尝试查找文件: {file_path}")
        if os.path.exists(file_path):
            print(f"找到甲基化数据文件: {file_path}")
            return file_path
    
    # 如果没有找到精确匹配，尝试部分匹配
    all_files = os.listdir(base_dir)
    for file in all_files:
        if cancer_type in file and ('methylation' in file.lower() or 'meth' in file.lower() or 'dma' in file.lower()):
            file_path = os.path.join(base_dir, file)
            print(f"找到可能的甲基化数据文件: {file_path}")
            return file_path
    
    print(f"未找到{cancer_type}的甲基化数据文件")
    return None

def find_probe_id_column(df):
    """
    在DataFrame中查找可能的探针ID列
    """
    # 可能的探针ID列名
    probe_cols = ['Probe_ID', 'ProbeID', 'Probe', 'ID', 'probe_id', 'probe']
    
    # 检查是否有完全匹配的列名
    for col in probe_cols:
        if col in df.columns:
            return col
    
    # 检查是否有部分匹配的列名
    for col in df.columns:
        if any(probe_term in col.lower() for probe_term in ['probe', 'cg']):
            return col
    
    # 如果第一列看起来像探针ID（例如，以cg开头），使用第一列
    first_col = df.columns[0]
    if isinstance(df[first_col].iloc[0], str) and (df[first_col].iloc[0].startswith('cg') or re.match(r'cg\d+', df[first_col].iloc[0])):
        return first_col
    
    return None

def load_methylation_data(cancer_type, promoter_probes):
    """
    加载甲基化数据，改进的探针匹配策略
    """
    # 查找甲基化数据文件
    file_path = find_methylation_file(cancer_type)
    if not file_path:
        return None
    
    try:
        # 读取甲基化数据
        meth_data = pd.read_csv(file_path)
        
        # 查找探针ID列
        probe_col = find_probe_id_column(meth_data)
        if not probe_col:
            print(f"错误：无法在{cancer_type}的甲基化数据中找到探针ID列")
            print(f"可用列：{meth_data.columns.tolist()[:5]}...")
            return None
        
        print(f"使用'{probe_col}'列作为探针ID")
        
        # 设置探针ID为索引
        meth_data.set_index(probe_col, inplace=True)
        
        # 打印一些索引信息以供参考
        print(f"甲基化数据索引样例：{list(meth_data.index)[:5]}...")
        print(f"启动子探针样例：{promoter_probes[:5]}...")
        
        # 统计数据框的形状
        original_shape = meth_data.shape
        print(f"原始甲基化数据形状: {original_shape}")
        
        # 更灵活的探针匹配方式
        # 1. 精确匹配
        valid_probes_exact = [p for p in promoter_probes if p in meth_data.index]
        
        # 2. 不区分大小写匹配
        if not valid_probes_exact:
            valid_probes_case = []
            meth_index_lower = [idx.lower() if isinstance(idx, str) else str(idx).lower() for idx in meth_data.index]
            for p in promoter_probes:
                p_lower = p.lower()
                if p_lower in meth_index_lower:
                    orig_idx = meth_data.index[meth_index_lower.index(p_lower)]
                    valid_probes_case.append(orig_idx)
            if valid_probes_case:
                print(f"通过不区分大小写匹配找到{len(valid_probes_case)}个探针")
                valid_probes_exact = valid_probes_case
        
        # 3. 部分匹配（如果精确匹配和不区分大小写都没有结果）
        if not valid_probes_exact:
            valid_probes_partial = []
            for p in promoter_probes:
                for idx in meth_data.index:
                    if isinstance(idx, str) and p in idx:
                        valid_probes_partial.append(idx)
            if valid_probes_partial:
                print(f"通过部分匹配找到{len(valid_probes_partial)}个探针")
                valid_probes_exact = valid_probes_partial
        
        if not valid_probes_exact:
            print(f"警告：{cancer_type}的甲基化数据中没有找到任何启动子区域探针")
            return None
        
        print(f"找到{len(valid_probes_exact)}个有效的启动子区域探针")
        meth_data = meth_data.loc[valid_probes_exact]
        
        # 检查甲基化数据的形状
        print(f"过滤后的甲基化数据形状: {meth_data.shape}")
        
        # 确保所有列都是样本数据，而不是元数据
        # 假设样本ID以TCGA_开头
        sample_cols = [col for col in meth_data.columns if isinstance(col, str) and ('TCGA_' in col or col.startswith('TCGA-'))]
        if not sample_cols:
            # 如果找不到TCGA样本，尝试使用所有看起来像数值的列
            sample_cols = [col for col in meth_data.columns if col != probe_col and pd.api.types.is_numeric_dtype(meth_data[col])]
        
        if len(sample_cols) < len(meth_data.columns):
            print(f"移除{len(meth_data.columns) - len(sample_cols)}个非样本列")
            meth_data = meth_data[sample_cols]
        
        # 确保所有值都是数值类型
        for col in meth_data.columns:
            meth_data[col] = pd.to_numeric(meth_data[col], errors='coerce')
        
        # 删除全是NaN的行和列
        meth_data.dropna(how='all', axis=0, inplace=True)  # 删除全是NaN的行
        meth_data.dropna(how='all', axis=1, inplace=True)  # 删除全是NaN的列
        
        print(f"最终甲基化数据形状: {meth_data.shape}")
        print(f"加载了{len(meth_data)}个探针和{len(meth_data.columns)}个样本的甲基化数据")
        return meth_data
        
    except Exception as e:
        print(f"加载{cancer_type}的甲基化数据时出错: {str(e)}")
        import traceback
        print(traceback.format_exc())
        return None

def get_promoter_probes():
    """
    获取MIR100HG启动子区域的探针列表
    """
    try:
        # 尝试从文件加载探针列表
        probes_df = pd.read_csv(promoter_probes_file)
        
        # 检查必要的列是否存在
        if 'Probe_ID' not in probes_df.columns:
            # 尝试查找类似的列
            probe_cols = [col for col in probes_df.columns if 'probe' in col.lower() or 'id' in col.lower() or 'cg' in col.lower()]
            if probe_cols:
                print(f"使用'{probe_cols[0]}'作为探针ID列")
                probes_df.rename(columns={probe_cols[0]: 'Probe_ID'}, inplace=True)
            else:
                print(f"错误：启动子探针文件中找不到探针ID列")
                print(f"可用列：{probes_df.columns.tolist()}")
                return None
        
        # 检查是否有重复的探针
        unique_probes = probes_df['Probe_ID'].unique()
        if len(unique_probes) < len(probes_df):
            print(f"警告：探针列表中有重复项，原始数量：{len(probes_df)}，去重后：{len(unique_probes)}")
        
        # 使用去重后的探针列表
        promoter_probes = unique_probes.tolist()
        print(f"从文件中加载了{len(promoter_probes)}个启动子探针")
        return promoter_probes
        
    except Exception as e:
        print(f"读取启动子探针文件时出错: {str(e)}")
        return None

def analyze_single_cancer(cancer_type):
    """
    为单个癌症类型进行MIR100HG启动子甲基化与表达相关性分析
    """
    print(f"\n{'='*50}")
    print(f"分析 {cancer_type} 中MIR100HG启动子甲基化与表达的相关性")
    print(f"{'='*50}\n")
    
    # 获取启动子探针列表
    promoter_probes = get_promoter_probes()
    
    if promoter_probes is None or len(promoter_probes) == 0:
        print(f"错误: 无法获取启动子探针列表，分析终止")
        return
    
    # 加载表达数据
    expression = load_expression_data(cancer_type)
    if expression is None or len(expression) == 0:
        print(f"无法进行{cancer_type}的分析: 缺少表达数据")
        return
    
    # 加载甲基化数据
    methylation = load_methylation_data(cancer_type, promoter_probes)
    if methylation is None or methylation.empty:
        print(f"无法进行{cancer_type}的分析: 缺少甲基化数据")
        return
    
    # 找出表达数据和甲基化数据中共同的样本
    common_samples = list(set(expression.index).intersection(set(methylation.columns)))
    
    if len(common_samples) == 0:
        print(f"无法进行{cancer_type}的分析: 表达数据和甲基化数据没有共同的样本")
        # 尝试部分匹配样本ID
        expr_samples = [s for s in expression.index if isinstance(s, str)]
        meth_samples = [s for s in methylation.columns if isinstance(s, str)]
        
        common_partial = []
        for e_sample in expr_samples:
            for m_sample in meth_samples:
                # 移除可能的前缀/后缀，只比较核心部分
                e_core = re.sub(r'[^a-zA-Z0-9]', '', e_sample)
                m_core = re.sub(r'[^a-zA-Z0-9]', '', m_sample)
                if e_core in m_core or m_core in e_core:
                    common_partial.append((e_sample, m_sample))
        
        if common_partial:
            print(f"通过部分匹配找到{len(common_partial)}个共同样本")
            # 创建映射并重命名甲基化数据列名以匹配表达数据
            rename_dict = {m: e for e, m in common_partial}
            methylation.rename(columns=rename_dict, inplace=True)
            # 重新查找共同样本
            common_samples = list(set(expression.index).intersection(set(methylation.columns)))
            print(f"重命名后找到{len(common_samples)}个共同样本")
        else:
            print("即使通过部分匹配也找不到共同样本")
            return
    
    print(f"找到{len(common_samples)}个共同样本")
    
    # 筛选只包含共同样本的数据
    expression_filtered = expression.loc[common_samples]
    methylation_filtered = methylation[common_samples]
    
    # 打印数据形状信息
    print(f"表达数据形状: {expression_filtered.shape}")
    print(f"甲基化数据形状: {methylation_filtered.shape}")
    
    # 计算每个探针与表达的相关性
    results = []
    
    for probe in methylation_filtered.index:
        # 获取探针的甲基化值
        probe_values = methylation_filtered.loc[probe]
        expr_values = expression_filtered['MIR100HG_Expression']
        
        # 移除任何NaN值
        valid_indices = ~(probe_values.isna() | expr_values.isna())
        valid_probe_values = probe_values[valid_indices]
        valid_expr_values = expr_values[valid_indices]
        
        # 确保有足够的有效样本进行相关性分析
        if len(valid_probe_values) < 10:
            print(f"警告: 探针{probe}只有{len(valid_probe_values)}个有效样本，跳过")
            continue
            
        # 计算Spearman相关系数
        try:
            rho, p_value = spearmanr(valid_probe_values, valid_expr_values)
            
            # 检查结果是否有效
            if np.isnan(rho) or np.isnan(p_value):
                print(f"警告: 探针{probe}的计算结果为NaN，跳过")
                continue
                
            results.append({
                'Probe_ID': probe,
                'Spearman_rho': rho,
                'P_value': p_value,
                'Valid_samples': len(valid_probe_values)
            })
        except Exception as e:
            print(f"计算探针{probe}的相关性时出错: {str(e)}")
    
    # 如果没有有效的结果，则终止
    if not results:
        print(f"未能为{cancer_type}计算出任何有效的相关性结果")
        return
    
    # 创建结果DataFrame
    results_df = pd.DataFrame(results)
    
    # 添加FDR校正
    results_df['FDR_adjusted_P'] = benjamini_hochberg_correction(results_df['P_value'].values)
    
    # 添加显著性标记
    results_df['Significant_P05'] = results_df['P_value'] < 0.05
    results_df['Significant_FDR05'] = results_df['FDR_adjusted_P'] < 0.05
    
    # 按相关系数排序（最负相关的在前）
    results_df = results_df.sort_values('Spearman_rho')
    
    # 保存结果
    output_file = os.path.join(base_dir, output_file_pattern.format(cancer_type))
    results_df.to_csv(output_file, index=False)
    print(f"结果已保存至: {output_file}")
    
    # 打印相关性摘要
    total_probes = len(results_df)
    neg_corr = sum(results_df['Spearman_rho'] < 0)
    pos_corr = sum(results_df['Spearman_rho'] > 0)
    
    sig_neg = sum((results_df['Spearman_rho'] < 0) & results_df['Significant_P05'])
    sig_pos = sum((results_df['Spearman_rho'] > 0) & results_df['Significant_P05'])
    
    fdr_sig_neg = sum((results_df['Spearman_rho'] < 0) & results_df['Significant_FDR05'])
    fdr_sig_pos = sum((results_df['Spearman_rho'] > 0) & results_df['Significant_FDR05'])
    
    print(f"\n{cancer_type}相关性分析摘要:")
    print(f"  总探针数: {total_probes}")
    print(f"  负相关探针: {neg_corr} ({neg_corr/total_probes*100:.1f}%)")
    print(f"    - p<0.05显著: {sig_neg} ({sig_neg/total_probes*100:.1f}%)")
    print(f"    - FDR<0.05显著: {fdr_sig_neg} ({fdr_sig_neg/total_probes*100:.1f}%)")
    print(f"  正相关探针: {pos_corr} ({pos_corr/total_probes*100:.1f}%)")
    print(f"    - p<0.05显著: {sig_pos} ({sig_pos/total_probes*100:.1f}%)")
    print(f"    - FDR<0.05显著: {fdr_sig_pos} ({fdr_sig_pos/total_probes*100:.1f}%)")
    
    # 打印最显著的负相关探针
    if sig_neg > 0:
        print("\n  显著负相关的探针 (p<0.05):")
        top_neg = results_df[(results_df['Spearman_rho'] < 0) & results_df['Significant_P05']].head(10)
        for i, (_, row) in enumerate(top_neg.iterrows(), 1):
            fdr_mark = "**" if row['Significant_FDR05'] else ""
            print(f"    {i}. {row['Probe_ID']}: rho = {row['Spearman_rho']:.4f}, p = {row['P_value']:.4e}{fdr_mark}")
    
    # 打印最显著的正相关探针
    if sig_pos > 0:
        print("\n  显著正相关的探针 (p<0.05):")
        top_pos = results_df[(results_df['Spearman_rho'] > 0) & results_df['Significant_P05']].nsmallest(10, 'P_value')
        for i, (_, row) in enumerate(top_pos.iterrows(), 1):
            fdr_mark = "**" if row['Significant_FDR05'] else ""
            print(f"    {i}. {row['Probe_ID']}: rho = {row['Spearman_rho']:.4f}, p = {row['P_value']:.4e}{fdr_mark}")
    
    print(f"\n{cancer_type}分析完成!\n")
    return results_df

def main():
    """
    主函数，为每个癌症类型单独进行分析
    """
    print("==== MIR100HG启动子甲基化与表达相关性分析 ====\n")
    
    # 分别分析每个癌症类型
    for cancer_type in cancer_types:
        # 为该癌症类型进行单独分析
        analyze_single_cancer(cancer_type)
    
    print("\n全部分析完成!")

if __name__ == "__main__":
    main()

