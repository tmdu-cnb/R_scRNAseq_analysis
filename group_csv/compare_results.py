#!/usr/bin/env python3
"""
2つのDEG結果ファイルを比較するスクリプト
"""

import csv

def load_idep_results(idep_file):
    """iDEPの結果を読み込む"""
    results = {}
    print(f"Loading iDEP results from {idep_file}...")
    with open(idep_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row['User_ID'].strip()
            if gene_id:
                results[gene_id] = {
                    'symbol': row.get('symbol', ''),
                    'baseMean': row.get('baseMean', ''),
                    'log2FC': row.get('CTL-DTA_log2FC', ''),
                    'adjPval': row.get('CTL-DTA_adjPval', '')
                }
    print(f"Loaded {len(results)} genes from iDEP")
    return results

def load_loupe_results(loupe_file):
    """ルーペブラウザーの結果を読み込む"""
    results = {}
    print(f"Loading Loupe results from {loupe_file}...")
    with open(loupe_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row['FeatureID'].strip()
            if gene_id:
                results[gene_id] = {
                    'symbol': row.get('FeatureName', ''),
                    'average': row.get('dtcl calca8 Average', ''),
                    'log2FC': row.get('dtcl calca8 Log2 Fold Change', ''),
                    'pval': row.get('dtcl calca8 P-Value', '')
                }
    print(f"Loaded {len(results)} genes from Loupe")
    return results

def compare_results(idep_file, loupe_file, output_file):
    """2つの結果を比較"""
    
    idep_results = load_idep_results(idep_file)
    loupe_results = load_loupe_results(loupe_file)
    
    # 共通の遺伝子を取得
    common_genes = set(idep_results.keys()) & set(loupe_results.keys())
    print(f"\nCommon genes: {len(common_genes)}")
    
    # 比較結果を書き込む
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        
        header = [
            'Gene_ID', 'Symbol',
            'iDEP_baseMean', 'iDEP_log2FC', 'iDEP_adjPval',
            'Loupe_Average', 'Loupe_log2FC', 'Loupe_Pval',
            'log2FC_diff', 'abs_log2FC_diff'
        ]
        writer.writerow(header)
        
        large_diff_count = 0
        for gene_id in sorted(common_genes):
            idep = idep_results[gene_id]
            loupe = loupe_results[gene_id]
            
            # log2FCの差を計算
            try:
                idep_fc = float(idep['log2FC']) if idep['log2FC'] and idep['log2FC'] != 'NA' else 0
                loupe_fc = float(loupe['log2FC']) if loupe['log2FC'] else 0
                fc_diff = abs(idep_fc - loupe_fc)
                
                if fc_diff > 1.0:  # log2FCの差が1以上
                    large_diff_count += 1
            except:
                fc_diff = 'N/A'
                abs_fc_diff = 'N/A'
            
            row = [
                gene_id,
                idep['symbol'] or loupe['symbol'],
                idep['baseMean'],
                idep['log2FC'],
                idep['adjPval'],
                loupe['average'],
                loupe['log2FC'],
                loupe['pval'],
                idep_fc - loupe_fc if 'idep_fc' in locals() else 'N/A',
                fc_diff if 'fc_diff' in locals() and fc_diff != 'N/A' else 'N/A'
            ]
            writer.writerow(row)
        
        print(f"\nGenes with large log2FC difference (>1.0): {large_diff_count}")
    
    print(f"\nComparison results written to {output_file}")

if __name__ == '__main__':
    idep_file = '/Users/godakyosuke/Desktop/group_csv/scRNAseqからidepでDEG.csv'
    loupe_file = '/Users/godakyosuke/Desktop/group_csv/dt calca8 Features.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DEG_comparison.csv'
    
    compare_results(idep_file, loupe_file, output_file)


