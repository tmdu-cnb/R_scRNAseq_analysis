#!/usr/bin/env python3
"""
元のResultsファイルからDESeq2の結果データ（baseMean, log2FC, adjPval）を復元するスクリプト
"""

import csv

def load_deseq2_results(original_file):
    """元のResultsファイルからDESeq2の結果データを読み込む"""
    deseq2_data = {}
    
    print(f"Loading DESeq2 results from {original_file}...")
    with open(original_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # 列のインデックスを取得
        user_id_idx = None
        base_mean_idx = None
        log2fc_idx = None
        adj_pval_idx = None
        
        for i, col in enumerate(header):
            col_clean = col.strip('"')
            if col_clean == 'User_ID':
                user_id_idx = i
            elif col_clean == 'baseMean':
                base_mean_idx = i
            elif col_clean == 'CTL-DTA_log2FC':
                log2fc_idx = i
            elif col_clean == 'CTL-DTA_adjPval':
                adj_pval_idx = i
        
        if user_id_idx is None:
            # ensembl_IDを試す
            for i, col in enumerate(header):
                col_clean = col.strip('"')
                if col_clean == 'ensembl_ID':
                    user_id_idx = i
                    break
        
        if user_id_idx is None:
            print("Error: Could not find User_ID or ensembl_ID column")
            return {}
        
        print(f"Found columns: User_ID={user_id_idx}, baseMean={base_mean_idx}, log2FC={log2fc_idx}, adjPval={adj_pval_idx}")
        
        for row in reader:
            if not row or user_id_idx >= len(row):
                continue
            
            gene_id = row[user_id_idx].strip('"')
            if gene_id:
                base_mean = row[base_mean_idx].strip('"') if base_mean_idx and base_mean_idx < len(row) else ''
                log2fc = row[log2fc_idx].strip('"') if log2fc_idx and log2fc_idx < len(row) else ''
                adj_pval = row[adj_pval_idx].strip('"') if adj_pval_idx and adj_pval_idx < len(row) else ''
                
                deseq2_data[gene_id] = {
                    'baseMean': base_mean,
                    'log2FC': log2fc,
                    'adjPval': adj_pval
                }
    
    print(f"Loaded DESeq2 results for {len(deseq2_data)} genes")
    return deseq2_data

def restore_deseq2_to_results(results_file, original_file, output_file):
    """ResultsファイルにDESeq2の結果データを復元"""
    
    # 元のファイルからDESeq2の結果を読み込む
    deseq2_data = load_deseq2_results(original_file)
    
    if not deseq2_data:
        print("Error: Could not load DESeq2 results")
        return
    
    print(f"\nProcessing {results_file}...")
    
    # Resultsファイルを読み込んで更新
    rows_data = []
    header = None
    
    with open(results_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows_data = [row for row in reader if row]
    
    print(f"Found {len(rows_data)} data rows")
    
    # User_ID列のインデックスを探す
    user_id_idx = None
    base_mean_idx = None
    log2fc_idx = None
    adj_pval_idx = None
    
    for i, col in enumerate(header):
        col_clean = col.strip('"')
        if col_clean == 'User_ID':
            user_id_idx = i
        elif col_clean == 'baseMean':
            base_mean_idx = i
        elif col_clean == 'CTL-DTA_log2FC':
            log2fc_idx = i
        elif col_clean == 'CTL-DTA_adjPval':
            adj_pval_idx = i
    
    if user_id_idx is None:
        print("Error: Could not find User_ID column in results file")
        return
    
    print(f"Updating columns: baseMean={base_mean_idx}, log2FC={log2fc_idx}, adjPval={adj_pval_idx}")
    
    # ファイルに書き込む
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        
        writer.writerow(header)
        
        matched_count = 0
        
        for row in rows_data:
            # User_IDを取得
            user_id = row[user_id_idx].strip('"') if user_id_idx < len(row) else ''
            
            # DESeq2の結果データを更新
            if user_id in deseq2_data:
                if base_mean_idx and base_mean_idx < len(row):
                    row[base_mean_idx] = deseq2_data[user_id]['baseMean']
                if log2fc_idx and log2fc_idx < len(row):
                    row[log2fc_idx] = deseq2_data[user_id]['log2FC']
                if adj_pval_idx and adj_pval_idx < len(row):
                    row[adj_pval_idx] = deseq2_data[user_id]['adjPval']
                matched_count += 1
            
            writer.writerow(row)
        
        print(f"Matched and updated {matched_count} genes with DESeq2 results")
    
    print(f"Done! Output written to {output_file}")

if __name__ == '__main__':
    results_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    original_file = '/Users/godakyosuke/Desktop/Results_LFC_Pval_DESeq2.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    
    restore_deseq2_to_results(results_file, original_file, output_file)




