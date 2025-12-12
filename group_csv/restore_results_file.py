#!/usr/bin/env python3
"""
Results_LFC_Pval_DESeq2.csvファイルを復元するスクリプト
DCとDTのカウントデータから構造を復元
"""

import csv

def load_count_data(count_file):
    """カウントファイルから遺伝子IDをキーとしてカウントデータを読み込む"""
    count_data = {}
    gene_names = {}
    
    print(f"Loading count data from {count_file}...")
    with open(count_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # カウントデータ列を取得
        count_cols = [col.strip('"') for col in header[2:] if col.strip()]
        
        for row in reader:
            if not row or not row[0]:
                continue
            
            gene_id = row[0].strip()
            if gene_id:
                if len(row) > 1:
                    gene_names[gene_id] = row[1].strip()
                counts = row[2:2+len(count_cols)]
                count_data[gene_id] = counts
    
    print(f"Loaded {len(count_data)} genes with {len(count_cols)} count columns")
    return count_data, count_cols, gene_names

def load_gene_info(id_gene_file):
    """idとgene.csvから遺伝子情報を読み込む"""
    gene_info = {}
    print(f"Loading gene information from {id_gene_file}...")
    with open(id_gene_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row['User_ID']
            gene_name = row['symbol']
            gene_info[gene_id] = gene_name
    
    print(f"Loaded {len(gene_info)} gene mappings")
    return gene_info

def restore_results_file(dc_file, dt_file, id_gene_file, output_file):
    """Resultsファイルを復元"""
    
    # DCとDTのカウントデータを読み込む
    dc_data, dc_cols, dc_gene_names = load_count_data(dc_file)
    dt_data, dt_cols, dt_gene_names = load_count_data(dt_file)
    
    # 遺伝子情報を読み込む
    gene_info = load_gene_info(id_gene_file)
    
    # すべての遺伝子IDを取得（ソート済み）
    all_gene_ids = sorted(set(dc_data.keys()) | set(dt_data.keys()))
    print(f"\nTotal unique genes: {len(all_gene_ids)}")
    
    # ヘッダーを作成（最初に見た構造に合わせる）
    header = [
        'symbol', 'ensembl_ID', 'User_ID',
        'baseMean', 'CTL-DTA_log2FC', 'CTL-DTA_adjPval',
        'Processed data:', 'search_label'
    ] + dc_cols + dt_cols
    
    # ファイルに書き込む
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        
        writer.writerow(header)
        
        for gene_id in all_gene_ids:
            # 遺伝子名を取得（優先順位: DC > DT > idとgene.csv）
            gene_name = dc_gene_names.get(gene_id) or dt_gene_names.get(gene_id) or gene_info.get(gene_id, '')
            
            # DCカウントデータを取得
            dc_counts = dc_data.get(gene_id, ['0'] * len(dc_cols))
            while len(dc_counts) < len(dc_cols):
                dc_counts.append('0')
            
            # DTカウントデータを取得
            dt_counts = dt_data.get(gene_id, ['0'] * len(dt_cols))
            while len(dt_counts) < len(dt_cols):
                dt_counts.append('0')
            
            # search_labelを作成
            search_label = f"{gene_name} | {gene_id}" if gene_name else gene_id
            
            # 行を作成
            # DESeq2の結果列は空にする（計算できないため）
            row = [
                gene_name,           # symbol
                gene_id,             # ensembl_ID
                gene_id,             # User_ID
                '',                  # baseMean (空)
                '',                  # CTL-DTA_log2FC (空)
                '',                  # CTL-DTA_adjPval (空)
                '',                  # Processed data: (空)
                search_label         # search_label
            ] + dc_counts[:len(dc_cols)] + dt_counts[:len(dt_cols)]
            
            writer.writerow(row)
    
    print(f"Done! Restored {len(all_gene_ids)} genes to {output_file}")
    print("Note: DESeq2 result columns (baseMean, log2FC, adjPval) are empty and need to be calculated separately")

if __name__ == '__main__':
    dc_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    dt_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    id_gene_file = '/Users/godakyosuke/Desktop/idとgene.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    
    restore_results_file(dc_file, dt_file, id_gene_file, output_file)




