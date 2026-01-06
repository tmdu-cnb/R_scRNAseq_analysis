#!/usr/bin/env python3
"""
DT_PURKINJE_counts.csvに遺伝子IDを追加するスクリプト
"""

import csv

def load_gene_mapping(mapping_file):
    """idとgene.csvから遺伝子名→ENSMUSG番号のマッピングを読み込む"""
    gene_to_id = {}
    
    print(f"Loading gene mapping from {mapping_file}...")
    with open(mapping_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_name = row['symbol']
            gene_id = row['User_ID']
            # 大文字小文字を区別せずにマッピング（念のため）
            gene_to_id[gene_name] = gene_id
            gene_to_id[gene_name.upper()] = gene_id
            gene_to_id[gene_name.lower()] = gene_id
    
    print(f"Loaded {len(gene_to_id)} gene mappings")
    return gene_to_id

def add_gene_ids(input_file, mapping_file, output_file):
    """DT_PURKINJE_counts.csvに遺伝子IDを追加"""
    
    # 遺伝子マッピングを読み込む
    gene_to_id = load_gene_mapping(mapping_file)
    
    print(f"\nProcessing {input_file}...")
    
    matched_count = 0
    unmatched_count = 0
    unmatched_genes = []
    
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # ヘッダー行を読み込む
        header = next(reader)
        
        # 新しいヘッダーを作成（最初にgene_idを追加）
        new_header = ['gene_id'] + header
        writer.writerow(new_header)
        
        # データ行を処理
        for row in reader:
            if not row:
                continue
            
            gene_name = row[0].strip('"')  # 遺伝子名（最初の列、クォートを除去）
            
            # 遺伝子IDを検索
            gene_id = gene_to_id.get(gene_name)
            
            if gene_id:
                matched_count += 1
            else:
                unmatched_count += 1
                unmatched_genes.append(gene_name)
                gene_id = ''  # 見つからない場合は空文字
            
            # 新しい行を作成（gene_idを最初に追加）
            new_row = [gene_id] + row
            writer.writerow(new_row)
    
    print(f"\nResults:")
    print(f"  Matched genes: {matched_count}")
    print(f"  Unmatched genes: {unmatched_count}")
    
    if unmatched_genes:
        print(f"\nFirst 10 unmatched genes: {unmatched_genes[:10]}")
    
    print(f"\nOutput written to {output_file}")

if __name__ == '__main__':
    input_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts.csv'
    mapping_file = '/Users/godakyosuke/Desktop/idとgene.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    
    add_gene_ids(input_file, mapping_file, output_file)




