#!/usr/bin/env python3
"""
DC_PURKINJE_counts.csvに遺伝子IDを追加し、カウント列名をctl1, ctl2, ... に変更するスクリプト
"""

import csv
import re

def load_gene_mapping(mapping_file):
    """idとgene.csvから遺伝子名→ENSMUSG番号のマッピングを読み込む"""
    gene_to_id = {}
    
    print(f"Loading gene mapping from {mapping_file}...")
    with open(mapping_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_name = row['symbol']
            gene_id = row['User_ID']
            gene_to_id[gene_name] = gene_id
            gene_to_id[gene_name.upper()] = gene_id
            gene_to_id[gene_name.lower()] = gene_id
    
    print(f"Loaded {len(gene_to_id)} gene mappings")
    return gene_to_id

def extract_gene_id_number(gene_id):
    """ENSMUSG番号から数値部分を抽出してソート用のキーを返す"""
    if not gene_id:
        return float('inf')
    match = re.search(r'ENSMUSG(\d+)', gene_id)
    if match:
        return int(match.group(1))
    return float('inf')

def process_file(input_file, mapping_file, output_file):
    """遺伝子IDを追加し、列名を変更し、ソートする"""
    
    # 遺伝子マッピングを読み込む
    gene_to_id = load_gene_mapping(mapping_file)
    
    print(f"\nProcessing {input_file}...")
    
    # データを読み込む
    rows = []
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        for row in reader:
            if not row:
                continue
            
            gene_name = row[0].strip('"')
            gene_id = gene_to_id.get(gene_name, '')
            
            # gene_idを最初に追加
            new_row = [gene_id] + row
            rows.append(new_row)
    
    print(f"Loaded {len(rows)} rows")
    
    # 遺伝子IDでソート
    print("Sorting by gene ID...")
    rows_sorted = sorted(rows, key=lambda x: extract_gene_id_number(x[0] if x else ''))
    
    # 新しいヘッダーを作成
    # 最初の2列: gene_id, gene
    # 残りの列: ctl1, ctl2, ...
    count_columns = header[1:]  # 元のgene列以降
    new_header = ['gene_id', 'gene']
    for i in range(len(count_columns)):
        new_header.append(f'ctl{i+1}')
    
    print(f"Renaming {len(count_columns)} count columns to ctl1, ctl2, ...")
    
    # ファイルに書き込む
    print(f"Writing to {output_file}...")
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(new_header)
        writer.writerows(rows_sorted)
    
    print(f"Done! Created {output_file} with {len(rows_sorted)} rows.")

if __name__ == '__main__':
    input_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts.csv'
    mapping_file = '/Users/godakyosuke/Desktop/idとgene.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    
    process_file(input_file, mapping_file, output_file)




