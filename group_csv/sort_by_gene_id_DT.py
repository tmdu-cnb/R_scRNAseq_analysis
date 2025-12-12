#!/usr/bin/env python3
"""
DT_PURKINJE_counts_with_ids.csvを遺伝子ID（ENSMUSG番号）の少ない順に並べ替えるスクリプト
"""

import csv
import re

def extract_gene_id_number(gene_id):
    """ENSMUSG番号から数値部分を抽出してソート用のキーを返す"""
    if not gene_id:
        return float('inf')  # 空の場合は最後に
    
    # ENSMUSG00000000759 から数値部分を抽出
    match = re.search(r'ENSMUSG(\d+)', gene_id)
    if match:
        return int(match.group(1))
    return float('inf')  # パターンに一致しない場合は最後に

def sort_csv_by_gene_id(input_file, output_file):
    """CSVファイルを遺伝子IDでソート"""
    
    print(f"Reading {input_file}...")
    
    # ファイルを読み込む
    rows = []
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # ヘッダーを取得
        
        for row in reader:
            if row:  # 空行をスキップ
                rows.append(row)
    
    print(f"Loaded {len(rows)} rows")
    
    # 遺伝子IDでソート（最初の列がgene_id）
    print("Sorting by gene ID...")
    
    # ヘッダー行を除いてソート
    # 空のgene_idを持つ行は最後に
    rows_sorted = sorted(rows, key=lambda x: extract_gene_id_number(x[0] if x else ''))
    
    # ソート済みファイルを書き込む
    print(f"Writing sorted data to {output_file}...")
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)  # ヘッダーを書き込み
        writer.writerows(rows_sorted)
    
    print(f"Done! Sorted {len(rows_sorted)} rows.")

if __name__ == '__main__':
    input_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    
    sort_csv_by_gene_id(input_file, output_file)




