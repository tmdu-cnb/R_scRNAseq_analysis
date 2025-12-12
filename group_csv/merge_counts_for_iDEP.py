#!/usr/bin/env python3
"""
DC_PURKINJE_counts_with_ids.csvとDT_PURKINJE_counts_with_ids.csvを連結して
counts_for_iDEP.csv形式のファイルを作成するスクリプト
"""

import csv

def merge_counts_files(dc_file, dt_file, output_file):
    """2つのカウントファイルを連結"""
    
    print(f"Reading {dc_file}...")
    dc_data = {}
    dc_header = None
    
    with open(dc_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        dc_header = next(reader)  # gene_id, gene, ctl1, ctl2, ...
        
        for row in reader:
            if not row or not row[0]:  # gene_idが空の行をスキップ
                continue
            gene_id = row[0]
            # gene_idをキーとして、カウントデータ（ctl1-ctl7）を保存
            dc_data[gene_id] = row[2:]  # ctl1-ctl7のデータ
    
    print(f"Loaded {len(dc_data)} genes from DC file")
    
    print(f"\nReading {dt_file}...")
    dt_data = {}
    dt_header = None
    
    with open(dt_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        dt_header = next(reader)  # gene_id, gene, dta1, dta2, ...
        
        for row in reader:
            if not row or not row[0]:  # gene_idが空の行をスキップ
                continue
            gene_id = row[0]
            # gene_idをキーとして、カウントデータ（dta1-dta46）を保存
            dt_data[gene_id] = row[2:]  # dta1-dta46のデータ
    
    print(f"Loaded {len(dt_data)} genes from DT file")
    
    # 両方のファイルに存在する遺伝子IDを取得（ソート済み）
    all_gene_ids = sorted(set(dc_data.keys()) | set(dt_data.keys()))
    print(f"\nTotal unique genes: {len(all_gene_ids)}")
    
    # 新しいヘッダーを作成
    # Gene列 + ctl1-ctl7 + dta1-dta46
    ctl_columns = [col for col in dc_header[2:] if col.strip()]  # ctl1-ctl7
    dta_columns = [col for col in dt_header[2:] if col.strip()]  # dta1-dta46
    new_header = ['Gene'] + ctl_columns + dta_columns
    
    print(f"Writing merged data to {output_file}...")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)  # クォートで囲む
        
        # ヘッダーを書き込み
        writer.writerow(new_header)
        
        # 各遺伝子についてデータを書き込み
        for gene_id in all_gene_ids:
            dc_counts = dc_data.get(gene_id, ['0'] * len(ctl_columns))  # 存在しない場合は0で埋める
            dt_counts = dt_data.get(gene_id, ['0'] * len(dta_columns))  # 存在しない場合は0で埋める
            
            # データが不足している場合は0で埋める
            while len(dc_counts) < len(ctl_columns):
                dc_counts.append('0')
            while len(dt_counts) < len(dta_columns):
                dt_counts.append('0')
            
            # Gene列 + DCカウント + DTカウント
            row = [gene_id] + dc_counts[:len(ctl_columns)] + dt_counts[:len(dta_columns)]
            writer.writerow(row)
    
    print(f"Done! Created {output_file} with {len(all_gene_ids)} genes.")

if __name__ == '__main__':
    dc_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    dt_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/counts_for_iDEP.csv'
    
    merge_counts_files(dc_file, dt_file, output_file)




