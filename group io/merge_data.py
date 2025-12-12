#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv
import sys

# ファイルパス
csv_file = '/Users/godakyosuke/Desktop/group io/C_ioとT_ioのDEG.csv'
txt_file = '/Users/godakyosuke/Desktop/group io/counts_for_iDEP　1118 remove p.txt'
output_file = '/Users/godakyosuke/Desktop/group io/merged_data.csv'

# txtファイルのデータを読み込む（Gene列をキーとして）
txt_data = {}
with open(txt_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    txt_columns = header[1:]  # Gene列以外のカラム名
    
    for row in reader:
        if row:  # 空行をスキップ
            gene_id = row[0]
            txt_data[gene_id] = row[1:]

# CSVファイルを読み込んで結合
merged_rows = []
with open(csv_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    csv_header = next(reader)
    
    # 出力ヘッダーを作成（CSVの全カラム + txtのカラム）
    output_header = csv_header + txt_columns
    merged_rows.append(output_header)
    
    # ensembl_ID列のインデックスを取得
    ensembl_id_idx = csv_header.index('ensembl_ID')
    
    # 各行を処理
    for row in reader:
        if row:
            ensembl_id = row[ensembl_id_idx].strip('"')  # クォートを除去
            
            # txtデータから対応するデータを取得
            if ensembl_id in txt_data:
                merged_row = row + txt_data[ensembl_id]
            else:
                # マッチしない場合は空値を追加
                merged_row = row + [''] * len(txt_columns)
            
            merged_rows.append(merged_row)

# 結果を出力
with open(output_file, 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(merged_rows)

print(f"結合が完了しました。出力ファイル: {output_file}")
print(f"結合された行数: {len(merged_rows) - 1}")  # ヘッダーを除く

