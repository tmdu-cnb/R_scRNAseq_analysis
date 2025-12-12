#!/usr/bin/env python3
"""
CSVファイルのカウントデータ列名をctl1, ctl2, ... に変更するスクリプト
"""

import csv

def rename_count_columns(input_file, output_file):
    """カウントデータ列名をctl1, ctl2, ... に変更"""
    
    print(f"Reading {input_file}...")
    
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # ヘッダー行を読み込む
        header = next(reader)
        
        # 新しいヘッダーを作成
        # 最初の2列（gene_id, gene）はそのまま、残りをctl1, ctl2, ... に変更
        new_header = header[:2]  # gene_id, gene
        
        # カウントデータ列をctl1, ctl2, ... に変更
        count_columns = header[2:]  # カウントデータ列
        for i, col in enumerate(count_columns, start=1):
            new_header.append(f'ctl{i}')
        
        print(f"Renaming {len(count_columns)} count columns to ctl1, ctl2, ...")
        
        # 新しいヘッダーを書き込み
        writer.writerow(new_header)
        
        # データ行をそのまま書き込み
        row_count = 0
        for row in reader:
            writer.writerow(row)
            row_count += 1
        
        print(f"Processed {row_count} data rows")
    
    print(f"Done! Output written to {output_file}")

if __name__ == '__main__':
    input_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    
    rename_count_columns(input_file, output_file)




