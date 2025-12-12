#!/usr/bin/env python3
"""
DT_PURKINJE_counts_with_ids.csvのカウントデータ列名をdta1, dta2, ... に変更するスクリプト
"""

import csv

def rename_count_columns(input_file, output_file):
    """カウントデータ列名をdta1, dta2, ... に変更"""
    
    print(f"Reading {input_file}...")
    
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # ヘッダー行を読み込む
        header = next(reader)
        
        # 新しいヘッダーを作成
        # 最初の2列（gene_id, gene）はそのまま、残りをdta1, dta2, ... に変更
        new_header = header[:2]  # gene_id, gene
        
        # カウントデータ列をdta1, dta2, ... に変更
        # 空の列は除外
        count_columns = [col for col in header[2:] if col.strip()]  # 空でない列のみ
        
        for i, col in enumerate(count_columns, start=1):
            new_header.append(f'dta{i}')
        
        print(f"Renaming {len(count_columns)} count columns to dta1, dta2, ...")
        
        # 新しいヘッダーを書き込み
        writer.writerow(new_header)
        
        # データ行を処理
        row_count = 0
        for row in reader:
            if not row:
                continue
            
            # 最初の2列（gene_id, gene）と、カウントデータ列を書き込み
            # 空でない列の数だけカウントデータを取得
            new_row = row[:2]  # gene_id, gene
            if len(row) > 2:
                # カウントデータ列の数だけ追加（空列は除外）
                new_row.extend(row[2:2+len(count_columns)])
            
            writer.writerow(new_row)
            row_count += 1
        
        print(f"Processed {row_count} data rows")
    
    print(f"Done! Output written to {output_file}")

if __name__ == '__main__':
    input_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    
    rename_count_columns(input_file, output_file)

