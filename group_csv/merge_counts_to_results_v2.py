#!/usr/bin/env python3
"""
Results_LFC_Pval_DESeq2.csvの右側にcounts_for_iDEP.csvを結合するスクリプト
"""

import csv
import os

def load_counts_data(counts_file):
    """counts_for_iDEP.csvから遺伝子IDをキーとしてカウントデータを読み込む"""
    counts_data = {}
    
    print(f"Loading count data from {counts_file}...")
    with open(counts_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # Gene列を除いたカウントデータ列を取得
        count_cols = [col.strip('"') for col in header[1:] if col.strip()]  # Gene列の後の列
        
        for row in reader:
            if not row or not row[0]:
                continue
            
            gene_id = row[0].strip('"')
            if gene_id:
                # カウントデータを取得（Gene列を除く）
                counts = row[1:1+len(count_cols)]
                counts_data[gene_id] = counts
    
    print(f"Loaded {len(counts_data)} genes with {len(count_cols)} count columns")
    return counts_data, count_cols

def merge_files(results_file, counts_file, output_file):
    """Resultsファイルの右側にcounts_for_iDEPのデータを結合"""
    
    # ファイルの存在確認
    if not os.path.exists(results_file) or os.path.getsize(results_file) == 0:
        print(f"Error: {results_file} is empty or does not exist")
        return
    
    # counts_for_iDEP.csvのデータを読み込む
    counts_data, count_cols = load_counts_data(counts_file)
    
    print(f"\nProcessing {results_file}...")
    
    # Resultsファイルを読み込んで結合
    rows_data = []
    header = None
    
    try:
        with open(results_file, 'r', encoding='utf-8') as infile:
            reader = csv.reader(infile)
            header = next(reader)
            rows_data = [row for row in reader if row]
    except (StopIteration, FileNotFoundError) as e:
        print(f"Error reading {results_file}: {e}")
        return
    
    if not rows_data:
        print(f"Error: No data rows found in {results_file}")
        return
    
    print(f"Found {len(rows_data)} data rows in Results file")
    
    # User_ID列のインデックスを探す
    user_id_idx = None
    for i, col in enumerate(header):
        col_clean = col.strip('"')
        if col_clean == 'User_ID':
            user_id_idx = i
            break
    
    if user_id_idx is None:
        # ensembl_ID列を試す
        for i, col in enumerate(header):
            col_clean = col.strip('"')
            if col_clean == 'ensembl_ID':
                user_id_idx = i
                break
    
    if user_id_idx is None:
        print("Error: Could not find User_ID or ensembl_ID column")
        return
    
    print(f"Using column index {user_id_idx} for gene ID matching")
    
    # 新しいヘッダーを作成（counts_for_iDEPの列を追加）
    # 既存の列名を確認
    existing_cols = {col.strip('"') for col in header}
    
    # counts_for_iDEPの列を追加（重複していない場合）
    new_header = header.copy()
    for col in count_cols:
        if col not in existing_cols:
            new_header.append(col)
        else:
            # 重複している場合は、_iDEPを追加
            new_header.append(f"{col}_iDEP")
    
    # ファイルに書き込む
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        
        writer.writerow(new_header)
        
        matched_count = 0
        
        for row in rows_data:
            # User_IDを取得
            user_id = None
            if user_id_idx < len(row):
                user_id = row[user_id_idx].strip('"')
            
            # 新しい行を作成
            new_row = row.copy()
            
            # counts_for_iDEPのデータを追加
            if user_id and user_id in counts_data:
                counts = counts_data[user_id]
                # データが不足している場合は0で埋める
                while len(counts) < len(count_cols):
                    counts.append('0')
                new_row.extend(counts[:len(count_cols)])
                matched_count += 1
            else:
                # マッチしない場合は0で埋める
                new_row.extend(['0'] * len(count_cols))
            
            writer.writerow(new_row)
        
        print(f"Matched {matched_count} genes with count data")
    
    print(f"Done! Output written to {output_file}")

if __name__ == '__main__':
    results_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    counts_file = '/Users/godakyosuke/Desktop/group_csv/counts_for_iDEP.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    
    merge_files(results_file, counts_file, output_file)




