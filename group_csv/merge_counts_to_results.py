#!/usr/bin/env python3
"""
Results_LFC_Pval_DESeq2.csvの右側にcounts_for_iDEP.csvを結合するスクリプト
"""

import csv

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
    
    # counts_for_iDEP.csvのデータを読み込む
    counts_data, count_cols = load_counts_data(counts_file)
    
    print(f"\nProcessing {results_file}...")
    
    # Resultsファイルを読み込んで結合
    with open(results_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        
        # ヘッダー行を読み込む
        header = next(reader)
        
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
        # 重複を避けるため、既存の列名をチェック
        new_header = header.copy()
        
        # 既存の列名を確認
        existing_cols = {col.strip('"') for col in header}
        
        # counts_for_iDEPの列を追加（重複していない場合）
        for col in count_cols:
            if col not in existing_cols:
                new_header.append(col)
            else:
                # 重複している場合は、_iDEPを追加
                new_header.append(f"{col}_iDEP")
        
        writer.writerow(new_header)
        
        # データ行を処理
        row_count = 0
        matched_count = 0
        
        for row in reader:
            if not row:
                continue
            
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
            row_count += 1
        
        print(f"Processed {row_count} rows")
        print(f"Matched {matched_count} genes with count data")
    
    print(f"\nDone! Output written to {output_file}")

if __name__ == '__main__':
    results_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    counts_file = '/Users/godakyosuke/Desktop/group_csv/counts_for_iDEP.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    
    merge_files(results_file, counts_file, output_file)
