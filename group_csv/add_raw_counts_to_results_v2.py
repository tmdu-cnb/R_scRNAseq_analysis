#!/usr/bin/env python3
"""
Results_LFC_Pval_DESeq2.csvにDCとDTの生カウントデータを追加するスクリプト
既存のカウント列（正規化済み）を生カウント（整数値）に置き換える
"""

import csv

def load_count_data(count_file):
    """カウントファイルから遺伝子IDをキーとしてカウントデータを読み込む"""
    count_data = {}
    
    print(f"Loading count data from {count_file}...")
    with open(count_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # カウントデータ列のインデックスを取得（gene_idとgene列を除く）
        count_cols = [col.strip('"') for col in header[2:] if col.strip()]  # gene_id, geneの後の列
        
        for row in reader:
            if not row or not row[0]:
                continue
            
            gene_id = row[0].strip()
            if gene_id:
                # カウントデータを取得
                counts = row[2:2+len(count_cols)]  # gene_id, geneの後のデータ
                count_data[gene_id] = counts
    
    print(f"Loaded {len(count_data)} genes with {len(count_cols)} count columns")
    return count_data, count_cols

def add_counts_to_results(results_file, dc_file, dt_file, output_file):
    """Resultsファイルにカウントデータを追加"""
    
    # DCとDTのカウントデータを読み込む
    dc_data, dc_cols = load_count_data(dc_file)
    dt_data, dt_cols = load_count_data(dt_file)
    
    print(f"\nProcessing {results_file}...")
    
    # Resultsファイルを読み込んで、カウントデータを追加
    rows_data = []
    
    try:
        with open(results_file, 'r', encoding='utf-8') as infile:
            reader = csv.reader(infile)
            header = next(reader)
            rows_data = [row for row in reader if row]
    except StopIteration:
        print("Warning: Results file appears to be empty or has no data rows")
        return
    
    if not rows_data:
        print("Error: No data rows found in Results file")
        return
    
    print(f"Found {len(rows_data)} data rows in Results file")
    
    # ヘッダーの構造を確認
    # User_IDまたはensembl_ID列を探す
    user_id_idx = None
    for i, col in enumerate(header):
        col_clean = col.strip('"')
        if col_clean == 'User_ID':
            user_id_idx = i
            break
        elif col_clean == 'ensembl_ID' and user_id_idx is None:
            user_id_idx = i
    
    if user_id_idx is None:
        print("Warning: Could not find User_ID or ensembl_ID column")
        # 2番目の列を試す
        if len(header) > 1:
            user_id_idx = 1
    
    print(f"Using column index {user_id_idx} for gene ID")
    
    # 既存のカウント列があるかチェック
    existing_ctl_start = None
    existing_dta_start = None
    
    for i, col in enumerate(header):
        col_clean = col.strip('"')
        if col_clean == 'ctl1' and existing_ctl_start is None:
            existing_ctl_start = i
        if col_clean == 'dta1' and existing_dta_start is None:
            existing_dta_start = i
    
    # 新しいヘッダーを作成
    if existing_ctl_start is not None and existing_dta_start is not None:
        # 既存のカウント列を置き換える
        new_header = header[:existing_ctl_start] + dc_cols + dt_cols
        if existing_dta_start + len(dt_cols) < len(header):
            new_header.extend(header[existing_dta_start + len(dt_cols):])
        replace_mode = True
    else:
        # カウント列がない場合は追加
        new_header = header + dc_cols + dt_cols
        replace_mode = False
    
    # ファイルに書き込む
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        
        writer.writerow(new_header)
        
        matched_count = 0
        
        for row in rows_data:
            # User_IDを取得
            user_id = None
            if user_id_idx is not None and user_id_idx < len(row):
                user_id = row[user_id_idx].strip('"')
            
            # 新しい行を作成
            if replace_mode:
                # 既存のカウント列を置き換える
                new_row = row[:existing_ctl_start]
                
                # DCカウントデータを追加
                dc_counts = dc_data.get(user_id, ['0'] * len(dc_cols)) if user_id else ['0'] * len(dc_cols)
                new_row.extend(dc_counts[:len(dc_cols)])
                
                # DTカウントデータを追加
                dt_counts = dt_data.get(user_id, ['0'] * len(dt_cols)) if user_id else ['0'] * len(dt_cols)
                new_row.extend(dt_counts[:len(dt_cols)])
                
                # 残りの列を追加
                if existing_dta_start + len(dt_cols) < len(row):
                    new_row.extend(row[existing_dta_start + len(dt_cols):])
            else:
                # カウント列がない場合は追加
                new_row = row.copy()
                
                # DCカウントデータを追加
                dc_counts = dc_data.get(user_id, ['0'] * len(dc_cols)) if user_id else ['0'] * len(dc_cols)
                new_row.extend(dc_counts[:len(dc_cols)])
                
                # DTカウントデータを追加
                dt_counts = dt_data.get(user_id, ['0'] * len(dt_cols)) if user_id else ['0'] * len(dt_cols)
                new_row.extend(dt_counts[:len(dt_cols)])
            
            writer.writerow(new_row)
            
            if user_id and (user_id in dc_data or user_id in dt_data):
                matched_count += 1
        
        print(f"Matched {matched_count} genes with count data")
    
    print(f"Done! Output written to {output_file}")

if __name__ == '__main__':
    results_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    dc_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    dt_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2_with_raw_counts.csv'
    
    add_counts_to_results(results_file, dc_file, dt_file, output_file)




