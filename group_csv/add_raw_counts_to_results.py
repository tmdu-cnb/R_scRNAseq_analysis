#!/usr/bin/env python3
"""
Results_LFC_Pval_DESeq2.csvにDCとDTの生カウントデータを追加するスクリプト
"""

import csv

def load_count_data(count_file, gene_id_col=0):
    """カウントファイルから遺伝子IDをキーとしてカウントデータを読み込む"""
    count_data = {}
    
    print(f"Loading count data from {count_file}...")
    with open(count_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # カウントデータ列のインデックスを取得（gene_idとgene列を除く）
        count_cols = header[2:]  # gene_id, geneの後の列
        
        for row in reader:
            if not row or not row[gene_id_col]:
                continue
            
            gene_id = row[gene_id_col].strip()
            if gene_id:
                # カウントデータを取得
                counts = row[2:2+len(count_cols)]  # gene_id, geneの後のデータ
                count_data[gene_id] = counts
    
    print(f"Loaded {len(count_data)} genes")
    return count_data, count_cols

def add_counts_to_results(results_file, dc_file, dt_file, output_file):
    """Resultsファイルにカウントデータを追加"""
    
    # DCとDTのカウントデータを読み込む
    dc_data, dc_cols = load_count_data(dc_file)
    dt_data, dt_cols = load_count_data(dt_file)
    
    print(f"\nProcessing {results_file}...")
    
    # Resultsファイルを読み込んで、カウントデータを追加
    with open(results_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        
        # ヘッダー行を読み込む
        header = next(reader)
        
        # 既存のカウント列があるかチェック
        # 既存のctl1-ctl7, dta1-dta46列を探す
        existing_ctl_start = None
        existing_dta_start = None
        
        for i, col in enumerate(header):
            if col == 'ctl1' or col == '"ctl1"':
                existing_ctl_start = i
            if col == 'dta1' or col == '"dta1"':
                existing_dta_start = i
        
        # 新しいヘッダーを作成
        if existing_ctl_start is not None and existing_dta_start is not None:
            # 既存のカウント列を置き換える
            new_header = header[:existing_ctl_start] + dc_cols + dt_cols
            if existing_dta_start + len(dt_cols) < len(header):
                new_header.extend(header[existing_dta_start + len(dt_cols):])
        else:
            # カウント列がない場合は追加
            new_header = header + dc_cols + dt_cols
        
        writer.writerow(new_header)
        
        # データ行を処理
        row_count = 0
        matched_count = 0
        
        for row in reader:
            if not row:
                continue
            
            # User_IDまたはensembl_IDを取得
            user_id = None
            if len(row) > 2:
                # User_ID列を探す
                for i, col in enumerate(header):
                    if col == 'User_ID' or col == '"User_ID"':
                        if i < len(row):
                            user_id = row[i].strip('"')
                            break
                    elif col == 'ensembl_ID' or col == '"ensembl_ID"':
                        if i < len(row):
                            user_id = row[i].strip('"')
                            break
            
            # 新しい行を作成
            if existing_ctl_start is not None and existing_dta_start is not None:
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
            row_count += 1
            
            if user_id and (user_id in dc_data or user_id in dt_data):
                matched_count += 1
        
        print(f"Processed {row_count} rows")
        print(f"Matched {matched_count} genes with count data")
    
    print(f"\nDone! Output written to {output_file}")

if __name__ == '__main__':
    results_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    dc_file = '/Users/godakyosuke/Desktop/group_csv/DC_PURKINJE_counts_with_ids.csv'
    dt_file = '/Users/godakyosuke/Desktop/group_csv/DT_PURKINJE_counts_with_ids.csv'
    output_file = '/Users/godakyosuke/Desktop/group_csv/Results_LFC_Pval_DESeq2.csv'
    
    add_counts_to_results(results_file, dc_file, dt_file, output_file)




