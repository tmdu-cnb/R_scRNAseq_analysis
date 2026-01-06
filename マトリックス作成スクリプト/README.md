# スクリプト

## 内容

- **add_gene_id.py** - DC用の遺伝子ID追加スクリプト
- **add_gene_id_DT.py** - DT用の遺伝子ID追加スクリプト
- **add_ids_and_rename_columns.py** - DC用のID追加＋カラム名変更スクリプト
- **add_ids_and_rename_columns_DT.py** - DT用のID追加＋カラム名変更スクリプト
- **ensmusg idとgenenameの対応表.csv** - 遺伝子名とIDのマッピングテーブル

## 説明

Purkinje細胞限定マトリックスに遺伝子IDを追加するためのPythonスクリプトです。

## 各スクリプトの詳細

### add_gene_id.py / add_gene_id_DT.py
- **機能**: Purkinje細胞限定マトリックス（DC_PURKINJE_counts.csv または DT_PURKINJE_counts.csv）に遺伝子ID（ENSMUSG番号）を追加
- **入力**: 
  - Purkinje細胞限定マトリックス（遺伝子名のみ）
  - 遺伝子名とIDのマッピングテーブル
- **出力**: `*_with_ids.csv` ファイル（最初の列にgene_idが追加される）
- **処理内容**: 
  1. マッピングテーブルから遺伝子名→ENSMUSG番号の対応を読み込み
  2. 各遺伝子行の先頭にgene_id列を追加
  3. マッチしない遺伝子は空欄のまま出力

### add_ids_and_rename_columns.py / add_ids_and_rename_columns_DT.py
- **機能**: 遺伝子ID追加に加えて、サンプル列名を統一形式（ctl1, ctl2, ... または dta1, dta2, ...）に変更し、遺伝子ID順にソート
- **入力**: 
  - Purkinje細胞限定マトリックス（遺伝子名のみ）
  - 遺伝子名とIDのマッピングテーブル
- **出力**: `*_with_ids.csv` ファイル（ID追加＋列名変更＋ソート済み）
- **処理内容**: 
  1. 遺伝子IDを追加
  2. サンプル列名を正規化（例: ATGGTGAGTGGGTGTGATACGTCA-9 → ctl1）
  3. 遺伝子IDの数値部分でソート（ENSMUSG00000000001 → ENSMUSG00000000002 → ...）
  4. iDEPなどの解析ツールで使用しやすい形式に整形

### ensmusg idとgenenameの対応表.csv
- **内容**: 遺伝子シンボル名（symbol）とENSMUSG番号（User_ID）の対応表
- **用途**: 上記スクリプトで遺伝子名をIDに変換する際の参照テーブル
- **形式**: CSV形式、列は `symbol` と `User_ID`

