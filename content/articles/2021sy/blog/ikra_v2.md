Title:RNA-seq pipeline ikra v2.0 リリース
Date: 2021.05.05
Tags: Bioinformatics
Authors: 山崎, 北畠, 梅津, 安水

阪大医学部Python会では、RNA-Seqのパイプラインである **ikra** ([https://github.com/yyoshiaki/ikra](https://github.com/yyoshiaki/ikra))をリリースしておりました (cf. [阪医Python会特製 RNA-seq pipeline ver. 1.0 リリース](https://oumpy.github.io/blog/2019/03/original_rnaseq_pipeline.html))。このたびikraの機能を大幅にアップデートし、ikra v2.0としてリリースしました。アップデートのポイントは以下の通りです。

## aligner optionの追加
- `-a, --align (hisat2 or star)`モードを付けて実行すると、salmonによる定量化と同時にリファレンスゲノムへのマッピングを行う機能を実装しました。alignerはhisat2とSTARのいずれかを選べます。
- マッピングではbamファイルとbwファイルが生成されますので、これらをIGVブラウザで直接閲覧することができます。IGVとはマッピング結果を可視化できるツールです。下の画像は、テストデータのマッピング結果のうち、Actb（βアクチン）遺伝子にマップされたtranscriptを可視化したものです。

![IGV]({attach}./images/ikra_v2_figs/IGV_Actb.png)

- リファレンスゲノムは、humanならGRCh38、mouseならGRCm38です。

## gencode optionの追加
- `-g, --gencode (version)`モードを付けて実行すると、salmonに利用するreference transcriptomeのバージョンを指定できるようにしました。
- デフォルトは、humanならv37、mouseならvM26です。

## ツールのアップデート
- ikraの実行に用いるツール群をアップデートし、新しいものに入れ替えました。


今後とも、阪医Python会およびikraをどうぞよろしくお願いいたします。