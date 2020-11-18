Title:阪医Python会特製 RNA-seq pipeline ver. 1.0 リリース
Date: 2019.03.19
Tags: Bioinformatics
Authors: 菅波, 西田, 大森, 安水, 小川, 山田, 川嶋, 川崎, 平岡, 廣瀬, 柳澤, 淡田, 中村, 依藤

阪医Python会のbioinformaticsチームの一つの成果として、RNA-seqのパイプラインのv1.0がリリースとなったので記事とさせていただきます。SRR idから遺伝子✕サンプルのテーブルにするまでには意外に大変ですが、それをすべて自動化しました。ダウンロード、詳細等は以下にあります。

**[https://github.com/yyoshiaki/auto_counttable_maker](https://github.com/yyoshiaki/auto_counttable_maker)**

なお、以下のイラストはikraのアイコンとなっています。

<img src="{attach}images/original_rnaseq_pipeline_figs/ikra.png" width="250px"> 

## 特徴

今回、自分たちの使いやすさを考えてツールの設計を行いました。他サンプルのファイル名の管理など、煩わしいところをすべて自動化することで、ミスも減らせると思います。

1. 準備するのは簡単なCSVファイル（テーブルだけ）で、極力パラメーターを減らした。
2. すべてdocker上で動くため、ツールを各々インストールする必要がないし、バージョンに苦しむこともない。
3. udockerにも対応しているため、ユーザー権限しかないサーバー上でも実行可能。
4. outputは[idep](http://bioinformatics.sdstate.edu/idep/)に対応。
5. もちろんマルチスレッド対応。

## 使い方

必要なテーブルは

| name | SRR or fastq | Layout | condition1 | ... |
| ---- | ---- | - | - | - |
| Treg_LN_1 | SRR5385247 | SE | Treg | ... |
| Treg_LN_2 | SRR5385248 | SE | Treg | ... |

のような形式で、前3列が必須です。簡単ですね！データの集め方は、論文についているaccession number等をたどるのでもいいし、新しくなって爆速になった[DDBJ Search](http://sra.dbcls.jp/)もおすすめ。

コマンドはオプションが指定でき、リード数を100000に絞ったテストモードやマルチプロセスにも対応。

```txt
Usage: bash MakeCountTable_Illumina_trimgalore_SRR.sh experiment_table.csv spiece [--test, --help, --without-docker, --udocker] [--threads [VALUE]]
args
1.experiment matrix(csv)
2.reference(human or mouse)

Options:
--test test mode(MAX_SPOT_ID=100000).(dafault : False)
-u, --udocker
-w, --without-docker
-t, --threads
-h, --help Show usage.
```

なお、自前のfastq filesからの実行はv1.1で載せようと思っています。また、出力はscaled TPMを採用。(Soneson, C., Love, M. I. &amp; Robinson, M. D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4, 1521 (2015).)。

## pipelineの構成

1. fasterq-dump : シーケンスデータの取得
2. fastqc : QC
3. trim-galore : トリミング
4. salmon : RNA定量
5. multiqc : QCログの回収、可視化
6. tximport : 遺伝子テーブルの生成

となっています。各ツールの説明は省きますが、今時のツールの選定になっていると思います。

## idep

本ツールは[idep](http://bioinformatics.sdstate.edu/idep/)を意識した設計になっています。idepはRNA-seqの解析をinteractiveに行えるプラットフォームで、Differential expressed genes(DEGs)の検出だけではなく、遺伝子、サンプルのクラスタリング、パスウェイ解析、可視化などが行えます。idepについては以下がとても参考になります。

[macでインフォマティクス : インタラクティブなRNA seq解析webアプリケーション iDEP](http://kazumaxneo.hatenablog.com/entry/2018/12/29/153838)


![2]({attach}images/original_rnaseq_pipeline_figs/screenshot-from-2019-03-19-23-14-31.png)

## githubを用いたチーム開発

今回、githubを用いてチーム開発を行いました。githubはエンジニアの間では当たり前のツールで、チームでのソフト開発によく用いられます。bioinformatics界隈でチーム開発を経験できることは意外に少なく、非常にいい経験になりました。

![3]({attach}images/original_rnaseq_pipeline_figs/ikra_git.png)

雑多にはなりましたが、阪医Python会bioinformaticsチームの成果をアナウンスさせていただきました。完成までには3ヶ月ほどを要し、各人のアイデアや努力が詰まっております。今後もどんどん開発を進めていこうと思います。また、皆様のissue, Pull Requestもお待ちしております。

*Enjoy bioinformatics life!*

![4]({attach}images/original_rnaseq_pipeline_figs/ios-e381aee794bbe5838f.jpg)
