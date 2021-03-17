Title: TCGAbiolinksとその使い方
Date: 2018.11.11
Modified: 2018.11.11
Tags: Bioinformatics
Author: 淡田

皆さん初めまして淡田公久と申します。
リレー投稿初週ということで、僕が現在基礎配属で使用していますデータやツールに関してのお話をさせて頂きます。

## TCGA
TCGAというサイトをご存じでしょうか？
TCGAは２００６年からアメリカが開始したガンゲノムプロジェクトで、２０種類以上の癌腫についての「ゲノム、メチル化異常、遺伝子発現」などの解析データをまとめている大型ゲノムバンクのような役割をしております。

[Home - The Cancer Genome Atlas - Cancer Genome - TCGA](https://cancergenome.nih.gov/)

実際に検索してみると分かりますが、データ数はほんとに膨大です（日本ではここまで大きのがなく、悲しいところ(´;ω;｀)）

中には「contorlled」となっており、アクセスできないデータもありますが、多くは「open」となってアクセスできるものがほとんどなので解析などの際に困ることは無いはず。

## TCGAbiolinks
では実際にそのデータを使ってどう解析するかといったときにどうするか？？

ここで本題の**TCGAbiolinks**のお話になります。
TCGAbiolinksとは、Rというツール上でTCGAのデータを解析するためのパッケージです。
> [Bioconductor - TCGAbiolinks](http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)

こちらにインストールの方法やらが詳しく説明されています。

先ほども述べた通り、TCGAにはメチル化データや、発現量データ（他にも臨床データなど）があるため、それぞれ解析方法は異なりますが、大枠のworkflowは

```
(1) queryという関数でデータを引用
    ↓
(2) downloadで引用データをダウンロード
    ↓
(3) データーをprepare
    ↓
(4) 目的に応じて解析
    ↓
(5) 結果を画像化
```
という流れです。
まぁ僕自身もまだこのツールにかんして半人前にもなっていないレベルなのですが（笑）、openソースの癌データをRというフリーツールで解析できるのはおもしろいのではないでしょうか？

TCGAのデータ解析には他にもsubio platform（[マイクロアレイ・ＮＧＳの無料解析ソフト \| Subio Platform](https://www.subioplatform.com/ja/products/subioplatform/)）というものあるようですが、TCGAbiolinksのがメジャーではあるでしょう。

このbiolinksに付随して**elmer**というパッケージもあるのでそれについてもお話しできればと思います。

ありがとうございました。
