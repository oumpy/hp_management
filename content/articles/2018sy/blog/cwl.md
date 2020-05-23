Title:Common Workflow Language入門
Date: 2018.12.07
Tags: Bioinformatics
Author: 安水

先日[workflow-meetup](https://github.com/manabuishii/workflow-meetup/wiki/20181126)にお誘いを頂いて参加してきました。そこでThe Common Workflow Language(CWL)というものを習ったので忘れないうちに**医学部生にもわかりやすく**まとめます。

## CWLって？

![cwl](https://github.com/common-workflow-language/cwl-website/blob/master/site/CWL-Logo-Header.png?raw=true)

一言で言うと、bioinformatics処理の自動化です。どういうことかもうちょっと詳しく見てみましょう。

&gt;ソフトウェアを組み合わせて構成される一連の作業手順 (ワークフロー) を記述するための仕組みや言語、GUI ソフトウェアは既に多く存在します。しかし、それらは特定の実行環境 (ハードウェア、ソフトウェア) に依存したものであり、異なる環境の間でワークフローを共有、再実行することは困難です。この問題を解決するために、異なるワークフロー実行ソフトウェア (実行エンジン) の間で共通してインポート/エクスポートできるフォーマットを目指して、Common Workflow Language (CWL) の開発が始まりました。

[素晴らしい日本語ドキュメント](https://github.com/pitagora-galaxy/cwl/wiki/CWL-Start-Guide-JP)から取ってきました。要はbioinformaticsなど、複数のソフトの連携が必要な処理を自動化、共有がスムーズになるような仕組みです。ローカルでテストして、サーバーやクラウドに持っていくのも簡単です。最近はSangerやBroad Institute(ハーバードとMITの共同研究施設)などでも使われているとのことです。

まずは**[60秒でわかるCWL - youtube](https://www.youtube.com/embed/86eY8xs-Vo8?cc_load_policy=1&amp;cc_lang_pref=ja&amp;autoplay=1)**をみてください。

`.cwl`で記述できるのはtoolとworkflowの二種類です。toolは例えばSTARやkallistoなど一つのツールを一つのcwlファイルに記述します。workflowはそれらのtoolをbindして実際のworkflowにします。更に、workflowをnestして新しいworkflowを作ったりも出来ます。

## 始め方

1. [素晴らしい日本語ドキュメント](https://github.com/pitagora-galaxy/cwl/wiki/CWL-Start-Guide-JP)を読みましょう。10分もあれば読めると思います。
2. [cwl-intro-gui-workshop スライド](https://docs.google.com/presentation/d/13K6BKQoimiaOIFSnComw8in3HlCxEm820OASP8_jnDI/edit?usp=sharing)をやってみましょう。1時間ほどです。
3. [user guide](http://www.commonwl.org/user_guide/)をやりましょう。3時間ほどで終わるように設計されています。(理論値)

## CWLのここがすごい

入門は各サイトにお任せするとして、今回は私がCWLをさわってみてすごいと思ったところをまとめてみます。素人目線です。

### rabix composerのGUIがきれいすぎる

![rabix gif](https://github.com/rabix/composer/raw/master/doc/images/workflows.gif)

画像はbioinformaticsでお馴染みのSTARをrabix composerで操っているところです。[rabix composer](http://rabix.io/)を使えば、GUIで直感的にCWLが触れます。更に出来たworkflowをクラウドやサーバーにそのまま持っていってCLI越しに使うことが出来ます。

### dockerをベースにして完全にreproducibleにできる

いかにも今どきなのですが、toolの記述はすべてdockerのコンテナを使って記述することが出来ます。つまり、環境が変わってインストールし直す手間がまったく無いということです。すばらしいですね。さらに`--user-space-docker-cmd=udocker`というオプション一つで、udockerというUser権限しか無い環境でも使えるdockerを使って動かすことができるようです。

### すでにいろいろなソフトのcwlが出来ている

Communityで管理されている[CWL Tools &amp; Workflows](https://github.com/common-workflow-language/workflows)があります。ほかにもDBCLS太田さんの管理されている[Pitagora Workflows in CWL](https://github.com/pitagora-galaxy/cwl)もいろいろ揃っています。

### cwltoolの安定感がすごい

出来たcwlファイルをCLI環境で使うにはcwltoolを使います。感動したのはヘルプが自動で生成されるところ。しかも、pythonで書かれていてpythonに組み込むことも出来ます。

```python
import cwltool.factory
fac = cwltool.factory.Factory()

echo = fac.make("echo.cwl")
result = echo(inp="foo")

# result["out"] == "foo"
```

### 日本人contributerが多い

困ったときに助けてくれるやさしいcontributerの方々が日本に何人かいます。また、日本語ドキュメントも多くてとても助かります。更に、開発者の[Michael R. Crusoe](https://twitter.com/biocrusoe)は来日しており、[#CommonWLjp](https://twitter.com/hashtag/CommonWLjp?src=hash)などでCWLの輪を広げています。

<blockquote class="twitter-tweet" data-lang="ja"><p lang="en" dir="ltr">On my way to the 11th <a href="https://twitter.com/NBDC_info?ref_src=twsrc%5Etfw">@NBDC_info</a> / <a href="https://twitter.com/dbcls?ref_src=twsrc%5Etfw">@dbcls</a> <a href="https://twitter.com/hashtag/biohack18?src=hash&amp;ref_src=twsrc%5Etfw">#biohack18</a> (with a week long stop in Tokyo for <a href="https://twitter.com/hashtag/usegalaxyjp?src=hash&amp;ref_src=twsrc%5Etfw">#usegalaxyjp</a> meeting and other <a href="https://twitter.com/hashtag/CommonWL?src=hash&amp;ref_src=twsrc%5Etfw">#CommonWL</a> presentations and a bilingual workshop) <a href="https://t.co/W4gsO9bnA1">pic.twitter.com/W4gsO9bnA1</a></p>&mdash; Michael R. Crusoe (@biocrusoe) <a href="https://twitter.com/biocrusoe/status/1069104245100748801?ref_src=twsrc%5Etfw">2018年12月2日</a></blockquote>


## kallistoをCWL + Rabix Composerで試してみる

まずは実際の作業工程を見てみましょう。

![gif](https://github.com/yyoshiaki/cwl_user_guide/blob/master/kallisto/kallisto.gif?raw=true)

ちゃんとkallistoのoutputが生成されていますね！データは[github](https://github.com/yyoshiaki/cwl_user_guide/tree/master/kallisto)にまとめておきました。テスト用のシーケンスデータも付けてあるので、クローンして遊んでみてください。`kallisto-index.cwl`と`kallisto-quant.cwl`は公式の[CWL Tools &amp; Workflows](https://github.com/common-workflow-language/workflows)から取ってきました。もちろんこれらはdocker imageを使っています。新しいworkflowを作ってみましょう。動画のように直感的にできると思います。

今回、errorが出て進まなかった場所があったのですが、石井さん、西田さんの応援や、cwlの作者のMichaelやRabixの作者のKaushikがerrorに対処してくれました。
> [issue : kallisto workflow is incompatible with Rabix composer #418](https://github.com/rabix/composer/issues/418#issuecomment-444257454)

この対応のフレンドリーさはbioinformaticsならではですね。

## workflow業界のあれこれ

ここで一旦workflowについて広く見てみましょう。ハーバードを始めとする海外の研究所ではサーバーを研究所ごとに管理するのに変えて、クラウドの利用が盛んになってきています。そういう時代感もあり、workflowやコンテナなど、reproducibleな環境の整備というのはますます重要視されてきています。workflow業界ではcwlの他に古き良き[Galaxy](http://wiki.pitagora-galaxy.org/wiki/index.php/Workflows)、Broad Instituteが開発している[WDL](https://software.broadinstitute.org/wdl/)や[nextflow](https://www.nextflow.io/)、Pythonで書かれている[Snakemake](https://snakemake.readthedocs.io/en/stable/)など、さまざまな選択肢があります。どのworkflowが好きか、Cambridgeの[@AlbertVilella](https://twitter.com/@AlbertVilella)が世界中のbioinformaticianにアンケートをとっていました。

<blockquote class="twitter-tweet" data-lang="ja"><p lang="en" dir="ltr">90 responses in <a href="https://t.co/hbaY7ShocR">https://t.co/hbaY7ShocR</a> <a href="https://twitter.com/nextflowio?ref_src=twsrc%5Etfw">@nextflowio</a> <a href="https://twitter.com/SBGenomics?ref_src=twsrc%5Etfw">@SBGenomics</a> <a href="https://twitter.com/commonwl?ref_src=twsrc%5Etfw">@commonwl</a> Snamemake <a href="https://t.co/1kc1DaXdBa">pic.twitter.com/1kc1DaXdBa</a></p>&mdash; Albert Vilella (@AlbertVilella) <a href="https://twitter.com/AlbertVilella/status/1070219898306084865?ref_src=twsrc%5Etfw">2018年12月5日</a></blockquote>


<blockquote class="twitter-tweet" data-lang="ja"><p lang="en" dir="ltr">No votes so far for <a href="https://twitter.com/ensembl?ref_src=twsrc%5Etfw">@ensembl</a> Hive, Apache Taverna or <a href="https://twitter.com/arvados?ref_src=twsrc%5Etfw">@arvados</a> , even though we know they have plenty of users <a href="https://t.co/Db2el8vPDt">https://t.co/Db2el8vPDt</a> <a href="https://twitter.com/hashtag/Bioinformatics?src=hash&amp;ref_src=twsrc%5Etfw">#Bioinformatics</a> <a href="https://twitter.com/hashtag/Workflows?src=hash&amp;ref_src=twsrc%5Etfw">#Workflows</a> <a href="https://t.co/hbaY7ShocR">https://t.co/hbaY7ShocR</a></p>&mdash; Albert Vilella (@AlbertVilella) <a href="https://twitter.com/AlbertVilella/status/1070579264192462849?ref_src=twsrc%5Etfw">2018年12月6日</a></blockquote>

どうやら一番twitter界隈のbioinformaticianはnextflowが好きなようです。ただ、下のspreadsheetを見ると、githubではGalaxy,bcbioについで三番目に盛んなようです。

他にも、cwlやnextflowを取り上げたgenome解析の論文も出ていたりします。

> [Baichoo, S. et al. Developing reproducible bioinformatics analysis workflows for heterogeneous computing environments to support African genomics. BMC Bioinformatics 19, 457 (2018).](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2446-1)

## 結論　みんなでCWLを使おう！

医学部生にできるだけわかりやすく書いてみました。みんなでcwlを使ってどんどん解析を楽にしましょう！

国試も自動化できたらなあ。。。

## 参考資料

- [CWL公式ページ](https://www.commonwl.org/) : リンク集がついている。
- [cwl-intro-gui-workshop スライド](https://docs.google.com/presentation/d/13K6BKQoimiaOIFSnComw8in3HlCxEm820OASP8_jnDI/edit?usp=sharing) : 二階堂研石井さんのスライド。わかりやすいです。
