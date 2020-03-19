Title:VELOCYTO
Date: 2019.04.14
Tags: bioinformatics
Slug: velocyto
Author: 廣瀬
Summary:

Velocyto は、RNAseq の結果に含まれるイントロンの割合からその細胞の分化指向性を
算出するという解析手法です。

<h3>RNAseq でイントロン??</h3>

RNAseq においては totalRNA のうち 99%ともいわれる rRNA を除き mRNA のみを効率
よく Sequence する目的で、Poly(A)で成熟 mRNA を濃縮することがしばしば行われてい
ます。しかし、そのようなサンプルにおいても、実際の Sequence 結果ではイントロンに
あたる配列が読まれてくることが指摘されており、この原因は Poly(A)類似モチーフの存
在であると推察されています。

結果の例
概日時間の検討ではとてもきれいな結果が得られています。関連遺伝子の unspliced(u) とspliced(s)を観察していくと、circadian timeの経過とともにu がs に変化していきます。
また、未分化細胞からの分化の方向性を見積もることも可能と報告されています。マウス
のオリゴデンドロサイト前駆体細胞がオリゴデンドロサイトへ分化する方向に Velocity を
持っていることが示されました。

実践編
Developer の WEB(http://velocyto.org/)の tutorial を参考に手持ちの BAM で挑戦しま
す。R 用と Python 用とあるようですが、もちろん Python です。二部構成となっており、
前半はCommand lineでBAMファイルから.loomファイルを作成、後半では作成した.loom
ファイルから解析を行うとのことです。

準備

[code lang="text"]
$ conda install numpy scipy cython numba matplotlib scikit-learn h5py click
$ pip install velocyto
$ pip install scanpy
$ pip install -U scvelo
[/code]

※pysam は Windows 環境下では仮想環境下であってもインストールできないそうです。
(https://qiita.com/chaoi/items/6d7702cd70430610f844)これに気付かず数日はまりまし
た。。
Genome annotation file のダウンロード
GENCODE のサイトから GTF ファイルをダウンロードし、適当な場所に保存して
gunzip。
実行
Sequence の手法ごとに Command を選べます。今回は smart-seq2 で行います。

[code lang="text"]
$ velocyto run-smartseq2 -e hogehoge /*.bam annotation.gtf
[/code]

オプション
-o, --outputfolder ¶:出力ディレクトリを作成・指定
-e, --sampleid ¶:出力ファイル名
-m, --repmask ¶:リピート配列をマスクする場合
-t, --dtype ¶:出力ファイルのデータ型(default は unit32)
-d, --dump ¶
-v, --verbose¶
サンプルごとに BAM があることを想定して複数ファイルを*で受け付けてくれます。
BAM と GTF ファイルは必須です。なお、BAM の代わりに SAM を入れてもエラーは吐か
ないようです。これにより、hogehoge.loom ファイルが得られました。この中には、各分
子のスプライシング状態が格納されています。
※TOPHAT でアラインメントしたデータはそのままではうまく行かないようです。
→TopHat-Recondition で unmap のデータを拾ってくればいけるかもしれません。今回は
HISAT2 でアラインメントし直してから velocyto をかけ、改めて hogehoge.loom ファイル
を得ました。
描画
Scanpy と組み合わせて使用するために scvelo を用いました。

[code lang="text"]
import scanpy as sc
import scvelo as scv
[/code]

データ読み込み

[code lang="text"]
adata = sc.read("hogehoge.loom", cache=True)
[/code]

次元圧縮

[code lang="text"]
sc.tl.pca(adata, n_comps=10)
[/code]

速度計算

[code lang="text"]
scv.pp.moments(adata)
scv.tl.velocity(adata)
[/code]

描画

[code lang="text"]
scv.pl.velocity_embedding(adata)
[/code]

<img src="https://pythonoum.files.wordpress.com/2019/03/picture1.png" alt="Picture1.png" width="359" height="351" class="alignnone size-full wp-image-506" />

※それぞれの細胞における特定遺伝子の発現レベルを色で表したり、PCA 以外の次元圧縮
法を用いることもできます。この場合は tSNE のほうがきれいに見えるようです。パラメ
ータはたくさん用意されているので、各データセットで最適な条件を検討する必要がある
でしょう。

<img src="https://pythonoum.files.wordpress.com/2019/03/picture2.png" alt="Picture2.png" width="415" height="361" class="alignnone size-full wp-image-507" />

<img src="https://pythonoum.files.wordpress.com/2019/03/picture3.png" alt="Picture3.png" width="416" height="362" class="alignnone size-full wp-image-508" />

https://www.nature.com/articles/s41586-018-0414-6
http://velocyto.org/velocyto.py/index.html
catway.jp/bioinformatics/qc/rmrepeat.html
https://github.com/theislab/scvelo
