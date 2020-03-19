Title:Kallistoを用いたRNA-seq解析パイプライン
Date: 2019.04.14
Tags: bioinformatics
Slug: kallisto_rnaseq_pipeline
Author: 平岡
Summary:

!!!!!!!!!!! 
今回はKallistoを用いたRNA-seq解析パイプラインを紹介しま
す。<a href="https://pythonoum.wordpress.com/2018/10/16/linux%E3%81%A7%E3%81%AEbioinformatics%E7%92%B0%E5%A2%83%E6%A7%8B%E7%AF%89_01/">LinuxでのBioinformatics環境構築_01</a>でこの記事への準備はすべて終了している流れになります。<a href="https://ja.wikipedia.org/wiki/%E3%82%AB%E3%83%AA%E3%82%B9%E3%83%88_(%E5%B0%8F%E6%83%91%E6%98%9F)">Kallisto</a>は小惑星の名前のようです。つっこみどころありましたら、コメントいただけると嬉しいです！それではいきましょう！

<strong><em>リファレンスのダウンロード</em></strong>
kallistoでは、transcriptにシュードアラインメントするので、リファレンスにはcDNAを用います。今回は<a href="https://www.gencodegenes.org/">GenCodeGenes</a>のヒトtranscript sequencesのデータを用いました。

[code lang="text"]
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz
[/code]

<hr />

<strong><em>FASTQファイルをダウンロードする場合</em></strong>
<a href="https://www.ebi.ac.uk/arrayexpress/">ArrayExpress</a>からFASTQファイルをダウンロード、解凍する。今回のデータは、ヒトES細胞と成熟膵島細胞のデータ。single-end readとなっている。

[code lang="text"]
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR266/ERR266335/ERR266335.fastq.gz
$ gunzip ERR266335.fastq.gz
[/code]

ERR266349 ERR266351 ERR266338 ERR266347についても同様に！

参照：次世代シークエンサーDRY解析教本 (細胞工学別冊)

<hr />

<strong><em>SRAファイルをダウンロードする場合</em></strong>
DDBJの<a href="https://ddbj.nig.ac.jp/DRASearch/">DRA search</a>からSRAファイルをダウンロード、SRAファイルをFASTQに変換する。pfast-dumpで .sraをペアエンド.fastqに変換。 (<strong>kallistoはsraファイルを扱えない</strong>ので、pfastq-dumpでfastqに変換する必要がある。）

[code lang="text"]

#!/bin/bash
# download sra files.
mkdir sra-fastq
id=(ERR266335 ERR266337 ERR266338 ERR266347 ERR266349 ERR266351)
for item in ${id[@]}
do
echo start download ${item}.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/ERX/ERX182/ERX182652/${item}/${item}.sra
pfastq-dump -s ${item} -t 8 -O sra-fastq
done
[/code]

pfastq-dumpのオプション
-s: SRAファイルのID
-t: スレッド数
-O: 出力ファイル

今回はSRAファイルのダウンロードとpfastq-dumpを使ってsraをfastqに変換する処理を<a>シェルスクリプト</a>を使って行いました。共通項のあるかつ時間のかかるterminalでの処理はシェルスクリプトを使うと便利です。

<hr />

SRAファイルとは？

Sequence Read Archiveの略。（かつてはNGSにリードが短い特徴があったのでShort Read Archiveと呼ばれていた。）NGSの登場により配列の品質情報を塩基配列とともに記述形式であるFASTQ形式が使用されるようになった時にできたバイナリ形式のデータフォーマット。よってpfastq-dumpなどのツールでFASTQに変換することができる。<a href="https://en.wikipedia.org/wiki/International_Nucleotide_Sequence_Database_Collaboration">INSDC</a>、<a href="https://en.wikipedia.org/wiki/European_Bioinformatics_Institute">EBI</a>、<a href="https://en.wikipedia.org/wiki/DNA_Data_Bank_of_Japan">DDBJ</a>が共同で運営しているデータベース<a href="https://en.wikipedia.org/wiki/Sequence_Read_Archive">SRA</a>に保存してある。

[code lang="text"]
$ less ERR266335.sra
[/code]

バイナリデータなので、SRAファイルの中をのぞいてみると下のようになる。
<img class="alignnone size-full wp-image-308" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-19-29-38.png" alt="Screen Shot 2018-10-31 at 19.29.38.png" width="488" height="250" />

<hr />

pfastq-dumpとは？
fastq-dumpを並列処理するbashスクリプト。Sequence Read Archive（wiki）からダウンロードされたシーケンスデータ（SRAフォーマット ）をfastq-dumpの並列処理で素早くfastqに変換することができる。<a href="https://github.com/inutanoOhta">Ohta</a>さんが公開されている。

<hr />

<strong><em>Fastqc</em></strong>
クオリティチェック。

[code lang="text"]
#!/bin/bash
id=(ERR266335 ERR266337 ERR266338 ERR266347 ERR266349 ERR266351)
mkdir fastqc
for item in ${id[@]}
do
echo start quality check ${item}
mkdir fastqc/fastqc_${item}
fastqc -t 8 -o fastqc/fastqc_${item} sra-fastq/${item}.fastq -f fastq
done

[/code]

sra-fastq/${item}.fastqがインプットファイル。
fastqcのオプションについて
-t: スレッド数。
-O: 解析結果の保存先のディレクトリを指定する。今回はfastqcというディレクトリを作ってそこに入れている。
-f: インプットファイルのフォーマット。bam, samにも対応。.fastq.gzもfastqで指定する。

実行結果は下記のようになる。
<img class="alignnone size-full wp-image-309" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-19-58-17.png" alt="Screen Shot 2018-10-31 at 19.58.17.png" width="355" height="253" />

fastqcディレクトリをのぞいてみる。fastqcのなかにfastqc_${item}というディレクトリが自動生成されている。

<img class="alignnone size-full wp-image-310" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-20-05-02.png" alt="Screen Shot 2018-10-31 at 20.05.02.png" width="829" height="32" />

fastqc_ERR266335の中をのぞいてみると。htmlファイルとzipファイルが生成されている。
<img class="alignnone size-full wp-image-311" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-20-07-03.png" alt="Screen Shot 2018-10-31 at 20.07.03.png" width="358" height="25" />

htmlファイルをブラウザでいることができる。

[code lang="text"]
$ open ERR266335_fastqc.html
[/code]

<img class="alignnone size-full wp-image-312" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-20-09-37.png" alt="Screen Shot 2018-10-31 at 20.09.37.png" width="1013" height="554" />

参照：次世代シークエンサーDRY解析教本

<hr />

<strong><em>Multiqc</em></strong>
クオリティチェックの結果、ログファイルなどをまとめていい感じにレポートにしてくれるツール。

[code lang="text"]
# After analysis, run Multiqc by commands below. You can create report.
$ multiqc .
$ open multiqc_report.html
[/code]

multiqcの実行により、関連ファイルが下記のように自動生成される。
<img class="alignnone size-full wp-image-314" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-22-at-21-20-18.png" alt="Screen Shot 2018-10-22 at 21.20.18.png" width="627" height="138" />
複数のリードのクオリティチェックの結果を同時に表示できる。
<img class="alignnone size-full wp-image-313" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-31-at-20-02-17.png" alt="Screen Shot 2018-10-31 at 20.02.17.png" width="1065" height="466" />
参照：https://multiqc.info/

<hr />

<strong><em>Trimmomatic</em></strong>
Java で書かれているアダプタートリミングツールである。 Trimmomatic はアダプターの除去のみならず、リードの末端から一定数の塩基をトリムしたりする、簡単なクオリティフィルタリングも行える。

[code lang="text"]
$ trimmomatic SE -phred33 ERR266335.fastq output_ERR266335.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
[/code]

SEオプションでsingle-end readを指定している。pair-end readでは、オプションでPEとかき、paired outputとunpaired outputの出力先２つを指定する必要がある。

ILLUMINACLIP: 除去するアダプター配列をFASTA形式で与える。そのあとにミスマッチ許容数、palindrome clip threshold、simple clip thresholdの順に指定していく。

参照
https://bi.biopapyrus.jp/rnaseq/qc/trimmomatic.html
http://www.usadellab.org/cms/?page=trimmomatic

<hr />

<strong><em>Coffee Break</em></strong>
<strong>Single-end</strong>, <strong>Pair-end</strong>ってなに？？
<img class="size-full wp-image-233" src="https://pythonoum.files.wordpress.com/2018/10/img_1051.jpg" width="638" height="479" /><img class="size-full wp-image-232" src="https://pythonoum.files.wordpress.com/2018/10/img_1050.jpg" width="689" height="445" />
シーケンスする機器によって、cDNAのかた方のみ読む(single-end read)方法と両端を読む(pair-end read)方法があります。トリミング、アラインメントにおいて、single-endなのか、pair-endなのかを指定しあげる必要があります。

<hr />

<strong><em>Kallisto</em></strong>
RNA-Seqデータ、またはより一般的にはハイスループットシーケンシングリードを用いて転写産物の量を定量化するためのプログラムである。

kallisto や Salmon を利用して定量したデータを使って、edgeR や DESeq2 などで発現量の群間比較を行うことができる。この際に、Bioconductor の tximport パッケージを利用することで、簡単に kallisto/Salmon の定量結果を edgeR/DESEq2 に渡すことができる。

[code lang="text"]
$ time kallisto/kallisto index -i hsGRCh38_kallisto Homo_sapiens.GRCh38.rna.fa.gz

$ time kallisto/kallisto quant -i hsGRCh38_kallisto sra_fastqc/ERR266335.fastq -o ERR266335exp_kallisto
[/code]

kallisto quantにおいて、<strong>-iと-oのオプションは強制</strong>である。
-i:作成したインデックスの指定
-o:出力結果の保存先
<strong>デフォルトではペアエンドを読もうとする</strong>ので、シングルリードの場合は--singleオプションをつける。シングルの時は
-s:Estimated standard deviation of fragment lengthシーケンシング用のライブラリー中のフラグメントの長さの偏差
-l:Estimated average fragment lengthシーケンシング用のライブラリー中のフラグメントの長さの平均
のオプションを追加するのが必須となる。
(kallistoはsraファイルを扱えないので、pfastq-dumpでfastqに変換する必要があった。)
ERR266337 ERR266349 ERR266351 ERR266338 ERR266347も同様に

kallisto_quant.shというシェルスクリプトを書き実行した。今回に限らず、時間がかかるかつ繰り返しの処理はシェルスクリプトを書くと良い（私もこれから練習します）。

[code lang="text"]
#!/bin/bash
id=(ERR266335 ERR266337 ERR266349 ERR266351 ERR266338 ERR266347)
for item in ${id[@]}
do
echo start mapping ${item} with Kallisto
result_dir=${item}_exp_kallisto
kallisto/kallisto quant -i hsGRCh38_kallisto -o ${item} --single -l 200 -s 20 -b 100 sra_fastqc/${item}.fastq
done
[/code]

abundance.tsv, target_id, length, eff_length, est_counts, tpm

参照：
https://scilifelab.github.io/courses/rnaseq/labs/kallisto
https://bi.biopapyrus.jp/rnaseq/mapping/kallisto/kallisto-single-end-reads.html
http://kazumaxneo.hatenablog.com/entry/2018/07/14/180503
https://scilifelab.github.io/courses/rnaseq/labs/kallisto

<hr />

<strong><em>Coffee Break</em></strong>

<em>FPKM</em>,<em>RPKM</em>,<em>TPM</em>とは？
転写産物にマッピングされるリードの数は、サンプル中の総リード数（sequence depth）と転写産物の長さに影響されるので、RNA-Seq データから得られたリードカウントデータは、そのまま転写産物（遺伝子）発現量を表すわけではない。そのため、RNA-Seq データから得られるリードカウントデータを転写産物発現量として利用するには、総リード数や転写産物長で補正する必要がある。

補正計算として、かつてはFPKM,RPKMが用いられてきたが、現在ではかわりにTPMが用いられている。TPMではサンプルごとの値の合計が同じになるので、比較する目的のためにはTPMの方が都合が良い。

FPKM/RPKM の計算

FPKM: Fragments Per Kilobase of exon per Killion reads Mapped

RPKM: Reads Per Kilobase of exon per Million mapped reads
<img class="alignnone size-full wp-image-265" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-23-at-0-40-05.png" alt="Screen Shot 2018-10-23 at 0.40.05.png" width="341" height="59" />

N: リファレンスにマッピングできた全リード数

Yi: そのうち転写産物 i の領域にマッピングされたリード数

Li: 転写産物 i の長さ

TPMの計算

TPMの計算

TPM: Transcripts Per Kilobase Million

<img class="alignnone  wp-image-297" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-28-at-21-48-15.png" alt="Screen Shot 2018-10-28 at 21.48.15.png" width="133" height="92" />

Yt : 転写産物 t にマッピングされたリードカウント

Lt:  転写産物 t の長さ

Tt: 転写産物 t の 1,000 bp あたりのリード数

<img class="alignnone  wp-image-296" src="https://pythonoum.files.wordpress.com/2018/10/screen-shot-2018-10-28-at-21-48-28.png" alt="Screen Shot 2018-10-28 at 21.48.28.png" width="197" height="66" />

転写産物長による補正後の総リードカウントが 100 万となるように補正

参照：
<a href="https://bi.biopapyrus.jp/rnaseq/analysis/normalizaiton/fpkm.html">bipapyrus fpkm</a>
<a href="https://bi.biopapyrus.jp/rnaseq/analysis/normalizaiton/tpm.html">bipapyrus tpm</a>

<hr />

* Tximport, DESeq2を用いた解析はコメントをいただき、ただ今編集中となっております。少々お待ちください。
