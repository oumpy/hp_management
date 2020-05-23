Title:LinuxでのBioinformatics環境構築
Date: 2018.10.16
Tags: Bioinformatics
Author: 平岡

ただいま、私のMacbook Proが入院しておりまして、古いWindows10をUbuntuとデュアルブートして作業しております。下級生でもWindowsしか持っていない、でもBioinformaticsに関心があるという人が、スムーズに環境構築できるようにと今回の記事を書きます。なお、前提として、Ubuntuのインストールが完了しているものとします。なお、筆者のインストールしたUbuntuは18.04.1 LTSです。

今回は[Kallistoを用いたRNA-seq解析パイプライン](https://oumpy.github.io/articles/2018/10/kallisto_rnaseq_pipeline.html)で使う、RNA-seq系のツールのインストールを行っていきますが、順次別の解析目的のツールインストールも紹介したいと考えております。今回の環境構築で、Kallistoを用いたRNA-seq解析パイプラインに進むことができます。

## Pythonのインストール
Ubuntuのバージョンによってはpython2.7がデフォルトとなっている場合もあるので、python3系をダウンロードし、デフォルトに設定しよう。

```bash
# check the version of python.
$ python --version
>>> Python 2.7.15rc1
```

## Biocondaのインストール
<a href="https://www.anaconda.com/download/#macos">Anaconda</a>のインストーラーをダウンロードします。

```bash
$ bash ~/Downloads/Anaconda3-4.1.0-Linux-x86_64.sh
$ echo 'export PATH=/home/user/anaconda3/bin:$PATH' &gt;&gt; ~/.bashrc
$ conda -V
```

```bash
$ python --version
>>> Python 3.6.6+
```

バイオインフォマティクスの分野では解析のために様々なツールを利用しますが、インストールのたびにパスを通すなどの作業をしていると大変煩雑ですし、バージョン管理もしにくくなります。パッケージマネージャーを使ってツールを一元管理するのが賢明です。
Homebrewなどのパッケージ管理システムは有名ですが、対応していないツールの多く、現状はBiocondaというパッケージマネージャーがおすすめです。

```bash
# Download a Miniconda package for linux python3
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

EnterやYesを入力していくとインストールが終了し、パスも通った状態になります。

```bash
# condaコマンドが正常に動作していれば成功です。（terminalを開き直しましょう。）
$ conda -h
```

biocondaのチャンネルを追加。biocondaチャンネルが最上位に来るように設定。

```bash
$ conda config --add channels conda-forge
$ conda config --add channels defaults
$ conda config --add channels r
$ conda config --add channels bioconda
```

それでは早速ツールをインストールしていきましょう。基本 `conda install ???`でインストールが完了し正常に動作していきます。非常に便利ですね。

```bash
$ conda install parallel-fastq-dump
$ conda install fastqc
$ conda install multiqc
$ conda install trimmomatic
$ conda install kallisto
```

### 参照
- <http://imamachi-n.hatenablog.com/entry/2017/01/14/212719>

## Rツールのインストール
Rのコンソールを開いて、Rのツールをインストールしていきます。
Tximportのインストール

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("tximport")
```

DESeq2のインストール

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("DESeq2")
```

### 参照
- <https://bioconductor.org/packages/release/bioc/html/tximport.html>
- <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>
