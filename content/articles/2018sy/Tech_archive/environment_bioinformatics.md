Title:LinuxでのBioinformatics環境構築_01
Date: 2019.04.14
Tags: bioinformatics
Slug: environment_bioinformatics
Author: 平岡
Summary:

お疲れ様です。ただいま、私のMacbook Proが入院しておりまして、古いWindows10をUbuntuとデュアルブートして作業しております。下級生でもWindowsしか持っていない、でもBioinformaticsに関心があるという人が、スムーズに環境構築できるようにと今回の記事を書きます。なお、前提として、Ubuntuのインストールが完了しているものとします。なお、筆者のインストールしたUbuntuは18.04.1 LTSです。

今回は<strong>＜<em>Kallisto</em>を用いた<em>RNA-seq</em>解析パイプライン＞</strong>で使う、RNA-seq系のツールのインストールを行っていきますが、順次別の解析目的のツールインストールも紹介したいと考えております。今回の環境構築で、<a href="https://pythonoum.wordpress.com/2018/10/16/kallisto%E3%82%92%E7%94%A8%E3%81%84%E3%81%9Frna-seq%E8%A7%A3%E6%9E%90%E3%83%91%E3%82%A4%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3/">Kallistoを用いたRNA-seq解析パイプライン</a>に進むことができます。

<strong>Ubuntuではpython2.7がデフォルト</strong>となっているので、python3をダウンロードし、デフォルトに設定しよう。

[code lang="text"]
# check the version of python.
$ python --version
>>> Python 2.7.15rc1
[/code]

<a href="https://www.anaconda.com/download/#macos">Anaconda</a>のインストーラーをダウンロードします。

```
$ bash ~/Downloads/Anaconda3-4.1.0-Linux-x86_64.sh
$ echo 'export PATH=/home/user/anaconda3/bin:$PATH' &gt;&gt; ~/.bashrc
$ conda -V
```

[code lang="text"]
$ python --version
>>> Python 3.6.6+
[/code]

バイオインフォマティクスの分野では解析のために様々なツールを利用しますが、インストールのたびにパスを通すなどの作業をしていると大変煩雑ですし、バージョン管理もしにくくなります。パッケージマネージャーを使ってツールを一元管理するのが賢明です。
Homebrewなどのパッケージ管理システムは有名ですが、対応していないツールの多く、現状はBiocondaというパッケージマネージャーがおすすめです。

[code lang="text"]
# Download a Miniconda package for linux python3
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
[/code]

EnterやYesを入力していくとインストールが終了し、パスも通った状態になります。

[code lang="text"]
# condaコマンドが正常に動作していれば成功です。（terminalを開き直しましょう。）
$ conda -h
[/code]

biocondaのチャンネルを追加。biocondaチャンネルが最上位に来るように設定。

[code lang="text"]
$ conda config --add channels conda-forge
$ conda config --add channels defaults
$ conda config --add channels r
$ conda config --add channels bioconda
[/code]

それでは早速ツールをインストールしていきましょう。基本 conda install ???でインストールが完了し正常に動作していきます。非常に便利ですね。

[code lang="text"]
$ conda install parallel-fastq-dump
$ conda install fastqc
$ conda install multiqc
$ conda install trimmomatic
$ conda install kallisto
[/code]

参照：http://imamachi-n.hatenablog.com/entry/2017/01/14/212719

<strong><em>R</em>ツールのインストール</strong>
Rのコンソールを開いて、Rのツールをインストールしていきます。
Tximportのインストール

[code lang="text"]
> source("https://bioconductor.org/biocLite.R")
> biocLite("tximport")
[/code]

DESeq2のインストール

[code lang="text"]
> source("https://bioconductor.org/biocLite.R")
> biocLite("DESeq2")
[/code]

参照：
https://bioconductor.org/packages/release/bioc/html/tximport.html
https://bioconductor.org/packages/release/bioc/html/DESeq2.html
