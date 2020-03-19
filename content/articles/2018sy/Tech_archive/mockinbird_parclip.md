Title:Mockinbirdを用いたPAR-CLIP解析
Date: 2019.04.14
Tags: bioinformatics
Slug: mockinbird_parclip
Author: 平岡
Summary:

***Mockinbirdを用いたPAR-CLIP解析***
５年　平岡　悠

今回はPAR-CLIP解析のAll-in-oneパイプラインソフトウェアである、<a href="//wwwuser.gwdg.de/~compbiol/mockinbird/doc/intro.html)”">Mockinbird</a>を紹介します。日本語での情報がほぼなく、説明が長くなってしまいそうなので、２回に分けて書きます。

**PAR-CLIP**って何？
PAR-CLIPは photoactivatable ribonucleoside-enhanced crosslinking and immunoprecipitationの略でRNAやmicroRNA結合タンパクやの結合サイトを特定するために使います。

実験法としては、まず4SU（4-チオウリジン）を細胞培地中に添加し、4SUによるラベル標識を行います。次にUVの照射を行い、RNA結合タンパクとRNAの間で架橋反応を起こします。その後、RBPと結合したリードの抽出を行い、サンプルとして用います。

サンプルRNAを逆転写するときに、4SUがシトシン( C )に置換されます(misread)。そのため、リードがピークを形成し、かつT -&gt; C置換が入っていた場合、RBPと結合していた可能性が高いということになります。これによりあるタンパクが、RNAのどの領域に結合するのか？結合領域間でどのような相互作用があるのかを知ることができます。特に一塩基という高い解像度での解析ができることが特徴です。

![](https://lh6.googleusercontent.com/554iMVEXLP3Hq6U6C8FeveBV1FeMT19zXQuj-728Db4UxDU7JDSnzpF-RvezYq0DW3z2kQjDiOlZkeQLJRZhwNm4SqkaeJnBYyOV7yUbmuV0peAaJ4TcKazkJaqvRMF65Rgldk2l)

**Mockinbird**とは？

<a href="//wwwuser.gwdg.de/~compbiol/mockinbird/doc/intro.html”">公式ドキュメント</a>、<a href="//github.com/soedinglab/mockinbird)”">GitHub</a>に詳細は記載されていますが、日本語での情報が少ない、というか全くない。のでまとめていきたいと思います。

PAR-CLIPの解析は生のFASTQデータから、Quality Check -&gt; Trimming -&gt; Mapping -&gt; Downstream Analysisという流れで行いますが、それらすべての解析をMockinbirdでは一気に行ってくれます。また、MockinbirdModule（詳細は後述）を使うことでmock PAR-CLIP experimentを考慮に入れた解析を行うことができるのも大きな特徴です。

このソフトでは、２つにパイプラインが分かれており、preprocessing phaseとpostprocessing phaseと呼ばれています。タンパク結合サイトの情報が書かれたテーブルを出力するまでが、preprocessing phaseでそれ以降のDownstream解析をpostprocessing phaseと呼んでいます。それぞれにYAMLファイルが用意されており、ユーザーはYAMLファイルで解析したいサンプルデータ、使いたいリファレンスデータ、使いたいモジュール（Mappingの時にSTARやBowtieを選べたりする。）閾値などのパラメータを指定します。最後にコマンドを一行実行するだけで、解析がすべて自動で流れていきます。あるモジュールで出力されたファイルのパスが次のモジュールに自動で設定されていくシステムになっているので、慣れると非常に便利です。

ただ、ドキュメントでも、既存の解析手法に完全に取って代わるものではなく、PAR-CLIP実験のコンディションのトラブルシューティングを再現可能性が高い状態で行なったり、研究の最初の仮説づくりを短時間で行うために使ってもらうことが目的と書かれていました。![](https://lh6.googleusercontent.com/rd3A2Nm8czl8cdPu_47SuddQXr-i-ac0ThA-ZupDlDIa67geEayilGTo2Bp3VnJXt6UXNe9b7Rq6EckkNigvg68OCJ9Wk2XrKfQ-FwJu5I_HjkwXXrpMUR7x-0pBojlg6E1EtlH2)

環境構築

前置きが長くなりましたが、実際に環境構築から進めていきます。

```conda create -n mockinbird -c bioconda -c conda-forge python=3.6 mockinbird```

Anacondaでmockinbirdという名前の仮想環境を作ります。mockinbirdが仮想環境にインストールされることになり、これにより、下記のライブラリとツールが環境内で使えるようになります。

```source activate mockinbird```

でmockinbirdの仮想環境を立ち上げ、

```source deactivate mockinbird```

で仮想環境を閉じます。

```conda remove --all -n mockinbird```

で環境を削除します。

```git clone [https://github.com/soedinglab/mockinbird.git](https://github.com/soedinglab/mockinbird.git)```

でgithubからクローンします。

```cd mockinbird/mockinbird/data```

とディレクトリを進めていくと、```preprocess.yaml```, ```postprocess.yaml```というファイルがあります。基本はこの二つのYAMLファイルを調整していくことになります。terminalでコマンドを出すときも ```mockinbird/mockinbird/data```ディレクトリで実行するのが個人的にはおすすめです（ファイルの構成が綺麗になる。パスの調整も楽）。

<a href="//wwwuser.gwdg.de/~compbiol/mockinbird/mockinbird_tutorial_nomock.tar.gz”">Tutorialデータ</a>がドキュメントの方に書かれていましたが、私が見たときはリンクがNotFoundになっていました。ので、私は<a href="//wwwuser.gwdg.de/~compbiol/mockinbird/doc/intro.html”">公式ドキュメント</a>、<a href="//github.com/soedinglab/mockinbird)”">GitHub</a>を参考にYAMLファイルの調整などを行いました。

```preprocess```, ```postprocess```の実行のコードは以下のようになっています。

```
# mockinbird preprocess [parclip_fastq] [output_dir] [prefix] [config_file]
$ mockinbird preprocess nrd1.fastq nrd1 nrd1 preprocess.yaml

# mockinbird postprocess {{genomefasta}} {{output_dir}} {{output_dir}} {{script_dir}}
$ mockinbird postprocess nrd1 nrd1_pp postprocess.yaml
```

ですが、コードの実行はYAMLの設定ができてからになります。ということでYAMLの中身について説明していきます。

***preprocess.yaml***

preprocess.yamlの４つの区画から構成されています。

① 変数設定。
② general. 必須情報の設定。
③ reads. リードについての情報。
④ pipeline. 使うモジュールと各引数の設定。

① 変数設定
下記のスクショはドキュメントから引用していますが、モジュールで使う、ディレクトリの変数化と、mock_processing をFalseに設定しています。（ここ重要！）

![](https://lh5.googleusercontent.com/axBM--SAveLQ5WpApwgiK3VmRPCEywGRY1ap3jXUNrp1ejgabu7k_3yMrnFT5dEA6xO5j2wBdjLh4UFNdJ5zbtwjq4nndYXoBVOBQyRRT7VskRhZUiruBhH8CK-AeFxu4Fv3bX3i)

② general. 必須情報の設定
アダプター配列の指定、リファレンスデータの指定、UMI(Unique molecular identifiers)の有無、スレッド数の指定。

③ reads. リードについての情報
リードの最短長、T -&gt; C mutationの指定

![](https://lh3.googleusercontent.com/STLDXYeF8NC-rqiRjaW4rArxJNHYNNcpDkrpj-OkX6W0Ltr7eo0FkqMd6VXWM6asWlMQMDhkg0l0pViKfzpFpXKFhN_KN1tyJkLDYeNt2fRqyNMVNQvp6rhwyF0LOZN7W_mJaSxy)

④ pipeline. 使うモジュールと各引数の設定
クオリティチェック -&gt; トリミング -&gt; マッピング
![](https://lh5.googleusercontent.com/KASBpB_VhdJ80ZHeIyOucHGZShtPan6lgSWzCcrKl6qGoQahRudn5Gq9n5gu--xl2o3bWuOh7u8y4huuc6adxcis-BLRJhW9iczYI2tREs8MRWkpzoWsx_WtQXIvtlxskNRaZK2h)

![](https://lh4.googleusercontent.com/7vAyDIsxenmYiwISI6p1k6mjbbw1HESriPeyla17MOvM8dPjJxYxD5cLF8OWPXtO7AHCt4AvDjXPZmndQ7p6sNUkME4BIdM-CsmRfKyf0L7ni1nQY3LXoHXDdv6cNY9Xg_EGQyWb)
MockinbirdModuleを使うことでmock PAR-CLIP experimentを考慮に入れた解析ができますと、すでに書きましたが、Mock Experimentの結果をバックグラウンドとして処理して、PAR-CLIP Experimentの結果を出すことができます。

具体的な方法としては、まず
```mock_processing = True```とし、
```mockinbird preprocess nrd1_mock.fastq nrd1_mock nrd1_mock preprocess.yaml```を実行します。これで```{% if not mock_processing %}```までで解析がストップします。

PileupModule, BamStatisticsModuleによって、Mock Experimentの.mpileupファイルと_stat.jsonファイルがそれぞれ生成されます。（この二つのファイルを次に使うことになります。）

次に```mock_processing = False```に戻し、
```$ mockinbird preprocess nrd1.fastq nrd1 nrd1 preprocess.yaml```を実行します。すると下記のコードも一緒に実行されます。

![](https://lh4.googleusercontent.com/hKqjVPmLLGnF0bhn637a_19_gZZ6KJBDceofpi8TVC_NfnhM1rZiZVb8MiBz-TmVv0FtStVpck_3IdL_CmYS52kpQqm0uzGbDWK5xyw1InCRVyLE8ne99FQIuAFm4VkuvpEkEtC-)

これによって、Mock Experimentの``.mpileupファイル``と``_stat.jsonファイル``を受け取りつつ、最終的なCLIP-seqのテーブルを作ることができ、preprocessは終了となります。下記のようなテーブルが出力され、これをもとにpostprocessに進んでいくことになります。
![](https://lh3.googleusercontent.com/yb1p70EaelwH-dEGIRhTLnjUoDLgKHSNzqD6dpG9IT9y5fM0LApdvtHKNgmKA73G6HgfflmWkEnoFvm_HEjpIxf4T6Uk15uG4ovIJsG2ILHRGGpXELLQ381BWMkFU6k8AxtSgocl)![](https://lh4.googleusercontent.com/3bwU07uRrmxRkOJxVPU8N-S11uYBRU3NrRFFrkoDyz_RlXjHgClXH2aDfJ2IU10yV-kURjSGBOR5nsdyaUBrh6C2lkMW3u8C128ktZXIWXniOzjGyIoD68ArK3awsfazmT4bRr2Z)

長くなってしまったので、モジュールとpostprocessの説明は次回に回したいと思います。長い間お付き合いいただきありがとうございました。
