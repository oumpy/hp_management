Title: WSI＊自己教師あり学習のメモ
Date: 2022.07.03
Tags: Machine Learning,
Author: 安部

最近病理の画像をバイトで扱うことが増えています。せっかくなので自己教師あり学習を使ってみたいと思い、病理画像に自己教師あり学習を用いた論文をいくつか調べました。

#### WSIとは

Whole slide images（WSI）とは，バーチャルスライドとも呼ばれる病理組織プレパラート標本全体をスキャンしたものです。画像サイズが巨大なので128*128などのパッチに分割して扱うことが多いです。パッチの作成時に組織以外の余白を除くことが多いです。


#### 自己教師あり学習とは

ラベルなしデータに対して、データ自身から独自のラベルを機械的に作成したものから画像の表現を獲得するタスクです。
大まかに以下の3種類に分けることができると著者は考えています。
1. 画像の一部をmaskして隠された部分を再構成できるように学習する
2. 画像Aとそれに摂動を加えた画像A'、別画像Bを用意して、それぞれを埋め込んだときに似ているAとA'が近く、異なるBは遠くになるように学習する
3. 画像Aに異なる摂動を加えたA'とA''を用意して別々に埋め込みを得た時にA'とA''が似ていることを学習する
例として、<br>
1. [MAE](https://arxiv.org/abs/2111.06377)
2. [SIMCLR](https://arxiv.org/abs/2002.05709)
3. [DINO](https://arxiv.org/abs/2104.14294)

などがあります。<br>
②の手法では正例と負例の類似度が高すぎるとうまくいかないこと、①は分類よりもセグメンテーションに強い特徴が得られること、③は分類には強いがセグメンテーションに適した特徴が得られなあい可能性があること、が[最近の論文](https://arxiv.org/pdf/2204.07141.pdf)にて指摘されています。

### 　Self supervised contrastive learning for digital histopathology

[原著](https://www.sciencedirect.com/science/article/pii/S2666827021000992)

#### 概要

WSIにSIMCLRを適応。TCGAやCPTAC、コンペデータから得た多臓器の57のデータセットで事前学習させた。
評価は分類、回帰、セグメンテーションなど複数の下流タスクを行なった。

#### 結果
少アノテーションにてimagenetで事前学習させたモデルに勝る結果が得られた。
<img src="{attach}./images/wsi_ssl_figs/path_simclr_fig1.jpg" alt="path_simclr_fig1" width="400px">
教師なしでクラスタリングしてみると、似ている画像が近傍にあることが確認された。<br>
高解像度で事前学習させた方が性能が上がり、異なる解像度を組み合わせることでも性能が上がった。<br>
前立腺画像で学習させたモデルが乳房画像で学習させたモデルよりも乳房画像での分類タスクの性能が良いことがあった。複数のデータセットを組み合わせるのが重要だった。<br>


### Contrastive learning-based computational histopathology predict differential expression of cancer driver genes

[原著](https://arxiv.org/abs/2204.11994)


#### 概要
WSIの扱いとして定番な「patch分割→patchの情報を集約」の手法にてpatchの埋め込み作成に対照学習(AdCo:https://arxiv.org/abs/2011.08435)を使用。集約にはattention poolingを使った。
下流タスクとして腫瘍診断と癌ドライバー遺伝子の差次的発現予測を行なった。

#### 結果
対照学習を使用しない従来手法より性能が向上した。

<img src="{attach}./images/wsi_ssl_figs/histcode_fig1.png" alt="histcode_fig1" width="400px">
attentionを可視化すると病理医と類似した領域に注目していることが確認された。

<img src="{attach}./images/wsi_ssl_figs/histcode_fig2.png" alt="histcode_fig2" width="400px">
T/Bリンパ球に得意的なマーカーについてのアノテーションを追加で作成してみると、腫瘍細胞の場所とリンパ球が集まっていた場所が一致し、免疫細胞が固形腫瘍に浸潤する事実に一致する結果が得られた。
attentionをみることによりどこの組織に何の遺伝子が発現している細胞があるのかを可視化することができる。





###　　Scaling Vision Transformers to Gigapixel Images via　Hierarchical Self-Supervised Learning

[原著](https://arxiv.org/pdf/2206.02647.pdf)

#### 概要

「patch分割→patchの情報を集約」において従来手法ではpatchの埋め込み作成だけに自己教師あり学習を使用していたが、埋め込みの集約の部分にも自己教師あり学習を使用した。自己教師あり学習の手法としてDINOを使用。

<img src="{attach}./images/wsi_ssl_figs/ireko_wsi_fig1.png" alt="ireko_wsi_fig1" width="400px">

ViTのpatch埋め込みにresnet50を使用したハイブリットversionよろしく上層のViTのpatch埋め込みに下層のViTを使用する。下層からDINOを学習していく。



#### 結果


少アノテーションにて従来手法よりも高い性能を得た。

<img src="{attach}./images/wsi_ssl_figs/ireko_wsi_fig2.png" alt="ireko_wsi_fig2" width="400px">
<img src="{attach}./images/wsi_ssl_figs/ireko_wsi_fig3.png" alt="ireko_wsi_fig3" width="400px">
最小patch単位のattentionではheadごと正常細胞、異型細胞、リンパ球、余白(内腔スペース、脂肪領域、エアポケット)といったミクロな特徴を捉えることができていた。また、集約層のattentionを見ると腫瘍の成長パターン、脂肪や間質領域への腫瘍浸潤、その他の組織-組織間の関係性を捉えることに成功していた。

制約として、最後の集約層はスライド単位の学習をしているので十分なデータを集めるのは困難である点や、事前学習を入れ子で行うにはマシンパワーが必要な点が挙げられる。





###　まとめ
patchレベルで自己教師あり学習を用いるとimagenetの事前学習済みモデルを使用したときよりも性能が向上したり、少アノテーションでも性能が担保されたりすることが確認されているようです。
上記論文中ではどの自己教師あり学習を使用すればよいかという比較はされていなかったので下流タスクに応じて幾つか試す必要がありそうです。














