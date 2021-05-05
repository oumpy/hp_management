Title: 【論文まとめ】DINO: Emerging Properties in Self-Supervised Vision Transformers
Date: 2021.05.07
Tags: Machine Learning, 論文まとめ
Author: 山本

本記事では主にFacebook AI Researchの研究者らによって提案された**DINO**という画像モデルにおける自己教師あり学習の解説を行います．

> Caron, Mathilde, Hugo Touvron, Ishan Misra, Hervé Jégou, Julien Mairal, Piotr Bojanowski, and Armand Joulin. 2021. “Emerging Properties in Self-Supervised Vision Transformers.” arXiv [cs.CV].  [http://arxiv.org/abs/2104.14294](http://arxiv.org/abs/2104.14294).

(cf.) [Facebook ブログ](https://ai.facebook.com/blog/dino-paws-computer-vision-with-self-supervised-transformers-and-10x-more-efficient-training), [GitHub](https://github.com/facebookresearch/dino), [Yannic Kilcherの解説動画](https://www.youtube.com/watch?v=h3ij3F3cPIk)

**要点**：画像モデル (e.g. ResNet, Vision transformers)における，ラベル無しの自己教師あり学習の新しい手法DINOを考案したよ．ImageNetの画像をDINOで学習されたモデルを用いて特徴空間に埋め込むと，物体の種類ごとにクラスターが生まれ，線形変換あるいはkNNを適応するだけで教師あり学習に匹敵する精度の物体認識ができたよ．特に画像モデルにVision Transformersを用いたとき，そのAttention mapは物体を識別し，segmentationのようなことができていたよ．

![fig1]({attach}./images/dino_figs/fig1_dino.jpg)


## 背景
### Self-supervised learning
**自己教師あり学習 (Self-supervised learning; SSL)** は機械学習モデルの訓練方法の一種です．自己教師あり学習では，欠損のないデータを用意し，データの一部分（欠損したデータ）だけをモデルに入力し，残りのデータを予測するようにモデルを訓練します．下図は自己教師あり学習の概念図です ([Self-supervised learning: The dark matter of intelligence](https://ai.facebook.com/blog/self-supervised-learning-the-dark-matter-of-intelligence/)から引用)．

![ssl]({attach}./images/dino_figs/ssl.png)

データ全体を大きな直方体に見立てたとき，入力されるデータは緑色の部分，予測するデータは残りの灰色の部分となります．自己教師あり学習の例として，言語モデルではBERT ([Devlin et al., 2018](https://arxiv.org/abs/1810.04805)) のように文の一部を隠して，残りの単語から隠された単語を予測するような学習法があります．画像モデルでは，画像の一部を塗りつぶして残りの部分から塗りつぶされた部分を予測するInpainting (e.g. DeepFill; [Yu et al., CVPR. 2018](https://arxiv.org/abs/1801.07892); [Yu et al., CVPR. 2019](https://arxiv.org/abs/1806.03589))や，現在のフレームから未来のフレームを予測する動画予測 (e.g. PredNet; [Lotter et al., ICLR. 2017](https://arxiv.org/abs/1605.08104))などがあります．計算論的神経科学における視覚モデルとしても，こうした自己教師あり学習/教師なし学習のモデルは生物学的妥当性 (biologically plausible)があるとされて研究が進められています ([Zhuang et al., PNAS. 2021](https://www.pnas.org/content/118/3/e2014196118.long))．

この記事で紹介するDINOは**自己教師あり表現学習 (Self-supervised representation learning)** の一種と言えます．自己教師あり表現学習には多数の研究がありますが，2つの系統について説明します．まず，**対照学習 (contrastive learning) **を用いた手法としては，MoCo ([He et al., 2019](https://arxiv.org/abs/1911.05722)), SimCLR ([Chen et al., 2020](https://arxiv.org/abs/2002.05709))などが代表的です．対照学習を用いた手法では，元画像とそれをaugmentした画像 (positive sample)，関係ない画像 (negative sample)の3つの画像を入力し，元画像とpositive sampleとの表現の類似度が小さく，negative sampleとの表現の類似度が大きくなるように学習を進めます．次に，2020年後半頃より，**非対称的な2つのネットワークを用いる手法**が複数考案されてきました．**BYOL** ([Grill et al., 2020](https://arxiv.org/abs/2006.07733)), SimSiam ([Chen & He, 2020](https://arxiv.org/abs/2011.10566)), Barlow Twins ([Zbontar et al., 2021](https://arxiv.org/abs/2103.03230)) などが該当します．これらの手法ではaugmentationを用いた非対称的な入力と，stop-gradientを用いた非対称的な重み更新による2つの非対称性を生み出し，2つのネットワークの出力の類似度を小さくするように学習を行います．この記事で紹介するDINOはBYOLに似た手法ですが，以下の3つが異なります．

1. BYOLとは異なる類似度損失を用いている．
2. 2つのネットワークが同じモデル構造をしている．
3. 手法をVision Transformerにも適応している．
  
### Vision Transformers
**Vision Transformer (ViT）** は言語モデルにおいてRNNを置き換えるモジュールとして考案されたTransformerを画像モデルに適応したモデルのことです ([Dosovitskiy et al., ICLR. 2021](https://openreview.net/forum?id=YicbFdNTTy))．ViTはビタミンの略称でないことに注意しましょう．

ViTでは，入力画像をパッチに分割（基本的に16×16のような格子状に分割）し，画像パッチをベクトルに埋め込んだ後，Transformerのencoderに入力するという操作をします（下図参照: Dosovitskiy et al., 2021; Fig.1）．

![vit]({attach}./images/dino_figs/vit.jpg)

この記事ではViTに関して詳しい解説はしませんが，以下の資料が参考になります．

- [Transformer メタサーベイ - Slideshare](https://www.slideshare.net/cvpaperchallenge/transformer-247407256) 
- [画像認識の大革命。AI界で話題爆発中の「Vision Transformer」を解説！ - Qiita](https://qiita.com/omiita/items/0049ade809c4817670d7)
  
## DINO : self-*di*stillation with *no* labels
### モデル構造
Caronらが提案するDINO (self-**di**stillation with **no** labels) とは「**ラベル無しでの自己蒸留**」を意味します．ここでの**蒸留 (distillation) **とは，枝刈り (pruning) や量子化 (quantization)に並ぶニューラルネットワークのモデル圧縮手法です．通常の蒸留では，ラベル付きデータセットとパラメータ数が多いモデル（**教師モデル; teacher model**）を用意し，ラベルと教師モデルの出力を教師信号 (hard & soft target) としてパラメータ数が少ないモデル（**生徒モデル; student model**）を訓練します ([Hinton, Vinyals & Dean, NIPS, 2014](https://arxiv.org/abs/1503.02531))．データセットだけでscratchから生徒モデルを訓練するより，教師モデルを用いた方が生徒モデルの性能は高くなるということが知られています．

DINOは教師モデルと生徒モデルを用いるところは通常の蒸留と同じですが，ラベル付きデータは用いず，**教師モデルは生徒モデルから**作られるという違いがあります．このことを念頭に置いて，モデルの構造を見てみましょう．以下はモデルの概略図です（[Facebook ブログ](https://ai.facebook.com/blog/dino-paws-computer-vision-with-self-supervised-transformers-and-10x-more-efficient-training)より引用および改変）．生徒モデルを$g _ {\theta s}$，教師モデルを $g _ {\theta t}$とし，入力画像を $\mathbf{x}$とします．

![model1]({attach}./images/dino_figs/model1.jpg)

論文中には次のようなPytorch-likeな擬似コードが掲載されています (Caron et al., 2021; Algorithm 1)．以下ではこの擬似コードを用いながら解説を行います．

```python
# gs, gt: student and teacher networks
# C: center (K)
# tps, tpt: student and teacher temperatures
# l, m: network and center momentum rates
gt.params = gs.params
for x in loader: # load a minibatch x with n samples
    x1, x2 = augment(x), augment(x) # random views
    s1, s2 = gs(x1), gs(x2) # student output n-by-K
    t1, t2 = gt(x1), gt(x2) # teacher output n-by-K
    loss = H(t1, s2)/2 + H(t2, s1)/2
    loss.backward() # back-propagate
    # student, teacher and center updates
    update(gs) # SGD
    gt.params = l*gt.params + (1-l)*gs.params
    C = m*C + (1-m)*cat([t1, t2]).mean(dim=0)

def H(t, s):
    t = t.detach() # stop gradient
    s = softmax(s / tps, dim=1)
    t = softmax((t - C) / tpt, dim=1) # center + sharpen
    return - (t * log(s)).sum(dim=1).mean()
```
  
### Augmentation
各モデルに対して$\mathbf{x}$はそのまま入力せず，画像 $\mathbf{x}$ を切り出す(crop)ようなaugmentationした画像を入力します．切り出し方は，入力画像 $\mathbf{x}$から2つのglobal 画像 $x_1^g, x_2^g$，および複数のlocal画像 $x_j^\ell\ (\ell = 1, 2, \ldots)$を生成するようにします．Global画像とlocal画像は切り出す範囲（解像度）が異なり，例えばglobal画像は元画像の50%以上，local画像は元画像の50%未満などとします．さらに切り出した画像全ての集合を $V = \left\{x_1^g, x_2^g, x_j^\ell\right\} \ (\ell = 1, 2, \ldots)$とします．ここで，生徒モデルには集合 $V$の全ての要素，すなわちglobal画像とlocal画像の両方を入力しますが，教師モデルにはglobal画像 $x_1^g, x_2^g$ のみを入力します．こうすることで，localからglobal (local-to-global)への対応関係が生み出されると，Caronらは述べています．

擬似コードで対応する部分は以下のようになっています．この擬似コードでは両方のモデルに`x1, x2`を入力しているので，`x1, x2`はglobal画像であると思われます（誤解している可能性があります）．

```python
for x in loader: # load a minibatch x with n samples
    x1, x2 = augment(x), augment(x) # random views
    s1, s2 = gs(x1), gs(x2) # student output n-by-K
    t1, t2 = gt(x1), gt(x2) # teacher output n-by-K
```
  
### Sharpeningとcentering
モデルの崩壊（e.g. 出力が一様分布化）を避けるためにモデルの出力に**sharpening**と**centering**の2つの操作を行うことが提案されています．生徒モデルにはsharpeningのみ，教師モデルには両方の操作を適応します．

まず，sharpeningは通常使われるSoftmax関数を修正することで実装されます．

$$
P(x)^{(i)}=\frac{\exp \left(g_{\theta}(x)^{(i)} / \tau\right)}{\sum_{k=1}^{K} \exp \left(g_{\theta}(x)^{(k)} / \tau \right)}
$$

ここで$\tau$は温度のパラメータ（Softmaxの元となったカノニカル分布に由来する名称）であり，分布の鋭さ (sharpness)に影響します．DINOでは生徒，教師モデルで異なる温度パラメータ $\tau_{s}, \tau_{t}$をそれぞれ用います．なお，$\tau>1$の場合，分布はむしろ緩やかになるので，$0<\tau<1$の値が用いられます．

次に，centeringは教師モデルの出力から$C$を除算することで出力の分布を平均に近づける操作です．ここで$C$はスカラーではなく，出力と同じサイズのベクトルであることに注意してください（softmaxの入力からスカラーを除算しても同じ値しか出力されません）．$C$はゼロベクトルで初期化され，次のように$g_{\theta_{t}}$の出力の指数移動平均で更新されます．

$$
c \leftarrow m c+(1-m) \frac{1}{B} \sum_{i=1}^{B} g_{\theta_{t}}\left(x_{i}\right)
$$

なお，$B$はバッチサイズです．擬似コードでは次のように実装されます．

```python
C = m*C + (1-m)*cat([t1, t2]).mean(dim=0)

s = softmax(s / tps, dim=1)
t = softmax((t - C) / tpt, dim=1) # center + sharpen
```

それぞれの効果は次のようになります．

![sharpening_centering]({attach}./images/dino_figs/sharpening_centering.png)

上図を描画するPythonコードは以下の通りです．

```python
import numpy as np
import matplotlib.pyplot as plt

def softmax(x):
    ex = np.exp(x - np.max(x))
    return ex / np.sum(ex)

K, B, tau = 10, 5, 0.07 # output dims, batch size, temp param
x, C = np.random.rand(K), np.mean(np.random.rand(5, K), axis=0) # input, center
pos = range(K) # for bar plot
plt.figure(figsize=(12, 3), dpi=100)
plt.subplot(1,4,1); plt.bar(pos, x); plt.title(r"$x  \in \mathbb{R}^K$")
plt.subplot(1,4,2); plt.bar(pos, softmax(x)); plt.title("softmax"+r"$(x)$")
plt.subplot(1,4,3); plt.bar(pos, softmax(x/tau)); plt.title("Sharpening: softmax"+r"$(x/\tau)$")
plt.subplot(1,4,4); plt.bar(pos, softmax(x-C)); plt.title("Centering: softmax"+r"$(x-C)$")
plt.tight_layout()
```
  
### 損失関数とパラメータの更新
以下は損失関数の計算とパラメータの更新の概略図です（[Facebook ブログ](https://ai.facebook.com/blog/dino-paws-computer-vision-with-self-supervised-transformers-and-10x-more-efficient-training)より引用および改変）．

![model2]({attach}./images/dino_figs/model2.jpg)

まず，生徒モデルの出力 $P_s$と教師モデルの出力 $P_t$を用い，損失関数 $H(P_s, P_t):=-P_t\log P_s$を計算します．次に損失関数を最小化するように生徒モデルのパラメータ $\theta_s$ をbackpropで更新します．なお，損失関数および生徒モデルの最適化問題は，Augmentationの節で述べたglobal画像$x_1^g, x_2^g$とlocal画像$x_j^\ell\ (\ell = 1, 2, \ldots)$，および全体の画像集合 $V$を用いると次のように表されます．

$$
\min _{\theta_{s}} \sum_{x \in\left\{x_{1}^{g}, x_{2}^{g}\right\}} \sum_{x^{\prime} \in V \atop x^{\prime} \neq x} H\left(P_{t}(x), P_{s}\left(x^{\prime}\right)\right)
$$


一方，教師モデルのパラメータ $\theta_t$ は$\theta_s$ を**指数移動平均 (exponential moving average; EMA)** することにより生成されます．

$$
\theta_t \leftarrow \lambda \theta_t + (1-\lambda) \theta_s
$$

なお，$\theta_t$の初期値は$\theta_s$とします．これに対応する擬似コードは以下の部分です．

```python
loss = H(t1, s2)/2 + H(t2, s1)/2
loss.backward() # back-propagate
update(gs) # SGD
gt.params = l*gt.params + (1-l)*gs.params
```
  
### 学習後のDINOモデルの特徴
学習後のDINOモデルの特徴としては，次の2点があります．

1. ImageNetの画像をDINOで学習されたモデルを用いて特徴空間に埋め込むと，物体の種類ごとにクラスターが生まれ，線形変換あるいはkNNを適応するだけで教師あり学習に匹敵する精度の物体認識ができた．
2. 画像モデルにViTを用いたとき，そのAttention mapは物体を識別し，segmentationのようなことができた．

![output_clusters]({attach}./images/dino_figs/output_clusters.jpg)
