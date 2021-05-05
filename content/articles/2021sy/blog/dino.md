Title: 【論文まとめ】DINO: Emerging Properties in Self-Supervised Vision Transformers
Date: 2021.05.05
Tags: Machine Learning, 論文まとめ
Author: 山本

本記事では主にFacebook AI Researchの研究者らによって提案された**DINO**という画像モデルにおける自己教師あり学習の解説を行います．

> Caron, Mathilde, Hugo Touvron, Ishan Misra, Hervé Jégou, Julien Mairal, Piotr Bojanowski, and Armand Joulin. 2021. “Emerging Properties in Self-Supervised Vision Transformers.” arXiv [cs.CV].  [http://arxiv.org/abs/2104.14294](http://arxiv.org/abs/2104.14294).

(cf.) [Facebook ブログ](https://ai.facebook.com/blog/dino-paws-computer-vision-with-self-supervised-transformers-and-10x-more-efficient-training), [GitHub](https://github.com/facebookresearch/dino), [Yannic Kilcherの解説動画](https://www.youtube.com/watch?v=h3ij3F3cPIk)

**要点**：画像モデルにおける，ラベル無しの自己教師あり学習の新しい手法DINOを考案したよ．画像モデルの最終層の出力は物体の種類ごとにクラスター化し，線形変換あるいはkNNを適応するだけで教師あり学習に匹敵する精度の物体認識ができたよ．特に画像モデルにVision Transformerを用いたとき，そのAttention mapは物体を識別し，segmentationのようなことができていたよ．

![fig1]({attach}./images/dino_figs/fig1.png)


## 背景
### Vision Transformers

Vision Transformer (ViT; ※ビタミンではない）

### Self-supervised learning
自己教師あり学習 (Self-supervised learning; SSL) 

言語モデルではBERT，ViTなら

**自己教師あり学習を用いた表現学習 (Self-supervised representation learning)** 手法は近年研究が盛んに行われています．**対照学習 (contrastive learning) **を用いた手法としては，MoCo ([He et al., 2019](https://arxiv.org/abs/1911.05722)), SimCLR ([Chen et al., 2020](https://arxiv.org/abs/2002.05709))などが代表的です．対照学習を用いた手法では，元画像とそれをaugmentした画像 (positive sample)，関係ない画像 (negative sample)の3つの画像を入力し，元画像とpositive sampleとの表現の類似度が小さく，negative sampleとの表現の類似度が大きくなるように学習を進めます．一方で，2020年後半頃より，**非対称的な2つのネットワークを用いる手法**が複数考案され，**BYOL** ([Grill et al., 2020](https://arxiv.org/abs/2006.07733)), SimSiam ([Chen & He, 2020](https://arxiv.org/abs/2011.10566)), Barlow Twins ([Zbontar et al., 2021](https://arxiv.org/abs/2103.03230)) などが該当します．これらの手法ではaugmentationを用いた非対称的な入力と，stop-gradientを用いた非対称的な重み更新による2つの非対称性を生み出し，2つのネットワークの出力の類似度を小さくするように学習を行います．この記事で紹介するDINOはBYOLに似た手法ですが，以下の3つが異なります．

1. BYOLとは異なる類似度損失を用いている．
2. 2つのネットワークが同じモデル構造をしている．
3. 手法をVision Transformerにも適応している．

## DINO : self-*di*stillation with *no* labels
### モデル構造
Caronらが提案するDINO (self-**di**stillation with **no** labels) とは「**ラベル無しでの自己蒸留**」を意味します．ここでの**蒸留 (distillation) **とは，枝刈り (pruning) や量子化 (quantization)に並ぶニューラルネットワークのモデル圧縮手法です．通常の蒸留では，ラベル付きデータセットとパラメータ数が多いモデル（**教師モデル; teacher model**）を用意し，ラベルと教師モデルの出力を教師信号 (hard & soft target) としてパラメータ数が少ないモデル（**生徒モデル; student model**）を訓練します ([Hinton, Vinyals & Dean, NIPS, 2014](https://arxiv.org/abs/1503.02531))．データセットだけでscratchから生徒モデルを訓練するより，教師モデルを用いた方が生徒モデルの性能は高くなるということが知られています．

DINOは教師モデルと生徒モデルを用いるところは通常の蒸留と同じですが，ラベル付きデータは用いず，**教師モデルは生徒モデルから**作られるという違いがあります．このことを念頭に置いて，モデルの構造を見てみましょう．

生徒モデルを$g _ {\theta s}$，教師モデルを $g _ {\theta t}$とします．入力画像を $\mathbf{x}$，画像 $\mathbf{x}$ を切り出す(crop)ようなaugmentationをした画像を$\mathbf{x}_1, \mathbf{x}_2$とします．globalな画像 (e.g. 元画像の50%以上)とlocalな画像 (e.g. 元画像の50%未満)を生成します．生徒モデルにはglobalとlocalの両方を入力するが，教師モデルにはglobalのみを入力します．こうすることで，“local-to-global”の対応が生成されます．


![model1]({attach}./images/dino_figs/model1.png)


モデルの崩壊を避けるために、モメンタムティーチャーの出力のセンタリングとシャープ化のみでも動作します。セクション5.3で実験的に示されているように、センタリングは1つの次元が支配的になるのを防ぎますが、一様分布への崩壊を促進します。両方の操作を適用することで、その効果がバランスされ、モメンタムティーチャーがある場合の崩壊を回避するのに十分な効果が得られます。崩壊を避けるためにこの方法を選択することは、安定性とバッチへの依存度の低さを交換することになります。センタリング操作は、1次のバッチ統計にのみ依存し、教師にバイアス項cを追加すると解釈できます：gt(x)←gt(x)+c。センタリングcは指数移動平均で更新され、


### パラメータの更新
![model2]({attach}./images/dino_figs/model2.png)


一方，教師モデルのパラメータ $\theta_t$ は生徒モデルのパラメータ $\theta_s$ を**指数移動平均 (exponential moving average; EMA)** することにより生成されます．
$$
\theta_t \leftarrow \lambda \theta_t + (1-\lambda) \theta_s
$$

### 実装

Pytorch-likeな擬似コードが掲載されています (Caron et al., 2021; Algorithm 1)．

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

