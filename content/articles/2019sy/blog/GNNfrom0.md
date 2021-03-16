Title: Graph Neural Network (GNN)をゼロから実装した話
Date: 2019.05.13
Modified: 2021.03.15
Tags: Machine Learning
Author: 小川

医学科3年の小川です。
最近どうも[某ネトゲ](https://atcoder.jp)の話しかしていない気がするので、たまには少し違う話題。

## 某社の試験課題
その筋で有名な某社の某イベント[[注1](#note1)]に応募したところ、選抜試験として次の課題が送られてきました。  

> **Graph Neural Network (GNN)を言語標準ライブラリのみで実装せよ。**
> ただし行列演算ライブラリ（Pythonならnumpy、C++ならEigenなど）だけは用いてもよい。

つまりTensorflowとか使ってはいけない。
広義のフルスクラッチ実装。

本記事はそれでやってみた（GWに丸2日くらいかけて実装・実験した）ことの記録です。
主要部分を紹介解説しますが、コード全体を見てみたい方は[こちら](https://github.com/pybach/2019GNNfrom0/)に置いてあります。

（注：その方面では全然高度な話じゃないけれど、**初心者には厳しめの内容**です。）

## 参考：オコゼ本
<img src="https://www.oreilly.co.jp/books/images/picture_large978-4-87311-758-4.jpeg" width=18%>

ニューラルネットワークのフルスクラッチ実装は珍しいネタではなく、[「ゼロから作るDeep Learning 〜Pythonで学ぶディープラーニングの理論と実装〜」（斎藤康毅著、O'Reilly Japan）](https://www.oreilly.co.jp/books/9784873117584/)は街の書店でも良く見掛けます。  

実はこの本、全然読んだことがなかったのですが、後から見ると僕が書いたものと**多くのコードがそっくり**で、**これくらいの安直な実装は誰が書いても大差ない**ことを示唆しているのかもしれない。
誰が書いてもそっくり、は**Pythonの売り**でもあります。  
（なお本記事で解説する**メモリ共有**はこの本ではやっていない。
また**クラスの継承**も使っていない。
どちらも概念的にやや複雑になるので避けた、のかも。）

ともあれ、この本は読みやすいうえ、扱っている内容は当然本記事より広く、基本的な事項をかなりカバーしていてお勧めできます。  
（一応表面上は、知識ゼロでも読めるように書かれている。
とはいえ、パラパラ見て知らない単語だらけの場合は、まだやめておいた方が無難かも。）


## Graph Neural Network (GNN)
画像データなどではなく、「ネットワークのデータ」を分類したりするのに用いられるニューラルネットワーク。
例えば化学物質の構造（どの原子と原子が結合している／していないか）を与えて物性や生理活性を予測するとか、そういうことをイメージしてもらえばよいです。
この例なら**創薬に使う**ことなどが期待できます。

最近流行していて、グラフ畳み込みネット(GCN)とか様々なバージョンが研究されていますが、ここでは最も単純なものを考えます。  

無向重み無しグラフ$G=(\mathrm{nodes},\mathrm{edges})$が与えられたとき、各頂点$i\in\mathrm{nodes}$に付随する$D$次元の初期状態ベクトル$x_i^{(0)}$を適当に与える（整数$D$はGNNの固定パラメータ）。
そこから以下の操作（**集約**=aggregation）を定められた回数（$T$回）反復する。

- **集約1**：$a_i^{(t)} = \sum_{j\; s.t.\, (i,j)\in\mathrm{edges}}x_j^{(t)}$
- **集約2**：$x_i^{(t+1)} = f(W\cdot a_i^{(t)})$

ここで$W$はGNNの持つ$D\times D$行列値の学習パラメータ。
$f(\cdot)$は成分ごとに適用される適当な活性化関数(activation function)。
今回はReLU関数（x>0でxを、x<=0で0を返す関数）とする。

集約を終えた後、全頂点ベクトルの総和として、グラフ$G$の特徴ベクトル$h_G$を返す（**読み出し**=readout）。

- **読み出し**：$h_G = \sum_{i\in\mathrm{nodes}} x_i^{(T)}$

最後、この$h_G$をロジスティック／ソフトマックス層などに突っ込んで、分類問題などを解くわけです。
ニューラルネットワークとしては再帰ニューラルネット(RNN)の特別な場合と見ることができます。

## 問題
実際に出された問題（要約）です。
4段階に分かれていて、順番に機能を追加実装していけ、という感じ。

- **課題1**：行列$W$を固定として、与えられたグラフ$G$に対して$h_G$を求めるGNN関数またはクラスを実装してテストせよ。
- **課題2**：グラフに対する二値分類問題（$y=0, 1$）を解くため、陽性確率を$p=\mathrm{Sigmoid}(s)$、$s=A\cdot h_G+b$、学習パラメータセットを$\theta=\{W,A,b\}$として（$A$は$D$次元ベクトル、$b$はスカラー）、単一データ$(G,y)$に対する勾配降下法を実装してテストせよ。  
（勾配は数値微分、損失関数はbinary crossentropyとする。)
- **課題3**：データセットに対する学習が行えるよう、確率的勾配降下法(SGD)およびmomentum SGDによるミニバッチ学習を実装し、与えられたN=2000の実データセットで性能評価せよ。
- **課題4**：勾配降下をさらにAdamにするとか、$W$を多層化するとか好きにやってみて、最後にラベルのないテストデータに対する予言値を提出せよ。

使用言語は7言語（C, C++, Python, Ruby, Go, Java, Rust）から選択。
**やっぱりPython**だよね、うん。  
（まあ実用的なプロジェクトならC++でしょうが、そもそもフルスクラッチとかしないし、Tensorflowのようなフレームワークを使う前提なら、余程でない限り**結局Python**になります。）

## 実装
**クラスとして実装**しました。というか他のやり方が事実上思いつかない。
（なので、言語選択肢の中でもCは勘弁願いたいです。）  
課題が進む毎に、**前課題のクラスを継承して機能追加**した**派生クラス**を作っていきます。

### 課題1 : 集約と読み出し
> **課題1**：行列$W$を固定として、与えられたグラフ$G$に対して$h_G$を求めるGNN関数またはクラスを実装してテストせよ。

`GNN1`クラスとして実装。
Gを隣接行列で与えることにすると、集約は2段階とも単なる行列積（と活性化関数）です。numpyの**np.dot()**で簡単・高速に行えます。

```python
# task1.py
import numpy as np

class GNN1:
    def __init__(self, D, T, W):
        ...
    ...
    @staticmethod
    def f(x):
        return np.frompyfunc(relu, 1, 1)(x)
        
    @staticmethod
    def aggregate1(G, x):
        a = np.dot(G, x)
        return a
    def aggregate2(self, a):
        x = self.f(np.dot(a, self.W))
        return x
        
    def readout(self, G):
        N = len(G)
        x = self.x0(N)
        for _ in range(self.T):
            a = self.aggregate1(G, x)
            x = self.aggregate2(a)
        return x.sum(axis=0)
```

### 課題2 : 勾配降下
> **課題2**：グラフに対する二値分類問題 $(y=0, 1)$ を解くため、陽性確率を $p=\mathrm{Sigmoid}(s)$、$s=A\cdot h_G+b$、学習パラメータセットを $\theta=\{W,A,b\}$ として（ $A$は$D$次元ベクトル、$b$はスカラー）、単一データ $(G,y)$ に対する勾配降下法を実装してテストせよ。  
（勾配は数値微分、損失関数はbinary crossentropyとする。)

前課題の`GNN1`クラスを継承し、**損失関数や勾配を計算する機能を追加**した派生クラス`GNN2`を作成します。

実装のポイントは初期化時に** `Theta` と `(W,A,b)` でメモリを共有**するようにしておくこと。
今の場合、`Theta` は(D^2 +D+1)次元ベクトルであると**同時に、**そのうち最初の$D^2$成分は$DxD$行列 `W`、次のD成分はD次元ベクトル `A`、最後の成分は `b`、としてもアクセスできるようになっています。

こうすることで、数値微分や勾配降下の際に**学習パラメータの構造を何も考えなくてよい**という利点があります。
実際、`GNN2.grad()` や `GNN2.shift_Theta()` 関数の中には `W` や `A` の文字すら登場しません。  
（そのため、クラスを継承・拡張してパラメータ構造に変化があっても、修正する必要がない。）  

```python
# task2.py
from task1 import *

class GNN2(GNN1):
    def __init__(self, D, T, sigma=0.4,
                 epsilon=1.0e-3): # 数値微分の微小変分
        ...
    ...
    
    # 学習パラメータの初期化。Thetaと(W,A,b)でメモリを共有する。
    # bも1要素numpy配列にしておくこと。
    def init_Theta(self, sigma):
        self.Theta = np.random.normal(0,sigma,self.D*self.D+self.D +1)
        self.W, self.A, self.b = self.decode_Theta()
    def decode_Theta(self):
        return (self.Theta[:self.D*self.D].reshape((self.D,self.D)),
                self.Theta[self.D*self.D:self.D*(self.D+1)],
                self.Theta[self.D*(self.D+1):self.D*(self.D+1)+1])
        # numpy配列のスライスやreshapeはもとのオブジェクトとメモリを共有している。
        # （メモリのコピーを作成しない。）
        # 不意の副作用を気にすることが多いが、ここでは積極的に利用。

    # GNN1で定義したreadout()を用いて陽性確率／損失関数を計算する
    def s(self, G):
        h = self.readout(G)
        return np.dot(h,self.A) + self.b[0]
    def p(self, G):
        return sigmoid(self.s(G))
    def loss(self,G,y):
        return binary_cross_entropy(y, self.s(G))

    # 損失関数loss()の数値微分
    def grad(self, G, y):
        L0 = self.loss(G,y)
        K = len(self.Theta)
        dL = np.empty(K)
        for i in range(K):
            self.Theta[i] += self.epsilon
            dL[i] = (self.loss(G,y) - L0)/self.epsilon
            self.Theta[i] -= self.epsilon
        return dL

    # 勾配降下
    def shift_Theta(self, dTheta):
        self.Theta += dTheta
    def descendant_update(self, G, y, alpha):
        self.shift_Theta(-alpha*self.grad(G,y))
```
 
### 課題3：データセットからのミニバッチ学習（SGD, momentum SGD）
> **課題3**：データセットに対する学習が行えるよう、確率的勾配降下法(SGD)およびmomentum SGDによるミニバッチ学習を実装し、与えられたN=2000の実データセットで性能評価せよ。

ここでは (momentum) SGDの実装を求められているわけですが、アルゴリズム個々に学習メソッドを実装するのはさすがに無駄が多過ぎます。  
実際は学習メソッドは `GNN3.fit()` ひとつだけでいいです。

```python
# task3.py
from task2 import *

class GNN3(GNN2):
    ...
    
    def fit(self,
            x, y,   # 学習データ。
            optimizer,      # Optimizer派生クラスのオブジェクトを与える。
            epochs=1, batchsize=16,
            validation=None     # 検証データ。タプル (vx,vy) で与える。
            ):
        n = len(y)
        ...
        
        for e in range(epochs):
            # epoch毎にデータをランダムに並べ替える。
            order = np.random.permutation(n)

            for b in range((n-1)//batchsize+1):
                start = batchsize * b
                end = min(start+batchsize, n)
                # 勾配ミニバッチ平均
                grad = np.mean([ self.grad(x[i],y[i])
                                 for i in order[start:end]],axis=0)
                # optimizerで変化量を求め、Thetaを動かす。
                self.shift_Theta(optimizer.update(grad))
                ...
        ...
        # loss, accの履歴を返す。
        return batch_losses, epoch_losses, epoch_accs, batch_vlosses, epoch_vlosses, epoch_vaccs

```

ポイントは `optimizer` 引数。
ここに `SGD()` やら `MomentumSGD()` やらを、`gnn3.fit(x,y,optimizer=SGD())` のように与えて処理を切り替えます。  
（普段ニューラルネットワーク系のライブラリを使っていると、このあたりの設計は当然ですが。）

また `order` を順に見ていくことでエポック内非復元抽出を行います（さほど重要というほどではありませんが。オコゼ本では復元抽出）。

オプティマイザの中身はこちら。

```python
# task3.py つづき
### ほぼダミーの基底クラス。（抽象クラスにした方がすっきりするが、さぼった）
class Optimizer():
    @staticmethod
    def update(grad):
        # 勾配を受け取り、内部状態を更新し、変化量を返すメソッド。
        # （派生クラスではそのような関数に再定義する。）
        return -grad

### 確率的勾配降下法
# Optimizerを継承してupdate()を上書き再定義。
class SGD(Optimizer):
    def __init__(self, alpha=1.0e-4):
        self.alpha = alpha

    def update(self, grad):
        # SGDは更新される内部状態を持たない。
        return -self.alpha * grad

### 運動量付き確率的勾配降下法
class MomentumSGD(Optimizer):
    def __init__(self, alpha=1.0e-4, eta=0.9):
        self.alpha = alpha
        self.eta = eta
        self.w = None

    def update(self, grad):
        # 内部のwを更新しながら変化量を返す。
        if self.w is None:
            self.w = np.zeros(len(grad))
        self.w *= self.eta
        self.w -= self.alpha * grad
        return self.w
```
これだけ。
処理が分離されているので、わかりやすいですよね？

### 課題4 : その他の工夫と未知データの推論
> **課題4**：勾配降下をさらにAdamにするとか、$W$を多層化するとか好きにやってみて、最後にラベルのないテストデータに対する予言値を提出せよ。

#### (a) Adamの実装
オプティマイザ派生クラスを作るだけ。GNNクラスには1ミリも触れなくてOK。

```python
# task4a.py
from task3 import *

# 原論文 https://arxiv.org/abs/1412.6980 の通り、愚直に実装する。
class Adam(Optimizer):
    def __init__(self,
                 alpha=1.0e-3, beta1=0.9, beta2=0.999,
                 epsilon=1.0e-8):
        # 内部状態ベクトル m,v を初期化
        ...
        
    def update(self, grad):
        # m,v を更新して変化量を返す
        ...
```
#### (b) 多層ニューラルネット化
祖先クラスである `GNN1` や `GNN2` で定義した関数を、`W` の多層化に合わせて上書き再定義した派生クラス `GNN4` を作る。
変更箇所はわりと少なく、実質これだけ。

```python
# task4b.py
from task3 import *

class GNN4(GNN3):
    def __init__(self, D, T, sigma=0.4,
                 Nw=2): # Wのlayer数。 W.shape=(Nw,D,D)
        ...

    ### パラメータ初期化の関数2つを上書きする。
    def init_Theta(self, sigma):
        self.Theta = np.random.normal(0,sigma,Nw*self.D*self.D+self.D+1)
        self.W, self.A, self.b = self.decode_Theta()

    def decode_Theta(self):
        Wsize, Asize = self.Nw*self.D*self.D, self.D
        return (self.Theta[:Wsize].reshape((self.Nw,self.D,self.D)),
                self.Theta[Wsize:Wsize+Asize],
                self.Theta[Wsize+Asize:Wsize+Asize+1])

    ### 集約2を多層ニューラルネットに変更。
    def aggregate2(self, a):
        x = a
        for i in range(self.Nw):
            x = self.f(np.dot(x, self.W[i]))
        return x
```

## 結果
### SGD, momentum SGD, Adamの学習曲線（単層ネットワーク）

与えられたN=2000のデータセットで、(学習,検証)=(1600,400) に分割してとりあえず描いてみた学習曲線。
とにかくAdamの性能が圧倒的にいいですね。
このグラフは1-shotだけど、反復しても傾向は同じ。

<a href="{attach}./images/GNNfrom0_figs/task4a_plot03.pdf"><img src="{attach}./images/GNNfrom0_figs/task4a_plot03.png"></a>
<a href="{attach}./images/GNNfrom0_figs/task4a_plot03a.pdf"><img src="{attach}./images/GNNfrom0_figs/task4a_plot03a.png"></a>

### 多層化GNNの性能評価
- 単層`GNN3`と`Nw=2`の`GNN4`で比較（Adamを使用）。
2層化で成績が向上しているとは言えないっぽい。。
どちらも正答率(vacc)は60%程度。
（図のaccはvaccの意味。）

<a href="{attach}./images/GNNfrom0_figs/task4b_plot04.pdf"><img src="{attach}./images/GNNfrom0_figs/task4b_plot04.png"></a>
<a href="{attach}./images/GNNfrom0_figs/task4b_plot04a.pdf"><img src="{attach}./images/GNNfrom0_figs/task4b_plot04a.png"></a>

- 補足と反省
	- 後で試したところ、実は `Nw=6` くらいまで深層化すると **vacc>=0.65** くらいにはなることがわかった。
    やはり多層化は無駄ではないらしい。
   - ちゃんと評価するなら、データセットを (train,valid,test) に3分割し、validのスコアが高いパラメータセットを選んでtestを評価、という手順を踏むべきである。
   今回はそこまでやっていない。
   - さらに本気で学習＆パラメータ探索をするなら、cupy版を作ってGPUを使う方がいい。

## まとめ
GNNを一通り、少なくとも明白なバグなく、ミニバッチ学習までちゃんと動かすことができました。
時間が十分確保できず、最後のネットワーク改良を全然詰めることができなかった点は残念。
でもフルスクラッチで手書きした学習がちゃんと動いてそれなりに非自明な結果を出す、というのはなかなか楽しく、一種の全能感（？）を味わうことができます。

またGNNって殆ど勉強したことがなかったのですが、実装・実験してみて感覚は何となくわかりました（小並感）。
折角なので拡張や理論をもう少しちゃんと勉強して、実際の研究にも使ってみたいと思います。

おしまい。

## 注記

### <span id="note1"></span>注1 : 某社の某イベント
[Preferred Networks社](https://www.preferred.jp) 2019年度夏季インターン、でした。

本稿に掲載したコードほぼそのままで提出。
後日談となりますが、面接を経て採用して頂きました。

インターン期間は8-9月で、医学科の夏休みは8月だけなのですが、ちょうど[基礎配属](http://www.edu.med.osaka-u.ac.jp/assignment/)期間の一部を利用してフル参加することができた、という事情もあったりします。
