Title: 行列演算の初歩の初歩とNumPy
Date: 2019.05.16
Modified: 2021.08.28
Tags: Python
Author: 小川

機械学習をはじめ、何をやるにも行列演算は必須です。
本稿ではその初歩と、ライブラリ [NumPy](https://numpy.org/) での基本的な実行方法を簡単に解説します。
(主に1年生向け。)

> 引用表記の部分

はやや発展的な内容なので、初めて読むときは飛ばしてOKです。

## ベクトルの積
### 定義
まず、**ベクトルの積**を以下のように定義します。  
(Python流に0-index記法を採用。)

$$
x=(x_0,x_1,...,x_{D-1}),
\quad
y=\left(\begin{array} 
\@
y_0 \\
y_1 \\
\vdots\\
y_{D-1} \\
\end{array}\right)
\quad\Rightarrow\quad
x\cdot y = x_0y_0+x_1y_1+\dots+x_{D-1}y_{D-1}
$$
高校で学ぶベクトルの内積と基本的に同じですが、**行ベクトル (横ベクトル) と列ベクトル (縦ベクトル) を区別**するのが大事なところです。
掛け算は、同じ次元の**「行ベクトル $\times$ 列ベクトル」の順序**でしか行えません。  
(本当は逆順でもできるが、結果は全く違うものになる。
どうなるか？最後まで読んだら考えてみましょう。)

> 行ベクトルと列ベクトルはそれぞれ互いに双対なベクトル空間の元 (の数ベクトル表現) です。
> 列ベクトル同士・行ベクトル同士で内積が定義できるのは、列ベクトルと行ベクトルを互いに変換する写像 (計量) が定義された**計量ベクトル空間**に限られます。  
> 高校で習うベクトル空間に内積がしれっと入っているのは、空間上に自然なユークリッド計量が暗黙に仮定されているからです。


## 行列：ベクトルを並べたもの
行列は、このように成分を2次元状 ($E\times D$) に並べたもの。
$$
X = \left(
\begin{array}\@
x_{0,0}   & x_{0,1} & \dots  & x_{0,D-1} \\
x_{1,0}   & x_{1,1} &        & \vdots \\
\vdots    &         & \ddots & \\
x_{E-1,0} & \dots   &        & x_{E-1,D-1}
\end{array}
\right)
$$
これは  
(1)「$D$次元行ベクトルを縦に$E$個並べたもの」  
または  
(2)「$E$次元列ベクトルを横に$D$個並べたもの」  
です。  
**どちらの解釈をしても構いません**。

### (1) 列ベクトルとの掛け算
前者(1)の解釈をすると、上の列ベクトル$y$と掛け算ができます。
$$
X\cdot y =
\left(\begin{array}
\@
x_{0,0}y_0+x_{0,1}y_1+\dots+x_{0,D-1}y_{D-1} \\
x_{1,0}y_0+x_{1,1}y_1+\dots+x_{1,D-1}y_{D-1} \\
\vdots \\
x_{E-1,0}y_0+x_{E-1,1}y_1+\dots+x_{E-1,D-1}y_{D-1}
\end{array}\right)
$$
結果は**E次元列ベクトル**になりました。

**列ベクトルに (正しい次元の) 行列を左から掛ければ列ベクトルに**なります。

### (2) 行ベクトルとの掛け算
一方で解釈(2) をとると、$E$ 次元行ベクトル $z=(z_0,z_1,\dots,z_{E-1})$ に対して、**右からの掛け算**  $z\cdot X$ を行うことができます。
結果は **D次元行ベクトル** になります。

つまり、**行ベクトルに（正しい次元の）行列を右から掛ければ行ベクトルになる**ことがわかります。

> ### 掛け算の結合則
行列に対して、(1)(2) どちらの解釈をしても構わない、と言いました。
この意味は例えば
$$
z\cdot X\cdot y
$$
という式が与えられた際に現れます。  
(1) の解釈をするなら、この式は$z\cdot(X\cdot y)$と計算、(2) ならば$(z\cdot X)\cdot y$です。
どちらもただ一つの数を与えますが、実は両者は**常に同じ**になることが、実際に計算すれば簡単に証明できます。  
これは後述の行列積にも拡張でき、行列積は全て順序によりません。
ただし**交換則は成立しない** ($X\cdot Y\ne Y\cdot X$) ので注意。

## 行列同士の掛け算
ここまでくると、掛け算をさらに拡張できる。
$D\times F$ 行列 $Y$ を
$$
Y = 
\left(\begin{array}\@
y_{0,0}           & y_{0,1} & \dots  & y_{0,F-1} \\
y_{1,0}           & y_{1,1} &        & \vdots \\
\vdots            &         & \ddots & \\
y_{D-1} & \dots   &        & y_{D-1,F-1}
\end{array}\right)
$$
とするとき、行列同士の掛け算は
$$
X\cdot Y = W = 
\left(\begin{array}\@
w_{0,0} & \dots & w_{0,F-1} \\
\vdots  & \ddots & \vdots \\
w_{E-1,0} & \dots & w_{E-1,F-1}
\end{array}\right)
$$

$$
w_{i,j} = x_{i,0}y_{0,j}+x_{i,1}y_{1,j}+\dots+x_{i,D-1}y_{D-1,j}
$$
となる。
左の $X$ から行ベクトル、右の $Y$ から列ベクトルを取り出して、それぞれ掛け算したものを順番に並べる。  

$E\times D$ と $D\times F$ を掛けると真ん中の $D$ が潰れて $E\times F$ 行列、と覚えればいいです。
もちろん $D$ が揃っていないと掛け算できない。

## ベクトル、行列、掛け算の意味
- ベクトルは通常、「データ」とか「状態」などの、比較的直接的な意味を持ちます。
特に数式として書く場合、普通は列ベクトルにその役割を持たせます。
- データや状態がたくさんある場合、それら列ベクトルを横に並べれば行列になり、その行列がデータ／状態の集合を表します。
別の行列を左から掛ければ、各列に独立・一括で演算ができます。
- 行列は基本的に、上で見た通り、**列ベクトルに左から作用して別の列ベクトルを作る演算子**です。
行ベクトルも、列ベクトルに作用する行列のうち、行が一つしかない特別な場合と見なせます。
- 例えば、画像を表すベクトルに行列を掛け算することで、**画像の持つ特徴 (線の向きとか) を表す別のベクトルを得る**、といったことなどができます。
(CNN = 畳み込みニューラルネットワーク)
- 行列の掛け算は合成関数のようなもの。掛け算の結合則があるので、ベクトル$x$に順番に変換を施す操作 $A\cdot(B\cdot x)$ は、$=(A\cdot B)\cdot x$ のように予め左側を先に計算しておけます。
これにより、$x$ を様々に取り替えて計算する場合に計算量を低減できます。

## NumPyによる行列演算
Python の NumPy パッケージは強力な行列演算機能を備えています。
行列の定義は `np.array()`、掛け算は `np.dot()` または `@`演算子 (Python3.5以降) で行えます。

```python
import numpy as np

A = np.array([[1,2,3],[4,5,6]]) # 2行3列(2x3)行列
x = np.array([[1],[2],[3]])     # 3次元列ベクトル
print(A @ x)
```
これで、

```python
array([[14],
       [32]])
```
と表示されれば成功。結果は2次元列ベクトル。

他にも色々試してみてください。

なお若干ややこしいことに、NumPy にはいくつも便利機能があり、数学的には不正な式でも適当に補完してやってくれたりします。
(例えば行列とベクトルの足し算などもできる。)  
ただ慣れるまでは、ひとまず数学的に忠実な記法を心掛けた方がいいでしょう。

とりあえず以上。
おしまい。
