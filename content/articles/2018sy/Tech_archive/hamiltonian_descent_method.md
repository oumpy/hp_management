Title:Hamiltonian Descent Methodの要点
Date: 2019.04.14
Tags: deeplearning
Slug: hamiltonian_descent_method
Author: 山本
Summary:

解説自体は[Hamiltonian Descent Methodsの実装についての解説]("https://omedstu.jimdo.com/2018/09/26/hamiltonian-descent-methods%E3%81%AE%E5%AE%9F%E8%A3%85%E3%81%AB%E3%81%A4%E3%81%84%E3%81%A6%E3%81%AE%E8%A7%A3%E8%AA%AC/")に書いたのですが、こちらには要点だけまとめておきます。

## 要点のまとめ

![1537941381]({attach}images/hamiltonian_descent_method_figs/1537941381.jpg)
- Hamiltonian Descent法を用いると凸関数最適化が高速かつ高精度で行える。
- 超収束(Super-Convergence, [arxiv](https://arxiv.org/pdf/1708.07120.pdf))とは関係ある？論文読んだ限りはなさそう。
- ニューラルネットワークの最適化に応用するには運動エネルギー関数の研究が必要。
- なぜ高速に学習できるのか、なぜSGDと比べてニューラルネットワークの学習が不能なのか不明。でも学習できるならもっと実装あるよなと思っています。

## アルゴリズム（1つ目の陽解法）

パラメータ$x$以外に運動量のパラメータ$p$を用意します。以下の更新式に従ってパラメータを更新します。
$$
p_{i+1} = \delta p_i - \epsilon\delta\nabla f(x_i)\\
x_{i+1} = \delta x_i + \epsilon\nabla k(p_{i+1})
$$
ここで$x_t$は時刻$t$における解、$p_t$は時刻tにおける運動量、$k$は運動エネルギー関数、$f$は最小化したい関数（位置エネルギー関数）です。2回微分が必要ですが、運動エネルギー関数$k(p)$を解析的に微分可能にしておけば、$\nabla k(p)$を定義すれば微分は1回でokです。
