Title:Hamiltonian Descent Methodの要点
Date: 2019.04.14
Category:
Tags: deeplearning
Slug: hamiltonian_descent_method
Author: 山本
Summary:

解説自体は<a href="https://omedstu.jimdo.com/2018/09/26/hamiltonian-descent-methods%E3%81%AE%E5%AE%9F%E8%A3%85%E3%81%AB%E3%81%A4%E3%81%84%E3%81%A6%E3%81%AE%E8%A7%A3%E8%AA%AC/" target="_blank" rel="noopener">Hamiltonian Descent Methodsの実装についての解説</a>に書いたのですが、こちらには要点だけまとめておきます。

<h3>要点のまとめ</h3>

<img class="aligncenter size-full wp-image-170" src="https://pythonoum.files.wordpress.com/2018/09/1537941381.jpg" alt="1537941381" width="800" height="269" />
- Hamiltonian Descent法を用いると凸関数最適化が高速かつ高精度で行える。
- 超収束(Super-Convergence, <a href="https://arxiv.org/pdf/1708.07120.pdf">arxiv</a>)とは関係ある？論文読んだ限りはなさそう。
- ニューラルネットワークの最適化に応用するには運動エネルギー関数の研究が必要。
- なぜ高速に学習できるのか、なぜSGDと比べてニューラルネットワークの学習が不能なのか不明。でも学習できるならもっと実装あるよなと思っています。

<h3>アルゴリズム（1つ目の陽解法）</h3>

パラメータx以外に運動量のパラメータpを用意します。以下の更新式に従ってパラメータを更新します。
<img class="aligncenter size-medium wp-image-168" src="https://pythonoum.files.wordpress.com/2018/09/texclip20180929003048.png?w=600" alt="texclip20180929003048" width="300" height="81" />
ここでx_tは時刻tにおける解、p_tは時刻tにおける運動量、kは運動エネルギー関数、fは最小化したい関数（位置エネルギー関数）です2回微分が必要ですが、運動エネルギー関数k(p)を解析的に微分可能にしておけば、∇k(p)を定義すれば微分は1回でokです。
