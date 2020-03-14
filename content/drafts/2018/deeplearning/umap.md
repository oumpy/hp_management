Title:超高速次元圧縮アルゴリズムUMAP
Date: 2019.04.14
Category: deeplearning
Tags: deeplearning
Slug: umap
Author: 安水
Summary:超高速次元圧縮アルゴリズムUMAP

<h1>超高速次元圧縮アルゴリズムUMAP</h1>

世間では、高性能な次元圧縮アルゴリズムとしてtSNEがよく使われています。 tSNEは便利ですが、少し遅いです。（パラメーターも意外に面倒。 <a href="https://deepage.net/machine_learning/2017/03/08/tsne.html">tsneのパラメータについて</a>）

そこで今回紹介するのはUMAP。 arxivで今月publishされたばかりのアルゴリズムです(記事執筆時 2018年2月時点)。 試しにMNIST(手書き数字画像。28*28=764次元)70000枚の次元圧縮をしてみました。

tSNEではちょうど1時間30分でしたが、UMAPではたったの1分でした。

詳細はこちら。 https://github.com/lmcinnes/umap

インストールは $ pip install umap

こちらが今回のソースコード。

<table>
<tbody>
<tr>
<td>from sklearn import datasets
import matplotlib.pyplot as plt
%matplotlib inline

digits = datasets.fetch_mldata('MNIST original')

print(digits['data'].shape)

import umap
embedding = umap.UMAP().fit_transform(digits.data)
plt.figure(figsize=(12,8))
plt.scatter(embedding[:,0],embedding[:,1],c=digits.target, s=0.1)
plt.title('UMAP')
plt.savefig('umap.png')

from sklearn.manifold import TSNE
model = TSNE(n_components=2)
tsne_result = model.fit_transform(digits.data)
plt.figure(figsize=(12,8))
plt.scatter(tsne_result[:,0],tsne_result[:,1],c=digits.target, s=0.1)
plt.title('tSNE')
plt.savefig('tsne.png')</td>
</tr>
</tbody>
</table>

<img title="" src="https://pythonoum.files.wordpress.com/2018/09/null1.png" alt="" width="334" height="222" /><img title="" src="https://pythonoum.files.wordpress.com/2018/09/null.png" alt="" width="333" height="221" />
