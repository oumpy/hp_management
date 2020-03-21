Title:超高速次元圧縮アルゴリズムUMAP
Date: 2019.04.14
Tags: deeplearning
Slug: umap
Author: 安水
Summary:

世間では、高性能な次元圧縮アルゴリズムとしてtSNEがよく使われています。
tSNEは便利ですが、少し遅いです。（パラメーターも意外に面倒。https://deepage.net/machine_learning/2017/03/08/tsne.html）

そこで今回紹介するのはUMAP。
arxivで今月publishされたばかりのアルゴリズムです(2018年2月時点)。
試しにMNIST(手書き数字画像。28*28=764次元)70000枚の次元圧縮をしてみました。

tSNEではちょうど1時間30分でしたが、UMAPではたったの1分でした。

詳細はこちら。
https://github.com/lmcinnes/umap

インストールは
`$ pip install umap`


こちらが今回のソースコード。

```
from sklearn import datasets
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
plt.savefig('tsne.png')
```
