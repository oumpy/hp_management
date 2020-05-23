Title:超高速次元圧縮アルゴリズムUMAP
Date: 2018.09.15
Tags: Machine Learning
Author: 安水

## tSNEとUMAP
世間では、高性能な次元圧縮アルゴリズムとしてtSNEがよく使われています。 tSNEは便利ですが、少し遅いです。（パラメーターも意外に面倒。 tSNEのパラメータについては以下のリンク参照）

<https://deepage.net/machine_learning/2017/03/08/tsne.html>

そこで今回紹介するのは**UMAP**。 arxivで今月publishされたばかりのアルゴリズムです(記事執筆時 2018年2月時点)。 試しにMNIST(手書き数字画像。28*28=764次元)70000枚の次元圧縮をしてみました。

tSNEではちょうど1時間30分でしたが、UMAPではたったの1分でした。詳細は以下のリンクより。

<https://github.com/lmcinnes/umap>

## コード
インストールは `pip install umap`でok。

こちらが今回のソースコード。
```python
from sklearn import datasets
import matplotlib.pyplot as plt
%matplotlib inline
digits = datasets.fetch_mldata(‘MNIST original’)

print(digits[‘data’].shape)

import umap
embedding = umap.UMAP().fit_transform(digits.data)
plt.figure(figsize=(12,8))
plt.scatter(embedding[:,0],embedding[:,1],c=digits.target, s=0.1)
plt.title(‘UMAP’)
plt.savefig(‘umap.png’)

from sklearn.manifold import TSNE
model = TSNE(n_components=2)
tsne_result = model.fit_transform(digits.data)
plt.figure(figsize=(12,8))
plt.scatter(tsne_result[:,0],tsne_result[:,1],c=digits.target, s=0.1)
plt.title(‘tSNE’)
plt.savefig(‘tsne.png’)
```

## 結果
### tSNE
![tsne]({attach}images/umap_figs/tsne.png)

### UMAP
![umap]({attach}images/umap_figs/umap.png)
