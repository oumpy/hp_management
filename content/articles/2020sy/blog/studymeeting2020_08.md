Title: 【勉強会資料 2020 第8回】 視覚とSparse coding
Date: 2020.08.21
Tags: 勉強会, Neuroscience
Author: 山本

2020.08.21に行われた第8回勉強会のオンライン勉強会の資料を公開します。視覚生理学の基礎と第一次視覚野(V1)のモデルであるSparse codingモデル (Olshausen & Field, 1996) の説明を行いました。Sparse codingモデルをPythonで実装し、ガボールフィルタが生じることを確認しました 。また、単純細胞の受容野をよく再現するガボールフィルタの描画と画像への適応や、自然画像に対するICA・PCAについても紹介しました。

## 使用したNotebook
[Jupyter Notebook]({attach}./attach/studymeeting2020_08/第8回勉強会_sparse_coding.ipynb)

## 紹介スライド
[% embed src="{attach}./attach/studymeeting2020_08/第8回勉強会_視覚とSparsecoding.pdf" height="500" %]

## 参考資料
発表者が作成しているサイトに、Sparse codingモデルの生成モデルとしての解説を掲載しています (Julia実装)。

- [11.2 Sparse coding (Olshausen & Field, 1996) モデル — Juliaで学ぶ計算論的神経科学](https://compneuro-julia.github.io/11-2_sparse-coding.html)