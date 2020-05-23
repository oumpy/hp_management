Title: 【論文まとめ】Improving Breast Cancer Detection using Symmetry Information with Deep Learning
Date: 2018.09.25
Tags: Machine Learning, 論文まとめ
Author: 秋山


> Improving Breast Cancer Detection using Symmetry Information with Deep Learning.
<http://arxiv.org/abs/1808.08273>

## 概要
- 放射線専門医がするようにマンモグラムの左右差を見て乳がんの判定をするCNNモデルを提案
- MICCAI (医療画像解析のトップ会議) 2018 採択
- 1st authorは放射線科医

## データセット
- マンモグラフィー7k症例, 正常4k, 悪性腫瘍3k
- 放射線専門医により腫瘍領域マスクをアノテーションした
- 正常乳房は2年フォローして乳がん発症が無かったことを確認
- 乳がん乳房は腫瘍領域を生検し悪性腫瘍であるとすべて診断されている
- 非公開データ

## 前処理
#### Fig. 1
![Fig. 1]({attach}images/breast_cancer_detection_figs/2352f582-60c8-486b-a9ec-c247787510f7.png)

- 従来法の画像特徴量を用いた手法で腫瘍候補点を計算 (Fig.1, 赤点)
- 候補点を中心に300x300 pixels (6cm x 6cm) を切り出し (Fig.1, 緑枠)
- 左右差を比較するために反対側の領域を切り出し (Fig.1, 青枠)
- 切り出した領域に腫瘍のマスク領域が含まれていれば悪性腫瘍のラベルを割り当てる

## モデル
- baseline model: 標準的なCNN (VGG like)
- symmetry model: 候補領域とその左右対称となる領域の2画像を入力とするCNN

#### Fig.2
![Fig. 2]({attach}images/breast_cancer_detection_figs/67ac70b1-0afe-478b-b046-a2bacb1450f8.png)

## 結果
![AUC table]({attach}images/breast_cancer_detection_figs/e382b9e382afe383aae383bce383b3e382b7e383a7e38383e38388-2018-09-25-17-28-58.png)

AUCがわずかに改善 (有意差なし)

#### Fig.3a
![Fig. 3a]({attach}images/breast_cancer_detection_figs/6aaac880-7ac2-4425-8be3-4d78367088e2.png)

FROC曲線(偽陽性率を横軸, 感受性を縦軸) で比較すると有意に提案手法がよかった

## 考察

#### Fig. 4
![Fig. 4]({attach}images/breast_cancer_detection_figs/ef07a0c8-ed1e-4127-b747-f9103aa6095e.png)

- a) baseline modelが正常と誤判定し, symmetry modelが正しく悪性腫瘍と判定した画像の例 (上下の画像が同じ患者の左右の乳房)
- b) baseline modelが悪性腫瘍と誤判定し, symmetry modelが正しく正常と判定した画像の例 (上下の画像が同じ患者の左右の乳房)

暗めの悪性腫瘍, 明るい正常像は誤判定しやすいがsymmetry modelは左右の比較によって正しく判定していることが読み取れる

## 読んだ感想
性能差がちょっと微妙だけど考察の納得感はある
