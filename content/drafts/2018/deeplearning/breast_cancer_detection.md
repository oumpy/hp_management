Title:Improving Breast Cancer Detection using Symmetry Information with Deep Learning
Date: 2019.04.14
Category:
Tags: deeplearning
Slug: breast_cancer_detection
Author: 秋山
Summary:Improving Breast Cancer Detection using Symmetry Information with Deep Learning

<a href="http://arxiv.org/abs/1808.08273">Improving Breast Cancer Detection using Symmetry Information with Deep Learning</a>
放射線専門医がするようにマンモグラムの左右差を見て乳がんの判定をするCNNモデルを提案
arxiv 2018.08
MICCAI (医療画像解析のトップ会議) 2018 採択
1st authorは放射線科医

データセット
マンモグラフィー7k症例, 正常4k, 悪性腫瘍3k
放射線専門医により腫瘍領域マスクをアノテーションした
正常乳房は2年フォローして乳がん発症が無かったことを確認
乳がん乳房は腫瘍領域を生検し悪性腫瘍であるとすべて診断されている
非公開データ

前処理
従来法の画像特徴量を用いた手法で腫瘍候補点を計算 (Fig.1, 赤点)
候補点を中心に300x300 pixels (6cm x 6cm) を切り出し (Fig.1, 緑枠)
左右差を比較するために反対側の領域を切り出し (Fig.1, 青枠)
切り出した領域に腫瘍のマスク領域が含まれていれば悪性腫瘍のラベルを割り当てる

<img class="alignnone size-full wp-image-153" src="https://pythonoum.files.wordpress.com/2018/09/2352f582-60c8-486b-a9ec-c247787510f7.png" alt="2352f582-60c8-486b-a9ec-c247787510f7" width="857" height="567" />

Fig. 1

モデル
baseline model: 標準的なCNN (VGG like)
symmetry model: 候補領域とその左右対称となる領域の2画像を入力とするCNN

<img class="alignnone size-full wp-image-152" src="https://pythonoum.files.wordpress.com/2018/09/67ac70b1-0afe-478b-b046-a2bacb1450f8.png" alt="67ac70b1-0afe-478b-b046-a2bacb1450f8" width="1487" height="470" />

結果

<img class="alignnone size-full wp-image-155" src="https://pythonoum.files.wordpress.com/2018/09/e382b9e382afe383aae383bce383b3e382b7e383a7e38383e38388-2018-09-25-17-28-58.png" alt="スクリーンショット 2018-09-25 17.28.58" width="1516" height="296" />
AUCがわずかに改善 (有意差なし)

<img class="alignnone size-full wp-image-151" src="https://pythonoum.files.wordpress.com/2018/09/6aaac880-7ac2-4425-8be3-4d78367088e2.png" alt="6aaac880-7ac2-4425-8be3-4d78367088e2" width="671" height="484" />

Fig.3a
FROC曲線(偽陽性率を横軸, 感受性を縦軸) で比較すると有意に提案手法がよかった

考察

<img class="alignnone size-full wp-image-154" src="https://pythonoum.files.wordpress.com/2018/09/ef07a0c8-ed1e-4127-b747-f9103aa6095e.png" alt="ef07a0c8-ed1e-4127-b747-f9103aa6095e" width="960" height="360" />

Fig. 4
a) baseline modelが正常と誤判定し, symmetry modelが正しく悪性腫瘍と判定した画像の例 (上下の画像が同じ患者の左右の乳房)
b) baseline modelが悪性腫瘍と誤判定し, symmetry modelが正しく正常と判定した画像の例 (上下の画像が同じ患者の左右の乳房)
暗めの悪性腫瘍, 明るい正常像は誤判定しやすいがsymmetry modelは左右の比較によって正しく判定していることが読み取れる

読んだ感想
性能差がちょっと微妙だけど考察の納得感はある

&nbsp;
