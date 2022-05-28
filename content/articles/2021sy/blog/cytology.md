Title: 甲状腺細胞診分類モデル＠日本臨床細胞学会
Date: 2021.06.24
Modified: 2022.03.05
Tags: Machine Learning
Author: 安部

自分が大阪大学データビリティフロンティア機構の技術補佐員としてモデル開発に携わっている研究が、6/4~6/6に開催された第６２回日本臨床細胞学会総会にて発表がありました。
テーマは『**甲状腺細胞診画像から癌の種類を予測する**』でした。

この[抄録](https://jscc2021.jp/files/jscc2021_abstractse.pdf)のp112にある
「s3-3AI を用いた甲状腺細胞診支援システム(ADDICT)の開発に向けて」です。


## 初学会参加

旅費を援助して頂いたので幕張まで行って学会に参加してきました。様々な研究を聞いて新たな着想を得たり、質問をする過程で議論が深まったりと良い経験となりました。

## ドメイン知識

モデル作成にあたってCAMを見る必要があったのですが、4月にはCAMが妥当か全くわかりませんでした。しかし、かの有名な病理総論のテスト勉強のおかげで~~癌は~~という特徴があるから~~な細胞に注目していると妥当だろう！というような判断が少しできました。深層学習モデル作成に必要なドメイン知識が辛い？テスト勉強から得られたのは貴重な経験でした。


<img src="{attach}./images/cytology_figs/cytology_photo.jpg" alt="cytology_img" width="400px">

## 追記(2022/03/05)

この研究が以下の雑誌に掲載されました！

新岡 宏彦, 廣川 満良, 鈴木 彩菜,[**安部 政俊**]({author}安部), 新井 悠介, 式見 彰浩, 長原 一, 宮内 昭
“深層学習を用いた甲状腺細胞診自動診断システム (AI differential diagnosis for cytology of the thyroid:ADDICT)の開発とその現状”
PHARM TECH JAPAN, Vol. 38, No. 2, pp. 87 - 94.

廣川満良, 新岡宏彦, 鈴木彩菜, [**安部 政俊**]({author}安部), 式見彰浩, 長原一, 宮内昭
AIを用いた甲状腺細胞診支援システム(AI differential diagnosis for cytology of the thyroid:ADDICT)の開発と利用 
Journal of the Japanese Society of Clinical Cytology 61(3) 200-207 2022年6月