Title: 自転車コンペ
Date: 2021.11.08
Modified: 2021.11.08
Tags: Data Science Competition,Machine Learning
Author: 梅津

# シェアサイクルコンペ
2021/10/01~2021/10/27にsignateで開催された[SIGNATE Student Cup 2021秋【予測部門】](https://signate.jp/competitions/550)で377人中26位になり銀メダルを獲得しました。
## コンペの概要
予測する日以前のデータを使用し、予測当日での70ヶ所の駐輪場での利用可能な自転車数を予測する回帰問題。時系列データで1時間ごとの値が与えられた。<br>
予測する日は複数あり、それらは全てが連続しているわけではなく前後に訓練に回せるデータがある場合もある。しかし、上記に書いた通り予測する日に対し翌日からのデータを使用するのは禁止。リークに注意する必要があった。<br>
利用可能な自転車数から気候など様々な特徴量が4つのテーブルから与えられた。<br>
評価指標はRMSE(root mean squared error)
## したこと
今回は単独での参加だったのでEDAからモデリングまで自分でやりました。<br>
1. feature engineering
平日、祝日、休日や春夏秋冬の季節を作成しました。
2. train&validation
予測する日以前のデータの全てをtrainにまわしました。欠損したところはその直前の予測日以前のデータでtrainしたモデルでの予測値を代入し、pseudo labelingをしました。120個のlgbmをtrainすることになりましたが、結局は月毎で12個つくった方が精度が出るみたいです。validationを上記の訓練する状況に似せようと1日だけにしたのが悪手なのかも...
1つ目のlgbmでoptunaを使いハイパーパラメータを最適化し、他にも流用しました。
seed値を変えた３つの120個のlgbmのセットをensambleに使いました。
## 他にやりたかったこと
LSTMやTabnetなどのNNも試してみたかった。<br>
ラグ特徴量のような時系列データで効果的な特徴量も作りたかった。<br>
RFのハイパーパラメータを最適化し使おうとしたが何時間たっても終わる気配がせず断念。最適化しなかった場合、lgbmとensambleするとlgbmのRandom Seed Averageよりも精度が低かった。<br>
station毎にモデルを訓練させたかった。
## 反省点
多くありますが、最初にEDAがあまく、予測フラグが立っていない行の目的変数がNaNであることがありえることに気付かず、精度が出ませんでした。EDAは大切！<br>
train用のデータの説明変数に目的変数を入れたままなのに気付かず、お手上げ状況になった。GBDTへの理解が足らず、リークした際の状況を推測できなかったのが原因。<br>
カテゴリ変数のencodingを一切関数を使わずに行ったため、次からは使いたい。
## とても参考になったもの
[Kaggleで勝つデータ分析の技術](https://www.amazon.co.jp/Kaggle%E3%81%A7%E5%8B%9D%E3%81%A4%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E3%81%AE%E6%8A%80%E8%A1%93-%E9%96%80%E8%84%87-%E5%A4%A7%E8%BC%94/dp/4297108437/ref=sr_1_1?adgrpid=109872001568&dchild=1&gclid=Cj0KCQjw5oiMBhDtARIsAJi0qk2fQg-LTS200pIygImVOdjg50nrlrsEpo_AZNS3GSV43BlV3ZLwrJgaArUyEALw_wcB&hvadid=553921514572&hvdev=c&hvlocphy=1009543&hvnetw=g&hvqmt=e&hvrand=11529920702029821113&hvtargid=kwd-818170082585&hydadcr=27492_14478797&jp-ad-ap=0&keywords=kaggle%E3%81%A7%E5%8B%9D%E3%81%A4%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E3%81%AE%E6%8A%80%E8%A1%93&qid=1635941789&sr=8-1)<br>
[How to Win a Data Science Competition: Learn from Top Kagglers](https://www.coursera.org/learn/competitive-data-science)<br>
使ったコードは[こちら](https://github.com/yumezu121/signate_bikes/tree/main) です。