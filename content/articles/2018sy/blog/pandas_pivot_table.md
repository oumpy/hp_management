Title:pandasのpivot_tableを用いた高速データ処理
Date: 2018.11.17
Tags: Data Science Competition
Author: 金子

## 概要
pandasのpivot_tableは強力な機能で、カテゴリごとの集計や計算を高速に行うことができます。

pivot_tableを使った計算で個人的によく使う処理をまとめたものを[kaggle のkernel](https://www.kaggle.com/nadare/feature-engenieering-with-pivot-table)で公開しました。

![1]({attach}images/pandas_pivot_table_figs/pandas_pivot_table.png)

このkernelでは簡単なダミーデータでpivot_tableに対する計算の仕方をまとめた後、実例として[PLAsTiCC コンペ](https://www.kaggle.com/c/PLAsTiCC-2018)の[Starter Kit](https://www.kaggle.com/michaelapers/the-plasticc-astronomy-starter-kit)にあった特徴量の計算をpandasのpivot_tableを用いて高速化しました。

## どんなことができるようになるの？
- カテゴリごとに組み込みの集計関数より高度な関数を適用できる
- カテゴリごとの移動平均をかけるようになる
- 未来のデータを含まないmean_encodingがかける

## どれくらい早くなるの？
上記のkernelはEDA中に実際に僕が書いたコードに少し修正を加えたものですが、

愚直なコード(1時間以上) → groupbyでの処理(2分半) → pivot_table(4秒)

という感じで早くなりました。
