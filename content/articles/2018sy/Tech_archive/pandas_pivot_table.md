Title:pandasのpivot_tableを用いた高速データ処理
Date: 2019.04.14
Tags: kaggle
Slug: pandas_pivot_table
Author: 金子
Summary:

<h1>pandasのpivot_tableを用いた高速データ処理</h1>
<h2>概要</h2>
pandasのpivot_tableは強力な機能で、カテゴリごとの集計や計算を高速に行うことができます。

pivot_tableを使った計算で個人的によく使う処理をまとめたものをkaggle の<strong><a href="https://www.kaggle.com/nadare/feature-engenieering-with-pivot-table">kernel</a></strong>で公開しました。

このkernelでは簡単なダミーデータでpivot_tableに対する計算の仕方をまとめた後、実例として<a href="https://www.kaggle.com/c/PLAsTiCC-2018">PLAsTiCC コンペ</a>の<a href="https://www.kaggle.com/michaelapers/the-plasticc-astronomy-starter-kit">Starter Kit</a>にあった特徴量の計算をpandasのpivot_tableを用いて高速化しました。
<h2>どんなことができるようになるの？</h2>
<ul>
	<li>カテゴリごとに組み込みの集計関数より高度な関数を適用できる</li>
	<li>カテゴリごとの移動平均をかけるようになる</li>
	<li>未来のデータを含まないmean_encodingがかける</li>
</ul>
<h2>どれくらい早くなるの？</h2>
上記のkernelはEDA中に実際に僕が書いたコードに少し修正を加えたものですが、

愚直なコード(1時間以上)→groupbyでの処理(2分半)→pivot_table(4秒)

という感じで早くなりました。
 
