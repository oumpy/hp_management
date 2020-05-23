Title:Pythonで可視化入門
Date: 2019.01.18
Tags: Python
Author: 宮崎

Pythonでいい感じのグラフを書いてみたい！でも面倒！よくわからない！

これを読めばそんなあなたも簡単にいい感じのグラフがかける！

この記事では、Pythonのライブラリである、定番のmatplotlibとseaborn、pandasを使った可視化、インタラクティブなplotly、複数グラフ表示が簡単にできるラッパーであるholoviewsの簡単な解説をします。

## Pythonの可視化ライブラリ
### matplotlib：<u><a href="https://seaborn.pydata.org/">https://matplotlib.org/</a></u>(公式)
定番のmatplotlibです。まずはギャラリー ( <u><a href="https://matplotlib.org/gallery/index.html">https://matplotlib.org/gallery/index.html</a></u>)でいろいろなグラフを見てみましょう。グラフがたくさん並んでいて楽しいですね！

### seaborn：<u><a href="https://seaborn.pydata.org/">https://seaborn.pydata.org/</a></u>(公式)
matplotlibをベースにしたseabornです。こちらもギャラリー(<u><a href="https://seaborn.pydata.org/examples/index.html">https://seaborn.pydata.org/examples/index.htm</a><a href="https://seaborn.pydata.org/examples/index.html">l</a></u>)を見てみましょう。ヒートマップがおしゃれですね！

### pandas：<u><a href="https://pandas.pydata.org/index.html">https://pandas.pydata.org/index.html</a></u>(公式)
表計算で便利なpandasです。DataFrame形式のデータをpandas.DataFrame.plot(<u><a href="http://pandas.pydata.org/pandas-docs/stable/visualization.html">http://pandas.pydata.org/pandas-docs/stable/visualization.html</a></u>)でグラフを書くのが一番手間がかからない気がします。bar(棒グラフ)やhistogramやscatter(散布図)を作るならこれで十分なことが多いです。DataFrameの.describe (<u><a href="https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.describe.html">https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.describe.html</a></u>)で統計量が簡単に計算できるので、グラフがあっているかの確認も一連の流れでできます。

### plotly：<u><a href="https://plot.ly/">https://plot.ly/</a></u>(公式)
インタラクティブに動かせるplotlyです。<u><a href="https://plot.ly/python/">https://plot.ly/python/</a></u>(Python用) を見てみると、なんとグラフを動かせます！楽しい！Kaggleでもよく見ます。データ構造を調べるのに便利だからかもしれないですね。(例<u><a href="https://www.kaggle.com/andresionek/what-makes-a-kaggler-valuable?utm_medium=social&amp;utm_source=twitter.com&amp;utm_campaign=Weekly-Kernel-Awards">https://www.kaggle.com/andresionek/what-makes-a-kaggler-valuable?utm_medium=social&amp;utm_source=twitter.com&amp;utm_campaign=Weekly-Kernel-Awards</a></u>) Dashと組み合わせると ( <u><a href="https://plot.ly/products/dash/">https://plot.ly/products/dash/</a></u>) 簡単にWebアプリにできます。

### holoview：<u><a href="http://holoviews.org/">http://holoviews.org/</a></u>(公式)
Pythonの可視化ツールはHoloViewsが標準になるかもしれない(<u><a href="https://qiita.com/driller/items/53be86cea3c3201e7e0f">https://qiita.com/driller/items/53be86cea3c3201e7e0f</a></u>)とまで言われるholoviewは、matplotlibやplotly、bokehを簡単に使えるようにするラッパーです。特筆すべき機能として、introduction(<u><a href="http://holoviews.org/getting_started/Introduction.html">http://holoviews.org/getting_started/Introduction.html</a></u>)をみるとわかるのですが、複数グラフ表示は足し算(例：`Compositional Layouts layout =scatter +hv.Histogram`)、オーバーレイは掛け算(例：`Compositional Overlays image +image*points`)で定義できとてもシンプルに書くことができます。

可視化のまとめ(kaggle：<u><a href="https://www.kaggle.com/maheshdadhich/strength-of-visualization-python-visuals-tutorial">https://www.kaggle.com/maheshdadhich/strength-of-visualization-python-visuals-tutorial</a></u>)もよくまとまっていて非常に勉強になります。

そのほかにも、tensorflowの embedding <u><a href="https://www.tensorflow.org/guide/embedding">https://www.tensorflow.org/guide/embedding</a></u>から Mnist や Word2Vecのデモ <u><a href="http://projector.tensorflow.org/">http://projector.tensorflow.org/</a></u>(重いかも) や <u><a href="https://distill.pub/">https://distill.pub/</a></u>も動かしてみると楽しいです！

## まとめ
- holoview最強説(plotlyも使えるラッパーなので)
- 楽にしたいならpandasのplot、動かしたければplotly！

keywordでググると詳しい記事がたくさん出てきます。上記以外の可視化手法や面白いものあれば教えてください！
