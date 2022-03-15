Title: Kaggle Tabular Playground Series - Dec 2021
Date: 2022.03.10
Modified: 2022.03.10
Tags: Data Science Competition
Author: 石本

## やったこと
- baseline作成（[参考記事](https://www.kaggle.com/chryzal/features-engineering-for-you)）　public score = 0.95102
    - **Sequential モデル**を使用
    - データフレームからあまり重要でなさそうな行や列を削除
    - Aspectはコンパスの方向のことなので0~359の範囲に入っていないものは±360して**正規化**
    - x_dist_hydrlgy, y_dist_hydrlgyをマンハッタン距離とユークリッド距離に変換
    - hilshadeが0~255に入っていない場合は外れ値とみなして0または255に**正規化**
    - RobustScalerでデータのスケール変換・移動
    - **メモリ節約**のためにデータフレームのデータ型を変換
    - 試運転のため、epoch数を200 -> 1に変更
- soil_typeの合計とwilderness_areaの合計を特徴量に加えた。（[参考記事](https://www.kaggle.com/c/tabular-playground-series-dec-2021/discussion/292823)）
- 進捗を見ながら学習を進めるために、tqdmを導入した。
- epoch数を1 -> 200に変更　public score = 0.95664

## 結果・感想
結果は、273位でした。baselineに手を加えられたのは初めてだったので成長を感じることができました。とはいえ、baselineのモデルを完全に理解できたわけではなかったので有効な改善ができたかは不明です。
これからも積極的にコンペに参加していきたいです。

## 参考にした本
今回のコンペでは[『kaggleで勝つデータ分析の技術』](https://www.amazon.co.jp/Kaggle%E3%81%A7%E5%8B%9D%E3%81%A4%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E3%81%AE%E6%8A%80%E8%A1%93-%E9%96%80%E8%84%87-%E5%A4%A7%E8%BC%94/dp/4297108437)を参考にし、モデルの改善に取り組みました。特に、「特徴量の作成」の項が参考になりました。部分的にしか読めていないので、今後コンペに参加していく中で完読を目指していきます。
