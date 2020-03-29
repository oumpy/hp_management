Title:寿司打タイピング自動化
Date: 2018.09.20
Tags: Python
Slug: sushi_typing
Author: 安水
Summary:

<iframe width="560" height="315" src="https://www.youtube.com/embed/SqOO9I1tFjk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

寿司打の自動化についてです。上の動画のようになりました。寿司打は結構昔からあるタイピング練習サイトです。寿司打は<a href="http://typing.sakura.ne.jp/sushida/">こちら</a>です。
1万円コースで1万円超えない人はタイピング遅いので特訓しましょう。

自動化の流れとしては

- スクリーンショットを取る(PIL)
- OCRで文字認識をする(pyocr)
- キーボード入力(pyautogui)

です。
ご認識を防ぐため、調整用のところでできるだけ正確にローマ字部分が入るようにスクリーンショットの座標をあわせないといけません。

今回はMBPをつかって一万円コース63120 円でした。OCRの部分が律速になっていそうなので、リソースかアルゴリズムでここを高速化するとプラトーになるのではと考えています。
一晩中自動でゲームを繰り返してスコアをメモし続けるように書き直したので、少し冗長なコードになっています。
ところで、ランキングをチェックすると、抜かれていたので、だれかリベンジしてください。

<https://gist.github.com/yyoshiaki/a9a2eea1105a654e6b9beaaf2b854871>
