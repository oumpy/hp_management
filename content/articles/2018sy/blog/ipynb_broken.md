Title: ipynb破損時の対応
Date: 2018.09.15
Modified: 2020.05.16
Tags: Python
Author: 安水

先日１ヶ月以上を費やしてまとめていたipynbファイルが破損して開けなくなってしまい、やさぐれていました。
なんとかならないか調べても、壊れたファイルは戻ってこないよ的な塩対応しかなかったので絶望していたのですが、pythonでなんとかできるのではないかとあがいてみました。
なんとか復旧したのでその方法を共有します。

![image01]({attach}./images/ipynb_broken_figs/image01.png)

まず開こうと思ったらこんな感じでした。
絶望感しかありません。
ipynbファイルの弱点を見たような気がしました。

そこで、新しいipynbを作り、jsonを読んでみました。

![image02]({attach}./images/ipynb_broken_figs/image02.png)

てことで677378行目付近を見てみる。

```bash
less -N BAP_note.ipynb
```

Nオプションで行数が表示されます。
そして、`677378G` と打つと677378行目に行けます。

![image03]({attach}./images/ipynb_broken_figs/image03.png)

たしかに表示が崩れているので、`v` で編集モード（vimモード）を立ち上げ、`dd` でけしていきました。
`:x` で保存できます。
（vimじゃなくてもemacs等でも行ける。）

復活しました。
うれしい!!

![image04]({attach}./images/ipynb_broken_figs/image04.png)
