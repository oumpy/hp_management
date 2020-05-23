Title:Particle filter
Date: 2018.10.05
Tags: Machine Learning
Author: 柳澤

今回は**パーティクルフィルター**の紹介をしたいと思います。
といってもやり始めたばっかなので、間違っていたらご指摘お願いします。

## 手順
基本は以下の４つのサイクルを繰り返すだけです。  

1. リサンプリング  
2. 予測  
3. 観測  
4. 尤度計算 重みの更新  

と言っても分かりにくかったんで、イメージで話すと

1. まず粒子（パーティクル）を全体に振りかけます  
2. 振りかけたなかで、あってそうな粒子だけ生き残ってもらい、それ以外は消えてもらいます。  
3. あってそうなものは、尤度（確からしさ）を計算して其れに（おおよそ）従い新たに粒子を撒きなおします（ちょっとランダムウォークさせます）  

ってかんじで対象の動きを推定してくれます。

なんでいきなりパーティクルフィルターの話をしたかというと、こいつは画像解析の分野ではノイズや予想外の動きによって影響を受けにくく、しっかりと標的のものを追ってくれるかなり有用な方法らしいからです。

まあとりあえずやってみよう。ということでpython3+opencvを使いました。
Opencvをpythonで使えるようになったのは結構最近なので、なかなか良い本がないのですが、[公式のチュートリアル](http://labs.eecs.tottori-u.ac.jp/sd/Member/oyamada/OpenCV/html/py_tutorials/py_tutorials.html)が結構役に立ちます。

でも、やってたら結構間違いもあるので、注意してください。あと、有料ですがUdemyでもopencv+pythonの講座があるので、試してみても良いかもしれません。
先ほど紹介したコンピュータビジョン最先端ガイドも理論がわからない時便利です。

うまくいった例はネットに大量にあるので、興味のある方は検索してみてください。

ちなみにパーティクルフィルターでの一番の肝は、↑の4．尤度関数の設定の仕方です。
今回は色の尤度の他に、距離でも尤度を設定して（つまりある時点で散らばっている粒子の重心からの距離を考えるということ）それらをかけたものを最終の尤度のしました。
ブラウン運動ではなく、今回のように細胞の動きを追う場合では、次のフレームで動きそうなところの尤度を大きくすればもっと正確に動きを追いかけられるみたいです。難しい。。。

## 実装
以下、[『Udemy 【Pythonで学ぶ】OpenCVでの画像処理入門』](https://www.udemy.com/pythonopencv/)のコードを参考にさせていただきました。

``` python
import cv2
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

for i in tqdm(range(100)):
    time.sleep(0.1) #プログレスバーの表示

args = sys.argv #コマンドライン引数

cap = cv2.VideoCapture("/movie/lps.avi")

if cap.isOpened() == False:
    sys.exit

ret, frame = cap.read()
h, w =frame.shape[:2]#大きさを取得

fourcc = cv2.VideoWriter_fourcc(*"mp4v")
output_dst = cv2.VideoWriter("/output/test[{0},{1}].m4v".format(args[1],args[2]),fourcc,5.0,(w,h))#動画出力の設定

np.random.seed(100)#乱数の初期化,毎回同じ乱数になる
Np = 50#粒子の数
obj = [int(float(args[1])),int(args[2])] #目的の（追いかける）標的の座標0~512
WD = 100

px = np.zeros((Np),dtype=np.int64)#粒子のx座標
py = np.zeros((Np),dtype=np.int64)#粒子のy座標
lc = np.zeros((Np))#粒子の色の尤度
ls = np.zeros((Np))#粒子の空間の尤度
lt = np.zeros((Np))#粒子の尤度total
index = np.arange(Np)

#objの周りに撒く
d = 10
px = np.random.normal(obj[0], d, Np).astype(np.int)
py = np.random.normal(obj[1], d, Np).astype(np.int)

j = 0
artists = []
while True:
    ret, frame = cap.read()#１枚読み込み
    if ret == False:
        break#最後になったらループから抜ける

    gx = np.average(px)
    gy = np.average(py)#１フレーム前の粒子の重心

    for i in range(Np):
        lc[i] = frame[py[i],px[i]][1] / 255.0#色の尤度
        ls[i] = np.exp(-((px[i] - gx) ** 2 + (py[i] - gy) ** 2)/(WD ** 2))
        lt[i] = lc[i] * ls[i]
    lt = lt / lt.sum()

    pnew_index = np.array(random.choices(population=index,weights=lt,k=Np))
    pxnew = px[pnew_index] + np.random.randint(-15,15,Np)
    pynew = py[pnew_index] + np.random.randint(-15,15,Np)

    plt.hist(lt)

    #リサンプリングした,ある程度ランダムウォーク
    px = np.where(pxnew > w-1, w-1, pxnew)
    py = np.where(pynew > h-1, h-1, pynew)
    px = np.where(px < 0, 0, px)
    py = np.where(py < 0, 0, py)#ランダムウォークで画面外に出る場合の処理

    for i in range(Np): #画像の中に粒子を描く
        cv2.circle(frame,(px[i],py[i]),1,(255,255,255),1)

    cv2.imwrite("/output/test_tiff/test" + str(j) + ".tif" ,frame) #tiffでも保存
    j = j + 1

    output_dst.write(frame)
```

## 結果
こんな感じです。

![1]({attach}images/particle_filter_figs/result.gif)

最初は尤度を計算するとき以下のようにやってました。
```python
pxnew = np.array(random.choices(population=px,weights=lt,k=Np)) + np.random.randint(-15,15,Np)
pynew = np.array(random.choices(population=py,weights=lt,k=Np)) + np.random.randint(-15,15,Np)
```

これではxとyを別々に計算してしまっているので、微妙に結果がおかしい感じになってました。（点が四角っぽくなる）
確かにそれはそうか。ランダムに選ぶのは１回でいいはず。。。。

こういうとき、適当にindexとかおいてやるとうまくいくんですね。
今回もpnew_indexをおいてしまうという感じでやってます。

```python
pnew_index = np.array(random.choices(population=index,weights=lt,k=Np))
pxnew = px[pnew_index] + np.random.randint(-15,15,Np)
pynew = py[pnew_index] + np.random.randint(-15,15,Np)
```

## 参考文献
- <a href="https://www.amazon.co.jp/%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%82%BF%E3%83%93%E3%82%B8%E3%83%A7%E3%83%B3%E6%9C%80%E5%85%88%E7%AB%AF%E3%82%AC%E3%82%A4%E3%83%891-CVIM%E3%83%81%E3%83%A5%E3%83%BC%E3%83%88%E3%83%AA%E3%82%A2%E3%83%AB%E3%82%B7%E3%83%AA%E3%83%BC%E3%82%BA-%E5%80%89%E7%88%AA-%E4%BA%AE/dp/4915851346">コンピュータビジョン最先端ガイド</a>
