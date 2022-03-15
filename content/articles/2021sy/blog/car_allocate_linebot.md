Title:”車割り”を最適化
Date:2022.03.13
Modified:2022.03.13
Tags:自動化,競技プログラミング
Authors:富本,安部

”車割り”という作業を競プロ力で自動化→lineのbotで運用　ということをやってみました。


### 車割り　とは

部活の帰りに親切な部員数名が他の部員を最寄駅付近などまで送り届けてあげるとき、どの部員が誰の車に乗るか割り振る　というシステムです。<br>
「自分の帰り道から大きく外れる場所を通るのはなるべく避けたいが、できるだけ多くの人数が乗車しているのが望ましい」という考えのもと毎回数人の知恵を結集して割り振りが決められています。<br>
この作業にはかなりの時間がかかるので自動化できれば嬉しいはずです。


### 車割りを最適化問題に落とす(富本part)

車割り問題は重み付き二部グラフ最大マッチング問題に帰着できます。<br>
まず各車と各部員に対応する頂点を用意し、次に各車に各部員を乗せたときの利益(たとえば送っていける距離)を事前に決めてその利益分の辺を張っていきます。<br>
すると、このグラフは二部グラフとなり、できるだけ利益が大きくなるようにマッチングさせる問題になります。<br>
このままだと各車に一人の部員しかマッチングできませんが、たとえばある車に運転手以外に人が４人乗れる場合は、その車に対応する頂点を４つ用意しておくことで、最大４人の部員をマッチングさせられます。<br>
(実際には頂点を４つ用意するのではなく、後で説明する最小費用流アルゴリズムにおいて張る辺の容量を4にします。)


### アルゴリズムの解説　(富本part)

重み付き二部グラフ最大マッチングは最小費用流アルゴリズムを用いると以下のようにして解くことができます。<br>
1. 各車と各部員の頂点のほかに２つの頂点SとTを用意し、十分大きい整数をINFと置く<br>
2. Sから各車に辺(容量:その車に乗れる人数、コスト:0)を張る<br>
3. 各車から各部員に辺(容量:1、コスト:INF-利益)を張る(ただし、利益＜0の場合は辺を張らない)<br>
4. 各部員からTに辺(容量:1、コスト:0)を張る<br>
5. SからTにフローを流せるだけ流す<br>
6. 車(Aとする)から部員(Bとする)に張られており(正の)フローが流れているようなすべての辺について、部員Bは車Aに乗ることにする<br>


### アルゴリズムをline botに載せる

line botの作成は[公式のドキュメント](https://developers.line.biz/ja/docs/messaging-api/overview/)が丁寧なのでそれに従いました。
"サンプルボットを作成する"の項目にある[pythonのtoolkit](https://github.com/line/line-bot-sdk-python)を使用します。
まず、[echobotの例](https://github.com/line/line-bot-sdk-python/tree/master/examples/flask-echo)をそのまま実行して動作確認をしました。


次に、実装されたアルゴリズムを組み込んでいきます。
1. 各部員に対して学年、最寄り駅の経度緯度、車に対して乗車可能人数を初期値として持っておきます。
```python
grades = {"name1":1,"name2":2....}
places = {"start1":(000,000),"name1":(000,000),,"name1":(000,000)...}
car = {"name1":4,"name5":3,...}

```

2. 入力に応じて出力を変えるように場合わけをします。
①割り振り②初期値の更新③初期値の一部確認④helpを想定される入力の様式としてそれ以外は"helpとうって期待する入力形式を確認するように"とエラーを出します。

①割り振り、では　車のスタート地点、割り振られる部員を入力として
```python
@app.route("/callback", methods=['POST'])
def callback():
    ~~~~
    input,has_error = check_input(input)
    if has_error:
        posst = "入力が違います..."
    else:
        members,start_place,absent,guest=input
        posst = assignment(members=members,start_place=start_place,absent=absent,guest=guest)#重み付き二部グラフ最大マッチング実行


    line_bot_api.reply_message(
            event.reply_token,
            TextSendMessage(text=posst)
        )
#>>posst
#out_of_car:[name1,name2,name3, ...]
#name4's car: [name5,name6,name7]
#name8's car: [name9,name10,name11]
#name12's car: [name13,name14,name15]
```
という出力を得ます。

②初期値の更新、では変数をglobalにすると動きました。
```python

X = ....

@app.route("/callback", methods=['POST'])
def callback():
    ~~~~
    global X
    X = new_X
```

実際には、membersはあらかじめ数個に分かれているgroups(X=[ ],Y=[ ],Z..)のlistを使用することで入力の文字列をできるだけ短くしています。
　

### 感想

いままで触ってみたら楽しいかもしれないけれどとっつきにくい、と感じていたflaskやbotに1mmでもタッチできたのは良い経験でした。<br>
競技プログラミングが実生活に役立つ瞬間がみれて感動しました。<br>
最適化の部分はcolabに貼り付けて自分で変数を書き換え→実行すれば動くようになっていましたがほとんど利用されなかったこともあり、休みを利用してlineのbotに載せました。
botの部分はサンプルをいじっただけでしたが、競プロerの力を借りることで巷によるあるbotよりも良さそうなものが作れたのではないかと思います。

#### 参考
https://www.youtube.com/watch?v=jBsvdgFMZtg
heroku login:https://qiita.com/WEByamori/items/4cab4be8ce336ed3e4d6


