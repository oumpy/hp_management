https://oumpy.github.io/
jupyterもそのまま出力できるようにしました。
写真ファイルのパス指定や内部リンク調整、（Dateは自動で今の日時打ち込まれるようにしました。）手打ちをした部分（Author, Category, Tagsなど）のチェックが必要なのですが、一人では流石に厳しいので、後輩に任せていきますねー！もう働きたくないよぅ。

基本はcontentディレクトリにmdを書いていきます(htmlで書いても大丈夫)。imagesフォルダ作ったので、そこに写真おいてパス通せば表示できる。

jupyter notebookに関してはcontent/jupyter_labに入れて、メタデータのファイル作って、notebookの最初のセルにメタデータの内容と同じものを書けば大丈夫です。

以下カスタムで書いたスクリプト
bash create.sh arg1 arg2 でarg1.mdという名前でarg2というカテゴリのマークダウン（必要なメタデータ付き） 生成して、カテゴリー名のディレクトリに送ります。
bash push.shを実行すれば、全てのファイルをhtmlにコンパイルして、/outputへ、そしてgithubにプッシュします。
カスタムファイルは自分の趣味で色々いじってくれたらいいと思う〜。効率的な方法思いついたら教えて。

やってほしい調整は、写真ファイルのパス指定や内部リンク調整（現状ネットワーク上のパスを参照しているので表示できていますが、wordpressを潰したら消えるので、github上にあげておきましょう。）、（Dateは自動で今の日時打ち込まれるようにしました。）手打ちをした部分（Author, Category, Tagsなど）、authorも


## Jupyter notebookをHTML出力できるようにする。

pluginを追加
```
mkdir plugins
git clone https://github.com/danielfrg/pelican-ipynb.git plugins/ipynb
```
[この記事](https://qiita.com/driller/items/49a990cbdfb51afed620)に従って、pelicanconf.pyを書き換える。
