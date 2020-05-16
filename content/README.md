# サイト内容（contentディレクトリ）の内容と書き方

サイトのソースはこのディレクトリ (`content`) に入っている。

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management/content
```

## サイト内容の更新

### 新しいページ・記事の作り方

`content` ディレクトリに移動し、Pythonスクリプト `create.py` で雛形を作ります。

あるいは、 (その方が楽であれば) `content/template.md` を手動でコピーして編集してもかまいません。

``` bash
$ cd hp_management/content
$ python create.py
```

質問に従って必要項目を順に入力します。
ブログ記事の場合、`date` (日付) と `slug` (ファイル名) は入力必須。
その他は全部デフォルトにして後から書き換えるのでも構いません。
ファイルは `content/articles/(schoolyear)/(category)/(slug).md` として作成されます。

- (category) はデフォルトで `blog`。
  他に `news`, `page` (サイトの固定ページ) を指定できる。
- (schoolyear) は `date` から自動で計算される。
  例えば `date` が `2020.03.15` なら (schoolyear) は `2019sy`。
- 質問に答える代わりに、`--slug`, `--date` などのコマンドラインオプションを使うことも可能。
  詳細は `$ python create.py --help` とすると出てくる。
- `Author:` の下を1行空け、その後にMarkdown原稿を貼り付ける。
  タイトル (`#` タグ) は自動生成されるので削り、分節タグは原則全て `##` 以下を使う。
  原稿ファイルがすでにある場合、`--content` オプションでファイルを渡すこともできる。

#### サマリーについて
記事一覧に表示される各記事のサマリーは、記事本文冒頭の140文字まで（タグ・空白・改行など除く）が使用されます。
ただし、

- 記事原稿ファイルのメタデータに `Summary:` を加えると、そちらがサマリーとして優先的に使用されます。
- サマリーを本文中の特定範囲 (140文字以内) にしたい場合は、
```
<!-- PELICAN_BEGIN_SUMMARY -->
```
```
<!-- PELICAN_END_SUMMARY -->
```
の2行を本文中に挿入してください。
(常にセットである必要はなく、いずれも省略可能です。)
例えば本文を
```
## はじめに
<!-- PELICAN_BEGIN_SUMMARY -->
今回はすごいことをやってみました。
<!-- PELICAN_END_SUMMARY -->
## すごいことの詳細
とてもすごい。どこがとてもすごいかというと....
```
のようにすると、サマリーは「今回はすごいことをやってみました。」となります。
`<!-- ... -->` は本文中ではコメントとして無視されます (表示されません)。

#### Jupyter Notebookの扱い

jupyter notebookに関しては他の記事（mdファイル）と同じ場所に入れ、さらに同じ場所にメタデータファイル（`myarticle.ipynb` の場合は`myarticle.nbdata`）を置いてmdファイルと同様のメタデータを書きます。
`articles/2018sy/blog` の `lorentz.ipynb` および `lorentz.nbdata` を参考にしてください。

#### タグの付け方
現在は10タグ設定している。基本的にはこの中から最も当てはまるタグを1つ選ぶ。ただし、新しい分野も歓迎します（相談してください）。また、大文字・小文字等の違いに注意。
- `Bioinformatics`：バイオインフォマティクス関連
- `Machine Learning`：機械学習関連 (Deep Learningなど)
- `Statistics`：統計学（機械学習との違いは若干曖昧）
- `Data Science Competition`：KaggleやSignateなど
- `Python`：Pythonに関すること（上記に当てはまらない場合）
- `Shell script`：シェルスクリプトの使い方など
- `GitHub`：Githubの使い方など
- `競技プログラミング`：Atcoderなど
- `論文関連`：論文の探し方、書き方、管理の仕方などのtips
- `海外留学`：海外留学の報告など

#### リンクなど

ソース内のmarkdownファイルなどへのリンクが自動的にhtmlへのリンクに変換される機能があります（画像などはそのまま）。
リンクの書き方としては、

- 通常通りにリンクを書く。
- その後、リンクの頭に「おまじない」を入れる。
  - 画像など、「その記事に付随して初めてアップするファイル」は`{attach}`
  - すでに存在している他記事などの場合は`{filename}`

という順序になります。
以下、具体例です。

##### 画像の設置と読み込み

記事に貼る画像は、

```markdown
content/articles/(schoolyear)/(category)/images/(filename)_figs/(imagefile)
```

として保存します。

例えば、`sugoikiji.md` に画像ファイル `sugoigazou.png` を読み込みたい場合、記事本体が

```markdown
content/articles/2020sy/blog/sugoikiji.md
```

であれば、画像ファイルは

```markdown
content/articles/2020sy/blog/images/sugoikiji_figs/sugoigazou.png
```

のようにします。

これで、 `sugoikiji.md` の中で

```markdown
![すごい画像]({attach}./images/sugoikiji_figs/sugoigazou.png)
```

と書けば読み込まれます。

##### 他記事等へのリンク

例えば同じ年度の他の関連記事などへのリンクを貼りたい場合は、

```markdown
[前回の記事]({filename}./my_previous_article.md)
```

のようにします。
`{filename}` がおまじないです。

過年度の記事であれば、少しややこしいですが、ディレクトリを辿って

```markdown
[以前の記事]({filename}../../2019sy/blog/old_article.md)
<!-- 相対パス指定 -->
```

または

```markdown
[以前の記事]({filename}/articles/2019sy/blog/old_article.md)
<!-- content/ をトップとする絶対パス指定 -->
```

のようにします。

### 記事を投稿する方法

投稿はGitHubのプルリクエストによりますが、方法は2つあります。

方法2の方が若干ですが初心者向けで、プレビューも使えるので推奨ですが、github/oumpy のメンバーに登録される必要があります。
未登録でも方法1は可能です。

#### 方法その1 (レポジトリのフォークを使う方法)

「フォークからのプルリクエスト」によって行います。

以下簡単に。各段階での詳しいことはネット上にもたくさんありますが、わからなければSlackで質問してください。

##### 事前に必要な設定

###### Git/GitHub一般の設定

（すでにGitHubを使っている人は飛ばして次へ）

- GitHub アカウント ( `hoge` とする。以下自分のものに読み替え) と Gitクライアント (SourceTreeなど) を設定する。
- 自分でGitHub上に作ったレポジトリにプッシュできることを確認する。
  認証関係はハマりどころなので資料を確認しながら慎重に。

###### 本レポジトリに関する設定

- `oumpy/hp_management` をクローンしてローカルレポジトリを作成。
  リモートレポジトリ (`oumpy/hp_management`) の名前はデフォルトで `origin` となる。

  （本README冒頭のコマンド。実行済みの場合はそのままでOK。）

- GitHubの `oumpy/hp_management` で、右上にある「fork」から自分のアカウントにフォーク (コピー) を作成する。

- `hoge/hp_management` を2つ目のリモートレポジトリとして設定。2つのリモートレポジトリには適切な名前をつけて区別する。
  
  - `hoge/ph_mangement` を `origin`、`oumpy/hp_management` を  `upstream` 。
  - `hoge/ph_mangement` を `hoge`、`oumpy/hp_management` を  `oumpy` 。こちらはGitHub上での呼称と整合する。
  
  どちらでもOK。以下では後者で説明する。

##### 記事投稿ごとにやること

- masterブランチを `oumpy/master` からプルして `master` と `oumpy/master` を一致させた後、この点に新しいブランチを作成。
  ここでは `sugoikiji` とする。
- 書いた記事ファイルをブランチ `sugoikiji` にコミットする。
  記事タイトルなどをコミットコメントとして書く。
- ブランチ `sugoikiji` を自分のレポジトリにプッシュする ( ブランチ `hoge/sugoikiji` ができる ) 。
- webブラウザで自分のレポジトリまたはoumpyの元レポジトリに行き、自分の `sugoikiji` から`oumpy/master` へのプルリクエストを作成。
  必要な説明などを同時にコメントとして書く。
  作成が完了すると、ラベル `article` が自動的に付けられる。

通常は以上でやることは終わりです。

- 問題なければ幹部が承認し、アップされます。
  また「マージされました」というお知らせがアカウントに紐付けられたメールアドレスにGitHubから送られます。
- エラーその他で修正が必要で、本人に依頼する必要がある場合は修正リクエストが発行されます。その場合もメールでお知らせが送られます。

#### 方法2 (oumpy上に直接プッシュする方法)

##### 事前に必要な設定

###### Git/GitHub一般の設定

方法1と同じです。すでにGitHubを使っている場合は必要ありません。

###### アカウント登録依頼

`github/oumpy` へのメンバー登録を管理者に依頼する。招待メールが送られてくるので受諾の手続きをする。

###### 本レポジトリに関する設定

`oumpy/hp_management` をクローンしてローカルレポジトリを作成する。
(本README冒頭のコマンドを実行していればそれでOK。)

リモートの名前は `origin` のままでも良いが、何かの時に方法1を使ったりHP開発に参加したりしたい場合は、方法1と同じように `oumpy` や `upstream` に改称。

以下では方法1と同様に `oumpy` と名づけたものとする。

##### 記事投稿ごとにやること

- masterブランチを `oumpy/master` からプルして `master` と `oumpy/master` を一致させた後、この点に新しいブランチを作成。
  **この時、名前の先頭に `article/` をつける。**
  ここでは `article/sugoikiji` とする。
- 書いた記事ファイルをブランチ `article/sugoikiji` にコミットする。
  記事タイトルなどをコミットコメントとして書く。
- ブランチ `article/sugoikiji` をリモートレポジトリ `oumpy` にプッシュする ( ブランチ `oumpy/article/sugoikiji` ができる ) 。
  (ここでブランチ名が正しくないと、権限がないと言われてプッシュに失敗します。)

- webブラウザで`https://github.com/oumpy/hp_management/`に行き、 `article/sugoikiji` から`master` へのプルリクエストを作成。
  必要な説明などを同時にコメントとして書く。
  作成が完了すると、ラベル `article` が自動的に付けられる。
  またプレビューへのリンクが自動投稿される。

通常は以上でやることは終わりです。

その後の承認・マージや修正依頼については方法1の場合と同じです。

##### 記事のプレビュー

この方法2では、プルリクエストの前、ブランチをリモートにプッシュした時点で、そのブランチの内容から、サイトのプレビューが自動的に生成される。
今の場合であれば、 `https://oumpy.github.io/previews/refs/heads/article/sugoikiji/` にブラウザでアクセスすると、プレビューを見ることができる。
上記の通り、プルリクエストを作成すると、このリンクは当該プルリクエストのページに自動で投稿される。

ブランチをプッシュしてからプレビューがアップされるまでは通常1分程度 (場合によっては10分かかることなどもある)。
同じブランチを更新してプッシュするたびにプレビューも更新される。
またリモートブランチを削除するとプレビューも削除され、現在のサイトへのリンクで置き換えられる。
(中身は全く同一だが、URLはリダイレクトされない。
 プレビューが削除されているかどうかは、サイトのタイトルに「 (テスト用ページ) 」とついているかどうかでわかる。)

## サイト全体の情報設定

`content/contentconf.py` に、サイト表題・副題やリンク先などの情報が書かれている。
**この設定は `pelicanconf.py` の一部であり、また `pelicanconf.py` 本体の記述よりも優先される。**

同様に `content/contentpublishconf.py` もあり、テスト時にはない方がよい設定を記述する。

