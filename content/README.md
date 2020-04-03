# サイト内容（contentディレクトリ）の内容と書き方

サイトのソースはこのディレクトリ (`content`) に入っている。

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management/content
```

## サイト内容の更新

### 新しいページ・記事の作り方

`content` ディレクトリに移動し、独自スクリプト `create.sh` で雛形を作ります。

あるいは、 (その方が楽であれば) `content/template.md` を手動でコピーして編集してもかまいません。

#### page (固定ページ) を作成する場合

```bash
$ cd hp_management/content
$ sh create.sh (filename)
```

`(filename)` に生成するファイル名を指定します (内容に応じて適切に)。
`content/pages/` の下にメタデータ入りのMarkdownファイルが生成されます。

これに内容を追記します。htmlも認識されます。

#### article (ニュースやブログ記事) を作成する場合

``` bash
$ cd hp_management/content
$ sh create.sh (filename) (category) (tags)
```

Pageの場合と同じコマンド。(category)を指定するとarticleになり、`content/articles/(schoolyear)/(category)/(filename).md`として作成されます。

- (category) はNewsかBlog。
  (category) にpageまたはpagesを指定するとpageになる。
- (tags) は省略可能、ただしarticleの場合は編集の際に必ず入れること。
  (この仕様は変更するかも。)
- (category) は各ファイルのメタデータには記載されない。
  mdファイルを入れるディレクトリの名前がカテゴリとして認識される。
  (記載した場合はそちらが優先されるが、記載しないでください。)
- 自動作成されたファイルは `Title` と `Author` が空欄なので、埋める。
- `Author` の下を1行空け、その後にMarkdown原稿を貼り付ける。
  タイトル (`#` タグ) は自動生成されるので削り、分節タグは原則全て `##` 以下を使う。

##### Jupyter Notebookの扱い

jupyter notebookに関しては他の記事（mdファイル）と同じ場所に入れ、さらに同じ場所にメタデータファイル（`myarticle.ipynb` の場合は`myarticle.nbdata`）を置いてmdファイルと同様のメタデータを書きます。
`articles/2018sy/Blog` の `lorentz.ipynb` および `lorentz.nbdata` を参考にしてください。

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
content/articles/2020sy/Blog/sugoikiji.md
```

であれば、画像ファイルは

```markdown
content/articles/2020sy/Blog/images/sugoikiji_figs/sugoigazou.png
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
[以前の記事]({filename}../../2019sy/Blog/old_article.md)
<!-- 相対パス指定 -->
```

または

```markdown
[以前の記事]({filename}/articles/2019sy/Blog/old_article.md)
<!-- content/ をトップとする絶対パス指定 -->
```

のようにします。

### 記事を投稿する方法

「フォークからのプルリクエスト」によって行います。

以下簡単に。各段階での詳しいことはネット上にもたくさんありますが、わからなければSlackで質問してください。

#### 事前に必要な設定

##### Git/GitHub一般の設定

（すでにGitHubを使っている人は飛ばして次へ）

- GitHub アカウント ( `hoge` とする。以下自分のものに読み替え) と Gitクライアント (Source Treeなど) を設定する。
- 自分でGitHub上に作ったレポジトリにプッシュできることを確認する。
  認証関係はハマりどころなので資料を確認しながら慎重に。

##### 本レポジトリに関する設定

- GitHubの `oumpy/hp_management` で、右上にある「fork」から自分のアカウントにフォーク (コピー) を作成する。
- `hoge/hp_management` をクローンしてローカルレポジトリを作成。
  リモートレポジトリ (`hoge/hp_management`) の名前はデフォルトで `origin` となる。
- `oumpy/hp_management` を2つ目のリモートレポジトリとして設定。名前は `upstream` とする。
  - 名前は `upstream` / `origin` の代わりに、oumpyの方を `oumpy` 、自分のものを `hoge` として区別してもよい。
    こちらはGitHub上での呼称と整合する。お好みで。

#### 記事投稿ごとにやること

- masterブランチを `upstream/master` からプルして `master` と `upstream/master` を一致させた後、この点に新しいブランチを作成。
  ここでは `sugoikiji` とする。
- 書いた記事ファイルをブランチ `sugoikiji` にコミットする。
  記事タイトルなどをコミットコメントとして書く。
- ブランチ `sugoikiji` を自分のレポジトリにプッシュする ( ブランチ `origin/sugoikiji` ができる ) 。
- webブラウザで自分のレポジトリまたはoumpyの元レポジトリに行き、自分の `sugoikiji` から`oumpy/master` へのプルリクエストを作成。
必要な説明などを同時にコメントとして書く。
記事の投稿・修正の場合、`article` ラベルをつける。

通常は以上でやることは終わりです。

- 問題なければ幹部が承認し、アップされます。
  また「マージされました」というお知らせがアカウントに紐付けられたメールアドレスにGitHubから送られます。
- エラーその他で修正が必要で、本人に依頼する必要がある場合は修正リクエストが発行されます。その場合もメールでお知らせが送られます。

## サイト全体の情報設定

`content/contentconf.py` に、サイト表題・副題やリンク先などの情報が書かれている。**この設定は `pelicanconf.py` の一部であり、また `pelicanconf.py` 本体の記述よりも優先される。**

同様に `content/contentpublishconf.py` もあり、テスト時にはない方がよい設定を記述する。

## ToDo

- Member,Calenderなどは必要か？
