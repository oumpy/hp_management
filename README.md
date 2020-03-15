# [阪医Python会HP](https://oumpy.github.io/)

原作：平岡悠  
現メンテナ：AtamaokaC

## 全体の仕組み

- ファイルが全て入った、[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)とoutputディレクトリのHTMLファイルのみが入った、[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)がある。

- 基本は管理用レポジトリをプルして、設定を変更したり、content内の新しいブログを追加したりする。完成したら、管理用レポジトリにプッシュするとともに、出力用レポジトリにもoutputの中身をpushする。この繰り返し。

- 会員は本レポジトリ（hp_management）を自分のGitHubアカウントにfork、`content/`下のMarkdownファイルを編集・追加してプルリクエストを送ることを推奨。

- 技術的説明は下に色々あるが、管理に携わる人以外は**「新しいページの作り方」以外は知らなくてもよい**。（「導入」等もやらなくてよい。）

  Markdownとgitの使い方は本HP関係なく必須科目です。勉強しましょう。

- ブログ記事は（少なくとも当面）はてなブログに掲載（当HPからリンク）。

## 新しいページの作り方

独自スクリプト`create.sh`で雛形を作ります。

あるいは、（その方が楽であれば）`template.md` を手動でコピーして編集してもかまいません。

### Pageを作成する場合

```bash
$ cd hp_management
$ bash create.sh (filename)
```

`(filename)` に生成するファイル名を指定します（内容に応じて適切に）。`./content/pages/`の下にメタデータ入りのMarkdownファイルが生成されます。

これに内容を追記します。htmlも認識されます。画像も `![]({attach}images/filename)` のように指定すると、 `./content/pages/images/filename` を読み込めます。

### Articleを作成する場合

``` bash
$ cd hp_management
$ bash create.sh (filename) (category) (tags)
```

Pageの場合と同じ。(category)を指定するとarticleになり、`content/articles/(schoolyear)/(category)/(filename).md`として作成されます。(category)はnewsかtechを推奨。(tags)は省略可能。なお(category)にpageまたはpagesを指定するとpageになります。

#### Jupyter Notebookの扱い

jupyter notebookに関してはmdファイルと同じ場所に入れて、メタデータのファイル作って、notebookの最初のセルにメタデータの内容と同じものを書けば大丈夫です。テストファイルがあるので参考にしてください。

## 導入した機能など

### Jupyter notebookをHTML出力できるようにする。

[この記事](https://qiita.com/driller/items/49a990cbdfb51afed620)に従って、インストール、pluginの導入、pelicanconf.pyの書き換えを行った。

### Themeの導入

[voidy-bootstrap](https://github.com/robulouski/voidy-bootstrap)というテーマを導入した。pelicanconf.pyの書き換えを行った。Twitterアカウントへのリンク設定などもpelicanconf.pyからできる。

### トップページの変更

デフォルトではarticlesのインデックスがトップページになる。これを変更する正式な方法が何かはよくわからないが、ひとまず `pelicanconf.py`中で

```python
INDEX_SAVE_AS = 'articles.html'
```

としてarticleのインデックスURLを`index.html`から`articles.html`に変更、次いで、pagesのうちの一つのメタデータで

```
Slug: index
```

と設定すると、このページがトップになる。この方法で `about.md` をトップに設定した。

## 使い方（管理者向け）

### 導入

#### ツールのインストール

環境：Python 3.6以降

```bash
$ pip install pelican Markdown
```

#### レポジトリのクローン

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management
$ git clone https://github.com/oumpy/oumpy.github.io.git ./output
$ git submodule update -i
$ cd themes && git submodule update voidy-bootstrap
```

2行目以降は  `./init.sh` で一括実行できる。

### 更新のアップ

更新を強制的に一括pushするスクリプトを用意しています。以下には書いていませんがGitHubの認証に関する設定も必要です。

#### [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へのpush

（スクリプト調整中。現在は非推奨です。）

```bash
$ cd hp_management
$ bash manage_push.sh
```
で、`./output/`以外の全てのファイルを[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へpushします。

#### [出力用レポジトリ](https://github.com/oumpy/oumpy.github.io)へのpush

こちらもGitの使い方としては乱暴ですが、サイト本体は真面目に履歴管理する対象ではない（ソースの方を管理すればOK）という思想に基づきます。

```bash
$ cd hp_management
$ bash blog_push.sh
```
を実行すれば、全てのファイルをhtmlにコンパイルして、`./output/`へ、そして[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)にプッシュします。


## 課題
### 管理システム

- 管理者以外も手を出しやすいよう、ディレクトリ構造やスクリプト類の整理。（かなり進展した。ひとまずOK？）

### サイト内容

- CSSの調整など。
- Home,Member,Calender,Contact, Activityのページを作る。

### 記事(article)類

- 本サイトで完結できるかの検討研究
- 写真ファイルのパス指定や内部リンク調整

- 手打ちをした部分（Author, Category, Tagsなど）のチェック

