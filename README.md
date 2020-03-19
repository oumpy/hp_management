# [阪医Python会HP](https://oumpy.github.io/)

原作：平岡悠  
現メンテナ：AtamaokaC

## 全体の仕組み

- ファイルが全て入った、[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)とoutputディレクトリのHTMLファイルのみが入った、[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)がある。

- 基本は管理用レポジトリをプルして、設定を変更したり、content内の新しいブログを追加したりする。完成したら、管理用レポジトリにプッシュするとともに、出力用レポジトリにもoutputの中身をpushする。この繰り返し。

- 会員は本レポジトリ（hp_management）を自分のGitHubアカウントにfork、`content/`下のMarkdownファイルを編集・追加してプルリクエストを送ることを推奨。

- 技術的説明は下に色々あるが、管理に携わる人以外は **「新しいページの作り方」以外は知らなくてもよい** 。（「導入」等もやらなくてよい。）

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

Pageの場合と同じコマンド。(category)を指定するとarticleになり、`content/articles/(schoolyear)/(category)/(filename).md`として作成されます。

- (category)はnewsかtech。(category)にpageまたはpagesを指定するとpageになる。
- (tags)は省略可能、ただしarticleの場合は編集の際に必ず入れること。（この仕様は変更するかも。）
- (category)は各ファイルのメタデータには記載されない。mdファイルを入れるディレクトリの名前がカテゴリとして認識される。
  （記載した場合はそちらが優先されるが、記載しないでください。）

#### 画像の設置と読み込み

画像については、`content/articles/(schoolyear)/(category)/images/(filename)_figs/(imagefile)`として保存し、`{attach}images/(filename)_figs/(imagefile)`で読み込むのを標準とします。

例えば、`sugoikiji.md`に画像ファイル `sugoigazou.png` を読み込みたい場合は `content/articles/2020sc/tech/images/sugoikiji_figs/sugoigazou.png` のように設置し、

```markdown
![Sugoi Gazou]({attach}images/sugoikiji_figs/sugoigazou.png)
```

と書けば読み込まれます。

#### Jupyter Notebookの扱い

jupyter notebookに関しては他の記事（mdファイル）と同じ場所に入れ、さらに同じ場所にメタデータファイル（`myarticle.ipynb`の場合は`myarticle.nbdata`）を置いてmdファイルと同様のメタデータを書きます。 `article/2018sy/tech_archive`の`lorentz.ipynb`および`lorentz.nbdata`を参考にしてください。

## 導入した機能など

### Jupyter notebookをHTML出力できるようにする。

プラグイン[pelican-ipynb](https://github.com/danielfrg/pelican-ipynb)を導入。pelicanconf.pyには以下の記述を追加。

```python
MARKUP = ('md', 'ipynb')
PLUGIN_PATH = './plugins'
PLUGINS = ['ipynb']
```

（初期導入時の参考記事：<https://qiita.com/driller/items/49a990cbdfb51afed620>）

### Themeの導入

[voidy-bootstrap](https://github.com/robulouski/voidy-bootstrap)というテーマを導入した。pelicanconf.pyの書き換えを行った。Twitterアカウントへのリンク設定などもpelicanconf.pyからできる。

さらにテーマ改変・ファイル追加によりLook & Feelの変更と機能追加を行っている。

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

#### レポジトリのクローンとテーマファイル(voidy-bootstrap)のコピー

```bash
$ cd anywhere_you_like
$ sh init.sh
```

テーマのファイルのみコピーし直したい場合は`sh init.sh -c`でOK。

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
- Recruit, Activity, Achievementsなどの内容作成。
- Member,Calenderなどは必要か？

### 記事(article)類

- 本サイトで完結できるかの検討研究
- 写真ファイルのパス指定や内部リンク調整
- 手打ちをした部分（Author, Tagsなど）のチェック
- 古い記事のいくつかはhtmlのコピペでMarkdownとしては酷い。

