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

## 導入した機能など

### 初期設定

初期導入時の参考記事：<https://qiita.com/driller/items/49a990cbdfb51afed620>

### Pluginの導入

プラグイン[pelican-plugins](https://github.com/getpelican/pelican-plugins)および[pelican-ipynb](https://github.com/danielfrg/pelican-ipynb)を導入。pelicanconf.pyには以下の記述を追加。

```python
MARKUP = ('md', 'ipynb')
PLUGIN_PATHS = ['./plugins']
PLUGINS = ['pelican-ipynb.markup', 'render_math']
```

これでjupyter notebookファイル(.ipynb)とLaTeX数式の使用がそれぞれ可能になる。

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
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management
$ sh tools/init.sh
```

テーマのファイルのみコピーし直したい場合は`sh init.sh -c`でOK。

### 更新のアップ

#### [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)

通常のレポジトリ管理を行います。

- サイト内容(`content`ディレクトリ内)：`content`ブランチへのプルリクエストを受け付けます。
- それ以外：開発参加者間で適切に管理します。

#### [出力用レポジトリ](https://github.com/oumpy/oumpy.github.io)へのpush

このレポジトリは出力にすぎないので、あまり真面目な変更履歴管理は行いません。
masterブランチに全て上書きしていく形でOKです。
そのために一括commit & pushするスクリプト`tools/pushblog.sh`を用意しています。
なおGitHubの認証に関する設定が事前に必要です。

```bash
$ cd hp_management
$ sh tools/pushblog.sh "Sugoi Kiji added."
```
を実行すれば、全てのファイルをhtmlにコンパイルして、`./output/`へ、そして[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)にプッシュします。
`tools/pushblog.sh` に引数として与えた文字列がコミットのコメントになります。
省略すると単に "Update" になります。

##### リモートブランチ（ソース）のpullと出力のpush

```bash
$ sh tools/updateblog.sh "Yabai update"
```

で、masterブランチをpull & checkout、コンパイルして出力用レポジトリに指定したコメント付きでcommit & pushします。第2引数としてmaster以外のソースブランチ、第3引数としてmaster以外のターゲットブランチを指定することも可能。

このスクリプト `tools/updateblog.sh` は主に自動push用に用意されていますが、動作を理解していれば手動で用いても問題ありません。

## 課題
### 管理システム

- 管理者以外も手を出しやすいよう、ディレクトリ構造やスクリプト類の整理。（かなり進展した。ひとまずOK？）

### サイト内容

- CSSの調整など。
- Member,Calenderなどは必要か？

### 記事(article)類

- 本サイトで完結できるかの検討研究
