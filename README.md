# [阪医Python会HP](https://oumpy.github.io/)

原作：平岡悠  
メンテナ (HP係)：AtamaokaC、takyamamoto

## 全体の仕組み

- 管理ツールおよびサイト原稿 (主にMarkdownファイル) が入った [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)と、`output/`ディレクトリ (主にHTMLファイル) が入った、[webサイト本体のレポジトリ (出力用レポジトリ)](https://github.com/oumpy/oumpy.github.io) がある。

- 基本は管理用レポジトリをプルして、設定を変更したり、`content/article/`内に新しい記事を追加したりする。
  完成したら、管理用レポジトリにプッシュするとともに、出力用レポジトリにもoutputの中身をpushする。この繰り返し。
   (後述の自動更新システムを使えば、日常更新作業で出力用レポジトリを手動で扱う必要はない。)

- 会員はアカウントを`oumpy`オーガニゼーションのメンバーに登録し、Markdown/Jupyter notebookファイルを`content/`下に追加してプルリクエストを送ることで、記事を投稿できる。

- 以下は主に管理者・開発者向けの技術的内容です。**記事投稿の方法などは`content/README.md`を参照** のこと。
以下の内容をやる必要はありません。

  (Markdownとgitの使い方は本HP関係なく必須科目です。勉強しましょう。)

## 導入した機能など

### 初期設定

初期導入時の参考記事：<https://qiita.com/driller/items/49a990cbdfb51afed620>

### Pluginの導入

プラグイン[pelican-plugins](https://github.com/getpelican/pelican-plugins)および[pelican-ipynb](https://github.com/danielfrg/pelican-ipynb)を導入。
pelicanconf.pyには以下の記述を追加。

```python
MARKUP = ('md', 'ipynb')
PLUGIN_PATHS = ['./plugins']
PLUGINS = ['pelican-ipynb.markup', 'render_math']
```

これでjupyter notebookファイル(.ipynb)とLaTeX数式の使用がそれぞれ可能になる。

### Themeの導入

[voidy-bootstrap](https://github.com/robulouski/voidy-bootstrap)というテーマを導入した。
pelicanconf.pyの書き換えを行った。
Twitterアカウントへのリンク設定などもpelicanconf.pyからできる。

さらにテーマ改変・ファイル追加によりLook & Feelの変更と機能追加を行っている。

### トップページの変更

デフォルトではarticlesのインデックスがトップページになる。
これを変更する正式な方法が何かはよくわからないが、ひとまず `pelicanconf.py`中で

```python
INDEX_SAVE_AS = 'articles.html'
```

としてarticleのインデックスURLを`index.html`から`articles.html`に変更、次いで、pagesのうちの一つのメタデータで

```
Slug: index
```

と設定すると、このページがトップになる。
この方法で `about.md` をトップに設定した。

### その他

- page/articleに関するファイル場所・URL等の基本設定。
- サイト固有文字列等の設定を `content/contentconf.py`, `content/contentpublishconf.py` に分離。

## 使い方（管理者向け）

### 導入

#### レポジトリのクローンとテーマファイル(voidy-bootstrap)のコピー

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management
$ sh tools/init.sh
```
(GitHub Pages のレポジトリは、 `content/contentpublishconf.py` の中で
```python
SITEREPOSITORY = 'https://github.com/oumpy/oumpy.github.io.git'
```
により定義。
`init.sh` の中で自動的に読み込まれる。）

テーマのファイルのみコピーし直したい場合も `sh init.sh` でOK。

#### Pythonパッケージのインストール

環境：Python 3.6以降

`hp_management/` 直下で

```bash
$ pip install -r requirements.txt
```

とすると、pelicanほか必要なパッケージを一括インストールできる。

### 更新のアップ

#### [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)

プルリクエストベースでのレポジトリ管理を行います。
- 投稿・開発の基本は `master` ブランチへのプルリクエスト。
  マージにはレビュー権限のある人の承認が必要。
- 記事投稿もシステム開発も基本は同じ。
  使い捨てブランチから `master` にプルリクを行う。
- 原則として記事投稿は幹部が、固定ページやニュースは正副代表が、システム設定はメンテナが、レビュー・承認して `master` にマージする。
  (`executives`, `hp_maintainers` などのチームと `CODEOWNERS` によりレビュー権限設定。)

#### [出力用レポジトリ](https://github.com/oumpy/oumpy.github.io)へのpush (手動で行う場合)

このレポジトリは出力にすぎないので、あまり真面目な変更履歴管理は行いません。
masterブランチに全て上書きしていく形でOKです。
そのために一括commit & pushするスクリプト`tools/pushsite.sh`を用意しています。
なおGitHubの認証に関する設定が事前に必要です。

```bash
$ cd hp_management
$ sh tools/pushsite.sh "Sugoi Kiji added."
```
を実行すれば、全てのファイルをhtmlにコンパイルして、`./output/`へ、そして[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)にプッシュします。
`tools/pushsite.sh` に引数として与えた文字列がコミットのコメントになります。
省略すると単に "Update" になります。

##### リモートブランチ（ソース）のpullと出力のpush

```bash
$ sh tools/updatesite.sh "Yabai update"
```

で、masterブランチをpull & checkout、コンパイルして出力用レポジトリに指定したコメント付きでcommit & pushします。
第2引数としてmaster以外のソースブランチ、第3引数としてmaster以外のターゲットブランチを指定することも可能。

このスクリプト `tools/updatesite.sh` は主に自動push用に用意されていますが、動作を理解していれば手動で用いても問題ありません。

## 自動更新 & プレビュー機構
### GitHub Actions (デフォルト)
#### 自動更新
本レポジトリのmasterブランチにプルリクエストがマージされると、GitHubのシステム上でソースが自動的にコンパイルされ、サイトにプッシュされます。
この指示の実体 (ワークフロー) は `.github/workflow/deploy_on_site.yml` です。

#### プレビュー自動作成
masterブランチ以外のブランチ (仮に `article/sugoikiji` とする) を更新 (プッシュ) すると、やはりGitHubのシステム上でソースが自動的ににコンパイルされて
https://oumpy.github.io/previews/refs/heads/article/sugoikiji/ にアップされます。
ブランチを削除するとプレビューも削除されます。
この指示の実体 (ワークフロー) は `.github/workflow/preview.yml` / `.github/workflow/delete_preview.yml` です。

### Webhook と外部サーバを用いた機構
webサーバ上にローカルレポジトリを設置することで、GitHubのwebhook機能を用いてサイトを自動更新/プレビューすることができます。
ただしプレビューは任意のブランチのプッシュで任意のコードを実行できるので、他用途と兼用しているサーバ上での使用は非推奨です。

#### 本機能を使用する方法

1. webサーバ上の適当なディレクトリ（webからアクセスできないところ）にhp_managementを正しく設置。

2. webhookを受け取るcgiの場所とファイル名を好きなように決める。
   （ローカルパスを以下仮に`/home/hoge/www/cgi-bin/deploy.cgi`、対応するURLを http://www.example.io/cgi-bin/deploy.cgi とする。）

3. hp_management/webhookconf.py`を以下のような内容で作成する。

```python
# webhookconf.py
# cgipath, secret の2つの変数を以下のように設定する。
# 上で決めたcgiのフルパスを記述 (webサーバの設定に依存)
cgipath = '/home/hoge/www/cgi-bin/deploy.cgi'
# githubとの間でのパスワードとなるバイト列を設定
secret = b'xyzabc....'
```

4. `sh tools/init.sh` を実行。設定したパスにcgiが設置される。

5. GitHubの本レポジトリでwebhookを設定する。
   cgiのURL (http://www.example.io/cgi-bin/deploy.cgi) とsecretを設定し、`content type` に `application/json`、また "Just the push event." を選択。
   自動更新のみの場合は対象を `master` ブランチ（推奨）、プレビューも使う場合は全てのブランチにする。

以上により、本レポジトリ (hp_management) のブランチ更新を自動検出して`tools/updatesite.sh` が実行されます。
