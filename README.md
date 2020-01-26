# [阪医Python会HP](https://oumpy.github.io/)

ファイルが全て入った、[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)とoutputディレクトリのHTMLファイルのみが入った、[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)がある。

基本は管理用レポジトリをプルして、設定を変更したり、content内の新しいブログを追加したりする。完成したら、管理用レポジトリにプッシュするとともに、出力用レポジトリにもoutputの中身をpushする。この繰り返し。

## 導入した機能など

### Jupyter notebookをHTML出力できるようにする。

[この記事](https://qiita.com/driller/items/49a990cbdfb51afed620)に従って、インストール、pluginの導入、pelicanconf.pyの書き換えを行った。

### Themeの導入

[MinimalXY](https://github.com/petrnohejl/MinimalXY/tree/87f0ebb57543b7810dffc9ebe05ed96bc897ffd1)という、シンプルなテーマを導入した。pelicanconf.pyの書き換えを行った。Twitterアカウントへのリンク設定などもpelicanconf.pyからできる。

## 使い方

### 導入

#### ツールのインストール

環境：Python 3.6以降

```bash
$ pip install pelican ghp-import Markdown
```

#### レポジトリのクローン

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ git clone https://github.com/oumpy/oumpy.github.io.git ./output/
```

### 新しいページの作り方

独自スクリプトを使っていくので、我流でも大丈夫です。
```bash
$ cd hp_management
$ bash create.sh arg1 arg2
```
で`arg1.md`という名前で`arg2`というカテゴリのマークダウン（必要なメタデータ付き） 生成して、カテゴリー名のディレクトリに送ります。
`./content/`内にマークダウンが入っているので、内容を書きます。htmlも認識されます。画像も適当なフォルダを作ってその中に置き、パスを通せば表示できます。

#### Jupyter Notebookの扱い

jupyter notebookに関しては`./content/notebook/`に入れて、メタデータのファイル作って、notebookの最初のセルにメタデータの内容と同じものを書けば大丈夫です。テストファイルがあるので参考にしてください。

### 更新のアップ（簡易・非推奨）

更新を強制的に一括pushするスクリプトを用意しています。Gitの使い方としては乱暴・危険で、非推奨です。少しのミスで容易に事故が起きるので、**何をしているか理解**していて、**一人で管理**していて、かつ**とにかく何でもいいからpushしたい**、とき用。できるだけgitを通常の手順でちゃんと使いましょう。

#### [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へのpush

```bash
$ cd hp_management
$ bash manage_push.sh
```
で、`./output/`以外の全てのファイルを[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へpushします。

#### [出力用レポジトリ](https://github.com/oumpy/oumpy.github.io)へのpush

```bash
$ cd hp_management
$ bash blog_push.sh
```
を実行すれば、全てのファイルをhtmlにコンパイルして、`./output/`へ、そして[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)にプッシュします。


## 課題
まずは、wordpressと同じ状況に近づける。4月中に完成させて、5月に完全移行するのが目標。

・Home,Member,Calender,Contact, Activityのページを作る。

・写真ファイルのパス指定や内部リンク調整

・手打ちをした部分（Author, Category, Tagsなど）のチェック

---平岡悠---