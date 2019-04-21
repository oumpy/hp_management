## [阪医Python会HP](https://oumpy.github.io/)
ファイルが全て入った、[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)とoutputディレクトリのHTMLファイルのみが入った、[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)がある。

基本は管理用レポジトリをプルして、設定を変更したり、content内の新しいブログを追加したりする。完成したら、管理用レポジトリにプッシュするとともに、出力用レポジトリにもoutputの中身をpushする。この繰り返し。

## Jupyter notebookをHTML出力できるようにする。
[この記事](https://qiita.com/driller/items/49a990cbdfb51afed620)に従って、インストール、pluginの導入、pelicanconf.pyの書き換えを行った。

## Themeの導入
[MinimalXY](https://github.com/petrnohejl/MinimalXY/tree/87f0ebb57543b7810dffc9ebe05ed96bc897ffd1)という、シンプルなテーマを導入した。pelicanconf.pyの書き換えを行った。Twitterアカウントへのリンク設定などもpelicanconf.pyからできる。

## 新しいブログの作り方
私が作ったスクリプト使っていくので、我流でも大丈夫です。
```
bash create.sh arg1 arg2
```
でarg1.mdという名前でarg2というカテゴリのマークダウン（必要なメタデータ付き） 生成して、カテゴリー名のディレクトリに送ります。
content内にマークダウンが入っているので、ブログを書きます。htmlも認識されます。imagesフォルダ作ったので、そこに写真おいてパス通せば表示できます。

## Jupyter Notebookの扱い
jupyter notebookに関してはcontent/notebookに入れて、メタデータのファイル作って、notebookの最初のセルにメタデータの内容と同じものを書けば大丈夫です。テストファイルがあるので参考にしてください。

## [ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へのpush
```
bash manage_push.sh
```
で、全てのファイルを[ブログ管理用のレポジトリ](https://github.com/oumpy/hp_management)へpushします。

## [出力用レポジトリ](https://github.com/oumpy/oumpy.github.io)へのpush
```
bash blog_push.sh
```
を実行すれば、全てのファイルをhtmlにコンパイルして、/outputへ、そして[HTMLをユーザーに出力するレポジトリ](https://github.com/oumpy/oumpy.github.io)にプッシュします。


## 課題
まずは、wordpressと同じ状況に近づける。4月中に完成させて、5月に完全移行するのが目標。
・Home,Member,Calender,Contact, Activityのページを作る。
・写真ファイルのパス指定や内部リンク調整
・手打ちをした部分（Author, Category, Tagsなど）のチェック
