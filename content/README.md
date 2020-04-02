# サイト内容（contentディレクトリ）の内容と書き方

サイトのソースはこのディレクトリ (`content`) に入っている。

```bash
$ cd anywhere_you_like
$ git clone https://github.com/oumpy/hp_management.git
$ cd hp_management/content
```

## サイト内容の更新

### 新しいページの作り方

`content` ディレクトリに移動し、独自スクリプト `create.sh` で雛形を作ります。

あるいは、 (その方が楽であれば) `content/template.md` を手動でコピーして編集してもかまいません。

####  Pageを作成する場合

```bash
$ cd hp_management/content
$ sh create.sh (filename)
```

`(filename)` に生成するファイル名を指定します (内容に応じて適切に)。
`content/pages/` の下にメタデータ入りのMarkdownファイルが生成されます。

これに内容を追記します。htmlも認識されます。

#### Articleを作成する場合

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

## サイト全体の情報設定

`content/contentconf.py` に、サイト表題・副題やリンク先などの情報が書かれている。**この設定は `pelicanconf.py` の一部であり、また `pelicanconf.py` 本体の記述よりも優先される。**

同様に `content/contentpublishconf.py` もあり、テスト時にはない方がよい設定を記述する。

## ToDo

- 写真ファイルのパス指定や内部リンク調整
- 手打ちをした部分（Author, Tagsなど）のチェック
- 古い記事のいくつかはhtmlのコピペでMarkdownとしては酷い。
- Member,Calenderなどは必要か？
