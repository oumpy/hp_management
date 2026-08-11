# oumpy.github.io システム技術仕様書

大阪大学医学部Python会Webサイト <https://oumpy.github.io> の生成・配信システム
(レポジトリ [`oumpy/hp_management`](https://github.com/oumpy/hp_management)) の技術仕様。

- 対象読者: HP係・システムメンテナ・将来の開発者
- 最終更新: 2026-08-11 (システム近代化アップデート時)
- 記事投稿の手順は `content/README.md` を、運用手順の概要はトップの `README.md` を参照。
  本書は「システムがどう作られているか」を記述する。

---

## 1. 全体アーキテクチャ

```
[hp_management レポジトリ]                | [oumpy.github.io レポジトリ]
  content/   : 記事・固定ページ原稿        |   master ブランチ = 公開HTML
  (外側)     : ビルドシステム一式           |   previews/ : ブランチプレビュー
       |                                  |        ^
       |  Pelican (静的サイトジェネレータ)  |        |
       +--> output/ --- git push ---------+--------+
```

- **ソースレポジトリ** `oumpy/hp_management`: Markdown / Jupyter notebook の原稿と、
  Pelican ベースのビルドシステム。プルリクエストベースで運用。
- **出力レポジトリ** `oumpy/oumpy.github.io`: ビルド成果物 (HTML) のみを置く
  GitHub Pages レポジトリ。履歴は「上書きコミットの積み重ね」であり管理対象ではない。
  URL は `content/contentpublishconf.py` の `SITEREPOSITORY` で定義。
- master へのマージ → GitHub Actions が自動ビルド → 出力レポジトリへ push、が基本経路 (§9)。

### 使用技術スタック (2026-08 現在)

| 層 | 技術 | バージョン方針 |
|---|---|---|
| ジェネレータ | Pelican | >= 4.11 (CI では常に最新) |
| Python | CPython | 3.9+ (CI は 3.13) |
| Markdown 変換 | Python-Markdown | >= 3.6 |
| ipynb 変換 | 自作 `myplugins/ipynb_reader` | nbconvert 等の Jupyter スタックに **非依存** |
| CSS フレームワーク | Bootstrap | 5.3.7 (jsDelivr CDN, SRI付き) |
| アイコン | Font Awesome | 6.7.2 (jsDelivr CDN) |
| 数式 | MathJax | 2.7.3 (pelican-render-math が注入) |
| JS | Vanilla JS + Bootstrap bundle | **jQuery 非依存** |

---

## 2. ディレクトリ構成

```
hp_management/
├── pelicanconf.py            # メイン設定 (開発ビルド用)
├── publishconf.py            # 本番ビルド用設定 (pelicanconf を継承)
├── requirements.txt          # Python 依存パッケージ
├── .github/workflows/        # CI/CD (§9)
├── content/                  # ★サイトの内容 (原稿・コンテンツ側設定)
│   ├── contentconf.py        #   サイト固有設定 (メニュー・SNS・文言など)
│   ├── contentpublishconf.py #   本番のみの設定 (SITEURL, SITEREPOSITORY)
│   ├── articles/<年度>sy/{blog,news}/   # 記事 (.md / .ipynb + .nbdata)
│   ├── pages/                #   固定ページ
│   └── create.py             #   記事スケルトン生成スクリプト
├── myplugins/                # ★自作/vendor プラグイン (§6)
├── theme/voidy-bootstrap/    # ★テーマ一式 (完全 vendor 済み, §7)
├── tools/                    # 初期化・デプロイスクリプト (§9, §10)
└── 3rdtools/misc-tools/      # (サブモジュール) レガシー webhook 用ロックツール
```

サブモジュールは `3rdtools/misc-tools` のみ (レガシー機構専用。通常のビルドには不要)。
かつて存在した `themes/` (pelican-themes) と `plugins/` (pelican-plugins) の
巨大サブモジュールは 2026-08 に廃止し、必要ファイルをレポジトリ内へ vendor した。

---

## 3. 設定ファイルの構造と読み込み順

「システム設定」と「コンテンツ設定」を分離するのが本システムの特徴。

```
pelicanconf.py                     (システム側・開発ビルド)
  └─ from content.contentconf import *      (コンテンツ側)
publishconf.py                     (本番ビルド; pelican -s publishconf.py で使用)
  ├─ from pelicanconf import *
  ├─ PLUGINS += [sitemap]; robots.txt 追加
  └─ from content.contentpublishconf import *
```

- `tools/lib/pelicanns.py` は**空のモジュール**で、名前空間の受け渡しに使う:
  `pelicanconf.py` が読み込み済みグローバルを `pelicanns` に代入してから
  `contentconf.py` を import し、`contentconf.py` 側は
  `from tools.lib.pelicanns import *` でシステム側の定義 (`datetime`,
  `PLUGIN_PATHS` 等) を参照する。**削除してはならない。**
- コンテンツ側 (`contentconf.py`) にはメニュー定義 `ADD_ON_MENU`、SNS リンク
  `SOCIAL`、Google CSE ID、タググループ、サイドバー文言などを置く。
- プレビュービルドでは CI が `contentconf.py` の末尾に `SITEURL` 等を追記する (§9.2)。

### 主要な Pelican 設定 (pelicanconf.py)

- `MARKUP = ['md', 'ipynb']`
- `ARTICLE_SAVE_AS = ARTICLE_URL = '{category}/{date:%Y}/{date:%m}/{slug}.html'`
- `PAGE_SAVE_AS = '{slug}.html'`, `CATEGORY_SAVE_AS = '{slug}.html'`
- `INDEX_SAVE_AS = 'articles.html'` — トップページは `about.md` (`Slug: index`) が担う
- `SLUGIFY_SOURCE = 'basename'` — スラグはファイル名由来
- ページネーション: 2ページ目以降は `{base_name}/latests/{number}/`
- フィード: `feeds/all.atom.xml`, `feeds/all.rss.xml`
- `ARTICLE_EXCLUDES_DIRNAMES = PAGE_EXCLUDES_DIRNAMES = ['attach', 'images']`
  (自作プラグイン excludes_dirnames が処理)

---

## 4. 依存パッケージ (requirements.txt)

```
pelican[markdown] >= 4.11
Markdown >= 3.6
beautifulsoup4                      # autosummary が使用
pelican-sitemap                     # 本番のみ有効
pelican-render-math                 # 数式 (MathJax)
pelican-tag-cloud
pelican-related-posts
pelican-simple-footnotes
pelican-neighbors
minchin.pelican.plugins.nojekyll    # .nojekyll 生成
```

すべて PyPI の現行 namespace plugin。バージョン上限固定なし
(旧システムの `nbconvert<6`, `jinja2<3.1`, `pelican-jupyter`, `ipython_genutils`,
`libsass`, `requests` は 2026-08 に廃止)。

---

## 5. ビルドパイプライン

### 5.1 コマンド

| 操作 | コマンド |
|---|---|
| 開発ビルド | `pelican` (pelicanconf.py を自動読込) |
| ローカルサーバ | `pelican -l [-p 8000]`、自動再生成付きは `pelican -r -l` |
| 本番ビルド | `pelican -s publishconf.py` |
| 出力先変更 | `pelican -o output.new` |

Makefile / tasks.py は 2026-08 に廃止し、Pelican 公式 CLI を直接使う
(現行 Pelican の標準的な使い方。旧 `make html`/`make publish` 相当は上表)。

### 5.2 処理の流れ

1. Pelican が `content/` を走査。リーダーで `.md` (標準) / `.ipynb` (自作) を HTML 化
2. 各プラグイン (§6) が signal 経由で内容を加工
3. テーマ (§7) の Jinja2 テンプレートでページ生成

### 5.3 数式 (MathJax) の扱い

- pelican-render-math が Markdown 拡張として `$...$` / `$$...$$` を保護し
  `<span class="math">\(...\)</span>` に変換、数式を含むページに MathJax 2.7.3
  (cdnjs) のローダスクリプトを注入する。
- ipynb 記事の Markdown セルも同じ Markdown 設定 (`MARKDOWN` セッティング) で
  変換されるため、同一の機構で数式が処理される。
- 注入スクリプトは `id` ガード付きのため重複実行はされない。

(かつて存在した `content/postprocess.sh` による旧 URL 互換シンボリックリンク
生成は 2026-08 に廃止。汎用の `postprocess` プラグインも使用箇所がなくなった
ため撤去した。必要なら git 履歴から復元できる。)

---

## 6. プラグイン (myplugins/)

`PLUGIN_PATHS = ['./myplugins']`。各プラグインは Pelican の signal API で接続する。

### 6.1 ipynb_reader (自作, 2026-08 新規)

Jupyter notebook 記事のリーダー。**依存は markdown + pygments のみ**。

- **登録**: `readers_init` signal で `readers.reader_classes['ipynb']` を差し替え。
- **メタデータ**: 同名の `.nbdata` ファイル (必須) から `Key: Value` 形式で読む。
  - 例: `Title:`, `Date:`, `Modified:`, `Tags:`, `Author:`, `Summary:`
  - `Subcells: [first, last]` — 表示セル範囲のスライス指定 (pelican-jupyter 互換。
    `ast.literal_eval` で解釈し `cells[first:last]`)。
  - `.nbdata` が無い notebook はエラー (記事として扱わない)。
- **対応 nbformat**: v4 のみ (全既存記事が v4)。
- **セル変換**:
  - markdown セル → サイト共通の `MARKDOWN` 設定で `markdown.Markdown` により変換
    (render-math 拡張が効く)。`attachment:名前` は data URI に置換。
  - code セル → Pygments (`python3` レクサ、他言語は kernelspec から判定) で
    ハイライト。フォーマッタの cssclass は `highlight hl-ipython3`。
  - raw セル → エスケープして `<pre class="raw_cell">`。
- **出力 (outputs) の変換**: MIME 優先順位
  `text/html > image/svg+xml > image/png > image/jpeg > text/latex > text/markdown > text/plain`。
  - stream → `<pre>` (`output_stdout` / `output_stderr`)、ANSI エスケープは除去
  - error → traceback を `<pre class="... output_error">`
  - 画像 → base64 data URI の `<img>`
  - text/latex → `<div class="... output_latex math">` (MathJax に委ねる)
- **生成 HTML 構造**: 旧 nbconvert 5 "basic" テンプレート互換のクラス名を維持:

  ```html
  <div class="jupyter-notebook">
    <div class="cell border-box-sizing code_cell rendered">
      <div class="input">
        <div class="prompt input_prompt">In [1]:</div>
        <div class="inner_cell"><div class="input_area">
          <div class="highlight hl-ipython3"><pre>...</pre></div>
        </div></div>
      </div>
      <div class="output_wrapper"><div class="output">
        <div class="output_area">
          <div class="prompt output_prompt">Out[1]:</div>
          <div class="output_subarea ...">...</div>
        </div>
      </div></div>
    </div>
  </div>
  ```

  この互換性は (a) `voidybootstrap-custom.css` の notebook 用スタイル、
  (b) Colab ボタンスクリプトの記事判定 (`highlight hl-ipython3` の存在で
  ipynb 記事と判定) のために**維持必須**。
- **レイアウト CSS**: `voidybootstrap-custom.css` 末尾の
  `.jupyter-notebook ...` ブロックが In/Out プロンプトを左側に置く flex
  レイアウトを提供 (768px 未満では縦積み)。

### 6.2 その他の自作プラグイン

| プラグイン | 機能 |
|---|---|
| `autosummary` | 記事冒頭から自動で要約生成 (bs4 使用)。`summary` と併存中 (要整理) |
| `summary` (vendor) | `<!-- PELICAN_BEGIN_SUMMARY -->` マーカーによる要約。旧 pelican-plugins 由来 |
| `shortcodes` (vendor) | `SHORTCODES` 設定によるショートコード展開。`youtube`, `embed` を定義済み |
| `category_names` | カテゴリ表示名の差し替え (`CATEGORYNAMES_ALTERNATIVES`: blog→技術ブログ 等) |
| `apply_jinja2` | 文字列を Jinja2 評価する `apply_jinja2` フィルタ。また `Jinja2: True` メタデータを持つ記事/ページの本文・タイトルを Jinja2 として描画 (例: `search.md` が `GOOGLE_CSE_ID` を埋め込むのに使用) |
| `path2obj` | URL からページ/記事オブジェクトを引く `url2obj` フィルタ |
| `subsections` | ページ階層 (サブセクション) 情報を構築。メニュー AUTO 展開に使用 |
| `makemenu` | 多階層ナビゲーションメニュー HTML の生成 (§7.3) |
| `excludes_dirnames` | `*_EXCLUDES_DIRNAMES` 設定でディレクトリ名単位の除外 |
| `skiptags` | 特定タグの除外処理 |

vendor プラグイン (`summary`, `shortcodes`) の由来はアーカイブ済み
[pelican-plugins](https://github.com/getpelican/pelican-plugins) commit `18a59c3`。

---

## 7. テーマ (theme/voidy-bootstrap)

[voidy-bootstrap](https://github.com/robulouski/voidy-bootstrap) (MIT, commit
`83f4d80`) を基に大幅カスタマイズしたもの。**上流ファイルも含め全て vendor 済み**
で、ビルド時のコピー処理はない (旧方式: `cp -an` で上流とマージしていた)。

### 7.1 CSS 構成 (読み込み順)

1. **Bootstrap 5.3.7** — jsDelivr CDN、SRI 付き (base.html にハードコード)
2. **Font Awesome 6.7.2** — jsDelivr CDN (`FONT_AWESOME_LINK` 設定で変更可能。
   `integrity` キーは省略可)
3. `theme/css/pygment.css` — コードハイライト基本色
4. `theme/css/theme-overrides.css` — **Bootstrap カスタマイズ層**。
   CSS 変数でプライマリ色 `#0f5889` 等を上書き。旧 SASS パイプライン
   (theme.scss + libsass + Bootstrap ソース DL) の代替
5. `theme/css/voidybootstrap-custom.css` — サイト固有スタイル本体
   (ナビバー配色 `#1d2113`、ジャンボトロン、見出し装飾、サイドバー box、
   notebook セル、シンタックスハイライト Monokai 風配色、寄付者フォント
   サイズ `.donator-E*` など)

Bootstrap のバージョンを上げる場合は `base.html` の `<link>`/`<script>` の
URL と SRI ハッシュを両方更新すること。

### 7.2 テンプレート構成

- `templates/base.html` — 全体骨格。ナビバー (BS5 `navbar-expand-lg sticky-top`)、
  CDN 読み込み、フッタ、`voidy-menu.js` の読み込み。
- 上書きテンプレート: `index.html`, `page.html`, `archives.html`, `category.html`,
  `tag.html`, `author(s).html` 等。
- `templates/includes/custom/` — 本サイト独自の部品群。主要なもの:

| ファイル | 役割 |
|---|---|
| `sidebar.html`, `sb_*.html` | サイドバー (検索, SNS, 新着, タグ, 支援, X timeline) |
| `open_in_colab_scripts.html` | ipynb 記事に GitHub/Colab バッジを自動付加 (§7.4) |
| `toc_scripts.html` | h2/h3 から目次を生成 (vanilla JS) |
| `scrolltop_scripts.html` | ページトップへ戻るボタン |
| `article_header_info.html` | 記事ヘッダ (日付・著者・タグバッジ shields.io) |
| `*_showmodified_scripts.html` | `Modified:` メタデータの表示 |
| `utterances.html` | コメント欄 (utteranc.es, issue-term: title) |
| `footer.html` | 著作権表示フッタ |

- `CUSTOM_*` 設定 (pelicanconf.py) でどの include を使うかを制御する
  (voidy-bootstrap のフック機構)。

### 7.3 ナビゲーションメニュー

3つの構成要素からなる:

1. **`contentconf.py` の `ADD_ON_MENU`** — `MenuItem(url, title=, subsections=,
   active_pages=, self_in_subsections=)` のリストでメニュー構造を宣言。
   - `subsections=MenuItem.AUTO` → subsections プラグインの情報から自動展開
   - `self_in_subsections=True` → サブメニュー先頭に親ページ自身+区切り線を挿入
     (モバイルで親ページへ到達するための導線。**重要**)
   - `active_pages` — 正規表現。現在ページがマッチしたら親項目を active 表示
2. **`makemenu` プラグイン** — `resolve_menu` Jinja2 フィルタを提供。
   `ADD_ON_MENU` を `{url, title, active, divider, children}` の素朴な dict 木に
   解決するだけで、**HTML は一切生成しない**。active 判定 (URL 一致 /
   `active_pages` 正規表現 / 子の active の伝播) もここで行う。
   トラバーサルは原実装のスタイルを保った非再帰 DFS (明示スタック +
   FORWARD/BACK の2フェーズ訪問)。
   マークアップはテーマ側の再帰マクロ `templates/includes/menu.html` が担当
   (トップ階層: `nav-item dropdown`, 下位階層: `dropdown dropend`。
   `data-bs-toggle` は**付けない** — Bootstrap の Dropdown JS を意図的に使わない)。
   相対 URL は Pelican が per-page に相対化する `{{ SITEURL }}` に任せる
   (かつての rooturl 自前計算は廃止)。
3. **`static/js/voidy-menu.js` + CSS** — 動作定義:
   - デスクトップ (>=992px): CSS `:hover` でサブメニュー表示
     (`#main-navbar .dropdown:hover > .dropdown-menu`)。親クリックはリンク遷移。
   - モバイル (<992px): ハンバーガー (BS5 collapse) 内で、親タップは遷移せず
     サブメニューをインライン展開 (`.show` 切替、aria-expanded 更新)。
     ページ遷移はサブメニュー内の親自身項目から。
   - ブレークポイント 992px は `navbar-expand-lg` と一致させてある。変更時は
     CSS メディアクエリと voidy-menu.js の `DESKTOP` 定数を両方変更のこと。

### 7.4 Colab / GitHub バッジ

`open_in_colab_scripts.html` (vanilla JS):

- ページ内に `highlight hl-ipython3` があれば ipynb 記事と判定し、
  meta タグ (`article:date`, `article:category:slug`, `article:slug`,
  `source-repository`) から
  `content/articles/<年度>sy/<カテゴリ>/<slug>.ipynb` のパスを組み立てて
  GitHub / Colab へのバッジリンクを `#colablink` に挿入する。
  - 年度は「4月始まり」: 記事月が 1〜3 月なら年-1 を年度とする。
- また、本文中の `.ipynb` へのローカルリンクすべての直後に Colab バッジを付加。

### 7.5 X (Twitter) タイムライン

`sb_twittertl.html`:

- 公式 widgets.js (`platform.twitter.com/widgets.js`) による timeline 埋め込みを
  まず試行 (`data-dnt`, `data-chrome="noheader nofooter noborders"`)。
- **フォールバック**: 6 秒後に iframe が描画されていなければ、box の中身を
  `@oumed_python をフォロー` ボタン (リンク) に置換する。
  X のタイムライン埋め込みは 2023 年以降未ログイン閲覧者にはほぼ表示されないため。
- 設定: `TWITTER_TIMELINE_URL`, `TWITTER_USERNAME`, `TWITTER_TIMELINE_HEIGHT`,
  `CUSTOM_TWITTERTL_TITLE` (contentconf.py / pelicanconf.py)。

### 7.6 外部サービス一覧

| サービス | 用途 | 場所 |
|---|---|---|
| jsDelivr | Bootstrap / Font Awesome 配信 | base.html, `FONT_AWESOME_LINK` |
| cdnjs | MathJax 2.7.3 | render-math が注入 |
| Google PSE | サイト内検索 (`GOOGLE_CSE_ID`)。検索窓はサイドバーのフォームが `/search.html?q=...` へ送信し、`search.md` (Jinja2: True) の `gcse-searchresults-only` 要素が現行の cse.js 埋め込みで結果を表示 | `sb_google_cse.html`, `search.md` |
| utteranc.es | 記事コメント欄 (GitHub Issues 連携) | `utterances.html` |
| shields.io | タグ/GitHub バッジ画像 | `article_header_info.html` ほか |
| platform.twitter.com | X タイムライン・共有ボタン | §7.5, `sharing_scripts.html` |
| connect.facebook.net | いいねボタン | `sharing_scripts.html` |
| colab.research.google.com | ノートブック実行リンク | §7.4 |

---

## 8. URL 設計

- 記事: `/{category}/{yyyy}/{mm}/{slug}.html` (例 `/blog/2018/10/bayesian_ttest.html`)
- 固定ページ: `/{slug}.html`。トップは `about.md` (Slug: index) → `/index.html`
- カテゴリ: `/{slug}.html` (blog.html, news.html) — 固定ページと同じ階層に置く設計
- 記事インデックス: `/articles.html`、ページネーション `/articles/latests/N/`
- タグ: `/tag/<tag>.html`、著者: `/author/<name>.html`
- フィード: `/feeds/all.atom.xml`, `/feeds/all.rss.xml`

---

## 9. CI/CD (GitHub Actions)

Secrets: `bot_identity` (デプロイ用 SSH 秘密鍵), `known_hosts`。
出力レポジトリへの push は SSH config の `Host github` + `url.github:.insteadOf`
書き換えで行う。

### 9.1 deploy_on_site.yml — 本番デプロイ

- トリガ: `master` への push
- 手順: checkout (submodules 含む) → Python 3.13 + pip cache →
  `pip install -r requirements.txt` → `tools/init.sh` →
  `pelican -s publishconf.py -o output.new` → 旧 output の `.git` と `previews/` を
  output.new に移植 → 鮮度チェック (master が最新か) → commit & push
- ビルド失敗はワークフロー失敗になる (`|| true` は 2026-08 に除去)

### 9.2 preview.yml — ブランチプレビュー

- トリガ: master 以外のブランチ push
- build job: `contentconf.py` 末尾に `SITEURL = '/previews/refs/heads/<branch>'`
  等を追記してから `pelican` を実行。成果物を artifact として保存
- push job: 出力レポジトリを clone し、`previews/refs/heads/<branch>/` に配置して
  push。プレビュー URL: `https://oumpy.github.io/previews/refs/heads/<branch>/`
- 鮮度チェック: ビルドしたコミットがまだブランチ先端のときのみ push

### 9.3 delete_preview.yml / post_preview_link.yml

- ブランチ削除時に対応するプレビューを削除
- PR 作成時にプレビュー URL をコメントとして自動投稿 (同一レポジトリの PR のみ)

### 9.4 label.yml

`.github/labeler.yml` に基づく PR 自動ラベル付け。

---

## 10. レガシー: webhook 自動更新機構

GitHub Actions 導入以前の自前サーバ用機構。現在は未使用だが残置してある。

- `tools/lib/deploy.py` — GitHub webhook を受ける CGI のテンプレート
  (`tools/setcgi.sh` が `webhookconf.py` の設定でインスタンス化)
- `tools/updatesite.sh` — pull → ビルド → 出力 push を一括実行。
  `3rdtools/misc-tools` の `pexlock`/`punlock` で排他制御
- `tools/pushsite.sh` — output/ を出力レポジトリへ一括 commit & push
- 完全撤去する場合: `tools/{setcgi.sh,updatesite.sh,gotobranch.sh}`,
  `tools/lib/deploy.py`, `webhook_requirements.txt`, `3rdtools/` サブモジュール,
  README の該当節を削除すればよい (Actions 経路には影響しない)

---

## 11. 2026-08 近代化アップデートの設計判断

背景と決定の記録 (2026-08 の一連のシステム更新コミット群)。

1. **サブモジュール廃止 / vendor 化** — pelican-themes / pelican-plugins は
   どちらもアーカイブ済み巨大モノレポで、クローンコストと供給リスクが大きい。
   実使用分 (テーマ1つ、プラグイン2つ) のみ取り込み。
2. **pelican-jupyter → 自作リーダー** — pelican-jupyter は開発停止で
   `nbconvert<6` を要求し、Jinja2 まで古い版に固定される諸悪の根源だった。
   ノートブックの「静的表示」に必要な機能は小さいため、markdown+pygments のみで
   再実装し、Jupyter スタック依存をゼロにした。旧 CSS/JS との互換のため
   nbconvert 5 のクラス構造を踏襲 (§6.1)。
3. **SASS パイプライン廃止** — 旧方式はビルドごとに Bootstrap ソース zip を
   GitHub から DL し libsass でコンパイルしていた (低速・脆弱、かつ libsass は
   非推奨で Bootstrap 5 サポート外)。カスタマイズ実態が「プライマリ色ほか数点」
   だったため、BS5 の CSS 変数による上書き (`theme-overrides.css`) で置換。
4. **Bootstrap 4.5 → 5.3** — 旧参照先 stackpath BootstrapCDN は 2024 年に
   サービス終了しており、**Bootstrap JS が実際にはロードされずモバイルメニューが
   機能していなかった**。BS5 + SRI 付き jsDelivr に移行し、ナビバーを BS5 作法で
   再構成。多階層ドロップダウンの外部 hack (raw.githubusercontent 直リンクで
   MIME 的にも無効だった) は自前の voidy-menu.js + CSS に置換。
5. **jQuery 廃止** — BS5 が jQuery 不要になったのに合わせ、theme 内の全スクリプト
   (TOC, scrolltop, modified 表示, Colab バッジ, ドロップダウン) を vanilla 化。
   ページ途中で jQuery 3.4.1 を二重ロードする問題も解消。
6. **X タイムライン** — 埋め込み再試行 + タイムアウトフォールバック方式 (§7.5)。
7. **CI 修正** — `[ $cur_hash=$new_hash ]` (常に真) の比較バグ修正、
   `|| true` によるビルド失敗の握り潰し除去、pip キャッシュ有効化。
8. **Bootstrap 5.3.7 の採用理由** — 更新時点の最新は 5.3.8 だが、公式ドキュメント
   で SRI ハッシュを検証できた 5.3.7 を採用 (差分は軽微)。更新する場合は §7.1。

### 既知の残課題

- `autosummary` と `summary` の併存 (`pelicanconf.py` にもコメントあり)
- Font Awesome の `<link>` に SRI 未設定 (`FONT_AWESOME_LINK` に `integrity`
  キーを足せば有効化される)
- MathJax が 2.7.3 (render-math プラグイン依存)。MathJax 3/4 への移行は
  render-math の対応待ちか自前注入への切り替えが必要
- レガシー webhook 機構 (§10) の撤去判断

---

## 12. 開発環境の作り方

```bash
git clone https://github.com/oumpy/hp_management.git
cd hp_management
python3 -m venv .venv && source .venv/bin/activate   # 任意
pip install -r requirements.txt
sh tools/init.sh          # 出力レポジトリの clone (push 権限がなければ省略可)
pelican                   # → output/
pelican -l                # http://localhost:8000 (自動再生成付きは pelican -r -l)
```

- Python 3.9 以降。OS 依存なし (シェルスクリプトは POSIX sh)。
- 本番相当を確認したいときは `pelican -s publishconf.py` (sitemap / robots.txt 付き)。

---

## 13. 変更時のチェックリスト

- [ ] `pelican` によるビルドが警告 ({attach} 関連の既知警告を除き) なしで通るか
- [ ] ipynb 記事 1 本 (例: `/blog/2018/10/bayesian_ttest.html`) の表示・数式・
      Colab バッジ
- [ ] デスクトップ: メニューのホバー展開・クリック遷移・active 表示
- [ ] モバイル幅 (<992px): ハンバーガー開閉・サブメニューのタップ展開
- [ ] サイドバー各 box (検索・SNS・新着・タグ・支援・X timeline フォールバック)
- [ ] `pelican -s publishconf.py` で sitemap.xml / robots.txt / .nojekyll が生成されるか
- [ ] CDN の URL と SRI ハッシュの対応 (Bootstrap 更新時)
