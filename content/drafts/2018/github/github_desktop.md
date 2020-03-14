Title:github desktopアプリの使い方
Date: 2019.04.14
Category: github
Tags: github
Slug: github_desktop
Author: 安水
Summary:github desktopアプリの使い方

githubとはなにかの説明は省略しています。githubについて知らない人は先にいろいろ調べてみましょう。今回は[github desktop](https://desktop.github.com/)についてです。

### レポジトリの読み込み（プロジェクト開始時の一回のみ。）

<img class="alignnone size-full wp-image-284" src="https://pythonoum.files.wordpress.com/2018/10/1.png" alt="1" width="550" height="202" />

- github上のレポジトリとLocal pathを指定。

<img class="alignnone size-full wp-image-285" src="https://pythonoum.files.wordpress.com/2018/10/2.png" alt="2" width="900" height="968" />

これでうまくcloneされていたらOK


---

### 毎回の流れ(基本編)

- 一日の作業の初めにFetchしてきて、最新の状態に同期する。

<img class="alignnone size-full wp-image-286" src="https://pythonoum.files.wordpress.com/2018/10/3.png" alt="3" width="928" height="398" />

ファイルを変更、加筆すると自動的にアプリに反映される

<img class="alignnone size-full wp-image-287" src="https://pythonoum.files.wordpress.com/2018/10/4.png" alt="4" width="2144" height="1544" />

- 作業の区切りがついたら左下のSummeryに適当にコメントを付け、commit to master。
コメントはもうちょっと丁寧に付けましょう。。。

<img class="alignnone size-full wp-image-288" src="https://pythonoum.files.wordpress.com/2018/10/5.png" alt="5" width="620" height="518" />

<img class="alignnone size-full wp-image-289" src="https://pythonoum.files.wordpress.com/2018/10/6.png" alt="6" width="2144" height="1544" />

- pushして変更を反映させる。

まとめると、**fetch(pull) -&gt; commit -&gt; push** のライフライクルです。

### CLIでの操作

デスクトップアプリはwin,macではあるが、Ubuntuでは無いらしいので、同じ作業をCLIでやる必要がある。ただし、基本は同じで

```bash
$ git pull origin master
$ git add --all
$ git commit -m "コメント"
$ git push origin master
```

って感じでやればOK。branchやmerge、プルリクエストについては触れていません。
