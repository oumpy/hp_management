Title:GitHub desktopアプリの使い方
Date: 2018.10.28
Tags: GitHub
Author: 安水

githubとはなにかの説明は省略しています。githubについて知らない人は先にいろいろ調べてみましょう。今回は[github desktop](https://desktop.github.com/)についてです。

## レポジトリの読み込み（プロジェクト開始時の一回のみ。）

![1]({attach}images/github_desktop_figs/2018_10_1.png)

github上のレポジトリとLocal pathを指定。

![2]({attach}images/github_desktop_figs/2018_10_2.png)

これでうまくcloneされていたらOK

---

## 毎回の流れ(基本編)

一日の作業の初めにFetchしてきて、最新の状態に同期する。

![3]({attach}images/github_desktop_figs/2018_10_3.png)

ファイルを変更、加筆すると自動的にアプリに反映される

![4]({attach}images/github_desktop_figs/2018_10_4.png)

作業の区切りがついたら左下のSummeryに適当にコメントを付け、commit to master。
コメントはもうちょっと丁寧に付けましょう。。。

![5]({attach}images/github_desktop_figs/2018_10_5.png)

pushして変更を反映させる。

![6]({attach}images/github_desktop_figs/2018_10_6.png)

まとめると、**fetch(pull) -&gt; commit -&gt; push** のライフライクルです。

## CLIでの操作

デスクトップアプリはwin,macではあるが、Ubuntuでは無いらしいので、同じ作業をCLIでやる必要がある。ただし、基本は同じで

```bash
$ git pull origin master
$ git add --all
$ git commit -m "コメント"
$ git push origin master
```

って感じでやればOK。branchやmerge、プルリクエストについては触れていません。
