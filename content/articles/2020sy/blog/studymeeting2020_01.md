Title: 【勉強会資料 2020 第1回】Python会とは? / Pythonの基礎
Date: 2020.05.19
Modified: 2020.05.23
Tags: 勉強会, Python
Author: 山本

2020.05.19に開催したオンライン勉強会の資料を公開します。

今回のみ、会外部の一般の方々にも参加頂く形で行われました。
次回からは会内部のみとなりますが、資料は公開してゆく予定です。

## スライド
[% embed src="{attach}./attach/studymeeting2020_01/200519_python会入門.pdf" %]

## Pythonの環境構築
- [Pythonの公式サイト](https://www.python.org/)からPythonだけをinstallしても使えるが、パッケージ管理や環境管理のためにAnacondaを使うと便利（ここで不要じゃという声が飛んでくる）。
- 少なくともWindowsだと使った方が楽
- Anacondaそのままは全部盛りに近い（＝要らないものなども自動でinstallされる）ので、最小限のパッケージに留めたMinicondaをinstallする方がよい
- Miniconda入れる前にpyenv入れた方が良いという話もある(MacOS, Linux)

## Pythonファイルの実行

```sh
python hoge.py
```
で実行可能。
IDEでは簡単に実行出来たり、他にはJupyter notebookを用いる手法も。

## Python Tutorial
[Jupyter Notebook]({attach}./attach/studymeeting2020_01/Python_tutorial.ipynb)

コードのみ以下に掲載します。
### NumPyで数値計算

```python
import numpy as np
x = np.arange(0, 1, 0.1)
print('x:', x)
y = 2*x + 10
print('y:', y)
```

### Matplotlibで図を描画する
Matplotlibは図を描画するためのライブラリ。
どんな図が描画できるかは公式ドキュメントの[ギャラリー](https://matplotlib.org/gallery/index.html)を参照。

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 1, 0.1)
y = 2*x + 10

plt.plot(x, y)
plt.show()
```

### Seabornによる回帰直線の描画
SeabornはMatplotlibがベースのライブラリ。
より綺麗で複雑な図を簡単に描画できる。
公式ドキュメントの[ギャラリー](https://seaborn.pydata.org/examples/index.html)を参照。
```python
import seaborn as sns
np.random.seed(0) # 乱数seedの設定

# Toyデータ(y=2*x+10+noise)の作成
n_data = 100 # データ数
x = np.random.randn(n_data)
y = 2*x + np.random.randn(n_data) + 10

# 回帰直線のplot
plt.figure(figsize=(4,4))
sns.regplot(x, y) # 回帰(regression)の実行
plt.show() # 画像表示
```

### Markdownについて
Markdownはプレーンテキスト形式で手軽に書いた文書からHTMLを生成するための言語
- 書き方は[Markdown記法サンプル集](https://qiita.com/tbpgr/items/989c6badefff69377da7)等を見るとよい
- 普段から使う場合は[Typora](https://typora.io/)を使うのがおススメ。
TeX形式の数式も書ける。
PandocをインストールすればWordやLaTeXなどに変換可能。

### Google Colabについて
[Google Colaboratory](https://colab.research.google.com)はJupyter Notebookをブラウザ上で使えるようにGoogleが提供しているサービス。
Pythonをインストールする必要はありません。
また、Python会の一部のブログは記事をそのままGoogle Colabで開くことができます(GitHubにアップロードしたJupyter Notebookファイル(.ipynb)はURLを修正するだけでColabで開くことができます)。
