Title: LaTeX の展開
Date: 2021.03.30
Modified: 2021.03.30
Tags: その他言語
Author: 河村

## TeXの処理動作
TeXは前からソースコードを認識し，`\textbackslash（英字）`というカタマリ[ref]
制御綴に使える文字はcatcodeが11の文字（英字・日本語）．
従って，日本語の混ざった制御綴も作れる．
`\relax␣あ`の場合は，`\relax`という制御綴と「あ」という文字として認識されるが，`\relaxあ`の場合は，`\relaxあ`という制御綴として認識され（定義しない限り）エラーを吐く。
[/ref]
を一つの制御綴（コントロールシークエンス）として処理します．
```LaTeX
\csone\cstwo            \csthree
```
この場合の制御綴は`\csone`, `\cstwo`, `\csthree`です（空白は制御綴に使われない）．

TeXは制御綴を認識した後，前から順番に展開（定義の代入）を行います．
例えば，次の制御綴を用意します．

- `\csone`の定義：one;
- `\cstwo`の定義：two;
- `\csthree`の定義：three;

この場合は，
```LaTeX
\csone\cstwo            \csthree
\csone\cstwo\csthree
one;\cstwo\csthree
one;two;\csthree
one;two;three;
```
と展開されます（複数の空白は1個の空白として処理され，制御綴の直後の空白は無視される）．

## 展開の制御
先程の例の通り，TeXは基本的に前から展開を行いますが，時には展開の順序を入れ替えたりしたい場合があります．
そのために，いくつかの命令が存在します：

- `\expandafter`．直後の制御綴の展開を一回遅らせる，
- `\edef`．制御綴の中身を全て展開して定義，
- `\noexpand`．`\edef`中で使用すると，定義の際に展開されない（実行時は展開される）．

例として，次の制御綴を用意します．

```LaTeX
\def\A#1#2{#1 is #2}
\def\B{\BB\BB}
\def\BB{foo}
\def\C{bar}
```

それぞれの制御綴の意味は，

- `\A`：引数は2個（\#1, \#2）で，展開すると「\#1 is \#2」（\#1, \#2が制御綴なら展開される）
- `\B`：展開すると「`\BB``\BB`」（`\BB`は展開される）
- `\BB`：展開すると「foo」（それ以上展開されない）
- `\C`：展開すると「bar」（それ以上展開されない）

### 例1
`\A\B\C`
```LaTeX
\A\B\C
% -> \B is \C
% -> \BB\BB is \C
% -> foo\BB is \C
% -> foofoo is \C
% -> foofoo is bar
```
1行目で`\A`を展開する際の，`\A`の引数は`\B`, `\C`です．

### 例2
`\expandafter\A\B\C`
```LaTeX
\expandafter\A\B\C
% -> \A\BB\BB\C
% -> \BB is \BB\C
% -> foo is \BB\C
% -> foo is foo\C
% -> foo is foobar
```
先ず，`\expandafter\A`によって，`\A`の展開は抑制され，先に`\B`が展開されます．

### 例3
`\expandafter\expandafter\A\B\C`
```LaTeX
\expandafter\expandafter\A\B\C
% -> \expandafter\B is \C
% -> \B is \C
% -> \BB\BB is \C
% -> foo\BB is \C
% -> foofoo is \C
% -> foofoo is bar
```
1行目で，先頭の`\expandafter`は直後の`\expandafter`を抑制するため，先に`\A`が展開されます．
2行目で，`\expandafter\B`の\textbf{直後}は制御綴でないため，特に効果はありません．

### 例4
`\expandafter\expandafter\expandafter\A\B\C`．`\expandafter`が3連続の場合は次のように，2つ後の制御綴が先ず展開されます．
```LaTeX
\expandafter\expandafter\expandafter\A\B\C
% -> \expandafter\A\BB\BB\C
% -> \Afoo\BB\C
% -> f is oo\BB\C
% -> f is oofoo\C
% -> f is oofoobar
```
3行目で`\A`の第1引数は`f`，第2引数は`o`です．

### 例5
さらに命令を追加してみます．
```LaTeX
\def\hoge{\B\C}
\A\hoge % -> error
\expandafter\A\hoge % -> \A\B\C -> foofoo is bar
```
`\A\hoge`は`\A`の第2引数が存在しないので，エラーとなります．

### 例6
`\def`の代わりに`\edef`で`\hoge`を定義します．
```LaTeX
\edef\hoge{\B\C} % -> \def\hoge{foofoobar}
\A\hoge % -> error
\expandafter\A\hoge % -> \Afoofoobar -> f is oofoobar
```

### 例7
`\edef`中で`\noexpand`をつければ，その直後の制御綴は（定義の際には）展開されません．
```LaTeX
\edef\hoge{\noexpand\B\C} % -> \def\hoge{\Bbar}
\A\hoge % -> error
\expandafter\A\hoge % -> \A\Bbar -> \B is bar -> \BB\BB is bar -> foofoo is bar
```
`\A\Bbar`で，第1引数は`\B`，第1引数は`b`です．

### 例8
別の例．
```LaTeX
\def\twice#1{#1#1}
\twice\twice hoge % -> error
\expandafter\twice\twice{ho}ge % -> \twicehohoge -> hhohoge
\expandafter\twice\twice{{ho}}ge % -> \twice{ho}hoge -> hohohoge
```
1つ目は無限ループになります．

## 脚注
