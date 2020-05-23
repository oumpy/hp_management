Title: JuliaとPythonと競技プログラミング
Date: 2019.03.16
Tags: 競技プログラミング
Author: 小川

最近友人と話していてJulia (<https://julialang.org>)が話題になったことがあったので、少しだけ調べて試してみた話。

## Juliaってなに？

Pythonのような動的型付けのスクリプト言語です。実行時にコンパイルを行いC言語にも迫る実行速度、科学技術系の数値計算もどんと来い、という触れ込みで、人気上昇中らしいです。Pythonを含む他言語のライブラリを読み込む仕組みを備えているのもすごいところ。  
<s>ただ、コードの見た目が激しくMatlab風味で思わず目を背けたくなります。</s>

## AtCoderでのJulia

ご多分に漏れず、実行速度に惹かれました。半年ほど前から参加している競技プログラミングサイト**AtCoder** (<https://atcoder.jp>)で、常用しているPython（← Python会だからね！）で実行時間切れ、C/C++に書き換えると正答、という経験を何度かしてきたので。

そういうわけで、ここでは最近の**第121回 AtCoder Beginner Contest (ABC121)** (<https://atcoder.jp/contests/abc121/>)の問題で、Juliaのパフォーマンスを実際に見てみます。  
（言語知識0から数十分調べて書いたコードのため、動きはするが思わぬところで非効率、ということはあるかもしれません。）

### AtCoder Beginner Contest 121 問題A "White Cells"

（問題はこちら → <https://atcoder.jp/contests/abc121/tasks/abc121_a>）

入力値4個を読み込んで簡単な演算結果を返すだけの問題です。

```julia
# ABC121-A "White Cells" in Julia
inpl() = map(parse,split(readline()))
(H,W) = inpl()
(h,w) = inpl()
println((H-h)*(W-w))
```

所要時間が入力にほぼ依存しない問題ですが、、各言語でのAtCoder上実行結果。

|   言語    |   実行時間 | 消費メモリ |
| :-------: | ---------: | ---------: |
| **Julia** | **360 ms** | **110 MB** |
|   PyPy3   |     180 ms |      40 MB |
|  Python3  |      17 ms |       3 MB |
|    C++    |       1 ms |     256 kB |

C++にもPythonにも、実行時間と消費メモリ双方で惨敗。   
おそらくですが、**実行ごとにまずコンパイルを行う**ので、簡単な問題だとそれが相対的に巨大なオーバーヘッドになってしまうようです。  
Pythonの半コンパイラ型実装であるPyPyはJuliaの半分程度でした。

<!---

### AtCoder Beginner Contest 121-B

問題はこちら → <https://atcoder.jp/contests/abc121/tasks/abc121_b>

```julia
# ABC121-B by Julia
inpl() = map(parse,split(readline()))
(N,M,C) = inpl()
B = input()
A = Array{Int}(N,M)
for n=1:N
	A[n,:] = inpl()
end
ans = 0
D = A*B+C
for d in D
  if d > 0
  	ans+=1
  end
end
println(ans)
```

--->

### AtCoder Beginner Contest 121 問題C "Energy Drink Collector"

（問題はこちら → <https://atcoder.jp/contests/abc121/tasks/abc121_c>）

読み込んだ値の列をソートして、条件判定をしながら順番に足し上げていく問題。

```julia
# ABC121-C "Energy Drink Collector" in Julia
inpl() = map(parse,split(readline()))
(N,M) = inpl()
A = Array{Int}(N,2)
for i in 1:N
  A[i,:] = inpl()
end
A = sortrows(A, by=x->x[1])
ans = 0
for i in 1:N
  if A[i,2] >= M
    ans += M*A[i,1]
    break
  else
    M -= A[i,2]
    ans += A[i,1]*A[i,2]
  end
end
println(ans)
```

なんと、テスト16個中15個でタイムアウト（2000ms以上）。。。

Pythonのコードはこちら：

```python
# ABC121-C "Energy Drink Collector" in Python3
inpl = lambda: list(map(int,input().split()))
N,M = inpl()
A = []
for i in range(N):
  A.append(inpl())
A.sort(key=lambda x: x[0])
ans = 0
for i in range(N):
  if A[i][1] >= M:
    ans += A[i][0]*M
    break
  else:
    M -= A[i][1]
    ans += A[i][0]*A[i][1]
print(ans)
```

こちらは最大466msでクリア（同一コードのPyPy3では733ms）。ここで134msのテストもJuliaではタイムアウト。悲しい。

今回Juliaでタイムアウトになったのは、言語の特性や正しいコーディングの仕方を知らないから、という可能性は高いです。ただ実際にこの問題でJuliaを使って提出されている答案は、最速クリアのものでも1724ms（そもそもJuliaでの提出数自体が少ないですが）。やはり上の安直なPythonコードが圧勝しています。

## Juliaは競技プログラミングに向かない？

Juliaの実行速度が速いこと自体は（今回検証していませんが、きっと）本当なんだと思います。しかしそれは時間のかかる複雑・大規模な処理の場合であって、競技プログラミングのような高々2-3秒の計算にはコンパイルのオーバーヘッドがやはり大きいのかな、という印象です。  

AtCoderなどでは、C++などのあからさまなコンパイル型言語はコンパイル時間が実行時間に算入されず、一方でJuliaに対しては算入されます。やや理不尽な感じはしなくもないけれど、このルール下でJuliaの高速性能を生かすことは（あくまでAtCoderのような短時間型競技プログラミングの話ですが）なかなか難しそう。残念ながら、普通にやるとPythonよりもずっと遅い。  

なので、当面の競技プログラミング用言語は**やっぱりPython**（とC/C++）、と個人的には結論づけたところです。おしまい。

**< History >**  
2019.03.21 ver1.0  
2019.03.21 ver1.1 微修正、PyPy3の成績を追記  
2019.03.21 ver1.2 JuliaとPythonのコードが対応するよう修正
