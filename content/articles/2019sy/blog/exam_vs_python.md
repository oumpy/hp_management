Title: 試験 vs Python
Date: 2019.09.03
Modified: 2019.09.03
Tags: 勉学支援
Author: 小川

## クイズ丸暗記選手権
あさっては某科目の試験です。
クイズみたいな過去問をそのまま出すから理屈抜きで丸暗記しろという、哲学的な出題だそうです。
そして僕は諸事情によりまったく勉強<s>して</s>できていなかったので、何とかしないといけません。
<s>僕が悪いんじゃない世のなかが悪いんだ。</s>

こういうときは、もうPythonでスクリプトを書くしかないですね。  
過去問からコピペした問題リストより、ランダムに順次出題します。  

## ソースコード
```python
###
# 無限クイズ出題システム quiz.py
# Usage: $ python quiz.py problems_datafile.tsv
#
# 問題→（キー）→解答→（キー）→問題→（無限に続く）
# 'q' or 'Q' or ESC で終了。
###

import sys
import random
from getch import getch
#from farbric import colors
BLACK = '\033[30m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
PURPLE = '\033[95m'
CYAN = '\033[96m'
WHITE = '\033[97m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
ESC = 27

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'problems.tsv'

quit_ch = {ord('q'),ord('Q'), ESC}
with open(filename,'r') as f:
    problem_set = [ line.rstrip('\n').split('\t') for line in f ]
N = len(problem_set)
n = 0
k = -1
while True:
    n += 1
    k_prev = k
    while k == k_prev or len(problem_set[k]) != 2:
        k = random.randint(0,N-1)
    prb, ans = problem_set[k]
    print('第{}問'.format(n))
    print(BLUE, prb, ENDC)
    key = ord(getch())
    if key in quit_ch:
        break
    else:
        print('解答')
        print(RED, ans, ENDC)

    key = ord(getch())
    if key in quit_ch:
        break
    else:
        print()
```

## 実行例
$ python quiz.py parasitology_problems.tsv  
第1問  
 <span style="color:blue">世界的に広く分布し、日本においては、STDの一種である。5億人の感染者の〜95%は無症状感染者 (asymptomatic cyst carrier) である。腸管に寄生し、ミトコンドリアが変化したマイトソームという細胞小器官を持つ。</span>  
解答  
 <span style="color:red">赤痢アメーバ</span>  

第2問  
 <span style="color:blue">成虫は消化管粘膜に寄生し幼虫を産出するが、幼虫は宿主の外に出ることなく同じ宿主の横紋筋で被嚢し、この宿主がほかの肉食動物に捕食されるのを待つ。</span>  
解答  
 <span style="color:red">旋毛虫</span> 

第3問  
 <span style="color:blue">奄美、沖縄になお多くの感染者が存在し、自家感染で虫体数が増加する。</span>  
解答  
 <span style="color:red">糞線虫</span>  

第4問  
...

## 結論
まじめに勉強するしかないようです。  
<s>（試験はクイズじゃなくて競プロにしてくれたらいいのに）</s>

おしまい。
