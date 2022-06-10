Title: 英語のカルテデータのコンペに参加してみた(際の失敗と反省)
Date: 2022.05.24
Modified: 2022.06.08
Tags: Data Science Competition
Author: 杉原

## この記事の趣旨
医学生かつ技術がある程度わかる人としてkaggleに関わる機会を得た。
その経験を通して学んだことなど今から書く。臨床を学んだ医師がコーディング以外の強みを活かしてテクノロジーに関わっていくとはどういった形のものなのか、イメージがつくようになってもらえたら嬉しい。なお、今回は自然言語処理周りの専門家の方と行ったのでモデル周りのことは全く把握していない。悪しからず。

## コンペのタスクの概要
アメリカの医師試験で行われる、模擬診察の時の学生が書いたメモを扱ったコンペになる。
この医学生が書いたメモの自動採点をするために、その症例でこぼしてはならない重要な情報や症状を固有表現抽出 (NER) タスクをして抽出することを行う。
[詳しくはこちら](https://www.kaggle.com/code/yufuin/nbme-japanese/notebook)

## やったことと考えたこと
1. 教師データと現状のモデルでの出力結果を眺める
    ![fig1]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig1.png)
    上記のデータは上が教師データ、下が[uth-bert](https://github.com/noroka/uth_bert)を元に学習させて出力したデータになる。
    そもそも教師データ自体に何の法則性もなく、年齢が含まれていたり性別が含まれていなかったりした。
    あまりにもランダムだと感じられたので、次の手に移行。

2. 出てくる単語を診療科の単語毎に区分け
    まず教師データ、テストデータ、評価データの性質を均等化しようと考えた。
    ここで取り入れた評価軸は２つ。
    患者が男性or女性
    患者の想定される疾患「1.循環器,2.消化器,3.婦人科,4.精神科,5.循環器と消化器,6.その他」
    ![fig2]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig2.png)
    上記のように頻出単語にマークし、それを主訴により5種類の疾患分野に分類した。

3. 教師データとテストデータを均等化した学習結果
    〜1.循環器,2.消化器,3.婦人科,4.精神科,5.循環器と消化器,6.その他が教師データとテストデータに占める割合を調整してみた〜

    分け方1:何も分けない(下記)
    ![分け方1]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig3.png)
    分け方2:男性or女性で分ける,疾患で分けない(下記)
    ![分け方2]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig4.png)
    分け方3:男性or女性で分けない,疾患で分ける(下記)
    ![分け方3]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig5.png)
    分け方4:男性or女性で分ける,疾患で分ける(下記)
    ![分け方4]({attach}./images/medicalrecord_competition_figs/sugiharacompe_fig6.png)
    →分け方4はf1値が最高！分けてよかったね。(めっちゃ微差だけど)

4. 顛末
    モデルが完全一致の数で精度を測るものになってしまっており、部分一致を評価できないモデルだったらしく全員意気消沈。今回は諦めることに。

5. 最後に：感想戦
    記録上位のチームの記事を見てもらってもわかる通り、試行回数を積むことが大事だということがわかるだろう。
    [一位のチームの記事](https://www.kaggle.com/competitions/nbme-score-clinical-patient-notes/discussion/323095)
    [二位のチームの記事](https://www.kaggle.com/competitions/nbme-score-clinical-patient-notes/discussion/323085)

実は取り組み始めたのは締切の3週間前、1.の出力の法則性があまりないことに気づいたのは2週間前だった。とにかく試行して提出して改善を繰り返せたチームが勝つことを思い知らされた。みなさん、早めに取り組み始めてくださいね。
