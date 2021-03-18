Title: attention勉強会2
Date: 2021.03.14
Modified: 2021.03.19
Tags: Machine Learning,勉強会,論文まとめ
Author: 安部

2021/03/14開催の勉強会の資料です。
最近attentionの流行を肌で感じだし、焦って勉強しています。

## スライド
[% embed src="{attach}./attach/studymeeting2021_03/transformer勉強会2.pdf" %]

## maskについての訂正
- softmax(QK)についてではなく、QKにmaskをかけます。

## Linformer
- D(Q,K,V)=softmax(QK)VのQKが低ランク
- L(Q,K,V)=softmax(Q(EK))(FV)と、E,Fをかけて次元を落とす
- maskはQKにかけるのではなく、K,V,(Q)にかける
- 高速化は未確認...
