Title:必要n数の決定～研究計画～
Date: 2019.03.22
Tags: Statistics
Author: 依藤

## シングルセル解析での細胞数
① **実現可能性**：金銭面、使用機器のスループットを考慮した細胞数  
② **検出力**：統計学的検討により細胞種間の差が充分に検出される症例数であるか  
③ **倫理面**：ヒトサンプル、あるいはモデル動物サンプルの取り扱い  

とくに②に関して：

- 分類するclusterに関して過去の文献から予測する。
- 想定されるcluster数
- 1clusterあたりの最低細胞数

上記の見積もりができれば、ウェブサイトでの計算が可能。

## 臨床研究での症例数
① **実現可能性**：それだけの症例数を集めることができるか  
② **検出力**：統計学的検討によりoutcomeの差が充分に検出される症例数であるか  
③ **倫理面**：不必要に多くの被験者を研究に参加させるべきではない  

とくに②に関して：

症例数の予測に必要なものを過去の文献やパイロット研究を行うことで予測する。

- outcomeの差
- データのばらつき(標準偏差)

上記の見積もりができれば、ソフトウェアでの計算が可能。 (cf.) [Power and Sample Size Calculation (PS)](http://biostat.mc.vanderbilt.edu/wiki/Main/PowerSampleSize)

## まとめ
シングルセル解析では、やはり値段や機器の問題が大きくなってくるでしょうか。特に細胞種がかなり希少な場合に、細胞数(n)をどんどん増やすしかないとなると、その点が一番気になります。

一方で臨床研究では、単施設で可能なのか、多施設共同でないと症例数(n)が集まらないのかという問題は切実です。そしてヒトを対象にする以上、症例数決定の明確な理由が必要です。[CONSORT2010 checklist](http://www.consort-statement.org/media/default/downloads/consort%202010%20checklist.pdf)にも記載があります。

いずれにせよ、事前の文献調査やプレ実験・パイロット研究によって着地点を明確にしておかないと、途中で頓挫してしまい、出るはずの結果がでない事態になりそうです。今後はシングルセル解析で前提となる統計手法に関しても勉強していこうと思っています。

## 参考文献
- 今日から使える医療統計　(新谷歩, 2015年, 医学書院)
- Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies　(Atefeh Lafzi et el., Nature Protocols volume 13, pages2742–2757 (2018))
