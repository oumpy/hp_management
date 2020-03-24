Title:必要n数の決定～研究計画～
Date: 2019.04.14
Tags: statistics
Slug: n_number
Author: 依藤
Summary:

<span style="font-weight:400;">実験計画(シングルセル解析)を考えるときに、臨床研究だとどうだったかなと思うことがあります。</span>

<b>〈シングルセル解析〉group Aの1000細胞　と　group Bの1000細胞をサブタイピング</b>

<b>〈臨床研究〉介入Aの結果、心血管イベントありの1000人　と　心血管イベントなしの1000人</b>

<span style="font-weight:400;">細胞はヒトの数ほど多彩であり、それを分類するのがシングルセル解析であり、(勿論、データ量の桁が異なりますし、結果を見ているのか背景を見ているのかという違いはありますが)臨床研究のデータ解析の考え方を活かせるように思います。今回はサンプル調整に必要な細胞数(n)の決定に関して、臨床研究での症例数(n)の決定と比較してみました。</span>
<ul>
	<li><span style="font-weight:400;">シングルセル解析での細胞数</span></li>
</ul>
<span style="font-weight:400;">①実現可能性：金銭面、使用機器のスループットを考慮した細胞数</span>

<span style="font-weight:400;">②検出力：統計学的検討により細胞種間の差が充分に検出される症例数であるか</span>

<span style="font-weight:400;">③倫理面：ヒトサンプル、あるいはモデル動物サンプルの取り扱い</span>

<span style="font-weight:400;">とくに②に関して：</span>

<span style="font-weight:400;">分類するclusterに関して過去の文献から予測する。</span>

<span style="font-weight:400;">想定されるcluster数</span>

<span style="font-weight:400;">1clusterあたりの最低細胞数</span>

<span style="font-weight:400;">上記の見積もりができれば、<a href="https://satijalab.org/howmanycells" target="_blank" rel="noopener noreferrer">ウェブサイト</a>での計算が可能。</span>
<ul>
	<li><span style="font-weight:400;">臨床研究での症例数</span></li>
</ul>
<span style="font-weight:400;">①実現可能性：それだけの症例数を集めることができるか</span>

<span style="font-weight:400;">②検出力：統計学的検討によりoutcomeの差が充分に検出される症例数であるか</span>

<span style="font-weight:400;">③倫理面：不必要に多くの被験者を研究に参加させるべきではない</span>

<span style="font-weight:400;">とくに②に関して：</span>

<span style="font-weight:400;">症例数の予測に必要なものを過去の文献やパイロット研究を行うことで予測する。</span>

<span style="font-weight:400;">　outcomeの差</span>

<span style="font-weight:400;">　データのばらつき(標準偏差)</span>

<span style="font-weight:400;">上記の見積もりができれば、ソフトウェアでの計算が可能。</span>

<a href="http://biostat.mc.vanderbilt.edu/wiki/Main/PowerSampleSize" target="_blank" rel="noopener noreferrer">Power and Sample Size Calculation (PS)</a>

<span style="font-weight:400;">　シングルセル解析では、やはり値段や機器の問題が大きくなってくるでしょうか。特に細胞種がかなり希少な場合に、細胞数(n)をどんどん増やすしかないとなると、その点が一番気になります。</span>

<span style="font-weight:400;">　一方で臨床研究では、単施設で可能なのか、多施設共同でないと症例数(n)が集まらないのかという問題は切実です。そしてヒトを対象にする以上、症例数決定の明確な理由が必要です。</span>

<a href="http://www.consort-statement.org/media/default/downloads/consort%202010%20checklist.pdf">CONSORT2010 checklist</a><span style="font-weight:400;">にも記載があります。</span>

<span style="font-weight:400;">いずれにせよ、事前の文献調査やプレ実験・パイロット研究によって着地点を明確にしておかないと、途中で頓挫してしまい、出るはずの結果がでない事態になりそうです。今後はシングルセル解析で前提となる統計手法に関しても勉強していこうと思っています。</span>

<span style="font-weight:400;">【参考文献】</span>

<span style="font-weight:400;">今日から使える医療統計　(新谷歩, 2015年, 医学書院)</span>

<span style="font-weight:400;">Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies　(Atefeh Lafzi et el., Nature Protocols volume 13, pages2742–2757 (2018))</span>

&nbsp;
