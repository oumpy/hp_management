Title:(論文まとめ)医師国家試験問題自動生成AI
Date: 2019.04.14
Category:
Tags: deeplearning
Slug: doctor_exam_ai
Author: 秋山
Summary:

On the Generation of Medical Question-Answer Pairs.
Shen S, Li Y, Du N, Wu X, Xie Y, Ge S, et al. arXiv. 2018.
http://arxiv.org/abs/1811.00681

上記論文を軽く紹介。
Tensent Mediacal AI lab Internからの論文。
Deep learning技術により質問文に対して回答するAI (question answering, QA)が発展している。しかし、QAを医療に応用するためにはAIを学習させるためのデータが不足している。そこで質問文と解答のペアを自動生成するモデルを提案した。

手法
Key Phrase Detector
質問文の各フレーズが解答の決め手となるキーフレーズであるかを評価する。キーフレーズであるかどうかは特定の解答に対して高頻度で質問文に出現するフレーズがキーフレーズであるとしてdetectorを学習させる。例えば「日本脳炎」が解答である場合「項部硬直」などがキーフレーズとなる。
Conditional Variational Autoencoder (CVAE)
キーフレーズは維持しつつ、それ以外のフレーズを生成モデルCVAEによって言い換える。これによって答が同じな新たな質問文が作られる。

&nbsp;

<img class="alignnone size-full wp-image-427" src="https://pythonoum.files.wordpress.com/2018/12/393E4E28-4B67-4C25-A11B-BE7CA70BB9AC.png" alt="393E4E28-4B67-4C25-A11B-BE7CA70BB9AC" width="2055" height="880" />

データ
中国医師国家試験18,798問
中国のWikipedia風医療サイト (http://xywy.com/)
医学辞書19冊
医学論文2,130,128本
医学専門書518冊

実験結果
アルゴリズムによる評価、人間による評価ともにベースラインを上回った。

読んだ感想
生成された問題文の例が載ってないのですごいのかよくわからなかった。データの量はすごい。
