Title: deeplabcutの使い方まとめ
Date: 2019.12.16
Modified: 2020.05.16
Tags: Python
Author: 味岡

## DeepLabCut とは
deeplabcutは、動物の行動実験などに使える動画のトラッキングツールです。
package化されているので簡単に使うことができます。
googlecolabでもできるそうです。
詳しい使い方は下記の資料が参考になります。

日本語資料：  
<https://qiita.com/riichirohira/items/b92723278ba92fb42db9>  
<https://qiita.com/Cony-san/items/62df5365878baaa5e275>

公式サイト(github)：  
<https://github.com/AlexEMG/DeepLabCut>


以下は資料にはない僕が使っていく中で気づいた点のメモ書きです。

## 各OSでのインストールの時に気をつけること
### Windows
- インストールはできる (yamlなしでも) 。
- プロジェクトを作るのに管理者権限が必要。
- ROI取りまでは出来た。
しかしその後はエラーが多発した。  
(hdf5ファイルが出来ない、numpy import errorなど)

(追記：2019/12/9時点 最新版のdeeplabcutではwindowsでも実装できた（研究室の先輩談）)

### Mac OS
- インストールは、`dlc-macOS-CPU.yaml` をダウンロードして、
```bash
conda env create -f dlc-macOS-CPU.yaml
```
で環境を作成 (`pip install deeplabcut` はNG)

- GUIでのROI取りをするため、`python` ではなく、「`pythonw`」で起動

### Linux (Ubuntu)
- pipでもインストール可能、普通に起動すればOK

## 補足情報
- 日本語資料には動画はAVIではないといけないとあったが、実際はMPEGでもMP4でも大丈夫だった。
また、tensorflow-gpuは1.4.0でも大丈夫だった。(cuda8.0、cuDNN6、python3.5)
- ROIとりでは、目的部位が隠れた場合には点を打たないこともできる
- 学習はlossが0.01以下で一定になってきたらOK
- ROIとり枚数の基準としては、20～50枚くらいが良さそう

## ソースコードの例 (Ubuntuの場合)
### トレーニング (train.py)
```python
import deeplabcut
deeplabcut.create_new_project('日付','名前',['ファイル名.avi'])
path_config = "生成されたconfig.yamlのパス"
deeplabcut.extract_frames(path_config, 'automatic', 'kmeans')
deeplabcut.label_frames(path_config)
deeplabcut.create_training_dataset(path_config_file)
deeplabcut.train_network(path_config_file, shuffle=1, displayiters=10,saveiters=500)
```
- ホームフォルダにモデル用のフォルダが生成されます。
- configファイルを編集することで設定を変えられます。
- `deeplabcut.label_frames` の部分でGUIが起動してROI取りできるようになります。
- 走らせると学習がずっと進むので、十分学習できたらCtrl+Cで止めましょう。

### 推論 (inference.py)
```python
import deeplabcut
path_config = "使いたいモデルのconfig.yamlのパス"
videofile_path = ['ファイル名.avi']
deeplabcut.analyze_videos(path_config_file,videofile_path)
deeplabcut.create_labeled_video(path_config_file,videofile_path)
deeplabcut.analyze_videos_converth5_to_csv("videoがあるフォルダのパス")
```
- 生成物として、h5ファイル、ラベルがついた動画ファイル、追跡の座標のcsvファイルが作られます。
