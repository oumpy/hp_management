Title: プレゼン形式の動画からスライドのpdfを自動作成する＋OCRをつける
Date: 2021.03.23
Modified: 2021.03.31
Tags: 自動化, 勉学支援
Author: 山田

昨年度はずっとオンラインの授業であったため、それに役立つようなプログラムを作ったりしていました。
中でも、
「プレゼン形式の動画からスライドのpdfを自動作成する＋OCRをつける」
はかなり役立っています。

zoomなどで授業がある際、教員によってはレジュメを配ってくれない場合もあります。
そんなときは授業動画を録画してmp4ファイルを準備すれば、自動でpdfを作ってくれる、というものです。

## 仕組み
流れとしては、

- mp4を一定時間（10秒など）ごとの画像に分割
- それぞれの画像に対して、opencvで1つ前の画像との差分を取る
- 差分を数値化してある値以上のものをはじく
- 残った画像をpdfに変換して、まとめる

という感じです。詳しくはコードを見てみて下さい。
実行時間は数分なので、待てないほどではないです。

## OCR
OCRを付けたい場合は、**ocrmypdf**が無料で使える中では精度がいいので、学習データを少しいじって、スライド向けのオプションにすればそれなりに認識してくれます。

<https://github.com/jbarlow83/OCRmyPDF>

## コーデックの変換
mp4についても、コーデックをH.265に変換してあげると、画質はほとんど変わらず容量がかなり削減されます。  
もしmp4も残しておきたいといった場合は、変換はffmpegで簡単にできるので、おすすめです。
（この変換はかなり時間がかかりますが...）

参考：<https://life.craftz.dog/entry/save-storage-with-h265-ffmpeg>

## プログラム
プログラムを載せておくので、使い方が分からない時は聞いて下さい。

### movie2pdf.py
```python
# 使い方
# python movie2pdf.py -dir 実行フォルダのパス

import argparse
import subprocess
import glob
import cv2
import shutil
import os
# pdfへ変換
from PIL import Image
import img2pdf
# pdfをまとめる
import re
import PyPDF2
from send2trash import send2trash

# コマンドライン引数の取得
parser = argparse.ArgumentParser()
parser.add_argument('-dir')
parser.add_argument('-num_f', default=0.5)
args = parser.parse_args()

base_dir = args.dir
num_frames = args.num_f
# floatに変換
num_frames = float(num_frames)

# avi→mp4へ変換
input_avi = base_dir + "input.avi"
if os.path.exists(input_avi):
    cmd = 'ffmpeg -i {0} -vcodec libx265 -tag:v hvc1 {1}input.mp4'.format(
        input_avi, base_dir)
    subprocess.run(cmd, shell=True)
    print("Done! - avi2mp4")
else:
    print("Not Avi file")

# mp4→連番のpng
input_mp4 = base_dir + "input.mp4"
os.makedirs(base_dir + "photo")
# 5秒に1枚のpngへ
# -rは1秒あたりのコマ数
cmd = 'ffmpeg -i {0} -vcodec png -r {1} {2}photo/%05d.png'.format(
    input_mp4, num_frames, base_dir)
subprocess.run(cmd, shell=True)
print("Done! - mp4topng")

# 連番png→重複除いてpdf作成
# ２値化


def binarization(img, threshold=100):
    ret, img = cv2.threshold(img, threshold, 255, cv2.THRESH_BINARY)
    return img

# 関数
# 差分を数値化


def getDiff(path1, path2):
    # 画像の読み込み
    img1 = cv2.imread(path1)
    img2 = cv2.imread(path2)
    # グレースケール変換
    img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)
    # 差分取得
    mask = cv2.absdiff(img1, img2)
    # ２値化
    mask = binarization(mask)
    return cv2.countNonZero(mask)  # 白の要素数


# 重複を除く
files = glob.glob(base_dir + "photo/*.png")
files = sorted(files)

list_frame_save = []  # 保存する画像のリスト
threshold = 10000  # しきい値

while True:
    if len(files) > 1:
        diff = getDiff(files[0], files[1])
        if diff == 0:
            files.remove(files[1])
        elif 0 < diff < threshold:
            files.remove(files[1])
        else:
            print(files[0])
            list_frame_save.append(files[0])
            files.remove(files[0])

    else:
        break

list_frame_save = list_frame_save + files

# ファイルのコピー
os.makedirs(base_dir + "save")
for i in list_frame_save:
    file_name = os.path.basename(i)
    shutil.copyfile(i, base_dir + "save/" + file_name)

# pdfへ変換
# png→jpg
# 保存先のディレクトリを作成
os.makedirs(base_dir + "jpg")

list_png = list_frame_save
for i in list_png:
    name, ext = os.path.splitext(i)
    name_png = name.split("/")[-1] + ".jpg"
    # 変換
    im = Image.open(i)
    rgb_im = im.convert('RGB')
    rgb_im.save(base_dir + "jpg/" + name_png)

# jpg→pdf
path_jpg = base_dir + "jpg/*.jpg"
list_jpg = glob.glob(path_jpg)
# 保存先のディレクトリを作成
os.makedirs(base_dir + "pdf")

for i in list_jpg:
    name, ext = os.path.splitext(i)
    name_pdf = name.split("/")[-1] + ".pdf"
    # 変換
    # Pillowモジュールを使用し画像の読み込み
    img = Image.open(i)
    # 画像→pdfファイルに変換
    cov_pdf = img2pdf.convert(i)
    # pdfファイルを読み込み（pdf_nameで指定したpdfがない場合、pdf_nameをファイル名として新規にpdfファイルを作成）
    file = open(base_dir + "pdf/" + name_pdf, "wb")
    # pdfファイルを書き込み
    file.write(cov_pdf)

    # 開いているファイルを閉じる
    img.close()
    file.close()

# 複数のpdfファイルを結合する
pdf_path = base_dir + "pdf/"

merge = PyPDF2.PdfFileMerger()
for j in sorted(os.listdir(pdf_path), key=lambda s: int(re.search(r'\d+', s).group())):
    merge.append(pdf_path + "/" + j)
merge.write(base_dir + 'slide.pdf')
merge.close()

# フォルダとファイルを削除
shutil.rmtree(base_dir + "jpg")
shutil.rmtree(base_dir + "pdf")
send2trash(base_dir + "photo")  # photoフォルダはごみ箱へ移動
if os.path.exists(input_avi):
    send2trash(base_dir + "input.avi")
else:
    pass

print("Done! All")
```

実行方法
```
python movie2pdf.py -dir 実行フォルダのパス
```

### movie2pdf_folder_ocr.sh
```bash
for file in *.mp4; do
# echo $file
name=$(basename $file .mp4)
echo $name
# タイトルのフォルダを作成
mkdir $name
# mp4を移動＆名前変更
mv $file $name/input.mp4

# 変換実行
python ~/python/program/201214_mp4pdf/script/movie2pdf.py -dir ./$name/ -num_f 0.5
# OCR
ocrmypdf -l jpn --output-type pdfa --tesseract-pagesegmode 6 --sidecar $name/output.txt $name/slide.pdf $name/slide_ocr.pdf
# outputのslide.pdfのタイトル変更
mv $name/slide.pdf $name/$name.pdf
mv $name/slide_ocr.pdf $name/$name_ocr.pdf
done
```

指定のフォルダ内にて、  
```
sh movie2pdf_folder_ocr.sh
```


## 注意事項
授業の録画は禁止されている場合もあると思うので、プログラムの利用は全て自己責任でお願いします。

オンライン授業はテックに強いと快適に過ごせるため、今後もずっとオンラインでいいなと思ってしまいます。
