Title: SwinTransformerで物体検出
Date: 2021.12.18
Modified: 2021.12.18
Tags: Machine Learning
Author: 安部


## なんの記事か

最近流行りのSwinTransformerを使って物体検出を触ってみます。[SIIM2021 3rd解法](https://www.kaggle.com/c/siim-covid19-detection/discussion/263654)では Swin Transformer with RepPoints(mmdetection)がアンサンブルの種として使用されています。自分で手を動かせばわかることがあるかもしれないと思い過去コンペで追試をしてみます。

行った追試は以下の通りです。

1. swin_tiny mask_r_cnnの公式configをいじって使う
2. augmentationを追加
3. さらにaugmentationを強める
4. backboneを大きいモデルに変更する
5. 他のモデルを試す
6. cocoの事前学習が有効か見てみる


### 追試に用いたデータ

kaggleの過去コンペ&csv提出ができる&自分のドメインに近い、という理由で[VinBig](https://www.kaggle.com/c/vinbigdata-chest-xray-abnormalities-detection/overview) を追試の対象としました。<br>
アノテーションはdiscussionに多く挙げられていた「複数のあのテーターのbboxはナイーブにWBFで統合する」という方法で前処理し、single foldで学習、subをしています。手元のmAP0.5とprivateでのスコアを確認していきます。(スコアの確認が歯切れの悪い形になっています。ご了承ください)<br>以下全てのsubに[2cls fillter](https://www.kaggle.com/c/vinbigdata-chest-xray-abnormalities-detection/discussion/210351) をつけています。コンペ特有の課題(複数のアノテーターへの対処)は本記事では扱いません。推論のコードは[これ](https://www.kaggle.com/kc3222/mmdet-pytorch-framework-infer-vfnet)を参考にしました。


### 実験

#### 1. swin_tiny mask_r_cnnの公式configをいじって使う


[mmdetectionの公式のレポジトリ](https://github.com/open-mmlab/mmdetection/tree/master/configs/swin)にはmask_Rcnn系列のbackboneにswin_tiny/smallを用いた場合のconfig/log/weightが存在します。<br>
まずこれをそのまま使用したいのですがmaskの情報はvinbigには無いのでmask headを消しておきます。


```python
   roi_head=dict(
       type='StandardRoIHead',
       bbox_roi_extractor=dict(
           type='SingleRoIExtractor',
           roi_layer=dict(type='RoIAlign', output_size=7, sampling_ratio=0),
           out_channels=256,
           featmap_strides=[4, 8, 16, 32]),
       bbox_head=dict(
           type='Shared2FCBBoxHead',
           in_channels=256,
           fc_out_channels=1024,
           roi_feat_size=7,
           num_classes=14,
           bbox_coder=dict(
               type='DeltaXYWHBBoxCoder',
               target_means=[0., 0., 0., 0.],
               target_stds=[0.1, 0.1, 0.2, 0.2]),
           reg_class_agnostic=False,
           loss_cls=dict(
               type='CrossEntropyLoss', use_sigmoid=False, loss_weight=1.0),
           loss_bbox=dict(type='L1Loss', loss_weight=1.0)),
       mask_roi_extractor=dict(
           type='SingleRoIExtractor',
           roi_layer=dict(type='RoIAlign', output_size=14, sampling_ratio=0),
           out_channels=256,
           featmap_strides=[4, 8, 16, 32]),
       #mask_head=dict(
       #    type='FCNMaskHead',
       #    num_convs=4,
       #    in_channels=256,
       #    conv_out_channels=256,
       #    num_classes=14,
       #    loss_mask=dict(
       #        type='CrossEntropyLoss', use_mask=True, loss_weight=1.0))
       ),

```
次にbackboneを変更します。

```python
   backbone=dict(
       type='SwinTransformer',
       embed_dims=96,
       depths=[2, 2, 6, 2],
       num_heads=[3, 6, 12, 24],
       window_size=7,
       mlp_ratio=4,
       qkv_bias=True,
       qk_scale=None,
       drop_rate=0.,
       attn_drop_rate=0.,
       drop_path_rate=0.2,
       patch_norm=True,
       out_indices=(0, 1, 2, 3),
       with_cp=False,
       convert_weights=True,
       init_cfg=dict(type='Pretrained', checkpoint="https://github.com/SwinTransformer/storage/releases/download/v1.0.0/swin_tiny_patch4_window7_224.pth'")),
 
   neck=dict(
       type='FPN',
       in_channels=[96, 192, 384, 768],
       out_channels=256,
       num_outs=5),
```

load_from = ”~/mmdetection/checkpoints/mask_rcnn_swin-t-p4-w7_fpn_ms-crop-3x_coco_20210906_131725-bacf6f7b.pth”のように指定すれば事前学習済みの重みを使用できます。<br>
結果は以下です。

|    |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  tiny |  0.3600  |  0.223 |
|  small |  0.3560 |  0.221  |




#### 2. augmentationを追加

[albumentationsの使用例](https://github.com/open-mmlab/mmdetection/blob/master/configs/albu_example/mask_rcnn_r50_fpn_albu_1x_coco.py)があったのでそのままコピペします。SIIM2021の3rd解法のものを借用したものも試しました。

|    |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  small |  0.3560 |  0.221  |
|  small+aug |  0.3790  |  0.240  |
|  small+aug_3rd  |  0.3840  |  0.234  |


augmentationの追加は気兼ねなくできる＆伸び代もありそうなので良さそうです。

#### 3. さらにaugmentationを強める

mmdetetionでのmosaic/mixupの使用には[yoloxのconfig](https://github.com/open-mmlab/mmdetection/blob/master/configs/yolox/yolox_s_8x8_300e_coco.py)をみれば良さそうです。datasetのtypeをMultiImageMixDatasetにするなどの修正が必要です。



```python

img_scale = (1024,1024) # (1333, 800)
train_pipeline = [
    #dict(type='LoadImageFromFile'),
    #dict(type='LoadAnnotations', with_bbox=True),
    dict(type='Mosaic', img_scale=img_scale, pad_val=114.0),
    #dict(
    #    type='RandomAffine',
    #    scaling_ratio_range=(0.1, 2),
    #    border=(-img_scale[0] // 2, -img_scale[1] // 2)),
    dict(
        type='Albu',
        transforms=albu_train_transforms,
        bbox_params=dict(
            type='BboxParams',
            format='pascal_voc',
            label_fields=['gt_labels'],
            min_visibility=0.0,
            filter_lost_elements=True),
        keymap={
            'img': 'image',
            #'gt_masks': 'masks',
            'gt_bboxes': 'bboxes'
        },
        update_pad_shape=False,
        skip_img_without_anno=True),
    #dict(
    #    type='MixUp',
    #    img_scale=img_scale,
    #    ratio_range=(0.8, 1.6),
    #    pad_val=114.0),
    dict(type='Resize', img_scale=img_scale, keep_ratio=True),
    dict(type='RandomFlip', flip_ratio=0.5),
    dict(type='Normalize', **img_norm_cfg),
    dict(type='Pad', size_divisor=32),
    dict(type='DefaultFormatBundle'),
    dict(type='Collect', keys=['img', 'gt_bboxes', 'gt_labels']),
]

train_dataset = dict(
    type='MultiImageMixDataset',
    dataset=dict(
        type=dataset_type,
        classes=classes,
        ann_file= '~/vinbig/dataset/fold0/annotations/train.json',
        img_prefix= '~/vinbig/dataset/fold0/train2017/',
        pipeline=[
            dict(type='LoadImageFromFile'),
            dict(type='LoadAnnotations', with_bbox=True)
        ],
        filter_empty_gt=False,
    ),
    pipeline=train_pipeline)

data = dict(
    samples_per_gpu=8,
    workers_per_gpu=12,
    #train=dict(
    #    type=dataset_type,
    #    classes=classes,
    #    ann_file='~/vinbig/dataset/fold0/annotations/train.json',
    #    img_prefix='~/vinbig/dataset/fold0/train2017/',
    #    pipeline=train_pipeline),
    train = train_dataset,
    val=dict(
        type=dataset_type,
        classes=classes,
        ann_file='~/vinbig/dataset/fold0/annotations/val.json',
        img_prefix= "~/vinbig/dataset/fold0/val2017/",
        pipeline=test_pipeline),
```

|    |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  small+aug |  0.3790  |  0.240  |
|  small+mosaic |  0.3310  |  0.208  |


適当に使って良いものでは無さそうです。



#### 4. backboneを大きいモデルに変更する

[SIIM2021の3rd解法のレポジトリ](https://github.com/yujiariyasu/siim_covid19_detection/blob/main/ian-siim/mmdetection/configs/swin)にはモデルの大きさを大きくしてみるという実験を行った痕跡が残されているので倣ってみます。(swin_rsna000.py swin_rsna001.py)<br>backboneとneckの値をいじればなんとかなりそうです。以下はtiny→baseの変更をおこなったものです。


```python
   backbone=dict(
       type='SwinTransformer',
       #embed_dims=96,
       embed_dim=128,
       #embed_dim=192, #large
       #depths=[2, 2, 6, 2],
       depths=[2, 2, 18, 2], #largeもおなじ
       #num_heads=[3, 6, 12, 24],
       num_heads=[4, 8, 16, 32],
       #num_heads=[6, 12, 24, 48], #large
       #window_size=7,
       window_size=12,
       mlp_ratio=4,
       qkv_bias=True,
       qk_scale=None,
       drop_rate=0.,
       attn_drop_rate=0.,
       drop_path_rate=0.2,
       patch_norm=True,
       out_indices=(0, 1, 2, 3),
       with_cp=False,
       convert_weights=True,
       init_cfg=dict(type='Pretrained', checkpoint="https://github.com/SwinTransformer/storage/releases/download/v1.0.0/swin_base_patch4_window12_384_22kto1k.pth'")),
 
   neck=dict(
       type='FPN',
       #in_channels=[96, 192, 384, 768],
       in_channels=[128, 256, 512, 1024],
       #in_channels=[192, 384, 768, 1536],, #large
       out_channels=256,
       num_outs=5),
```

base large smallの結果は以下です。モデルを大きくすると収束が遅くなるのでは無いかと思ったのでlargeを長めに学習させるものも試しました。augmentationはSIIM2021の3rd解法のものを使用しています。


|  TH  |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  small |  0.3840 |  0.234  |
|  base |  0.3770  |  0.246 |
|  lagre  |  0.3330  |  0.230  |
|  lagre(3x)  |  0.3640  |  0.234  |


largeは収束が遅い&脳死で使って良いわけでは無さそうです。

#### 5. 他のモデルを試す

CascadeRcnn,retinanet,HTCもsubしてみます。各々backboneとneckを書き換えれば使用可能です。

|  TH  |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  Cascade |  0.3840 | 0.226 |
|  retinanet |  0.3790 |  0.232 |
|  HTC |  0.3640 | 0.229 |

使いこなせないですね...


#### 6. cocoの事前学習が有効か見てみる

mask R cnnでのcoco2017の事前学習済みの重みはmmdetectionに公開されていますが、他のモデルについては公式のmmdetectionには見当たりません。ネットを探しても見つからない場合、自分でcoco2017の事前学習をすることになります。A100１枚でswin smallのretinanetを12epoch学習させるのに2.5日ほどかかるので非常にしんどいです。<br>
事前学習でのmAPが高ければ高いほどよいのか確認するため8epoch終了時(学習途中)と12epoch終了時(最終epoch)の2つの重みを用いました。2つのcocoでのmAP_valは以下の通りです。


|  TH  |  val_mAP@0.5:0.95 |  val_mAP@0.5 |  val_mAP@0.75 |  val_mAP small |  val_mAP medium |  val_mAP large |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|  8epoch |  0.407 |  0.619  | 0.438 | 0.245 | 0.450 | 0.542 |
|  12epoch |   0.457  |  0.670 | 0.492 | 0.300 | 0.498 | 0.605 |

結果は以下です。

|  TH  |  val_mAP@0.5 |  private LB |
| ---- | ---- | ---- |
|  without pretraining |  0.3790 |  0.232  |
|  8epoch |  0.3820  | 0.239 |
|  12epoch  |  0.3940  |  0.238  |


cocoの事前学習を用いればLBが上昇しました。事前学習でのmAPが高ければ高いほどよいとも言いきれ無さそうです。


### まとめ


mmdetectionを使用してSwinTransformerを物体検出に適用しました。勘所は不明ですがmmdetectionを気兼ねなく使えるようになったのはよかったと思います。追試はモチベがなかなか出ない問題がありますが一日のsub数など気にせず実験のフィードバックが得られるので楽しかったです。


### おまけ

mmdetection以外でもswinで物体検出が可能です。

1. [yolov5](https://github.com/yl305237731/flexible-yolov5)

yolov5のbackboneなどを自分で変えて実験できます。coco2017で事前学習した重みは載っていません。[SIIM2021の2ndのレポジトリ](https://github.com/nvnnghia/siim2021/tree/main/detection)に似たようなものがありました。
<br>
これをvinbigに試しても散々でした。
nvnn氏曰く「SIIM2021の検出タスクは解像度が低くても&cocoの精度が低いものでもcv/lbがそんなに変わらなかったのでcocoで事前学習されていないモデルでも構わず使用した。有効かどうかはデータセット依存だよ」とのことです。”そんなに変わらない”の肌感は不明です｡｡｡

2. [detectron2](https://github.com/xiaohu2015/SwinT_detectron2)

detectron2でswin_Tinyのretinanetを36epoch事前学習させた重みなどがあります。詳細は未確認です。