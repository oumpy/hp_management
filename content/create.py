#!/usr/bin/env python3
from datetime import datetime
from collections import defaultdict
import os
import argparse
import shutil

template = "./template.md"

today = datetime.now().strftime("%Y.%m.%d")
fields = [  # (field, default, comment)
    ('slug',      'my_new_article',   'ファイル名。英数字、拡張子なし'),
    ('date',      today,              '記事日付。YYYY.mm.DD または YYYY.mm.DD HH.MM 形式'),
    ('modified',  None,               '更新日。None の場合は date と同じ'),
    ('title',     None,               '記事タイトル。content 内に # 要素がある場合は指定しない。'),
    ('category',  'blog',             'blog, news または page'),
    ('tags',      None,               'Bioinformatics など、サイト参照。複数の場合カンマ区切り'),
    ('author',    None,               '記事執筆者の名前。姓またはHN、複数の場合カンマ区切り'),
    ('content',   None,               '記事内容ファイル (.md/.ipynb) がすでにある場合はパスを指定'),
]

if __name__ == '__main__':
    # Read command-line arguments
    parser = argparse.ArgumentParser()
    shortused = defaultdict(bool)
    for field, default, comment in fields:
        if shortused[field[0].lower()]:
            parser.add_argument('--{}'.format(field.lower()), default=None, help=comment)
        else:
            parser.add_argument('--{}'.format(field.lower()), '-{}'.format(field[0].lower()), default=None, help=comment)
            shortused[field[0].lower()] = True
    parser.add_argument('--filename', '-f', default=None, help='--slug と同じ。両方指定すると slug が優先')
    parser.add_argument('--default', action='store_true', default=None, help='未指定の項目を質問せずデフォルト値にする')

    args = parser.parse_args()
    args.slug = args.slug or args.filename

    print('新しい記事ファイルを作成します。', end=' ')
    if not args.default:
        print('必要項目を順に入力してください。')
        print('[ ] 内がデフォルト値です。', end=' ')
        print('各項目はファイル作成後に手動で編集することもできます。')
    print('---')

    # Input information
    values = dict()
    for field, default, comment in fields:
        Field = field.capitalize()
        if args.default or getattr(args, field):
            print('{0}: {1}'.format(Field, getattr(args, field) or fields[field][1]))
            values[field] = args[field]
        else:
            print('{0} ({2}) [ {1} ]: '.format(Field, default, comment), end='')
            given = input().strip()
            if given:
                values[field] = given
            else:
                values[field] = default
    print('---')

    # Set school year of Japan
    year, month = map(int, values['date'].split('.')[:2])
    if month < 4:
        schoolyear = year - 1
    else:
        schoolyear = year

    # Set directory to save
    if values['category'].lower() in {'page', 'pages'}:
        directory = './pages'
    else:
        directory = "./articles/{0}sy/{1}".format(schoolyear, values['category'])
    filepath = "{0}/{1}.md".format(directory, values['slug'])
    imagedirpath = "{0}/images/{1}_figs".format(directory, values['slug'])

    # set modified date
    if not values['modified']:
        values['modified'] = values['date']

    # read content
    content = []
    if values['content']:
        contentfile = values['content']
        if os.path.isfile(contentfile):
            suffix = contentfile.split('.')[-1]
            if suffix == 'md':
                with open(contentfile, 'r') as cf:
                    for line in cf.readlines():
                        line = line.rstrip('\n')
                        if not values['title']:
                            first, second = line.lstrip()[:2]
                            if first == '#' and second != '#':
                                values['title'] = line.strip()[1:].lstrip()
                                continue
                        content.append(line)
            elif suffix == 'ipynb':
                filepath_body = filepath.rstrip('.md')
                filepath = filepath_body + '.nbdata'
                filepath_ipynb = filepath_body + '.ipynb'
                if os.path.exists(filepath_ipynb):
                    print("'{}' はすでに存在します。".format(filepath))
                else:
                    shutil.copy(contentfile, filepath_ipynb)
                    print("ファイル '{}' を作成 (コピー) しました。".format(filepath_ipynb))
            else:
                print("ファイル '{}' は対応するファイル形式ではありません。指定を無視します。".format(contentfile))
        else:
            print("ファイル '{}' は存在しません。指定を無視します。".format(contentfile))

    if os.path.exists(filepath):
    # file already exists ?
        print("Error: '{}' はすでに存在します。".format(filepath))
        exit()

    # make the file
    if not os.path.exists(directory):
        os.makedirs(directory)
    suffix = filepath.split('.')[-1]
    with open(filepath, 'w') as tf:
        with open(template, 'r') as sf:
            for line in sf.readlines():
                line = line.rstrip('\n')
                Field = line.split(':')[0]
                field = Field.lower()
                if field in values.keys() and values[field]:
                    if field in {'author', 'tags'}:
                        values[field] = values[field].strip().strip(',')
                    if field == 'author' and len(values[field].split(',')) > 1:
                        print('{0}s: {1}'.format(Field, values[field]), file=tf)
                    else:
                        print('{0}: {1}'.format(Field, values[field]), file=tf)
                elif suffix == '.md' or line.strip():
                    print(line, file=tf)
        for line in content:
            print(line, file=tf)

    print("ファイル '{}' を作成しました。".format(filepath))
    if suffix == '.md':
        try:
            os.makedirs(imagedirpath)
            print("画像保存用ディレクトリ'{}'を作成しました。".format(imagedirpath))
        except FileExistsError:
            if os.path.isdir(imagedirpath):
                print("画像保存用ディレクトリ'{}'はすでに存在します。".format(imagedirpath))
            else:
                print("画像保存用ディレクトリ'{}'を作成できませんでした。".format(imagedirpath))
