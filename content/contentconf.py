# -*- coding: utf-8 -*- #
# site-specific settings
import sys, os
sys.path.append(os.curdir + '/..')
from tools.lib.pelicanns import *
for pluginpath in PLUGIN_PATHS:
    sys.path.append(os.curdir + '/' + pluginpath)
from makemenu import MenuItem
import re
striptags = lambda html: re.compile(r'<[^>]*?>').sub('', html)

LANG = 'ja'

GITHUB_ACCOUNT = 'oumpy'
SOURCEREPOSITORY_NAME = 'hp_management'

# Author
AUTHOR = 'Python会'
SITENAME = '大阪大学医学部 Python会 <span style="font-size:smaller;">(情報医科学研究会)</span>'
SITETAG = striptags(SITENAME)
SITEURL = ''
SITESUBTITLE ='Now is better than never.'
COPYRIGHT_YEAR = datetime.date.today().year
COPYRIGHT_AUTHOR = '大阪大学医学部 情報医科学研究会 (Python会)'

LOGOIMG = 'images/logo.jpg'
# TOP_LOGO_IMAGE = LOGOIMG
OPEN_GRAPH_IMAGE = LOGOIMG

# Social
SOCIAL = [ # (name, URL, icon, color, size)
    # ('facebook', '#3B5998'),
    # ('技術Blog (はてな)','https://oumedpython.hatenablog.com/'),
    ('Twitter', 'https://twitter.com/oumed_python', '<i class="fab fa-twitter"></i>', '#1DA1F2', 'larger'),
    ('E-mail', 'mailto:handai.python@gmail.com', '<i class="far fa-envelope"></i>', '#0078D4', 'larger'),
    ('GitHub Organization', 'https://github.com/oumpy', '<i class="fab fa-github"></i>', '#211F1F', 'larger'),
    ('YouTube Channel', 'https://www.youtube.com/channel/UCh1eAeDCpsZeOh0Z9paNfHQ', '<i class="fab fa-youtube"></i>', '#c4302b', 'larger'),
    ('Connpass', 'https://oum-python.connpass.com', '<img width="35px" src="https://connpass.com/static/img/72_72.png" style="display: inline;"/>', '#000000', 'normal'),
    ('Atom Feed', './feeds/all.atom.xml', '<i class="fa fa-rss fa-fw fa-lg"></i>', '#00008b', 'smaller'),
    ('RSS Feed', './feeds/all.rss.xml', '<i class="fas fa-rss-square fa-fw fa-lg"></i>', '#f26522', 'normal'),
]

# Blogroll
LINKS = [
        # ('Python会ブログ','https://oumedpython.hatenablog.com/'),
        #('大阪大学医学部', 'http://www.med.osaka-u.ac.jp/'),
        #('Python.org', 'http://python.org/'),
        #('Pelican（本サイトで使用）', 'http://getpelican.com/'),
        # ('Jinja2', 'http://jinja.pocoo.org/'),
        # ('You can modify those links in your config file', '#'),
]

# How many recent posts appear in sidebar?
RECENT_POST_COUNT = 8

# Settings for Google Custom Search
GOOGLE_CSE_ID = '012109292501676780101:lytlaxeswfy'
SEARCHBOX_MESSAGE = 'サイト内を検索'

# Settings for Twitter Timeline
TWITTER_TIMELINE_URL = "https://twitter.com/oumed_python?ref_src=twsrc%5Etfw"
TWITTER_USERNAME = 'oumed_python'

CUSTOM_SOCIAL_TITLE = "ソーシャル"
CUSTOM_RECENTPOSTS_TITLE = '新着記事'
CUSTOM_CATEGORIES_TITLE = "記事カテゴリ"
CUSTOM_TAGS_TITLE = "タグ"
CUSTOM_LINKS_TITLE = "リンク"
CUSTOM_TWITTERTL_TITLE = "Timeline"
CUSTOM_RELATED_ARTICLES_TITLE = "関連記事"
CUSTOM_ARCHIVES_TITLE = "過去記事一覧"
CUSTOM_AUTHORS_TITLE = "著者一覧"
CUSTOM_TAGS_TITLE = "タグ一覧"
CUSTOM_SIDEBAR_SUPPORT_TITLE = "会へのご支援"

SIDEBAR_SUPPORTLINK_LIST = [
    ('https://paypal.me/oumpy/',
     '<span style="font-size:larger;"><span class="text-primary" style="margin-right:5px;"><i class="fab fa-paypal"></i></span><font color="#003087">Pay</font><font color="#009cde">Pal</font><font color="#012169">.Me</font></span>',
     '<br><span style="font-size:smaller; margin-left:0px;">(<a href="{{ SITEURL }}/donations.html">ご寄付を頂いた方々</a>)</span>',
    ),
    # '<img src="https://www.paypalobjects.com/webstatic/paypalme/images/social/pplogo384.png"/>'
]

DISPLAY_PAGES_ON_MENU = False
CATEGORYNAMES_ALTERNATIVES = {
    'news': ('お知らせ', 'Python会からのお知らせ'),
    'blog': ('技術ブログ',),
}

ADD_ON_MENU = [
    MenuItem('index.html',
             self_in_subsections=True,
             subsections=[
                 MenuItem('constitution.html', title='会規約', subsections=MenuItem.AUTO),
             ]),
    MenuItem('activities.html', self_in_subsections=True, subsections=MenuItem.AUTO),
    MenuItem('achievements.html', self_in_subsections=True, subsections=MenuItem.AUTO),
    MenuItem('blog.html', title=CATEGORYNAMES_ALTERNATIVES['blog'][0],
             self_in_subsections=True,
             active_pages=r'(blog/|tag/(?!news).+\.html$|author/(?!pythonhui).+\.html$)',
             subsections = [
                 # MenuItem('tags.html', title='タグ一覧'),
                 MenuItem('authors.html', title='著者一覧'),
                 MenuItem('archives.html', title='過去記事一覧'),
             ]),
    MenuItem('news.html', title=CATEGORYNAMES_ALTERNATIVES['news'][0],
             self_in_subsections=True,
             active_pages=r'(news/|tag/news\.html$|author/pythonhui\.html$)'),
    MenuItem('recruit.html', self_in_subsections=True, subsections=MenuItem.AUTO),
    MenuItem('student_server.html',
             title='学生用計算機',
             self_in_subsections=True, subsections=MenuItem.AUTO),
    MenuItem('contact.html', self_in_subsections=True),
]
MENU_STEPS = 1

HIDE_ARCHIVES_ON_MENU = True
SHOW_FEED_ATOM_ON_MENU = SHOW_FEED_RSS_ON_MENU = False
SIDEBAR_HIDE_CATEGORIES = True

SHOW_CATEGORIES_ON_LIST = False
SHOW_CATEGORY_TITLE = True
def readfile(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    return ''.join(content)
PAGE_EXCLUDES = ['pages/includes']
CATEGORY_CONTENTS = {
    'blog' : readfile('content/pages/includes/blog_content.html'),
    'news' : readfile('content/pages/includes/news_content.html'),
}

CUSTOM_TAG_BADGE_COLOR = 'blue'
TAG_GROUPS = [ # (groupname, [articles,...,], badge_color )
    ('Research tools & topics',
     ['Bioinformatics', 'Machine Learning', 'Statistics', 'Data Science Competition', 'ハードウェア',
      'Simulation', 'Neuroscience'],
     'darkorange'),
    ('Programming', ['Python', 'Shell script', '競技プログラミング', 'GitHub', 'Unix', '自動化'], 'green'),
    ('その他', ['論文まとめ', '論文関連', '検定試験', '海外留学', '勉強会', '勉学支援'], CUSTOM_TAG_BADGE_COLOR),
]
CUSTOM_TAG_BADGE_COLORS = {'News' : 'hotpink'}
for group in TAG_GROUPS:
    for tag in group[1]:
        CUSTOM_TAG_BADGE_COLORS[tag] = group[2]
TAG_CLOUD_BADGE = True

PREVIEW_SITENAME_APPEND = ' (テスト用ページ)'

PLUGINS += [
    'postprocess',
]
POSTPROCESS_COMMAND = 'sh content/postprocess.sh'
