# -*- coding: utf-8 -*- #
# site-specific settings
GITHUB_ACCOUNT = 'oumpy'
SOURCEREPOSITORY_NAME = 'hp_management'

# Author
AUTHOR = 'Python会'
SITENAME = '大阪大学医学部 Python会'
SITEURL = ''
AUTHOR_INTRO = u'大阪大学医学部所属のPython職人集団です'
AUTHOR_DESCRIPTION = u'Now is better than never'
# AUTHOR_AVATAR = '../images/'
# AUTHOR_WEB = 'https://twitter.com/oumed_python'
SITESUBTITLE ='Now is better than never.'
SITETAG = SITENAME

# Social
SOCIAL = (
    # ('facebook', ''),
    # ('技術Blog (はてな)','https://oumedpython.hatenablog.com/'),
    ('Twitter', 'https://twitter.com/oumed_python'),
    ('E-mail', 'mailto:handai.python@gmail.com'),
    ('GitHub', 'https://github.com/oumpy'),
    ('YouTube', 'https://www.youtube.com/channel/UCh1eAeDCpsZeOh0Z9paNfHQ'),
)

# Blogroll
LINKS = (
        # ('Python会ブログ','https://oumedpython.hatenablog.com/'),
        #('大阪大学医学部', 'http://www.med.osaka-u.ac.jp/'),
        #('Python.org', 'http://python.org/'),
        #('Pelican（本サイトで使用）', 'http://getpelican.com/'),
        # ('Jinja2', 'http://jinja.pocoo.org/'),
        # ('You can modify those links in your config file', '#'),
        )

# Settings for Twitter Timeline
TWITTER_TIMELINE_URL = "https://twitter.com/oumed_python?ref_src=twsrc%5Etfw"
TWITTER_USERNAME = 'oumed_python'

CUSTOM_SOCIAL_TITLE = "ソーシャル"
CUSTOM_CATEGORIES_TITLE = "記事カテゴリ"
CUSTOM_TAGS_TITLE = "タグ"
CUSTOM_LINKS_TITLE = "リンク"
CUSTOM_TWITTERTL_TITLE = "Timeline"
CUSTOM_RELATED_ARTICLES_TITLE = "関連記事"

OPEN_GRAPH_IMAGE = 'logo.jpg'

DISPLAY_PAGES_ON_MENU = False
blog_title = '技術ブログ'
news_title = 'Python会からのお知らせ'
ADD_ON_MENU = (
    ('Python会について', 'index.html'),
    ('活動内容', 'activities.html'),
    ('実績', 'achievements.html'),
    (blog_title, 'blog.html'),
    ('お知らせ', 'news.html'),
    ('会員募集', 'recruit.html'),
    ('Contact', 'contact.html'),
)
HIDE_ARCHIVES_ON_MENU = True
SIDEBAR_HIDE_CATEGORIES = True

SHOW_CATEGORY_TITLE = True
CUSTOM_CATEGORY_TITLES = {'Blog': blog_title, 'News': news_title}
def readfile(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    return ''.join(content)
PAGE_EXCLUDES = ['pages/includes']
CATEGORY_CONTENTS = {
    'Blog' : readfile('content/pages/includes/blog_content.html'),
    'News' : readfile('content/pages/includes/news_content.html'),
}
DEFAULT_PAGINATION = 10

TAG_GROUPS = [
    ('Research tools & techniques', ['Bioinformatics', 'Machine Learning', 'Statistics', 'Data Science Competition']),
    ('Programming', ['Python', 'Shell script', 'GitHub', '競技プログラミング']),
    ('その他', ['論文関連', '検定試験', '海外留学']),
]

PREVIEW_SITENAME_APPEND = ' (テスト用ページ)'
