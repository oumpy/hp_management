#!/bin/sh
if [ "$1" != "-c" ]; then
    git clone https://github.com/oumpy/oumpy.github.io.git ./output
    git submodule update -i
    git submodule update --init pelican-ipynb
fi

themename="voidy-bootstrap"
cd themes
git submodule update --init "$themename"
cd ..
mkdir -p "theme/$themename"
cp -an "themes/$themename"/* "theme/$themename/"

# put webhook-cgi.
cgitemplate=tools/deploy.py
if [ -e "webhook.py" ]; then
    cgipath=`python3 -c "import os,sys; sys.path.append(os.curdir); import webhookconf; print(webhookconf.cgipath)"`
    if [ -n "$cgipath" ]; then
        if [ ${cgipath: -1} == "/" ]; then
            cgipath="$cgipath`basename "$cgitemplate"`"
        fi
        pythonpath=`which python3`
        hpmanagement_path=`pwd`
        mkdir -p `dirname "$cgipath"`
        sed -e "s|\$pythonpath|$pythonpath|g; s|\$hpmanagement_path|$hpmanagement_path|g" \
            < $cgitemplate \
            > "$cgipath"
        chmod +x "$cgipath"
    fi
fi
