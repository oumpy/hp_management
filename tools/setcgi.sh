#!/bin/sh
# put webhook-cgi.
cgitemplate="tools/lib/deploy.py"
if [ -e "webhookconf.py" ]; then
    cgipath=`python3 -c "import os,sys; sys.path.append(os.curdir); import webhookconf; print(webhookconf.cgipath)"`
    if [ -n "$cgipath" ]; then
        if [ `echo "$cgipath" | rev | cut -c 1` == "/" ]; then
            cgipath="$cgipath`basename "$cgitemplate"`"
        fi
        pythonpath=`which python3`
        hpmanagement_path=`pwd`
        path="$PATH"
        mkdir -p `dirname "$cgipath"`
        cat "$cgitemplate" | \
        sed -e "s|\$pythonpath|$pythonpath|g" | \
        sed -e "s|\$hpmanagement_path|$hpmanagement_path|g" | \
        sed -e "s|\$path|$path|g" \
            > "$cgipath"
        chmod +x "$cgipath"
    fi
fi
