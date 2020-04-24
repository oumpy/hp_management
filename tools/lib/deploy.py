#!$pythonpath
# -*- coding: utf-8 -*-
import os, sys
import subprocess
import json
import hashlib, hmac
from datetime import datetime
import fasteners

hpmanagement_path = "$hpmanagement_path"
sys.path.append(hpmanagement_path)
from webhookconf import secret
webhooklogfilepath = hpmanagement_path + '/webhook.log'
updatelogfilepath = hpmanagement_path + '/update.log' 
lockfilepath = hpmanagement_path + '/deploycgi.lock'

os.environ['PATH'] = '$path'

print('Content-Type:text/html\n\n')

data = sys.stdin.read()
payload = json.loads(data)
signature_header = 'HTTP_X_HUB_SIGNATURE'
if signature_header in os.environ.keys():
    signature = os.environ[signature_header]
else:
    signature = ''

event_header = 'HTTP_X_GITHUB_EVENT'
if event_header in os.environ.keys():
    event = os.environ[event_header]
else:
    print('No event was given.')
    exit()

if secret:
    hasher = 'sha1=' + hmac.new(secret, bytes(data, 'utf-8'), hashlib.sha1).hexdigest()
    # print('Signature : {}'.format(signature))
    # print('Calculated: {}'.format(hasher))
    if signature != hasher:
        print('Signature could not be verified.')
        exit()
    else:
        print('Signature verified.')

branch = payload['ref']

lock = fasteners.InterProcessLock(lockfilepath)
lock.acquire()
f = open(webhooklogfilepath, 'a')
uf = open(updatelogfilepath, 'a')

now = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.now())
print(now, 'Branch', branch, 'pushed: ', end='', file=f)
print(now, file=uf)
os.chdir(hpmanagement_path)

if not event in {'push','create','delete'}:
    print('nothing done.', file=f)
else:
    if event == 'push':
        bs = branch.split('/')
        if len(bs) == 3 and bs[0] == 'refs' and bs[1] == 'heads':
            branch = bs[2]
        else:
            branch = ''
    if branch:
        subprocess.call([
            'sh', 'tools/updatesite.sh',
            'Auto-push by hp_management/%s update.' % branch,
            branch],
            stdout=uf, stderr=uf)
        print('website compiled and (possibly) pushed.', file=f)
    else:
        print('nothing done.', file=f)

f.close()
uf.close()
lock.release()
