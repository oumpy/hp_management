#!/bin/sh
# Prune branch previews in the output repository (previews/refs/heads/...).
#
# Run from the root of the source repository, with the output repository
# cloned at ./output (tools/init_output.sh).  Deletions are staged in
# ./output but neither committed nor pushed: the caller commits, so that
# the preview and deploy workflows can fold the pruning into their own
# commit (one commit means one GitHub Pages build).
#
# Environment:
#   PRUNE_MODE      orphans  : delete previews whose branch no longer exists
#                              in the source repository (default)
#                   branches : delete the previews listed in PRUNE_BRANCHES
#                   all      : delete every preview
#                   none     : delete nothing (list only)
#   PRUNE_BRANCHES  branch names for mode "branches", whitespace-separated
#   PRUNE_KEEP      if set to a number N, additionally keep only the N most
#                   recently updated previews (any mode)
#   PRUNE_DRY_RUN   "true" to report what would be deleted without deleting
#   SOURCE_REPO_URL repository to list live branches from (URL or remote
#                   name; default "origin" of the current checkout, whose
#                   credentials are then used)
#
# A Markdown report (all previews, their size and last update, and what
# was deleted) is appended to $GITHUB_STEP_SUMMARY when that is set, and
# printed to stdout otherwise.

set -e

MODE=${PRUNE_MODE:-orphans}
KEEP=${PRUNE_KEEP:-}
DRY=${PRUNE_DRY_RUN:-false}
SRC=${SOURCE_REPO_URL:-origin}
OUT=./output
ROOT=previews/refs/heads

[ -d "$OUT/$ROOT" ] || { echo "No previews directory ($OUT/$ROOT)."; exit 0; }

report() { if [ -n "$GITHUB_STEP_SUMMARY" ]; then printf '%s\n' "$*" >> "$GITHUB_STEP_SUMMARY"; fi; printf '%s\n' "$*"; }

# --- inventory: a preview root is a directory that contains articles.html
cd "$OUT"
previews=$(find "$ROOT" -name articles.html -type f | sed "s#^$ROOT/##; s#/articles.html\$##" | sort)
cd ..

# --- live branches of the source repository
live=$(git ls-remote --heads "$SRC" 2>/dev/null | sed 's#.*refs/heads/##')

# --- collect facts: name, size (MB), last update (epoch), status
tmp=$(mktemp)
for p in $previews; do
  size=$(du -sm "$OUT/$ROOT/$p" | cut -f1)
  ts=$(cd "$OUT" && git log -1 --format=%ct -- "$ROOT/$p" 2>/dev/null || echo 0)
  [ -n "$ts" ] || ts=0
  if printf '%s\n' "$live" | grep -qx "$p"; then st=live; else st=orphan; fi
  printf '%s\t%s\t%s\t%s\n' "$ts" "$p" "$size" "$st" >> "$tmp"
done

# --- decide what to delete
delete=""
case "$MODE" in
  orphans)  delete=$(awk -F'\t' '$4=="orphan"{print $2}' "$tmp") ;;
  branches) for b in $PRUNE_BRANCHES; do
              if awk -F'\t' '{print $2}' "$tmp" | grep -qx "$b"; then delete="$delete $b"; else echo "No preview for branch '$b'."; fi
            done ;;
  all)      delete=$(awk -F'\t' '{print $2}' "$tmp") ;;
  none)     delete="" ;;
  *)        echo "Unknown PRUNE_MODE '$MODE'." >&2; exit 2 ;;
esac
if [ -n "$KEEP" ]; then
  # newest first; everything beyond the first KEEP entries goes as well
  extra=$(sort -t"$(printf '\t')" -k1,1nr "$tmp" | awk -F'\t' -v k="$KEEP" 'NR>k{print $2}')
  delete="$delete $extra"
fi
delete=$(printf '%s\n' $delete | sort -u)

# --- report
title="### Branch previews (mode: $MODE"
[ -n "$KEEP" ] && title="$title, keep newest $KEEP"
[ "$DRY" = true ] && title="$title, dry run"
report "$title)"
report ""
report "| preview | size | last update | status | action |"
report "|---|---:|---|---|---|"
sort -t"$(printf '\t')" -k1,1nr "$tmp" | while IFS="$(printf '\t')" read -r ts p size st; do
  when=$(date -u -d "@$ts" +%Y-%m-%d 2>/dev/null || date -u -r "$ts" +%Y-%m-%d 2>/dev/null || echo "?")
  if printf '%s\n' "$delete" | grep -qx "$p"; then act=delete; else act=keep; fi
  report "| \`$p\` | $size MB | $when | $st | $act |"
done
n=$(awk 'END{print NR}' "$tmp"); sz=$(awk -F'\t' '{s+=$3} END{print s+0}' "$tmp")
d=$(printf '%s\n' $delete | grep -c . || true)
report ""
report "$n previews, $sz MB in total; deleting $d."
rm -f "$tmp"

# --- delete
[ -n "$delete" ] || exit 0
if [ "$DRY" = true ]; then echo "Dry run: nothing deleted."; exit 0; fi
cd "$OUT"
for p in $delete; do
  rm -rf "$ROOT/$p"
  # remove now-empty parent directories (branch name prefixes)
  d=$(dirname "$ROOT/$p")
  while [ "$d" != "$ROOT" ] && [ -d "$d" ] && [ -z "$(ls -A "$d")" ]; do rmdir "$d"; d=$(dirname "$d"); done
done
git add -A previews
echo "Deleted (staged in $OUT): $(printf '%s ' $delete)"
