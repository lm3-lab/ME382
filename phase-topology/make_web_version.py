#!/usr/bin/env python3
"""Derive the publishable Artifact fragment from index.html.

index.html is a complete standalone document so it opens from file:// on any
machine. The Artifact platform supplies its own <!doctype>/<head>/<body>
skeleton, so the published version must be a fragment: <title>, the font
<link>s and <style> at the top, then the body content. Keeping one source of
truth means the two versions cannot drift.

    python3 make_web_version.py [out.html]
"""
import re, sys, pathlib

src = pathlib.Path(__file__).with_name("index.html").read_text()
out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "web.html")

head = re.search(r"<head>(.*?)</head>", src, re.S)
body = re.search(r"<body>(.*?)</body>", src, re.S)
if not head or not body:
    sys.exit("could not find <head> / <body> in index.html")

# Keep title, font links and styles. Drop charset/viewport: the skeleton sets
# both, and its viewport carries viewport-fit=cover for the safe-area padding.
keep = []
for m in re.finditer(r"<title>.*?</title>|<link\b[^>]*>|<style>.*?</style>", head.group(1), re.S):
    keep.append(m.group(0))

frag = "\n".join(keep) + "\n" + body.group(1).strip() + "\n"
# \b so <header> does not trip the <head> check
stray = re.search(r"<\s*(!doctype|html|head|body)\b", frag, re.I)
if stray:
    sys.exit("fragment still contains %s" % stray.group(0))
out.write_text(frag)
print("wrote %s  (%d bytes, %d lines)" % (out, len(frag), frag.count("\n") + 1))
