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
from html.parser import HTMLParser

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


class _Text(HTMLParser):
    """Reader-visible text. A real parser, not a tag-stripping regex: a bare '<'
    in prose (as in '0 <= x <= 1') makes a regex swallow the rest of the line as
    if it were a tag, which would hide exactly the content this check exists to
    protect. Entities are decoded, so &lt; and &middot; compare as what is read."""
    SKIP = {"script", "style", "title"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts, self.depth = [], 0

    def handle_starttag(self, tag, attrs):
        if tag in self.SKIP:
            self.depth += 1

    def handle_endtag(self, tag):
        if tag in self.SKIP and self.depth:
            self.depth -= 1

    def handle_data(self, data):
        if not self.depth:
            self.parts.append(data)


def visible(h):
    t = _Text()
    t.feed(h)
    return re.sub(r"\s+", " ", "".join(t.parts)).strip()


# Guard against a silent content loss: every word visible in index.html must
# survive into the fragment. A regex that quietly eats a paragraph would
# otherwise publish a page missing text nobody notices until a student reads it.
want = visible(body.group(1))
got = visible(frag)
if want != got:
    n = next((i for i, (a, b) in enumerate(zip(want, got)) if a != b), min(len(want), len(got)))
    sys.exit("visible text differs from index.html at character %d\n"
             "  index.html: %r\n  generated : %r" % (n, want[n - 70:n + 70], got[n - 70:n + 70]))

for tag in ("style", "script"):
    a, b = len(re.findall(r"<%s\b" % tag, src)), len(re.findall(r"<%s\b" % tag, frag))
    if a != b:
        sys.exit("lost a <%s> block: index.html has %d, generated has %d" % (tag, a, b))

out.write_text(frag)
print("wrote %s  (%d bytes, %d lines)" % (out, len(frag), frag.count("\n") + 1))
print("checked: visible text identical to index.html (%d chars)" % len(want))
