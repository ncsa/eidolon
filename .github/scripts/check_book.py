#!/usr/bin/env python3
"""Guards for the mdBook site. mdBook's own exit code does not cover these.

Measured on mdbook 0.5.4:

* A SUMMARY.md entry pointing at a file that does not exist makes mdbook **create** the
  file and exit 0, publishing a blank chapter with no warning at all. So the SUMMARY is
  checked against the working tree BEFORE the build, while the file is still absent.
* A `{{#include}}` whose source is missing prints ERROR and WARN lines and still exits 0,
  embedding the error text in the page. The caller greps the build log for those.
* Nothing checks links. A chapter can move without its SUMMARY entry following, and a
  file included from `docs/` can carry a relative link that resolved there and not here.

Usage:
    check_book.py summary <book-dir>    # before `mdbook build`
    check_book.py built   <build-dir>   # after
"""
import os
import re
import sys

# A page shorter than this has no content worth publishing; the blank chapter mdbook
# creates for a missing SUMMARY entry is empty, and a stub that renders to nothing is a
# hole in the site either way.
MIN_PAGE_CHARS = 40


def fail(msg):
    print(f"::error::{msg}")


def check_summary(book_dir):
    src = os.path.join(book_dir, "src")
    summary = os.path.join(src, "SUMMARY.md")
    if not os.path.exists(summary):
        sys.exit(f"no SUMMARY.md at {summary}")
    entries = re.findall(r"\]\(([^)#]+\.md)\)", open(summary, encoding="utf-8").read())
    if not entries:
        sys.exit(f"{summary} lists no chapters; refusing to build an empty book")
    bad = 0
    for rel in entries:
        path = os.path.join(src, rel)
        if not os.path.exists(path):
            fail(f"SUMMARY.md lists {rel}, which does not exist. mdbook would create it "
                 f"and publish a blank chapter.")
            bad += 1
        elif os.path.getsize(path) == 0:
            fail(f"SUMMARY.md chapter {rel} is empty")
            bad += 1
    print(f"checked {len(entries)} SUMMARY entries, {bad} missing or empty")
    return bad


def check_built(build_dir):
    broken, thin, pages = [], [], 0
    for dirpath, _, files in os.walk(build_dir):
        for name in files:
            if not name.endswith(".html"):
                continue
            path = os.path.join(dirpath, name)
            html = open(path, encoding="utf-8").read()
            body = re.search(r"<main>(.*?)</main>", html, re.S)
            if not body:
                continue
            rel = os.path.relpath(path, build_dir)
            pages += 1
            text = re.sub(r"<[^>]+>", "", body.group(1)).strip()
            if len(text) < MIN_PAGE_CHARS and rel not in ("toc.html",):
                thin.append(f"{rel} ({len(text)} chars)")
            for href in re.findall(r'href="([^"]+)"', body.group(1)):
                if href.startswith(("http://", "https://", "#", "mailto:")):
                    continue
                target = os.path.normpath(os.path.join(dirpath, href.split("#")[0]))
                if not os.path.exists(target):
                    broken.append(f"{rel} -> {href}")
    if pages == 0:
        sys.exit(f"no pages found under {build_dir}: the build produced nothing to check")
    for b in broken:
        fail(f"broken internal link: {b}")
    for t in thin:
        fail(f"page has essentially no content: {t}")
    print(f"checked {pages} pages, {len(broken)} broken links, {len(thin)} empty pages")
    return len(broken) + len(thin)


mode, target = sys.argv[1], sys.argv[2]
sys.exit(1 if (check_summary(target) if mode == "summary" else check_built(target)) else 0)
