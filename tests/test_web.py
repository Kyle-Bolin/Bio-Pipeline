"""The browser front end: its files ship with the package and it only calls tools that exist."""

import re

from cellfactory import tools
from cellfactory.web import STATIC_DIR, TOOLS_USED


def test_static_files_present():
    for name in ("index.html", "app.js", "style.css"):
        assert (STATIC_DIR / name).is_file(), name


def test_page_only_calls_existing_tools():
    names = {f.__name__ for f in tools.TOOLS}
    assert set(TOOLS_USED) <= names
    js = (STATIC_DIR / "app.js").read_text()
    called = set(re.findall(r'(?:timed|call)\("([a-z_]+)"', js)) | set(re.findall(r'"/tools/([a-z_]+)"', js))
    assert called, "no tool calls found in app.js"
    assert called <= set(TOOLS_USED), called - set(TOOLS_USED)


def test_page_is_self_contained():
    html = (STATIC_DIR / "index.html").read_text()
    assert not re.search(r'(src|href)="https?://', html), "the page must work offline"
