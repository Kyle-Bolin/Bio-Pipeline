"""The browser front end: its built files ship with the package and it only calls tools that exist."""

import re
from pathlib import Path

import pytest

from cellfactory import tools
from cellfactory.web import STATIC_DIR, TOOLS_USED

FRONTEND_SRC = Path(__file__).resolve().parents[1] / "frontend" / "src"


def test_built_files_present():
    assert (STATIC_DIR / "index.html").is_file()
    html = (STATIC_DIR / "index.html").read_text()
    assets = re.findall(r'(?:src|href)="\./(assets/[^"]+)"', html)
    assert any(a.endswith(".js") for a in assets), "index.html doesn't load the app bundle"
    for a in assets:
        assert (STATIC_DIR / a).is_file(), a


def test_page_is_self_contained():
    html = (STATIC_DIR / "index.html").read_text()
    assert not re.search(r'(src|href)="https?://', html), "the page must work offline"


def test_tools_used_exist():
    assert set(TOOLS_USED) <= {f.__name__ for f in tools.TOOLS}


@pytest.mark.skipif(not FRONTEND_SRC.is_dir(), reason="front end source not in this checkout")
def test_source_only_calls_listed_tools():
    names = {f.__name__ for f in tools.TOOLS}
    api_ts = (FRONTEND_SRC / "api.ts").read_text()
    listed = re.search(r"TOOL_NAMES = \[([^\]]*)\]", api_ts)
    assert listed and set(re.findall(r'"([a-z_]+)"', listed.group(1))) == set(TOOLS_USED)
    # The SLOW set names tools by string; a typo there would silently skip the job queue.
    slow = re.search(r"SLOW = new Set<string>\(\[([^\]]*)\]\)", api_ts)
    assert slow and set(re.findall(r'"([a-z_]+)"', slow.group(1))) <= names
    called = set()
    for f in FRONTEND_SRC.rglob("*.ts*"):
        text = f.read_text()
        called |= set(re.findall(r'\b(?:timed|call)(?:<[^>(]*>)?\(\s*"([a-z_]+)"', text))
        called |= set(re.findall(r'"/tools/([a-z_]+)"', text))
    assert called, "no tool calls found in frontend/src"
    assert called <= set(TOOLS_USED), called - set(TOOLS_USED)
    # The All tools tab builds its tool names from GET /api/tools, so those exist by construction.
