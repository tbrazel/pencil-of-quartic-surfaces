"""Build the static site: inline code/data/pencil.json into code/template.html -> site/index.html."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SITE = HERE.parent / "site"
SITE.mkdir(exist_ok=True)
template = (HERE / "template.html").read_text(encoding="utf-8")
data = (HERE / "data" / "pencil.json").read_text(encoding="utf-8").strip()
assert "__DATA__" in template
(SITE / "index.html").write_text(template.replace("__DATA__", data), encoding="utf-8")
print(f"wrote {SITE / 'index.html'}")
