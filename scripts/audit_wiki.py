import os
import re
import difflib

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIKI_DIR = os.path.join(BASE_DIR, "wiki")
MODEL_MD = os.path.join(WIKI_DIR, "model.md")
INDEX_HTML = os.path.join(WIKI_DIR, "index.html")
MODEL_TSX = os.path.join(WIKI_DIR, "src_gitlab", "model.tsx")

with open(MODEL_MD, 'r', encoding='utf-8') as f:
    md = f.read().strip()

with open(INDEX_HTML, 'r', encoding='utf-8') as f:
    html = f.read()

m = re.search(r'<script type="text/markdown" id="raw-markdown-content">([\s\S]*?)</script>', html)
if m:
    emb_md = m.group(1).strip()
    print(f"Length model.md: {len(md)}, Length embedded: {len(emb_md)}")
    if md == emb_md:
        print("EXACT MATCH between model.md and raw-markdown-content in index.html!")
    else:
        print("DIFFERENCES DETECTED between model.md and embedded markdown:")
        diff = list(difflib.unified_diff(md.splitlines(), emb_md.splitlines(), lineterm='', fromfile='model.md', tofile='embedded_in_index.html'))
        for line in diff[:40]:
            print(line)
else:
    print("Could not find raw-markdown-content in index.html")
