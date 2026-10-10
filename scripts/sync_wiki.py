import os
import re

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIKI_DIR = os.path.join(BASE_DIR, "wiki")
MODEL_MD = os.path.join(WIKI_DIR, "model.md")
INDEX_HTML = os.path.join(WIKI_DIR, "index.html")

with open(MODEL_MD, 'r', encoding='utf-8') as f:
    md_content = f.read().strip()

with open(INDEX_HTML, 'r', encoding='utf-8') as f:
    html_content = f.read()

start_tag = '<script type="text/markdown" id="raw-markdown-content">\n'
end_tag = '\n</script>'

start_pos = html_content.find(start_tag)
if start_pos == -1:
    raise ValueError("Start tag not found in index.html")

end_pos = html_content.find(end_tag, start_pos)
if end_pos == -1:
    raise ValueError("End tag not found in index.html")

new_html = html_content[:start_pos + len(start_tag)] + md_content + html_content[end_pos:]

with open(INDEX_HTML, 'w', encoding='utf-8') as f:
    f.write(new_html)

print("wiki/index.html successfully updated and synchronized!")
