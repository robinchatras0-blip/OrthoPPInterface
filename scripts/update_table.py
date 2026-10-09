#!/usr/bin/env python3
"""
update_table.py - Updates the candidate evaluation table in wiki/model.md,
wiki/index.html, and wiki/src_gitlab/model.tsx with the streamlined, mathematically
coherent 8-column layout requested by the user:
- Pair (Pair 01 to Pair 19)
- iPTM(s1) (A' · B')
- iPTM(s2) (A' · B_WT)
- iPTM(s3) (A_WT · B')
- iPTM(s4) (A_WT · B_WT, reference ceiling = 0.545)
- Delta iPTM (Selectivity margin)
- DockQ(A' · B') (Design pose)
- DockQ(A · B) (Wild-type reference = 1.000)
"""

import os
import re

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_MD = os.path.join(BASE_DIR, "wiki", "model.md")
INDEX_HTML = os.path.join(BASE_DIR, "wiki", "index.html")
MODEL_TSX = os.path.join(BASE_DIR, "wiki", "src_gitlab", "model.tsx")

table_md = """### 4.1 Evaluation Table of the 19 Champions: 4-Way iPTM Matrix & CAPRI DockQ

The table below summarizes the core biophysical evaluation of the 19 strict orthogonal champions. Each column strictly corresponds to our state indices from Section 2.1 ($s_1$ to $s_4$), the resulting selectivity margin $\\Delta\\text{iPTM}$, and the CAPRI binding pose fidelity $\\text{DockQ}(A' \\cdot B')$ compared to native reference $\\text{DockQ}(A \\cdot B) = 1.000$:

| Pair | $\\text{iPTM}(s_1)$<br>$(A' \\cdot B')$ | $\\text{iPTM}(s_2)$<br>$(A' \\cdot B_{\\text{WT}})$ | $\\text{iPTM}(s_3)$<br>$(A_{\\text{WT}} \\cdot B')$ | $\\text{iPTM}(s_4)$<br>$(A_{\\text{WT}} \\cdot B_{\\text{WT}})$ | $\\Delta\\text{iPTM}$ | $\\text{DockQ}(A' \\cdot B')$ | $\\text{DockQ}(A \\cdot B)$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Pair 01** | **0.745** | 0.295 | 0.251 | 0.545 | <strong style="color: #2ed578;">+0.494</strong> | 0.481 | 1.000 |
| **Pair 02** | **0.781** | 0.237 | 0.299 | 0.545 | <strong style="color: #2ed578;">+0.483</strong> | 0.547 | 1.000 |
| **Pair 03** | **0.701** | 0.295 | 0.259 | 0.545 | <strong style="color: #2ed578;">+0.442</strong> | 0.582 | 1.000 |
| **Pair 04** | 0.602 | 0.246 | 0.221 | 0.545 | <strong style="color: #2ed578;">+0.381</strong> | 0.497 | 1.000 |
| **Pair 05** | 0.626 | 0.259 | 0.247 | 0.545 | <strong style="color: #2ed578;">+0.379</strong> | 0.507 | 1.000 |
| **Pair 06** | **0.741** | 0.295 | 0.384 | 0.545 | <strong style="color: #2ed578;">+0.357</strong> | 0.540 | 1.000 |
| **Pair 07** | 0.600 | 0.259 | 0.243 | 0.545 | <strong style="color: #2ed578;">+0.357</strong> | 0.473 | 1.000 |
| **Pair 08** | 0.625 | 0.259 | 0.270 | 0.545 | <strong style="color: #2ed578;">+0.355</strong> | 0.537 | 1.000 |
| **Pair 09** | 0.632 | 0.295 | 0.282 | 0.545 | <strong style="color: #2ed578;">+0.351</strong> | 0.603 | 1.000 |
| **Pair 10** | 0.629 | 0.255 | 0.279 | 0.545 | <strong style="color: #2ed578;">+0.350</strong> | 0.584 | 1.000 |
| **Pair 11** | **0.741** | 0.240 | 0.402 | 0.545 | <strong style="color: #2ed578;">+0.340</strong> | 0.478 | 1.000 |
| **Pair 12** | 0.639 | 0.246 | 0.313 | 0.545 | <strong style="color: #2ed578;">+0.326</strong> | 0.509 | 1.000 |
| **Pair 13** | 0.608 | 0.227 | 0.296 | 0.545 | <strong style="color: #2ed578;">+0.312</strong> | 0.561 | 1.000 |
| **Pair 14** | 0.620 | 0.240 | 0.331 | 0.545 | <strong style="color: #2ed578;">+0.289</strong> | 0.534 | 1.000 |
| **Pair 15** | 0.636 | 0.227 | 0.356 | 0.545 | <strong style="color: #2ed578;">+0.281</strong> | 0.464 | 1.000 |
| **Pair 16** | 0.687 | 0.207 | 0.412 | 0.545 | <strong style="color: #2ed578;">+0.275</strong> | 0.571 | 1.000 |
| **Pair 17** | 0.629 | 0.255 | 0.356 | 0.545 | <strong style="color: #2ed578;">+0.273</strong> | 0.564 | 1.000 |
| **Pair 18** | 0.625 | 0.259 | 0.400 | 0.545 | <strong style="color: #2ed578;">+0.225</strong> | 0.486 | 1.000 |
| **Pair 19** | 0.612 | 0.208 | 0.419 | 0.545 | <strong style="color: #2ed578;">+0.194</strong> | 0.460 | 1.000 |"""

# 1. Update model.md
with open(MODEL_MD, 'r', encoding='utf-8') as f:
    content_md = f.read()

pattern_md = r'### 4\.1 Evaluation Table.*?(?=\n<div style="text-align: center; margin: 2rem 0;")'
match_md = re.search(pattern_md, content_md, re.DOTALL)
if match_md:
    new_md = content_md[:match_md.start()] + table_md + "\n" + content_md[match_md.end():]
    with open(MODEL_MD, 'w', encoding='utf-8') as f:
        f.write(new_md)
    print("Updated wiki/model.md successfully!")
else:
    print("ERROR: Could not match table in wiki/model.md")

# 2. Update index.html
with open(INDEX_HTML, 'r', encoding='utf-8') as f:
    content_html = f.read()

match_html = re.search(pattern_md, content_html, re.DOTALL)
if match_html:
    new_html = content_html[:match_html.start()] + table_md + "\n" + content_html[match_html.end():]
    # Also add cache busting to fetch('model.md') in index.html to guarantee no stale cache
    new_html = new_html.replace("fetch('model.md')", "fetch('model.md?v=' + Date.now(), { cache: 'no-store' })")
    with open(INDEX_HTML, 'w', encoding='utf-8') as f:
        f.write(new_html)
    print("Updated wiki/index.html successfully!")
else:
    print("ERROR: Could not match table in wiki/index.html")

# 3. Build React TSX Table for model.tsx
pairs_data = [
    ("01", "0.745", "0.295", "0.251", "0.545", "+0.494", "0.481", "1.000"),
    ("02", "0.781", "0.237", "0.299", "0.545", "+0.483", "0.547", "1.000"),
    ("03", "0.701", "0.295", "0.259", "0.545", "+0.442", "0.582", "1.000"),
    ("04", "0.602", "0.246", "0.221", "0.545", "+0.381", "0.497", "1.000"),
    ("05", "0.626", "0.259", "0.247", "0.545", "+0.379", "0.507", "1.000"),
    ("06", "0.741", "0.295", "0.384", "0.545", "+0.357", "0.540", "1.000"),
    ("07", "0.600", "0.259", "0.243", "0.545", "+0.357", "0.473", "1.000"),
    ("08", "0.625", "0.259", "0.270", "0.545", "+0.355", "0.537", "1.000"),
    ("09", "0.632", "0.295", "0.282", "0.545", "+0.351", "0.603", "1.000"),
    ("10", "0.629", "0.255", "0.279", "0.545", "+0.350", "0.584", "1.000"),
    ("11", "0.741", "0.240", "0.402", "0.545", "+0.340", "0.478", "1.000"),
    ("12", "0.639", "0.246", "0.313", "0.545", "+0.326", "0.509", "1.000"),
    ("13", "0.608", "0.227", "0.296", "0.545", "+0.312", "0.561", "1.000"),
    ("14", "0.620", "0.240", "0.331", "0.545", "+0.289", "0.534", "1.000"),
    ("15", "0.636", "0.227", "0.356", "0.545", "+0.281", "0.464", "1.000"),
    ("16", "0.687", "0.207", "0.412", "0.545", "+0.275", "0.571", "1.000"),
    ("17", "0.629", "0.255", "0.356", "0.545", "+0.273", "0.564", "1.000"),
    ("18", "0.625", "0.259", "0.400", "0.545", "+0.225", "0.486", "1.000"),
    ("19", "0.612", "0.208", "0.419", "0.545", "+0.194", "0.460", "1.000"),
]

tsx_rows = []
for p, s1, s2, s3, s4, d_iptm, dq, dq_wt in pairs_data:
    highlight_rescue = 'fw-bold text-success' if float(s1) >= 0.70 else 'text-light'
    tsx_rows.append(f"""              <tr>
                <td className="text-center fw-bold">Pair {p}</td>
                <td className="text-center {highlight_rescue}">{s1}</td>
                <td className="text-center text-secondary">{s2}</td>
                <td className="text-center text-secondary">{s3}</td>
                <td className="text-center text-muted">{s4}</td>
                <td className="text-center fw-bold text-success fs-6">{d_iptm}</td>
                <td className="text-center fw-bold">{dq}</td>
                <td className="text-center text-info">{dq_wt}</td>
              </tr>""")

tsx_table = f"""        <h3 className="h4 fw-bold mt-5 mb-3 text-info" id="champions-table">
          4.1 Evaluation Table of the 19 Orthogonal Champions: 4-Way iPTM Matrix &amp; CAPRI DockQ
        </h3>
        <p className="small text-secondary mb-3">
          Core biophysical evaluation of the 19 strict orthogonal champions. Each column strictly corresponds to our state indices from Section 2.1 ($s_1$ to $s_4$), the resulting selectivity margin $\\Delta\\text{{iPTM}}$, and the CAPRI binding pose fidelity $\\text{{DockQ}}(A' \\cdot B')$ compared to native reference $\\text{{DockQ}}(A \\cdot B) = 1.000$:
        </p>

        <div className="table-responsive my-3">
          <table className="table table-dark table-hover table-bordered small align-middle">
            <thead className="table-success text-dark text-nowrap">
              <tr>
                <th className="text-center">Pair</th>
                <th className="text-center bg-primary text-light">iPTM(s₁)<br/><span className="fw-normal small">A' &middot; B'</span></th>
                <th className="text-center bg-danger text-light">iPTM(s₂)<br/><span className="fw-normal small">A' &middot; B_WT</span></th>
                <th className="text-center bg-danger text-light">iPTM(s₃)<br/><span className="fw-normal small">A_WT &middot; B'</span></th>
                <th className="text-center bg-secondary text-light">iPTM(s₄)<br/><span className="fw-normal small">A_WT &middot; B_WT</span></th>
                <th className="text-center bg-success text-light">&Delta;iPTM</th>
                <th className="text-center">DockQ(A'&middot;B')</th>
                <th className="text-center text-info">DockQ(AB)</th>
              </tr>
            </thead>
            <tbody>
{chr(10).join(tsx_rows)}
            </tbody>
          </table>
        </div>"""

with open(MODEL_TSX, 'r', encoding='utf-8') as f:
    content_tsx = f.read()

pattern_tsx = r'<h3 className="h4 fw-bold mt-5 mb-3 text-info" id="champions-table">.*?</table>\s*</div>'
match_tsx = re.search(pattern_tsx, content_tsx, re.DOTALL)
if match_tsx:
    new_tsx = content_tsx[:match_tsx.start()] + tsx_table + content_tsx[match_tsx.end():]
    with open(MODEL_TSX, 'w', encoding='utf-8') as f:
        f.write(new_tsx)
    print("Updated wiki/src_gitlab/model.tsx successfully!")
else:
    print("ERROR: Could not match table in wiki/src_gitlab/model.tsx")

print("All files updated successfully!")
