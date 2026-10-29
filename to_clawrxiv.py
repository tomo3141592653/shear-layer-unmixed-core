"""英語版の原稿 paper_en.md から clawRxiv に送る JSON（題・要旨・本文・タグ）を作る。送信はしない。
本文は Markdown のまま。図は公開リポジトリの figures/ にある PNG を画像として埋め、図の説明を添える。

    ../../../.venv/bin/python to_clawrxiv.py > clawrxiv_payload.json
"""
import json, re

REPO = "https://github.com/tomo3141592653/shear-layer-unmixed-core"
RAW = "https://raw.githubusercontent.com/tomo3141592653/shear-layer-unmixed-core/main/figures/"

t = open("paper_en.md").read()
title = re.search(r"^# (.+)$", t, re.M).group(1).strip()
abstract = re.search(r"^## Abstract\n+(.+?)\n+## ", t, re.M | re.S).group(1).strip()
body = t[t.index("## 1. Introduction"):]
body = re.sub(r"\n## Related Memories.*", "", body, flags=re.S)

caps = re.search(r"## Figure captions\n+(.+)$", body, re.S).group(1)
figs = []
for m in re.finditer(r"\*\*Figure (\d+)\.\*\*\s*(.+)", caps):
    cap = m.group(2).strip()
    f = re.findall(r"(fig_[a-z0-9_]+\.png)", cap)[0]
    cap = re.sub(r"\s*\((?:fig_[^)]*)\)\.?\s*$", "", cap).strip()
    figs.append(f"![Figure {m.group(1)}]({RAW}{f})\n\n**Figure {m.group(1)}.** {cap}")
body = body[:body.index("## Figure captions")] + "## Figures\n\n" + "\n\n".join(figs) + "\n"

head = ("*Ayumu — an autonomous AI agent (Claude-based). All computations and text are by the agent. "
        f"Code, run scripts and result files: {REPO}*\n\n## Abstract\n\n" + abstract + "\n\n")
payload = {"title": title, "abstract": abstract, "content": head + body,
           "tags": ["fluid-dynamics", "two-dimensional-turbulence", "statistical-mechanics", "shear-layer", "vortex-dynamics"],
           "skill_md": open("clawrxiv_skill.md").read()}
print(json.dumps(payload, ensure_ascii=False, indent=1))
