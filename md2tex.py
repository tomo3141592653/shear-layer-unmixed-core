"""英語版の原稿 paper_en.md を LaTeX（paper_en.tex）に起こす。pandoc を使わず、この原稿で使っている書き方だけを扱う。
記号は TeX の数式に置き換え、付属の標準書体だけで組めるようにする（tectonic paper_en.tex で PDF）。

    python md2tex.py && tectonic paper_en.tex
"""
import re

GREEK = {"α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta", "ε": r"\varepsilon", "θ": r"\theta", "λ": r"\lambda",
         "μ": r"\mu", "ν": r"\nu", "π": r"\pi", "ρ": r"\rho", "σ": r"\sigma", "φ": r"\phi", "ψ": r"\psi", "ω": r"\omega",
         "Γ": r"\Gamma", "Δ": r"\Delta", "Σ": r"\Sigma", "Ω": r"\Omega", "Ψ": r"\Psi"}
MATHSYM = {"×": r"\times", "≈": r"\approx", "≤": r"\le", "≥": r"\ge", "±": r"\pm", "−": "-", "∇": r"\nabla", "∂": r"\partial",
           "∈": r"\in", "∝": r"\propto", "→": r"\to", "′": r"\prime", "·": r"\cdot", "°": r"^\circ"}
SUP = dict(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻ˣ", "0123456789-x"))
SUP.update({"²": "2", "³": "3", "¹": "1"})
SUB = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
TEXT = {"–": "--", "—": "---", "…": r"\ldots{}", "é": r"\'e", "ü": r'\"u'}
SPECIAL = {"&": r"\&", "%": r"\%", "#": r"\#", "$": r"\$", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\^{}", "\\": r"\textbackslash{}"}


def esc_plain(s):
    """数式でない平文の特殊文字だけを逃がす（_ は呼ぶ側で数式にしてから来る）。"""
    return "".join(SPECIAL.get(c, c) for c in s)


def conv_text(s):
    """一続きの平文を、記号を数式にしながら LaTeX にする。"""
    s = s.replace("√(νt)", "\0SQRT\0")
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if s.startswith("\0SQRT\0", i):
            out.append(r"$\sqrt{\nu t}$"); i += 6; continue
        nxt = s[i + 1] if i + 1 < len(s) else ""
        if nxt == "̄":                                   # 結合用の上線: x̄
            base = GREEK.get(c, c)
            out.append(f"$\\bar{{{base}}}$"); i += 2; continue
        if c in GREEK or c in MATHSYM:
            m = GREEK.get(c) or MATHSYM[c]
            j = i + 1
            sub = ""
            while j < len(s) and (s[j] in SUB or (s[j] in "0123456789" and c in GREEK and s[i + 1:j + 1].isdigit() and False)):
                sub += SUB[s[j]]; j += 1
            out.append(f"${m}{'_{' + sub + '}' if sub else ''}{'{}' if m.isalpha() or m.startswith(chr(92)) else ''}$"); i = j; continue
        if c in SUP:
            j = i
            sup = ""
            while j < len(s) and s[j] in SUP:
                sup += SUP[s[j]]; j += 1
            out.append(f"$^{{{sup}}}$"); i = j; continue
        if c in SUB:
            j = i
            sub = ""
            while j < len(s) and s[j] in SUB:
                sub += SUB[s[j]]; j += 1
            out.append(f"$_{{{sub}}}$"); i = j; continue
        if c in TEXT:
            out.append(TEXT[c]); i += 1; continue
        if c == "_":
            out.append(r"\_"); i += 1; continue
        out.append(SPECIAL.get(c, c)); i += 1
    t = "".join(out)
    t = re.sub(r"\$\$", "", t)                                # 隣り合った数式をつなぐ
    return t


def sub_body(s):
    """添字の中身: 1 語なら立体、i+1 のような式はそのまま数式に。"""
    s = s.replace("\\{", "").replace("\\}", "")
    return s if re.search(r"[+\-]", s) else r"\mathrm{" + s + "}"


def inline(s):
    """行の中の `code`・**太字**・*斜体*・添字 X_c を処理する。"""
    parts = re.split(r"(`[^`]+`)", s)
    out = []
    for p in parts:
        if p.startswith("`") and p.endswith("`"):
            out.append(r"\texttt{" + esc_plain(p[1:-1]).replace("_", r"\_") + "}")
            continue
        p = re.sub(r"(?<![A-Za-z])([A-Za-zα-ωΓΔΣΩΨ])_(?:\{([^}]+)\}|([A-Za-z0-9]+))", lambda m: "\x01" + m.group(1) + "\x02" + (m.group(2) or m.group(3)) + "\x03", p)
        p = re.sub(r"\*\*(.+?)\*\*", "\x04\\1\x05", p)
        p = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", "\x06\\1\x07", p)
        t = conv_text(p)
        t = re.sub("\x01(.+?)\x02(.+?)\x03", lambda m: f"${m.group(1).strip('$')}_{{{sub_body(m.group(2))}}}$", t)
        t = t.replace("\x04", r"\textbf{").replace("\x05", "}").replace("\x06", r"\emph{").replace("\x07", "}")
        out.append(re.sub(r"\$\$", "", t))
    return "".join(out)


def main():
    lines = open("paper_en.md").read().split("\n")
    body, title, figs = [], "", []
    i, mode, tabcap = 0, None, ""
    while i < len(lines):
        l = lines[i]
        if l.startswith("(English manuscript"):
            i += 1; continue
        if l.startswith("# "):
            title = inline(l[2:]); i += 1; continue
        if l.startswith("## Figure captions"):
            if mode == "ref":
                body.append(r"\end{list}\end{small}")
            mode = "fig"; i += 1; continue
        if mode == "fig":
            m = re.match(r"\*\*Figure (\d+)\.\*\*\s*(.*)", l)
            if m:
                cap = m.group(2)
                files = re.findall(r"(fig_[a-z0-9_]+\.png)", cap)
                cap = re.sub(r"\s*\((?:fig_[^)]*)\)\.?\s*$", "", cap).strip()
                figs.append((files[0] if files else None, inline(cap)))
            i += 1; continue
        if l.startswith("## "):
            h = re.sub(r"^\d+\.\s*", "", l[3:])
            star = "" if re.match(r"^\d", l[3:]) else "*"
            if h == "Abstract":
                body.append(r"\begin{abstract}"); mode = "abs"; i += 1; continue
            if mode == "abs":
                body.append(r"\end{abstract}"); mode = None
            if h == "References":
                body.append(r"\section*{References}" + "\n" + r"\begin{small}\begin{list}{}{\leftmargin=1.5em\itemindent=-1.5em\itemsep=1pt}")
                mode = "ref"; i += 1; continue
            if mode == "ref":
                body.append(r"\end{list}\end{small}"); mode = None
            body.append(f"\\section{star}{{{inline(h)}}}"); i += 1; continue
        if l.startswith("### "):
            h = re.sub(r"^\d+\.\d+\s*", "", l[4:])
            body.append(f"\\subsection{{{inline(h)}}}"); i += 1; continue
        if mode == "ref" and l.startswith("- "):
            body.append(r"\item " + inline(l[2:])); i += 1; continue
        if re.match(r"^Table \d+\.", l):
            tabcap = inline(re.sub(r"^Table \d+\.\s*", "", l)); i += 1; continue
        if l.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                if not re.match(r"^\|[-| ]+\|$", lines[i]):
                    rows.append([c.strip() for c in lines[i].strip("|").split("|")])
                i += 1
            n = len(rows[0])
            cap = (r"\caption{" + tabcap + r"}\label{tab1}") if tabcap else ""
            tabcap = ""
            sup = lambda m: m.group(1) + "×10⁻" + "".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(d)] for d in m.group(2))
            body.append(r"\begin{table}[t]\centering" + cap + r"\scriptsize\setlength{\tabcolsep}{3pt}\begin{tabular}{l" + "c" * (n - 1) + r"}\hline")
            for k, r in enumerate(rows):
                body.append(" & ".join(inline(re.sub(r"(\d)e-(\d+)", sup, c)) for c in r) + r" \\" + (r"\hline" if k == 0 else ""))
            body.append(r"\hline\end{tabular}\end{table}")
            continue
        if re.match(r"^(- |\d+\. )", l):
            env = "itemize" if l.startswith("- ") else "enumerate"
            body.append(f"\\begin{{{env}}}")
            while i < len(lines) and re.match(r"^(- |\d+\. )", lines[i]):
                body.append(r"\item " + inline(re.sub(r"^(- |\d+\. )", "", lines[i]))); i += 1
            body.append(f"\\end{{{env}}}")
            continue
        body.append(inline(l)); i += 1
    for k, (f, cap) in enumerate(figs, 1):
        if f:
            body.append(r"\begin{figure}[p]\centering\includegraphics[width=0.85\textwidth]{" + f + r"}\caption{" + cap + r"}\label{fig" + str(k) + r"}\end{figure}")
    head = r"""\documentclass[11pt]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage{amsmath,amssymb,graphicx}
\usepackage[hidelinks]{hyperref}
\title{""" + title + r"""}
\author{Ayumu}
\date{}
\begin{document}
\maketitle
"""
    open("paper_en.tex", "w").write(head + "\n".join(body) + "\n\\end{document}\n")


main()
