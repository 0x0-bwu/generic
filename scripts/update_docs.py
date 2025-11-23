#!/usr/bin/env python3
"""
Simple doc updater:
- Adds or updates Doxygen-style comments for C/C++/JS-like files (basic heuristics).
- Adds or updates Python docstrings for functions.
- Regenerates an API list in README.md between markers  and .

This is intentionally conservative and template-driven: it inserts structured placeholders
which you (or a more advanced tool) can refine later.
"""

import re
import os
import sys
import subprocess
from pathlib import Path

SRC_EXTS_C_LIKE = {'.c', '.cpp', '.cc', '.cxx', '.h', '.hpp', '.hh', '.hpp', '.js', '.ts'}
SRC_EXTS_PY = {'.py'}
REPO_ROOT = Path('.').resolve()
README = REPO_ROOT / 'README.md'
DOCS_START = '<!-- AUTO_DOCS_START -->'
DOCS_END = '<!-- AUTO_DOCS_END -->'

def find_source_files():
    """
    Brief description of find_source_files.
    :returns:
    """
    files = []
    for root, dirs, filenames in os.walk(REPO_ROOT):
        if '.git' in root.split(os.sep):
            continue
        for fn in filenames:
            ext = Path(fn).suffix
            if ext in SRC_EXTS_C_LIKE or ext in SRC_EXTS_PY:
                files.append(Path(root) / fn)
    return files

def ensure_py_docstring(text):
    """
    Brief description of ensure_py_docstring.
    :param text:
    :returns:
    """
    lines = text.splitlines()
    out = []
    i = 0
    sigs = []
    while i < len(lines):
        line = lines[i]
        m = re.match(r'^(\s*)def\s+([A-Za-z_]\w*)\s*\((.*?)\)\s*:', line)
        out.append(line)
        if m:
            indent, name, params = m.groups()
            sigs.append(f'{name}({params.strip()})')
            j = i + 1
            while j < len(lines) and re.match(r'^\s*($|#)', lines[j]):
                out.append(lines[j])
                j += 1
            has_doc = False
            if j < len(lines):
                if re.match(r'^\s*(?:[ruRU]{0,2}["\']{3})', lines[j]):
                    has_doc = True
            if not has_doc:
                doc = [
                    f'{indent}    """',
                    f'{indent}    Brief description of {name}.',
                ]
                params_list = [p.strip() for p in params.split(',')] if params.strip() else []
                for p in params_list:
                    if p:
                        pname = p.split('=')[0].strip()
                        if pname:
                            doc.append(f'{indent}    :param {pname}:')
                doc.append(f'{indent}    :returns:')
                doc.append(f'{indent}    """')
                for d in doc:
                    out.append(d)
            i = j
            continue
        i += 1
    return '\n'.join(out) + '\n', sigs

def insert_doxygen_for_c_like(text):
    """
    Brief description of insert_doxygen_for_c_like.
    :param text:
    :returns:
    """
    lines = text.splitlines()
    out = []
    sigs = []
    i = 0
    func_pattern = re.compile(r'^\s*([A-Za-z_][\w:><\*\s,&]*?)\s+([A-Za-z_]\w*)\s*\((.*?)\)\s*\{')
    while i < len(lines):
        line = lines[i]
        m = func_pattern.match(line)
        if m:
            ret, name, params = m.groups()
            sigs.append(f'{name}({params.strip()})')
            k = len(out) - 1
            has_doc = False
            if k >= 0:
                lookback = 3
                snippet = '\n'.join(out[max(0, k - lookback + 1):k + 1])
                if '/**' in snippet or '///' in snippet or '/*!' in snippet:
                    has_doc = True
            if not has_doc:
                params_list = [p.strip() for p in params.split(',')] if params.strip() else []
                doc = []
                doc.append('/**')
                doc.append(f' * @brief Brief description of {name}.')
                for p in params_list:
                    if p:
                        pname = p.split()[-1] if p.split() else p
                        pname = pname.replace('*', '').replace('&', '').split('=')[0].strip()
                        if pname:
                            doc.append(f' * @param {pname}')
                doc.append(f' * @return {ret.strip()}')
                doc.append(' */')
                out.extend(doc)
        out.append(line)
        i += 1
    return '\n'.join(out) + '\n', sigs

def update_readme(api_entries):
    """
    Brief description of update_readme.
    :param api_entries:
    :returns:
    """
    if not README.exists():
        print("README.md not found; skipping README update.")
        return False
    text = README.read_text(encoding='utf-8')
    if DOCS_START not in text or DOCS_END not in text:
        print("README markers not found; insert the markers first:")
        print(DOCS_START)
        print(DOCS_END)
        return False
    start = text.index(DOCS_START) + len(DOCS_START)
    end = text.index(DOCS_END)
    header = '\n\nGenerated API summary\n\n'
    body_lines = []
    for s in sorted(set(api_entries)):
        body_lines.append(f'- `{s}` — short description.')
    new_section = header + '\n'.join(body_lines) + '\n\n'
    new_text = text[:start] + '\n' + new_section + text[end:]
    if new_text != text:
        README.write_text(new_text, encoding='utf-8')
        return True
    return False


def git_commit_and_push():
    """
    Brief description of git_commit_and_push.
    :returns:
    """
    try:
        res = subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True, check=True)
        if not res.stdout.strip():
            print("No changes to commit.")
            return False
        subprocess.run(['git', 'config', 'user.name', os.environ.get('GITHUB_ACTOR', 'doc-agent')], check=True)
        subprocess.run(['git', 'config', 'user.email', f"{os.environ.get('GITHUB_ACTOR','doc-agent')}@users.noreply.github.com"], check=True)
        subprocess.run(['git', 'add', '-A'], check=True)
        subprocess.run(['git', 'commit', '-m', 'docs: update doxygen-style comments and README'], check=True)
        subprocess.run(['git', 'push'], check=True)
        print("Changes committed and pushed.")
        return True
    except subprocess.CalledProcessError as e:
        print("Git operation failed:", e)
        return False


def main():
    """
    Brief description of main.
    :returns:
    """
    files = find_source_files()
    print(f"Found {len(files)} source files to scan.")
    all_sigs = []
    changed = False
    for f in files:
        try:
            text = f.read_text(encoding='utf-8')
        except Exception as e:
            print(f"Skipping {f}: could not read ({e})")
            continue
        new_text = text
        sigs = []
        if f.suffix in SRC_EXTS_PY:
            new_text, sigs = ensure_py_docstring(text)
        elif f.suffix in SRC_EXTS_C_LIKE:
            new_text, sigs = insert_doxygen_for_c_like(text)
        if new_text != text:
            f.write_text(new_text, encoding='utf-8')
            print(f"Updated {f}")
            changed = True
        all_sigs.extend(sigs)
    try:
        readme_changed = update_readme(all_sigs)
        if readme_changed:
            print("Updated README.md API section.")
            changed = True
    except Exception as e:
        print("Failed to update README.md:", e)

    if changed:
        ok = git_commit_and_push()
        if not ok:
            print("Failed to commit/push changes or no permission. If running locally, run 'git status' and push manually.")
    else:
        print("No updates were necessary.")

if __name__ == '__main__':
    main()
