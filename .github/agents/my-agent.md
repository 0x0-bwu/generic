---
name: doc
description: add or update doc
---

# Doc Agent

This agent scans source files and:
- Adds or updates Doxygen-style comments for C/C++/JS-like files.
- Adds or updates Python function docstrings.
- Regenerates an API section in README.md between the markers:
  <!-- AUTO_DOCS_START --> and <!-- AUTO_DOCS_END -->

Usage (workflow): a GitHub Actions workflow runs the doc-updater script and will commit & push changes back to the branch when updates are made.


Requirements/Behavior:
- README.md must include the markers <!-- AUTO_DOCS_START --> and <!-- AUTO_DOCS_END --> where the API list will be inserted.
- Workflow commits via GITHUB_TOKEN provided to the Actions runner.
- The script is intentionally conservative and inserts placeholder comments for later refinement.
