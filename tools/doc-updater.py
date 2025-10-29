#!/usr/bin/env python3
"""
Documentation updater script for the generic library.
Scans source files and updates README.md with API documentation.
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict


def extract_file_description(filepath):
    """Extract the @brief description from a file's Doxygen comment."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # Look for @file and @brief in the file header
        file_match = re.search(r'@file\s+(\S+)', content)
        brief_match = re.search(r'@brief\s+(.+?)(?:\n\s*\*\s*@|\n\s*\*/)', content, re.DOTALL)
        
        if file_match and brief_match:
            filename = file_match.group(1)
            brief = brief_match.group(1).strip()
            # Clean up the brief description
            brief = re.sub(r'\s*\n\s*\*\s*', ' ', brief)
            brief = brief.strip()
            return filename, brief
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
    
    return None, None


def find_header_files(root_dir):
    """Find all .hpp and .h files in the repository."""
    header_files = []
    root_path = Path(root_dir)
    
    for pattern in ['**/*.hpp', '**/*.h']:
        for filepath in root_path.glob(pattern):
            # Skip test files and hidden directories
            if '.git' in filepath.parts or 'test' in filepath.parts:
                continue
            header_files.append(filepath)
    
    return sorted(header_files)


def organize_by_directory(files_info):
    """Organize files by their directory."""
    organized = defaultdict(list)
    
    for filepath, filename, brief in files_info:
        # Get the directory relative to root
        dir_name = filepath.parent.name
        if dir_name:
            organized[dir_name].append((filename, brief))
    
    return organized


def generate_api_docs(root_dir):
    """Generate API documentation from source files."""
    header_files = find_header_files(root_dir)
    
    files_info = []
    for filepath in header_files:
        filename, brief = extract_file_description(filepath)
        if filename and brief:
            files_info.append((filepath, filename, brief))
    
    # Organize by directory
    organized = organize_by_directory(files_info)
    
    # Generate markdown
    lines = []
    lines.append("### Header Files Overview\n")
    
    # Sort directories
    for dir_name in sorted(organized.keys()):
        lines.append(f"\n#### {dir_name}/\n")
        for filename, brief in sorted(organized[dir_name]):
            lines.append(f"- **{filename}**: {brief}\n")
    
    return ''.join(lines)


def update_readme(readme_path, api_docs):
    """Update README.md with the generated API documentation."""
    try:
        with open(readme_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Find the AUTO_DOCS markers
        start_marker = '<!-- AUTO_DOCS_START -->'
        end_marker = '<!-- AUTO_DOCS_END -->'
        
        if start_marker not in content or end_marker not in content:
            print(f"Error: AUTO_DOCS markers not found in {readme_path}", file=sys.stderr)
            return False
        
        # Replace content between markers
        pattern = f'{re.escape(start_marker)}.*?{re.escape(end_marker)}'
        new_content = f'{start_marker}\n{api_docs}{end_marker}'
        updated = re.sub(pattern, new_content, content, flags=re.DOTALL)
        
        # Write back
        with open(readme_path, 'w', encoding='utf-8') as f:
            f.write(updated)
        
        print(f"Successfully updated {readme_path}")
        return True
        
    except Exception as e:
        print(f"Error updating README: {e}", file=sys.stderr)
        return False


def main():
    """Main function."""
    # Get the repository root directory
    script_dir = Path(__file__).parent.absolute()
    repo_root = script_dir.parent
    
    print(f"Scanning repository at: {repo_root}")
    
    # Generate API documentation
    api_docs = generate_api_docs(repo_root)
    
    # Update README.md
    readme_path = repo_root / 'README.md'
    if update_readme(readme_path, api_docs):
        print("Documentation update completed successfully!")
        return 0
    else:
        print("Documentation update failed!", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
