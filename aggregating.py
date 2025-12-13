#!/usr/bin/env python3
"""
File Collector Script
Collects all files with relative path comments into one big text file.
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Set, List, Dict
from datetime import datetime

# Directories to ignore (Added 'test' and 'public')
DEFAULT_IGNORE_DIRS = {
    '.git', '.svn', '.hg', '__pycache__', 'node_modules',
    '.idea', '.vscode', '.vs', 'build', 'dist', 'target',
    'venv', 'env', '.env', 'virtualenv', 'coverage', '.mypy_cache',
    'bin', 'obj', 'test', 'tests', 'public', '.next'
}

# Specific files to ignore (Added configs, locks, and script files)
DEFAULT_IGNORE_FILES = {
    # Lock files
    'package-lock.json', 'yarn.lock', 'pnpm-lock.yaml', 'bun.lockb',
    'cargo.lock', 'Gemfile.lock', 'composer.lock', 'mix.lock',
    'poetry.lock', 'Pipfile.lock', 'go.sum', 'migration_lock.toml',

    # System/IDE
    '.DS_Store', 'Thumbs.db', 'LICENSE', 'LICENSE.txt',

    # Configs (Boilerplate noise)
    'tsconfig.json', 'tsconfig.app.json', 'tsconfig.node.json', 'tsconfig.build.json',
    'eslint.config.js', '.eslintrc.js', '.eslintrc.json',
    'postcss.config.js', 'tailwind.config.js', 'vite.config.ts', 'vite.config.js',
    'nest-cli.json', '.prettierrc', 'jest-e2e.json', 'jest.config.js',

    # Scripts themselves (to prevent collecting them)
    'aggregating.py', 'indexing.py', 'cleanup.py',

    # Web root
    'index.html', 'favicon.ico'
}

# File extensions to process
DEFAULT_EXTENSIONS = {
    '.py', '.js', '.ts', '.jsx', '.tsx', '.java', '.c', '.cpp', '.h', '.hpp',
    '.cs', '.go', '.rs', '.php', '.rb', '.swift', '.kt', '.scala', '.m', '.mm',
    '.css', '.scss', '.sass', '.less', '.xml', '.json', '.yaml', '.yml',
    '.sql', '.graphql', '.gql', '.dockerfile', 'dockerfile',
    '.toml', '.ini', '.cfg', '.conf', '.env'
}

# Patterns to identify relative path comments
COMMENT_PATTERNS = {
    '# {}', '// {}', '/* {} */', '<!-- {} -->', '-- {}', '.. {}'
}

def has_relative_path_comment(file_path: Path, relative_path: str) -> bool:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for i in range(10):
                line = f.readline()
                if not line: break
                for pattern in COMMENT_PATTERNS:
                    if pattern.format(relative_path) in line:
                        return True
    except:
        return False
    return False

def read_file_content(file_path: Path, max_lines: int = None) -> str:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            if max_lines:
                lines = [f.readline() for _ in range(max_lines)]
                return ''.join(line for line in lines if line)
            return f.read()
    except UnicodeDecodeError:
        return f"[Binary file or encoding issue: {file_path}]\n"
    except Exception as e:
        return f"[Error reading file: {e}]\n"

def collect_files_with_comments(
    root_dir: Path, base_dir: Path, ignore_dirs: Set[str],
    ignore_files: Set[str], extensions: Set[str], verbose: bool = False
) -> List[Dict]:
    collected_files = []

    for item in root_dir.iterdir():
        if item.name in ignore_files:
            continue

        if item.is_dir():
            if item.name in ignore_dirs:
                continue
            collected_files.extend(collect_files_with_comments(
                item, base_dir, ignore_dirs, ignore_files, extensions, verbose
            ))
            continue

        # Ignore test files explicitly
        if item.name.endswith('.spec.ts') or item.name.endswith('.test.ts'):
            continue

        file_ext = item.suffix.lower()
        if file_ext == '' and item.name.lower() == 'dockerfile': file_ext = 'dockerfile'

        if file_ext not in extensions:
            continue

        try:
            relative_path = str(item.relative_to(base_dir)).replace('\\', '/')
        except ValueError:
            relative_path = str(item)

        if has_relative_path_comment(item, relative_path):
            collected_files.append({
                'path': relative_path,
                'full_path': item,
                'size': item.stat().st_size
            })

    return sorted(collected_files, key=lambda x: x['path'])

def create_collection_file(output_file: Path, collected_files: List[Dict], verbose: bool = False) -> Dict:
    stats = {'total_files': len(collected_files), 'total_size': 0}

    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"Project Collection\nTotal files: {len(collected_files)}\n{'='*80}\n\n")

            # File Index
            f.write("INDEX:\n")
            for file_info in collected_files:
                f.write(f"{file_info['path']}\n")
            f.write(f"{'-'*80}\n\n")

            for idx, file_info in enumerate(collected_files):
                stats['total_size'] += file_info['size']
                if verbose: print(f"Adding: {file_info['path']}")

                f.write(f"\n{'='*80}\nFILE: {file_info['path']}\n{'='*80}\n\n")
                f.write(read_file_content(file_info['full_path']))
                f.write(f"\n")

    except Exception as e:
        print(f"Error: {e}")
        if output_file.exists(): output_file.unlink()
        raise
    return stats

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', nargs='?', default='.')
    parser.add_argument('--output', '-o', default='project_code.txt')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    base_dir = Path(args.directory).resolve()
    output_file = Path(args.output).resolve()

    # Dynamic ignores
    ignore_files = DEFAULT_IGNORE_FILES.copy()
    ignore_files.add(output_file.name) # Ignore the output file itself

    # Also ignore other common output names
    ignore_files.update({'testing.txt', 'output.txt', 'code_dump.txt'})

    print(f"Scanning: {base_dir}")
    print(f"Ignoring configs and system files...")

    collected = collect_files_with_comments(
        base_dir, base_dir, DEFAULT_IGNORE_DIRS, ignore_files, DEFAULT_EXTENSIONS, args.verbose
    )

    if not collected:
        print("No files found. Did you run indexing.py first?")
        sys.exit(0)

    stats = create_collection_file(output_file, collected, args.verbose)
    print(f"\nComplete! Collected {stats['total_files']} files into {output_file.name}")

if __name__ == '__main__':
    main()
