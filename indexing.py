#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
from typing import Set

# Strict ignores for boilerplate/config
DEFAULT_IGNORE_DIRS = {
    ".git",
    "__pycache__",
    "node_modules",
    ".vscode",
    ".idea",
    "dist",
    "build",
    "coverage",
    "test",
    "tests",
    "public",
    ".next",
}

DEFAULT_IGNORE_FILES = {
    "package-lock.json",
    "yarn.lock",
    "pnpm-lock.yaml",
    "bun.lockb",
    "cargo.lock",
    "migration_lock.toml",
    ".DS_Store",
    "tsconfig.json",
    "tsconfig.app.json",
    "tsconfig.node.json",
    "tsconfig.build.json",
    "eslint.config.js",
    ".eslintrc.js",
    "postcss.config.js",
    "tailwind.config.js",
    "vite.config.ts",
    "nest-cli.json",
    ".prettierrc",
    "jest-e2e.json",
    "aggregating.py",
    "indexing.py",
    "cleanup.py",
    "output.txt",
    "testing.txt",
}

DEFAULT_EXTENSIONS = {
    ".py",
    ".js",
    ".ts",
    ".jsx",
    ".tsx",
    ".java",
    ".c",
    ".cpp",
    ".h",
    ".cs",
    ".go",
    ".rs",
    ".php",
    ".rb",
    ".css",
    ".scss",
    ".html",
    ".sql",
    ".graphql",
    ".dockerfile",
    "dockerfile",
    ".yaml",
    ".yml",
    ".env",
    ".xml",
    ".fxml",
}

COMMENT_SYNTAX = {
    "default": "# {}",
    ".js": "// {}",
    ".ts": "// {}",
    ".jsx": "// {}",
    ".tsx": "// {}",
    ".java": "// {}",
    ".c": "// {}",
    ".cpp": "// {}",
    ".h": "// {}",
    ".cs": "// {}",
    ".go": "// {}",
    ".rs": "// {}",
    ".php": "// {}",
    ".swift": "// {}",
    ".kt": "// {}",
    ".css": "/* {} */",
    ".scss": "/* {} */",
    ".html": "<!-- {} -->",
    ".xml": "<!-- {} -->",
    ".md": "<!-- {} -->",
    ".sql": "-- {}",
}


def get_comment_syntax(ext: str, path: str) -> str:
    return COMMENT_SYNTAX.get(ext.lower(), COMMENT_SYNTAX["default"]).format(path)


def add_comment(file_path: Path, rel_path: str, dry_run: bool) -> bool:
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
    except:
        return False  # Skip binary/encoding errors

    comment = get_comment_syntax(file_path.suffix, rel_path)

    # Check if comment exists in first 5 lines
    if any(comment in line for line in content.split("\n", 5)[:5]):
        return False

    if not dry_run:
        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(f"{comment}\n{content}")
            return True
        except:
            return False
    return True


def process_dir(root: Path, base: Path, dry_run: bool, verbose: bool) -> dict:
    stats = {"processed": 0, "skipped": 0}

    for item in root.iterdir():
        if item.name in DEFAULT_IGNORE_FILES:
            continue

        if item.is_dir():
            if item.name in DEFAULT_IGNORE_DIRS:
                continue
            sub = process_dir(item, base, dry_run, verbose)
            stats["processed"] += sub["processed"]
            stats["skipped"] += sub["skipped"]
            continue

        # Ignore tests and configs manually
        if item.name.endswith(".spec.ts") or item.name.endswith(".test.ts"):
            continue

        ext = item.suffix.lower()
        if ext == "" and item.name.lower() == "dockerfile":
            ext = "dockerfile"

        if ext not in DEFAULT_EXTENSIONS:
            stats["skipped"] += 1
            continue

        rel_path = str(item.relative_to(base)).replace("\\", "/")
        if verbose:
            print(f"Checking: {rel_path}")

        if add_comment(item, rel_path, dry_run):
            stats["processed"] += 1
            print(f"{'[DRY]' if dry_run else '✓'} Added: {rel_path}")

    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", nargs="?", default=".")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    base = Path(args.directory).resolve()
    print(f"Indexing: {base}")

    stats = process_dir(base, base, args.dry_run, args.verbose)
    print(f"\nDone. Processed: {stats['processed']}, Skipped: {stats['skipped']}")


if __name__ == "__main__":
    main()
