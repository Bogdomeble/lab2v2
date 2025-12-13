#!/usr/bin/env python3
"""
File Comment Cleanup Script
Removes relative path comments added by the indexing script.
Restores files to their original state.
"""

import argparse
import os
from pathlib import Path
from typing import Set

# Strict ignores (Same as before)
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
    "eslint.config.js",
    ".eslintrc.js",
    "postcss.config.js",
    "tailwind.config.js",
    "vite.config.ts",
    "nest-cli.json",
    ".prettierrc",
    "aggregating.py",
    "indexing.py",
    "cleanup.py",
    "output.txt",
    "testing.txt",
}

# Added .json explicitly to ensure it matches the # {} syntax
COMMENT_SYNTAX = {
    "default": "# {}",
    ".js": "// {}",
    ".ts": "// {}",
    ".jsx": "// {}",
    ".tsx": "// {}",
    ".java": "// {}",
    ".c": "// {}",
    ".cpp": "// {}",
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
    ".fxml": "<!-- {} -->",
    ".sql": "-- {}",
    ".json": "# {}",  # Explicitly defined to match indexing.py output
}


def get_expected_comment(ext: str, path: str) -> str:
    """Get the comment string that we expect to find at the top of the file."""
    # Fallback to default if extension not found
    syntax = COMMENT_SYNTAX.get(ext.lower(), COMMENT_SYNTAX["default"])
    return syntax.format(path)


def clean_file(file_path: Path, rel_path: str, dry_run: bool) -> bool:
    """
    Remove the first line if it matches the expected relative path comment.
    """
    try:
        # Read file into memory
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except:
        return False

    if not lines:
        return False

    ext = file_path.suffix.lower()
    if ext == "" and file_path.name.lower() == "dockerfile":
        ext = "dockerfile"

    # Generate what the comment should look like
    expected = get_expected_comment(ext, rel_path)

    # Check if the first line is the comment
    if lines[0].strip() == expected:
        if not dry_run:
            try:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.writelines(lines[1:])  # Write everything except the first line
                return True
            except Exception as e:
                print(f"Error writing {file_path}: {e}")
                return False
        return True
    return False


def process_dir(root: Path, base: Path, dry_run: bool, verbose: bool) -> int:
    cleaned_count = 0
    for item in root.iterdir():
        # Skip specific ignore files
        if item.name in DEFAULT_IGNORE_FILES:
            continue

        # Recursively process directories
        if item.is_dir():
            if item.name in DEFAULT_IGNORE_DIRS:
                continue
            cleaned_count += process_dir(item, base, dry_run, verbose)
            continue

        # Process the file (JSON check removed)
        try:
            rel_path = str(item.relative_to(base)).replace("\\", "/")
        except ValueError:
            rel_path = str(item)

        if clean_file(item, rel_path, dry_run):
            cleaned_count += 1
            if verbose or dry_run:
                print(f"{'[DRY] ' if dry_run else ''}Cleaned: {rel_path}")

    return cleaned_count


def main():
    parser = argparse.ArgumentParser(
        description="Remove relative path comments from files"
    )
    parser.add_argument("directory", nargs="?", default=".", help="Project directory")
    parser.add_argument("--dry-run", action="store_true", help="Simulate cleanup")
    parser.add_argument("--verbose", "-v", action="store_true", help="Detailed output")
    args = parser.parse_args()

    base = Path(args.directory).resolve()
    if not base.exists():
        print(f"Error: {base} does not exist.")
        return

    print(f"Cleaning directory: {base}")
    if args.dry_run:
        print("DRY RUN MODE - No files will be modified")

    count = process_dir(base, base, args.dry_run, args.verbose)

    print(f"\nCleanup complete. Total files cleaned: {count}")
    if count > 0 and not args.dry_run:
        print("Your JSON files should now be valid again.")


if __name__ == "__main__":
    main()
