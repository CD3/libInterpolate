import argparse
import pathlib
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="A script to help with publishing releases.")
    parser.add_argument('version', help='The version number for the release.')
    args = parser.parse_args()

    if not re.fullmatch(r"\d+\.\d+\.\d+(-([a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*))?", args.version):
        print(f"Error: '{args.version}' is not a valid semantic version string (e.g., 1.2.3 or 1.2.3-rc.1).", file=sys.stderr)
        sys.exit(1)

    if subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True).stdout:
        print("Error: Git working tree is not clean. Please commit or stash your changes.", file=sys.stderr)
        sys.exit(1)

    print(f"Starting release process for version {args.version}")

    version_cmake_path = pathlib.Path('version.cmake')
    print(f"Writing version to '{version_cmake_path}'")
    version_cmake_path.write_text(f'set( libInterpolate_FULL_VERSION "{args.version}")')

    print("\nPlease review the changes below:")
    subprocess.run(['git', 'status'])

    response = input("\nCommit and tag? (yes/no): ").lower().strip()
    if response not in ['yes', 'y']:
        print("Aborting. You can revert the changes with 'git restore {version_cmake_path}'.")
        sys.exit(0)

    print("Committing version bump...")
    try:
        subprocess.run(['git', 'add', 'version.cmake'], check=True)
        subprocess.run(['git', 'commit', '-m', 'chore: version bump'], check=True)
        print("Committed.")

        tag_name = f"v{args.version}"
        print(f"Tagging commit with {tag_name}...")
        subprocess.run(['git', 'tag', tag_name], check=True)
        print("Tagged.")
    except subprocess.CalledProcessError as e:
        print(f"Error during git operation: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Removing version from '{version_cmake_path}'")
    version_cmake_path.write_text(
f'''# this file is used to embed a build version number into th elibrary source.
# when we make a release to tag, this file will set the libInterpolate_FULL_VERSION
# variable here.
    ''')

    try:
        subprocess.run(['git', 'add', 'version.cmake'], check=True)
        subprocess.run(['git', 'commit', '-m', 'chore: removing release version'], check=True)
        print("Committed.")

    except subprocess.CalledProcessError as e:
        print(f"Error during git operation: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
