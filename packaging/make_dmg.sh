#!/bin/sh
# Wraps dist/SHMtrace.app in a drag-to-Applications disk image, the way
# SciCanvas ships (electron-builder's dmg target). Run after PyInstaller:
#   sh packaging/make_dmg.sh
# Output: dist/SHMtrace-<version>-<arch>.dmg
set -e
cd "$(dirname "$0")/.."
APP=dist/SHMtrace.app
[ -d "$APP" ] || { echo "build $APP first: pyinstaller packaging/shmtrace.spec"; exit 1; }
VER=$(/usr/libexec/PlistBuddy -c "Print CFBundleShortVersionString" "$APP/Contents/Info.plist")
ARCH=$(uname -m)
STAGE=$(mktemp -d)
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
OUT="dist/SHMtrace-$VER-$ARCH.dmg"
rm -f "$OUT"
hdiutil create -volname "SHMtrace $VER" -srcfolder "$STAGE" -fs HFS+ -format UDZO -ov "$OUT"
rm -rf "$STAGE"
echo "$OUT"
