#!/bin/sh
# Wraps dist/nanoshm.app in a drag-to-Applications disk image, the way
# SciCanvas ships (electron-builder's dmg target). Run after PyInstaller:
#   sh packaging/make_dmg.sh
# Output: dist/nanoshm-<version>-<arch>.dmg
set -e
cd "$(dirname "$0")/.."
APP=dist/nanoshm.app
[ -d "$APP" ] || { echo "build $APP first: pyinstaller packaging/nanoshm.spec"; exit 1; }
VER=$(/usr/libexec/PlistBuddy -c "Print CFBundleShortVersionString" "$APP/Contents/Info.plist")
ARCH=$(uname -m)
STAGE=$(mktemp -d)
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
OUT="dist/nanoshm-$VER-$ARCH.dmg"
rm -f "$OUT"
hdiutil create -volname "nanoshm $VER" -srcfolder "$STAGE" -fs HFS+ -format UDZO -ov "$OUT"
rm -rf "$STAGE"
echo "$OUT"
