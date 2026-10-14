# Building the SHMtrace desktop app

The app is a Python program with a local web UI in a native window. It is built
on the machine it will run on: macOS builds the `.app`, Windows builds the
`.exe`. There is no cross-building.

## Both platforms

```bash
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[app,bam,sanger]" pyinstaller
pyinstaller packaging/shmtrace.spec --noconfirm
```

Output:
- **macOS:** `dist/SHMtrace.app`. Drag it to Applications. The first launch needs
  right-click then Open, because the app is not signed by Apple. To avoid that for
  other people, sign and notarise it: `codesign --deep --force --options runtime
  --sign "Developer ID Application: <you>" dist/SHMtrace.app`, then
  `xcrun notarytool submit`.
- **Windows:** `dist/shmtrace/shmtrace.exe` with its folder. Zip the folder, or wrap
  it with Inno Setup for an installer.

## macOS installer (like SciCanvas)

The spec gives the app its icon (`packaging/icon/icon.icns`, redrawn with
`python packaging/icon/make_icon.py`) and its name, version and category. After
PyInstaller, wrap it in a drag-to-Applications disk image:

```bash
sh packaging/make_dmg.sh        # writes dist/SHMtrace-0.2.1-arm64.dmg
```

## What it does at run time

It serves on `127.0.0.1` on a free port and opens a window on it. Nothing listens
on the network and no data leaves the machine. `shmtrace-app --no-window` serves
without a window and prints the address, which is how to use it over SSH.

## Testing a build

```bash
dist/shmtrace/shmtrace --no-window        # then open the printed address
```

Check: the Set-up tab lists presets and built-in references, a small project runs
to the end, figures appear in Results, and the ZIP export downloads.

## Size

The bundle is roughly 300–400 MB, mostly SciPy, NumPy and matplotlib. Dropping
matplotlib would save about 60 MB but lose the static figure exports, which the
HTML report and the PNG/SVG/PDF buttons use.

## A warning you can ignore

PyInstaller prints `hidden import "scipy.special._cdflib" not found`. That
module was removed in SciPy 1.14 and the warning comes from PyInstaller's own
bundled SciPy hook, not from this spec. The build and the app are fine.
