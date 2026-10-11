# Running nanoshm on your Mac

Two ways. The first takes a few minutes and gives you the app; the second
builds the double-clickable version.

## 1. Run it now

```bash
cd ~/Nanoshm/nanoshm_v0.2
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[app,bam,sanger]"
nanoshm-app
```

A window opens. Everything runs on your own machine and nothing is uploaded.
If the window does not appear, run `nanoshm-app --no-window` and open the
address it prints in a browser.

Python 3.10 or newer. If `python3` is missing, install it from python.org or
with Homebrew (`brew install python`).

## 2. Build the double-clickable app

```bash
source .venv/bin/activate
pip install pyinstaller
pyinstaller packaging/nanoshm.spec --noconfirm
```

`dist/nanoshm.app` appears. Drag it to Applications. The first launch needs
right-click then Open, because the app is not signed by Apple. The same spec
builds `dist/nanoshm/nanoshm.exe` on Windows. `packaging/BUILD.md` has the
details, including signing if you want to hand it to other people.

I built and ran this on Linux to check the packaging works end to end: it
starts, runs a full analysis and exports. The Mac build should behave the same,
but it is the one step I could not test from here.

## Where to start

`docs/QUICKSTART.md` is one page. `docs/USER_GUIDE.md` is the full account,
including what every number means and what to check before believing it.
