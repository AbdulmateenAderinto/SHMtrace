"""Entry point for the frozen app.

PyInstaller runs its entry script as __main__ with no package, so the app
cannot start from shmtrace/desktop.py directly: its relative imports would
fail. This imports the package properly and hands over.
"""
import multiprocessing
import os
import sys
from pathlib import Path

if __name__ == "__main__":
    multiprocessing.freeze_support()
    # PyInstaller's matplotlib hook points MPLCONFIGDIR at a fresh temporary
    # folder on every launch, so the font cache is rebuilt each time and startup
    # takes several seconds longer. A user cache folder keeps it between launches.
    cache = Path.home() / ("Library/Caches" if sys.platform == "darwin" else ".cache") \
        / "shmtrace" / "matplotlib"
    try:
        cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(cache)
    except OSError:
        pass
    from shmtrace.desktop import main
    sys.exit(main())
