"""Entry point for the frozen app.

PyInstaller runs its entry script as __main__ with no package, so the app
cannot start from nanoshm/desktop.py directly: its relative imports would
fail. This imports the package properly and hands over.
"""
import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from nanoshm.desktop import main
    sys.exit(main())
