# PyInstaller spec for the shmtrace desktop app.
#   pip install -e ".[app,bam,sanger]" pyinstaller
#   pyinstaller packaging/shmtrace.spec --noconfirm
# Produces dist/shmtrace (a folder) and, on macOS, dist/SHMtrace.app.
import sys
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = [
    ("../src/shmtrace/webapp/static", "shmtrace/webapp/static"),
    ("../src/shmtrace/presets", "shmtrace/presets"),
    ("../src/shmtrace/references", "shmtrace/references"),
]
datas += collect_data_files("mappy", include_py_files=False)

hiddenimports = (
    collect_submodules("uvicorn")
    + collect_submodules("shmtrace")
    + ["anyio", "h11", "click", "multipart", "openpyxl", "scipy.special.cython_special"]
    # matplotlib picks its PDF and SVG writers by name at save time, so
    # PyInstaller cannot see them and the figure exports fail without these
    + ["matplotlib.backends.backend_agg", "matplotlib.backends.backend_pdf",
       "matplotlib.backends.backend_svg"]
)

a = Analysis(
    ["launch.py"],
    pathex=["../src"],
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["tkinter", "pytest", "IPython", "notebook", "PyQt5", "PySide6"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="shmtrace",
          console=False, disable_windowed_traceback=False, argv_emulation=sys.platform == "darwin")
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="shmtrace")

if sys.platform == "darwin":
    app = BUNDLE(coll, name="SHMtrace.app", icon="icon/icon.icns",
                 bundle_identifier="org.martinlab.shmtrace",
                 info_plist={"CFBundleName": "SHMtrace",
                             "CFBundleDisplayName": "SHMtrace",
                             "CFBundleShortVersionString": "0.2.1",
                             "CFBundleVersion": "0.2.1",
                             "LSApplicationCategoryType": "public.app-category.education",
                             "NSHighResolutionCapable": True,
                             "LSMinimumSystemVersion": "11.0",
                             "NSHumanReadableCopyright": "Martin Lab, MIT licence"})
