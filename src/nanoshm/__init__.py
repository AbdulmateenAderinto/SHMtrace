"""The old name of SHMtrace, kept so existing scripts that `import nanoshm` still run.

`import nanoshm` and `from nanoshm.x import y` resolve to `shmtrace` and its modules,
with a DeprecationWarning. New code should import shmtrace.
"""
import importlib
import pkgutil
import sys
import warnings

import shmtrace

warnings.warn("nanoshm is now SHMtrace: import shmtrace instead", DeprecationWarning, stacklevel=2)

for _m in pkgutil.walk_packages(shmtrace.__path__, "shmtrace."):
    _old = "nanoshm" + _m.name[len("shmtrace"):]
    if _old not in sys.modules:
        try:
            sys.modules[_old] = importlib.import_module(_m.name)
        except Exception:            # optional extras (webapp, pysam) may be missing
            pass

sys.modules[__name__] = shmtrace
