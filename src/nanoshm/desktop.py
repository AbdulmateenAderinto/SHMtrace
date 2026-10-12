"""The desktop app: a local server plus a native window.

`nanoshm-app` starts uvicorn on a free port on localhost only, waits for it,
and opens a window on it. With no GUI toolkit available it prints the address
and keeps serving, so the same command works over SSH or in a container.
"""
from __future__ import annotations

import argparse
import socket
import sys
import threading
import time
import urllib.request

from . import __version__

HOST = "127.0.0.1"


def _free_port() -> int:
    with socket.socket() as s:
        s.bind((HOST, 0))
        return s.getsockname()[1]


def _port_taken(port: int) -> bool:
    with socket.socket() as s:
        try:
            s.bind((HOST, port))
        except OSError:
            return True
        return False


def _serve(port: int):
    import uvicorn

    from .webapp.main import app
    uvicorn.run(app, host=HOST, port=port, log_level="warning")


def _wait(url: str, timeout: float = 40.0) -> bool:
    end = time.time() + timeout
    while time.time() < end:
        try:
            with urllib.request.urlopen(url, timeout=1):
                return True
        except Exception:
            time.sleep(0.2)
    return False


def main():
    ap = argparse.ArgumentParser(description=f"nanoshm {__version__} desktop app")
    ap.add_argument("--port", type=int, default=0, help="0 picks a free port")
    ap.add_argument("--no-window", action="store_true", help="serve only, print the address")
    args = ap.parse_args()

    port = args.port or _free_port()
    if _port_taken(port):
        # otherwise _wait() would find whatever else is listening and open that
        print(f"port {port} is already in use; pick another or leave --port out",
              file=sys.stderr)
        return 1
    url = f"http://{HOST}:{port}/"
    threading.Thread(target=_serve, args=(port,), daemon=True).start()
    if not _wait(url):
        print("nanoshm could not start its local server", file=sys.stderr)
        return 1
    print(f"nanoshm {__version__} is running at {url}", flush=True)

    if args.no_window:
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            return 0

    try:
        import webview
    except ImportError:
        print("No window toolkit (pywebview) — open the address above in a browser.", flush=True)
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            return 0
    # pywebview ignores downloads unless told otherwise, which would leave every
    # button on the Export tab doing nothing
    webview.settings["ALLOW_DOWNLOADS"] = True
    # ?desktop=1 tells the page it is the desktop app, so it can offer "save to a
    # folder"; a query survives reloads, and is there before the page's script runs
    webview.create_window(f"nanoshm {__version__}", url + "?desktop=1", width=1280,
                          height=880, min_size=(900, 640))
    webview.start()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
