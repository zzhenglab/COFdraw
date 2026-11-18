"""Windows desktop launcher; local drawing only, no OCR/model imports."""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import re
import socket
import sys
import threading
import time
from urllib.request import urlopen

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.exception_handlers import request_validation_exception_handler
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .drawing_api import create_drawing_router

VERSION = "0.8.6"
STATIC = Path(__file__).with_name("web")


def create_desktop_app(data_root):
    app = FastAPI(title="COFdraw", version=VERSION, docs_url=None, redoc_url=None)

    @app.exception_handler(RequestValidationError)
    async def validation_error(request: Request, exc: RequestValidationError):
        if request.url.path.rstrip('/') in {'/api/drawing/describe', '/api/drawing/describe/chat'}:
            # Pydantic's default errors include rejected input values. This
            # endpoint may contain provider credentials even in malformed input.
            return JSONResponse({'detail': 'Invalid description request. Check the text, mode, catalog selections and API settings.'}, status_code=422)
        return await request_validation_exception_handler(request, exc)

    @app.middleware("http")
    async def local_only(request: Request, call_next):
        if request.url.hostname not in {"127.0.0.1", "localhost", "testserver"}:
            return JSONResponse({"detail": "Local connections only"}, status_code=403)
        if request.method not in {"GET", "HEAD", "OPTIONS"}:
            origin = request.headers.get("origin")
            if origin and origin.rstrip("/") != str(request.base_url).rstrip("/"):
                return JSONResponse({"detail": "Cross-origin writes are disabled"}, status_code=403)
            if request.headers.get("sec-fetch-site") == "cross-site":
                return JSONResponse({"detail": "Cross-site writes are disabled"}, status_code=403)
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Content-Security-Policy"] = "default-src 'self'; img-src 'self' blob: data:; style-src 'self' 'unsafe-inline'; script-src 'self'; connect-src 'self'; object-src 'none'; frame-ancestors 'none'"
        return response

    @app.get("/")
    def index():
        return FileResponse(STATIC / "cofdraw.html", media_type="text/html")

    @app.get("/api/desktop")
    def about():
        return {"version": VERSION, "offline": True, "optional_llm_api": True,
                "default_mode": "offline_rules", "data_directory": str(Path(data_root).resolve()), "training": False}

    app.include_router(create_drawing_router(Path(data_root)))
    app.mount("/", StaticFiles(directory=STATIC), name="assets")
    return app


class DesktopBridge:
    def __init__(self, origin):
        self._origin = origin
        self._window = None

    def _check(self):
        if self._window is None or not str(self._window.get_current_url()).startswith(self._origin + "/"):
            raise ValueError("Save is available only inside the local COFdraw window")

    def _save(self, suggested, data):
        import webview
        self._check()
        name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(suggested))[:160]
        if Path(name).suffix.lower() not in {".png", ".svg", ".mol", ".sdf", ".cdxml", ".json", ".cif", ".xyz", ".pdb", ".csv", ".zip"}:
            raise ValueError("Unsupported export extension")
        selected = self._window.create_file_dialog(webview.FileDialog.SAVE, save_filename=name)
        if not selected:
            return {"saved": False, "cancelled": True}
        path = Path(selected if isinstance(selected, str) else selected[0])
        path.write_bytes(data)
        return {"saved": True, "path": str(path)}

    def save_export(self, relative_path, suggested_filename):
        self._check()
        if not re.fullmatch(r"/api/drawing/exports/[0-9a-f]{32}/(png|svg|mol|sdf|cdxml|json|cif|xyz|pdb|csv|zip)", str(relative_path)):
            raise ValueError("Not a generated COFdraw export")
        limit = 100_000_000 if str(relative_path).endswith('/zip') else 20_000_000
        with urlopen(self._origin + relative_path, timeout=30) as response:
            data = response.read(limit+1)
        if len(data) > limit:
            raise ValueError("Export exceeds the desktop save limit")
        return self._save(suggested_filename, data)

    def save_text(self, suggested_filename, content, mime_type=None):
        self._check()
        if not isinstance(content, str) or len(content) > 2_000_000:
            raise ValueError("Project is too large")
        return self._save(suggested_filename, content.encode("utf-8"))

    def save_image(self, suggested_filename, data_url):
        import base64
        self._check()
        if not isinstance(data_url,str) or len(data_url)>28_000_000 or not data_url.startswith('data:image/png;base64,'):
            raise ValueError('Expected a bounded PNG image')
        data=base64.b64decode(data_url.split(',',1)[1],validate=True)
        if not data.startswith(b'\x89PNG\r\n\x1a\n') or Path(suggested_filename).suffix.lower()!='.png':
            raise ValueError('Expected a PNG filename and image')
        return self._save(suggested_filename,data)


def main(argv=None):
    parser = argparse.ArgumentParser(description="COFdraw offline desktop")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--status-file", type=Path)
    parser.add_argument("--smoke-window", action="store_true")
    args = parser.parse_args(argv)
    home = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parents[3]
    data = (args.data_dir or home / "data").resolve()
    data.mkdir(parents=True, exist_ok=True)
    log = data / "cofdraw.log"
    logging.basicConfig(filename=log, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    # Windowed PyInstaller has no stdout/stderr; libraries expect writable streams.
    if sys.stdout is None:
        sys.stdout = log.open("a", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = log.open("a", encoding="utf-8")
    import uvicorn
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", args.port))
    port = sock.getsockname()[1]
    origin = f"http://127.0.0.1:{port}"
    server = uvicorn.Server(uvicorn.Config(create_desktop_app(data), host="127.0.0.1", port=port, log_config=None, access_log=False, loop="asyncio", http="h11", ws="none"))
    status = {"version": VERSION, "url": origin + "/", "data_directory": str(data), "pid": os.getpid(), "status": "starting"}
    def write_status():
        if args.status_file:
            args.status_file.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    write_status()
    if args.headless:
        status["status"] = "headless_server"; write_status()
        server.run(sockets=[sock])
        return
    worker = threading.Thread(target=server.run, kwargs={"sockets": [sock]}, daemon=True)
    worker.start()
    try:
        for _ in range(100):
            if server.started:
                break
            if not worker.is_alive():
                raise RuntimeError("The local drawing service failed to start")
            time.sleep(.1)
        if not server.started:
            raise RuntimeError("The local drawing service did not become ready")
        import webview
        bridge = DesktopBridge(origin)
        window = webview.create_window(f"COFdraw {VERSION}", origin + "/", js_api=bridge, width=1400, height=960, min_size=(980, 700), hidden=args.smoke_window, text_select=True, background_color="#f4f3ee")
        bridge._window = window
        status["status"] = "desktop_window"; write_status()
        def smoke():
            deadline = time.monotonic() + 40
            try:
                while time.monotonic() < deadline:
                    ready = window.evaluate_js("document.readyState === 'complete' && !!document.querySelector('main') && !!(window.pywebview && window.pywebview.api && window.pywebview.api.save_export)")
                    if ready:
                        status.update(status="desktop_smoke_passed", title=window.evaluate_js("document.title"), native_save_bridge=True)
                        write_status(); return
                    time.sleep(.2)
                raise RuntimeError("Desktop document did not become ready")
            except Exception as exc:
                status.update(status="desktop_smoke_failed", error=str(exc)); write_status()
            finally:
                window.destroy()
        webview.start(smoke if args.smoke_window else None, gui="edgechromium", debug=False, private_mode=True)
        if args.smoke_window and status.get('status') != 'desktop_smoke_passed':
            return 1
    except Exception as exc:
        logging.exception("Desktop startup failed")
        status.update(status="failed", error=str(exc)); write_status()
        if not args.smoke_window:
            import ctypes
            ctypes.windll.user32.MessageBoxW(0, "COFdraw could not open its desktop window.\n\n"
                "Extract the complete current release and keep COFdraw.exe.config, _internal, and data beside COFdraw.exe.\n\n"
                "The desktop also requires Microsoft WebView2 Runtime.\n\nDetails: " + str(log), "COFdraw", 0x10)
        return 1
    finally:
        server.should_exit = True
        worker.join(timeout=5)
        sock.close()


if __name__ == "__main__":
    main()
