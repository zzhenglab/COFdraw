from pathlib import Path
import pytest
from fastapi.testclient import TestClient
from cof_vision.cofdraw_desktop import create_desktop_app, DesktopBridge

def test_desktop_serves_offline_assets_and_rejects_foreign_origin(tmp_path):
    with TestClient(create_desktop_app(tmp_path)) as client:
        page = client.get("/")
        assert page.status_code == 200 and "cofdraw.js" in page.text
        assert "script-src 'self'" in page.headers["content-security-policy"]
        assert client.get("/api/desktop").json()["training"] is False
        assert client.get("/cofdraw.js").status_code == 200
        assert client.get("/", headers={"host": "untrusted.example"}).status_code == 403
        assert client.post("/api/drawing/monomer", json={"smiles": "CC"}, headers={"origin": "https://untrusted.example"}).status_code == 403
        assert client.post("/api/drawing/monomer", json={"smiles": "CC"}, headers={"origin": "http://testserver"}).status_code == 200

def test_desktop_bridge_rejects_untrusted_window_and_non_export_paths():
    bridge = DesktopBridge("http://127.0.0.1:12345")
    with pytest.raises(ValueError):
        bridge.save_export("/api/drawing/exports/" + "a" * 32 + "/svg", "drawing.svg")
    class Window:
        def get_current_url(self): return "http://127.0.0.1:12345/"
    bridge._window = Window()
    for path in ("https://example.org/drawing.svg", "/../private.txt", "/api/drawing/exports/" + "a" * 32 + "/exe"):
        with pytest.raises(ValueError):
            bridge.save_export(path, "drawing.svg")
