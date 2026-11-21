"""Freeze, ad-hoc sign and archive a native COFdraw.app on a Mac."""
import platform
import plistlib
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from PIL import Image, ImageOps
from build_cofdraw_desktop import (
    REPO, VERSION, copy_licenses, freeze_desktop, sha, stage_scientific_data,
    validate_catalog, validate_raspa_resources, validate_synthesis_resources,
)


def main():
    if sys.platform != "darwin":
        raise SystemExit("Build on a Mac of the target architecture")
    arch = platform.machine()
    destination = REPO / "dist"
    archive = destination / f"COFdraw-macOS-{arch}.zip"
    if archive.exists():
        raise SystemExit(f"Output already exists: {archive}")
    destination.mkdir(exist_ok=True)
    validate_raspa_resources()
    validate_synthesis_resources()
    data = REPO / "data"
    with TemporaryDirectory(prefix="cofdraw-mac-") as temporary:
        work = Path(temporary)
        stage_scientific_data(data, work, validate_catalog(data))
        notices = work / "notices"
        copy_licenses(notices)
        shutil.copytree(REPO / "docs/third-party-licenses", notices, dirs_exist_ok=True)
        shutil.copy2(REPO / "LICENSE", notices / "COFdraw-LICENSE.txt")
        shutil.copy2(REPO / "src/cof_vision/cofdraw_xray_factors.py", notices)
        icon = work / "COFdraw.icns"
        with Image.open(REPO / "src/cof_vision/web/cofdraw-logo.png") as original:
            ImageOps.pad(original.convert("RGBA"), (1024, 1024), color=(0, 0, 0, 0)).save(icon)
        freeze_desktop(work, icon, [
            "--osx-bundle-identifier", "org.zzhenglab.cofdraw",
            "--hidden-import", "webview.platforms.cocoa",
            "--add-data", f"{work / 'data'}:data", "--add-data", f"{notices}:notices",
        ])
        app = work / "dist/COFdraw.app"
        # WKWebView enforces the page's script policy on the injected bridge.
        # pywebview 6.2.1 builds its wrappers with new Function. Install equivalent
        # closures in the bundled JS so the page can retain its strict CSP.
        bridge = app / "Contents/Resources/webview/js/api.js"
        js = bridge.read_text("utf-8")
        start = js.index("      // Define the function body")
        end = js.index("      window.pywebview._returnValuesCallbacks[funcName]", start)
        js = js[:start] + """      nestedObject[functionName] = (function (name) {
        return function () {
          var id = (Math.random() + "").substring(2);
          var args = Array.prototype.slice.call(arguments);
          return new Promise(function (resolve, reject) {
            window.pywebview._checkValue(name, resolve, reject, id);
            window.pywebview._jsApiCallback(name, args, id);
          });
        };
      })(funcName);
""" + js[end:]
        bridge.write_text(js, "utf-8")
        info = app / "Contents/Info.plist"
        metadata = plistlib.loads(info.read_bytes())
        metadata.update(CFBundleShortVersionString=VERSION, CFBundleVersion=VERSION)
        info.write_bytes(plistlib.dumps(metadata))
        subprocess.run(["codesign", "--force", "--sign", "-", str(app)], check=True)
        subprocess.run(["codesign", "--verify", "--deep", "--strict", str(app)], check=True)
        subprocess.run(["ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", str(app), str(archive)], check=True)
    archive.with_suffix(".sha256").write_text(f"{sha(archive)}  {archive.name}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
