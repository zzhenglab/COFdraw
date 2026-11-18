"""Verify the real measured diameter overlay and annotated exports.

Uses a local backend (or --base-url), with no external network/provider calls.
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import socket
import sys
import threading
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--base-url')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    from playwright.sync_api import sync_playwright
    from PIL import Image
    server = worker = None
    if args.base_url:
        base = args.base_url.rstrip('/')
    else:
        from cof_vision.cofdraw_desktop import create_desktop_app
        import uvicorn
        runtime = output / 'runtime-data'
        shutil.copytree(ROOT / 'data', runtime, dirs_exist_ok=True, ignore=shutil.ignore_patterns('exports', '*.log'))
        with socket.socket() as sock:
            sock.bind(('127.0.0.1', 0))
            port = sock.getsockname()[1]
        server = uvicorn.Server(uvicorn.Config(create_desktop_app(runtime), host='127.0.0.1', port=port, log_level='warning'))
        worker = threading.Thread(target=server.run, daemon=True)
        worker.start()
        for _ in range(200):
            if server.started:
                break
            time.sleep(.05)
        if not server.started:
            raise RuntimeError('Local test service did not start')
        base = f'http://127.0.0.1:{port}'
    checks, errors, failure = [], [], None
    try:
        with sync_playwright() as p:
            executable = next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),
                Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.is_file()), None)
            browser = p.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={'width': 1280, 'height': 1000}, bypass_csp=True)
            page.set_default_timeout(30000)
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(base)
            page.wait_for_function('window.COFdrawStudio?.state.online && !document.getElementById("build").disabled')
            page.locator('#mode-manual').click()
            page.locator('#build').click()
            page.wait_for_function('window.COFdrawStudio.state.result && !window.COFdrawStudio.state.busy')
            page.locator('#measure-pore').click()
            page.locator('#pore-measurement-overlay').wait_for(state='visible')
            record = page.evaluate('''() => {
              const r=COFdrawMeasure.active;
              return {result:r.result,layout:r.render.audit.drawing_style.layout,
                line:r.group.querySelector('line').outerHTML,canvas:r.canvasSize};
            }''')
            line = ET.fromstring(record['line'])
            pixels = [[float(line.get('x1')), float(line.get('y1'))], [float(line.get('x2')), float(line.get('y2'))]]
            result, layout = record['result'], record['layout']
            scale = layout['pixels_per_unit']
            for actual, expected in zip(pixels, result['annotation']['endpoints']):
                expected_pixels = [layout['atom_pixels'][0][0] + scale*(expected[0]-layout['atom_coordinates'][0][0]),
                                   layout['atom_pixels'][0][1] - scale*(expected[1]-layout['atom_coordinates'][0][1])]
                assert math.dist(actual, expected_pixels) < 1e-7
            diameter = math.dist(*pixels) / scale * result['annotation']['angstrom_per_drawing_unit']
            assert abs(diameter - result['diameter_angstrom']) < 1e-8
            assert f"{result['diameter_nm']:.2f} nm" in page.locator('#pore-measurement-overlay').text_content()
            checks.append('Visible double-ended arrow matches exact backend diameter and nm label')
            page.locator('#canvas-viewport').scroll_into_view_if_needed()
            page.screenshot(path=str(output / 'pore-arrow-1280.png'))
            for width, height in [(1024, 900), (390, 844)]:
                page.set_viewport_size({'width': width, 'height': height})
                page.wait_for_timeout(100)
                assert page.evaluate('''() => {
                  const a=document.getElementById('pore-measurement-overlay').getBoundingClientRect();
                  const b=document.getElementById('framework-image').getBoundingClientRect();
                  return ['x','y','width','height'].every(k=>Math.abs(a[k]-b[k])<1);
                }''')
            checks.append('Overlay preserves the image viewport at 1024 px and phone widths')
            page.set_viewport_size({'width': 1280, 'height': 1000})
            with page.expect_download() as download:
                page.locator('#save-pore-svg').click()
            svg_path = output / 'pore-annotated.svg'
            download.value.save_as(svg_path)
            svg = ET.fromstring(svg_path.read_text('utf-8'))
            group = next(element for element in svg.iter() if element.get('data-cofdraw-pore-measurement') == '1')
            exported_line = next(element for element in group if element.tag.endswith('line'))
            assert all(exported_line.get(key) == line.get(key) for key in ('x1', 'y1', 'x2', 'y2'))
            metadata = next(element for element in svg.iter() if element.get('data-cofdraw-measurement') == '1')
            provenance = json.loads(metadata.text)
            assert provenance['scope'] == result['scope'] and provenance['sampling_lower_bound'] is True
            checks.append('Annotated SVG preserves exact arrow endpoints plus method and scope metadata')
            with page.expect_download() as download:
                page.locator('#save-pore-png').click()
            png_path = output / 'pore-annotated.png'
            download.value.save_as(png_path)
            image = Image.open(png_path).convert('RGB')
            assert list(image.size) == record['canvas']
            import numpy as np
            rgb = np.asarray(image)
            assert ((rgb[:,:,0] < 80) & (rgb[:,:,1] > 70) & (rgb[:,:,1] < 130) & (rgb[:,:,2] > 80) & (rgb[:,:,2] < 150)).sum() > 100
            checks.append('Annotated PNG contains the measured arrow at the exported canvas size')
            original = page.evaluate('''async () => await (await fetch('.'+COFdrawStudio.state.result.downloads.svg)).text()''')
            assert 'data-cofdraw-pore-measurement' not in original
            checks.append('Annotation leaves ordinary figure and chemical exports unchanged')
            page.locator('#appearance-options').evaluate('element=>element.open=true')
            page.locator('[data-style="perspective"]').click()
            page.wait_for_function('!COFdrawStudio.state.busy && COFdrawStudio.state.result.audit.drawing_style.id==="perspective"')
            assert page.locator('#pore-measurement-overlay').is_hidden()
            page.locator('#measure-pore').click()
            page.wait_for_function('document.getElementById("pore-annotation-note")?.textContent.includes("flat chemical style")')
            assert page.locator('#pore-measurement-overlay').is_hidden()
            checks.append('Projected style clears the arrow and explains why a flat style is required')
            page.locator('[data-style="standard"]').click()
            page.wait_for_function('!COFdrawStudio.state.busy && COFdrawStudio.state.result.audit.drawing_style.id==="standard"')
            page.locator('#measure-pore').click()
            page.locator('#pore-measurement-overlay').wait_for(state='visible')
            page.locator('#recipe-json').evaluate('element=>{element.value+=" ";element.dispatchEvent(new Event("input",{bubbles:true}));}')
            assert page.locator('#pore-measurement-overlay').is_hidden()
            assert page.locator('#save-pore-svg').count() == 0
            checks.append('Pending JSON geometry edits remove annotation and annotated-save actions')
            # Real backend response, made stale only while in flight.
            page.evaluate('COFdrawStudio.state.jsonPending=false;COFdrawStudio.state.dirty=false;COFdrawStudio.sync()')
            def late_result(route):
                response = route.fetch()
                page.evaluate('COFdrawStudio.state.dirty=true;COFdrawStudio.sync()')
                route.fulfill(response=response)
            page.route('**/api/drawing/science/measure', late_result)
            page.locator('#measure-pore').click()
            page.wait_for_timeout(1000)
            assert page.locator('#pore-measurement-overlay').is_hidden(), 'Late response restored the arrow'
            assert page.locator('#pore-measurement').is_hidden(), 'Late response restored the measurement card'
            checks.append('A real measurement response arriving after an edit cannot restore stale output')
            assert not errors, errors
            browser.close()
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
    finally:
        if server:
            server.should_exit = True
            worker.join(10)
        report = {'checks': checks, 'javascript_errors': errors, 'failure': failure,
                  'scope': 'Local backend and browser regression; no external network or paid provider calls'}
        (output / 'measure-smoke-report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(report, indent=2))
    return 1 if failure is not None else 0


if __name__ == '__main__':
    raise SystemExit(main())
