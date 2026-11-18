"""Exercise Describe in a browser, including an explicitly local model fixture.

The fixture tests HTTP transport/validation, not the quality of a real model.
No external model service is contacted. Works with source or --base-url frozen app.
"""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import argparse
import json
from pathlib import Path
import shutil
import socket
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--base-url')
    args = parser.parse_args()
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    from playwright.sync_api import sync_playwright
    server = worker = None
    if args.base_url:
        base = args.base_url.rstrip('/')
    else:
        import uvicorn
        from cof_vision.cofdraw_desktop import create_desktop_app
        data = output/'runtime-data'
        shutil.copytree(ROOT/'data', data, dirs_exist_ok=True, ignore=shutil.ignore_patterns('exports', '*.log'))
        with socket.socket() as sock:
            sock.bind(('127.0.0.1', 0)); port = sock.getsockname()[1]
        server = uvicorn.Server(uvicorn.Config(create_desktop_app(data), host='127.0.0.1', port=port, log_level='warning'))
        worker = threading.Thread(target=server.run, daemon=True); worker.start()
        for _ in range(200):
            if server.started: break
            time.sleep(.05)
        if not server.started: raise RuntimeError('Description test service failed to start')
        base = f'http://127.0.0.1:{port}'
    token = 'cofdraw-fixture-token-not-a-real-credential'
    fixture = {'count': 0, 'auth_ok': False, 'compact': False, 'selection': {}, 'fail': False}

    class Provider(BaseHTTPRequestHandler):
        def log_message(self, *_): pass

        def do_POST(self):
            fixture['count'] += 1
            fixture['auth_ok'] = self.headers.get('Authorization') == 'Bearer '+token
            raw = self.rfile.read(int(self.headers.get('Content-Length', '0')))
            payload = json.loads(raw)
            encoded = json.dumps(payload)
            fixture['compact'] = (len(raw) < 196608 and token not in encoded
                                  and 'PRIVATE_PROJECT_NOT_FOR_PROVIDER' not in encoded
                                  and 'O=Cc1cc(C=O)cc(C=O)c1' not in encoded
                                  and self.path == '/v1/chat/completions')
            if fixture['fail']:
                self.send_response(503); self.end_headers(); self.wfile.write(b'Local fixture failure'); return
            selection = {key: value for key, value in fixture['selection'].items() if value is not None}
            proposal = {'selection': selection, 'explanation': 'Local test fixture selects existing catalog IDs.', 'unresolved': []}
            response = json.dumps({'choices': [{'message': {'content': json.dumps(proposal)}}]}).encode()
            self.send_response(200); self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(response))); self.end_headers(); self.wfile.write(response)

    provider = ThreadingHTTPServer(('127.0.0.1', 0), Provider)
    provider_thread = threading.Thread(target=provider.serve_forever, daemon=True); provider_thread.start()
    checks, errors = [], []
    failure = None
    try:
        with sync_playwright() as playwright:
            executable = next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'), Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.is_file()), None)
            browser = playwright.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={'width': 1280, 'height': 900}, bypass_csp=True)
            page.set_default_timeout(30000); page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(base)
            assert page.locator('#mode-manual').get_attribute('aria-selected') == 'true'
            page.locator('#mode-describe').click()
            page.wait_for_function('window.COFdrawDescribe?.state.catalog && window.COFdrawStudio.state.online && !window.COFdrawDescribe.state.busy')
            assert page.locator('#describe-panel').is_visible() and not page.locator('#manual-panel').is_visible()
            assert page.locator('#describe-create').is_disabled()
            page.screenshot(path=str(output/'describe-start-1280.png'))
            checks.append('Manual builder starts selected and Describe opens without generating a drawing')

            def show_inputs():
                if page.locator('#describe-edit-text').is_visible():
                    page.locator('#describe-edit-text').click()

            def demo(name):
                show_inputs()
                page.locator(f'#describe-demos [data-demo="{name}"]').click()
                page.locator('#describe-analyze').click()
                page.wait_for_function('!window.COFdrawDescribe.state.busy')

            demo('imine')
            assert page.evaluate('window.COFdrawDescribe.state.plan.status') == 'ready'
            fixture['selection'] = page.evaluate('window.COFdrawDescribe.state.plan.selection')
            assert page.locator('#describe-monomer-a').input_value() and page.locator('#describe-monomer-b').input_value()
            page.wait_for_function('document.querySelectorAll("#describe-monomer-previews img").length===2')
            page.locator('.building-panel').evaluate('el=>el.scrollTop=0')
            create_box=page.locator('#describe-create').bounding_box()
            assert create_box and create_box['y']+create_box['height']<=900, 'Ready drawing action must fit1280x900'
            page.screenshot(path=str(output/'describe-review-1280.png'), full_page=True)
            page.locator('#describe-create').click()
            page.wait_for_function('!window.COFdrawDescribe.state.busy && window.COFdrawStudio.state.result && !window.COFdrawStudio.state.busy')
            assert page.evaluate('window.COFdrawStudio.state.sourceId') is None
            assert page.evaluate('window.COFdrawStudio.state.result.audit.atoms') > 20
            assert page.locator('#exports button').count() >= 5
            page.evaluate('window.scrollTo(0,0)'); page.screenshot(path=str(output/'describe-drawing-1280.png'))
            checks.append('Local prose resolves exact monomers/template and renders real editable chemistry')

            show_inputs();page.locator('#describe-text').fill('Draw a COF')
            assert page.locator('#describe-create').is_disabled() and page.evaluate('window.COFdrawDescribe.state.plan') is None
            show_inputs();page.locator('#describe-analyze').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy')
            assert page.evaluate('window.COFdrawDescribe.state.plan.status') == 'needs_input'
            assert page.locator('#describe-create').is_disabled()
            demo('unsupported')
            assert page.evaluate('window.COFdrawDescribe.state.plan.status') == 'unsupported'
            assert page.locator('#describe-create').is_disabled()
            checks.append('Edited text invalidates readiness; ambiguous and unsupported requests cannot draw')

            demo('fcof'); page.locator('#describe-create').click()
            page.wait_for_function('!window.COFdrawDescribe.state.busy && !window.COFdrawStudio.state.busy')
            assert page.evaluate('window.COFdrawStudio.state.sourceId') == 'f-cof-1-source-v1'
            assert page.evaluate('window.COFdrawStudio.state.result.audit.atoms') == 90
            checks.append('Named source description preserves registered graph identity and source style')

            show_inputs();page.locator('#describe-demos [data-demo="imine"]').click()
            page.locator('#describe-provider summary').click()
            page.locator('#describe-mode').select_option('llm')
            page.locator('#describe-service').select_option('compatible')
            page.locator('#describe-base-url').fill(f'http://127.0.0.1:{provider.server_port}/v1')
            page.locator('#describe-model').fill('local-fixture')
            page.locator('#describe-api-key').fill(token)
            page.evaluate('window.COFdrawStudio.state.recipe.private_project_marker="PRIVATE_PROJECT_NOT_FOR_PROVIDER"')
            assert fixture['count'] == 0  # Settings, mode and prompt changes send nothing.
            page.locator('#describe-analyze').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy')
            assert page.evaluate('window.COFdrawDescribe.state.plan?.status') == 'ready', page.locator('#describe-status').inner_text()
            assert fixture['count'] == 1 and fixture['auth_ok'] and fixture['compact']
            assert 'Local test fixture selects existing catalog IDs.' in page.locator('#describe-evidence-content').text_content()
            show_inputs();page.locator('#describe-analyze').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy')
            assert fixture['count'] == 2 and page.evaluate('window.COFdrawDescribe.state.plan?.status') == 'ready'
            page.locator('#describe-adjust summary').click();page.locator('#describe-style').select_option('mono_wave')
            assert page.locator('#describe-create').is_disabled()
            assert page.locator('#describe-analyze').inner_text() == 'Check revised choices'
            page.locator('#describe-analyze').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy')
            assert fixture['count'] == 2, 'Manual review must not contact the provider again'
            assert page.evaluate('window.COFdrawDescribe.state.plan?.status') == 'ready', page.locator('#describe-status').inner_text()
            page.locator('#describe-create').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy && !window.COFdrawStudio.state.busy')
            assert page.evaluate('window.COFdrawStudio.state.result.audit.drawing_style.id') == 'mono_wave'
            checks.append('Local model fixture verifies real browser/backend/provider transport, ID validation, rendering and free local selection review')
            with page.expect_download() as download:
                page.locator('#save-project').click()
            project_path = output/'describe-project.json'; download.value.save_as(project_path)
            assert token not in project_path.read_text('utf-8')
            assert token not in page.evaluate('JSON.stringify(localStorage)+JSON.stringify(sessionStorage)')
            checks.append('Ephemeral API key is absent from project exports and browser storage')

            fixture['fail'] = True
            show_inputs();page.locator('#describe-analyze').click(); page.wait_for_function('!window.COFdrawDescribe.state.busy')
            assert page.locator('#describe-create').is_disabled()
            assert page.evaluate('window.COFdrawDescribe.state.plan') is None
            assert page.locator('#describe-status').evaluate('el=>el.classList.contains("error")')
            assert token not in page.locator('#describe-status').inner_text()
            checks.append('Provider errors clear the ready plan and show a credential-safe failure')

            page.locator('#describe-api-key').fill(''); page.locator('#describe-mode').select_option('rules')
            page.locator('#describe-provider summary').click()
            page.locator('#mode-manual').click(); assert page.locator('#manual-panel').is_visible()
            page.locator('#mode-manual').focus(); page.keyboard.press('ArrowLeft')
            assert page.locator('#describe-panel').is_visible() and page.locator('#mode-describe').get_attribute('aria-selected') == 'true'
            checks.append('Manual builder remains accessible; starting modes support keyboard navigation')
            for width, height in [(1280, 1000), (1024, 900), (390, 844)]:
                page.set_viewport_size({'width': width, 'height': height}); page.evaluate('window.scrollTo(0,0);document.querySelector(".building-panel").scrollTop=0;document.querySelector(".canvas-panel").scrollTop=0')
                page.screenshot(path=str(output/f'describe-responsive-{width}.png'), full_page=True)
                assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1'), width
            checks.append('Describe, canvas and accessible settings fit1280/1024/390px without horizontal overflow')
            assert not errors, errors
            browser.close()
    except Exception as exc:
        failure = f'{type(exc).__name__}: {exc}'; raise
    finally:
        provider.shutdown(); provider.server_close(); provider_thread.join(timeout=5)
        if server:
            server.should_exit = True; worker.join(timeout=10)
        report = {'checks': checks, 'javascript_errors': errors, 'failure': failure,
                  'provider_test': {'kind': 'local_deterministic_http_fixture_not_real_model', 'request_count': fixture['count'],
                                    'authorization_received': fixture['auth_ok'], 'compact_catalog_without_project_or_credentials': fixture['compact']}}
        (output/'describe-smoke-report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps({'passed': len(checks), 'failure': failure, 'javascript_errors': errors}))


if __name__ == '__main__':
    main()
