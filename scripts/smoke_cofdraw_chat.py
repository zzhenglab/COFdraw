"""Real local browser conversation, custom drawing and privacy regression.

The optional provider is an explicit local HTTP fixture, not a real model.
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
import traceback
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--base-url')
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    from playwright.sync_api import sync_playwright
    server = worker = None
    if args.base_url:
        parsed = urlsplit(args.base_url)
        if parsed.scheme != 'http' or parsed.hostname not in ('127.0.0.1', 'localhost', '::1') or parsed.username or parsed.password or parsed.path not in ('', '/') or parsed.query or parsed.fragment:
            parser.error('Use an HTTP loopback origin without credentials or a path.')
        base = args.base_url.rstrip('/')
    else:
        import uvicorn
        from cof_vision.cofdraw_desktop import create_desktop_app
        runtime = out / 'runtime-data'
        if runtime.exists():
            raise ValueError('Use a fresh output directory.')
        shutil.copytree(ROOT / 'data', runtime, ignore=shutil.ignore_patterns('exports', '*.log'))
        with socket.socket() as sock:
            sock.bind(('127.0.0.1', 0)); port = sock.getsockname()[1]
        server = uvicorn.Server(uvicorn.Config(create_desktop_app(runtime), host='127.0.0.1', port=port, log_level='warning'))
        worker = threading.Thread(target=server.run, daemon=True); worker.start()
        for _ in range(200):
            if server.started: break
            time.sleep(.05)
        assert server.started
        base = f'http://127.0.0.1:{port}'
    fixture = {'calls': 0, 'safe': False}
    token = 'chat-fixture-not-a-real-key'

    class Provider(BaseHTTPRequestHandler):
        def log_message(self, *_): pass

        def do_POST(self):
            raw = self.rfile.read(int(self.headers['Content-Length']))
            data = json.loads(raw)
            fixture['calls'] += 1
            fixture['safe'] = (token not in raw.decode() and 'PRIVATE_LOCAL_DRAWING' not in raw.decode()
                and 'layout_position' not in raw.decode() and 'O=Cc1cc(C=O)cc(C=O)c1' not in raw.decode()
                and self.headers.get('Authorization') == 'Bearer ' + token)
            prompt = json.loads(data['messages'][1]['content'])
            turn = json.loads(prompt['description'])
            assert set(turn) == {'current_request', 'current_catalog_selection', 'instruction'}
            proposal = {'selection': {'style_id': 'pastel_blocks', 'palette': 'pastel'}, 'explanation': 'Local fixture appearance selection.', 'unresolved': []}
            response = json.dumps({'choices': [{'message': {'content': json.dumps(proposal)}}]}).encode()
            self.send_response(200); self.send_header('Content-Type', 'application/json'); self.send_header('Content-Length', str(len(response))); self.end_headers(); self.wfile.write(response)

    provider = ThreadingHTTPServer(('127.0.0.1', 0), Provider)
    provider_thread = threading.Thread(target=provider.serve_forever, daemon=True); provider_thread.start()
    checks, errors, requests = [], [], []
    failure = None; page = None
    try:
        with sync_playwright() as p:
            executable = next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'), Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.exists()), None)
            browser = p.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={'width': 1280, 'height': 900}, bypass_csp=True)
            page.set_default_timeout(30000)
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.on('request', lambda request: requests.append(request.url) if request.url.endswith('/prepare') else None)
            page.goto(base)
            assert page.locator('#mode-manual').get_attribute('aria-selected') == 'true'
            page.locator('#mode-describe').click()
            page.wait_for_function('window.COFdrawDescribe?.state.catalog && window.COFdrawStudio.state.online && !window.COFdrawDescribe.state.busy')

            def idle():
                page.wait_for_function('!window.COFdrawDescribe.state.busy && !window.COFdrawStudio.state.busy')

            def prompt(text):
                page.locator('#describe-text').fill(text); page.locator('#describe-analyze').click(); idle()

            def follow(text):
                page.locator('#chat-text').fill(text); page.locator('#chat-send').click(); idle()

            def status():
                return page.evaluate('window.COFdrawDescribe.state.plan?.status')

            prompt('BD + 13C-Tfppy')
            assert status() == 'needs_input' and page.locator('#describe-monomer-b').input_value() == 'tfppy_13c'
            page.locator('button[data-monomer-choice="recovered_bd"]').click(); idle()
            assert status() == 'ready'
            assert {page.locator('#describe-monomer-'+slot).input_value() for slot in ('a','b')} == {'tfppy_13c','recovered_bd'}
            checks.append('Choosing an ambiguous BD by full name preserves its labelled partner and validates the pair')
            page.locator('#chat-new').click()
            prompt('Draw TFB and TAPB')
            assert status() == 'ready'
            follow('Make it pastel')
            assert status() == 'ready'
            assert page.locator('#describe-style').input_value() == 'pastel_blocks'
            assert page.locator('#chat-history .chat-turn').count() == 1
            follow('Change B to DETH')
            assert status() == 'ready'
            assert page.locator('#describe-linkage').input_value() == 'acylhydrazone'
            page.locator('#describe-create').click(); idle()
            assert page.evaluate('window.COFdrawStudio.state.result.audit.atoms') > 20
            assert 'successfully' in page.locator('#describe-status').inner_text()
            assert page.locator('#describe-create').is_disabled()
            checks.append('Session turns retain the pair, change appearance, replace B with DETH and render real acylhydrazone chemistry after review')

            # Actual local coordinate edit plus private metadata. A style turn
            # must use this existing graph, not reconstruct a catalog template.
            page.evaluate("""() => {const s=window.COFdrawStudio.state;s.recipe.instances[0].layout_position[0]+=.17;s.recipe.private_local_note='PRIVATE_LOCAL_DRAWING';s.sourceId=null;}""")
            page.locator('#chat-use-current').click()
            recipe = page.evaluate('JSON.stringify(window.COFdrawStudio.state.recipe)')
            before_prepares = len(requests)
            follow('Use wavy bonds')
            assert status() == 'ready'
            assert page.evaluate('window.COFdrawDescribe.state.plan.operation') == 'style_only'
            page.locator('#describe-create').click(); idle()
            assert page.evaluate('JSON.stringify(window.COFdrawStudio.state.recipe)') == recipe
            assert len(requests) == before_prepares
            assert page.evaluate('window.COFdrawStudio.state.result.audit.drawing_style.id') == 'mono_wave'
            checks.append('Current custom graph and edited coordinates remain byte-identical through a real style-only render with no preparation request')

            page.locator('#chat-use-current').click()
            def stale_response(route):
                actual = route.fetch()
                page.evaluate('window.COFdrawStudio.state.recipe.instances[0].layout_position[0]+=.09')
                route.fulfill(response=actual)
            page.route('**/api/drawing/describe/chat', stale_response)
            follow('Make it pastel')
            page.unroute('**/api/drawing/describe/chat', stale_response)
            assert page.evaluate('window.COFdrawDescribe.state.plan') is None
            assert page.locator('#describe-create').is_disabled()
            assert 'drawing changed' in page.locator('#describe-status').inner_text().lower()
            checks.append('A real response arriving after an in-flight drawing edit cannot enable Apply')

            page.locator('#chat-new').click()
            prompt('Draw TFB and TAPB with a 2 nm pore')
            assert status() == 'unsupported'
            follow('Make it pastel')
            assert status() != 'ready' and page.locator('#describe-create').is_disabled()
            assert 'physical pore' in page.locator('#describe-feedback').inner_text().lower()
            checks.append('A physical pore target remains unresolved across a later appearance request')

            page.locator('#chat-new').click()
            page.locator('#describe-provider > summary').click()
            page.locator('#describe-mode').select_option('llm')
            assert page.locator('#describe-service').input_value() == 'openai'
            assert page.locator('#describe-model').input_value().startswith('gpt-')
            page.locator('#describe-api-key').fill(token)
            page.locator('#describe-service').select_option('anthropic')
            assert page.locator('#describe-model').input_value().startswith('claude-')
            assert page.locator('#describe-api-key').input_value() == '' and page.locator('#describe-base-url').is_hidden()
            checks.append('Named GPT and Claude select their model defaults and clear the key when switching providers')
            page.locator('#describe-service').select_option('compatible')
            page.locator('#describe-base-url').fill(f'http://127.0.0.1:{provider.server_port}/v1')
            page.locator('#describe-model').fill('local-fixture')
            page.locator('#describe-api-key').fill(token)
            assert fixture['calls'] == 0
            page.locator('#chat-use-current').click()
            follow('Make it pastel')
            assert status() == 'ready', page.locator('#describe-status').inner_text()
            assert fixture['calls'] == 1 and fixture['safe']
            assert token not in page.evaluate('JSON.stringify(window.COFdrawDescribe.state.history)+JSON.stringify(localStorage)+JSON.stringify(sessionStorage)')
            page.locator('#describe-create').click(); idle()
            with page.expect_download() as saved:
                page.locator('#save-project').click()
            project = out / 'chat-project.json'; saved.value.save_as(project)
            assert token not in project.read_text('utf-8') and 'chat-history' not in project.read_text('utf-8')
            checks.append('Explicit optional API call sends only the new message/catalog IDs; private graph, credentials and session history stay out of provider content and project exports')

            follow('How do I check SA scores?')
            assert fixture['calls'] == 1 and page.evaluate('window.COFdrawDescribe.state.plan.operation') == 'help'
            page.locator('#chat-actions button').click()
            assert page.locator('#monomer-library-dialog').is_visible()
            assert page.locator('.library-sa').count() > 50
            page.keyboard.press('Escape')
            page.locator('#mode-describe').click()
            follow('How do I measure pores?')
            page.locator('#chat-actions button').click()
            assert page.locator('#science-porosity').is_visible()
            page.keyboard.press('Escape')
            assert fixture['calls'] == 1
            checks.append('Offline SA and pore guidance opens the actual library/geometry tools without a provider call or invented result')

            canvas_sizes = []
            for width, height in [(1280, 900), (1366, 768), (1920, 1080), (390, 844)]:
                page.set_viewport_size({'width': width, 'height': height})
                page.locator('.building-panel').evaluate('el => el.scrollTop=0')
                assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth + 1')
                box = page.locator('#canvas-viewport').bounding_box()
                assert abs(box['width'] - box['height']) <= 2, box
                assert not page.locator('.identity-more').evaluate('element => element.open')
                assert page.locator('#identity-string').is_visible()
                canvas_sizes.append({'viewport': [width, height], 'canvas': [round(box['width']), round(box['height'])]})
                page.screenshot(path=str(out / f'chat-{width}.png'), full_page=True)
            checks.append('Conversation fits1280/1366/1920/390px; canvas stays square, full notation stays visible and additional identity details stay collapsed: ' + json.dumps(canvas_sizes))
            page.locator('#chat-new').click()
            assert not page.evaluate('window.COFdrawDescribe.state.history.length')
            assert page.evaluate('window.COFdrawDescribe.state.context') is None
            assert page.evaluate('Boolean(window.COFdrawStudio.state.recipe)')
            checks.append('Start new clears transient chat context/history while preserving the drawing')
            assert not errors, errors
            browser.close()
    except Exception:
        failure = traceback.format_exc()
        if page:
            try: page.screenshot(path=str(out / 'failure.png'), full_page=True)
            except Exception: pass
    finally:
        provider.shutdown(); provider.server_close(); provider_thread.join(timeout=5)
        if server:
            server.should_exit = True; worker.join(timeout=10)
        report = {'passed': failure is None and not errors, 'checks': checks, 'javascript_errors': errors, 'failure': failure,
                  'provider_fixture': {'kind': 'local_http_fixture_not_a_model', **fixture}}
        (out / 'chat-smoke-report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps(report, indent=2))
    if not report['passed']: raise SystemExit(1)


if __name__ == '__main__':
    main()
