"""Optional real-Chromium interaction smoke test. Requires Python Playwright.

The API responses below are explicitly synthetic test doubles, not Rust physics.
Normal mode exercises served assets; --in-memory loads identical assets with only
module URLs rewritten when browser navigation is prohibited by the test host.
That mode does not establish browser/HTTP/CSP integration or persistent storage.
The Node server tests separately exercise real loopback HTTP and child processes.
"""
import argparse
import base64
import copy
import json
import math
import os
import re
from pathlib import Path
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[1]
args = argparse.ArgumentParser(description=__doc__)
args.add_argument('--url', default='http://127.0.0.1:7070')
args.add_argument('--in-memory', action='store_true')
args.add_argument('--browser', default=os.environ.get('CHROMIUM_PATH'))
args.add_argument('--screenshots', type=Path)
options = args.parse_args()
times = [i / 16000 for i in range(800)]
fixture = {'schema': 'frankensim.observatory/v1', 'title': 'UI test fixture · not native physics',
    'source': {'model': 'Browser integration fixture', 'material': 'Synthetic test values only',
        'notes': ['This input tests the viewer. It is not a Rust solver result.']},
    'inputs': {'voltage': .2, 'duration': .1, 'resistivity': 2, 'support': 'strain', 'strain': .0001, 'limit': .01},
    'section': {'kind': 'concentric-section', 'length_m': .65, 'radii_m': [.00025, .0005],
        'cold_temperature_k': [300, 300], 'hot_temperature_k': [307.3, 302.9], 'hot_stress_pa': [-4e6, 1.8e7]},
    'metrics': {'cold_tension_n': 15.7, 'hot_tension_n': 9.8, 'cold_axial_n': 157000, 'hot_axial_n': 153000,
        'stored_heat_j': .53, 'electrical_energy_j': .530000001, 'thermal_residual_j': 1e-9},
    'traces': {'time_s': times, 'series': [{'label': name, 'unit': 'Pa',
        'values': [.002 * math.sin(2 * math.pi * frequency * t) * math.exp(-t * 18) for t in times]}
        for name, frequency in [('Cold pressure', 330), ('Heated pressure', 270)]]}}

def module(path):
    return 'data:text/javascript;base64,' + base64.b64encode(path.read_bytes()).decode()

with sync_playwright() as playwright:
    browser = playwright.chromium.launch(headless=True, executable_path=options.browser)
    page = browser.new_page(viewport={'width': 1600, 'height': 1200}, reduced_motion='reduce')
    page.set_default_timeout(4000)
    errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    setup = """data => {
      const original = window.fetch.bind(window);
      window.qa = {ready:false, delay:20, reject:false, data};
      window.fetch = async (url, opts = {}) => {
        if (url === '/api/status') return Response.json({ready:qa.ready,material:'synthetic',token:'test-only'});
        if (url === '/api/run') return new Promise((resolve,reject) => {
          const timer=setTimeout(()=>resolve(Response.json(qa.reject ? {error:'Test-double refusal'} : {snapshot:qa.data}, {status:qa.reject?400:200})),qa.delay);
          opts.signal.addEventListener('abort',()=>{clearTimeout(timer);reject(new DOMException('Aborted','AbortError'));},{once:true});
        });
        return original(url,opts);
      };
    }"""
    if options.in_memory:
        html = (ROOT / 'index.html').read_text()
        html = re.sub(r'<link[^>]+href="/?styles.css"[^>]*>', '', html)
        html = re.sub(r'<script[^>]+src="/?app.mjs"[^>]*></script>', '', html)
        page.set_content(html)
        page.add_style_tag(content=(ROOT / 'styles.css').read_text())
        page.evaluate(setup, fixture)
        # about:blank is not a secure context. Production loopback has native UUIDs.
        page.evaluate("if(!crypto.randomUUID)crypto.randomUUID=()=>Math.random().toString(36).slice(2)")
        app = (ROOT / 'app.mjs').read_text().replace("'./protocol.mjs'", repr(module(ROOT / 'protocol.mjs')))
        app = app.replace("'./visuals.mjs'", repr(module(ROOT / 'visuals.mjs')))
        page.add_script_tag(content=app, type='module')
    else:
        page.add_init_script('(' + setup + ')(' + json.dumps(fixture) + ')')
        page.goto(options.url)
    expect(page.locator('#engine-label')).to_have_text('Local viewer mode')
    expect(page.locator('#run')).to_be_disabled()
    expect(page.locator('#metric-temperature')).to_have_text('—')

    def load(value, name='fixture.json'):
        contents = value if isinstance(value, str) else json.dumps(value)
        page.locator('#file-input').set_input_files({'name': name, 'mimeType': 'text/plain', 'buffer': contents.encode()})

    load(fixture)
    expect(page.locator('#scene-state')).to_have_text('IMPORTED · HOT PHYSICAL STATE')
    expect(page.locator('#metric-temperature')).to_have_text('307.3K')
    page.locator('#cold').click()
    expect(page.locator('#metric-temperature')).to_have_text('300K')
    page.locator('#hot').click()
    page.locator('#field').select_option('stress')
    expect(page.locator('#core-value')).to_contain_text('-4.00e+6 Pa')
    page.locator('#sleeve').click()
    expect(page.locator('#selected-layer')).to_have_text('Sleeve')
    page.locator('#scene').focus()
    old = page.locator('#scene').evaluate('c=>c.toDataURL()')
    page.keyboard.press('ArrowLeft')
    page.wait_for_timeout(40)
    assert page.locator('#scene').evaluate('c=>c.toDataURL()') != old
    page.locator('#reset-view').click()
    page.locator('#field').select_option('temperature')
    page.locator('#core').click()
    page.locator('#trace').focus()
    page.keyboard.press('ArrowRight')
    expect(page.locator('#cursor-value')).to_have_text('6.25e-5 s')
    page.locator('#zoom').select_option('4')
    page.locator('#pin-baseline').click()
    shifted = copy.deepcopy(fixture)
    shifted['title'] = 'Shifted UI test fixture'
    shifted['traces']['time_s'] = [t + 1 for t in times]
    load(shifted)
    expect(page.locator('#channels')).to_contain_text('outside recorded time range')
    page.locator('#clear-baseline').click()
    count = page.locator('.run-card').count()
    load('time_s,p_pa\n0,1\n0,2', 'invalid.csv')
    expect(page.locator('#notification')).to_contain_text('Import refused')
    assert page.locator('.run-card').count() == count
    load('time_s,pressure_pa,temperature_k\n0,1,300\n1,2,301', 'mixed.csv')
    expect(page.locator('#plot-tabs button')).to_have_count(2)
    page.locator('#plot-tabs').get_by_role('button', name='K', exact=True).click()
    expect(page.locator('#channels')).to_contain_text('300 K')
    assert 'pressure' not in page.locator('#channels').inner_text()
    expect(page.locator('#listen-cold')).to_be_disabled()
    load({**fixture, 'title': '<img src=x onerror=alert(1)>'})
    expect(page.locator('img')).to_have_count(0)
    page.locator('#command').click()
    page.get_by_placeholder('Search actions…').fill('evidence')
    page.keyboard.press('Enter')
    expect(page.locator('#dialog-title')).to_have_text('Evidence, not decoration')
    expect(page.locator('#dialog-content')).to_contain_text('not authenticated')
    page.keyboard.press('Escape')
    load(fixture)
    with page.expect_download() as exported:
        page.locator('#export-json').click()
    assert json.loads(Path(exported.value.path()).read_text())['title'] == fixture['title']
    with page.expect_download() as exported:
        page.locator('#export-csv').click()
    assert len(Path(exported.value.path()).read_text().splitlines()) == 801
    with page.expect_download() as exported:
        page.locator('#snapshot-image').click()
    assert Path(exported.value.path()).read_bytes().startswith(b'\x89PNG')
    page.locator('#listen-cold').click()
    expect(page.locator('#audio-level')).to_contain_text('0 clipped samples')
    expect(page.locator('#stop-audio')).to_be_disabled()
    for width in [390, 768, 1024, 1440, 1600]:
        page.set_viewport_size({'width': width, 'height': 1200})
        assert page.evaluate('document.documentElement.scrollWidth') <= width
    page.locator('#theme').click()
    expect(page.locator('html')).to_have_attribute('data-theme', 'light')
    page.locator('#theme').click()
    if options.screenshots:
        options.screenshots.mkdir(parents=True, exist_ok=True)
        page.wait_for_timeout(5100)  # Let the transient import notification expire.
        page.screenshot(path=str(options.screenshots / 'observatory-test-input.png'), full_page=True)
    page.evaluate('qa.ready=true')
    page.locator('#engine-status').click()
    expect(page.locator('#run')).to_be_enabled()
    page.locator('#run').click()
    expect(page.locator('#notification')).to_contain_text('Accepted native result')
    # This only verifies the UI state transition with an explicit API test double.
    count = page.locator('.run-card').count()
    page.evaluate('qa.reject=true')
    page.locator('#run').click()
    expect(page.locator('#notification')).to_contain_text('Test-double refusal')
    assert page.locator('.run-card').count() == count
    page.evaluate('qa.reject=false;qa.delay=10000')
    page.locator('#run').click()
    expect(page.locator('#cancel')).to_be_visible()
    page.locator('#cancel').click()
    expect(page.locator('#notification')).to_contain_text('cancelled')
    assert page.locator('.run-card').count() == count
    if not options.in_memory:
        page.reload()
        expect(page.locator('.run-card').first).to_contain_text('RESTORED')
    assert not errors, errors
    browser.close()
print('Browser interaction checks passed (' + ('in-memory assets; HTTP/CSP and persistence not covered' if options.in_memory else 'served assets') + '). API fixture is not a Rust execution.')
