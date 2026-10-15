"""Actual-browser mesh workspace checks. All field values are explicit UI fixtures.
Requires Python Playwright and Chromium. --in-memory rewrites only module/style
URLs to load the same assets on hosts that prohibit navigation; it does not test
served-asset CSP integration or durable localStorage. No security flags are disabled.
"""
import argparse
import base64
import copy
import json
import math
import re
from functools import lru_cache
from pathlib import Path
from playwright.sync_api import sync_playwright, expect
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--in-memory',action='store_true')
parser.add_argument('--url',default='http://127.0.0.1:7070')
parser.add_argument('--browser',default='/usr/bin/chromium')
parser.add_argument('--screenshots',type=Path)
args=parser.parse_args()

def make_fixture():
    # A synthetic curved sheet for UI geometry/field mapping only. Not a simulation.
    nx,ny=30,20
    p=[[i/nx*1.4-.7,j/ny-.5,.06*math.sin(i/nx*math.pi)*math.sin(j/ny*math.pi)] for j in range(ny+1) for i in range(nx+1)]
    tri=[]
    for j in range(ny):
        for i in range(nx):
            a=j*(nx+1)+i; tri.extend([[a,a+1,a+nx+2],[a,a+nx+2,a+nx+1]])
    time=[0,.15,.4,.8]
    temperatures=[[300+f*12+65*(i/nx)*(1-j/ny) for j in range(ny+1) for i in range(nx+1)] for f in range(4)]
    vectors=[[v for _ in p for v in [3,4+f,0]] for f in range(4)]
    cells=[[100*frame+(idx%2)*35 for idx in range(len(tri))] for frame in range(4)]
    displacement=[[v for x,y,z in p for v in [0,0,frame*.02*(x+.7)]] for frame in range(4)]
    return {'schema':'frankensim.observatory/mesh-v1','title':'UI fixture · not native physics',
      'source':{'model':'Explicit browser test fixture','material':'Synthetic UI values; no material prediction','notes':['These frames exercise browser rendering only. They are not a Rust solver result.']},
      'mesh':{'positions_m':p,'triangles':tri,'triangle_cells':list(range(len(tri))),'cell_count':len(tri),
      'source_cells':copy.deepcopy(tri),'source_types':[5]*len(tri),'time_s':time,'displacement_field':3,'fields':[
      {'label':'Temperature','unit':'K','association':'point','components':1,'values':temperatures},
      {'label':'Heat flux','unit':'W/m²','association':'point','components':3,'values':vectors},
      {'label':'Cell stress','unit':'Pa','association':'cell','components':1,'values':cells},
      {'label':'Displacement','unit':'m','association':'point','components':3,'values':displacement}]}}

@lru_cache(None)
def module(name):
    text=(ROOT/name).read_text()
    for dep in re.findall(r"from\s+['\"]\./([^'\"]+)['\"]",text):
        text=text.replace("'./"+dep+"'",repr(module(dep)))
    if name=='mesh-ui.mjs':
        css=base64.b64encode((ROOT/'mesh.css').read_bytes()).decode()
        text=text.replace("link.href='/mesh.css'", "link.href='data:text/css;base64,"+css+"'")
    return 'data:text/javascript;base64,'+base64.b64encode(text.encode()).decode()

fixture=make_fixture()
with sync_playwright() as p:
    browser=p.chromium.launch(executable_path=args.browser,headless=True)
    page=browser.new_page(viewport={'width':1600,'height':1250},reduced_motion='reduce')
    page.set_default_timeout(5000)
    errors=[];page.on('pageerror',lambda e: errors.append(str(e)))
    if args.in_memory:
        html=(ROOT/'index.html').read_text()
        html=re.sub(r'<link[^>]+href="/styles.css"[^>]*>','',html)
        html=re.sub(r'<script[^>]+src="/app.mjs"[^>]*></script>','',html)
        page.set_content(html)
        page.add_style_tag(content=(ROOT/'styles.css').read_text())
        # The in-memory document lacks loopback origin and secure-context UUIDs.
        # This is test plumbing, not production metadata or a simulated solver.
        page.evaluate("""() => {window.fetch=async()=>Response.json({ready:false,material:'synthetic'});
          if(!crypto.randomUUID)crypto.randomUUID=()=>Math.random().toString(36).slice(2);} """)
        page.add_script_tag(url=module('app.mjs'),type='module')
    else:
        page.goto(args.url)
    expect(page.locator('#mesh-inspector')).to_have_count(1)
    expect(page.locator('#engine-label')).to_have_text(re.compile('Local viewer mode|Native engine ready'))
    def load(value,name='mesh.json'):
        text=value if isinstance(value,str) else json.dumps(value)
        page.locator('#file-input').set_input_files({'name':name,'mimeType':'text/plain','buffer':text.encode()})
    load(fixture)
    expect(page.locator('body')).to_have_class(re.compile('mesh-mode'))
    expect(page.locator('#mesh-clock')).to_have_text('0 s')
    expect(page.locator('#metric-temperature')).to_have_text('300K')
    expect(page.locator('#metric-tension')).to_have_text('365K')
    expect(page.locator('#mesh-max')).to_contain_text('401')
    graphics=page.evaluate("document.getElementById('mesh-scene').getContext('webgl') ? 'WebGL' : 'software'")
    print('Renderer:',graphics)
    # Probe the original point on the final frame; colors keep the full-history range.
    page.locator('#mesh-entity').fill('0');page.locator('#mesh-entity').dispatch_event('change')
    page.locator('#mesh-frame').fill('3');page.locator('#mesh-frame').dispatch_event('input')
    expect(page.locator('#metric-heat')).to_have_text('336K')
    expect(page.locator('#mesh-max')).to_contain_text('401')
    page.locator('#mesh-field').select_option('1')
    expect(page.locator('#metric-heat')).to_contain_text('7.6158')
    page.locator('#mesh-component').select_option('1')
    expect(page.locator('#metric-heat')).to_have_text('7W/m²')
    page.locator('#mesh-field').select_option('2')
    page.locator('#mesh-entity').fill('1');page.locator('#mesh-entity').dispatch_event('change')
    expect(page.locator('#selected-layer')).to_have_text('Cell 1')
    expect(page.locator('#metric-heat')).to_have_text('335Pa')
    with page.expect_download() as event:page.locator('#export-csv').click()
    rows=Path(event.value.path()).read_text().splitlines()
    assert rows[1:]==['0,35','0.15,135','0.4,235','0.8,335'],rows
    # Display scaling changes GPU positions, not the source data or probe values.
    page.locator('#mesh-deformation').fill('5');page.locator('#mesh-deformation').dispatch_event('change')
    expect(page.locator('#scene-subtitle')).to_contain_text('displacement ×5')
    with page.expect_download() as event:page.locator('#export-json').click()
    saved=json.loads(Path(event.value.path()).read_text())
    assert saved['mesh']==fixture['mesh']
    page.locator('#mesh-deformation').fill('0');page.locator('#mesh-deformation').dispatch_event('change')
    page.locator('#mesh-field').select_option('0')
    page.locator('#mesh-frame').fill('0');page.locator('#mesh-frame').dispatch_event('input')
    # A surface hit changes the exact original probe index.
    page.locator('#mesh-scene').click(position={'x':page.locator('#mesh-scene').bounding_box()['width']/2,'y':210})
    assert page.locator('#mesh-entity').input_value()!='0'
    page.locator('#mesh-scene').focus();page.keyboard.press('ArrowLeft')
    page.locator('#reset-view').click()
    page.locator('#mesh-play').click()
    expect(page.locator('#mesh-clock')).to_have_text('0.8 s')
    expect(page.locator('#mesh-play')).to_have_text('▶ Play')
    page.locator('#trace').focus();page.keyboard.press('ArrowLeft')
    expect(page.locator('#mesh-clock')).to_have_text('0.4 s')
    page.locator('#pin-baseline').click()
    changed=copy.deepcopy(fixture);changed['title']='Second explicit UI fixture'
    for values in changed['mesh']['fields'][0]['values']:
        for i in range(len(values)): values[i]+=10
    load(changed)
    expect(page.locator('#channels')).to_contain_text('same source geometry')
    bad=copy.deepcopy(fixture);bad['mesh']['triangles'][0][0]=999999
    count=page.locator('.run-card').count();load(bad)
    expect(page.locator('#notification')).to_contain_text('Import refused')
    assert page.locator('.run-card').count()==count
    page.locator('#clear-baseline').click()
    # Static VTU: coordinates must be declared explicitly, Unit attrs remain unchanged.
    vtu='''<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="3" NumberOfCells="1"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1000 0 0 0 1000 0</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 2</DataArray><DataArray type="Int64" Name="offsets" format="ascii">3</DataArray><DataArray type="UInt8" Name="types" format="ascii">5</DataArray></Cells><PointData><DataArray type="Float64" Name="Temperature" Unit="K" format="ascii">310 320 330</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>'''
    load(vtu,'native-writer-format.vtu')
    expect(page.locator('#confirm-mesh-import')).to_be_disabled()
    page.locator('#import-length-unit').select_option('mm');page.locator('#confirm-mesh-import').click()
    expect(page.locator('#mesh-clock')).to_have_text('Static snapshot')
    expect(page.locator('#mesh-play')).to_be_disabled()
    expect(page.locator('#layer-details')).to_contain_text('1 × 1 × 0')
    # Signal/native and spatial data remain one library, not separate apps.
    load('time_s,pressure_pa\n0,1\n1,2','signal.csv')
    expect(page.locator('body')).not_to_have_class(re.compile('mesh-mode'))
    expect(page.locator('#signal-heading')).to_have_text('The physical signal')
    expect(page.locator('#metric-temperature')).to_have_text('—')
    load(fixture)
    page.locator('#mesh-frame').fill('2');page.locator('#mesh-frame').dispatch_event('input')
    page.locator('#mesh-entity').fill('220');page.locator('#mesh-entity').dispatch_event('change')
    with page.expect_download() as event:page.locator('#snapshot-image').click()
    assert Path(event.value.path()).read_bytes().startswith(b'\x89PNG')
    for width in [390,768,1024,1440,1600]:
        page.set_viewport_size({'width':width,'height':1250})
        assert page.evaluate('document.documentElement.scrollWidth')<=width,(width,page.evaluate('document.documentElement.scrollWidth'))
    page.locator('#theme').click();expect(page.locator('html')).to_have_attribute('data-theme','light')
    page.locator('#theme').click()
    if args.screenshots:
        args.screenshots.mkdir(parents=True,exist_ok=True)
        page.locator('#mesh-scene').focus()
        page.wait_for_timeout(5100)
        page.locator('#workspace').screenshot(path=str(args.screenshots/'mesh-workspace-ui-fixture.png'))
    page.locator('#mesh-new-experiment').click()
    expect(page.locator('body')).not_to_have_class(re.compile('mesh-mode'))
    expect(page.locator('#run')).to_be_visible()
    assert not errors,errors
    browser.close()
print('Mesh browser checks passed. In-memory assets: '+str(args.in_memory)+'. Inputs are explicit UI fixtures, not native physics.')
