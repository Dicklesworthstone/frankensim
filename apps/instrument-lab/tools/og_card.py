"""Share cards (og.png 1200x630, twitter.png 1200x600) painted with real fs-vkplate fields.

usage:
  ogfield DIR/fields.json            # crates/fs-vkplate example (OG_TC/OG_TS = snapshot times, s)
  ffmpeg -i public/audio/kit-groove-ours.m4a -f f32le -ac 1 -ar 8000 DIR/groove.f32
  uv run --with numpy python tools/og_card.py DIR apps/instrument-lab
  then screenshot DIR/og.html and DIR/twitter.html at 1200x630 / 1200x600 (headless Chrome,
  --allow-file-access-from-files) into public/og.png and public/twitter.png.
"""
import json, sys, math
from pathlib import Path
import numpy as np
S = Path(sys.argv[1]); L = Path(sys.argv[2])
fields = json.loads((S / "fields.json").read_text())
# show the elastic waves: remove each surface's rigid motion (best-fit plane a + b x + c y of the velocity)
for f in fields.values():
    X, Y, V = (np.array(f[k]) for k in ("x", "y", "v"))
    A = np.c_[np.ones_like(X), X, Y]
    coef, *_ = np.linalg.lstsq(A, V, rcond=None)
    f["v"] = [float(v) for v in (V - A @ coef)]
x = np.frombuffer((S / "groove.f32").read_bytes(), "<f4")
bins = 300
seg = len(x) // bins
env = [float(np.abs(x[i * seg:(i + 1) * seg]).max()) for i in range(bins)]
m = max(env) or 1
env = [round(e / m, 3) for e in env]
F = L / "node_modules/@fontsource-variable/fraunces/files"
I = L / "node_modules/@fontsource-variable/inter/files"
J = L / "node_modules/@fontsource/jetbrains-mono/files"
html = r"""<!doctype html><html><head><meta charset="utf-8"><style>
@font-face { font-family: Fraunces; src: url("file://%FR%/fraunces-latin-full-normal.woff2"); font-style: normal; font-weight: 100 900; }
@font-face { font-family: Fraunces; src: url("file://%FR%/fraunces-latin-full-italic.woff2"); font-style: italic; font-weight: 100 900; }
@font-face { font-family: Inter; src: url("file://%IN%/inter-latin-wght-normal.woff2"); font-weight: 100 900; }
@font-face { font-family: Mono; src: url("file://%JB%/jetbrains-mono-latin-500-normal.woff2"); font-weight: 500; }
html, body { margin: 0; background: #07080b; }
#card { position: relative; width: 1200px; height: %H%px; overflow: hidden; background:
  radial-gradient(900px 520px at 78% 38%, rgba(217,169,79,0.16), rgba(0,0,0,0) 62%),
  radial-gradient(700px 420px at 92% 92%, rgba(99,210,192,0.10), rgba(0,0,0,0) 60%),
  radial-gradient(800px 600px at 0% 0%, rgba(140,110,60,0.10), rgba(0,0,0,0) 55%), #07080b; }
#art { position: absolute; inset: 0; }
.txt { position: absolute; left: 64px; top: 58px; width: 560px; color: #f1eadb; }
.eyebrow { font: 500 15px/1 Mono, monospace; letter-spacing: 0.22em; color: #d9a94f; text-transform: uppercase; display: flex; align-items: center; gap: 12px; }
.eyebrow i { display: inline-block; width: 34px; height: 1px; background: #d9a94f; }
h1 { font-family: Fraunces, serif; font-weight: 560; font-size: 92px; line-height: 0.95; letter-spacing: -0.02em; margin: 30px 0 0; font-variation-settings: "SOFT" 30, "opsz" 144; }
h1 em { font-style: italic; color: #e2b45f; font-weight: 520; }
.sub { font: 400 25px/1.38 Inter, sans-serif; color: rgba(241,234,219,0.78); margin-top: 26px; width: 520px; }
.chips { position: absolute; left: 64px; bottom: %CHIPB%px; display: flex; gap: 10px; flex-wrap: wrap; width: 640px; }
.chip { font: 500 12.5px/1 Mono, monospace; letter-spacing: 0.06em; color: rgba(241,234,219,0.86); border: 1px solid rgba(241,234,219,0.18); border-radius: 999px; padding: 8px 13px; background: rgba(255,255,255,0.02); }
.chip.b { color: #d9a94f; border-color: rgba(217,169,79,0.45); }
.chip.g { color: #63d2c0; border-color: rgba(99,210,192,0.45); }
.url { position: absolute; right: 56px; bottom: %URLB%px; font: 500 16px/1 Mono, monospace; letter-spacing: 0.08em; color: rgba(241,234,219,0.7); }
.url b { color: #f1eadb; font-weight: 500; }
.cap { position: absolute; font: 500 11.5px/1.3 Mono, monospace; letter-spacing: 0.08em; color: rgba(241,234,219,0.55); text-transform: uppercase; }
</style></head><body><div id="card">
<canvas id="art" width="2400" height="%H2%" style="width:1200px;height:%H%px"></canvas>
<div class="txt"><div class="eyebrow"><i></i>FrankenSim · Struck</div>
<h1>Nothing here was <em>recorded.</em></h1>
<div class="sub">A concert grand and a full drum kit, computed from strings, spruce, film and bronze.</div></div>
<div class="chips"><span class="chip b">Grand piano</span><span class="chip">Kick</span><span class="chip">Snare</span><span class="chip">Toms</span><span class="chip">Hi-hat</span><span class="chip">Crash</span><span class="chip">Ride</span></div>
<div class="url"><b>frankensim.org</b>/music</div>
<div class="cap" style="left:%CAPX%px; top:%CAPY%px">18-inch crash · velocity field 2.2 ms after the stick</div>
</div>
<script>
const F = %FIELDS%; const ENV = %ENV%; const H = %H%;
const cv = document.getElementById('art'); const g = cv.getContext('2d'); g.scale(2, 2);
function grid(f) { // node lookup on the regular grid
  const xs = f.x, ys = f.y; let h = Infinity;
  for (let i = 1; i < xs.length; i++) { const d = Math.abs(xs[i] - xs[i-1]); if (d > 1e-6 && d < h) h = d; }
  const map = new Map(); let vmax = 0;
  for (let i = 0; i < xs.length; i++) { map.set(Math.round(xs[i]/h) + ',' + Math.round(ys[i]/h), f.v[i]); vmax = Math.max(vmax, Math.abs(f.v[i])); }
  return { h, map, vmax, r: f.r };
}
function disc(f, cx, cy, rx, ry, base, opts) {
  const G = grid(f); const W = Math.ceil(rx * 2), Hh = Math.ceil(ry * 2);
  const img = g.createImageData(W * 2, Hh * 2);
  for (let py = 0; py < Hh * 2; py++) for (let px = 0; px < W * 2; px++) {
    const u = (px / 2 - rx) / rx, v = (py / 2 - ry) / ry; const rr = u*u + v*v;
    if (rr > 1) continue;
    const X = u * G.r / G.h, Y = -v * G.r / G.h;
    const i0 = Math.floor(X), j0 = Math.floor(Y), fx = X - i0, fy = Y - j0;
    const at = (i, j) => G.map.get(i + ',' + j) ?? 0;
    let val = (1-fx)*(1-fy)*at(i0,j0) + fx*(1-fy)*at(i0+1,j0) + (1-fx)*fy*at(i0,j0+1) + fx*fy*at(i0+1,j0+1);
    if (opts.hole && rr < opts.hole) continue;
    const s = Math.sign(val) * Math.pow(Math.min(1, Math.abs(val) / (G.vmax * opts.sat)), 0.55);
    // base material with radial shading and lathing rings
    const shade = 0.62 + 0.38 * (1 - rr) + 0.10 * (-u * 0.6 - v * 0.8);
    const ring = opts.lathe ? 0.06 * Math.sin(Math.sqrt(rr) * 260) : 0;
    let [r0, g0, b0] = base.map(c => c * (shade + ring));
    const a = Math.abs(s);
    const [cr, cg, cb] = s > 0 ? [255, 176, 72] : [80, 220, 200];
    const k = Math.min(1, a * opts.gain);
    const o = 4 * (py * W * 2 + px);
    img.data[o] = r0 * (1 - k) + cr * k; img.data[o+1] = g0 * (1 - k) + cg * k; img.data[o+2] = b0 * (1 - k) + cb * k;
    const edge = Math.min(1, (1 - Math.sqrt(rr)) * 60);
    img.data[o+3] = 255 * edge;
  }
  const off = document.createElement('canvas'); off.width = W * 2; off.height = Hh * 2;
  off.getContext('2d').putImageData(img, 0, 0);
  g.save(); g.shadowColor = 'rgba(0,0,0,0.6)'; g.shadowBlur = 40; g.shadowOffsetY = 18;
  g.drawImage(off, cx - rx, cy - ry, W, Hh); g.restore();
  // specular sheen
  const sh = g.createRadialGradient(cx - rx * 0.35, cy - ry * 0.45, 0, cx - rx * 0.35, cy - ry * 0.45, rx * 0.9);
  sh.addColorStop(0, 'rgba(255,240,210,' + opts.sheen + ')'); sh.addColorStop(1, 'rgba(255,240,210,0)');
  g.save(); g.beginPath(); g.ellipse(cx, cy, rx, ry, 0, 0, Math.PI * 2); g.clip(); g.fillStyle = sh; g.fillRect(cx - rx, cy - ry, 2 * rx, 2 * ry); g.restore();
}
// faint measurement grid
g.strokeStyle = 'rgba(241,234,219,0.035)'; g.lineWidth = 1;
for (let x = 0; x < 1200; x += 40) { g.beginPath(); g.moveTo(x + 0.5, 0); g.lineTo(x + 0.5, H); g.stroke(); }
for (let y = 0; y < H; y += 40) { g.beginPath(); g.moveTo(0, y + 0.5); g.lineTo(1200, y + 0.5); g.stroke(); }
// crash cymbal on its stand, then the snare in front
const CY = %CRASHY%;
g.strokeStyle = 'rgba(200,200,210,0.35)'; g.lineWidth = 3; g.beginPath(); g.moveTo(905, CY + 10); g.lineTo(905, H); g.stroke();
disc(F.crash, 905, CY, 255, 112, [176, 128, 62], { sat: 0.45, gain: 0.9, lathe: true, sheen: 0.22, hole: 0.0009 });
// bell
const bell = g.createRadialGradient(895, CY - 6, 2, 905, CY, 40); bell.addColorStop(0, 'rgba(255,226,160,0.55)'); bell.addColorStop(1, 'rgba(160,110,50,0)');
g.fillStyle = bell; g.beginPath(); g.ellipse(905, CY, 44, 19, 0, 0, Math.PI * 2); g.fill();
g.fillStyle = '#2a2a30'; g.beginPath(); g.ellipse(905, CY, 6, 3, 0, 0, Math.PI * 2); g.fill();
// snare: shell, hoop, head
const SX = 735, SY = %SNAREY%, SR = 150, SRY = 58, DEPTH = 62;
const shell = g.createLinearGradient(SX - SR, 0, SX + SR, 0);
shell.addColorStop(0, '#3a1c10'); shell.addColorStop(0.35, '#8a4422'); shell.addColorStop(0.55, '#b5652e'); shell.addColorStop(1, '#2e160c');
g.fillStyle = shell; g.beginPath(); g.ellipse(SX, SY + DEPTH, SR, SRY, 0, 0, Math.PI); g.lineTo(SX - SR, SY); g.ellipse(SX, SY, SR, SRY, 0, Math.PI, 0, true); g.closePath(); g.fill();
for (let i = 0; i < 8; i++) { const ang = Math.PI * (0.08 + 0.84 * i / 7); const lx = SX - SR * Math.cos(ang); g.fillStyle = 'rgba(225,228,235,0.85)'; g.fillRect(lx - 3, SY + 14 + Math.sin(ang) * SRY * 0.9, 6, 30); }
disc(F.snare, SX, SY, SR - 6, SRY - 3, [228, 224, 214], { sat: 0.6, gain: 0.8, lathe: false, sheen: 0.18 });
g.strokeStyle = '#e4e7ee'; g.lineWidth = 5; g.beginPath(); g.ellipse(SX, SY, SR - 2, SRY - 1, 0, 0, Math.PI * 2); g.stroke();
g.strokeStyle = 'rgba(255,255,255,0.6)'; g.lineWidth = 1.2; g.beginPath(); g.ellipse(SX, SY - 1.5, SR - 2, SRY - 1, 0, Math.PI * 1.05, Math.PI * 1.95); g.stroke();
// the groove waveform along the bottom
const WY = H - 34, WH = 22;
g.fillStyle = 'rgba(217,169,79,0.32)';
for (let i = 0; i < ENV.length; i++) { const x = 640 + i * (520 / ENV.length); const h = Math.max(1, ENV[i] * WH); g.fillRect(x, WY - h, 520 / ENV.length * 0.6, 2 * h); }
</script></body></html>"""
for H, name in ((630, "og.html"), (600, "twitter.html")):
    s = html.replace("%FR%", str(F)).replace("%IN%", str(I)).replace("%JB%", str(J))
    s = s.replace("%FIELDS%", json.dumps(fields)).replace("%ENV%", json.dumps(env))
    s = s.replace("%H2%", str(2 * H)).replace("%H%", str(H))
    s = s.replace("%CRASHY%", str(205 if H == 630 else 196)).replace("%SNAREY%", str(400 if H == 630 else 382))
    s = s.replace("%CHIPB%", str(52 if H == 630 else 46)).replace("%URLB%", str(64 if H == 630 else 58))
    s = s.replace("%CAPX%", "760").replace("%CAPY%", str(70 if H == 630 else 64))
    (S / name).write_text(s)
print("ok")
