// Instrument Explorer: data loading, DOM helpers and the hash router shared by the views.
import { stopAll, mountNowPlaying } from './audio.js';

export const DATA = 'data/';

const cache = new Map();

/** Fetch a JSON file of the built site once; later calls share the promise. A failure is not kept: one dropped
 *  connection would otherwise break every later page that needs the file until a reload. A network error, a cut-off
 *  body or a server error is retried once after a short pause before the page shows its error. */
export function getJSON(path) {
  if (!cache.has(path)) {
    const load = () => fetch(DATA + path).then(r => {
      if (!r.ok) throw Object.assign(new Error(`${path}: HTTP ${r.status}`), { transient: r.status >= 500 });
      return r.json();
    });
    const p = load().catch(e => e instanceof TypeError || e instanceof SyntaxError || e.transient
      ? new Promise(res => setTimeout(res, 600)).then(load)
      : Promise.reject(e));
    cache.set(path, p);
    p.catch(() => { if (cache.get(path) === p) cache.delete(path); });
  }
  return cache.get(path);
}

function append(node, kids, html = true) {
  for (const k of kids) {
    if (k === null || k === undefined || k === false) continue;
    if (Array.isArray(k)) append(node, k, html);
    else if (k instanceof Node) node.append(k);
    else if (html) node.append(...keepRanges(String(k)));
    else node.append(String(k));
  }
}

// A line may break after an en dash, which strands the "9" of "bars 5–9" at the start of the next line. In HTML text an
// en-dashed token goes into a nowrap span so it moves as one word; spaced dashes (" – ") stay ordinary breaks.
function keepRanges(t) {
  if (!t.includes('–')) return [t];
  const out = [];
  let at = 0;
  for (const m of t.matchAll(/\S+–\S+/g)) {
    if (m.index > at) out.push(t.slice(at, m.index));
    const span = document.createElement('span');
    span.className = 'nw';
    span.textContent = m[0];
    out.push(span);
    at = m.index + m[0].length;
  }
  if (at < t.length) out.push(t.slice(at));
  return out;
}

function setAttrs(node, attrs) {
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k.startsWith('on') && typeof v === 'function') node.addEventListener(k.slice(2), v);
    else if (k === 'style' && typeof v === 'object') Object.assign(node.style, v);
    else node.setAttribute(k, v === true ? '' : v);
  }
}

/** h('div.card.wide', {title: 'x'}, child, [children], 'text') */
export function h(tag, attrs, ...kids) {
  const [name, ...classes] = tag.split('.');
  const node = document.createElement(name || 'div');
  if (classes.length) node.className = classes.join(' ');
  setAttrs(node, attrs);
  append(node, kids);
  return node;
}

const SVG_NS = 'http://www.w3.org/2000/svg';

/** The same for SVG elements. */
export function s(tag, attrs, ...kids) {
  const [name, ...classes] = tag.split('.');
  const node = document.createElementNS(SVG_NS, name);
  if (classes.length) node.setAttribute('class', classes.join(' '));
  setAttrs(node, attrs);
  append(node, kids, false);
  return node;
}

export const fmt = {
  num(v, d = 1) {
    if (v === null || v === undefined || Number.isNaN(v)) return '–';
    const a = Math.abs(v);
    if (a !== 0 && (a < 10 ** -d || a >= 1e5)) return v.toExponential(1);
    return v.toFixed(d);
  },
  pct(f, d = 0) {
    return f === null || f === undefined ? '–' : `${(100 * f).toFixed(d)} %`;
  },
  date(iso) {
    if (!iso) return '–';
    const d = new Date(iso.length <= 10 ? `${iso}T12:00:00Z` : iso.replace(' ', 'T').replace(/Z?$/, 'Z'));
    if (Number.isNaN(d.getTime())) return iso;
    return d.toLocaleDateString(undefined, { month: 'short', day: 'numeric' });
  },
  datetime(iso) {
    if (!iso) return '–';
    const d = new Date(iso.replace(' ', 'T').replace(/Z?$/, 'Z'));
    if (Number.isNaN(d.getTime())) return iso;
    return d.toLocaleString(undefined, { month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' });
  },
  time(t) {
    if (!Number.isFinite(t)) return '0:00.0';
    const m = Math.floor(t / 60);
    return `${m}:${(t - 60 * m).toFixed(1).padStart(4, '0')}`;
  },
  plural(n, one, many = `${one}s`) {
    return `${n} ${n === 1 ? one : many}`;
  },
};

export const BADGES = {
  close: 'Inside the real spread',
  mixed: 'Partly inside the real spread',
  off: 'Outside the real spread',
  listening: 'Rendered; no recording to grade against',
  ungraded: 'Rendered, not graded',
  none: 'No render yet',
};

export function badge(b, label) {
  return h(`span.badge.b-${b}`, { title: label || BADGES[b] }, h('span.dot'), label || BADGES[b]);
}

/** A small "inside / outside" proportion bar. */
export function splitBar(inside, outside, title) {
  const total = inside + outside;
  if (!total) return h('div.split.empty', { title: title || 'no graded measures' });
  return h('div.split', { title: title || `${inside} inside, ${outside} outside` },
    h('span.in', { style: { flex: String(inside) } }),
    h('span.out', { style: { flex: String(outside) } }));
}

export function chip(text, cls = '') {
  return h(`span.chip${cls ? `.${cls}` : ''}`, null, text);
}

export function empty(text) {
  return h('p.empty', null, text);
}

/** On a phone the rows of a `stack` table become cards, one cell to a line (css/app.css), so no column sits past the
 *  right edge. Each cell is named by its column's heading, or by `labels[i]` where that heading is too long or should
 *  not show ('' for none). Returns the table. */
export function stack(table, labels = null) {
  const heads = [...table.querySelectorAll('thead th')].map((th, i) => (labels && labels[i] !== undefined ? labels[i] : th.textContent));
  table.classList.add('stack');
  for (const tr of table.querySelectorAll('tbody tr')) [...tr.children].forEach((td, i) => { if (heads[i]) td.dataset.label = heads[i]; });
  return table;
}

/** Routes: #/  #/instruments  #/f/<family>[/<section>[/<item>]]  #/p/<phrase>  #/method  #/decisions  #/process  #/about */
export function parseRoute() {
  const decode = x => {
    try { return decodeURIComponent(x); } catch (e) { return x; }
  };
  const parts = location.hash.replace(/^#\/?/, '').split('/').map(decode);
  if (parts[0] === 'p' && parts[1]) return { view: 'phrase', id: parts[1] };
  // the old per-family phrase links (#/f/<family>/hear/<phrase>) now open the phrase's own page
  if (parts[0] === 'f' && parts[2] === 'hear' && parts[3]) return { view: 'phrase', id: parts[3] };
  if (parts[0] === 'f' && parts[1]) return { view: 'family', id: parts[1], tab: parts[2] || null, item: parts[3] || null };
  if (['instruments', 'method', 'decisions', 'process', 'about'].includes(parts[0])) return { view: parts[0] };
  // anything else (an old or mistyped link) shows the start page, saying so
  return { view: 'home', unknown: parts[0] ? location.hash : null };
}

/** How a programmatic scroll over `distance` pixels should move: smoothly only for a short hop. A smooth flight over a
 *  long page (a tab 40 000 px down takes Chromium about 1.5 s) ignores the reader's own wheel and touch until it lands,
 *  so anything further than two screens jumps at once, and so does everything under reduced motion. */
export function scrollMotion(distance) {
  if (matchMedia('(prefers-reduced-motion: reduce)').matches) return 'auto';
  return Math.abs(distance) > 2 * innerHeight ? 'auto' : 'smooth';
}

export function link(...parts) {
  return `#/${parts.filter(p => p !== null && p !== undefined).map(encodeURIComponent).join('/')}`;
}

const VIEWS = {
  home: () => import('./views/home.js'),
  instruments: () => import('./views/instruments.js'),
  family: () => import('./views/family.js'),
  phrase: () => import('./views/analysis.js'),
  method: () => import('./views/method.js'),
  decisions: () => import('./views/decisions.js'),
  process: () => import('./views/process.js'),
  about: () => import('./views/about.js'),
};

let lastView = null;
let lastFamily = null;
let navigation = 0;

// Back, forward and reload return to where the reader was. Each history entry carries a key, the scroll position is
// remembered per key, and it is kept in sessionStorage over a reload. Tracking pauses while a route swaps the page,
// because the shrinking page scrolls by itself and would overwrite the position being returned to.
if ('scrollRestoration' in history) history.scrollRestoration = 'manual';
let positions = {};
try { positions = JSON.parse(sessionStorage.getItem('ie-scroll')) || {}; } catch (e) { positions = {}; }
let tracking = false;
addEventListener('scroll', () => {
  const k = tracking && history.state && history.state.k;
  if (k) positions[k] = Math.round(scrollY);
}, { passive: true });
addEventListener('pagehide', () => {
  try { sessionStorage.setItem('ie-scroll', JSON.stringify(positions)); } catch (e) { /* private mode: memory only */ }
});

/** A phone's header nav scrolls sideways: bring the current page's link into view (About sits past the edge). The
 * nav fades out over its last 36 px (app.css), so a link under the fade counts as hidden. */
function revealNavLink(a) {
  const nav = a.parentElement;
  const ar = a.getBoundingClientRect();
  const nr = nav.getBoundingClientRect();
  const fade = nav.scrollWidth > nav.clientWidth ? 36 : 0;
  if (ar.left < nr.left || ar.right > nr.right - fade) nav.scrollLeft += ar.left - nr.left - (nr.width - ar.width) / 2;
}
// the web fonts widen the links after the first route has placed them
document.fonts.ready.then(() => { const a = document.querySelector('header nav a.active'); if (a) revealNavLink(a); });

/** `fresh` treats the current entry as newly arrived at (a click on a link to the page already shown). */
async function route({ fresh = false } = {}) {
  const seq = ++navigation;
  const r = parseRoute();
  const state = history.state;
  const returning = !fresh && Boolean(state && state.k);
  if (!state || !state.k) history.replaceState({ ...(state || {}), k: Math.random().toString(36).slice(2, 10) }, '');
  if (returning && positions[state.k] !== undefined) r.restoreY = positions[state.k];
  tracking = false;
  const main = document.getElementById('main');
  for (const a of document.querySelectorAll('header nav a')) {
    const on = a.dataset.view === r.view || (r.view === 'family' && a.dataset.view === 'instruments')
      || (r.view === 'phrase' && a.dataset.view === 'home');
    a.classList.toggle('active', on);
    if (on) {
      a.setAttribute('aria-current', 'page');
      revealNavLink(a);
    } else a.removeAttribute('aria-current');
  }
  // A tab change inside one family keeps the page (and its header) and only swaps the tab.
  const sameFamily = r.view === 'family' && lastView === 'family' && lastFamily === r.id;
  if (!sameFamily) {
    stopAll();
    main.replaceChildren(h('p.loading', null, 'Loading…'));
    window.scrollTo(0, 0);
    if (lastView !== null) main.focus({ preventScroll: true });
  }
  lastView = r.view;
  lastFamily = r.view === 'family' ? r.id : null;
  try {
    const mod = await VIEWS[r.view]();
    // A new page is built off-screen and goes in only if no later navigation started while its data loaded.
    const host = sameFamily ? main : document.createElement('div');
    await mod.render(host, r, sameFamily);
    if (seq !== navigation) return;
    if (host !== main) main.replaceChildren(...host.childNodes);
    document.title = `${mod.title ? `${mod.title(r)} · ` : ''}Instrument tests · FrankenSim`;
    if (r.restoreY !== undefined) window.scrollTo(0, r.restoreY);
    tracking = true;
  } catch (err) {
    console.error(err);
    if (seq !== navigation) return;
    tracking = true;
    document.title = 'Could not show this page · Instrument tests · FrankenSim';
    main.replaceChildren(h('div.error', null, h('h2', null, 'Could not show this page'),
      h('p', null, String((err && err.message) || err)), h('p', null, h('a', { href: '#/' }, 'Back to the start'))));
  }
}

// Keyboard focus the browser scrolled to can still sit under a bar stuck to the top of the screen (the header, a family's
// section tabs, an analysis page's transport): move it out from under them.
document.addEventListener('focusin', e => {
  const el = e.target;
  if (!(el instanceof Element) || !el.matches(':focus-visible') || el.closest('header.site, nav.sectionbar, .anatransport')) return;
  let bottom = 0;
  for (const bar of document.querySelectorAll('header.site, nav.sectionbar, .anatransport')) {
    const r = bar.getBoundingClientRect();
    if (getComputedStyle(bar).position === 'sticky' && r.top <= 1 + bottom && r.bottom > bottom) bottom = r.bottom;
  }
  const r = el.getBoundingClientRect();
  if (bottom && r.top < bottom + 6 && r.bottom > 0) window.scrollBy(0, r.top - bottom - 12);
});

// The skip link moves focus; following its #main would be read as a route and land on the home page.
document.querySelector('a.skip').addEventListener('click', e => {
  e.preventDefault();
  document.getElementById('main').focus();
});

async function footer() {
  try {
    const idx = await getJSON('index.json');
    document.getElementById('built').textContent =
      `Data built ${fmt.datetime(idx.built_utc)} from FrankenSim ${idx.ref_sha}. Every simulated sound is computed from physical descriptions; recordings are credited in each instrument's References section.`;
  } catch (err) {
    document.getElementById('built').textContent =
      'No built data found. Run tools/sync.sh, then tools/build.py, then tools/serve.py (see the README).';
  }
}

window.addEventListener('hashchange', route);

// A link to the page already shown changes nothing in the address bar, so no route would run: on a family page it goes
// to the section it names again, anywhere else back to the top.
document.addEventListener('click', e => {
  const a = e.target.closest && e.target.closest('a[href^="#/"]');
  if (!a || e.defaultPrevented || a.getAttribute('href') !== location.hash || e.button || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
  e.preventDefault();
  const r = parseRoute();
  if (r.view === 'family' && r.tab) route({ fresh: true });
  else window.scrollTo({ top: 0, behavior: scrollMotion(scrollY) });
});
mountNowPlaying(document.getElementById('nowplaying'));
route();
footer();
