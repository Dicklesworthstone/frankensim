// Instrument Explorer: data loading, DOM helpers and the hash router shared by the views.
import { stopAll } from './audio.js';

export const DATA = 'data/';

const cache = new Map();

/** Fetch a JSON file of the built site once; later calls share the promise. */
export function getJSON(path) {
  if (!cache.has(path)) {
    cache.set(path, fetch(DATA + path).then(r => {
      if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
      return r.json();
    }));
  }
  return cache.get(path);
}

function append(node, kids) {
  for (const k of kids) {
    if (k === null || k === undefined || k === false) continue;
    if (Array.isArray(k)) append(node, k);
    else node.append(k instanceof Node ? k : String(k));
  }
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
  append(node, kids);
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

/** Routes: #/  #/f/<family>[/<tab>[/<item>]]  #/process  #/about */
export function parseRoute() {
  const parts = location.hash.replace(/^#\/?/, '').split('/').map(decodeURIComponent);
  if (parts[0] === 'f' && parts[1]) return { view: 'family', id: parts[1], tab: parts[2] || 'listen', item: parts[3] || null };
  if (parts[0] === 'process') return { view: 'process' };
  if (parts[0] === 'about') return { view: 'about' };
  return { view: 'home' };
}

export function link(...parts) {
  return `#/${parts.filter(p => p !== null && p !== undefined).map(encodeURIComponent).join('/')}`;
}

const VIEWS = {
  home: () => import('./views/home.js'),
  family: () => import('./views/family.js'),
  process: () => import('./views/process.js'),
  about: () => import('./views/about.js'),
};

let lastView = null;
let lastFamily = null;

async function route() {
  const r = parseRoute();
  const main = document.getElementById('main');
  for (const a of document.querySelectorAll('header nav a')) {
    a.classList.toggle('active', a.dataset.view === r.view || (r.view === 'family' && a.dataset.view === 'home'));
  }
  // A tab change inside one family keeps the page (and its header) and only swaps the tab.
  const sameFamily = r.view === 'family' && lastView === 'family' && lastFamily === r.id;
  if (!sameFamily) {
    stopAll();
    main.replaceChildren(h('p.loading', null, 'Loading…'));
    window.scrollTo(0, 0);
  }
  lastView = r.view;
  lastFamily = r.view === 'family' ? r.id : null;
  try {
    const mod = await VIEWS[r.view]();
    await mod.render(main, r, sameFamily);
  } catch (err) {
    console.error(err);
    main.replaceChildren(h('div.error', null, h('h2', null, 'Could not show this page'), h('pre', null, String(err && err.stack || err))));
  }
}

async function footer() {
  try {
    const idx = await getJSON('index.json');
    document.getElementById('built').textContent =
      `Built ${fmt.datetime(idx.built_utc)} from ${idx.ref} ${idx.ref_sha} and ${idx.mirror}.`;
  } catch (err) {
    document.getElementById('built').textContent =
      'No built data found. Run tools/sync.sh, then tools/build.py, then tools/serve.py (see the README).';
  }
}

window.addEventListener('hashchange', route);
route();
footer();
