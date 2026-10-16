// Loopback-only companion. No packages, shell, network fetches, or general file API.
import http from 'node:http';
import {readFile, mkdtemp, rm, access, stat} from 'node:fs/promises';
import {constants} from 'node:fs';
import {spawn} from 'node:child_process';
import {randomBytes} from 'node:crypto';
import {tmpdir} from 'node:os';
import {dirname, join, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';
import {MAX_BYTES, parseNative, validateInputs} from './protocol.mjs';
const here = dirname(fileURLToPath(import.meta.url));
export const BUILD = 'cargo build --release -p fs-couple --example heated_section_sound';
const assets = {'/': 'index.html', '/index.html': 'index.html', '/styles.css': 'styles.css',
  '/app.mjs': 'app.mjs', '/visuals.mjs': 'visuals.mjs', '/protocol.mjs': 'protocol.mjs',
  '/mesh-data.mjs': 'mesh-data.mjs', '/import-vtu.mjs': 'import-vtu.mjs',
  '/mesh-ui.mjs': 'mesh-ui.mjs', '/mesh-view.mjs': 'mesh-view.mjs', '/mesh.css': 'mesh.css', '/mesh-raster.mjs': 'mesh-raster.mjs'};
const types = {html: 'text/html', css: 'text/css', mjs: 'text/javascript'};
async function ready(engine) { try { await access(engine, constants.X_OK); return (await stat(engine)).isFile(); } catch { return false; } }
async function body(req) {
  let size = 0; const chunks = [];
  for await (const chunk of req) {
    size += chunk.length;
    if (size > 4096) throw Object.assign(new Error('Request exceeds 4 KiB.'), {status: 413});
    chunks.push(chunk);
  }
  return JSON.parse(Buffer.concat(chunks).toString('utf8'));
}
export function createObservatory({engine, packs = null, deadline = 60000} = {}) {
  engine = resolve(engine || join(here, '../../target/release/examples/heated_section_sound' + (process.platform === 'win32' ? '.exe' : '')));
  const token = randomBytes(32).toString('hex'); let active = false; let cancelActive = () => {};
  const server = http.createServer(async (req, res) => {
    const port = server.address().port;
    const hosts = [`127.0.0.1:${port}`, `localhost:${port}`];
    const send = (status, value) => { if (!res.destroyed && !res.writableEnded) {
      res.writeHead(status, {'Content-Type': 'application/json; charset=utf-8'}); res.end(JSON.stringify(value));
    }};
    res.setHeader('Cache-Control', 'no-store'); res.setHeader('X-Content-Type-Options', 'nosniff');
    res.setHeader('Referrer-Policy', 'no-referrer');
    res.setHeader('Content-Security-Policy', "default-src 'none'; script-src 'self'; style-src 'self'; img-src 'self' blob: data:; connect-src 'self'; media-src 'self' blob:; base-uri 'none'; frame-ancestors 'none'; form-action 'none'");
    if (!hosts.includes(req.headers.host)) return send(403, {error: 'Only loopback hosts are accepted.'});
    if (req.headers.origin && req.headers.origin !== `http://${req.headers.host}`) return send(403, {error: 'Cross-origin requests are refused.'});
    if (req.headers['sec-fetch-site'] && !['same-origin', 'none'].includes(req.headers['sec-fetch-site'])) return send(403, {error: 'Cross-site requests are refused.'});
    const path = req.url;
    try {
      if (req.method === 'GET' && path === '/api/status') return send(200, {ready: await ready(engine),
        busy: active, material: packs ? 'packs' : 'synthetic', token, build: BUILD});
      if (req.method === 'GET' && Object.hasOwn(assets, path)) {
        const file = assets[path]; const bytes = await readFile(join(here, file));
        res.writeHead(200, {'Content-Type': `${types[file.split('.').at(-1)]}; charset=utf-8`}); return res.end(bytes);
      }
      if (path !== '/api/run' || req.method !== 'POST') return send(404, {error: 'Not found.'});
      if (req.headers['x-observatory-token'] !== token) return send(403, {error: 'Missing local session token. Reload the page.'});
      if (req.headers['content-type']?.split(';')[0] !== 'application/json') return send(415, {error: 'Use application/json.'});
      if (active) return send(409, {error: 'Another experiment is still running. Cancel it or wait for completion.'});
      const input = validateInputs(await body(req));
      if (packs && input.resistivity !== 1) throw new Error('Pack materials cannot be rescaled. Set resistivity scale to 1.');
      if (!await ready(engine)) return send(503, {error: `Native engine is unavailable. Run: ${BUILD}`});
      if (active) return send(409, {error: 'Another experiment started while this request was being admitted.'});
      active = true;
      let dir;
      try {
        dir = await mkdtemp(join(tmpdir(), 'frankensim-observatory-'));
        if (res.destroyed) return;
        const output = join(dir, 'result');
        const args = packs ? ['--packs', ...packs, String(input.voltage), String(input.duration), '1', output]
          : ['--demo', String(input.resistivity), String(input.voltage), String(input.duration), '1', output];
        if (input.support === 'strain') args.push('--fixed-strain', String(input.strain), String(input.limit));
        const log = await new Promise((accept, reject) => {
          const child = spawn(engine, args, {cwd: dir, shell: false, windowsHide: true, stdio: ['ignore', 'pipe', 'pipe']});
          const chunks = []; let size = 0, reason = null, killTimer;
          const stop = message => {
            if (reason) return;
            reason = new Error(message); child.kill('SIGTERM');
            killTimer = setTimeout(() => child.kill('SIGKILL'), 1000); killTimer.unref();
          };
          cancelActive = () => stop('Experiment cancelled.');
          const disconnected = () => { if (!res.writableEnded) stop('Browser disconnected; experiment cancelled.'); };
          res.once('close', disconnected);
          const timer = setTimeout(() => stop('Experiment exceeded its 60-second execution budget.'), deadline);
          const read = chunk => { size += chunk.length; if (size > 256 * 1024) stop('Native diagnostics exceeded 256 KiB.'); else chunks.push(chunk); };
          child.stdout.on('data', read); child.stderr.on('data', read);
          child.once('error', error => { reason = error; });
          child.once('close', code => {
            clearTimeout(timer); clearTimeout(killTimer); res.off('close', disconnected); cancelActive = () => {};
            const log = Buffer.concat(chunks).toString('utf8');
            if (reason) reject(reason);
            else if (code !== 0) reject(new Error(`Native solver refused the experiment. ${log.slice(-3000)}`));
            else accept(log);
          });
        });
        if (res.destroyed) return;
        if ((await stat(output + '.csv')).size > MAX_BYTES) throw new Error('Native trace exceeds 8 MiB.');
        const snapshot = parseNative(await readFile(output + '.csv', 'utf8'), log, input, Boolean(packs));
        send(200, {snapshot});
      } finally {
        try { if (dir) await rm(dir, {recursive: true, force: true}); } // Only our uniquely created scratch directory.
        finally { active = false; }
      }
    } catch (error) { send(error.status || 400, {error: error.message || 'Experiment failed.'}); }
  });
  server.headersTimeout = 10000; server.requestTimeout = 15000;
  return {server, shutdown: () => {cancelActive(); server.close(); server.closeIdleConnections();}};
}
function cli() {
  const args = process.argv.slice(2); const options = {}; let port = 7070, open = false;
  for (let i = 0; i < args.length; i++) {
    const a = args[i];
    if (a === '--help') { console.log('node apps/observatory/server.mjs [--port 7070] [--open] [--engine PATH] [--core-pack PATH --sleeve-pack PATH]'); return; }
    if (a === '--open') {open = true; continue;}
    if (!['--port', '--engine', '--core-pack', '--sleeve-pack'].includes(a) || !args[i + 1]) throw new Error(`Unknown or incomplete option: ${a}`);
    const value = args[++i];
    if (a === '--port') port = Number(value); else options[a.slice(2)] = resolve(value);
  }
  if (!Number.isInteger(port) || port < 1024 || port > 65535) throw new Error('Port must be 1024–65535.');
  if (Boolean(options['core-pack']) !== Boolean(options['sleeve-pack'])) throw new Error('Supply both material packs.');
  const app = createObservatory({engine: options.engine,
    packs: options['core-pack'] ? [options['core-pack'], options['sleeve-pack']] : null});
  app.server.on('error', error => { console.error(`Observatory: ${error.message}`); process.exitCode = 1; });
  app.server.listen(port, '127.0.0.1', () => {
    const url = `http://127.0.0.1:${port}`;
    console.log(`\n  FRANKENSIM / OBSERVATORY\n  ${url}\n\n  Local only. Ctrl+C to stop.\n  Native engine build: ${BUILD}\n`);
    if (open) {
      // URL is generated above, never supplied by a browser request.
      const command = process.platform === 'darwin' ? ['open', [url]] : process.platform === 'win32'
        ? ['rundll32.exe', ['url.dll,FileProtocolHandler', url]] : ['xdg-open', [url]];
      spawn(command[0], command[1], {stdio: 'ignore', shell: false}).on('error', () => console.log(`Open ${url} in a browser.`));
    }
  });
  process.once('SIGINT', app.shutdown); process.once('SIGTERM', app.shutdown);
}
if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  try { cli(); } catch (e) {console.error(e.message); process.exitCode = 1;}
}
