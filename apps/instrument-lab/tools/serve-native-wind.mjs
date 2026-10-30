// Optional loopback companion for the original native physical performers.
import http from 'node:http';
import {access, lstat, mkdir, mkdtemp, open, realpath, rm, stat, writeFile} from 'node:fs/promises';
import {constants, createReadStream} from 'node:fs';
import {spawn} from 'node:child_process';
import {randomBytes} from 'node:crypto';
import {tmpdir} from 'node:os';
import {dirname, extname, isAbsolute, join, relative, resolve, sep} from 'node:path';
import {fileURLToPath} from 'node:url';
import {DEFAULT_TIMEOUT_MS, MAX_LOG_BYTES, MAX_OUTPUT_BYTES, MAX_SOURCE_BYTES,
  MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES, MAX_REQUEST_BYTES,
  MAX_HISTORY_BYTES, FINITE_GAS_STRING_RESULT_FORMAT, FINITE_GAS_STRING_SCHEMAS,
  NativeWindError, SUPPORTED_SCHEMAS, STRING_SCHEMAS, RESOLVED_SCHEMAS, validateInput, validateWave} from './native-wind-protocol.mjs';

import {nativePressureReceipt} from './native-pressure-receipt.mjs';
import {PRESSURE_RECEIPT_HEADER} from '../src/native-wave.mjs';

const here = dirname(fileURLToPath(import.meta.url));
const types = {'.html':'text/html; charset=utf-8', '.js':'text/javascript; charset=utf-8',
  '.mjs':'text/javascript; charset=utf-8', '.css':'text/css; charset=utf-8', '.json':'application/json',
  '.wasm':'application/wasm', '.svg':'image/svg+xml', '.png':'image/png', '.webp':'image/webp',
  '.jpg':'image/jpeg', '.jpeg':'image/jpeg', '.ico':'image/x-icon', '.wav':'audio/wav',
  '.flac':'audio/flac', '.m4a':'audio/mp4', '.mp4':'video/mp4', '.woff':'font/woff', '.woff2':'font/woff2',
  '.mid':'audio/midi', '.midi':'audio/midi', '.csv':'text/csv; charset=utf-8',
  '.tsv':'text/tab-separated-values; charset=utf-8', '.txt':'text/plain; charset=utf-8',
  '.performance':'text/plain; charset=utf-8'};
const inside = (root, file) => { const r = relative(root, file); return r !== '..' && !r.startsWith('..' + sep) && !isAbsolute(r); };
async function available(engine) {
  if (!engine) return false;
  try { await access(engine, constants.X_OK); return (await stat(engine)).isFile(); } catch { return false; }
}
function sourceBody(req, signal, bundled) {
  const limit = bundled ? MAX_REQUEST_BYTES : MAX_SOURCE_BYTES;
  return new Promise((accept, reject) => {
    let size = 0; const chunks = [];
    const cleanup = () => { req.off('data', data); req.off('end', end); req.off('error', fail);
      req.off('aborted', aborted); signal.removeEventListener('abort', aborted); };
    const fail = error => { cleanup(); req.resume(); reject(error); };
    const aborted = () => fail(signal.reason || new NativeWindError(400, 'Performance upload disconnected.'));
    const data = chunk => { size += chunk.length;
      if (size > limit) fail(new NativeWindError(413, 'Experiment exceeds the complete request limit.'));
      else chunks.push(chunk); };
    const end = () => { cleanup(); try { accept(validateInput(Buffer.concat(chunks), bundled)); } catch (error) { reject(error); } };
    req.on('data', data); req.once('end', end); req.once('error', fail); req.once('aborted', aborted);
    signal.addEventListener('abort', aborted, {once:true});
    if (signal.aborted || req.aborted) aborted();
  });
}
async function serveAsset(req, res, dist, pathname) {
  let decoded;
  try { decoded = decodeURIComponent(pathname); } catch { throw new NativeWindError(400, 'Invalid asset path.'); }
  if (decoded.includes('\0') || decoded.split('/').some(p => p.startsWith('.'))) return false;
  const file = resolve(dist, '.' + (decoded === '/' ? '/index.html' : decoded));
  if (!inside(dist, file)) return false;
  let actual, info;
  try { const root = await realpath(dist); actual = await realpath(file);
    if (!inside(root, actual)) return false;
    info = await stat(actual); if (!info.isFile()) return false;
  } catch { return false; }
  const headers = {'Content-Type':types[extname(file).toLowerCase()] || 'application/octet-stream',
    'Content-Length':info.size, 'Accept-Ranges':'bytes'};
  let start = 0, end = info.size - 1, status = 200;
  // Range applies to GET only. Unsupported/malformed or multipart ranges are
  // ignored; without representation validators, If-Range also receives the full file.
  const range = req.method === 'GET' && !req.headers['if-range']
    ? /^bytes=(\d*)-(\d*)$/.exec(req.headers.range || '') : null;
  if (range && (range[1] || range[2])) {
    if (range[1]) {
      start = Number(range[1]); end = range[2] ? Math.min(Number(range[2]), end) : end;
    } else start = Math.max(0, info.size - Number(range[2]));
    if (start >= info.size || start > end) {
      res.writeHead(416, {...headers, 'Content-Range':`bytes */${info.size}`, 'Content-Length':0});
      res.end(); return true;
    }
    status = 206; headers['Content-Range'] = `bytes ${start}-${end}/${info.size}`;
    headers['Content-Length'] = end - start + 1;
  }
  res.writeHead(status, headers);
  if (req.method === 'HEAD' || info.size === 0) res.end();
  else { const stream = createReadStream(actual, {start, end}); stream.on('error', () => res.destroy());
    res.once('close', () => stream.destroy()); stream.pipe(res); }
  return true;
}
function runNative(engine, input, output, dir, signal, outputLimit, experiment, historyLimit) {
    const {attachments, decimate, gasSources, stringProfiles, coatingStrainGuard, finiteGasString} = experiment;
    const sourceArgs = attachments.flatMap(file => ['--attached-source', file.path]);
    const gasArgs = gasSources ? ['--gas-water', gasSources.waterSpecies ?? 'none',
      '--gas-phase', gasSources.phase, ...gasSources.components.flatMap(component =>
        ['--gas-component', component.speciesPath, component.nasa9Path, component.transportPath])] : [];
    const profileArgs = [['bodyPath', '--body-profile'], ['stopPath', '--stop-profile'], ['contactThermalPath', '--contact-thermal']]
      .flatMap(([key, flag]) => stringProfiles?.[key] ? [flag, stringProfiles[key]] : []);
    const coatingArgs = coatingStrainGuard ? ['--coating-strain-guard',
      String(coatingStrainGuard.maximumBytes), String(coatingStrainGuard.maximumWork)] : [];
    const history = join(dir, 'output.history.csv');
    const command = experiment.renderer === 'strings' ? [] : [experiment.renderer ?? 'wind'];
    const args = experiment.renderer === 'finite-gas-string'
      ? [input, finiteGasString.solidPath, output, '--gas-record', finiteGasString.gasRecordPath,
        '--self-contained', '--history-output', history, '--require-sample-rate', '48000', ...sourceArgs, ...gasArgs]
      : [...command, input, output,
        '--self-contained', ...(decimate ? ['--decimate'] : []), ...sourceArgs, ...gasArgs, ...profileArgs, ...coatingArgs];
    const outputs = [[output, outputLimit, 'WAV'],
      ...(experiment.renderer === 'finite-gas-string' ? [[history, historyLimit, 'history']] : [])];
    return runNativeCommand(engine, args, dir, signal, outputs, 'renderer');
}
function runNativeCommand(engine, args, dir, signal, outputs, label) {
  return new Promise((accept, reject) => {
    if (signal.aborted) return reject(signal.reason);
    const child = spawn(engine, args, {cwd:dir, shell:false,
      windowsHide:true, detached:process.platform !== 'win32', stdio:['ignore', 'pipe', 'pipe']});
    let reason, size = 0, forceKill, finished = false; const chunks = [], stdout = [], stderr = [];
    const kill = name => { if (!child.pid) return;
      try { if (process.platform === 'win32') child.kill(name); else process.kill(-child.pid, name); }
      catch (error) { if (error.code !== 'ESRCH') child.kill(name); } };
    const stop = error => { if (reason || finished) return; reason = error; kill('SIGTERM');
      forceKill = setTimeout(() => kill('SIGKILL'), 1000); forceKill.unref(); };
    const abort = () => stop(signal.reason);
    signal.addEventListener('abort', abort, {once:true});
    const read = (bytes, stream) => { size += bytes.length;
      if (size > MAX_LOG_BYTES) stop(new NativeWindError(502, 'Native diagnostics exceeded 256 KiB.'));
      else { chunks.push(bytes); stream.push(bytes); } };
    child.stdout.on('data', bytes => read(bytes, stdout));
    child.stderr.on('data', bytes => read(bytes, stderr));
    // Bound in-progress native output as well as the completed result.
    let checking = false;
    const check = setInterval(async () => {
      if (checking) return; checking = true;
      try {
        for (const [path, limit, name] of outputs) {
          try { const info = await lstat(path);
            if (!info.isFile() || info.size > limit) stop(new NativeWindError(502, `Native ${name} exceeds the output limit or is not a regular file.`));
          } catch (error) { if (error.code !== 'ENOENT') stop(new NativeWindError(502, 'Could not inspect native output.')); }
        }
      }
      finally { checking = false; }
    }, 100); check.unref();
    child.once('error', () => stop(new NativeWindError(503, 'Selected native executable could not be started.')));
    // The configured native executable needs no surviving helper processes.
    child.once('exit', () => kill('SIGKILL'));
    child.once('close', code => {
      finished = true;
      clearInterval(check); clearTimeout(forceKill); signal.removeEventListener('abort', abort);
      if (reason) reject(reason);
      else if (code !== 0) reject(new NativeWindError(400,
        `Native ${label} refused the source (exit ${code}). ${Buffer.concat(chunks).toString('utf8').slice(-3000)}`));
      else accept({stdout:Buffer.concat(stdout).toString('utf8'), stderr:Buffer.concat(stderr).toString('utf8')});
    });
    if (signal.aborted) abort();
  });
}
async function readWave(output, limit) {
  let file;
  try {
    file = await open(output, constants.O_RDONLY | (constants.O_NOFOLLOW || 0));
    const info = await file.stat();
    if (!info.isFile() || info.size > limit) throw new NativeWindError(502, 'Native WAV exceeds the output limit or is not a regular file.');
    return validateWave(await file.readFile(), limit);
  } catch (error) {
    if (error instanceof NativeWindError) throw error;
    throw new NativeWindError(502, 'Native renderer did not leave a readable WAV.');
  } finally { await file?.close(); }
}

async function readBoundedOutput(output, limit, name, utf8 = false) {
  let file;
  try {
    file = await open(output, constants.O_RDONLY | (constants.O_NOFOLLOW || 0) | (constants.O_NONBLOCK || 0));
    const info = await file.stat();
    if (!info.isFile() || !info.size || info.size > limit)
      throw new NativeWindError(502, `Native ${name} exceeds the output limit or is not a nonempty regular file.`);
    const bytes = Buffer.alloc(Math.min(info.size + 1, limit + 1));
    let length = 0;
    while (length < bytes.length) {
      const read = await file.read(bytes, length, bytes.length - length, length);
      if (!read.bytesRead) break;
      length += read.bytesRead;
    }
    if (length !== info.size) throw new NativeWindError(502, `Native ${name} changed while being read.`);
    const result = bytes.subarray(0, length);
    if (utf8) new TextDecoder('utf-8', {fatal:true}).decode(result);
    return result;
  } catch (error) {
    if (error instanceof NativeWindError) throw error;
    throw new NativeWindError(502, `Native executable did not leave readable ${utf8 ? 'UTF-8 ' : ''}${name}.`);
  } finally { await file?.close(); }
}

async function compileMaterialManifests(engine, experiment, dir, signal) {
  if (!experiment.materialManifests) return experiment;
  const attachments = [...experiment.attachments];
  let total = attachments.reduce((sum, file) => sum + file.bytes.length, 0);
  for (const {packPath, manifestPath, materialState} of experiment.materialManifests) {
    if (signal.aborted) throw signal.reason;
    const limit = Math.min(MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES - total);
    if (limit <= 0) throw new NativeWindError(413, 'Compiled packs and original sources exceed the combined attachment limit.');
    const output = join(dir, packPath);
    await mkdir(dirname(output), {recursive:true, mode:0o700});
    const args = ['matdb-pack', '--manifest', join(dir, manifestPath), '--out', output,
      ...(materialState ? ['--material-state', materialState.chemistry, materialState.phase,
        materialState.process, String(materialState.revision)] : [])];
    // The native compiler owns source-path confinement, material identity,
    // normalization and pack encoding. All inputs are original attachments.
    await runNativeCommand(engine, args, dir, signal, [[output, limit, 'compiled material pack']], 'material compiler');
    const bytes = await readBoundedOutput(output, limit, 'compiled material pack');
    total += bytes.length;
    attachments.push({path:packPath, bytes});
  }
  return {...experiment, attachments};
}

export function createNativeWindServer({engine = process.env.FSIM_MUSIC_RENDER || null,
  stringsEngine = process.env.FSIM_MATERIAL_STRINGS || null,
  finiteGasStringEngine = process.env.FSIM_FINITE_GAS_STRING || null,
  materialCompiler = process.env.FSIM_MATDB_COMPILER || null,
  dist = resolve(here, '../dist'), timeoutMs = DEFAULT_TIMEOUT_MS,
  maximumOutputBytes = MAX_OUTPUT_BYTES, maximumHistoryBytes = MAX_HISTORY_BYTES, tempRoot = tmpdir()} = {}) {
  if ([engine, stringsEngine, finiteGasStringEngine, materialCompiler].some(path => path && !isAbsolute(path))) throw new Error('Native executables must be absolute paths.');
  if (!Number.isInteger(timeoutMs) || timeoutMs < 1 || timeoutMs > 600000
      || !Number.isInteger(maximumOutputBytes) || maximumOutputBytes < 44 || maximumOutputBytes > MAX_OUTPUT_BYTES
      || !Number.isInteger(maximumHistoryBytes) || maximumHistoryBytes < 1 || maximumHistoryBytes > MAX_HISTORY_BYTES)
    throw new Error('Invalid local rendering budget.');
  dist = resolve(dist); tempRoot = resolve(tempRoot);
  const token = randomBytes(32).toString('hex'); let active = null, closing = false;
  const server = http.createServer(async (req, res) => {
    const send = (status, value) => { if (res.destroyed || res.writableEnded) return;
      res.writeHead(status, {'Content-Type':'application/json; charset=utf-8', ...(status === 413 ? {'Connection':'close'} : {})});
      res.end(JSON.stringify(value)); };
    res.setHeader('Cache-Control', 'no-store'); res.setHeader('X-Content-Type-Options', 'nosniff');
    res.setHeader('Referrer-Policy', 'no-referrer');
    res.setHeader('Content-Security-Policy', "default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; media-src 'self' blob:; font-src 'self' data:; worker-src 'self' blob:; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'");
    const port = server.address()?.port, host = req.headers.host;
    if (![`127.0.0.1:${port}`, `localhost:${port}`, `[::1]:${port}`].includes(host)
        || (req.headers.origin && req.headers.origin !== `http://${host}`)
        || (req.headers['sec-fetch-site'] && !['same-origin', 'none'].includes(req.headers['sec-fetch-site'])))
      return send(403, {error:'Only same-origin loopback requests are accepted.'});
    let pathname;
    try { pathname = new URL(req.url, `http://${host}`).pathname; }
    catch { return send(400, {error:'Invalid request path.'}); }
    try {
      if (req.method === 'GET' && pathname === '/api/native-wind/status') {
        const [wind, strings, finiteGasString, compiler] = await Promise.all([available(engine), available(stringsEngine), available(finiteGasStringEngine), available(materialCompiler)]);
        return send(200,
        {available:!closing && (wind || strings || finiteGasString),
          renderers:{wind:!closing && wind, strings:!closing && strings, resolved:!closing && wind,
            'finite-gas-string':!closing && finiteGasString}, busy:active !== null, token,
          maximumSourceBytes:MAX_SOURCE_BYTES, maximumOutputBytes, maximumHistoryBytes, timeoutMs,
          supportedSchemas:[...(wind ? [...SUPPORTED_SCHEMAS, ...RESOLVED_SCHEMAS] : []), ...(strings ? STRING_SCHEMAS : []),
            ...(finiteGasString ? FINITE_GAS_STRING_SCHEMAS : [])],
          maximumAttachments:MAX_ATTACHMENTS, maximumAttachmentFileBytes:MAX_ATTACHMENT_FILE_BYTES,
          maximumAttachmentBytes:MAX_ATTACHMENT_BYTES, maximumRequestBytes:MAX_REQUEST_BYTES, supportsGasSources:true, supportsPressureReceipts:true,
          supportsCoatingStrainGuard:true, supportsMaterialManifests:!closing && compiler});
      }
      if (pathname !== '/api/native-wind/render' || req.method !== 'POST') {
        if (['GET', 'HEAD'].includes(req.method) && !pathname.startsWith('/api/') && await serveAsset(req, res, dist, pathname)) return;
        return send(404, {error:'Not found. Build Instrument Lab and pass its dist directory to this server.'});
      }
      if (closing) return send(503, {error:'Local renderer is shutting down.'});
      if (req.headers['x-native-wind-token'] !== token) return send(403, {error:'Missing local session token. Refresh renderer availability.'});
      const media = /^(text\/plain|application\/json)(?:\s*;\s*charset\s*=\s*utf-8)?\s*$/i.exec(req.headers['content-type'] || '');
      if (!media || (req.headers['content-encoding'] && req.headers['content-encoding'] !== 'identity'))
        return send(415, {error:'Send UTF-8 performance text or an attachment experiment as JSON.'});
      const bundled = media[1].toLowerCase() === 'application/json';
      if (Number(req.headers['content-length'] || 0) > (bundled ? MAX_REQUEST_BYTES : MAX_SOURCE_BYTES))
        return send(413, {error:'Experiment exceeds the complete request limit.'});
      if (active) return send(409, {error:'Another native render is running. Cancel it or wait.'});
      const controller = new AbortController(); let complete;
      const job = {controller, finished:new Promise(resolve => { complete = resolve; })}; active = job;
      const disconnected = () => { if (!res.writableEnded) controller.abort(new NativeWindError(499, 'Browser disconnected; render cancelled.')); };
      res.once('close', disconnected);
      const timer = setTimeout(() => controller.abort(new NativeWindError(504, 'Native render exceeded its execution budget.')), timeoutMs);
      let dir;
      try {
        let experiment = await sourceBody(req, controller.signal, bundled);
        const selectedEngine = experiment.renderer === 'strings' ? stringsEngine
          : experiment.renderer === 'finite-gas-string' ? finiteGasStringEngine : engine;
        if (!await available(selectedEngine)) throw new NativeWindError(503, experiment.renderer === 'strings'
          ? 'Native material-string executable is unavailable. Build material_strings and configure its absolute path.'
          : experiment.renderer === 'finite-gas-string'
            ? 'Native finite-gas string executable is unavailable. Build finite_gas_string and configure its absolute path.'
          : 'Native wind/resolved-flow executable is unavailable. Build music_render and configure its absolute path.');
        if (experiment.materialManifests && !await available(materialCompiler))
          throw new NativeWindError(503, 'Native material compiler is unavailable. Build xtask and configure --material-compiler with its absolute path.');
        if (controller.signal.aborted) throw controller.signal.reason;
        dir = await mkdtemp(join(tempRoot, 'frankensim-native-wind-'));
        const input = join(dir, 'input.performance'), output = join(dir, 'output.wav');
        await writeFile(input, experiment.source, {mode:0o600, flag:'wx'});
        for (const attachment of experiment.attachments) {
          if (controller.signal.aborted) throw controller.signal.reason;
          const destination = join(dir, attachment.path);
          await mkdir(dirname(destination), {recursive:true, mode:0o700});
          await writeFile(destination, attachment.bytes, {mode:0o600, flag:'wx'});
        }
        experiment = await compileMaterialManifests(materialCompiler, experiment, dir, controller.signal);
        const diagnostics = await runNative(selectedEngine, input, output, dir, controller.signal, maximumOutputBytes, experiment, maximumHistoryBytes);
        const wave = await readWave(output, maximumOutputBytes);
        if (controller.signal.aborted) throw controller.signal.reason;
        let receipt;
        try { receipt = await nativePressureReceipt(experiment.renderer ?? 'wind', output, diagnostics, wave); }
        catch (error) { throw new NativeWindError(502, error.message); }
        if (controller.signal.aborted) throw controller.signal.reason;
        if (experiment.renderer === 'finite-gas-string') {
          if (!receipt) throw new NativeWindError(502, 'Native finite-gas string pressure calibration is missing.');
          const history = await readBoundedOutput(join(dir, 'output.history.csv'), maximumHistoryBytes, 'history', true);
          if (controller.signal.aborted) throw controller.signal.reason;
          send(200, {format:FINITE_GAS_STRING_RESULT_FORMAT,
            wave:wave.toString('base64'), history:history.toString('base64'), receipt});
          return;
        }
        res.writeHead(200, {'Content-Type':'audio/wav', 'Content-Length':wave.length,
          ...(receipt ? {[PRESSURE_RECEIPT_HEADER]:JSON.stringify(receipt)} : {})}); res.end(wave);
      } finally {
        clearTimeout(timer); res.off('close', disconnected);
        try { if (dir) await rm(dir, {recursive:true, force:true}); } // Only this request's newly created private directory.
        finally { if (active === job) active = null; complete(); }
      }
    } catch (error) { send(error.status || 500, {error:error instanceof NativeWindError ? error.message : 'Local rendering failed.'}); }
  });
  server.headersTimeout = 10000; server.requestTimeout = 15000;
  return {server, async shutdown() {
    closing = true; const job = active;
    job?.controller.abort(new NativeWindError(503, 'Local renderer is shutting down.'));
    const closed = new Promise(resolve => server.close(() => resolve()));
    server.closeAllConnections(); await job?.finished; await closed;
  }};
}

async function cli() {
  let engine = process.env.FSIM_MUSIC_RENDER || null, dist = resolve(here, '../dist'), port = 8787;
  let timeoutMs = DEFAULT_TIMEOUT_MS;
  let stringsEngine = process.env.FSIM_MATERIAL_STRINGS || null;
  let finiteGasStringEngine = process.env.FSIM_FINITE_GAS_STRING || null;
  let materialCompiler = process.env.FSIM_MATDB_COMPILER || null;
  const args = process.argv.slice(2);
  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--help') {
      console.log('node tools/serve-native-wind.mjs [--engine /absolute/path/music_render] [--strings-engine /absolute/path/material_strings] [--finite-gas-string-engine /absolute/path/finite_gas_string] [--material-compiler /absolute/path/xtask] [--dist DIST] [--port 8787] [--timeout-ms 60000]\nFSIM_MUSIC_RENDER, FSIM_MATERIAL_STRINGS, FSIM_FINITE_GAS_STRING and FSIM_MATDB_COMPILER select the local executables. --timeout-ms admits 1..600000 milliseconds per render. Build the app first with npm run build.'); return;
    }
    if (!['--engine', '--strings-engine', '--finite-gas-string-engine', '--material-compiler', '--dist', '--port', '--timeout-ms'].includes(args[i]) || !args[i + 1]) throw new Error('Unknown or incomplete local-server option.');
    const option = args[i], value = args[++i];
    if (option === '--engine') engine = value; else if (option === '--strings-engine') stringsEngine = value;
    else if (option === '--finite-gas-string-engine') finiteGasStringEngine = value;
    else if (option === '--material-compiler') materialCompiler = value;
    else if (option === '--dist') dist = resolve(value);
    else if (option === '--timeout-ms') timeoutMs = Number(value);
    else port = Number(value);
  }
  if (!Number.isInteger(port) || port < 1024 || port > 65535) throw new Error('Port must be 1024–65535.');
  const app = createNativeWindServer({engine, stringsEngine, finiteGasStringEngine, materialCompiler, dist, timeoutMs});
  app.server.on('error', error => { console.error(error.message); process.exitCode = 1; });
  app.server.listen(port, '127.0.0.1', () => console.log(`Instrument Lab + native physics: http://127.0.0.1:${port}\nLocal only; Ctrl+C stops rendering and the server.`));
  process.once('SIGINT', () => void app.shutdown()); process.once('SIGTERM', () => void app.shutdown());
}
if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  cli().catch(error => { console.error(error.message); process.exitCode = 1; });
}
