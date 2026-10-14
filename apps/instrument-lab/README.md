# Struck — FrankenSim instrument lab

A public demo of FrankenSim's physically simulated drum, splash cymbal and concert grand.
Live at https://frankensim-struck.vercel.app.

Two tracks, always labelled:

- **Offline reference** — clips rendered natively by the `percussion` and `grand_piano`
  examples in `crates/fs-couple/examples`, with every clip's command, render time and energy
  closure parsed from its own log. One fixed playback gain per comparison; nothing is
  normalized, equalized or reverberated.
- **Live** — two WebAssembly engines that compute sound in the browser:
  `engine-drum` (FrankenSim's own FE heads, eigen-solver and prepared linear impact engine;
  its mechanics output is byte-identical to the reference binary at the same clock) and
  `engine-piano` (a lean kernel that steps the grand piano's own equations for the sounding
  keys only, validated against the reference `Instrument`).

## Layout

| Path | What |
|---|---|
| `index.html`, `src/` | the page (Vite, TypeScript, three.js; no framework) |
| `src/live/` | live section UIs and typed clients for the engine workers |
| `public/live/` | engine wasm + module workers + the AudioWorklet sink (static, unbundled) |
| `engine-drum/`, `engine-piano/` | the live engines (standalone cargo workspaces) |
| `tools/build_data.py` | turns a render campaign into `public/{audio,img,data,viz}` |
| `tools/percussion-viz-export.patch` | inert-unless-`FSIM_VIZ_DIR` export of mode shapes and modal coordinates used by the solver-state replays |

`public/{audio,img,data,viz}` and `public/live/*/engine*` are generated and not committed.

## Rebuilding

The workspace needs its sibling repos at the `constellation.lock` pins. If the shared sibling
checkouts have drifted, build from an export with pinned siblings beside it
(`git archive HEAD` for frankensim, `git archive <pin>` for each sibling).

1. Build the examples (release) with the viz patch applied in the export, then render the
   campaign; the exact command of every clip is in its "Run evidence" drawer on the page.
2. `npm run data -- <campaign-dir> public`
3. `engine-drum/build-wasm.sh` and the engine-piano build (see `engine-piano/`), which install
   into `public/live/`.
4. `npm run build`, then `vercel deploy --prod --yes` (prebuilt `dist/`; see `vercel.json`).
