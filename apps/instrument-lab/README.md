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
| `tools/verisim.py` | grades renders against real recordings: `public/data/verisim.json` plus `vs-*` spectrograms and level-matched clips |
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

## Grading against real recordings

The "Real vs ours" section compares renders with University of Iowa Musical Instrument Samples
(<https://theremin.music.uiowa.edu/MIS.html>, free to use without restriction): a Steinway Model B
(`Piano.{pp,mf,ff}.A3`, `Piano.mf.C6`) and an anechoic splash plus 13/17/18-inch crashes
(`splash.stick.normal.*`, `*crash.stick.normal.mf`). Convert each to 48 kHz float WAV, then:

    uv run --with numpy --with pillow python tools/verisim.py <campaign-dir> <wav48-dir> public

Both signals of a pair get the same 40–70 Hz zero-phase high-pass and onset alignment. The script
asserts that both notes start within 12 ms of their windows. Every gap is reported beside a ruler
built only from real-vs-real distances: the Iowa piano across dynamics, and the other real cymbals
against the splash. The page's findings list is computed from `verisim.json`, so it cannot drift
from the numbers.
