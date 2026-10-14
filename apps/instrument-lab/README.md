# Struck — FrankenSim instrument lab

A public demo of FrankenSim's physically simulated drum kit and concert grand.
Live at https://frankensim.org/music and https://frankensim-struck.vercel.app.

Two tracks, always labelled:

- **Offline reference** — clips rendered natively by the `percussion` and `grand_piano`
  examples in `crates/fs-couple/examples`, with every clip's command, render time and energy
  closure parsed from its own log. One fixed playback gain per comparison; nothing is
  normalized, equalized or reverberated.
- **Live** — WebAssembly engines that compute sound in the browser: `engine-kit` (the
  `fs-vkplate` drum-kit solver itself on a coarser grid, one Worker per piece) and
  `engine-piano` (a lean kernel that steps the grand piano's own equations for the sounding
  keys only, validated against the reference `Instrument`). `engine-drum` is the earlier modal
  drum (FE heads, eigen-solver, prepared linear impact engine), no longer on the page.

## Layout

| Path | What |
|---|---|
| `index.html`, `src/` | the page (Vite, TypeScript, three.js; no framework) |
| `src/live/` | live section UIs and typed clients for the engine workers |
| `public/live/` | engine wasm + module workers + the AudioWorklet sink (static, unbundled) |
| `engine-kit/`, `engine-piano/`, `engine-drum/` | the live engines (standalone cargo workspaces) |
| `public/live/mix-worklet.js`, `lane.js` | `fs-live-mix`: sums frame-stamped lanes from several engine Workers in lock-step with the AudioContext clock |
| `tools/strandbar.py` | builds the Strandbar studio loop and the live engines' score and levels (`data/strandbar.json`) |
| `tools/build_data.py` | turns a render campaign into `public/{audio,img,data,viz}` |
| `tools/diagnose.py` | traces each real-vs-ours gap to a link of the signal chain and scores ablation renders: `public/data/diagnosis.json` |
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

## Diagnosing the gaps

`tools/diagnose.py` (same arguments as `verisim.py`, whose measurement code it imports) asks
WHY each gap exists. It uses the A3 bridge trace (`--bridge-trace-csv`) to separate what the
board does from what radiation to the microphone does. It also scores ablation renders, which
switch one mechanism on at a time, whenever they are present in the campaign:

    # receivers near the Iowa microphone positions (stereo)
    grand_piano ...a3-mf flags... --microphone 0.529,0.498,0.233 --microphone-right 0.640,0.473,0.233 --render a3-mf-close.wav
    # longitudinal string waves on A3 only; every other key `linear,<key>` (88 rows)
    #   waves,57,166505,0.2,10000,8,0.01,0.01   EA = E*A, E = 200 GPa, A from the scale's own EI
    grand_piano ...a3-{mf,ff} flags... --string-stretching axial-a3.fspx --render a3-{mf,ff}-waves.wav

A card counts as confirmed only when its experiment moves the measured gap. Air loading
(`--rayleigh-loading`) cannot yet be admitted on the production board; see bead
frankensim-s59w0.

## The drum kit

Every kit piece comes from one solver, `crates/fs-vkplate` (its own workspace, std only; see its
README for the physics). Render each piece, then grade them against the Salamander Drumkit
(Alexander Holm, CC BY-SA 3.0; `archive.org/details/SalamanderDrumkit`, extract `OH/`):

    cd crates/fs-vkplate && cargo build --release --example kit
    for p in kick snare snare-off tom-12 tom-14; do kit $p $KIT/$p.wav; done        # 1-5 min each
    kit hihat-closed $KIT/hihat-closed.wav 0 0.6; kit hihat-open $KIT/hihat-open.wav 0 2.5
    kit crash $KIT/crash.wav 0 3; kit ride $KIT/ride.wav 0 3                       # ~20 min each
    uv run --with numpy --with pillow python tools/kitgrade.py $KIT $SALAMANDER/OH public

`kitgrade.py` writes `data/kit.json`, level-matched excerpts, spectrograms and the two grooves
(the same two bars on the real kit and on FrankenSim's). Its ruler is real-vs-real: another stroke
of the same real drum at the same dynamic, and the same drum one dynamic away. Octave-band balance
is floored at -50 dB, because the real overheads carry room noise and cymbal bleed where a render
is silent. `diagnose.py` takes `$KIT $SALAMANDER/OH` as optional fourth and fifth arguments to run
the cymbal-shimmer experiment on the kit's crash.

## The live kit

`engine-kit/build-wasm.sh` compiles `fs_vkplate::live` (wasm-bindgen) into `public/live/kit/`. Each
piece gets its own Worker (`public/live/kit/worker.js`), which times a stroke at grid scales 4-8
and keeps the finest that runs 1.7x faster than real time on the device. A silent piece is not
stepped. Strokes land where the head or cymbal is touched; the Worker sends the surface velocity
at ~45 Hz and, after each stroke, the first 8 ms every 0.25 ms, which the page replays 60x slower.
Levels: one fixed gain per piece, its 50 ms loudness matched to the Salamander overheads, capped so
a standard stroke peaks at 0.55 (the live cymbals' attacks are too sharp to match their loudness).

## Strandbar

The studio loop is three native renders that start 0.15 s before the first downbeat:
`grand_piano --performance strandbar-piano.csv` (Steinway D preset, 128 board coordinates up to
2.4 kHz, one substep) and fs-vkplate `perform` for the kick and the closed hi-hat. Long fs-vkplate
parts can be rendered as bar-sized simulations, each with a pre-roll of real strokes, and stitched
(the stitched kick matched a continuous render to 0.9 % RMS). The hi-hat's banded Airy solve is
memory-bound: run one or two such processes per L3 cache (`taskset`, `-C target-cpu=native`).

    uv run --with numpy --with soundfile python tools/strandbar.py ARRANGEMENT.json STEMS_DIR RECORDING.wav public

It balances the three parts against the recording's low/middle/high band energies, folds the
tail onto the loop's start, writes `audio/strandbar-studio.flac` and `data/strandbar.json` (the
same score for the live piano and kit lanes), and gives each live part the RMS of its studio part
by running the live engines' WebAssembly in Node (`tools/strandbar_live_rms.mjs`).
