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
| `public/live/` | engine wasm + module workers + the AudioWorklet sink (static, unbundled; a build serves them under `live/<hash>/`, see Rebuilding) |
| `engine-kit/`, `engine-piano/`, `engine-drum/` | the live engines (standalone cargo workspaces) |
| `render-gpu/` | offline kit renders on an NVIDIA GPU when one is present (one persistent CUDA kernel stepping the exported fs-vkplate instrument; batches, chunks, both GPUs), on the CPU otherwise; see its README for the measured agreement and speed |
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
3. `engine-kit/build-wasm.sh` (it also writes `engine_levels.json`), `engine-piano/build-wasm.sh`
   (run it in the export: it links fs-couple by path) and, for the older drum section,
   `engine-drum/build-wasm.sh`; they install into `public/live/`.
4. `npm run build`, then `vercel deploy --prod --yes` (prebuilt `dist/`; see `vercel.json`).

The live files are never bundled, and each Worker loads its glue, wasm and data relative to its own
URL. So a build (`vite.config.ts`) serves the whole of `public/live` under one directory named by a
hash of its contents, `live/<hash>/`, and the page's `LIVE_BASE_URL` points there: a page always gets
Workers, glue and wasm from one build, and `vercel.json` and `cloudflare/_headers` cache `live/`
immutably. The build fails if what it copied differs from what it hashed. A page left open across a
deploy asks for its old `live/<hash>/`, which is gone, and says to reload. The dev server serves
`public/live` at `live/`.

The page's data, audio, images, replays and MIDI rolls keep their paths and a one-hour cache, but a build
hands the page each file's content hash and `util.asset()` asks for `path?v=<hash>`: a file a deploy
changed gets a new URL, so no visitor pairs new code with an hour-old `strandbar.json` or recording.

## Grading against real recordings

The "Real vs ours" section compares renders with University of Iowa Musical Instrument Samples
(<https://theremin.music.uiowa.edu/MIS.html>, free to use without restriction): a Steinway Model B
(`Piano.{pp,mf,ff}.A3`, `Piano.mf.C6`) and an anechoic splash plus 13/17/18-inch crashes
(`splash.stick.normal.*`, `*crash.stick.normal.mf`). Convert each to 48 kHz float WAV. The campaign's
splash renders come from the kit's solver, `perc/vk-splash.wav` and `perc/vk-splash-soft.wav`
(`kit splash OUT.wav 2 2.6` and `kit splash OUT.wav 0.6 2.6`, see "The drum kit"); a pair's label
names what its render simulated. Then:

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

    # the same A3 mf heard where the Iowa microphones sat (stereo; see "Listening point" below)
    grand_piano ...a3-mf flags... --microphone 0.275,0.556,0.25 --microphone-right 0.867,0.147,0.233 --render a3-mf-close.wav
    # ...and by stereo pairs 0.1 m either side of the bass-side point, along x and along y
    grand_piano ...a3-mf flags... --microphone 0.175,0.556,0.25 --microphone-right 0.375,0.556,0.25 --render a3-mf-px.wav
    grand_piano ...a3-mf flags... --microphone 0.275,0.456,0.25 --microphone-right 0.275,0.656,0.25 --render a3-mf-py.wav
    # the former estimated bridge stations on the same mesh, with the same air and receiver
    piano_geometry --divisions 16 --rib-normal-step-m 0.03 --bridge-map estimated --output board-est16.fsb
    grand_piano ...a3-mf flags... --board-geometry board-est16.fsb --render a3-mf-estmap.wav
    # bare-board modes and per-key bridge displacements of the graded board (A3's driving-point mobility)
    grand_piano --preset steinway-d --board-geometry board-rt16.fsb --board-band-hz 1600 --dump-board board-modes.csv
    # longitudinal string waves on A3 only; every other key `linear,<key>` (88 rows)
    #   waves,57,166505,0.2,10000,8,0.01,0.01   EA = E*A, E = 200 GPa, A from the scale's own EI
    grand_piano ...a3-{mf,ff} flags... --string-stretching axial-a3.fspx --render a3-{mf,ff}-waves.wav

A card counts as confirmed only when its experiment moves the measured gap. A mechanism that has
moved into the graded renders is reported as fixed or partly fixed: the A3 renders
(`a3-{pp,mf,ff,fff}`) carry `--string-stretching axial-a3.fspx`, whose waves drive the bridge axially,
and when the campaign also has the same strokes on linear strings (`a3-{mf,ff,fff}-linear`) the card
compares the two. Every piano render the page grades (A3, C6), and every ablation above except in
the one thing it changes, couples each string to its bridge at the point RT-0425 publishes for it,
carries two-sided air loading and listens at the default point, 1 m above the middle of the board:

    piano_geometry --divisions 16 --rib-normal-step-m 0.03 --bridge-map rt0425 --output board-rt16.fsb  # 11,772 area samples
    grand_piano --preset steinway-d --board-geometry board-rt16.fsb --board-band-hz 2400 \
      --board-reduction 128,64,200,400,800,1600,2300 --rayleigh-loading 2,90,2620,213 \
      --rt0425-hammer-stiffness --rt0425-hammer-dissipation --rt0425-string-damping ...

The hammers follow RT-0425's published law: its stiffness K_H acts on each unison string (the preset's
former felt conversion divided it across the unison, which its own comment calls too soft) with the
published dissipation R_H, and the strings take RT-0425's published losses. The same mf/ff/8 m/s strokes
with the former felt law, `a3-{mf,ff,fff}-legacyfelt`, let the ff-brightness card show what changed.

The former estimated stations (`--bridge-map estimated`) were spaced evenly along each bridge. With
the published string lengths they put most agraffes over the soundboard itself (A3's 0.77 m string,
running straight back from its contact, would start 0.38 m inside the board's front edge) and
reversed the bass order. The reduced board keeps a coordinate at 2.60 kHz, so the air band runs to
2.62 kHz. Of the meshes tried, 12 to 15 divisions left the air cells under six per wavelength there;
16 divisions with a 0.03 m step is the coarsest that passes. The polarization ablation adds bridge
directions that reach 2.62 kHz, so it runs its air to 2.64 kHz (215 samples).

Listening point: the Iowa recordings used two Neumann KM 84 cardioids, "8 in above center bass
strings" (left channel) and "8 in above center treble strings" (right). In the board frame those are
the centroids of the bass and treble strings' speaking lengths, 0.203 m above the strings: the
published contacts, the agraffe line the treble strings imply, and the bass agraffes continuing that
line. The points are estimates, and ours are omnidirectional. The page does not grade there. Real
microphones at those points sit above the strings and the cast-iron plate, which covers much of the
board, usually under a lid; the model's receivers look at a bare board in an infinite baffle. Heard
there, the A3 partial pattern lies 14 dB from the recording on the matching channel, against about
8 dB at the default point, and averaging both channels' power does not change that. The real
piano's own two microphones differ by 8.4 dB. `a3-mf-close` keeps those positions as the receiver
experiment.

When the campaign also has the same A3 mf stroke without air on that board (`a3-mf-noair`), the
low-resonance and ringing cards compare the two. Each render first samples and fits the board's
complete 128-mode air reaction (15 minutes or more on the 16-division board; see the example's
RAYLEIGH_LOADING.md).

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

## Larger live soundboards

The live piano embeds a 32-coordinate board; `public/live/piano/boards/` holds 48/64/96/128-coordinate
boards (not committed, ~6 MB) that the browser fetches only when its timing probe says it can afford
them (`boards.js`). The score lanes (Strandbar, the player) choose board and mode cap together on a
sustained probe (`chooseQuality`); the keyboard section chooses the largest board whose 12-key probe runs
1.4x faster than real time (`chooseBoard`) and leaves the mode budget to its load governor. For the lanes,
a Chromium desktop (8+ logical cores, `userAgentData.mobile` false, no data saving) starts fetching board 96,
their usual pick there, while the embedded board is timed, so the fetch is off the start-up path. Each board the probe times is built once: the Worker keeps a copy
(`PianoEngine.copy()`) from before the timing and plays the chosen board's copy, which renders the same
samples as a new engine (building board 96 takes ~0.9 s, a copy under a millisecond); an engine that cannot
copy (a wasm from an older build) builds the chosen board again. On a device too slow for even the
embedded board at the 700-mode start floor, the piano falls back to its 16 exact low modes, started at
the mode cap its own timing predicts (not below 700) rather than the full 2,400; the load governor takes
it lower from there. The boards are the offline
binary's reductions of the same FE slice; generate them natively with the engine's `gen_board` example:

    gen_board public/live/piano/boards/model-d-board-2400-r64.bin 2400 64,32,200,400,800,1600,2300
    # r48: 48,24,...  r96: 96,48,...  r128: 128,64,... (the studio render's reduction)

Measured on the Strandbar piano part against the studio stem (RMS octave-balance error, 63 Hz-8 kHz):
32 coordinates 2.95 dB, 48 2.66, 64 1.28, 96 0.77, 128 0.36.

## The live kit

`engine-kit/build-wasm.sh` compiles `fs_vkplate::live` (wasm-bindgen) into `public/live/kit/`. Each
piece gets its own Worker (`public/live/kit/worker.js`), which times a stroke at grid scales 2.5-8
(the snare's 6-8: its 5x and 4x sound no closer to the 1x grid than 6x; the hats' also 4.5x: on 5x
they ring a tonal line near 558 Hz, which 4.5x weakens or drops; the kick's 6-8, the rack tom's 5-6
and the floor tom's 4-6, the grids past which they sound no closer to the 1x grid, so their CPU
goes to the pieces that do) and keeps the finest that
runs 1.7x faster than real time on the device. It times from the coarse
end and stops at the first grid that is too slow, so the costly fine grids are never built. It builds
each grid once and keeps a copy from before the timing stroke, and the chosen grid's copy plays: a
snare's or hi-hat's build is almost all its settle (wire preload, pedal; 0.05-0.9 s), a copy takes
under a millisecond and plays the same samples. An engine that cannot copy (a wasm from an older
build) gets the chosen grid built again, which plays the same samples. Each
grid is timed over the first 112 ms of a stroke, the snare's over 450 ms: its wires touch the head
from ~85 ms on and its steps then cost ~1.6x as much, so the short test promised ~1.5x the speed it
kept under rapid tapping, and the snare dropped audio. The player (`src/player/`) times its snare
the same way, at 6x or else 8x, against 1.5 (its lanes render a second ahead). A silent piece is not stepped. Strokes land where the head or cymbal is touched; the Worker sends
the surface velocity at ~45 Hz and, after each stroke, the first 8 ms every 0.25 ms, which the
page replays 60x slower.
Levels: one gain per piece, set on the 5x grid: its 50 ms loudness matched to the Salamander
overheads, capped so a standard stroke peaks at 0.55 (the live cymbals' attacks are too sharp to
match their loudness). A stroke's loudness moves with the grid, by piece and not monotonically (the
crash is 5 dB louder at 4x than at 5x, the closed hi-hat 5 dB quieter; the drums move by under
2 dB, the snare by up to 2.9 dB at 2x), so `build-wasm.sh` also writes `engine_levels.json`: each
piece's standard stroke per grid 2-10x, its loudness (largest 400 ms RMS, so a cymbal's wash counts)
and peak, measured with that wasm by `tools/kit_levels.mjs` and stamped with its SHA-256. The
Worker scales the page's gain so the stroke on the grid it picked is as loud as on 5x
(`public/live/kit/levels.js`): a cut always, a boost only while the stroke still peaks at or under
0.55. So the closed hi-hat stays up to 3 dB quiet at 3.5-4x, and the hats and cymbals at 6x and
coarser, whose strokes already peak over 0.55, are not raised (up to 4.3 dB quiet). Beside any
other wasm the table is ignored and the gain is used as given. The same applies wherever a kit Worker plays
(Strandbar's gains are also measured on 5x).

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
by running the live engines' WebAssembly in Node (`tools/strandbar_live_rms.mjs`, board 32 and the 5x
grid). One master level serves both versions: the studio loop peaks at -1 dBFS at most and the live
parts' sum at -2 dBFS, so the boards and grids a device picks stay under full scale.

On the page the live lanes render two seconds ahead in 512-frame messages. Each drum keeps the finest
grid whose timing reaches 1.4x real time (`src/strandbar.ts`): the kick from 2-8x, the hi-hat from
2.5-4.5x (4.5x drops the 558 Hz line 5x rings), 6x and then 8x; 5x is skipped, where a closed stroke
keeps a 250 Hz component ringing ~4 s even with the rims' friction. Strandbar does not restart
a lane on a coarser grid while it plays, so the lists end at 8x for slow devices.

## Player piano

`src/player/` plays any MIDI file on the live piano and kit: a paper roll over a tracker bar, the
88 keys, the soundboard's modal field in 3D and each drum piece's surface velocity in 3D.
`src/midi/` reads the file: `midi-salvage.ts` is copied from the owner's jazz chord editor
(`jazz_chord_progression_editor_html`), `smf.ts` reads the result under that editor's M0 laws plus
the controllers a player needs, `chords.ts` names chords with its reverse-T1 templates, and
`arrange.ts` maps the file (pitched channels to the piano, GM channel 10 to kit pieces, velocity to
hammer or stick speed, CC64 to the pedal) and counts every fold, merge, thinning and skipped drum. `orchestrate.ts` turns a General MIDI
arrangement into one pianist's part first: the file's own channel volume and expression (CC7, CC11)
shape every strike, parts that double or echo another part are left out, chiptune ornaments are
thinned, parts written at one flat velocity get dynamics by role, and the page lists what happened
to every part. The shelf's second pack (`public/midi/game-pack*`) holds video-game files of that
kind, supplied by the owner.

The shelf's rolls are that editor's chord-chart library played through its own bass and comping
generators, with drum parts written per groove (from the frankensim root; the output in
`public/midi/` is committed):

    bun apps/instrument-lab/tools/standards/build-standards.ts ../jazz_chord_progression_editor_html apps/instrument-lab/public
    uv run --with mido python apps/instrument-lab/tools/standards/verify_standards.py apps/instrument-lab/public

The lanes take a whole score at once (`{type:'score', events}` with absolute frames, sorted once in
the Worker and fed to the queue only as far as the look-ahead), and a stop or a jump flushes what
they rendered ahead (`mix.flush({epoch})` + `lane.resync(epoch)`; audio, snapshots and fields carry
the epoch they were computed in). The piano lane sends `engine.snapshot()` per 1/30 s of audio
when asked (`snapshotHz`), and a kit lane its fields per `fieldFrames` of audio.
