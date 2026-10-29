# Instrument Explorer

A web dashboard for FrankenSim's instrument program, public at
<https://frankensim.org/instrument-tests/>. It plays what the simulated instruments can play
(scales, double stops, melodies, grooves, rows of notes) with their spectrograms and the measured
pitch of every note. For every family it puts our renders next to real recordings, shows how far
apart they are, and follows how each family was built and qualified.

**No audio is committed to this repository, and none ever should be.** The audio and maps live
under `~/.local/share/instrument-explorer/` and in the ignored `dist-cf/`. The public build plays
every recording whose licence allows it, credited next to the clip: University of Iowa MIS (free
to use without restriction), CC0, CC BY, CC BY-SA, YouTube Creative Commons and GPL FreePats. The
49 YouTube clips under the standard licence (all rights reserved) were used only to grade. They are
listed with their source link but never published.

## Running it

```sh
apps/instrument-explorer/tools/sync.sh          # mirror renders, grades and recordings from the hosts
python3 apps/instrument-explorer/tools/build.py # turn the newest mirror into the site data
python3 apps/instrument-explorer/tools/serve.py --open   # http://localhost:8742
```

`sync.sh` and `build.py` are only needed when the program has moved; `serve.py` alone serves the
last build. Reload the page after a rebuild.

The phrases come from the newest `phrases-<UTC>/` directory next to the mirrors (`build.py
--phrases`). Each holds `phrases.json` plus a WAV and a spectrogram per phrase. They are made on
a render host from the committed instrument descriptions:

1. `tools/phrases/phrases.py` writes one `music_render` job per note and a plan for joining them.
2. The jobs run.
3. `tools/phrases/phrasekit.py` joins the renders (detaché, tongued, legato or held, depending on
   the instrument), draws the spectrograms and measures each note's pitch on the joined audio
   (`pitchlib.py`).

The piano, the kit and the trumpet phrases are single continuous renders.

## Publishing

```sh
python3 apps/instrument-explorer/tools/publish.py   # dist-cf/instrument-tests/: AAC audio, licence filter
npx wrangler deploy --config apps/instrument-explorer/cloudflare/wrangler.jsonc
```

The share cards (`share/og.png` 1200x630, `share/twitter.png` 1200x600) are drawn from one phrase's
real analysis by `python3 apps/instrument-explorer/tools/share/card.py` (default the trumpet fanfare,
`--phrase`/`--title` for another; it needs Google Chrome). Re-run it after a phrase re-render so the
card shows the current sound; `publish.py` refuses to run without both images.

`publish.py` converts every FLAC the data points at to AAC (`.m4a`, Apple's encoder, 160 kbit/s).
It drops the audio of recordings that may not be redistributed, rewrites licence texts to say what
each licence allows, and strips build-host paths. It refuses to finish if any text names the
private product. The names are regular expressions in `~/.config/instrument-explorer/private-names`, one per line, kept
out of this repository, and it refuses to start without that file. Only files the data points at are copied, and outputs newer than their source
are reused. Files an earlier build left behind are never deleted; they are listed in
`dist-cf/.assetsignore` so they are not uploaded. The assets-only Worker `frankensim-instrument-tests` serves
`frankensim.org/instrument-tests/*`; everything else on the domain is untouched. New assets can
take about 30 s to appear after a deploy.

## What is in it

- **Listen** (home): a featured player and every phrase as a card you can play in place, filtered
  by scales, chords, melodies and grooves, or rows of notes. Each phrase page has the spectrogram
  with the notes it asked for drawn on it (dashed where asked, solid where it actually sounds when
  more than 50 cents off), a note strip with the measured pitch, and how it was made.
- **Instruments**: one card per family, with an honest badge. A measure is *inside* when our gap to
  the recording is no larger than the spread between real recordings of the same note. *Close*
  means every measure is inside, *mixed* at least half, and *off* less than half. The badge shows
  the weaker of the loose and same-note rulers.
- Family pages:
  - **Hear it**: the family's phrases.
  - **Compare**: the simulation and the recording play in lock-step, and switching between them
    (keys 1 and 2) keeps the playhead. It has a blind test, the grader's spectrogram and scalogram
    maps with a difference map, and the provenance of both clips.
  - **Grades**, **History**, **Qualification**, **Timeline** and **References**.
- **How we grade**: the gate ladder with where each family stands and how one note is graded (a
  real example). It also shows every measure with how often the simulations pass it, the rulers
  checked against a second library, the score over time, and the measured pitch of every note
  each family can play.
- **Decisions**: the owner decisions, rulings, findings and grading-method changes that steer the
  program, from `decisions.json` (hand-written), plus the owner-decision beads.
- **Activity**: commits and closes per day, and the milestone matrix across families.

## Layout

| Path | What |
|---|---|
| `index.html`, `css/`, `js/` | the app: plain ES modules, no build step, no external libraries |
| `tools/sync.sh` | mirrors the hosts into a fresh `mirror-<UTC>/`, hard-linking files that did not change |
| `decisions.json` | the decision log shown on the Decisions page; edit by hand |
| `tools/build.py` | stdlib Python plus `ffmpeg`: writes `site/{index,process,decisions,phrases}.json`, `site/family/*.json`, and FLAC audio and webp maps |
| `tools/serve.py` | a local server with byte ranges for seeking; the app at `/`, the data at `/data/` |
| `tools/publish.py`, `cloudflare/` | the public build and its Worker config |

Hosts default to `trj-ts` (grading tree and recordings) and `ts1` (lane renders); override them
with `IE_TRJ` and `IE_TS1`. The work-in-progress render paths are listed at the top of `sync.sh`.
Status, grades history, claims and milestones are read from `origin/main` (`build.py --ref`).

## Credits

University of Iowa Electronic Music Studios, Musical Instrument Samples. The kit groove recording
is the Salamander Drumkit by Alexander Holm (CC BY-SA 3.0). The other sources' credits and
licences are shown next to each clip on the References tab.
