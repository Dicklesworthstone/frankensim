# Instrument Explorer

A local web dashboard for FrankenSim's instrument program. For every instrument family it lets you
listen to our renders against real recordings, see how far apart they are, and follow how each
family was built and qualified.

It runs only on your own machine. **No audio is committed to this repository, and none ever
should be.** Several of the reference recordings are not ours to redistribute: University of Iowa
MIS clips are free for the project's grading use only, and others are all-rights-reserved,
CC BY-NC, CC BY-SA or GPL. The audio and maps live under `~/.local/share/instrument-explorer/`.

## Running it

```sh
apps/instrument-explorer/tools/sync.sh          # mirror renders, grades and recordings from the hosts
python3 apps/instrument-explorer/tools/build.py # turn the newest mirror into the site data
python3 apps/instrument-explorer/tools/serve.py --open   # http://localhost:8742
```

`sync.sh` and `build.py` are only needed when the program has moved; `serve.py` alone serves the
last build. Reload the page after a rebuild.

## What is in it

- **Instruments**: one card per family, with an honest badge. A measure is *inside* when our gap to
  the recording is no larger than the spread between real recordings of the same note. *Close*
  means every measure is inside, *mixed* at least half, and *off* less than half. The badge shows
  the weaker of the loose and same-note rulers.
- **Listen**: the simulation and the recording play in lock-step, and switching between them keeps
  the playhead. There is a blind test, the grader's spectrogram and scalogram maps with a
  difference map, and the provenance of both clips. Renders the grader never saw get a spectrogram
  computed in the browser.
- **Grades**: every measure of every graded note, set against the real-versus-real ruler.
- **History**: how each note's grade moved across regradings.
- **Qualification**: the QP1–Q5 milestones with their prerequisite beads, plus evidence charts
  (violin body modes, bowed-string regime sweeps, kit band decays, bar before/after).
- **Timeline** and **References**: commits and bead closes for the family, and the screened
  reference recordings with their credits and licences.
- **Process**: commits and closes per day, and the milestone matrix across families.

## Layout

| Path | What |
|---|---|
| `index.html`, `css/`, `js/` | the app: plain ES modules, no build step, no external libraries |
| `tools/sync.sh` | mirrors the hosts into a fresh `mirror-<UTC>/`, hard-linking files that did not change |
| `tools/build.py` | stdlib Python plus `ffmpeg`: writes `site/{index,process}.json`, `site/family/*.json`, and FLAC audio and webp maps |
| `tools/serve.py` | a local server with byte ranges for seeking; the app at `/`, the data at `/data/` |

Hosts default to `trj-ts` (grading tree and recordings) and `ts1` (lane renders); override them
with `IE_TRJ` and `IE_TS1`. The work-in-progress render paths are listed at the top of `sync.sh`.
Status, grades history, claims and milestones are read from `origin/main` (`build.py --ref`).

## Credits

University of Iowa Electronic Music Studios, Musical Instrument Samples. The other sources'
credits and licences are shown next to each clip on the References tab.
