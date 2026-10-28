# Real instrument reference recordings

`references-supplemental.json`, together with its `baseline_ids` from
`references.json`, supplies actual recordings for all 22 instrument family
cards, plus piano, organ, steel-string guitar, timpani, splash cymbal and the nine
drum-kit comparison states. Originals remain outside Git. The existing Iowa manifest has
1,176 files (about 12.3 GB); its selected baseline requires 21 files
(110,212,120 bytes), including an ambient-noise control.

These recordings test predictions from independently specified geometry,
materials, contact and atmosphere. They are not synthesis inputs or targets for
adjusting physical constants until a spectrum matches. Acquisition is not a
claim that any instrument now sounds realistic.

## Selective acquisition

Set `FRANKENSIM_REFERENCE_DIR` to an existing directory on the reference volume.
The fetcher verifies SHA-256 and byte length, refuses changed files and preserves
the originals. It does not resample, normalize or retune.

```sh
# Listing needs neither a destination nor network access.
bash data/instruments/fetch-references.sh --list --instrument Violin
bash data/instruments/fetch-references.sh --list \
  --manifest data/instruments/references-supplemental.json

# One existing Iowa note sequence.
bash data/instruments/fetch-references.sh \
  --id uiowa-violin2012-violin-arco-mf-sula-a4b4-mono

# One supplemental family, without fetching any archives.
bash data/instruments/fetch-references.sh \
  --manifest data/instruments/references-supplemental.json --instrument Recorder
```

`--instrument` is case-insensitive and exact; both selectors are repeatable.
Combining `--instrument` and `--id` intersects them. Unknown selectors and empty
intersections fail before downloading. No selector means the entire chosen
manifest: about 2.2 GB for the supplement, including the 1.61 GB TU-Note violin
ZIP and the 388 MB drum-kit archive.

Fetch the smaller Iowa selection without copying its IDs:

```sh
python3 - <<'PY'
import json, subprocess
from pathlib import Path
cohort = json.loads(Path("data/instruments/references-supplemental.json").read_text())
args = ["bash", "data/instruments/fetch-references.sh"]
for recording in cohort["baseline_ids"]:
    args.extend(["--id", recording])
subprocess.run(args, check=True)
PY
```

The pan-flute, steel-string guitar and Salamander entries are tar archives. Each
lists the exact selected member names, output paths, SHA-256 and byte length.
Use `--members` to stream these archives and retain only their selected audio:

```sh
bash data/instruments/fetch-references.sh --members \
  --manifest data/instruments/references-supplemental.json \
  --instrument 'Pan pipes' --instrument 'Steel-string guitar' --instrument 'Drum kit'
```

This retains 33,517,346 bytes of original recordings instead of 431,599,967
bytes of source archives. Both the complete archive and each selected member
must match their pins before any new audio file is published. Failed downloads
leave `.part` files for inspection; existing recordings and partial files are
never overwritten. A repeat run verifies the extracted recordings without
downloading the archive again. `--pin --members` reports the archive and member
digests without retaining downloads. Other files selected with `--members`
follow the ordinary acquisition path.

The 27 Salamander members
supply every reference used by `apps/instrument-lab/tools/kitgrade.py`: a primary
hit, a repeat at the same dynamic and a second dynamic for each of nine kit
states. Point its `SALAMANDER_OH_DIR` argument at the extracted
`Salamander/OH` directory. An archive is a
source container, not an audio file.

## Coverage

Pitch/dynamic labels below come from the source. Chromatic-range files contain
individual notes and silence and need segmentation before note-level
comparison. Instrument keys such as “D whistle” or “C duduk” do not identify
every performed note.

| Instrument or state | Acquired recording | Source |
|---|---|---|
| Violin | arco and pizzicato A4–B4 on A string; arco G3–B3 on G string | Iowa 2012 |
| Cello | arco A3–Ab4 on A string; arco and pizzicato C2–B2 on C string | Iowa 2012 |
| Viola | arco A4–B4 on A string | Iowa 2012 |
| Double bass | arco E1–B1 on E string | Iowa 2012 |
| Guitar | B3 on B string, mf | Iowa |
| Steel-string guitar | E4, HV sample, SFZ key 64 | FreePats / FlameStudios |
| Flute | C5–B5, mf, no vibrato | Iowa |
| Bb clarinet | C4–B4, mf; ambient silence | Iowa |
| Oboe | C4–B4, mf | Iowa |
| Bassoon | C3–B3, mf | Iowa |
| Alto saxophone | C4–B4, mf, no vibrato | Iowa |
| Bb trumpet | C4–B4, mf, no vibrato | Iowa |
| French horn | C4–B4, mf | Iowa |
| Tenor trombone | C3–B3, mf | Iowa |
| Tuba | C2–B2, mf | Iowa |
| Marimba | C4–B4, mf, yarn mallet | Iowa |
| Xylophone | C5–B5, mf, hard-rubber mallet | Iowa |
| Vibraphone | sustained C4–B4, mf | Iowa |
| Recorder | baroque alto recorder, sustained C4/E4 labels, Main microphone | VCSL |
| Tin whistle | phrase on metal D whistle with wooden stop | TimRoss / Wikimedia Commons |
| Pan pipes | D4/G4 labels, native untuned recordings | pattingthestar / OpenGameArt |
| Duduk | C-duduk performance by Steen, public MP3 preview | juskiddink / Freesound |
| Voice | Mirko “ah”, file 60 (C4 label) | SpectMorph |
| Piano | Steinway B A4, pp/mf/ff | Iowa |
| Organ flue pipe | Renaissance organ 8-foot register C3, Room microphone | VCSL |
| Timpani / kettle | two Timpani 1 strikes, v2/v4 source layers, summed microphones | VCSL |
| Kick, snare wires on/off, high/low toms, closed/open hi-hat, crash, ride | all 27 overhead-microphone inputs to `kitgrade.py`, including repeat hits and second dynamics | Alexander Holm / Salamander |
| Splash cymbal | normal stick strikes, pp/mf/ff | Iowa 2013 |

On 2026-10-23, the selected existing pins matched, supplemental originals were
hashed, and all **45 retained audio files (1,569.812 seconds)** passed a complete
FFmpeg decode. The Salamander archive was fully hashed while streaming and only
its nine selected members were retained. These checks establish availability
and integrity, not acoustic fidelity.

On 2026-10-24, the complete Salamander archive matched its existing digest and
all **27 kit members** passed full FFmpeg decoding, including the 18 previously
missing repeat-hit and second-dynamic inputs. The three original Iowa splash
AIFFs also passed complete decoding: 24-bit, 44.1 kHz stereo. This adds 21
checksum-pinned recordings to the selection. These 30 files were verified in
memory because the current filesystem was full; the manifest records how to
retrieve them and does not assert that an audio cache is present.

On 2026-10-27, all 39 sources in the selected Iowa baseline and supplement were
downloaded again, including all three complete archive streams. Their 66
selected audio files (262,548,365 bytes) matched their SHA-256 and byte-length
pins and passed full FFmpeg decoding. This covers the 22 family cards and the
additional piano, organ, steel-string guitar and percussion comparisons above.
These checks streamed the actual recordings through memory; they do not claim
a persistent audio cache or improved simulation fidelity. The new `--members`
path was also exercised against the real pan-pipe and steel-string archives.

## Bowed notes, transitions and measured player motion

The supplement now also contains two independent bowed-string recording sets:

- [TU-Note Violin](https://doi.org/10.14279/depositonce-6747), by Henrik von
  Coler, Jonas Margraf and Paul Schuladen, recorded with violinist Michiko
  Feuerlein and engineer Philipp Pawlowski in the SIM Berlin anechoic chamber.
  The selection contains **20 original WAVs**: open G3/D4/A4/E5 at mf, A4 at
  pp/ff, a D4-to-A3 legato transition without vibrato, plain ascending/descending
  scales, and a Massenet *Meditation* phrase, each in both original microphone
  positions. The **12 original label files** retain attack, sustain, release
  and transition boundaries and the source's note tables. The two microphones
  observe the same performances; they are not independent takes. The source
  uses A4 = 443 Hz. Its insert describes 24-bit/96 kHz capture, while all selected
  distributed WAVs were measured as **16-bit/96 kHz mono**. The CC BY-ND 4.0
  recordings and labels remain unchanged.
- [Pàmies-Vilà and colleagues' cello player/robot recordings](https://zenodo.org/records/10696680)
  provide **three synchronized 68-second phrases** on open A3/D3/G2, adapting
  Vivaldi's *Spring*: one professional player and two robot performances. This
  is the source of the already-ingested
  `descriptions/players/cello-open-strings-bow-kinematics-pamies-vila-2024.tsv`.
  The original README specifies a ROGA RG50 microphone 45 cm in front of the
  cello and 240 fps motion capture; the delivered audio is PCM16/44.1 kHz mono.
  The source's CC BY 4.0 license, README and score are retained. The kinematics
  supply motion, **not normal force or a hair-contact measurement**; they do
  not characterize the instrument's materials and geometry.

Fetch the new selections with the existing command:

```sh
bash data/instruments/fetch-references.sh --members \
  --manifest data/instruments/references-supplemental.json --instrument Violin
bash data/instruments/fetch-references.sh \
  --manifest data/instruments/references-supplemental.json --instrument Cello
```

ZIP needs random access to its directory, so `--members` retains the verified
compressed ZIP at `TU-Note_Violin.zip` on the reference volume and extracts only
the 32 pinned members under `TU-Note/`. It does not unpack the complete library
or buffer the archive in memory. `--pin --members` for ZIP requires that local
source archive; it then reports member digests without writing extracted files.
Tar sources retain their existing streaming behavior. A partial earlier
download is preserved and causes refusal instead of being overwritten.

The 23 new recordings and all 66 previously selected recordings were retained,
hash-verified and fully decoded in this session: **89 audio files, 348,320,151
bytes and 2,225.17089 seconds**. The new violin archive's complete SHA-256 and
each selected member's hash are in the supplement. The three cello recordings
also matched the source-published MD5 checksums before SHA-256 pins were added.
Annotations and compressed archives are excluded from the audio counts.

These additions expose onset and transition errors and permit a comparison
with a recorded player motion. They do not authorize changing physical
parameters to fit the recordings, or establish same-specimen validation of the
current simulated bodies.

## Recording limits and attribution

The supplement has `gradable: false` throughout. It extends listening coverage;
an objective comparison still needs specimen, pitch, articulation and
recording-chain controls.

- [Iowa](https://theremin.music.uiowa.edu/MIS.html) documents equipment and
  anechoic recording for many families, including vibraphone. The piano is a
  Steinway B recorded in room 2017 of Voxman, with microphones over bass/treble
  strings; it is not a free-field Steinway D measurement. The
  [splash recordings](https://theremin.music.uiowa.edu/MIScymbals.html) are
  anechoic, but the source lists no cymbal model, diameter, thickness or alloy.
  They do not establish equivalence to the simulator's eight-inch specimen.
- [VCSL](https://github.com/sgossner/VCSL) is CC0 and supplies real recordings,
  but permits tuning, mild EQ and noise reduction. Its recorder, organ and
  timpani files do not establish measured specimen geometry or excitation.
- [FreePats](https://freepats.zenvoid.org/Guitar/steel-acoustic-guitar.html)
  explicitly documents filtering, noise suppression, edits and trimming of the
  Gary Campion / FlameStudios guitar recordings. Its license is GPL-3.0-or-later
  with the archive's FreePats composition exception.
- [The pan flute](https://opengameart.org/content/wip-raw-pan-flute-audio) is
  pattingthestar's CC0 amateur, untuned, one-row instrument.
  [The whistle](https://commons.wikimedia.org/wiki/File:MetalTinWhistleWoodenStop.ogg)
  is TimRoss's public-domain lossy phrase.
  [The duduk](https://freesound.org/people/juskiddink/sounds/180004/) is the
  openly supplied CC BY 4.0 MP3 preview by juskiddink, played by Steen, not the
  login-protected original WAV. These are listening references, not
  characterized specimens.
- [SpectMorph](https://github.com/swesterfeld/spectmorph-instruments) supplies
  a real voice recording under CC0. It does not replace the independent
  VocalSet cohort or establish glottal flow and vocal-tract geometry.
- [The archived Salamander release](https://archive.org/details/SalamanderDrumkit)
  by Alexander Holm states CC BY-SA 3.0. Its garage room and microphone response
  remain, and dynamic labels do not measure strike forces. This license pin
  refers to that release, not a later relicensed one.

Licenses, credit, source pages, exact original URLs and hashes accompany the
manifest entries. No audio is redistributed by these metadata files.

## From mechanics to the listening check

For bowed strings, test independently predicted open-string pitch,
inharmonicity and decay against pizzicato first. With the same string and body
inputs, inspect arco onset, sustained motion and release. Changing the contact
implement must affect that same coupled string/body system. Keep recording
room, microphone, player technique and unmeasured parameters explicit when
interpreting a discrepancy.

Existing geometry sources include `geometry/violin-mclennan2009`,
`geometry/violin-cannone-guarneri-1743-esrf2026` and
`geometry/loc-stradivari-family-bouts`. These are different specimens from the
Iowa recordings; a shared family name does not establish geometric equivalence.

[Woodhouse's friction/rosin treatment](https://euphonics.org/9-6-friction-and-rosin-a-sticky-problem/)
and [thermal-model discussion](https://euphonics.org/9-6-3-exploring-the-thermal-bowing-models/)
explain why rosin state and frictional heating matter and why a velocity-only
friction curve is an approximation. [Yamaha's bow construction guide](https://www.yamaha.com/en/musical_instrument_guide/violin/mechanism/mechanism003.html)
describes roughly 160–180 hairs in a ribbon and demonstrates the effect of rosin.
Neither supplies measured coefficients for the exact recorded bow. Finite
ribbon geometry, contact load distribution, material measurements and reciprocal
body/air response must constrain the model; reference audio alone cannot
identify them.

Across instruments, specified geometry and material laws produce excitation,
coupled motion and radiated pressure; independently acquired recordings then
challenge that prediction. Extra-reference beads remain open where their
characterized-recording requirements are not met.
