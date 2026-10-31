# Real instrument reference recordings

`references-supplemental.json`, together with its `baseline_ids` from
`references.json`, supplies actual recordings for all 27 families in
`apps/instrument-explorer/tools/build.py`, plus a separate tenor saxophone,
splash cymbal and all nine drum-kit comparison states. Originals remain
outside Git. The existing Iowa manifest has
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

# Independently identified tenor saxophone, not an alto substitute.
bash data/instruments/fetch-references.sh \
  --manifest data/instruments/references-supplemental.json \
  --id vcsl-tenor-saxophone-c3-nv-vl2
```

`--instrument` is case-insensitive and exact; both selectors are repeatable.
Combining `--instrument` and `--id` intersects them. Unknown selectors and empty
intersections fail before downloading. No selector chooses every entry.
The supplement's complete listening archives include the 1.61 GB TU-Note ZIP
and the 388 MB drum-kit archive. The separate 22 GB force experiment is
range-only: it refuses a complete download and requires `--members --zip-ranges`
as described below.

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
digests without retaining downloads. ZIP archives use a verified local archive
for random access, as described for TU-Note below. Other files selected with
`--members` follow the ordinary acquisition path.

The 27 Salamander members
supply every reference used by `apps/instrument-lab/tools/kitgrade.py`: a primary
hit, a repeat at the same dynamic and a second dynamic for each of nine kit
states. Point its `SALAMANDER_OH_DIR` argument at the extracted
`Salamander/OH` directory. An archive is a
source container, not an audio file.

### Duduk source availability

On 2026-10-30 the public Freesound preview returned HTTP 403. For listening
coverage, the existing screened manifest also pins the independent
[Sound of Armenian Duduk](https://commons.wikimedia.org/wiki/File:Sound_of_Armenian_Duduk.ogg)
recording, which was downloaded, hash-verified and fully decoded:

```sh
bash data/instruments/fetch-references.sh \
  --manifest data/instruments/references-screened.json \
  --id commons-sound-of-armenian-duduk
```

Its source license is CC BY-SA 4.0. This is a different performance and remains
a family-level listening reference: it is not the Freesound specimen, a
matched note, or a characterized excitation/recording-chain measurement.
Both the original Freesound preview and this Commons recording were available
and matched their existing pins in the later 2026-10-30 acquisition below.

## Coverage

Pitch/dynamic labels below come from the source. Chromatic-range files contain
individual notes and silence and need segmentation before note-level
comparison. Instrument keys such as “D whistle” or “C duduk” do not identify
every performed note.

| Instrument or state | Acquired recording | Source |
|---|---|---|
| Violin | arco and pizzicato A4–B4 on A string; arco G3–B3 on G string | Iowa 2012 |
| Violin | open-string notes at multiple dynamics, legato transition, scales and Massenet phrase; two microphone views | TU-Note 2018 |
| Cello | arco A3–Ab4 on A string; arco and pizzicato C2–B2 on C string | Iowa 2012 |
| Cello | open-string Vivaldi phrase, human player and two robot reproductions, with synchronized motion | mdw 2024 |
| Viola | arco A4–B4 on A string | Iowa 2012 |
| Double bass | arco E1–B1 on E string | Iowa 2012 |
| Guitar | B3 on B string, mf | Iowa |
| Steel-string guitar | E4, HV sample, SFZ key 64 | FreePats / FlameStudios |
| Flute | C5–B5, mf, no vibrato | Iowa |
| Bb clarinet | C4–B4, mf; ambient silence | Iowa |
| Oboe | C4–B4, mf | Iowa |
| Bassoon | C3–B3, mf | Iowa |
| Alto saxophone | C4–B4, mf, no vibrato | Iowa |
| Tenor saxophone | C3 filename label, velocity layer 2, no vibrato, Main microphone | VCSL |
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
| Duduk | independent phrase recording, source identity retained; fallback for unavailable preview | Wikimedia Commons |
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

## Synchronized cello motion and sound

The [mdw cello dataset](https://zenodo.org/records/10696680) pairs the professional
player's recorded sound with the motion capture already used by
`descriptions/players/derive_cello_bow_kinematics.py`. It also supplies two robot
reproductions, with source position controls of 16 mm and 15 mm. All three
original WAVs are 68.008367 seconds, mono, 16-bit, 44.1 kHz; their synchronized
stick, nut and bridge markers were recorded at 240 frames per second.

The supplement links each of the three audio recordings to its original motion
CSV by ID, and retains the original README and score in the same cohort. The
fetcher acquires all eight files and verifies their independent hashes:

```sh
bash data/instruments/fetch-references.sh \
  --manifest data/instruments/references-supplemental.json --instrument Cello
```

The common recording chain names the ROGA RG50 microphone 45 cm in front of the
instrument. The source's CC BY 4.0 license, attribution, readme and byte hashes
are pinned in the same manifest. These paired measurements can challenge
predicted responses to real bow trajectories; they do not supply bow force,
hair contact state, exact specimen geometry or material properties. The
existing kinematic derivation uses stick markers, so its line separation is
not a measured hair-string gap. The recordings remain ungraded, and no audio
has been used to fit an excitation or alter a material constant.

On 2026-10-28, all 66 previously selected audio files were reacquired, checked
against their pins and fully decoded. The three new cello WAVs and their three
motion CSVs were independently hashed on two downloads, and all three WAVs
passed full FFmpeg decode. All **69 audio files (280,543,511 bytes)** and the
motion sources were retained outside Git in the session reference directory;
the manifest does not depend on that transient cache.

## Violin notes, transitions and phrases

[TU-Note Violin](https://doi.org/10.14279/depositonce-6747), by Henrik von Coler,
Jonas Margraf and Paul Schuladen, supplies unchanged recordings of Michiko
Feuerlein, recorded by Philipp Pawlowski in the SIM Berlin anechoic chamber.
The selected 20 WAVs cover open G3, D4, A4 and E5 at mf; A4 at
pp and ff; a D4-to-A3 legato transition without vibrato; plain major ascending
and minor descending scales; and a phrase from Massenet's *Méditation*. Each
performance has both a DPA 4099 and a Brüel & Kjær 4006 microphone recording.
These are two microphone views of the same performance, not independent takes.

Twelve original annotation files accompany the WAVs, including note tables and
attack, sustain, release and transition labels. The source's note table uses
A4 = 443 Hz. The package insert describes 24-bit/96 kHz capture, while all 20
selected distributed WAVs decode as mono PCM16 at 96 kHz. The manifest records
both facts; acquisition does not retune or change the delivered encoding.

```sh
bash data/instruments/fetch-references.sh --members \
  --manifest data/instruments/references-supplemental.json --instrument Violin
```

ZIP extraction caches and verifies the complete 1,605,101,158-byte compressed
archive as `TU-Note_Violin.zip` for random access, then retains only the 32 selected members
(67,807,308 bytes) under `TU-Note`. It does not unpack the complete collection.
An existing verified member is reused. `--pin --members` for ZIP requires the
archive already in the reference directory and reports member digests without
extracting them. The source's CC BY-ND 4.0 license, attribution, original member
paths and all byte hashes are recorded in the manifest.

The combined selection contains **89 unique audio files**, totaling
**348,320,151 bytes and 2,225.17089 seconds**. All matched their SHA-256 and
byte-length pins and passed full FFmpeg decoding during acquisition. This
counts the three mdw cello WAVs once and includes the 20 TU-Note WAVs in addition
to the 66 existing recordings. Archives, motion CSVs and annotation files are
excluded from the audio count. The three cello originals also matched the
source-published MD5 checksums before their SHA-256 pins were added. These
checks establish available, intact
recordings; the TU-Note source does not provide measured bow force, bow speed or
specimen geometry, and no material or excitation parameter was fitted to it.

## Rechecked family coverage (2026-10-30)

The baseline, listening supplement and independent Commons duduk were reacquired
as **71 original audio files, 283,804,824 bytes and 1,962.947266 seconds**.
Every file matched its
SHA-256 and byte-length pins and passed a complete FFmpeg decode. This selection
covers all 27 explorer families, separate alto and tenor saxophones,
splash cymbal and all nine drum-kit comparison states. It includes
the three synchronized mdw cello performances; their three motion CSVs, original
README and score were also retained and verified. The full source tar archives
were checked before their selected members were published. Three truncated Iowa
transfers were preserved and recovered or retried against the original pins;
no checksum was changed to accept an incomplete recording.

This acquisition omitted the 20 TU-Note violin WAVs and their large ZIP;
the three Iowa violin sequences already supply arco and pizzicato coverage.
The larger TU-Note selection above remains separately reproducible. Sensor WAVs
are counted separately from these listening recordings. To fetch the 71-file
listening selection with the existing fetcher and an existing
`FRANKENSIM_REFERENCE_DIR`:

```sh
python3 - <<'PY'
import json, subprocess
from pathlib import Path
supplement = "data/instruments/references-supplemental.json"
cohort = json.loads(Path(supplement).read_text())
selections = [
    ("data/instruments/references.json", cohort["baseline_ids"]),
    (supplement, [f["id"] for f in cohort["files"]
                  if f["id"] not in {"tu-note-violin-2018",
                                      "lampis-2024-bowed-string-transient"}]),
    ("data/instruments/references-screened.json",
     ["commons-sound-of-armenian-duduk"]),
]
for manifest, ids in selections:
    command = ["bash", "data/instruments/fetch-references.sh",
               "--members", "--manifest", manifest]
    for recording in ids:
        command.extend(["--id", recording])
    subprocess.run(command, check=True)
PY
```

## Measured bow/string transient

The [Lampis, Mayer and Chatziioannou experiment](https://zenodo.org/records/10946413)
provides an independent component measurement for testing bow/string dynamics.
The selected original recording is the first sweep point in the 2023-10-18
session: **4.84 seconds, 242,000 frames, 50 kHz, 11 channels, PCM32**. The authors'
`main.m` identifies a Dominant G2 string on a rigid monochord, with a 0.7 m
playing length and 0.055 m bow-bridge distance. It is a component experiment;
there is no cello body or microphone-pressure channel in this WAV.

The original README specifies the channels as follows:

| One-based channels | Measurement | Stored units |
|---|---|---|
| 1–2 | Left and right bridge-force sensors | V |
| 3–4 | Left and right nut-force sensors | V |
| 5–6 | Two bow-force sensors | V |
| 7 | Unused by the source analysis | — |
| 8 | String-tension sensor | V |
| 9–11 | Robot X, Y and Z coordinates | m |

The six selected members comprise this WAV, the README and four original
analysis files documenting calibration and offsets. These total 10,715,565
uncompressed bytes. The fetcher can obtain them with about 3.84 MB of HTTP
ranges, without acquiring the 21,950,506,374-byte ZIP:

```sh
bash data/instruments/fetch-references.sh --members --zip-ranges \
  --manifest data/instruments/references-supplemental.json \
  --id lampis-2024-bowed-string-transient
```

Range mode requires exact HTTP 206 responses, matching source size and byte
ranges, then checks every member's SHA-256 and length before publishing any
new member. Its `archive-ranges-read` event reports transferred bytes and no
archive digest. The manifest retains Zenodo's published MD5 as a source fact;
the whole archive has **not** been independently hashed. The normal whole-file
path refuses this source because it has no verified complete SHA-256 pin.

Raw sensor voltage must not be treated as newtons, and robot position must not
be treated as measured velocity. The retained source code specifies calibration,
offset removal, filtering and trajectory processing; its force-correction helper
still uses earlier directory names, so that mapping needs review before a
calibrated force comparison. The unchanged trace can support independently
specified excitation and response checks once those controls are resolved.
One transient does not validate the complete playability map, friction law,
instrument sound, or material-substitution behavior.

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
Iowa and TU-Note recordings; a shared family name does not establish geometric
equivalence.

For cello, the existing Library of Congress collection supplies back length and
bout widths only. Its `unknowns.tsv` explicitly records the missing outline,
arching and plate thickness. A physical cello shell assembly still needs those
measurements and its internal parts; the acquired recordings cannot supply them.

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
