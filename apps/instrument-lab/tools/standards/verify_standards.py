"""Independent check of the standards MIDI files written by build-standards.ts.

Decodes every file listed in <public>/midi/standards.json with mido (no code shared with
the writer) and checks it against the catalog and the file contract: SMF format 1, PPQ 960,
four tracks (conductor, Bass ch0, Piano ch1, Drums ch9), one tempo and 4/4, markers,
balanced note-on/off with no same-pitch overlap on a channel, a closing CC64=0 when the
pedal is used, every track ending on the same tick, and the catalog's counts and duration.

    uv run --with mido python apps/instrument-lab/tools/standards/verify_standards.py \
        apps/instrument-lab/public
"""

import json
import sys
from pathlib import Path

import mido


def check(public: Path, entry: dict) -> tuple[dict, list[str]]:
    path = public / entry["file"]
    problems: list[str] = []
    raw = path.read_bytes()
    mid = mido.MidiFile(path, charset="utf-8")  # titles and chord markers are UTF-8
    if mid.type != 1 or mid.ticks_per_beat != 960 or len(mid.tracks) != 4:
        problems.append(f"header type={mid.type} ppq={mid.ticks_per_beat} tracks={len(mid.tracks)}")
    names = [next((m.name for m in t if m.type == "track_name"), None) for t in mid.tracks]
    if names[1:] != ["Bass", "Piano", "Drums"] or names[0] != entry["title"]:
        problems.append(f"track names {names}")

    tempos, sigs, markers, ends = [], [], [], []
    notes = {0: 0, 1: 0, 9: 0}
    pedal_events = 0
    piano_live = 0
    piano_peak = 0
    for index, track in enumerate(mid.tracks):
        tick = 0
        held: dict[tuple[int, int], int] = {}
        last_pedal = None
        for msg in track:
            tick += msg.time
            if msg.is_meta:
                if msg.type == "set_tempo":
                    tempos.append(msg.tempo)
                elif msg.type == "time_signature":
                    sigs.append((msg.numerator, msg.denominator))
                elif msg.type == "marker":
                    markers.append((tick, msg.text))
                elif msg.type == "end_of_track":
                    ends.append(tick)
                continue
            if index == 0:
                problems.append("channel message in conductor")
                continue
            expected = {1: 0, 2: 1, 3: 9}[index]
            if msg.channel != expected:
                problems.append(f"track {index} has channel {msg.channel}")
            if msg.type == "control_change":
                if msg.control == 64:
                    pedal_events += 1
                    last_pedal = msg.value
                continue
            if msg.type == "note_on" and msg.velocity > 0:
                key = (msg.channel, msg.note)
                if key in held:
                    problems.append(f"overlap ch{msg.channel} note {msg.note} at {tick}")
                held[key] = tick
                notes[msg.channel] += 1
                if msg.channel == 1:
                    piano_live += 1
                    piano_peak = max(piano_peak, piano_live)
            elif msg.type == "note_off" or (msg.type == "note_on" and msg.velocity == 0):
                key = (msg.channel, msg.note)
                if key not in held:
                    problems.append(f"orphan off ch{msg.channel} note {msg.note} at {tick}")
                    continue
                if held.pop(key) == tick:
                    problems.append(f"zero-length note ch{msg.channel} {msg.note} at {tick}")
                if msg.channel == 1:
                    piano_live -= 1
        if held:
            problems.append(f"track {index}: {len(held)} notes never released")
        if last_pedal not in (None, 0):
            problems.append("sustain pedal left down")

    if len(tempos) != 1 or sigs != [(4, 4)]:
        problems.append(f"tempo/meter {tempos} {sigs}")
    if len(set(ends)) != 1:
        problems.append(f"tracks end at different ticks {ends}")
    if sorted(markers, key=lambda m: m[0]) != markers:
        problems.append("markers out of order")
    chorus_marks = [m for m in markers if m[1].startswith("Chorus ")]
    if len(chorus_marks) != entry["choruses"]:
        problems.append(f"{len(chorus_marks)} chorus markers for {entry['choruses']} choruses")

    end_tick = ends[0] if ends else 0
    seconds = mido.tick2second(end_tick, 960, tempos[0]) if tempos else 0.0
    if end_tick != entry["bars"] * 3840:
        problems.append(f"length {end_tick} ticks for {entry['bars']} bars")
    if abs(seconds - entry["durationS"]) > 0.01:
        problems.append(f"duration {seconds:.2f}s vs catalog {entry['durationS']}")
    if notes[1] != entry["pianoNotes"] or notes[9] != entry["drumHits"]:
        problems.append(f"counts piano {notes[1]}/{entry['pianoNotes']} drums {notes[9]}/{entry['drumHits']}")
    if (pedal_events > 0) != (entry["groove"] == "ballad-comp@1"):
        problems.append(f"pedal events {pedal_events} on {entry['groove']}")
    if len(raw) > 25 * 1024:
        problems.append(f"{len(raw)} bytes exceeds 25 KB")

    row = {
        "id": entry["id"],
        "tracks": len(mid.tracks),
        "bass": notes[0],
        "piano": notes[1],
        "drums": notes[9],
        "pedal": pedal_events,
        "chords": len(markers) - len(chorus_marks),
        "seconds": seconds,
        "peak": piano_peak,
        "kb": len(raw) / 1024,
    }
    return row, problems


def main() -> int:
    public = Path(sys.argv[1] if len(sys.argv) > 1 else "apps/instrument-lab/public")
    catalog = json.loads((public / "midi" / "standards.json").read_text())
    print(f"source {catalog['source']}  entries {len(catalog['entries'])}  refused {catalog['refused']}")
    print(f"{'id':22} trk  bass piano drums pedal chords  seconds  peak    KB")
    failures = 0
    for entry in catalog["entries"]:
        row, problems = check(public, entry)
        print(
            f"{row['id']:22} {row['tracks']:3} {row['bass']:5} {row['piano']:5} {row['drums']:5} "
            f"{row['pedal']:5} {row['chords']:6} {row['seconds']:8.2f} {row['peak']:5} {row['kb']:5.1f}"
        )
        for problem in problems:
            failures += 1
            print(f"    PROBLEM {problem}")
    print("OK" if failures == 0 else f"{failures} problems")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
