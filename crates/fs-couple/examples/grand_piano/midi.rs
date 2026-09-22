//! Bounded Standard MIDI File import into the existing physical performance.
//! No synthesizer, sound bank, frequency retuning, envelope or pressure gain.
//! Format 0/1, PPQN tempo maps and SMPTE clocks; see MIDI.md for the admission.
use super::{Control, Event};
use std::io::Read;

pub const MAX_BYTES: usize = 32 * 1024 * 1024;
const MAX_EVENTS: usize = 1_000_000;
const MAX_TRACKS: u16 = 256;

/// This is an explicit controller-to-mechanics assumption, not a calibration.
#[derive(Debug, Clone, Copy)]
pub struct Mapping {
    /// Zero-based MIDI channel, chosen by the caller (CLI uses 1..=16).
    pub channel: u8,
    /// Velocity byte 127 maps to this POST-ESCAPEMENT hammer speed.
    pub maximum_velocity_m_s: f64,
    /// Opt into CC64/127 as physical pedal travel instead of its on/off switch.
    pub continuous_sustain: bool,
}
impl Default for Mapping {
    fn default() -> Self {
        Self { channel: 0, maximum_velocity_m_s: 4.5, continuous_sustain: false }
    }
}
#[derive(Debug, Default)]
pub struct Report {
    pub tracks: u16,
    pub selected_note_ons: usize,
    pub other_channel_messages: usize,
    pub ignored_channel_messages: usize,
    pub ignored_sysex_events: usize,
    pub end_releases: usize,
    pub end_sample: u64,
}
pub struct Parsed {
    pub events: Vec<Event>,
    pub report: Report,
}

struct Reader<'a> { bytes: &'a [u8], pos: usize }
impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self { Self { bytes, pos: 0 } }
    fn empty(&self) -> bool { self.pos == self.bytes.len() }
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        let end = self.pos.checked_add(n).ok_or("MIDI byte address overflow")?;
        let out = self.bytes.get(self.pos..end)
            .ok_or_else(|| format!("truncated MIDI data at byte {} (need {n})", self.pos))?;
        self.pos = end;
        Ok(out)
    }
    fn byte(&mut self) -> Result<u8, String> { Ok(self.take(1)?[0]) }
    fn data(&mut self) -> Result<u8, String> {
        let b = self.byte()?;
        if b >= 128 { return Err("MIDI channel data byte has its status bit set".into()); }
        Ok(b)
    }
    fn u16(&mut self) -> Result<u16, String> {
        let b = self.take(2)?; Ok(u16::from_be_bytes([b[0], b[1]]))
    }
    fn u32(&mut self) -> Result<u32, String> {
        let b = self.take(4)?; Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn vlq(&mut self) -> Result<u32, String> {
        let mut value = 0;
        for _ in 0..4 {
            let b = self.byte()?;
            value = (value << 7) | u32::from(b & 127);
            if b & 128 == 0 { return Ok(value); }
        }
        Err("MIDI variable-length quantity exceeds four bytes".into())
    }
}

#[derive(Clone, Copy)]
enum Message { Tempo(u32), On(u8, u8), Off(u8), Cc(u8, u8), End }
struct Timed { tick: u64, message: Message }

fn track(bytes: &[u8], index: u16, format: u16, mapping: Mapping,
    raw: &mut Vec<Timed>, report: &mut Report, budget: &mut usize) -> Result<(), String>
{
    let mut input = Reader::new(bytes);
    let mut tick = 0_u64;
    let mut running = None;
    while !input.empty() {
        if *budget == MAX_EVENTS { return Err("MIDI exceeds one million raw events".into()); }
        *budget += 1;
        tick = tick.checked_add(u64::from(input.vlq()?)).ok_or("MIDI tick overflow")?;
        let first = input.byte()?;
        let status = if first < 128 {
            input.pos -= 1; // Reuse this byte as the first channel data byte.
            running.ok_or("MIDI running status has no preceding channel status")?
        } else { first };
        if (0x80..=0xef).contains(&status) {
            running = Some(status);
            let a = input.data()?;
            let family = status & 0xf0;
            let b = if family == 0xc0 || family == 0xd0 { 0 } else { input.data()? };
            if status & 15 != mapping.channel {
                report.other_channel_messages += 1;
                continue;
            }
            let message = match family {
                0x80 => Message::Off(a),
                0x90 if b == 0 => Message::Off(a),
                0x90 => Message::On(a, b),
                0xb0 => Message::Cc(a, b),
                0xe0 if a != 0 || b != 64 => {
                    return Err("selected MIDI channel contains pitch bend; string tension is not a MIDI pitch wheel".into());
                }
                _ => { report.ignored_channel_messages += 1; continue; }
            };
            raw.push(Timed { tick, message });
        } else {
            // Meta and SysEx events cancel running status even when skipped.
            running = None;
            match status {
                0xff => {
                    let kind = input.data()?;
                    let len = usize::try_from(input.vlq()?).map_err(|_| "MIDI meta length overflow")?;
                    let data = input.take(len)?;
                    match kind {
                        0x2f => {
                            if len != 0 || !input.empty() {
                                return Err("MIDI end-of-track must be empty and the final track event".into());
                            }
                            raw.push(Timed { tick, message: Message::End });
                            return Ok(());
                        }
                        0x51 => {
                            if len < 3 || (format == 1 && index != 0) {
                                return Err("MIDI tempo needs three bytes and the conductor track in format 1".into());
                            }
                            let tempo = (u32::from(data[0]) << 16) | (u32::from(data[1]) << 8) | u32::from(data[2]);
                            if tempo == 0 { return Err("MIDI tempo must be positive".into()); }
                            raw.push(Timed { tick, message: Message::Tempo(tempo) });
                        }
                        _ => {} // Length-delimited metadata, not physical input.
                    }
                }
                0xf0 | 0xf7 => {
                    let len = usize::try_from(input.vlq()?).map_err(|_| "MIDI SysEx length overflow")?;
                    input.take(len)?;
                    report.ignored_sysex_events += 1;
                }
                _ => return Err("unsupported system status in MIDI track; expected channel, meta or SysEx event".into()),
            }
        }
    }
    Err("MIDI track is missing its end-of-track event".into())
}

fn push(events: &mut Vec<Event>, sample: u64, control: Control, frames: u64) -> Result<(), String> {
    if events.len() == MAX_EVENTS { return Err("MIDI expands beyond one million physical controls".into()); }
    if sample >= frames {
        return Err(format!("MIDI control at sample {sample} is outside the {frames}-frame render; increase --duration, no silent truncation"));
    }
    events.push(Event { sample, control });
    Ok(())
}

/// Read the actual bytes with a bound before constructing a mechanical image.
pub fn load(path: &str, keys: &[u8], rate: u32, frames: u64, mapping: Mapping) -> Result<Parsed, String> {
    let mut bytes = Vec::new();
    std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?
        .take(MAX_BYTES as u64 + 1).read_to_end(&mut bytes).map_err(|e| format!("{path}: {e}"))?;
    read(&bytes, keys, rate, frames, mapping).map_err(|e| format!("{path}: {e}"))
}

pub fn read(bytes: &[u8], keys: &[u8], rate: u32, frames: u64, mapping: Mapping) -> Result<Parsed, String> {
    if bytes.len() > MAX_BYTES || rate == 0 || frames == 0 || keys.is_empty()
        || mapping.channel > 15 || !mapping.maximum_velocity_m_s.is_finite()
        || mapping.maximum_velocity_m_s <= 0.0 || mapping.maximum_velocity_m_s > 8.0
        || mapping.maximum_velocity_m_s / 127.0 == 0.0
    { return Err("MIDI requires bounded bytes, nonempty scale/clock, channel 1..16 and finite hammer speed in (0,8] m/s".into()); }
    let mut input = Reader::new(bytes);
    if input.take(4)? != b"MThd" { return Err("expected a Standard MIDI File MThd header".into()); }
    let len = usize::try_from(input.u32()?).map_err(|_| "MIDI header length overflow")?;
    let mut header = Reader::new(input.take(len)?);
    let format = header.u16()?;
    let tracks = header.u16()?;
    let division = header.u16()?;
    if format > 1 || tracks == 0 || tracks > MAX_TRACKS || (format == 0 && tracks != 1) {
        return Err("MIDI requires format 0 (one track) or format 1 (1..256 simultaneous tracks); independent format-2 patterns are not merged".into());
    }
    // Exact rational clock: PPQN accumulates tick*microseconds_per_quarter.
    // Timecode accumulates tick*fps_denominator. -29 is 30000/1001, not 29 Hz.
    let (denominator, fixed_tick) = if division & 0x8000 == 0 {
        if division == 0 { return Err("MIDI ticks per quarter must be positive".into()); }
        (u128::from(division) * 1_000_000, None)
    } else {
        let ticks = division & 255;
        let (fps, fraction) = match (division >> 8) as u8 as i8 {
            -24 => (24_u32, 1_u32), -25 => (25, 1), -29 => (30_000, 1_001), -30 => (30, 1),
            _ => return Err("unsupported MIDI SMPTE frame rate".into()),
        };
        if ticks == 0 { return Err("MIDI ticks per SMPTE frame must be positive".into()); }
        (u128::from(fps) * u128::from(ticks), Some(u128::from(fraction)))
    };
    let mut report = Report { tracks, ..Report::default() };
    let mut raw = Vec::new();
    let mut found = 0;
    let mut budget = 0;
    let mut chunks = 0;
    while !input.empty() {
        chunks += 1;
        if chunks > 4096 { return Err("MIDI exceeds 4096 chunks".into()); }
        let kind = input.take(4)?;
        let len = usize::try_from(input.u32()?).map_err(|_| "MIDI chunk length overflow")?;
        let bytes = input.take(len)?;
        if kind == b"MTrk" {
            if found == tracks { return Err("more MIDI tracks than declared in the header".into()); }
            track(bytes, found, format, mapping, &mut raw, &mut report, &mut budget)
                .map_err(|e| format!("MIDI track {}: {e}", found + 1))?;
            found += 1;
        } else if kind == b"MThd" { return Err("duplicate MIDI header".into()); }
        // Unknown chunks are deliberately skipped by their bounded length.
    }
    if found != tracks { return Err("fewer MIDI tracks than declared in the header".into()); }
    // Stable ordering: time, then file track order, then original event order.
    // In particular, never sort a same-track sostenuto edge ahead of a note.
    raw.sort_by_key(|e| e.tick);
    let mut events = Vec::new();
    let mut held = [false; 128];
    let mut sustain = 0.0;
    let mut sostenuto = false;
    let mut una_corda = false;
    let mut tempo = 500_000_u32;
    let mut previous = 0;
    let mut elapsed = 0_u128;
    for event in raw {
        let interval = u128::from(event.tick - previous)
            .checked_mul(fixed_tick.unwrap_or(u128::from(tempo))).ok_or("MIDI time overflow")?;
        elapsed = elapsed.checked_add(interval).ok_or("MIDI time overflow")?;
        previous = event.tick;
        let scaled = elapsed.checked_mul(u128::from(rate)).ok_or("MIDI sample time overflow")?;
        // Ceiling of the ABSOLUTE timestamp: never early, no accumulated
        // rounding drift from tempo segments or from per-event delta times.
        let sample = u64::try_from(scaled / denominator + u128::from(scaled % denominator != 0))
            .map_err(|_| "MIDI sample address overflow")?;
        report.end_sample = sample;
        match event.message {
            Message::Tempo(value) => tempo = value,
            Message::End => {},
            Message::On(key, velocity) => {
                if !keys.contains(&key) { return Err(format!("MIDI key {key} is absent from the supplied physical scale")); }
                if held[usize::from(key)] { return Err(format!("overlapping MIDI note-ons for key {key}; one physical key cannot spawn independent voices")); }
                held[usize::from(key)] = true;
                report.selected_note_ons += 1;
                push(&mut events, sample, Control::NoteOn { key,
                    velocity_m_s: mapping.maximum_velocity_m_s * (f64::from(velocity) / 127.0) }, frames)?;
            }
            Message::Off(key) => {
                if held[usize::from(key)] {
                    held[usize::from(key)] = false;
                    push(&mut events, sample, Control::NoteOff { key }, frames)?;
                }
            }
            Message::Cc(controller, value) => match controller {
                64 => {
                    sustain = if mapping.continuous_sustain { f64::from(value) / 127.0 }
                        else if value >= 64 { 1.0 } else { 0.0 };
                    push(&mut events, sample, Control::Sustain(sustain), frames)?;
                }
                66 => { sostenuto = value >= 64; push(&mut events, sample, Control::Sostenuto(sostenuto), frames)?; }
                67 => { una_corda = value >= 64; push(&mut events, sample, Control::UnaCorda(una_corda), frames)?; }
                121 => {
                    sustain = 0.0; sostenuto = false; una_corda = false;
                    for control in [Control::Sustain(0.0), Control::Sostenuto(false), Control::UnaCorda(false)] {
                        push(&mut events, sample, control, frames)?;
                    }
                }
                123 => {
                    for (key, down) in held.iter_mut().enumerate() {
                        if *down { *down = false; push(&mut events, sample, Control::NoteOff { key: key as u8 }, frames)?; }
                    }
                }
                120 => return Err("MIDI All Sound Off is not a physical state reset; use key/pedal releases and allow ringdown".into()),
                _ => report.ignored_channel_messages += 1,
            },
        }
    }
    if report.selected_note_ons == 0 {
        return Err(format!("no notes on selected MIDI channel {}; select --midi-channel explicitly", mapping.channel + 1));
    }
    if report.end_sample >= frames {
        return Err(format!("MIDI ends at sample {}; increase --duration to include that sample and the desired acoustic tail", report.end_sample));
    }
    // End the host's performance by releasing still-held physical controls,
    // not clearing resonators, muting PCM or overwriting contact histories.
    let end = report.end_sample;
    for (key, down) in held.into_iter().enumerate() {
        if down { push(&mut events, end, Control::NoteOff { key: key as u8 }, frames)?; report.end_releases += 1; }
    }
    for (active, control) in [(sustain != 0.0, Control::Sustain(0.0)),
        (sostenuto, Control::Sostenuto(false)), (una_corda, Control::UnaCorda(false))]
    {
        if active { push(&mut events, end, control, frames)?; report.end_releases += 1; }
    }
    Ok(Parsed { events, report })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn vlq(mut n: u32) -> Vec<u8> {
        let mut out = vec![(n & 127) as u8]; n >>= 7;
        while n != 0 { out.push((n & 127) as u8 | 128); n >>= 7; }
        out.reverse(); out
    }
    fn events(rows: &[(u32, &[u8])]) -> Vec<u8> {
        let mut out = Vec::new();
        for (delta, bytes) in rows { out.extend(vlq(*delta)); out.extend_from_slice(bytes); }
        out
    }
    fn smf(format: u16, division: u16, tracks: &[Vec<u8>]) -> Vec<u8> {
        let mut out = b"MThd\0\0\0\x06".to_vec();
        out.extend(format.to_be_bytes()); out.extend((tracks.len() as u16).to_be_bytes()); out.extend(division.to_be_bytes());
        for t in tracks { out.extend(b"MTrk"); out.extend((t.len() as u32).to_be_bytes()); out.extend(t); }
        out
    }
    fn parse(bytes: &[u8]) -> Result<Parsed, String> { read(bytes, &[60, 64, 69], 48_000, 100_000, Mapping::default()) }
    fn simple() -> Vec<u8> { smf(0, 480, &[events(&[(0, &[0x90, 69, 127]), (480, &[69, 0]), (0, &[0xff, 0x2f, 0])])]) }

    #[test]
    fn tempo_map_merges_tracks_without_delta_rounding_or_losing_order() {
        let tempo = events(&[(0, &[0xff,0x51,3,7,0xa1,0x20]), (480, &[0xff,0x51,3,15,0x42,0x40]), (480, &[0xff,0x2f,0])]);
        let notes = events(&[(0, &[0x90,60,127]), (480, &[60,0]), (0, &[64,64]), (480, &[64,0]), (0, &[0xff,0x2f,0])]);
        let p = parse(&smf(1, 480, &[tempo, notes])).unwrap();
        assert_eq!(p.events.iter().map(|e|e.sample).collect::<Vec<_>>(), [0,24_000,24_000,72_000]);
        assert_eq!(p.events[0].control, Control::NoteOn { key:60, velocity_m_s:4.5 });
        assert_eq!(p.events[1].control, Control::NoteOff { key:60 });
        assert_eq!(p.report.end_sample,72_000);
        assert_eq!(p.report.tracks,2);
        let fractional = smf(0, 7, &[events(&[(0,&[0x90,69,1]), (7,&[69,0]), (0,&[0xff,0x2f,0])])]);
        assert_eq!(parse(&fractional).unwrap().events[1].sample,24_000);
    }
    #[test]
    fn many_fractional_tempo_segments_round_the_absolute_clock_only_once() {
        let mut conductor = Vec::new();
        let mut total = 0_u128;
        for i in 0..1000_u32 {
            let tempo = 500_001 + i;
            conductor.extend(events(&[(u32::from(i != 0),
                &[0xff, 0x51, 3, (tempo >> 16) as u8, (tempo >> 8) as u8, tempo as u8])]));
            total += u128::from(tempo);
        }
        conductor.extend(events(&[(1, &[0xff, 0x2f, 0])]));
        let notes = events(&[(0, &[0x90, 69, 127]), (1000, &[69, 0]), (0, &[0xff, 0x2f, 0])]);
        let bytes = smf(1, 997, &[conductor, notes]);
        for rate in [8000, 44100, 48000, 96000] {
            let p = read(&bytes, &[69], rate, 100_000, Mapping::default()).unwrap();
            let expected = (total * u128::from(rate)).div_ceil(997_000_000) as u64;
            assert_eq!(p.events[1].sample, expected);
            assert_eq!(p.report.end_sample, expected);
        }
    }
    #[test]
    fn length_delimited_extensions_and_sysex_are_skipped_not_played() {
        let track = events(&[(0, &[0xf0, 3, 0x7d, 1, 0xf7]), (0, &[0x90, 69, 127]),
            (1, &[69, 0]), (0, &[0xff, 0x2f, 0])]);
        let base = smf(0, 480, &[track]);
        let mut bytes = b"MThd\0\0\0\x08".to_vec();
        bytes.extend_from_slice(&base[8..14]); bytes.extend([9, 9]);
        bytes.extend_from_slice(b"JUNK\0\0\0\x03abc"); bytes.extend_from_slice(&base[14..]);
        let p = parse(&bytes).unwrap();
        assert_eq!(p.events.len(), 2); assert_eq!(p.report.ignored_sysex_events, 1);
        let mut wrong = simple(); wrong[9] = 2;
        assert!(parse(&wrong).is_err());
        wrong = simple(); wrong[11] = 2;
        assert!(parse(&wrong).is_err());
        let zero_tempo = smf(0, 480, &[events(&[(0, &[0xff, 0x51, 3, 0, 0, 0]),
            (0, &[0x90, 69, 1]), (1, &[0xff, 0x2f, 0])])]);
        assert!(parse(&zero_tempo).is_err());
    }
    #[test]
    fn pedal_edges_keep_file_order_and_half_pedal_is_an_explicit_mapping() {
        let bytes=smf(0,480,&[events(&[(0,&[0x90,69,127]),(0,&[0xb0,66,127]),
            (0,&[64,63]),(0,&[67,127]),(480,&[0x90,69,0]),(0,&[0xff,0x2f,0])])]);
        let p=parse(&bytes).unwrap();
        assert!(matches!(p.events[0].control,Control::NoteOn{..}));
        assert_eq!(p.events[1].control,Control::Sostenuto(true));
        assert_eq!(p.events[2].control,Control::Sustain(0.0));
        assert_eq!(p.report.end_releases,2);
        let p=read(&bytes,&[69],48_000,100_000,Mapping{continuous_sustain:true,..Mapping::default()}).unwrap();
        assert_eq!(p.events[2].control,Control::Sustain(63.0/127.0));
        assert_eq!(p.report.end_releases,3);
    }
    #[test]
    fn nonselected_channels_cannot_move_this_pianos_keys_or_pedals() {
        let bytes=smf(0,480,&[events(&[(0,&[0x99,1,127]),(0,&[0xb9,64,127]),
            (0,&[0x90,69,127]),(480,&[69,0]),(0,&[0xff,0x2f,0])])]);
        let p=parse(&bytes).unwrap();
        assert_eq!(p.events.len(),2); assert_eq!(p.report.other_channel_messages,2);
        assert!(read(&simple(),&[69],48_000,100_000,Mapping{channel:1,..Mapping::default()}).is_err());
    }
    #[test]
    fn smpte_drop_frame_is_exact_and_tempo_does_not_retime_timecode() {
        let bytes=smf(0,0xe364,&[events(&[(0,&[0xff,0x51,3,15,0x42,0x40]),
            (0,&[0x90,69,127]),(3000,&[69,0]),(0,&[0xff,0x2f,0])])]);
        assert_eq!(parse(&bytes).unwrap().events[1].sample,48_048);
        for division in [0xe700,0xe564,0] { assert!(parse(&smf(0,division,&[events(&[(0,&[0xff,0x2f,0])])])).is_err()); }
    }
    #[test]
    fn every_truncation_and_bad_running_status_refuses_before_publication() {
        let bytes=simple();
        for n in 0..bytes.len() { assert!(parse(&bytes[..n]).is_err(),"accepted prefix {n}"); }
        assert!(parse(&bytes).is_ok());
        for reset in [&[0xff,1,0][..],&[0xf0,1,0xf7][..],&[0xf7,0][..]] {
            let bad=smf(0,480,&[events(&[(0,&[0x90,69,127]),(0,reset),(1,&[69,0]),(0,&[0xff,0x2f,0])])]);
            assert!(parse(&bad).is_err());
        }
        assert!(parse(&smf(0,480,&[vec![0x80,0x80,0x80,0x80,0]])).is_err());
        assert!(parse(&smf(0,480,&[events(&[(0,&[0x90,69,0x80])])])).is_err());
    }
    #[test]
    fn no_implicit_key_substitution_voice_reset_or_render_truncation() {
        assert!(read(&simple(),&[60],48_000,100_000,Mapping::default()).is_err());
        assert!(read(&simple(),&[69],48_000,24_000,Mapping::default()).is_err());
        for bad in [&[0x90,69,2][..],&[0xe0,0,65][..],&[0xb0,120,0][..]] {
            let bytes=smf(0,480,&[events(&[(0,&[0x90,69,127]),(1,bad),(1,&[0xff,0x2f,0])])]);
            assert!(parse(&bytes).is_err());
        }
        for velocity in [0.0,-1.0,8.1,f64::NAN,f64::INFINITY] {
            assert!(read(&simple(),&[69],48_000,100_000,Mapping{maximum_velocity_m_s:velocity,..Mapping::default()}).is_err());
        }
    }
    #[test]
    fn end_of_file_releases_controls_and_all_notes_off_respects_pedals() {
        let bytes=smf(0,480,&[events(&[(0,&[0xb0,64,127]),(0,&[0x90,69,127]),
            (240,&[0xb0,123,0]),(240,&[0xff,0x2f,0])])]);
        let p=parse(&bytes).unwrap();
        assert_eq!(p.events[2],Event{sample:12_000,control:Control::NoteOff{key:69}});
        assert_eq!(p.events[3],Event{sample:24_000,control:Control::Sustain(0.0)});
        let held=smf(0,480,&[events(&[(0,&[0x90,69,127]),(480,&[0xff,0x2f,0])])]);
        assert_eq!(parse(&held).unwrap().report.end_releases,1);
        let reset=smf(0,480,&[events(&[(0,&[0x90,69,127]),(240,&[0xb0,121,0]),(240,&[0xff,0x2f,0])])]);
        let p=parse(&reset).unwrap();
        assert_eq!(p.events[1].control,Control::Sustain(0.0));
        assert_eq!(p.events[2].control,Control::Sostenuto(false));
        assert_eq!(p.events[3].control,Control::UnaCorda(false));
    }
}
