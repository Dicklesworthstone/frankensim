//! Bounded Standard MIDI File 0/1 decoding. Timing stays integer/rational.
//!
//! Reference: MIDI Association RP-001 (Standard MIDI Files 1.0), sections
//! 2/3. Four-byte VLQs, per-track running status, conductor tempo and explicit
//! end-of-track are checked. SMPTE, SysEx and unbound channel messages refuse.
use super::{PerformanceError as Error, SmfLimits, cap, invalid};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MidiNote {
    pub track: u16,
    pub channel: u8,
    pub note: u8,
    pub velocity: u8,
    pub on: u64,
    pub off: u64,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MidiController {
    pub channel: u8,
    pub controller: u8,
    pub value: u8,
    pub tick: u64,
}
#[derive(Clone, Debug)]
pub struct MidiFile {
    pub ppq: u16,
    pub notes: Vec<MidiNote>,
    pub controllers: Vec<MidiController>,
    pub tempos: Vec<(u64, u32)>,
    pub end_tick: u64,
}

pub(crate) struct Reader<'a> {
    pub(crate) bytes: &'a [u8],
    pub(crate) at: usize,
}
impl<'a> Reader<'a> {
    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8], Error> {
        let end = self
            .at
            .checked_add(n)
            .ok_or_else(|| invalid("byte extent"))?;
        let result = self
            .bytes
            .get(self.at..end)
            .ok_or(Error::Truncated { offset: self.at })?;
        self.at = end;
        Ok(result)
    }
    pub(crate) fn byte(&mut self) -> Result<u8, Error> {
        Ok(self.take(1)?[0])
    }
    pub(crate) fn u16(&mut self) -> Result<u16, Error> {
        Ok(u16::from_be_bytes(
            self.take(2)?.try_into().map_err(|_| invalid("u16"))?,
        ))
    }
    pub(crate) fn u32(&mut self) -> Result<u32, Error> {
        Ok(u32::from_be_bytes(
            self.take(4)?.try_into().map_err(|_| invalid("u32"))?,
        ))
    }
    pub(crate) fn u64(&mut self) -> Result<u64, Error> {
        Ok(u64::from_be_bytes(
            self.take(8)?.try_into().map_err(|_| invalid("u64"))?,
        ))
    }
    fn vlq(&mut self) -> Result<u32, Error> {
        let mut value = 0;
        for _ in 0..4 {
            let b = self.byte()?;
            value = (value << 7) | u32::from(b & 127);
            if b & 128 == 0 {
                return Ok(value);
            }
        }
        Err(invalid("MIDI VLQ exceeds four bytes"))
    }
    fn data(&mut self) -> Result<u8, Error> {
        let b = self.byte()?;
        if b >= 128 {
            Err(invalid("MIDI data byte has status bit"))
        } else {
            Ok(b)
        }
    }
}

/// Decode without guessing instrument programs, controller semantics or tempo.
pub fn decode(bytes: &[u8], limits: SmfLimits) -> Result<MidiFile, Error> {
    cap(bytes.len(), limits.max_input_bytes, "MIDI input bytes")?;
    let mut r = Reader { bytes, at: 0 };
    if r.take(4)? != b"MThd" {
        return Err(invalid("MThd"));
    }
    let header_len = r.u32()? as usize;
    if header_len < 6 {
        return Err(invalid("MIDI header length"));
    }
    let format = r.u16()?;
    let tracks = r.u16()?;
    if format > 1 {
        return Err(Error::UnsupportedFormat { format });
    }
    if tracks == 0 || (format == 0 && tracks != 1) {
        return Err(invalid("MIDI format/track count"));
    }
    cap(usize::from(tracks), limits.max_tracks, "MIDI tracks")?;
    let ppq = r.u16()?;
    if ppq & 0x8000 != 0 {
        return Err(Error::SmpteDivision { division: ppq });
    }
    if ppq == 0 {
        return Err(invalid("zero PPQ division"));
    }
    r.take(header_len - 6)?;
    let mut result = MidiFile {
        ppq,
        notes: Vec::new(),
        controllers: Vec::new(),
        tempos: Vec::new(),
        end_tick: 0,
    };
    let mut events = 0usize;
    for track in 0..tracks {
        // Alien chunks are length-delimited by RP-001 and have no MIDI events.
        let data = loop {
            let tag = r.take(4)?;
            let len = r.u32()? as usize;
            let data = r.take(len)?;
            if tag == b"MTrk" {
                break data;
            }
        };
        let mut t = Reader { bytes: data, at: 0 };
        let mut tick = 0u64;
        let mut running = None;
        let mut held = BTreeMap::<(u8, u8), (u64, u8)>::new();
        let mut ended = false;
        while t.at < data.len() {
            events = events
                .checked_add(1)
                .ok_or_else(|| invalid("event count overflow"))?;
            cap(events, limits.max_events, "MIDI events")?;
            tick = tick
                .checked_add(u64::from(t.vlq()?))
                .ok_or_else(|| invalid("MIDI tick overflow"))?;
            let first = t.byte()?;
            let status = if first < 128 {
                t.at -= 1;
                running.ok_or_else(|| invalid("running status without channel status"))?
            } else {
                first
            };
            if status < 0xf0 {
                running = Some(status);
            } else {
                running = None;
            }
            match status {
                0x80..=0x9f => {
                    let channel = status & 15;
                    let note = t.data()?;
                    let velocity = t.data()?;
                    if status & 0xf0 == 0x90 && velocity != 0 {
                        if held.insert((channel, note), (tick, velocity)).is_some() {
                            return Err(Error::OverlappingNote { channel, note });
                        }
                    } else {
                        let (on, velocity) = held
                            .remove(&(channel, note))
                            .ok_or(Error::UnmatchedNoteOff { channel, note })?;
                        if tick <= on {
                            return Err(Error::ZeroLengthNote { channel, note });
                        }
                        result.notes.push(MidiNote {
                            track,
                            channel,
                            note,
                            velocity,
                            on,
                            off: tick,
                        });
                    }
                }
                0xb0..=0xbf => result.controllers.push(MidiController {
                    channel: status & 15,
                    controller: t.data()?,
                    value: t.data()?,
                    tick,
                }),
                0xff => {
                    let kind = t.byte()?;
                    let len = t.vlq()? as usize;
                    let payload = t.take(len)?;
                    match kind {
                        0x2f => {
                            if len != 0 || t.at != data.len() {
                                return Err(invalid("end-of-track must end its chunk"));
                            }
                            ended = true;
                        }
                        0x51 => {
                            if track != 0 {
                                return Err(Error::TempoOutsideConductor { track });
                            }
                            if len != 3 {
                                return Err(invalid("Set Tempo length"));
                            }
                            let tempo = (u32::from(payload[0]) << 16)
                                | (u32::from(payload[1]) << 8)
                                | u32::from(payload[2]);
                            if tempo == 0 {
                                return Err(invalid("zero Set Tempo"));
                            }
                            result.tempos.push((tick, tempo));
                        }
                        // Descriptive metadata does not issue a player command.
                        0x00..=0x09 | 0x20 | 0x21 | 0x54 | 0x58 | 0x59 => {}
                        _ => return Err(Error::UnsupportedMeta { kind }),
                    }
                }
                _ => return Err(Error::UnsupportedMessage { track, status }),
            }
        }
        if !ended {
            return Err(invalid("missing end-of-track"));
        }
        if let Some((&(channel, note), _)) = held.first_key_value() {
            return Err(Error::UnterminatedNote { channel, note });
        }
        result.end_tick = result.end_tick.max(tick);
    }
    if r.at != bytes.len() {
        return Err(invalid("trailing MIDI chunks/bytes"));
    }
    result
        .notes
        .sort_by_key(|n| (n.on, n.track, n.channel, n.note));
    result.controllers.sort_by_key(|c| c.tick);
    // This importer represents one performer. Polyphonic/legato routing needs
    // an explicitly different ownership rule; it must never become last-wins.
    for pair in result.notes.windows(2) {
        if pair[0].off > pair[1].on {
            return Err(Error::OverlappingNote {
                channel: pair[1].channel,
                note: pair[1].note,
            });
        }
    }
    Ok(result)
}

impl MidiFile {
    /// Ceiling of the COMPLETE elapsed rational; never round each tempo span.
    pub fn control_tick(
        &self,
        tick: u64,
        rate: u32,
        initial_tempo_us: Option<u32>,
        limit: u64,
    ) -> Result<u64, Error> {
        if rate == 0 || self.ppq == 0 || self.ppq & 0x8000 != 0 || initial_tempo_us == Some(0) {
            return Err(invalid("control rate, PPQ division or initial tempo"));
        }
        let mut tempo = initial_tempo_us;
        let mut previous = 0u64;
        let mut numerator = 0u128;
        for &(at, next) in &self.tempos {
            if at < previous || next == 0 {
                return Err(invalid("unordered or zero tempo"));
            }
            if at > tick {
                break;
            }
            if at > previous {
                numerator = numerator
                    .checked_add(
                        u128::from(at - previous) * u128::from(tempo.ok_or(Error::MissingTempo)?),
                    )
                    .ok_or_else(|| invalid("tempo integral overflow"))?;
            }
            previous = at;
            tempo = Some(next);
        }
        let tempo = tempo.ok_or(Error::MissingTempo)?;
        numerator = numerator
            .checked_add(u128::from(tick - previous) * u128::from(tempo))
            .ok_or_else(|| invalid("tempo integral overflow"))?;
        let denominator = u128::from(self.ppq) * 1_000_000;
        let scaled = numerator
            .checked_mul(u128::from(rate))
            .ok_or_else(|| invalid("control tick overflow"))?;
        let rounded = scaled.div_ceil(denominator);
        if rounded > u128::from(limit) {
            return Err(Error::Limit {
                resource: "control ticks",
                cap: limit,
            });
        }
        u64::try_from(rounded).map_err(|_| invalid("control tick overflow"))
    }
}
