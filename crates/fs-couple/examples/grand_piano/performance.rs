//! Prepared sample-accurate control stream. Parsing/allocation is cold; dispatch
//! is a monotone cursor over a prevalidated schedule, not a scan per audio sample.
//! note_on values are post-escapement hammer velocity in m/s, not MIDI gain.
//! jack_staccato/jack_legato values are PEAK JACK FORCE IN NEWTONS, with
//! 7/100 ms sin-squared pulses (Chabassier/Durufle JSV 2014 Table 3). The
//! mechanical engine, not this schedule, determines let-off and strike velocity.

use super::engine::Instrument;

pub const HEADER: &str = "sample,event,key,value";
const MAX_EVENTS: usize = 1_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Control {
    NoteOn { key: u8, velocity_m_s: f64 },
    JackOn { key: u8, peak_n: f64, duration_s: f64 },
    NoteOff { key: u8 },
    Sustain(f64),
    Sostenuto(bool),
    UnaCorda(bool),
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Event {
    pub sample: u64,
    pub control: Control,
}
#[derive(Debug)]
pub struct Performance {
    events: Vec<Event>,
    cursor: usize,
}

impl Performance {
    pub fn read(text: &str, keys: &[u8], frames: u64) -> Result<Self, String> {
        let mut events = Vec::new();
        let mut header = false;
        let mut previous = 0;
        for (line, raw) in text.lines().enumerate() {
            let raw = raw.trim();
            if raw.is_empty() || raw.starts_with('#') { continue; }
            if !header {
                if raw != HEADER { return Err(format!("line {}: expected {HEADER}", line + 1)); }
                header = true;
                continue;
            }
            let parse = (|| -> Result<Event, String> {
                let f: Vec<_> = raw.split(',').map(str::trim).collect();
                if f.len() != 4 { return Err("expected four control fields".into()); }
                let sample = f[0].parse::<u64>().map_err(|_| "invalid sample index")?;
                if sample >= frames || sample < previous {
                    return Err("events must be in sample order and inside the render interval".into());
                }
                let key = f[2].parse::<u8>().map_err(|_| "invalid key")?;
                let value = f[3].parse::<f64>().map_err(|_| "invalid control value")?;
                if !value.is_finite() { return Err("nonfinite control".into()); }
                let control = match f[1] {
                    "note_on" if keys.contains(&key) && value > 0.0 && value <= 8.0 =>
                        Control::NoteOn { key, velocity_m_s: value },
                    "jack_staccato" | "jack_legato" if keys.contains(&key) && value > 0.0 && value <= 200.0 =>
                        Control::JackOn { key, peak_n:value,
                            duration_s:if f[1]=="jack_staccato" {0.007}else{0.100} },
                    "note_off" if keys.contains(&key) && value == 0.0 => Control::NoteOff { key },
                    "sustain" if key == 0 && (0.0..=1.0).contains(&value) => Control::Sustain(value),
                    "sostenuto" if key == 0 && (value == 0.0 || value == 1.0) => Control::Sostenuto(value == 1.0),
                    "una_corda" if key == 0 && (value == 0.0 || value == 1.0) => Control::UnaCorda(value == 1.0),
                    _ => return Err("unknown event, missing scale key, or out-of-range control".into()),
                };
                Ok(Event { sample, control })
            })();
            let event = parse.map_err(|e| format!("line {}: {e}", line + 1))?;
            previous = event.sample;
            if events.len() >= MAX_EVENTS { return Err("performance event budget exhausted".into()); }
            events.push(event);
        }
        if !header { return Err("missing performance header".into()); }
        Ok(Self { events, cursor: 0 })
    }

    /// Dispatch BEFORE the sample's mechanical step. Equal-time events retain
    /// file order, including note-off/pedal transitions. A skipped event refuses
    /// rather than applying it at the wrong time or silently discarding it.
    pub fn dispatch(&mut self, sample: u64, piano: &mut Instrument) -> Result<(), String> {
        while let Some(event) = self.events.get(self.cursor).copied() {
            if event.sample > sample { break; }
            if event.sample < sample { return Err("performance cursor skipped a scheduled sample".into()); }
            match event.control {
                Control::NoteOn { key, velocity_m_s } => piano.note_on(key, velocity_m_s),
                Control::JackOn { key, peak_n, duration_s } => piano.jack_on(key, peak_n, duration_s),
                Control::NoteOff { key } => piano.note_off(key),
                Control::Sustain(value) => piano.set_sustain(value),
                Control::Sostenuto(on) => { piano.set_sostenuto(on); Ok(()) }
                Control::UnaCorda(on) => { piano.set_una_corda(on); Ok(()) }
            }.map_err(|e| format!("sample {sample}, event {}: {e}", self.cursor))?;
            self.cursor += 1;
        }
        Ok(())
    }

    /// A demonstration that only plays keys admitted by this scale, so measured
    /// single-note and partial-keyboard studies no longer fail on hardcoded A4.
    pub fn demonstration(keys: &[u8], rate: u32, frames: u64,
        requested: Option<u8>, velocity: Option<f64>) -> Result<Self, String> {
        if rate == 0 || frames == 0 {
            return Err("demonstration needs a positive sample rate and duration".into());
        }
        let key = match requested {
            Some(key) if keys.contains(&key) => key,
            Some(_) => return Err("requested key is absent from the input scale".into()),
            None => *keys.iter().min_by_key(|&&k| k.abs_diff(69)).ok_or("empty scale")?,
        };
        if velocity.is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 8.0) {
            return Err("hammer velocity must be finite and within (0,8] m/s".into());
        }
        let mut text = format!("{HEADER}\n0,sustain,0,1\n");
        let rate = u64::from(rate);
        for (s, v) in [(0, 0.6), (rate, 2.0), (rate * 2, 4.5)] {
            if s < frames { text.push_str(&format!("{s},note_on,{key},{}\n", velocity.unwrap_or(v))); }
            if s + rate / 2 < frames { text.push_str(&format!("{},note_off,{key},0\n", s + rate / 2)); }
        }
        let chord: Vec<u8> = if requested.is_some() { Vec::new() } else {
            [48, 60, 64, 67].into_iter().filter(|k| keys.contains(k)).collect()
        };
        for (s, event) in [(rate * 3, "note_on"), (rate * 4, "note_off")] {
            if s < frames {
                for k in &chord {
                    let v = if event == "note_on" { velocity.unwrap_or(2.5) } else { 0.0 };
                    text.push_str(&format!("{s},{event},{k},{v}\n"));
                }
            }
        }
        for (s, v) in [(rate * 4 + rate / 2, 0.5), (rate * 5, 0.0)] {
            if s < frames { text.push_str(&format!("{s},sustain,0,{v}\n")); }
        }
        Self::read(&text, keys, frames)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_sample_controls_preserve_input_order() {
        let p = Performance::read(&format!("{HEADER}\n0,note_off,60,0\n0,sostenuto,0,1\n0,note_on,60,2.0\n"), &[60], 100).unwrap();
        assert_eq!(p.events[0].control, Control::NoteOff { key: 60 });
        assert_eq!(p.events[1].control, Control::Sostenuto(true));
        assert_eq!(p.events[2].control, Control::NoteOn { key: 60, velocity_m_s: 2.0 });
    }
    #[test]
    fn invalid_controls_and_missing_measurements_refuse() {
        for row in ["0,note_on,61,2", "0,note_on,60,NaN", "0,note_on,60,0",
            "0,sustain,60,1", "0,una_corda,0,0.5", "100,note_off,60,0",
            "9,note_off,60,0\n8,note_on,60,2", "0,unknown,0,0"] {
            assert!(Performance::read(&format!("{HEADER}\n{row}"), &[60], 100).is_err());
        }
    }
    #[test]
    fn measured_single_note_demo_never_plays_missing_keys() {
        for rate in [8_000, 44_100, 48_000, 96_000] {
            let p = Performance::demonstration(&[36], rate, u64::from(rate) * 6, None, None).unwrap();
            assert!(p.events.iter().all(|e| match e.control {
                Control::NoteOn { key, .. } | Control::NoteOff { key } => key == 36,
                _ => true,
            }));
        }
    }
    #[test]
    fn single_key_and_velocity_controls_preserve_the_existing_demo_contract() {
        let p = Performance::demonstration(&[36, 48, 60, 69], 48_000, 288_000, Some(36), Some(1.25)).unwrap();
        let mut strikes = 0;
        for event in &p.events {
            if let Control::NoteOn { key, velocity_m_s } = event.control {
                assert_eq!(key, 36);
                assert_eq!(velocity_m_s, 1.25);
                strikes += 1;
            }
        }
        assert_eq!(strikes, 3);
        assert!(Performance::demonstration(&[36], 48_000, 288_000, Some(69), None).is_err());
        assert!(Performance::demonstration(&[36], 0, 288_000, None, None).is_err());
        assert!(Performance::demonstration(&[36], 48_000, 0, None, None).is_err());
    }
    #[test]
    fn jack_force_is_not_parsed_as_hammer_velocity() {
        let p=Performance::read(&format!("{HEADER}\n0,jack_staccato,27,70\n50,jack_legato,69,30\n"),&[27,69],100).unwrap();
        assert_eq!(p.events[0].control,Control::JackOn{key:27,peak_n:70.0,duration_s:0.007});
        assert_eq!(p.events[1].control,Control::JackOn{key:69,peak_n:30.0,duration_s:0.1});
        for row in ["0,jack_staccato,28,70","0,jack_legato,27,201","0,jack_legato,27,NaN","0,jack_legato,27,0"] {
            assert!(Performance::read(&format!("{HEADER}\n{row}"),&[27],100).is_err());
        }
    }
}
