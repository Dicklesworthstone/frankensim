//! CSV force gestures and explicit MIDI mappings reach the same performance
//! owner as the normal grand piano. No alternate scheduler or note envelope.
use super::Options;
use super::super::{performance::{midi, Performance}, exterior_geometry::RATE};

pub struct Score {
    pub performance: Performance,
    pub report: String,
}
impl Options {
    /// All file decoding and complete-score admission happen before BEM work or
    /// any excitation. Equal-sample ordering is left with the existing owner.
    pub fn score(&self, keys: &[u8], frames: u64) -> Result<Score, String> {
        self.validate()?;
        if let Some(path) = &self.performance {
            let text = super::super::read_bounded(path, 32 * 1024 * 1024)?;
            return Score::csv(&text, keys, frames);
        }
        if let Some(path) = &self.midi {
            let parsed = midi::load(path, keys, RATE, frames, self.midi_mapping)?;
            return Score::midi(parsed, keys, frames, self.midi_mapping);
        }
        // Keep the old A4-only study when available, but never invent an absent
        // key for a supplied partial scale. An explicit --note still must exist.
        let key = self.note.or_else(|| keys.iter().min_by_key(|&&k| k.abs_diff(69)).copied())
            .ok_or("exterior demonstration needs a nonempty scale")?;
        let velocity = self.velocity.unwrap_or(2.);
        Ok(Score {
            performance: Performance::demonstration(keys, RATE, frames, Some(key), Some(velocity))?,
            report: format!("Demonstration key {key}, post-escapement hammer velocity {velocity} m/s; no force-driven key-action claim."),
        })
    }
}
impl Score {
    pub fn csv(text: &str, keys: &[u8], frames: u64) -> Result<Self, String> {
        if text.len() > 32 * 1024 * 1024 || frames == 0 {
            return Err("CSV performance needs a positive duration and at most 32 MiB".into());
        }
        Ok(Self {
            performance: Performance::read(text, keys, frames)?,
            report: String::from("CSV performance: output-sample clock; note_on in post-escapement m/s; jack_staccato/jack_legato in peak N; pedals use the existing physical controls. Not a full keyboard-action reconstruction."),
        })
    }
    fn midi(parsed: midi::Parsed, keys: &[u8], frames: u64, mapping: midi::Mapping) -> Result<Self, String> {
        let r = &parsed.report;
        let report = format!("MIDI: channel {}; velocity 127 -> {} m/s (uncalibrated); CC64 {}. {} selected note-ons, {} tracks, end at output sample {}; {} other-channel messages, {} ignored channel messages, {} skipped SysEx, {} synthetic end releases. No pitch-wheel retuning or audio-gain mapping.",
            mapping.channel + 1, mapping.maximum_velocity_m_s,
            if mapping.continuous_sustain {"continuous travel"} else {"switch"},
            r.selected_note_ons, r.tracks, r.end_sample, r.other_channel_messages,
            r.ignored_channel_messages, r.ignored_sysex_events, r.end_releases);
        Ok(Self { performance: Performance::from_events(parsed.events, keys, frames)?, report })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::Controls;
    use super::super::super::{board, engine, steinway_scale};
    fn parse(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| String::from(*s)).collect::<Vec<_>>())
    }
    fn piano(key: u8) -> engine::Instrument {
        let c = steinway_scale::courses().unwrap().into_iter().find(|c|c.midi==key).unwrap();
        Controls::from_texts(&[c],None,None,Some("estimated")).unwrap()
            .instrument(vec![c],&board::demonstration(),&Options::default()).unwrap()
    }
    #[test]
    fn explicit_mapping_and_score_selection_cannot_be_silently_ignored() {
        let p = parse(&["--midi", "score.mid", "--midi-channel", "2", "--midi-velocity-max-m-s", "2",
            "--midi-half-pedal"]).unwrap();
        assert_eq!(p.midi_mapping.channel,1);assert_eq!(p.midi_mapping.maximum_velocity_m_s,2.);
        assert!(p.midi_mapping.continuous_sustain);
        for args in [vec!["--midi", "one.mid", "two.mid"], vec!["one.mid", "--midi", "two.mid"],
            vec!["--performance", "p.csv", "one.mid"], vec!["--midi-half-pedal"],
            vec!["one.mid", "--midi-channel", "0"], vec!["one.mid", "--midi-channel", "17"],
            vec!["one.mid", "--midi-velocity-max-m-s", "NaN"], vec!["one.mid", "--midi-velocity-max-m-s", "9"],
            vec!["--performance", "p.csv", "--velocity", "2"], vec!["one.mid", "--note", "69"],
            vec!["--performance", "p.csv", "--midi-channel", "2"], vec!["--velocity", "inf"]] {
            assert!(parse(&args).is_err(),"accepted {args:?}");
        }
        assert!(Options::harmonic(&["--performance".into(),"p.csv".into()]).is_err());
        assert!(Options::harmonic(&["--note".into(),"69".into()]).is_err());
    }
    #[test]
    fn partial_scale_demo_uses_a_real_key_and_explicit_missing_keys_refuse() {
        let mut p=piano(21);let mut s=Options::default().score(&[21],2400).unwrap();
        assert!(s.report.contains("key 21"));s.performance.dispatch(0,&mut p).unwrap();
        assert!(p.accounting.input_work_j>0.);
        assert!(parse(&["--note","69"]).unwrap().score(&[21],2400).is_err());
        assert!(Options::default().score(&[],2400).is_err());
        for text in ["sample,event,key,value\n0,jack_staccato,69,NaN\n",
            "sample,event,key,value\n2400,note_off,69,0\n",
            "sample,event,key,value\n1,note_off,69,0\n0,note_on,69,1\n"] {
            assert!(Score::csv(text,&[69],2400).is_err());
        }
    }
    fn midi_bytes() -> Vec<u8> {
        // 480 ticks/quarter at 500000 us/quarter -> 50 output samples/tick.
        let track = [0,0x90,60,100, 0,0x91,69,127, 0,0xb1,64,63,
            0,0xb1,66,127, 0,0xb1,67,127, 12,0x81,69,0,
            0,0xb1,66,0, 6,0xb1,64,0, 0,0xb1,67,0, 0,0xff,0x2f,0];
        let mut b=b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk".to_vec();
        b.extend((track.len() as u32).to_be_bytes());b.extend(track);b
    }
    #[test]
    fn mapped_midi_drives_the_exact_physical_pedal_and_hammer_trajectory() {
        let options=parse(&["score.mid","--midi-channel","2","--midi-velocity-max-m-s","2",
            "--midi-half-pedal"]).unwrap();
        let bytes=midi_bytes();
        let parsed=midi::read(&bytes,&[69],RATE,2400,options.midi_mapping).unwrap();
        assert_eq!(parsed.report.selected_note_ons,1);assert_eq!(parsed.report.other_channel_messages,1);
        let mut score=Score::midi(parsed,&[69],2400,options.midi_mapping).unwrap();
        assert!(score.report.contains("channel 2"));assert!(score.report.contains("continuous travel"));
        let mut switched=options.midi_mapping;switched.continuous_sustain=false;
        let parsed=midi::read(&bytes,&[69],RATE,2400,switched).unwrap();
        let mut switched_score=Score::midi(parsed,&[69],2400,switched).unwrap();
        let mut actual=piano(69);let mut direct=piano(69);let mut switch=piano(69);
        for n in 0..1200 {
            score.performance.dispatch(n,&mut actual).unwrap();
            switched_score.performance.dispatch(n,&mut switch).unwrap();
            match n {
                0=>{direct.note_on(69,2.).unwrap();direct.set_sustain(63./127.).unwrap();
                    direct.set_sostenuto(true);direct.set_una_corda(true);},
                600=>{direct.note_off(69).unwrap();direct.set_sostenuto(false);},
                900=>{direct.set_sustain(0.).unwrap();direct.set_una_corda(false);},
                _=>{},
            }
            actual.step().unwrap();direct.step().unwrap();switch.step().unwrap();
        }
        assert_eq!(actual.bank.q,direct.bank.q);assert_eq!(actual.bank.v,direct.bank.v);
        assert_eq!(actual.accounting.damper_loss_j,direct.accounting.damper_loss_j);
        assert_ne!(actual.bank.q,switch.bank.q,"half travel must change physical release damping");
        assert!(actual.accounting.felt_loss_j>0.);assert!(actual.accounting.damper_loss_j>0.);
    }
}
