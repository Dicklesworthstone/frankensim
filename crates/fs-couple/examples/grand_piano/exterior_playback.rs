//! Physical preparation for finite-body piano playback. No second engine.
//! All controls are admitted before structural/BEM work and before score input
//! can change state. Supplied cards cover the complete scale, including silence.
use super::{engine, geometry::Course, hammer_materials, linear, steinway_scale};
use super::exterior_geometry::RATE;
use std::collections::BTreeSet;
use super::performance::midi;
#[path = "exterior_score.rs"]
mod score;
pub use score::Score;

#[derive(Clone, Debug)]
pub struct Options {
    pub substeps: usize,
    pub modes: usize,
    pub midi: Option<String>,
    pub performance: Option<String>,
    pub midi_mapping: midi::Mapping,
    pub note: Option<u8>,
    pub velocity: Option<f64>,
    pub(super) mapping_explicit: bool,
    pub hammers: Option<String>,
    pub hammer_footprints: Option<String>,
    pub dampers: Option<String>,
}
impl Default for Options {
    fn default() -> Self {
        Self { substeps: 4, modes: 24, midi: None, performance: None,
            midi_mapping: midi::Mapping::default(), note: None, velocity: None,
            mapping_explicit: false, hammers: None,
            hammer_footprints: None, dampers: None }
    }
}
impl Options {
    /// The legacy positional MIDI path remains legal, before or after flags.
    /// Unknown, repeated, missing and out-of-budget controls fail before I/O.
    pub fn parse(args: &[String]) -> Result<Self, String> {
        let mut result = Self::default();
        let mut seen = BTreeSet::new();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if !flag.starts_with('-') {
                if flag.is_empty() || result.midi.replace(flag.clone()).is_some() {
                    return Err("exterior playback accepts only one positional MIDI path".into());
                }
                continue;
            }
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate playback option {flag}")); }
            if flag == "--midi-half-pedal" {
                result.midi_mapping.continuous_sustain = true;
                result.mapping_explicit = true;
                continue;
            }
            if !["--modes", "--substeps", "--hammers", "--hammer-footprints", "--dampers",
                "--midi", "--performance", "--midi-channel", "--midi-velocity-max-m-s", "--note", "--velocity"].contains(&flag.as_str()) {
                return Err(format!("unknown exterior playback option {flag}"));
            }
            let value = args.next().filter(|v| !v.is_empty() && !v.starts_with("--"))
                .ok_or_else(|| format!("missing value for {flag}"))?;
            let invalid = || format!("invalid value for {flag}: {value}");
            match flag.as_str() {
                "--modes" => result.modes = value.parse().map_err(|_| invalid())?,
                "--substeps" => result.substeps = value.parse().map_err(|_| invalid())?,
                "--hammers" => result.hammers = Some(value.clone()),
                "--hammer-footprints" => result.hammer_footprints = Some(value.clone()),
                "--dampers" => result.dampers = Some(value.clone()),
                "--performance" => result.performance = Some(value.clone()),
                "--midi" => {
                    if result.midi.replace(value.clone()).is_some() { return Err("duplicate MIDI score".into()); }
                }
                "--midi-channel" => {
                    result.midi_mapping.channel = value.parse::<u8>().map_err(|_| invalid())?
                        .checked_sub(1).ok_or_else(invalid)?;
                    result.mapping_explicit = true;
                }
                "--midi-velocity-max-m-s" => {
                    result.midi_mapping.maximum_velocity_m_s = value.parse().map_err(|_| invalid())?;
                    result.mapping_explicit = true;
                }
                "--note" => result.note = Some(value.parse().map_err(|_| invalid())?),
                "--velocity" => result.velocity = Some(value.parse().map_err(|_| invalid())?),
                _ => unreachable!(),
            }
        }
        result.validate()?;
        Ok(result)
    }
    /// Frequency-domain analysis has no hammer or pedal state. It admits only
    /// resolution controls, so material/gesture options cannot be silently ignored.
    pub fn harmonic(args: &[String]) -> Result<Self, String> {
        let options = Self::parse(args)?;
        if options.midi.is_some() || options.performance.is_some() || options.note.is_some()
            || options.velocity.is_some() || options.hammers.is_some()
            || options.hammer_footprints.is_some() || options.dampers.is_some() {
            return Err("response/admittance accept only --modes and --substeps, not playback controls".into());
        }
        Ok(options)
    }
    pub fn validate(&self) -> Result<(), String> {
        if !(1..=linear::MAX_STRING_MODES).contains(&self.modes)
            || !(1..=16).contains(&self.substeps) {
            return Err("exterior playback requires --modes 1..512 and --substeps 1..16".into());
        }
        if self.midi.is_some() && self.performance.is_some() {
            return Err("select exactly one MIDI or CSV performance, not both".into());
        }
        if (self.midi.is_some() || self.performance.is_some()) && (self.note.is_some() || self.velocity.is_some()) {
            return Err("--note/--velocity are demonstration controls, not score overrides".into());
        }
        if self.mapping_explicit && self.midi.is_none() {
            return Err("MIDI mapping controls require a MIDI score".into());
        }
        if self.midi_mapping.channel >= 16 || !self.midi_mapping.maximum_velocity_m_s.is_finite()
            || self.midi_mapping.maximum_velocity_m_s <= 0. || self.midi_mapping.maximum_velocity_m_s > 8.
            || self.note.is_some_and(|k| !(21..=108).contains(&k))
            || self.velocity.is_some_and(|v| !v.is_finite() || v <= 0. || v > 8.) {
            return Err("invalid MIDI channel or finite post-escapement velocity/key range".into());
        }
        Ok(())
    }
    pub fn report(&self, piano: &engine::Instrument) -> String {
        format!("Mechanical rate {} Hz ({} substeps/output frame); string partial ceiling {}, retained {} coordinates including duplex; {} contact sites. Hammer cards: {}; hammer faces: {}; dampers: {}. These are retention/work budgets, not convergence or real-time certificates.",
            piano.bank.rate, self.substeps, self.modes, piano.bank.modes.len(),
            piano.bank.contact_strings.len(), self.hammers.as_deref().unwrap_or("source defaults"),
            self.hammer_footprints.as_deref().unwrap_or("point"),
            self.dampers.as_deref().unwrap_or("point"))
    }
}

/// Complete, already-admitted data. Keeping this separate from file names
/// prevents a failed supplied material/geometry file from selecting defaults.
pub struct Controls {
    materials: Vec<hammer_materials::Material>,
    footprints: Option<linear::hammer_footprint::Specification>,
    dampers: Option<linear::dampers::Specification>,
}
impl Controls {
    pub fn load(options: &Options, courses: &[Course]) -> Result<Self, String> {
        options.validate()?;
        let hammers = options.hammers.as_deref()
            .map(|p| super::read_bounded(p, 1024 * 1024)).transpose()?;
        let footprints = options.hammer_footprints.as_deref()
            .map(|p| super::read_bounded(p, 64 * 1024)).transpose()?;
        let dampers = match options.dampers.as_deref() {
            Some("estimated") => Some(String::from("estimated")),
            Some(path) => Some(super::read_bounded(path, 1024 * 1024)?),
            None => None,
        };
        Self::from_texts(courses, hammers.as_deref(), footprints.as_deref(), dampers.as_deref())
    }
    /// The owners perform geometry, constitutive, duplicate and coverage checks.
    /// 'estimated' is explicit damper selection, never a missing-file fallback.
    pub fn from_texts(courses: &[Course], hammers: Option<&str>, footprints: Option<&str>,
        dampers: Option<&str>) -> Result<Self, String> {
        let keys: Vec<_> = courses.iter().map(|c| c.midi).collect();
        let materials = match hammers {
            Some(text) => hammer_materials::read(text, &keys)?,
            None => courses.iter().map(steinway_scale::hammer_material).collect::<Result<_,_>>()?,
        };
        let footprints = footprints.map(|text|
            linear::hammer_footprint::Specification::read(text, courses)).transpose()?;
        let dampers = match dampers {
            Some("estimated") => Some(linear::dampers::Specification::estimated(courses)?),
            Some(text) => Some(linear::dampers::Specification::read(text, courses)?),
            None => None,
        };
        Ok(Self { materials, footprints, dampers })
    }
    /// Reuse the source shank and original nonlinear felt/bridge engine. Larger
    /// mode budgets never relax its output-frequency ceiling; changing substeps
    /// changes its clock, not its retained acoustic/structural frequency band.
    pub fn instrument(self, courses: Vec<Course>, board: &[linear::BoardMode], options: &Options)
        -> Result<engine::Instrument, String> {
        options.validate()?;
        let mut piano = engine::Instrument::new_with_contact_geometry(courses, board, RATE,
            options.substeps, options.modes, true, self.materials,
            Some(engine::ShankGeometry::published()), self.footprints.as_ref())?;
        if let Some(dampers) = &self.dampers { piano.configure_dampers(dampers)?; }
        Ok(piano)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| String::from(*s)).collect::<Vec<_>>())
    }
    fn course() -> Course { steinway_scale::courses().unwrap()[48] }
    fn felt(stress: f64) -> String {
        format!("frankensim-hammer-materials-v1\nfelt,69,{stress},0.2,2.5,3.2,0.25,0.8,2500000\nbranch,69,2000000,0.0002\n")
    }
    #[test]
    fn legacy_score_and_new_complete_controls_have_unambiguous_admission() {
        let old = parse(&["score.mid"]).unwrap();
        assert_eq!((old.substeps, old.modes), (4, 24));
        let new = parse(&["--modes", "512", "score.mid", "--substeps", "16", "--hammers", "cards.fsh",
            "--hammer-footprints", "faces.fshp", "--dampers", "estimated"]).unwrap();
        assert_eq!((new.substeps, new.modes), (16, 512));
        assert_eq!(new.midi.as_deref(), Some("score.mid"));
        assert!(Options::harmonic(&["--modes".into(),"128".into()]).is_ok());
        assert!(Options::harmonic(&["--dampers".into(),"estimated".into()]).is_err());
        for args in [vec!["--modes", "0"], vec!["--modes", "513"], vec!["--modes", "NaN"],
            vec!["--substeps", "0"], vec!["--substeps", "17"], vec!["--modes"],
            vec!["--modes", "24", "--modes", "48"], vec!["one.mid", "two.mid"],
            vec!["--hammers", "--dampers", "estimated"], vec!["--unknown", "value"]] {
            assert!(parse(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn bass_retention_can_reach_the_output_band_without_retuning_or_aliasing() {
        let course = steinway_scale::courses().unwrap()[0];
        let build = |modes, substeps| {
            let options = Options { modes, substeps, ..Options::default() };
            Controls::from_texts(&[course], None, None, None).unwrap()
                .instrument(vec![course], &super::super::board::demonstration(), &options).unwrap()
        };
        let reference = build(24, 4); let wide = build(512, 4); let fine = build(512, 8);
        let speaking = |p: &engine::Instrument| p.bank.strings.iter().find(|s| s.contact.is_some()).unwrap().modes.clone();
        let a = speaking(&reference); let b = speaking(&wide);
        assert!(b.len() > a.len());
        for (i, j) in a.zip(b) { assert_eq!(reference.bank.modes[i].omega, wide.bank.modes[j].omega); }
        assert_eq!(wide.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>(),
            fine.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>());
        assert!(wide.bank.modes.iter().all(|m| m.omega <= std::f64::consts::TAU * 0.45 * f64::from(RATE)));
        assert_eq!(fine.bank.rate, 8 * RATE);
        assert_eq!(fine.board_trace_len(), 8 * fine.bank.board_count);
    }
    #[test]
    fn supplied_felt_and_faces_change_real_contact_while_the_scale_is_preserved() {
        let c = course(); let options = Options { modes: 48, substeps: 8, ..Options::default() };
        let build = |stress, faces| Controls::from_texts(&[c], Some(&felt(stress)), faces, Some("estimated"))
            .unwrap().instrument(vec![c], &super::super::board::demonstration(), &options).unwrap();
        let faces = "frankensim-hammer-footprints-v1\nspan,69,0.008,2\n";
        let mut soft = build(300000., Some(faces)); let mut hard = build(600000., Some(faces));
        let point = build(300000., None);
        assert_eq!(soft.bank.contact_strings.len(), 2 * point.bank.contact_strings.len());
        assert!(soft.damper_resolution().is_some());
        assert_eq!(soft.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>(),
            hard.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>());
        soft.note_on(69, 1.).unwrap(); hard.note_on(69, 1.).unwrap();
        for _ in 0..2400 { soft.step().unwrap(); hard.step().unwrap(); }
        assert_ne!(soft.bank.q, hard.bank.q);
        for p in [&soft, &hard] {
            assert!(p.accounting.felt_loss_j > 0.);
            assert!(p.accounting.felt_relaxation_loss_j > 0.);
            assert!((p.accounting.input_work_j - p.energy_j() - p.accounting.dissipated_j()).abs() < 1e-7);
        }
    }
    #[test]
    fn incomplete_supplied_physics_never_selects_source_defaults() {
        let c = course(); let missing = "frankensim-hammer-materials-v1\n";
        assert!(Controls::from_texts(&[c], Some(missing), None, None).is_err());
        assert!(Controls::from_texts(&[c], None, Some("frankensim-hammer-footprints-v1\n"), None).is_err());
        assert!(Controls::from_texts(&[c], None, None, Some("frankensim-piano-dampers-v1\n")).is_err());
        let options = Options { hammers: Some("missing-cards.fsh".into()), ..Options::default() };
        assert!(Controls::load(&options, &[c]).is_err());
    }
}
