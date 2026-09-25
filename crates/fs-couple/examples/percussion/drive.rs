//! Authored playing forces, not fitted contact pulses or replacement audio.
//!
//! CSV rows are `time_s,force_n`. Positive force pushes the existing stick
//! toward the head; negative force lifts it. Linear interpolation is integrated
//! over each mechanical interval, so a knot between ticks is not discarded.
//! The first and last forces must be zero, with zero force outside the file.
//! A frankensim-stick-score-v1 file adds tempo, reusable SI force shapes and
//! beat-based strokes/rolls; it compiles before playback, including overlaps.
//! See SCORE.md. The original CSV path and hot-loop arithmetic are unchanged.
//! These are external player inputs. Contact, rebound, head/air/snare storage,
//! radiation and the energy gate remain owned by the existing mechanics.
use super::super::Error;
use fs_couple::render::plate::impact::ImpactError;
use std::io::Read;

#[path = "drive_score.rs"]
mod score;
#[path = "drive_spatial.rs"]
mod spatial;
pub use spatial::SpatialInput;

const MAX_KNOTS: usize = 65_536;
const MAX_BYTES: u64 = 4 * 1024 * 1024;

#[derive(Clone, Copy, Debug)]
struct Knot {
    time_s: f64,
    force_n: f64,
}

#[derive(Debug)]
pub struct Program {
    knots: Vec<Knot>,
}

impl Program {
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len() as u64 > MAX_BYTES {
            return Err("stick-force file exceeds the 4 MiB admission budget".into());
        }
        if score::selected(text) { return score::parse(text); }
        let mut knots: Vec<Knot> = Vec::new();
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() {
                continue;
            }
            let bad = || format!("stick-force line {}: expected finite time_s,force_n with strictly increasing nonnegative times", line + 1);
            let mut fields = row.split(',').map(str::trim);
            let time_s = fields.next().ok_or_else(bad)?.parse::<f64>().map_err(|_| bad())?;
            let force_n = fields.next().ok_or_else(bad)?.parse::<f64>().map_err(|_| bad())?;
            if fields.next().is_some() || !time_s.is_finite() || time_s < 0.0
                || !force_n.is_finite()
                || knots.last().is_some_and(|previous| time_s <= previous.time_s)
            {
                return Err(bad().into());
            }
            if knots.len() == MAX_KNOTS {
                return Err("stick-force file exceeds 65536 knots".into());
            }
            knots.push(Knot { time_s, force_n });
        }
        if knots.len() < 2 || knots[0].force_n != 0.0 || knots.last().unwrap().force_n != 0.0 {
            return Err("stick-force performance needs at least two knots and zero force at both endpoints".into());
        }
        Ok(Self { knots })
    }

    pub fn load(path: &str) -> Result<Self, Error> {
        // Check actual bytes read, rather than trusting metadata for a changing
        // file or allocating an unbounded string before admission.
        let mut text = String::new();
        std::fs::File::open(path)?.take(MAX_BYTES + 1).read_to_string(&mut text)?;
        Self::parse(&text)
    }

    fn admit(&self, dt_s: f64, steps: u64, tip_weight: f64) -> Result<(), Error> {
        if !dt_s.is_finite() || dt_s <= 0.0 || steps == 0 || steps > (1_u64 << 53)
            || !tip_weight.is_finite() || tip_weight <= 0.0
        {
            return Err("stick drive requires a finite positive clock, bounded steps and the physical inverse-root-mass tip weight".into());
        }
        let end = steps as f64 * dt_s;
        if !end.is_finite() || self.knots.last().unwrap().time_s > end {
            return Err("stick-force performance extends beyond the requested mechanical duration".into());
        }
        // Same numerical work envelope as the host. This is a refusal, not a
        // force clamp, material coefficient, or sound-level normalization.
        if self.knots.iter().any(|k| {
            let generalized = k.force_n * tip_weight;
            !generalized.is_finite() || generalized.abs() > super::super::config(steps, dt_s).maximum_generalized_force
        }) {
            return Err("stick-force performance exceeds the generalized-force envelope".into());
        }
        Ok(())
    }

    fn average(&self, begin: f64, end: f64) -> f64 {
        let duration = end - begin;
        let mut i = self.knots.partition_point(|k| k.time_s <= begin).saturating_sub(1);
        let mut mean = 0.0;
        while i + 1 < self.knots.len() {
            let a = self.knots[i];
            let b = self.knots[i + 1];
            if a.time_s >= end {
                break;
            }
            let left = begin.max(a.time_s);
            let right = end.min(b.time_s);
            if right > left {
                let width = b.time_s - a.time_s;
                let x = (left - a.time_s) / width;
                let y = (right - a.time_s) / width;
                let f_left = (1.0 - x) * a.force_n + x * b.force_n;
                let f_right = (1.0 - y) * a.force_n + y * b.force_n;
                // Weighted local trapezoids avoid subtracting two large
                // cumulative impulses late in a long performance.
                mean += (0.5 * f_left + 0.5 * f_right) * ((right - left) / duration);
            }
            i += 1;
        }
        mean
    }
}

/// Remove only this option. The normal playing parser still owns launch speed
/// and strike position; a drive never silently replaces either initial input.
pub fn option(args: &mut Vec<String>) -> Result<Option<Program>, Error> {
    file_option(args, "--stick-force-file")
}

pub fn second_option(args: &mut Vec<String>) -> Result<Option<Program>, Error> {
    file_option(args, "--second-stick-force-file")
}

fn file_option(args: &mut Vec<String>, flag: &str) -> Result<Option<Program>, Error> {
    let mut selected = None;
    let mut i = 0;
    while i < args.len() {
        if args[i] != flag {
            i += 1;
            continue;
        }
        if selected.is_some() {
            return Err(format!("{flag} may be supplied only once").into());
        }
        let path = args.get(i + 1).ok_or_else(|| format!("{flag} needs a force CSV or score path"))?;
        if path.starts_with("--") {
            return Err(format!("{flag} needs a force CSV or score path, not another option").into());
        }
        selected = Some(Program::load(path)?);
        args.drain(i..i + 2);
    }
    Ok(selected)
}

/// A physical tip port supplied by the instrument composition, never an audio
/// channel or an inferred modal address. Each stick or compliant mute jaw has
/// its own physical force program; the established tip_weight field is 1/sqrt(m).
pub struct Input {
    pub program: Program,
    pub coordinate: usize,
    pub tip_weight: f64,
}

/// Persistent, preallocated force staging. An owner refusal must not consume a
/// tick. All inputs share one clock and one mechanical commit; neither hand
/// can consume a scheduled force while the other hand's contact step refuses.
pub struct StickDrive {
    inputs: Vec<Input>,
    spatial_inputs: Vec<SpatialInput>,
    dt_s: f64,
    steps: u64,
    accepted: u64,
    force: Vec<f64>,
}

impl StickDrive {
    pub fn new(program: Program, dt_s: f64, steps: u64, tip_weight: f64, modes: usize) -> Result<Self, Error> {
        Self::new_inputs(vec![Input { program, coordinate: 0, tip_weight }], dt_s, steps, modes)
    }

    pub fn new_inputs(inputs: Vec<Input>, dt_s: f64, steps: u64, modes: usize) -> Result<Self, Error> {
        if modes == 0 || !(1..=4).contains(&inputs.len()) {
            return Err("player drive needs one to four distinct physical inputs (two sticks plus two mute jaws) and a nonempty force basis".into());
        }
        for (i, input) in inputs.iter().enumerate() {
            if input.coordinate >= modes || inputs[..i].iter().any(|p| p.coordinate == input.coordinate) {
                return Err("stick drive ports must be distinct coordinates inside the mechanical force basis".into());
            }
            input.program.admit(dt_s, steps, input.tip_weight)?;
        }
        Ok(Self { inputs, spatial_inputs: Vec::new(), dt_s, steps, accepted: 0, force: vec![0.0; modes] })
    }

    pub fn forces(&mut self, external: &[f64]) -> Result<&[f64], ImpactError> {
        if self.accepted == self.steps {
            return Err(ImpactError::Budget);
        }
        if external.len() != self.force.len() || external.iter().any(|f| !f.is_finite()) {
            return Err(ImpactError::Invalid("stick drive external-force dimensions or finiteness"));
        }
        let begin = self.accepted as f64 * self.dt_s;
        let end = (self.accepted + 1) as f64 * self.dt_s;
        if end <= begin {
            return Err(ImpactError::Invalid("stick drive clock lost representable resolution"));
        }
        self.force.copy_from_slice(external);
        // F/sqrt(m) is work-conjugate to q=sqrt(m)*x. The composition supplies
        // the actual stick addresses; no force enters a head or air directly.
        for input in &self.inputs {
            self.force[input.coordinate] += input.tip_weight * input.program.average(begin, end);
            if !self.force[input.coordinate].is_finite() {
                return Err(ImpactError::Invalid("stick drive force overflow"));
            }
        }
        // One physical hand force acts through its signed geometry-derived
        // row. No per-mode player clocks, normalization, or allocations here.
        for input in &self.spatial_inputs {
            let force = input.program.average(begin, end);
            for (f, weight) in self.force.iter_mut().zip(&input.weights) {
                if *weight != 0.0 { *f += weight * force; }
                if !f.is_finite() {
                    return Err(ImpactError::Invalid("spatial player force overflow"));
                }
            }
        }
        Ok(&self.force)
    }

    pub fn accept(&mut self) {
        self.accepted += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn separate_hand_programs_share_one_retryable_clock_not_one_amplitude() {
        let first = Input { program: Program::parse("0,0\n0.5,4\n1,0").unwrap(),
            coordinate: 0, tip_weight: 2.0 };
        let second = Input { program: Program::parse("0.5,0\n0.75,-8\n1,0").unwrap(),
            coordinate: 2, tip_weight: 0.5 };
        let mut d = StickDrive::new_inputs(vec![first,second],0.5,2,3).unwrap();
        let external = [1.0,3.0,-2.0];
        assert_eq!(d.forces(&external).unwrap(), &[5.0,3.0,-2.0]);
        assert!(d.forces(&[0.0]).is_err());
        assert_eq!(d.forces(&external).unwrap(), &[5.0,3.0,-2.0]);
        d.accept();
        assert_eq!(d.forces(&external).unwrap(), &[5.0,3.0,-4.0]);
        assert_eq!(d.accepted,1);
        d.accept();
        assert!(matches!(d.forces(&external),Err(ImpactError::Budget)));
    }

    #[test]
    fn two_hand_admission_rejects_duplicate_missing_and_out_of_range_ports() {
        let input = |coordinate| Input { program: Program::parse("0,0\n1,0").unwrap(),
            coordinate, tip_weight: 1.0 };
        assert!(StickDrive::new_inputs(vec![],0.5,2,3).is_err());
        assert!(StickDrive::new_inputs(vec![input(0),input(0)],0.5,2,3).is_err());
        assert!(StickDrive::new_inputs(vec![input(0),input(3)],0.5,2,3).is_err());
        assert!(StickDrive::new_inputs(vec![input(0),input(1),input(2),input(3),input(4)],0.5,2,5).is_err());
        assert!(StickDrive::new_inputs(vec![input(0),input(2)],0.5,1,3).is_err());
        // A force program for the second hand alone must not invent a first.
        assert!(StickDrive::new_inputs(vec![input(2)],0.5,2,3).is_ok());
        for flag in ["--stick-force-file","--second-stick-force-file"] {
            assert!(file_option(&mut vec![flag.into()],flag).is_err());
            assert!(file_option(&mut vec![flag.into(),"--strike-speed-m-s".into()],flag).is_err());
        }
    }

    #[test]
    fn integrates_signed_pulses_and_knots_between_mechanical_ticks() {
        let p = Program::parse("# external player force in SI\n0,0\n0.25,4\n0.5,0\n0.75,-4\n1,0\n").unwrap();
        assert!((p.average(0.0, 0.5) - 2.0).abs() < 1e-14);
        assert!((p.average(0.5, 1.0) + 2.0).abs() < 1e-14);
        assert!(p.average(0.0, 1.0).abs() < 1e-14);
        assert!((p.average(0.125, 0.375) - 3.0).abs() < 1e-14);
        assert_eq!(p.average(1.0, 2.0), 0.0);
        let delayed = Program::parse("1,0\n2,2\n3,0").unwrap();
        assert_eq!(delayed.average(0.0, 1.0), 0.0);
        assert!((delayed.average(0.0, 4.0) - 0.5).abs() < 1e-14);
    }

    #[test]
    fn impulse_is_independent_of_tick_partition() {
        let p = Program::parse("0,0\n0.123,7\n0.432,-3\n0.999,0").unwrap();
        let exact = 0.5 * 7.0 * 0.123 + 0.5 * 4.0 * (0.432 - 0.123)
            - 0.5 * 3.0 * (0.999 - 0.432);
        for n in [1, 3, 17, 1000] {
            let mut impulse = 0.0;
            for i in 0..n {
                let a = f64::from(i) / f64::from(n);
                let b = f64::from(i + 1) / f64::from(n);
                impulse += p.average(a, b) * (b - a);
            }
            assert!((impulse - exact).abs() < 1e-12);
        }
    }

    #[test]
    fn refuses_bad_rows_discontinuous_endpoints_and_truncated_performances() {
        for text in ["", "0,0", "0,1\n1,0", "0,0\n1,1", "0,0\n0,0",
            "1,0\n0,0", "-1,0\n0,0", "0,NaN\n1,0", "0,0\ninf,0",
            "0,0,1\n1,0", "0,0\n1", "time_s,force_n\n0,0\n1,0"]
        {
            assert!(Program::parse(text).is_err(), "accepted {text:?}");
        }
        let p = Program::parse("0,0\n1,0").unwrap();
        assert!(p.admit(0.1, 9, 2.0).is_err());
        assert!(p.admit(0.1, 10, 2.0).is_ok());
        assert!(p.admit(0.0, 10, 2.0).is_err());
        assert!(p.admit(0.1, 10, f64::NAN).is_err());
    }

    #[test]
    fn retry_does_not_consume_force_or_modify_other_participants() {
        let p = Program::parse("0,0\n0.5,4\n1,0").unwrap();
        let mut d = StickDrive::new(p, 0.5, 2, 2.0, 3).unwrap();
        let external = [1.0, 3.0, -2.0];
        assert_eq!(d.forces(&external).unwrap(), &[5.0, 3.0, -2.0]);
        assert!(d.forces(&[0.0]).is_err());
        assert_eq!(d.forces(&external).unwrap(), &[5.0, 3.0, -2.0]);
        assert_eq!(d.accepted, 0);
        d.accept();
        assert_eq!(d.forces(&external).unwrap(), &[5.0, 3.0, -2.0]);
        d.accept();
        assert!(matches!(d.forces(&external), Err(ImpactError::Budget)));
    }
}

#[cfg(test)]
#[path = "drive_score_tests.rs"]
mod score_tests;
