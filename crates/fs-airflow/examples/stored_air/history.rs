//! Bounded CSV input for the existing conduction DutyCycle, not a new time law.
//! Durations tile the window; power_scale multiplies the example's 2 W source.
use fs_conduction::{ConductionError, duty::{DutyCycle, DutySegment}};
use fs_couple::iqn_ils::driver::StepInterval;
use std::io::Read;
use std::path::Path;

const MAX_BYTES: u64 = 65_536;
const MAX_ROWS: usize = 256;
const HEADER: &str = "duration_s,power_scale";

pub(super) fn default_pulse() -> DutyCycle {
    DutyCycle::new(vec![DutySegment::constant(1.0, 1.0).unwrap(),
        DutySegment::constant(4.0, 0.0).unwrap()]).expect("fixed valid pulse")
}

pub(super) fn read(path: &Path) -> Result<DutyCycle, String> {
    let mut text = String::new();
    std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?
        .take(MAX_BYTES + 1).read_to_string(&mut text)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text)
}

fn parse(text: &str) -> Result<DutyCycle, String> {
    if text.len() as u64 > MAX_BYTES { return Err("duty CSV exceeds 64 KiB".into()); }
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some(HEADER) { return Err(format!("expected header {HEADER}")); }
    let mut segments = Vec::new();
    for (index, line) in lines.enumerate() {
        if line.trim().is_empty() { continue; }
        if segments.len() == MAX_ROWS { return Err("duty CSV exceeds 256 segments".into()); }
        let mut values = line.split(',').map(str::trim);
        let duration = values.next().ok_or("missing duration")?.parse::<f64>()
            .map_err(|_| format!("invalid duration on row {}", index + 2))?;
        let scale = values.next().ok_or("missing scale")?.parse::<f64>()
            .map_err(|_| format!("invalid scale on row {}", index + 2))?;
        if values.next().is_some() { return Err(format!("extra column on row {}", index + 2)); }
        if !(2.0e6 * scale).is_finite() { return Err("scaled volumetric source is not finite".into()); }
        segments.push(DutySegment::constant(duration, scale).map_err(|e| e.to_string())?);
    }
    let cycle = DutyCycle::new(segments).map_err(|e| e.to_string())?;
    if cycle.boundaries_s().windows(2).any(|p| p[0] >= p[1]) {
        return Err("a duration is too small to advance the accumulated physical time".into());
    }
    if !(2.0 * cycle.energy_scale_seconds()).is_finite() {
        return Err("declared source energy is not finite".into());
    }
    Ok(cycle)
}

pub(super) fn scale(cycle: &DutyCycle, interval: StepInterval) -> Result<f64, ConductionError> {
    let refuse = || ConductionError::Config { parameter: "stored-air duty interval",
        what: "a positive finite interval must lie wholly in one declared duty segment".into() };
    let start = interval.start_s;
    let end = interval.end_s;
    if !(start.is_finite() && end.is_finite() && start >= 0.0 && end > start && end <= cycle.window_s()) {
        return Err(refuse());
    }
    // DutyCycle::scale_at is left-continuous at shared endpoints. A new step
    // starts in the following segment, whereas the step ending there retains
    // the previous held load. Select by interval rather than adding an epsilon.
    let slot = cycle.boundaries_s().partition_point(|&t| t <= start) - 1;
    if end > cycle.boundaries_s()[slot + 1] { return Err(refuse()); }
    Ok(cycle.segments()[slot].start_scale())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn csv_reuses_declared_segments_and_exact_window_energy() {
        let cycle = parse("duration_s,power_scale\n0.25,1\n0.75,0\n0.5,2\n").unwrap();
        assert_eq!(cycle.boundaries_s(), &[0.0, 0.25, 1.0, 1.5]);
        assert_eq!(cycle.energy_scale_seconds(), 1.25);
    }
    #[test]
    fn discontinuity_uses_interval_ownership_without_time_epsilon() {
        let cycle = default_pulse();
        let step = |a, b| StepInterval { start_s: a, end_s: b };
        assert_eq!(scale(&cycle, step(0.5, 1.0)).unwrap(), 1.0);
        assert_eq!(scale(&cycle, step(1.0, 1.5)).unwrap(), 0.0);
        assert!(scale(&cycle, step(0.5, 1.5)).is_err());
        assert!(scale(&cycle, step(5.0, 5.1)).is_err());
        assert!(scale(&cycle, step(1.0, 1.0)).is_err());
    }
    #[test]
    fn malformed_nonfinite_oversized_and_unrepresentable_history_refuses() {
        for text in ["", HEADER, "duration_s,power_scale\n0,1\n", "duration_s,power_scale\n1,-1\n",
            "duration_s,power_scale\nNaN,1\n", "duration_s,power_scale\n1,inf\n",
            "duration_s,power_scale\n1,1,2\n", "duration_s,power_scale\n1e100,1\n1,1\n"] {
            assert!(parse(text).is_err(), "{text}");
        }
        assert!(parse(&format!("{HEADER}\n{}", "1,1\n".repeat(MAX_ROWS + 1))).is_err());
        assert!(parse(&" ".repeat(MAX_BYTES as usize + 1)).is_err());
    }
}
