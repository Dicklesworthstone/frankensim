//! SI-valued piano geometry. Manufacturer specifications are NOT a full scale.
//!
//! Source: https://www.steinway.com/pianos/steinway/grand/model-d
//! The D-274 envelope, longest speaking string and board thickness endpoints
//! below are published dimensions. Every other number in `demonstration_scale`
//! is an explicitly authored estimate, not a measurement of a Steinway.

use std::f64::consts::{PI, TAU};
use fs_math::det;

pub const D_LENGTH_M: f64 = 2.74;
pub const D_WIDTH_M: f64 = 1.56;
pub const D_LONGEST_STRING_M: f64 = 2.01;
pub const D_BOARD_CENTER_M: f64 = 0.009;
pub const D_BOARD_EDGE_M: f64 = 0.006;
pub const CSV_HEADER: &str = "midi,unison,length_m,linear_density_kg_m,tension_n,flexural_rigidity_nm2,strike_fraction,hammer_mass_kg,felt_thickness_m,felt_area_m2,duplex_length_m,detune_cents";

/// Per-string work ceiling. Retention stops earlier at the output audio band;
/// a fixed small count would remove most of the bass attack spectrum.
pub const MAX_STRING_MODES: usize = 512;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PartialRetention {
    pub count: usize,
    /// At least one physical partial inside the requested band was omitted.
    pub capped: bool,
}

/// `auto` follows the physical frequency band up to the admitted work ceiling.
/// A numeric selection remains an explicit, potentially narrower approximation.
pub fn parse_partial_budget(value: &str) -> Result<usize, String> {
    let count = if value == "auto" { MAX_STRING_MODES } else {
        value.parse::<usize>().map_err(|_| "string modes require auto or 1..512")?
    };
    if !(1..=MAX_STRING_MODES).contains(&count) {
        return Err("string modes require auto or 1..512".into());
    }
    Ok(count)
}

/// One uniform speaking course; individual unison members have tension offsets.
/// Measured effective mass and bending rigidity allow wound strings to be
/// represented without falsely treating the copper winding as a bonded rod.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Course {
    pub midi: u8,
    pub unison: usize,
    pub length_m: f64,
    pub linear_density_kg_m: f64,
    pub tension_n: f64,
    pub flexural_rigidity_nm2: f64,
    pub strike_fraction: f64,
    pub hammer_mass_kg: f64,
    pub felt_thickness_m: f64,
    /// Total active patch area, NOT an area to replicate on each string.
    pub felt_area_m2: f64,
    /// Zero explicitly disables a duplex segment.
    pub duplex_length_m: f64,
    pub detune_cents: f64,
}

impl Course {
    pub fn validate(&self) -> Result<(), String> {
        if !(21..=108).contains(&self.midi) || !(1..=3).contains(&self.unison) {
            return Err("course needs a piano MIDI key (21..108) and 1..3 strings".into());
        }
        for value in [self.length_m, self.linear_density_kg_m, self.tension_n,
            self.hammer_mass_kg, self.felt_thickness_m, self.felt_area_m2] {
            if !value.is_finite() || value <= 0.0 {
                return Err("course dimensions, mass, density and tension must be finite and positive".into());
            }
        }
        if !self.flexural_rigidity_nm2.is_finite() || self.flexural_rigidity_nm2 < 0.0
            || !self.strike_fraction.is_finite() || !(0.0..1.0).contains(&self.strike_fraction)
            || self.strike_fraction == 0.0 || !self.duplex_length_m.is_finite()
            || self.duplex_length_m < 0.0 || !self.detune_cents.is_finite()
            || self.detune_cents.abs() > 100.0 {
            return Err("invalid rigidity, strike station, duplex length or unison spread".into());
        }
        for n in [1, MAX_STRING_MODES] {
            if !self.partial_hz(n, self.tension_n).is_finite() {
                return Err("derived string frequency overflow".into());
            }
        }
        if !(self.modal_mass_kg() > 0.0 && self.modal_mass_kg().is_finite()) {
            return Err("derived modal mass overflow/underflow".into());
        }
        Ok(())
    }

    /// Pinned stiff-string eigenfrequency, not a frequency-assigned oscillator.
    pub fn partial_hz(&self, n: usize, tension_n: f64) -> f64 {
        let k = n as f64 * PI / self.length_m;
        det::sqrt((tension_n * k * k + self.flexural_rigidity_nm2 * k.powi(4))
            / self.linear_density_kg_m) / TAU
    }

    /// Select this segment's physical partials using its actual length, EI,
    /// mass density and member tension. The caller supplies the OUTPUT band,
    /// not the oversampled mechanical Nyquist frequency. A capped selection
    /// remains explicit rather than being reported as full-band retention.
    pub fn partial_retention(&self, tension_n: f64, band_hz: f64, cap: usize)
        -> Result<PartialRetention, String> {
        self.validate()?;
        if !tension_n.is_finite() || tension_n <= 0.0
            || !band_hz.is_finite() || band_hz <= 0.0
            || !(1..=MAX_STRING_MODES).contains(&cap) {
            return Err("invalid string tension, output band or partial budget".into());
        }
        // Positive tension and nonnegative EI make this sequence increasing.
        // Inspect one more partial to distinguish the work cap from the band.
        for n in 1..=cap + 1 {
            let hz = self.partial_hz(n, tension_n);
            if !hz.is_finite() {
                return Err("derived string frequency overflow".into());
            }
            if hz > band_hz {
                return Ok(PartialRetention { count: n - 1, capped: false });
            }
        }
        Ok(PartialRetention { count: cap, capped: true })
    }

    pub fn modal_mass_kg(&self) -> f64 {
        0.5 * self.linear_density_kg_m * self.length_m
    }

    /// Work-conjugate projection: F -> generalized force and q -> displacement
    /// use the SAME coefficient, in 1/sqrt(kg).
    pub fn strike_shape(&self, n: usize) -> f64 {
        det::sin(n as f64 * PI * self.strike_fraction) / det::sqrt(self.modal_mass_kg())
    }

    /// Retuning changes tensile state, leaving EI fixed. In particular the
    /// bending contribution must NOT be multiplied by the detune factor.
    pub fn tension_at_cents(&self, cents: f64) -> Result<f64, String> {
        if !cents.is_finite() || cents.abs() > 100.0 {
            return Err("tension offset must be finite and within 100 cents".into());
        }
        let bending = self.flexural_rigidity_nm2 * (PI / self.length_m).powi(2);
        let tension = (self.tension_n + bending) * det::exp(2.0 * cents * std::f64::consts::LN_2 / 1200.0) - bending;
        if !(tension > 0.0 && tension.is_finite()) {
            return Err("retuning would require nonpositive or nonfinite tension".into());
        }
        Ok(tension)
    }
}

/// Exact mass of an ideal uniform single-layer circular helix per axial metre.
/// This is a geometry calculation, NOT a copper/steel effective-EI rule.
pub fn helical_linear_density(core_d: f64, wrap_d: f64, pitch: f64,
    core_density: f64, wrap_density: f64) -> Result<f64, String> {
    if [core_d, pitch, core_density, wrap_density].iter().any(|x| !x.is_finite() || *x <= 0.0)
        || !wrap_d.is_finite() || wrap_d < 0.0 {
        return Err("invalid wire/helix geometry or density".into());
    }
    let helix = det::sqrt(1.0 + (PI * (core_d + wrap_d) / pitch).powi(2));
    let mu = PI * 0.25 * (core_density * core_d * core_d + wrap_density * wrap_d * wrap_d * helix);
    if !mu.is_finite() || mu <= 0.0 { return Err("wire mass overflow".into()); }
    Ok(mu)
}

/// An immediately playable, explicitly ESTIMATED 88-key scale. Only the A0
/// length is a manufacturer dimension. Neither intermediate lengths, break
/// points, winding, hammer voicing nor tensions are a Steinway scale drawing.
pub fn demonstration_scale() -> Result<Vec<Course>, String> {
    let mut courses = Vec::with_capacity(88);
    for midi in 21u8..=108 {
        let f = 440.0 * det::exp((f64::from(midi) - 69.0) * std::f64::consts::LN_2 / 12.0);
        let bass = midi <= 40;
        let length = if bass {
            D_LONGEST_STRING_M * det::pow(27.5 / f, 0.30)
        } else {
            D_LONGEST_STRING_M * det::pow(27.5 / 82.406_889_228_217_5, 0.30)
                * det::pow(82.406_889_228_217_5 / f, 0.85)
        };
        let fraction = (f64::from(midi) - 21.0) / 87.0;
        let core_d: f64 = if bass { 0.0012 } else { 0.00105 - 0.00030 * fraction };
        // Authored nominal material constants, NOT specimen measurements.
        let ei = 200.0e9 * PI * core_d.powi(4) / 64.0;
        let k = PI / length;
        let (mu, tension) = if bass {
            // Solve winding geometry for an authored 750 N tension. EI is the
            // unbonded steel-core estimate; measured EI can replace it in CSV.
            let desired = (750.0 * k * k + ei * k.powi(4)) / (TAU * f).powi(2);
            let mut lo = 1.0e-8;
            let mut hi = 0.004;
            for _ in 0..64 {
                let d = 0.5 * (lo + hi);
                if helical_linear_density(core_d, d, d, 7850.0, 8960.0)? < desired { lo = d; } else { hi = d; }
            }
            let d = 0.5 * (lo + hi);
            let mu = helical_linear_density(core_d, d, d, 7850.0, 8960.0)?;
            (mu, (TAU * f).powi(2) * mu / (k * k) - ei * k * k)
        } else {
            let mu = 7850.0 * PI * core_d * core_d / 4.0;
            (mu, (TAU * f).powi(2) * mu / (k * k) - ei * k * k)
        };
        let course = Course {
            midi, unison: if midi < 31 { 1 } else if bass { 2 } else { 3 },
            length_m: length, linear_density_kg_m: mu, tension_n: tension,
            flexural_rigidity_nm2: ei, strike_fraction: 0.12,
            hammer_mass_kg: 0.011 - 0.006 * fraction,
            felt_thickness_m: 0.009 - 0.003 * fraction,
            felt_area_m2: 0.00010, duplex_length_m: if bass { 0.0 } else { length / 3.03 },
            detune_cents: 0.8,
        };
        course.validate()?;
        courses.push(course);
    }
    Ok(courses)
}

/// Strict, unit-explicit input. No inferred dimensions, substituted rows or
/// silently repaired data. Subsets are legal for measured single-note studies.
pub fn read_scale(text: &str) -> Result<Vec<Course>, String> {
    let mut courses = Vec::new();
    let mut seen = [false; 88];
    let mut header = false;
    for (line, raw) in text.lines().enumerate() {
        let raw = raw.trim();
        if raw.is_empty() || raw.starts_with('#') { continue; }
        if !header {
            if raw != CSV_HEADER { return Err(format!("line {}: expected SI scale header", line + 1)); }
            header = true;
            continue;
        }
        let values: Vec<&str> = raw.split(',').map(str::trim).collect();
        if values.len() != 12 { return Err(format!("line {}: expected 12 fields", line + 1)); }
        let parse = |i: usize| values[i].parse::<f64>().map_err(|_| format!("line {}: invalid field {}", line + 1, i + 1));
        let midi = values[0].parse::<u8>().map_err(|_| format!("line {}: invalid MIDI key", line + 1))?;
        let unison = values[1].parse::<usize>().map_err(|_| format!("line {}: invalid unison count", line + 1))?;
        let c = Course { midi, unison, length_m: parse(2)?, linear_density_kg_m: parse(3)?,
            tension_n: parse(4)?, flexural_rigidity_nm2: parse(5)?, strike_fraction: parse(6)?,
            hammer_mass_kg: parse(7)?, felt_thickness_m: parse(8)?, felt_area_m2: parse(9)?,
            duplex_length_m: parse(10)?, detune_cents: parse(11)? };
        c.validate().map_err(|e| format!("line {}: {e}", line + 1))?;
        let index = usize::from(midi - 21);
        if seen[index] { return Err(format!("line {}: duplicate key {midi}", line + 1)); }
        seen[index] = true;
        courses.push(c);
    }
    if courses.is_empty() { return Err("scale contains no courses".into()); }
    courses.sort_by_key(|c| c.midi);
    Ok(courses)
}

pub fn write_scale(courses: &[Course]) -> String {
    let mut text = format!("{CSV_HEADER}\n");
    for c in courses {
        text.push_str(&format!("{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}\n",
            c.midi, c.unison, c.length_m, c.linear_density_kg_m, c.tension_n,
            c.flexural_rigidity_nm2, c.strike_fraction, c.hammer_mass_kg,
            c.felt_thickness_m, c.felt_area_m2, c.duplex_length_m, c.detune_cents));
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    // G0: independent stiff-string frequency inequalities, including the
    // former bass truncation and an intentionally insufficient work budget.
    #[test]
    fn frequency_retention_preserves_in_band_bass_partials_and_bounds_every_register() {
        let scale = demonstration_scale().unwrap();
        for c in &scale {
            for rate in [8_000, 44_100, 48_000, 96_000, 192_000] {
                let band = 0.45 * f64::from(rate);
                for cents in [-c.detune_cents, 0.0, c.detune_cents] {
                    let tension = c.tension_at_cents(cents).unwrap();
                    let kept = c.partial_retention(tension, band, MAX_STRING_MODES).unwrap();
                    assert!(kept.count <= MAX_STRING_MODES);
                    if kept.count != 0 {
                        assert!(c.partial_hz(kept.count, tension) <= band);
                    }
                    let next_in_band = c.partial_hz(kept.count + 1, tension) <= band;
                    assert_eq!(kept.capped, next_in_band);
                    if kept.capped { assert_eq!(kept.count, MAX_STRING_MODES); }
                }
            }
        }
        let bass = scale[0];
        let preview = bass.partial_retention(bass.tension_n, 21_600.0, 24).unwrap();
        let audible = bass.partial_retention(bass.tension_n, 21_600.0, MAX_STRING_MODES).unwrap();
        assert_eq!(preview, PartialRetention { count: 24, capped: true });
        assert!(bass.partial_hz(24, bass.tension_n) < 1_000.0);
        assert!(audible.count > 128 && !audible.capped);
        assert!(bass.partial_hz(audible.count, bass.tension_n) > 21_000.0);
    }

    #[test]
    fn cutoff_is_inclusive_and_duplex_uses_its_own_speaking_length() {
        let c = demonstration_scale().unwrap()[48];
        let cutoff = c.partial_hz(25, c.tension_n);
        assert_eq!(c.partial_retention(c.tension_n, cutoff, 512).unwrap().count, 25);
        assert_eq!(c.partial_retention(c.tension_n, cutoff.next_down(), 512).unwrap().count, 24);
        let below_fundamental = c.partial_hz(1, c.tension_n).next_down();
        assert_eq!(c.partial_retention(c.tension_n, below_fundamental, 512).unwrap(),
            PartialRetention { count: 0, capped: false });
        let duplex = Course { length_m: c.duplex_length_m, ..c };
        let speaking = c.partial_retention(c.tension_n, 21_600.0, 512).unwrap();
        let tail = duplex.partial_retention(c.tension_n, 21_600.0, 512).unwrap();
        assert!(tail.count > 0 && tail.count < speaking.count);
    }

    #[test]
    fn auto_budget_is_bounded_and_never_claims_a_capped_scale_is_full_band() {
        assert_eq!(parse_partial_budget("auto").unwrap(), MAX_STRING_MODES);
        assert_eq!(parse_partial_budget("24").unwrap(), 24);
        for text in ["0", "513", "-1", "NaN", "", "2.5"] {
            assert!(parse_partial_budget(text).is_err());
        }
        let c = Course { flexural_rigidity_nm2: 0.0, ..demonstration_scale().unwrap()[0] };
        let kept = c.partial_retention(c.tension_n, 21_600.0, MAX_STRING_MODES).unwrap();
        assert_eq!(kept, PartialRetention { count: MAX_STRING_MODES, capped: true });
        for (tension, band, cap) in [(0.0, 21_600.0, 24), (f64::NAN, 21_600.0, 24),
            (c.tension_n, 0.0, 24), (c.tension_n, f64::INFINITY, 24),
            (c.tension_n, 21_600.0, 0), (c.tension_n, 21_600.0, 513)] {
            assert!(c.partial_retention(tension, band, cap).is_err());
        }
    }

    #[test]
    fn geometry_roundtrips_all_keys_and_tunes_the_actual_first_partial() {
        let scale = demonstration_scale().unwrap();
        assert_eq!(scale.len(), 88);
        assert_eq!(read_scale(&write_scale(&scale)).unwrap(), scale);
        for c in scale {
            let target = 440.0 * 2.0f64.powf((f64::from(c.midi) - 69.0) / 12.0);
            assert!((c.partial_hz(1, c.tension_n) / target - 1.0).abs() < 1e-12);
            let shifted = c.tension_at_cents(12.0).unwrap();
            assert!((1200.0 * (c.partial_hz(1, shifted) / target).log2() - 12.0).abs() < 1e-9);
        }
    }
    #[test]
    fn mass_normalization_and_winding_are_not_unit_mass_or_solid_copper() {
        let c = demonstration_scale().unwrap()[48];
        assert!((c.strike_shape(1).powi(2) * c.modal_mass_kg() - det::sin(PI * c.strike_fraction).powi(2)).abs() < 1e-14);
        let tight = helical_linear_density(0.001, 0.001, 0.001, 7850.0, 8960.0).unwrap();
        let loose = helical_linear_density(0.001, 0.001, 0.002, 7850.0, 8960.0).unwrap();
        assert!(tight > loose);
    }
    #[test]
    fn malformed_or_duplicate_measurements_refuse() {
        let c = demonstration_scale().unwrap()[0];
        assert!(read_scale(&write_scale(&[c, c])).is_err());
        assert!(read_scale(&write_scale(&[Course { tension_n: f64::NAN, ..c }])).is_err());
        assert!(read_scale(&write_scale(&[Course { strike_fraction: 1.0, ..c }])).is_err());
        assert!(read_scale(CSV_HEADER).is_err());
    }
}
