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
    /// Coefficients integrated together by the shared bonded-section owner.
    /// Effective scale CSVs do not identify axial rigidity and leave this None.
    pub bonded_section: Option<BondedSection>,
    pub strike_fraction: f64,
    pub hammer_mass_kg: f64,
    pub felt_thickness_m: f64,
    /// Total active patch area, NOT an area to replicate on each string.
    pub felt_area_m2: f64,
    /// Zero explicitly disables a duplex segment.
    pub duplex_length_m: f64,
    pub detune_cents: f64,
}

/// Frozen numerical section binding. Private fields prevent an independent EA
/// override from impersonating an integral of the supplied material layers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondedSection {
    linear_density_kg_m: f64,
    flexural_rigidity_nm2: f64,
    axial_rigidity_n: f64,
    bending_loss: Option<BondedBendingLoss>,
    /// Part of the line density carried by unbonded helical wraps (zero for
    /// a bonded rod). Wraps add inertia only; EA, EI and loss are the core's.
    wrap_line_density_kg_m: f64,
}
impl BondedSection {
    pub const fn axial_rigidity_n(self) -> f64 { self.axial_rigidity_n }
    pub const fn bending_loss(self) -> Option<BondedBendingLoss> { self.bending_loss }
    pub const fn wrap_line_density_kg_m(self) -> f64 { self.wrap_line_density_kg_m }
}

/// Frozen bending-loss integral from the SAME bonded section, with the
/// intersection of its source-supported angular-frequency intervals: a
/// Kelvin-Voigt viscous rigidity C_b, or a loss-angle hysteretic rigidity H
/// (the imaginary part of the complex bending stiffness EI + iH).
/// Private fields prevent a detached loss coefficient or support override.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondedBendingLoss {
    viscous_rigidity_nm2s: f64,
    hysteretic_rigidity_nm2: f64,
    omega_band_rad_s: (f64, f64),
}
impl BondedBendingLoss {
    pub const fn viscous_rigidity_nm2s(self) -> f64 { self.viscous_rigidity_nm2s }
    pub const fn hysteretic_rigidity_nm2(self) -> f64 { self.hysteretic_rigidity_nm2 }
    pub const fn omega_band_rad_s(self) -> (f64, f64) { self.omega_band_rad_s }

    /// Rate multiplying relative modal velocity, (C_b + H / omega) k^4 / mu
    /// [1/s]. For a loss angle this is the viscous rate that dissipates the
    /// same energy per cycle at the mode's own frequency, so it is evaluated
    /// per mode. This is bending strain-rate loss; an uncurved moving chord adds none.
    pub fn strain_rate(self, omega: f64, k: f64, mu: f64) -> Result<f64, String> {
        let (lo, hi) = self.omega_band_rad_s;
        if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi {
            return Err(format!("string mode {omega} rad/s outside material bending-loss support [{lo}, {hi}]"));
        }
        if !k.is_finite() || k <= 0.0 || !mu.is_finite() || mu <= 0.0 {
            return Err("invalid material bending-loss wavenumber or line density".into());
        }
        let rate = (self.viscous_rigidity_nm2s / mu + self.hysteretic_rigidity_nm2 / (mu * omega)) * k.powi(4);
        if !rate.is_finite() || rate < 0.0
            || ((self.viscous_rigidity_nm2s > 0.0 || self.hysteretic_rigidity_nm2 > 0.0) && rate == 0.0) {
            return Err("material bending-loss rate is unrepresentable".into());
        }
        Ok(rate)
    }
}

/// A string core together with the source-bound bending loss of the SAME
/// section, under either constitutive law.
#[derive(Debug, Clone, Copy)]
pub enum CoreLoss<'a> {
    KelvinVoigt(&'a fs_material::section::ConcentricKelvinVoigtSection),
    LossAngle(&'a fs_material::section::loss_angle::ConcentricLossAngleSection),
}
impl<'a> CoreLoss<'a> {
    pub fn elastic(self) -> &'a fs_material::section::ConcentricElasticSection {
        match self {
            Self::KelvinVoigt(section) => section.elastic(),
            Self::LossAngle(section) => section.elastic(),
        }
    }
    fn loss(self) -> BondedBendingLoss {
        match self {
            Self::KelvinVoigt(section) => BondedBendingLoss {
                viscous_rigidity_nm2s: section.viscous_stiffness_n_m2_s(),
                hysteretic_rigidity_nm2: 0.0, omega_band_rad_s: section.omega_band_rad_s(),
            },
            Self::LossAngle(section) => BondedBendingLoss {
                viscous_rigidity_nm2s: 0.0,
                hysteretic_rigidity_nm2: section.hysteretic_stiffness_n_m2(),
                omega_band_rad_s: section.omega_band_rad_s(),
            },
        }
    }
}

impl Course {
    /// Replace a speaking string's inertial, bending and axial coefficients with those
    /// of an independently resolved, perfectly bonded concentric material section.
    /// The string PDE and its acoustic/bridge/contact coupling are unchanged:
    /// they consume the updated physical coefficients at their normal call sites.
    ///
    /// This is for solid or *bonded* concentric rods, NOT a helically wound
    /// piano bass string. The latter requires winding pitch/slip mechanics;
    /// its effective EI must not be computed by bonding the wrap to the core.
    /// Tension is an externally imposed prestress and is deliberately retained.
    /// No retuning to a previous nominal pitch hides the material substitution.
    ///
    /// # Errors
    /// Refuses any section/course whose derived frequency or modal mass is
    /// unrepresentable; the original course is not mutated on failure.
    pub fn with_bonded_material_section(
        &self,
        section: &fs_material::section::ConcentricElasticSection,
    ) -> Result<Self, String> {
        if section.inner_radius_m() != 0.0 {
            return Err("a speaking string needs a solid concentric section, not a hollow tube".into());
        }
        let mut changed = *self;
        changed.linear_density_kg_m = section.line_density_kg_m();
        changed.flexural_rigidity_nm2 = section.bending_stiffness_n_m2();
        changed.bonded_section = Some(BondedSection {
            linear_density_kg_m: section.line_density_kg_m(),
            flexural_rigidity_nm2: section.bending_stiffness_n_m2(),
            axial_rigidity_n: section.axial_stiffness_n(),
            bending_loss: None,
            wrap_line_density_kg_m: 0.0,
        });
        changed.validate()?;
        Ok(changed)
    }

    /// Bind a helically wound string: the solid bonded `core` carries the
    /// tension, EA, EI and Kelvin-Voigt bending loss, and `wrapping` adds only
    /// its wire mass. No wrap stiffness, slip friction or wrap loss is
    /// invented. Tension is retained, as for a bonded rod.
    ///
    /// # Errors
    /// Refuses a hollow core, wraps wound on a different core radius, and any
    /// unrepresentable derived course; the original course is not mutated.
    pub fn with_wound_kelvin_voigt_section(
        &self,
        core: &fs_material::section::ConcentricKelvinVoigtSection,
        wrapping: &fs_material::section::wound::HelicalWrapping,
    ) -> Result<Self, String> {
        self.with_loss_section(CoreLoss::KelvinVoigt(core), Some(wrapping))
    }

    /// Bind a bonded core with its own source-bound bending loss (either law),
    /// optionally under helical wraps that add only their wire mass. Tension
    /// is retained. Every other binder is this one with a fixed choice.
    ///
    /// # Errors
    /// Refuses a hollow core, wraps wound on a different core radius, and any
    /// unrepresentable derived course; the original course is not mutated.
    pub fn with_loss_section(
        &self,
        core: CoreLoss<'_>,
        wrapping: Option<&fs_material::section::wound::HelicalWrapping>,
    ) -> Result<Self, String> {
        let elastic = core.elastic();
        let mut changed = self.with_bonded_material_section(elastic)?;
        let section = changed.bonded_section.as_mut().expect("bound elastic section");
        section.bending_loss = Some(core.loss());
        if let Some(wrapping) = wrapping {
            if wrapping.core_outer_radius_m() != elastic.outer_radius_m() {
                return Err("string wraps were wound on a different core radius".into());
            }
            let line_density = elastic.line_density_kg_m() + wrapping.line_density_kg_m();
            section.linear_density_kg_m = line_density;
            section.wrap_line_density_kg_m = wrapping.line_density_kg_m();
            changed.linear_density_kg_m = line_density;
        }
        changed.validate()?;
        Ok(changed)
    }

    /// Bind elastic and dissipative section integrals together. The loss
    /// replaces the demonstration bending law in every retained span/plane.
    /// A later elastic-only substitution clears this loss binding.
    pub fn with_bonded_kelvin_voigt_section(
        &self,
        section: &fs_material::section::ConcentricKelvinVoigtSection,
    ) -> Result<Self, String> {
        self.with_loss_section(CoreLoss::KelvinVoigt(section), None)
    }

    pub fn validate(&self) -> Result<(), String> {
        if let Some(section) = self.bonded_section {
            if self.linear_density_kg_m != section.linear_density_kg_m
                || self.flexural_rigidity_nm2 != section.flexural_rigidity_nm2
                || !section.axial_rigidity_n.is_finite() || section.axial_rigidity_n <= 0.0 {
                return Err("bound string mass, EI and EA must come from the same material section".into());
            }
            if !section.wrap_line_density_kg_m.is_finite() || section.wrap_line_density_kg_m < 0.0
                || section.wrap_line_density_kg_m >= section.linear_density_kg_m {
                return Err("string wrap mass must be a finite nonnegative part of the line density".into());
            }
            if let Some(loss) = section.bending_loss {
                let (lo, hi) = loss.omega_band_rad_s;
                if !loss.viscous_rigidity_nm2s.is_finite() || loss.viscous_rigidity_nm2s < 0.0
                    || !loss.hysteretic_rigidity_nm2.is_finite() || loss.hysteretic_rigidity_nm2 < 0.0
                    || !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi < lo {
                    return Err("invalid bound material bending-loss integral or frequency support".into());
                }
            }
        }
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
            flexural_rigidity_nm2: ei, bonded_section: None, strike_fraction: 0.12,
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
            tension_n: parse(4)?, flexural_rigidity_nm2: parse(5)?, bonded_section: None, strike_fraction: parse(6)?,
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

/// Export the effective linear scale. A material-bound study must also retain
/// its original section file: this twelve-column format carries neither EA nor
/// source claims and cannot reconstruct a bonded material binding by itself.
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
    fn sourced_material_substitution_changes_string_spectrum_without_solver_edits() {
        use fs_evidence::ValidityDomain;
        use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
            PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
        use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
        use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point};
        use fs_qty::{Density, Pressure};

        fn section(density: f64, young: f64, bore: f64) -> ConcentricElasticSection {
            let mut claims = ClaimSet::new();
            for (name, dims, value) in [
                (DENSITY_PROPERTY, Density::DIMS, density),
                (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
            ] {
                claims.insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic string material integration fixture".into(),
                        license: "CC0-1.0".into(), artifact: None,
                    },
                }).unwrap();
            }
            let card = MaterialCard::assemble(
                MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
                    process: "test".into(), revision: 0 },
                claims, vec![],
            ).unwrap();
            let material = resolve_material_state_point(
                &card, &QueryPoint::new(), &[
                    ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
                        ScalarAdmissibility::StrictlyPositive).unwrap(),
                    ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive).unwrap(),
                ], MaterialPropertySelection::SingleClaimOnly,
            ).unwrap();
            ConcentricElasticSection::try_new(bore, vec![ConcentricMaterialLayer {
                outer_radius_m: 0.0005, material,
            }]).unwrap()
        }

        let base = demonstration_scale().unwrap()[48];
        let steel = section(7800.0, 200e9, 0.0);
        let silver_like = section(10500.0, 80e9, 0.0);
        let a = base.with_bonded_material_section(&steel).unwrap();
        let b = base.with_bonded_material_section(&silver_like).unwrap();
        assert_eq!(a.tension_n, base.tension_n);
        assert_eq!(b.tension_n, base.tension_n);
        assert_eq!(a.linear_density_kg_m, steel.line_density_kg_m());
        assert_eq!(b.flexural_rigidity_nm2, silver_like.bending_stiffness_n_m2());
        assert!(b.linear_density_kg_m > a.linear_density_kg_m);
        assert!(b.flexural_rigidity_nm2 < a.flexural_rigidity_nm2);
        assert!(b.partial_hz(1, b.tension_n) < a.partial_hz(1, a.tension_n));
        assert!(b.modal_mass_kg() > a.modal_mass_kg());
        assert!(base.with_bonded_material_section(&section(7800.0, 200e9, 0.0001)).is_err());
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
