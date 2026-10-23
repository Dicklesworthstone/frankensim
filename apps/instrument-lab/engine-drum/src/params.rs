//! Declared SI inputs of the live drum and their typed refusals.
//!
//! Every field is a physical or numerical-resolution input. There is no pitch,
//! gain, normalization or "tone" parameter: frequencies come only from the
//! eigen-solve of the heads these numbers describe.

use crate::json::{self, Value};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{ThermochemicalPhase, read_gas_state};
use fs_plate::shell::head::TensionVariation;
use fs_plate::shell::head::tension::TENSION_HARMONIC_COUNT;

/// A typed refusal that crosses the browser boundary as `{code, message}`.
#[derive(Debug, Clone, PartialEq)]
pub struct Refusal {
    /// Stable machine-readable class.
    pub code: &'static str,
    /// Human-readable cause, including the owner's own words where available.
    pub message: String,
}
impl Refusal {
    /// Construct a refusal.
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self { code, message: message.into() }
    }
    /// `{"code":..,"message":..}`.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!("{{\"code\":{},\"message\":{}}}", json::string(self.code), json::string(&self.message))
    }
}
impl core::fmt::Display for Refusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}
impl std::error::Error for Refusal {}

/// One immutable ambient declaration and the shared material owner's result.
/// Both the sealed cavity and exterior receivers consume this exact state.
#[derive(Debug, Clone, PartialEq)]
pub struct Atmosphere {
    declaration: String,
    state: GasState,
}

impl Atmosphere {
    /// Resolve `T P RH` or `mixture T P mole|mass N SPECIES FRACTION ...`.
    /// Temperatures are kelvin, pressure is pascals and fractions are unitless.
    /// The existing material parser owns species, composition and phase limits.
    pub fn from_declaration(text: &str) -> Result<Self, Refusal> {
        if text.len() > 8192 {
            return Err(Refusal::new("gas_refused", "gas declaration exceeds 8 KiB"));
        }
        let mut fields = text.split_ascii_whitespace();
        let state = read_gas_state(&mut fields, None)
            .map_err(|e| Refusal::new("gas_refused", e.to_string()))?;
        if fields.next().is_some() {
            return Err(Refusal::new("gas_refused", "unexpected fields after the gas declaration"));
        }
        Ok(Self { declaration: text.split_ascii_whitespace().collect::<Vec<_>>().join(" "), state })
    }

    /// The same declaration resolved through its source packs, framed exactly
    /// as the live kit's browser packet: each consecutive species, NASA-9 and
    /// transport triple of `lengths` names unchanged canonical bytes. A
    /// `thermochemical ...` record resolves only this way. An empty
    /// `water_species` means none; the phase is `frozen` or `liquid-water`.
    /// Nothing is defaulted.
    pub fn from_source_packet(text: &str, bytes: &[u8], lengths: &[u32], water_species: &str,
        phase: &str) -> Result<Self, Refusal> {
        if text.len() > 8192 {
            return Err(Refusal::new("gas_refused", "gas declaration exceeds 8 KiB"));
        }
        let phase = match phase {
            "frozen" => ThermochemicalPhase::FrozenGas,
            "liquid-water" => ThermochemicalPhase::LiquidWaterLimit,
            _ => return Err(Refusal::new("gas_refused", "gas source phase must be frozen or liquid-water")),
        };
        let water = (!water_species.is_empty()).then_some(water_species);
        let record = fs_material::gas_species::input::packs::read_source_packet(text, bytes, lengths, water, phase)
            .map_err(|e| Refusal::new("gas_refused", e.to_string()))?;
        Ok(Self { declaration: text.split_ascii_whitespace().collect::<Vec<_>>().join(" "), state: *record.state() })
    }

    /// Exact resolved coefficients, without reevaluation or editable overrides.
    #[must_use]
    pub const fn state(&self) -> GasState { self.state }

    /// The reproducible input accepted by the shared material parser.
    #[must_use]
    pub fn declaration(&self) -> &str { &self.declaration }
}

impl Default for Atmosphere {
    fn default() -> Self {
        Self::from_declaration("293.15 101325 0").expect("declared room air is admissible")
    }
}

/// One drumhead film: geometry-independent material, thickness and installation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeadParams {
    /// Film thickness [m].
    pub thickness_m: f64,
    /// Installed isotropic membrane tension [N/m].
    pub tension_n_m: f64,
    /// Viscous modal damping ratio (q'' + 2 zeta w q' + w^2 q = f).
    pub damping_ratio: f64,
    /// Equilibrated affine and harmonic tension variation (zero = uniform).
    pub tension_variation: TensionVariation,
}

/// Complete declared instrument, numerics and observation.
#[derive(Debug, Clone, PartialEq)]
pub struct DrumParams {
    /// Clear vibrating head radius [m].
    pub radius_m: f64,
    /// Shell depth = head separation [m]; sets the sealed air volume.
    pub depth_m: f64,
    /// Batter (struck, upper) head.
    pub batter: HeadParams,
    /// Resonant (lower) head.
    pub resonant: HeadParams,
    /// Film Young modulus [Pa] (both heads).
    pub young_pa: f64,
    /// Film Poisson ratio (both heads).
    pub poisson: f64,
    /// Film density [kg/m^3] (both heads).
    pub density_kg_m3: f64,
    /// FE mesh: radial intervals (1..=32).
    pub radial_intervals: usize,
    /// FE mesh: vertices per ring (8..=128).
    pub azimuths: usize,
    /// Retained head eigen-window [Hz], lower edge.
    pub band_low_hz: f64,
    /// Retained head eigen-window [Hz], upper edge.
    pub band_high_hz: f64,
    /// Mechanical steps per 48 kHz audio frame (1, 2, 4, 8 or 16).
    pub substeps: usize,
    /// Shared material gas for cavity stiffness and exterior radiation.
    pub atmosphere: Atmosphere,
    /// Microphone above the batter head, head-centred frame, z up from the head plane [m].
    pub mic_batter_m: [f64; 3],
    /// Microphone below the resonant head, same x/y frame, z = distance below that head [m].
    pub mic_resonant_m: [f64; 3],
    /// Fixed pressure that maps to float sample 1.0 [Pa]. Not adaptive.
    pub full_scale_pa: f64,
    /// Hunt-Crossley internal contact loss chi [s/m]; 0 = the reference's elastic Hertz tip.
    pub contact_loss_s_m: f64,
    /// Separation at which the departing stick is caught by the player [m].
    pub catch_gap_m: f64,
}

/// Hard ceiling on retained head coordinates (both heads together): the
/// reference's own 63-coordinate host budget (drum_spec.rs `Spec::prepare`).
pub const MAX_HEAD_MODES: usize = 63;
/// Fixed audio output rate [Hz].
pub const OUTPUT_RATE_HZ: u32 = 48_000;

impl Default for DrumParams {
    /// The offline reference's estimated 14 x 6.5 inch drum (drum_spec.rs
    /// `Spec::reference`), its 5 x 32 mesh and its 80-500 Hz head window.
    fn default() -> Self {
        let batter = HeadParams {
            thickness_m: 0.000_254,
            tension_n_m: 3000.0,
            damping_ratio: 0.001,
            tension_variation: TensionVariation::default(),
        };
        Self {
            radius_m: 0.1778 - 0.0075,
            depth_m: 0.1651,
            batter,
            resonant: HeadParams { thickness_m: 0.000_076_2, tension_n_m: 1500.0, ..batter },
            young_pa: 4e9,
            poisson: 0.38,
            density_kg_m3: 1390.0,
            radial_intervals: 5,
            azimuths: 32,
            band_low_hz: 80.0,
            band_high_hz: 500.0,
            substeps: 16,
            atmosphere: Atmosphere::default(),
            mic_batter_m: [0.08, 0.05, 0.35],
            mic_resonant_m: [0.08, 0.05, 0.35],
            // Measured with these defaults: a 2 m/s stroke at (0.06, 0.01) m
            // peaks at 191 Pa at this microphone and a 4 m/s stroke at 441 Pa.
            // Fixed, never adaptive; overs are counted, not limited.
            full_scale_pa: 400.0,
            contact_loss_s_m: 0.0,
            catch_gap_m: 0.005,
        }
    }
}

fn bad(message: impl Into<String>) -> Refusal {
    Refusal::new("bad_parameter", message)
}

fn scalar(key: &str, value: &Value) -> Result<f64, Refusal> {
    match value {
        Value::Number(v) => Ok(*v),
        _ => Err(bad(format!("\"{key}\" must be a finite number in SI units"))),
    }
}
fn count(key: &str, value: &Value) -> Result<usize, Refusal> {
    let v = scalar(key, value)?;
    if v < 0.0 || v.fract() != 0.0 || v > 1e6 {
        return Err(bad(format!("\"{key}\" must be a non-negative integer")));
    }
    Ok(v as usize)
}
fn vector<const N: usize>(key: &str, value: &Value) -> Result<[f64; N], Refusal> {
    match value {
        Value::Numbers(v) if v.len() == N => {
            let mut out = [0.0; N];
            out.copy_from_slice(v);
            Ok(out)
        }
        _ => Err(bad(format!("\"{key}\" must be an array of exactly {N} finite numbers"))),
    }
}
fn variation(
    key: &str,
    value: &Value,
    existing: TensionVariation,
) -> Result<TensionVariation, Refusal> {
    let v: [f64; 7] = vector(key, value)?;
    // The affine control owns seven coefficients, not the independent lug
    // harmonics. Preserve those on partial live edits and in either JSON order.
    Ok(TensionVariation {
        constant_n_m: [v[0], v[1], v[2]],
        gradient_n_m2: [v[3], v[4], v[5], v[6]],
        ..existing
    })
}
fn harmonics(key: &str, value: &Value) -> Result<[[f64;2]; TENSION_HARMONIC_COUNT], Refusal> {
    match value {
        Value::Numbers(v) if v.len() == 2*TENSION_HARMONIC_COUNT => {
            let mut out=[[0.0;2];TENSION_HARMONIC_COUNT];
            for (slot,pair) in out.iter_mut().zip(v.chunks_exact(2)) {
                slot.copy_from_slice(pair);
            }
            Ok(out)
        }
        _ => Err(bad(format!("\"{key}\" must contain cosine/sine N/m pairs for orders 2..12 (22 finite numbers)"))),
    }
}

impl DrumParams {
    /// Overlay a flat JSON object on `self`. Unknown keys are refused, so a
    /// misspelled slider cannot silently leave the instrument unchanged.
    ///
    /// # Errors
    /// `bad_json`, `unknown_parameter` or `bad_parameter`.
    pub fn overlay_json(&self, text: &str) -> Result<Self, Refusal> {
        self.clone().apply_json(text)
    }

    fn apply_json(mut self, text: &str) -> Result<Self, Refusal> {
        let entries = json::parse_flat_object(text).map_err(|e| Refusal::new("bad_json", e))?;
        for (key, value) in &entries {
            let k = key.as_str();
            match k {
                "radius_m" => self.radius_m = scalar(k, value)?,
                "depth_m" => self.depth_m = scalar(k, value)?,
                "batter_thickness_m" => self.batter.thickness_m = scalar(k, value)?,
                "resonant_thickness_m" => self.resonant.thickness_m = scalar(k, value)?,
                "batter_tension_n_m" => self.batter.tension_n_m = scalar(k, value)?,
                "resonant_tension_n_m" => self.resonant.tension_n_m = scalar(k, value)?,
                "damping_ratio" => {
                    let v = scalar(k, value)?;
                    self.batter.damping_ratio = v;
                    self.resonant.damping_ratio = v;
                }
                "batter_damping_ratio" => self.batter.damping_ratio = scalar(k, value)?,
                "resonant_damping_ratio" => self.resonant.damping_ratio = scalar(k, value)?,
                "batter_tension_variation" => {
                    self.batter.tension_variation = variation(k, value, self.batter.tension_variation)?;
                }
                "resonant_tension_variation" => {
                    self.resonant.tension_variation = variation(k, value, self.resonant.tension_variation)?;
                }
                "batter_tension_harmonics" => self.batter.tension_variation.harmonics_n_m = harmonics(k,value)?,
                "resonant_tension_harmonics" => self.resonant.tension_variation.harmonics_n_m = harmonics(k,value)?,
                "young_pa" => self.young_pa = scalar(k, value)?,
                "poisson" => self.poisson = scalar(k, value)?,
                "density_kg_m3" => self.density_kg_m3 = scalar(k, value)?,
                "radial_intervals" => self.radial_intervals = count(k, value)?,
                "azimuths" => self.azimuths = count(k, value)?,
                "band_low_hz" => self.band_low_hz = scalar(k, value)?,
                "band_high_hz" => self.band_high_hz = scalar(k, value)?,
                "substeps" => self.substeps = count(k, value)?,
                "gas" => match value {
                    Value::Text(text) => self.atmosphere = Atmosphere::from_declaration(text)?,
                    _ => return Err(bad("gas must be a material gas declaration string")),
                },
                "mic_batter_m" => self.mic_batter_m = vector(k, value)?,
                "mic_resonant_m" => self.mic_resonant_m = vector(k, value)?,
                "full_scale_pa" => self.full_scale_pa = scalar(k, value)?,
                "contact_loss_s_m" => self.contact_loss_s_m = scalar(k, value)?,
                "catch_gap_m" => self.catch_gap_m = scalar(k, value)?,
                _ => {
                    return Err(Refusal::new(
                        "unknown_parameter",
                        format!("\"{k}\" is not a parameter of the live drum"),
                    ));
                }
            }
        }
        Ok(self)
    }

    /// Admission of everything that does not need the FE owners. The plate,
    /// modal and radiation owners keep their own (stricter) refusals.
    ///
    /// # Errors
    /// `bad_parameter` naming the first violated physical or numerical bound.
    pub fn validate(&self) -> Result<(), Refusal> {
        let positive = |name: &str, v: f64| {
            if v.is_finite() && v > 0.0 { Ok(()) } else { Err(bad(format!("{name} must be finite and positive"))) }
        };
        positive("radius_m", self.radius_m)?;
        positive("depth_m", self.depth_m)?;
        if !(0.02..=1.0).contains(&self.radius_m) || !(0.01..=2.0).contains(&self.depth_m) {
            return Err(bad("radius_m must lie in 0.02..=1 m and depth_m in 0.01..=2 m"));
        }
        for (name, head) in [("batter", &self.batter), ("resonant", &self.resonant)] {
            positive(&format!("{name}_thickness_m"), head.thickness_m)?;
            positive(&format!("{name}_tension_n_m"), head.tension_n_m)?;
            head.tension_variation.validate().map_err(|e|bad(format!("{name} tension field: {e}")))?;
            if !head.damping_ratio.is_finite() || head.damping_ratio < 0.0 || head.damping_ratio >= 1.0 {
                return Err(bad(format!("{name}_damping_ratio must be a passive ratio in 0..1")));
            }
        }
        positive("young_pa", self.young_pa)?;
        positive("density_kg_m3", self.density_kg_m3)?;
        if !self.poisson.is_finite() || self.poisson <= -1.0 || self.poisson >= 0.5 {
            return Err(bad("poisson must lie in (-1, 0.5)"));
        }
        if !(1..=32).contains(&self.radial_intervals) || !(8..=128).contains(&self.azimuths) {
            return Err(bad("mesh needs 1..=32 radial intervals and 8..=128 azimuths"));
        }
        if !self.band_low_hz.is_finite() || !self.band_high_hz.is_finite()
            || self.band_low_hz < 0.0 || self.band_high_hz <= self.band_low_hz
        {
            return Err(bad("head window needs finite 0 <= band_low_hz < band_high_hz"));
        }
        if !matches!(self.substeps, 1 | 2 | 4 | 8 | 16) {
            return Err(bad("substeps must be 1, 2, 4, 8 or 16 (mechanical rate = 48 kHz x substeps)"));
        }
        // The reference's clock admission (drum_spec.rs admit_clock): the
        // declared window must stay below 0.45 of the mechanical rate.
        let rate = f64::from(OUTPUT_RATE_HZ) * self.substeps as f64;
        if self.band_high_hz / rate >= 0.45 {
            return Err(Refusal::new("nyquist_guard", "declared head window crosses the mechanical Nyquist guard"));
        }
        self.atmosphere.state().validate_coherence()
            .map_err(|e| Refusal::new("gas_refused", e.to_string()))?;
        positive("full_scale_pa", self.full_scale_pa)?;
        for (name, mic) in [("mic_batter_m", self.mic_batter_m), ("mic_resonant_m", self.mic_resonant_m)] {
            if mic.iter().any(|v| !v.is_finite()) || mic[2] < 0.05 || mic.iter().any(|v| v.abs() > 20.0) {
                return Err(bad(format!("{name} must be finite, within 20 m, and at least 0.05 m off its head plane")));
            }
        }
        if !self.contact_loss_s_m.is_finite() || self.contact_loss_s_m < 0.0 {
            return Err(bad("contact_loss_s_m (Hunt-Crossley chi) must be finite and non-negative"));
        }
        if !self.catch_gap_m.is_finite() || !(0.0005..=0.05).contains(&self.catch_gap_m) {
            return Err(bad("catch_gap_m must lie in 0.0005..=0.05 m"));
        }
        Ok(())
    }

    /// Weyl-law estimate of how many fixed-rim modes of one head fall inside
    /// the declared window: N(k) = (ka)^2/4 - ka/2, with k from the tensioned
    /// plate dispersion rho h w^2 = T k^2 + D k^4, D = E h^3 / (12 (1 - nu^2)).
    /// Used ONLY to refuse a hopeless window before the eigen-solve spends
    /// seconds on it; the eigen-solve remains the authority on every mode.
    #[must_use]
    pub fn estimated_head_modes(&self, head: &HeadParams) -> f64 {
        let areal = self.density_kg_m3 * head.thickness_m;
        let bending = self.young_pa * head.thickness_m.powi(3) / (12.0 * (1.0 - self.poisson * self.poisson));
        let tension = head.tension_n_m;
        let count = |frequency_hz: f64| {
            let w = core::f64::consts::TAU * frequency_hz;
            // k^2 = 2 sigma w^2 / (T + sqrt(T^2 + 4 D sigma w^2)): stable for D -> 0.
            let k2 = 2.0 * areal * w * w / (tension + (tension * tension + 4.0 * bending * areal * w * w).sqrt());
            let ka = k2.max(0.0).sqrt() * self.radius_m;
            (0.25 * ka * ka - 0.5 * ka).max(0.0)
        };
        (count(self.band_high_hz) - count(self.band_low_hz)).max(0.0)
    }

    /// Sealed air volume [m^3] of the rigid cylindrical shell.
    #[must_use]
    pub fn volume_m3(&self) -> f64 {
        core::f64::consts::PI * self.radius_m * self.radius_m * self.depth_m
    }
    /// Adiabatic small-signal bulk modulus rho c^2 [Pa].
    #[must_use]
    pub fn bulk_modulus_pa(&self) -> f64 {
        let gas = self.atmosphere.state();
        gas.density * gas.sound_speed * gas.sound_speed
    }
    /// Mechanical rate [Hz].
    #[must_use]
    pub fn mechanical_rate_hz(&self) -> u32 {
        OUTPUT_RATE_HZ * self.substeps as u32
    }
    /// True when `other` describes the same heads, mesh and window, so the
    /// existing eigen-solve can be reused (only numerics/observation changed).
    #[must_use]
    pub fn same_structure(&self, other: &Self) -> bool {
        self.radius_m == other.radius_m
            && self.batter == other.batter
            && self.resonant == other.resonant
            && self.young_pa == other.young_pa
            && self.poisson == other.poisson
            && self.density_kg_m3 == other.density_kg_m3
            && self.radial_intervals == other.radial_intervals
            && self.azimuths == other.azimuths
            && self.band_low_hz == other.band_low_hz
            && self.band_high_hz == other.band_high_hz
    }

    /// Every admitted value, as the flat object `overlay_json` accepts.
    #[must_use]
    pub fn to_json(&self) -> String {
        let tv = |v: &TensionVariation| {
            let mut all = v.constant_n_m.to_vec();
            all.extend_from_slice(&v.gradient_n_m2);
            json::numbers(&all)
        };
        let th = |v: &TensionVariation| {
            let all=v.harmonics_n_m.iter().flatten().copied().collect::<Vec<_>>();
            json::numbers(&all)
        };
        let n = json::number;
        format!(
            concat!(
                "{{\"radius_m\":{},\"depth_m\":{},",
                "\"batter_thickness_m\":{},\"resonant_thickness_m\":{},",
                "\"batter_tension_n_m\":{},\"resonant_tension_n_m\":{},",
                "\"batter_damping_ratio\":{},\"resonant_damping_ratio\":{},",
                "\"batter_tension_variation\":{},\"resonant_tension_variation\":{},",
                "\"batter_tension_harmonics\":{},\"resonant_tension_harmonics\":{},",
                "\"young_pa\":{},\"poisson\":{},\"density_kg_m3\":{},",
                "\"radial_intervals\":{},\"azimuths\":{},",
                "\"band_low_hz\":{},\"band_high_hz\":{},\"substeps\":{},",
                "\"gas\":{},",
                "\"mic_batter_m\":{},\"mic_resonant_m\":{},\"full_scale_pa\":{},",
                "\"contact_loss_s_m\":{},\"catch_gap_m\":{}}}"
            ),
            n(self.radius_m), n(self.depth_m),
            n(self.batter.thickness_m), n(self.resonant.thickness_m),
            n(self.batter.tension_n_m), n(self.resonant.tension_n_m),
            n(self.batter.damping_ratio), n(self.resonant.damping_ratio),
            tv(&self.batter.tension_variation), tv(&self.resonant.tension_variation),
            th(&self.batter.tension_variation), th(&self.resonant.tension_variation),
            n(self.young_pa), n(self.poisson), n(self.density_kg_m3),
            self.radial_intervals, self.azimuths,
            n(self.band_low_hz), n(self.band_high_hz), self.substeps,
            json::string(self.atmosphere.declaration()),
            json::numbers(&self.mic_batter_m), json::numbers(&self.mic_resonant_m), n(self.full_scale_pa),
            n(self.contact_loss_s_m), n(self.catch_gap_m),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_the_reference_drum_and_round_trip_through_json() {
        let d = DrumParams::default();
        d.validate().unwrap();
        assert!((d.radius_m - 0.1703).abs() < 1e-12);
        let again = DrumParams::default().overlay_json(&d.to_json()).unwrap();
        assert_eq!(d, again);
    }

    #[test]
    fn live_json_round_trips_and_builds_multi_lug_prestress() {
        let d=DrumParams::default();
        let mut values=vec![0.0;2*TENSION_HARMONIC_COUNT];
        values[2*(10-2)]=120.0; values[2*(10-2)+1]=-35.0;
        let json=format!(
            "{{\"batter_tension_harmonics\":{},\"batter_tension_variation\":[20,10,2,1,-2,3,-4]}}",
            json::numbers(&values)
        );
        let tuned=d.overlay_json(&json).unwrap();
        assert_eq!(tuned.batter.tension_variation.harmonics_n_m[8],[120.0,-35.0]);
        assert_eq!(tuned.batter.tension_variation.constant_n_m, [20.0, 10.0, 2.0]);
        assert_eq!(tuned.batter.tension_variation.gradient_n_m2, [1.0, -2.0, 3.0, -4.0]);
        tuned.validate().unwrap();
        let again=DrumParams::default().overlay_json(&tuned.to_json()).unwrap();
        assert_eq!(tuned,again);
        let plain=crate::instrument::Instrument::build(d).unwrap();
        let built=crate::instrument::Instrument::build(tuned).unwrap();
        assert_eq!(plain.heads[0].film.model.m,built.heads[0].film.model.m);
        assert!(plain.heads[0].modes.iter().zip(&built.heads[0].modes)
            .any(|(a,b)|(a.lambda-b.lambda).abs()>1e-8*a.lambda.max(1.0)));
    }


    #[test]
    fn affine_and_harmonic_controls_commute_for_both_heads() {
        let initial = DrumParams::default();
        let values: Vec<_> = (0..2 * TENSION_HARMONIC_COUNT)
            .map(|i| i as f64 - 10.0)
            .collect();
        for head in ["batter", "resonant"] {
            let affine = format!("\"{head}_tension_variation\":[20,10,2,1,-2,3,-4]");
            let harmonic = format!("\"{head}_tension_harmonics\":{}", json::numbers(&values));
            let affine_first = initial.overlay_json(&format!("{{{affine},{harmonic}}}")).unwrap();
            let harmonic_first = initial.overlay_json(&format!("{{{harmonic},{affine}}}")).unwrap();
            let sequential = initial
                .overlay_json(&format!("{{{harmonic}}}")).unwrap()
                .overlay_json(&format!("{{{affine}}}")).unwrap();
            assert_eq!(affine_first, harmonic_first);
            assert_eq!(affine_first, sequential);
            affine_first.validate().unwrap();
            assert!(!initial.same_structure(&affine_first));
            let round_trip = initial.overlay_json(&sequential.to_json()).unwrap();
            assert_eq!(round_trip, sequential);
            assert!(round_trip.same_structure(&sequential));
            let (changed, unchanged, original) = if head == "batter" {
                (sequential.batter, sequential.resonant, initial.resonant)
            } else {
                (sequential.resonant, sequential.batter, initial.batter)
            };
            assert_eq!(unchanged, original);
            assert_eq!(changed.tension_variation.constant_n_m, [20.0, 10.0, 2.0]);
            assert_eq!(changed.tension_variation.gradient_n_m2, [1.0, -2.0, 3.0, -4.0]);
            for (pair, expected) in changed.tension_variation.harmonics_n_m.iter()
                .zip(values.chunks_exact(2)) {
                assert_eq!(pair.as_slice(), expected);
            }
        }
    }

    #[test]
    fn explicit_zero_edits_clear_only_the_requested_tension_component() {
        let mut initial = DrumParams::default();
        for head in [&mut initial.batter, &mut initial.resonant] {
            head.tension_variation.constant_n_m = [20.0, 10.0, 2.0];
            head.tension_variation.gradient_n_m2 = [1.0, -2.0, 3.0, -4.0];
            head.tension_variation.harmonics_n_m = [[12.0, -3.0]; TENSION_HARMONIC_COUNT];
        }
        let zeros = json::numbers(&[0.0; 2 * TENSION_HARMONIC_COUNT]);
        for name in ["batter", "resonant"] {
            let mut expected = initial.clone();
            let head = if name == "batter" { &mut expected.batter } else { &mut expected.resonant };
            head.tension_variation.constant_n_m = [0.0; 3];
            head.tension_variation.gradient_n_m2 = [0.0; 4];
            let affine_reset = initial.overlay_json(&format!(
                "{{\"{name}_tension_variation\":[0,0,0,0,0,0,0]}}"
            )).unwrap();
            assert_eq!(affine_reset, expected);
            let mut expected = initial.clone();
            let head = if name == "batter" { &mut expected.batter } else { &mut expected.resonant };
            head.tension_variation.harmonics_n_m = [[0.0; 2]; TENSION_HARMONIC_COUNT];
            let harmonic_reset = initial.overlay_json(&format!(
                "{{\"{name}_tension_harmonics\":{zeros}}}"
            )).unwrap();
            assert_eq!(harmonic_reset, expected);
        }
    }

    #[test]
    fn typed_refusals_name_their_class() {
        let d = DrumParams::default();
        assert_eq!(d.overlay_json("{\"pitch_hz\":200}").unwrap_err().code, "unknown_parameter");
        assert_eq!(d.overlay_json("{\"radius_m\":").unwrap_err().code, "bad_json");
        assert_eq!(d.overlay_json("{\"radius_m\":\"big\"}").unwrap_err().code, "bad_parameter");
        for bad in [
            "{\"batter_tension_n_m\":-5}", "{\"batter_tension_n_m\":0}", "{\"radius_m\":5}",
            "{\"substeps\":3}", "{\"damping_ratio\":-0.1}", "{\"mic_batter_m\":[0,0,0.01]}",
            "{\"band_low_hz\":600}", "{\"azimuths\":4}", "{\"full_scale_pa\":0}",
        ] {
            let p = d.overlay_json(bad).unwrap();
            assert_eq!(p.validate().unwrap_err().code, "bad_parameter", "{bad}");
        }
        let p = d.overlay_json("{\"substeps\":1,\"band_high_hz\":22000}").unwrap();
        assert_eq!(p.validate().unwrap_err().code, "nyquist_guard");
    }

    #[test]
    fn material_gas_round_trip_and_pressure_temperature_composition_laws() {
        let base = DrumParams::default();
        let low = base.overlay_json(r#"{"gas":"293.15 50000 0"}"#).unwrap();
        let high = base.overlay_json(r#"{"gas":"293.15 100000 0"}"#).unwrap();
        let hot = base.overlay_json(r#"{"gas":"313.15 100000 0"}"#).unwrap();
        assert_eq!(low.atmosphere.state().sound_speed, high.atmosphere.state().sound_speed);
        assert_eq!(high.atmosphere.state().density, 2.0 * low.atmosphere.state().density);
        assert_eq!(high.bulk_modulus_pa(), 2.0 * low.bulk_modulus_pa());
        assert!(hot.atmosphere.state().sound_speed > high.atmosphere.state().sound_speed);
        assert!(hot.atmosphere.state().density < high.atmosphere.state().density);
        // Ideal-gas K = gamma p: heating at fixed pressure does not stiffen the cavity.
        assert!((hot.bulk_modulus_pa() / high.bulk_modulus_pa() - 1.0).abs() < 1e-14);
        let mixture = base.overlay_json(
            r#"{"gas":"mixture 293.15 100000 mole 3 N2 0.5 O2 0.2 CO2 0.3"}"#,
        ).unwrap();
        assert!(mixture.atmosphere.state().sound_speed < high.atmosphere.state().sound_speed);
        assert!(mixture.atmosphere.state().density > high.atmosphere.state().density);
        assert!(mixture.bulk_modulus_pa() < high.bulk_modulus_pa());
        for p in [low, high, hot, mixture] {
            p.validate().unwrap();
            assert!(p.same_structure(&base));
            assert_eq!(base.overlay_json(&p.to_json()).unwrap(), p);
        }
        assert!(base.overlay_json(r#"{"gas_density_kg_m3":2}"#).is_err());
        assert!(base.overlay_json(r#"{"gas_sound_speed_m_s":400}"#).is_err());
    }

    #[test]
    fn gas_admission_preserves_the_material_owners_refusals() {
        for gas in [
            "293.15 0 0", "293.15 101325 2", "293.15 101325 0 trailing",
            "mixture 293.15 100000 mole 2 N2 0.8 O2 0.8",
            "mixture 293.15 100000 mole 2 N2 0.9 H2O 0.1",
            "mixture 293.15 100000 mole 1 invented 1",
            "thermochemical 293.15 100000 mole 1 N2 1",
        ] {
            assert_eq!(Atmosphere::from_declaration(gas).unwrap_err().code, "gas_refused", "{gas}");
        }
    }
}
