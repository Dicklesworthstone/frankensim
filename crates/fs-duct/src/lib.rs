//! 1D viscothermal duct/horn acoustics
//! (bead frankensim-fsim-duct-acoustics-zdmm1, musical-acoustics
//! program): transfer matrices over cylinder and cone segments with
//! Zwikker–Kosten wide-tube losses, low-`ka` radiation loads, and
//! input impedance. A clarinet bore, a muffler, and an HVAC run are
//! the same object here — per the program doctrine there is no
//! instrument-specific code path, and every medium property is DERIVED
//! from an `fs_material::gas::GasState` (temperature, pressure, gas),
//! never hardcoded.
//!
//! TIME CONVENTION (pinned, matching `fs_bem::helmholtz`): fields vary
//! as `e^{-i omega t}`, travelling waves as `e^{i(kx - omega t)}`, so
//! attenuation requires `Im k > 0` and mass-like radiation reactance
//! is NEGATIVE imaginary. Both signs are pinned by executable tests
//! (passivity of `Re Z_in`, resonance flattening, decay), not trusted
//! from convention labels.
//!
//! PHYSICS — wide-tube (large shear number) Zwikker–Kosten first
//! order. With shear number `rv = r sqrt(rho omega / mu)` (the tube
//! radius over the viscous boundary-layer scale):
//!
//! `k  = k0 [1 + (1 + i)/(sqrt2 rv) (1 + (gamma - 1)/sqrt(Pr))]`
//! `Zc = (rho c / S) [1 + (1 + i)/(sqrt2 rv) (1 - (gamma - 1)/sqrt(Pr))]`
//!
//! equivalently `Im k = alpha = (1/(r c)) sqrt(nu omega / 2)
//! (1 + (gamma - 1)/sqrt(Pr))` — the classic Kirchhoff wall-loss
//! attenuation. Validity is a NAMED refusal, not degradation: below
//! [`MIN_SHEAR_NUMBER`] the wide-tube expansion is invalid (narrow
//! tubes need the full Bessel/Kelvin ZK solution — a recorded
//! follow-up), and the loss model states it covers straight rigid
//! smooth walls only.
//!
//! TRANSFER MATRICES are never transcribed from a book: each segment's
//! 2-port `[p1, U1] = M [p2, U2]` is built NUMERICALLY from the two
//! exact analytic basis solutions of its 1D Helmholtz problem — plane
//! waves `e^{+-ikx}` for cylinders, spherical waves `e^{+-ikx}/x` for
//! cones — with the transmission-line relation
//! `U = S p' / (i k z_specific)` under the pinned convention (the
//! lossless limit reduces to `U = S p'/(i omega rho)`; the lossy form
//! is what lets the ZK impedance correction enter the propagation).
//! Transcription errors are structurally impossible; the oracles below
//! (closed-form cotangent, cylinder limit of the cone, composition
//! exactness, complete-cone harmonics) arbitrate the construction.
//!
//! RADIATION LOADS are the classic low-`ka` fits: unflanged
//! `Z = (rho c / S)[(ka)^2/4 - i 0.6133 ka]` (Levine–Schwinger) and
//! flanged `Z = (rho c / S)[(ka)^2/2 - i 0.8216 ka]`, refused above
//! [`MAX_RADIATION_KA`] by name. The Helmholtz-BEM facility is the
//! recorded successor for computed loads.
//!
//! Determinism: sequential fixed-order arithmetic through
//! `fs_math::det` and `fs_la::eigen_complex::lu_complex`; repeat
//! evaluations are bitwise identical.
//!
//! Deferred with recorded triggers (see CONTRACT): full Bessel/Kelvin
//! narrow-tube ZK (needs fs-math Bessel functions), tone-hole
//! T-junction lattices and fingering tables (slice 3 of the bead),
//! multimodal horn expansion (trigger: bell mismatch beyond authored
//! tolerance), and BEM-computed radiation loads.

use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;

/// Wide-tube validity floor: below this shear number the first-order
/// Zwikker–Kosten expansion is refused by name.
pub const MIN_SHEAR_NUMBER: f64 = 10.0;

/// Low-`ka` radiation-fit validity ceiling.
pub const MAX_RADIATION_KA: f64 = 1.0;

/// Typed refusals with stable `FS-DUCT-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum DuctError {
    /// A parameter is non-finite or out of range.
    BadParameter {
        /// Which parameter refused.
        what: &'static str,
    },
    /// The wide-tube loss model is invalid for this tube/frequency.
    TooNarrow {
        /// Measured shear number.
        shear_number: f64,
    },
    /// The low-`ka` radiation fit is invalid at this frequency.
    RadiationKaTooLarge {
        /// Measured `ka` at the radiating mouth.
        ka: f64,
    },
    /// The duct has no segments.
    EmptyDuct,
    /// A dense solve refused (degenerate basis).
    Singular,
}

impl core::fmt::Display for DuctError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DuctError::BadParameter { what } => write!(f, "FS-DUCT-BAD-PARAMETER: {what}"),
            DuctError::TooNarrow { shear_number } => write!(
                f,
                "FS-DUCT-TOO-NARROW: shear number {shear_number:.2} below the wide-tube floor \
                 {MIN_SHEAR_NUMBER}"
            ),
            DuctError::RadiationKaTooLarge { ka } => write!(
                f,
                "FS-DUCT-RADIATION-KA: ka = {ka:.3} beyond the low-ka fit ceiling \
                 {MAX_RADIATION_KA}"
            ),
            DuctError::EmptyDuct => write!(f, "FS-DUCT-EMPTY: a duct needs at least one segment"),
            DuctError::Singular => write!(f, "FS-DUCT-SINGULAR: basis solve refused"),
        }
    }
}

impl std::error::Error for DuctError {}

/// Axial loss model for a segment chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LossModel {
    /// Ideal lossless propagation (`k` and `Zc` real) — the exact
    /// closed-form arm and the dropped-losses mutation baseline.
    Lossless,
    /// First-order wide-tube Zwikker–Kosten viscothermal losses.
    WideTube,
}

/// One duct segment. Lengths and radii in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Segment {
    /// Straight cylinder.
    Cylinder {
        /// Bore radius [m].
        radius: f64,
        /// Axial length [m].
        length: f64,
    },
    /// Truncated cone (linear radius taper; exact spherical-wave
    /// 1D propagation; viscothermal correction evaluated at the mean
    /// radius — the standard engineering treatment, documented, with
    /// refinement by subdivision available to callers).
    Cone {
        /// Inlet radius [m].
        inlet_radius: f64,
        /// Outlet radius [m].
        outlet_radius: f64,
        /// Axial length [m].
        length: f64,
    },
}

impl Segment {
    fn validate(&self) -> Result<(), DuctError> {
        let ok = |v: f64| v > 0.0 && v.is_finite();
        let valid = match *self {
            Segment::Cylinder { radius, length } => ok(radius) && ok(length),
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } => ok(inlet_radius) && ok(outlet_radius) && ok(length),
        };
        if valid {
            Ok(())
        } else {
            Err(DuctError::BadParameter {
                what: "segment radii and lengths must be positive and finite",
            })
        }
    }

    fn mean_radius(&self) -> f64 {
        match *self {
            Segment::Cylinder { radius, .. } => radius,
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                ..
            } => f64::midpoint(inlet_radius, outlet_radius),
        }
    }

    fn outlet_radius(&self) -> f64 {
        match *self {
            Segment::Cylinder { radius, .. } => radius,
            Segment::Cone { outlet_radius, .. } => outlet_radius,
        }
    }
}

/// Termination at the far (outlet) end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Termination {
    /// Rigid closed end (`Z -> infinity`).
    Closed,
    /// Ideal pressure-release open end (`Z = 0`) — the textbook limit
    /// used by the closed-form oracles.
    IdealOpen,
    /// Unflanged-pipe low-`ka` radiation load (Levine–Schwinger fit,
    /// end correction 0.6133 a).
    UnflangedOpen,
    /// Flanged-pipe low-`ka` radiation load (end correction 0.8216 a).
    FlangedOpen,
}

/// A duct: an ordered chain of segments from input to termination.
#[derive(Debug, Clone, PartialEq)]
pub struct Duct {
    /// Segments, input first.
    pub segments: Vec<Segment>,
}

/// The complex wavenumber and characteristic impedance of a segment at
/// one frequency, with the validity diagnostic that admitted them.
#[derive(Debug, Clone, Copy)]
pub struct SegmentWave {
    /// Complex axial wavenumber `k` [1/m] (`Im k >= 0` is decay).
    pub wavenumber: C64,
    /// Plane-wave characteristic impedance `Zc = rho c / S` with the
    /// viscothermal correction [Pa s/m^3].
    pub characteristic_impedance: C64,
    /// Area-free specific impedance `rho c (1 + eps_z)` [Pa s/m] — the
    /// transmission-line series relation the basis solutions use
    /// (`U = S p' / (i k z_specific)`), so the ZK impedance correction
    /// (equivalently the complex effective density) actually enters
    /// the propagation instead of being a decorative field.
    pub specific_impedance: C64,
    /// Shear number `rv` that admitted the wide-tube model (infinite
    /// for the lossless arm).
    pub shear_number: f64,
}

/// A solved input impedance with its validity diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct DuctResponse {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Input impedance `Z_in = p/U` at the inlet plane [Pa s/m^3].
    pub impedance: C64,
    /// Smallest shear number over the chain (the loss-model margin).
    pub min_shear_number: f64,
    /// `ka` at the radiating mouth (0 for closed/ideal terminations).
    pub mouth_ka: f64,
}

/// Viscothermal wavenumber and characteristic impedance for a tube of
/// the given radius: the Zwikker–Kosten wide-tube first order, derived
/// entirely from the first-principles [`GasState`].
///
/// # Errors
/// [`DuctError`] on bad inputs or a shear number below
/// [`MIN_SHEAR_NUMBER`] (narrow-tube regime, refused by name).
pub fn segment_wave(
    state: &GasState,
    radius: f64,
    omega: f64,
    loss: LossModel,
) -> Result<SegmentWave, DuctError> {
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "radius must be positive and finite",
        });
    }
    if !(omega > 0.0 && omega.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "frequency must be positive and finite",
        });
    }
    let area = core::f64::consts::PI * radius * radius;
    let k0 = omega / state.sound_speed;
    let z0 = state.density * state.sound_speed / area;
    match loss {
        LossModel::Lossless => Ok(SegmentWave {
            wavenumber: C64::from_re(k0),
            characteristic_impedance: C64::from_re(z0),
            specific_impedance: C64::from_re(state.density * state.sound_speed),
            shear_number: f64::INFINITY,
        }),
        LossModel::WideTube => {
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            if rv < MIN_SHEAR_NUMBER {
                return Err(DuctError::TooNarrow { shear_number: rv });
            }
            let thermal = (state.gamma - 1.0) / state.prandtl.sqrt();
            let scale = 1.0 / (core::f64::consts::SQRT_2 * rv);
            // (1 + i) scale (1 +- thermal): Im k > 0 decays under
            // e^{i(kx - omega t)}; the signs are pinned by the
            // passivity, flattening, and Q oracles below.
            let eps_k = C64::new(1.0, 1.0).scale(scale * (1.0 + thermal));
            let eps_z = C64::new(1.0, 1.0).scale(scale * (1.0 - thermal));
            Ok(SegmentWave {
                wavenumber: (C64::ONE + eps_k).scale(k0),
                characteristic_impedance: (C64::ONE + eps_z).scale(z0),
                specific_impedance: (C64::ONE + eps_z).scale(state.density * state.sound_speed),
                shear_number: rv,
            })
        }
    }
}

fn cis(theta: f64) -> C64 {
    C64::new(det::cos(theta), det::sin(theta))
}

/// `e^{i z}` for complex `z = a + b i`: `e^{-b} (cos a + i sin a)`.
fn exp_i(z: C64) -> C64 {
    cis(z.re).scale(det::exp(-z.im))
}

/// Basis solution pair `(p, U)` for one segment at axial position `t`
/// in `[0, length]`, under `e^{-i omega t}` with `U = S p'/(i omega
/// rho)`. Cylinder: plane waves `e^{+-ikt}`. Cone: spherical waves
/// `e^{+-ikx}/x` with `x` the distance from the (virtual) apex and the
/// area growing as `x^2`.
fn basis_pair(segment: &Segment, wave: &SegmentWave, t: f64, forward: bool) -> (C64, C64) {
    let k = wave.wavenumber;
    let sign = if forward { 1.0 } else { -1.0 };
    // Transmission-line series relation U = S p' / (i k z_specific):
    // for the forward cylinder wave this gives U = p S / z_specific
    // exactly, so the ZK impedance correction (complex effective
    // density) enters the propagation. Under e^{-i omega t} the
    // lossless limit reduces to U = S p'/(i omega rho).
    let series = (C64::new(0.0, 1.0) * k * wave.specific_impedance).recip();
    match *segment {
        Segment::Cylinder { radius, .. } => {
            let area = core::f64::consts::PI * radius * radius;
            // p = e^{+- i k t}; p' = +- i k p.
            let p = exp_i(k.scale(sign * t));
            let dp = C64::new(0.0, 1.0) * k.scale(sign) * p;
            (p, dp.scale(area) * series)
        }
        Segment::Cone {
            inlet_radius,
            outlet_radius,
            length,
        } => {
            // Linear taper r(t) = r1 + (r2 - r1) t / L; apex at
            // x = 0 with x(t) = x1 + t, x1 = r1 L / (r2 - r1) for an
            // expanding cone (contracting handled by the same algebra
            // with x decreasing; a zero-taper cone is delegated to the
            // cylinder basis by the caller).
            let slope = (outlet_radius - inlet_radius) / length;
            let x1 = inlet_radius / slope;
            let x = x1 + t;
            let r = inlet_radius + slope * t;
            let area = core::f64::consts::PI * r * r;
            // p = e^{+- i k x}/x; p' = (+- i k - 1/x) p.
            let p = exp_i(k.scale(sign * x)).scale(1.0 / x);
            let dp = (C64::new(0.0, sign) * k - C64::from_re(1.0 / x)) * p;
            (p, dp.scale(area) * series)
        }
    }
}

/// The segment 2-port `[p_in, U_in] = M [p_out, U_out]`, built
/// numerically from the analytic basis (no transcribed matrices).
fn segment_matrix(segment: &Segment, wave: &SegmentWave) -> Result<[C64; 4], DuctError> {
    let length = match *segment {
        Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => length,
    };
    // Delegate zero-taper cones to the cylinder basis (the apex
    // distance diverges).
    if let Segment::Cone {
        inlet_radius,
        outlet_radius,
        ..
    } = *segment
        && ((outlet_radius - inlet_radius).abs() < 1e-12 * inlet_radius)
    {
        let cyl = Segment::Cylinder {
            radius: inlet_radius,
            length,
        };
        return segment_matrix(&cyl, wave);
    }
    // Columns of the end-state matrices for the two basis solutions.
    let (p1f, u1f) = basis_pair(segment, wave, 0.0, true);
    let (p1b, u1b) = basis_pair(segment, wave, 0.0, false);
    let (p2f, u2f) = basis_pair(segment, wave, length, true);
    let (p2b, u2b) = basis_pair(segment, wave, length, false);
    // M = E1 * E2^{-1} with E = [[p_f, p_b], [U_f, U_b]].
    let e2 = [p2f, p2b, u2f, u2b];
    let lu = lu_complex(&e2, 2).map_err(|_| DuctError::Singular)?;
    // Solve E2^T? lu solves E2 x = b column-wise; we need E1 E2^{-1}:
    // compute columns of E2^{-1} by solving with unit vectors, then
    // multiply.
    let mut inv = [C64::ZERO; 4];
    for col in 0..2 {
        let mut b = [C64::ZERO, C64::ZERO];
        b[col] = C64::ONE;
        lu.solve(&mut b);
        inv[col] = b[0];
        inv[2 + col] = b[1];
    }
    let e1 = [p1f, p1b, u1f, u1b];
    let mut m = [C64::ZERO; 4];
    for row in 0..2 {
        for col in 0..2 {
            m[row * 2 + col] = e1[row * 2] * inv[col] + e1[row * 2 + 1] * inv[2 + col];
        }
    }
    Ok(m)
}

/// Termination impedance at the outlet.
///
/// # Errors
/// [`DuctError::RadiationKaTooLarge`] beyond the low-`ka` fit ceiling.
fn termination_impedance(
    termination: Termination,
    state: &GasState,
    radius: f64,
    omega: f64,
) -> Result<(Option<C64>, f64), DuctError> {
    let ka = omega / state.sound_speed * radius;
    let area = core::f64::consts::PI * radius * radius;
    let z0 = state.density * state.sound_speed / area;
    match termination {
        Termination::Closed => Ok((None, 0.0)),
        Termination::IdealOpen => Ok((Some(C64::ZERO), 0.0)),
        Termination::UnflangedOpen => {
            if ka > MAX_RADIATION_KA {
                return Err(DuctError::RadiationKaTooLarge { ka });
            }
            // Mass-like reactance is NEGATIVE imaginary under
            // e^{-i omega t} (the fs-bem pinned convention).
            Ok((Some(C64::new(0.25 * ka * ka, -0.6133 * ka).scale(z0)), ka))
        }
        Termination::FlangedOpen => {
            if ka > MAX_RADIATION_KA {
                return Err(DuctError::RadiationKaTooLarge { ka });
            }
            Ok((Some(C64::new(0.5 * ka * ka, -0.8216 * ka).scale(z0)), ka))
        }
    }
}

/// Input impedance of the duct at angular frequency `omega`, derived
/// entirely from the ambient [`GasState`].
///
/// # Errors
/// [`DuctError`] on empty/invalid geometry, the narrow-tube refusal,
/// or the radiation-fit ceiling.
pub fn input_impedance(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
) -> Result<DuctResponse, DuctError> {
    if duct.segments.is_empty() {
        return Err(DuctError::EmptyDuct);
    }
    for segment in &duct.segments {
        segment.validate()?;
    }
    // Chain matrix input -> outlet.
    let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let mut min_rv = f64::INFINITY;
    for segment in &duct.segments {
        let wave = segment_wave(state, segment.mean_radius(), omega, loss)?;
        min_rv = min_rv.min(wave.shear_number);
        let s = segment_matrix(segment, &wave)?;
        let mut next = [C64::ZERO; 4];
        for row in 0..2 {
            for col in 0..2 {
                next[row * 2 + col] = m[row * 2] * s[col] + m[row * 2 + 1] * s[2 + col];
            }
        }
        m = next;
    }
    let mouth_radius = duct
        .segments
        .last()
        .expect("non-empty checked above")
        .outlet_radius();
    let (z_load, mouth_ka) = termination_impedance(termination, state, mouth_radius, omega)?;
    let impedance = match z_load {
        // Z_in = (A Z_L + B)/(C Z_L + D); closed end is the C, D limit.
        Some(zl) => (m[0] * zl + m[1]) * (m[2] * zl + m[3]).recip(),
        None => m[0] * m[2].recip(),
    };
    Ok(DuctResponse {
        omega,
        impedance,
        min_shear_number: min_rv,
        mouth_ka,
    })
}

/// Uniformly sampled input-impedance sweep (inclusive endpoints).
///
/// # Errors
/// As for [`input_impedance`].
pub fn impedance_sweep(
    duct: &Duct,
    state: &GasState,
    omega_lo: f64,
    omega_hi: f64,
    count: usize,
    loss: LossModel,
    termination: Termination,
) -> Result<Vec<DuctResponse>, DuctError> {
    if !(omega_hi > omega_lo && omega_lo > 0.0) || count < 2 {
        return Err(DuctError::BadParameter {
            what: "sweep needs omega_hi > omega_lo > 0 and count >= 2",
        });
    }
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let omega = omega_lo + (omega_hi - omega_lo) * i as f64 / (count - 1) as f64;
        out.push(input_impedance(duct, state, omega, loss, termination)?);
    }
    Ok(out)
}

/// Indices of local maxima of `|Z_in|` in a sweep (interior points
/// strictly above both neighbours) — the resonance-peak helper the
/// oracles use.
#[must_use]
pub fn impedance_peaks(sweep: &[DuctResponse]) -> Vec<usize> {
    let mut peaks = Vec::new();
    for i in 1..sweep.len().saturating_sub(1) {
        let a = sweep[i - 1].impedance.abs();
        let b = sweep[i].impedance.abs();
        let c = sweep[i + 1].impedance.abs();
        if b > a && b > c {
            peaks.push(i);
        }
    }
    peaks
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn cot(x: f64) -> f64 {
        det::cos(x) / det::sin(x)
    }

    #[test]
    fn lossless_cylinder_matches_closed_forms() {
        // Under e^{-i omega t}: closed pipe Z_in = +i Zc cot(kL),
        // ideal-open pipe Z_in = -i Zc tan(kL) — both exact, computed
        // in-test from the scalar formulas, arbitrating the numeric
        // basis construction to near machine precision.
        let state = air();
        let (radius, length) = (0.0075, 0.5);
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let zc = state.density * state.sound_speed / (core::f64::consts::PI * radius * radius);
        for &f in &[80.0, 233.0, 617.0, 1201.0] {
            let omega = 2.0 * core::f64::consts::PI * f;
            let k = omega / state.sound_speed;
            let closed = input_impedance(
                &duct,
                &state,
                omega,
                LossModel::Lossless,
                Termination::Closed,
            )
            .expect("closed");
            let expected = C64::new(0.0, zc * cot(k * length));
            assert!(
                (closed.impedance - expected).abs() < 1e-9 * expected.abs().max(zc),
                "closed at {f} Hz: {:?} vs {expected:?}",
                closed.impedance
            );
            let open = input_impedance(
                &duct,
                &state,
                omega,
                LossModel::Lossless,
                Termination::IdealOpen,
            )
            .expect("open");
            let expected = C64::new(0.0, -zc * det::sin(k * length) / det::cos(k * length));
            assert!(
                (open.impedance - expected).abs() < 1e-9 * expected.abs().max(zc),
                "open at {f} Hz: {:?} vs {expected:?}",
                open.impedance
            );
        }
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"lossless-cylinder-closed-forms\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn quarter_wave_peaks_carry_the_end_correction() {
        use core::fmt::Write as _;
        // Unflanged open cylinder driven at the other end: |Z_in|
        // peaks at f_n ~ (2n - 1) c / (4 (L + 0.6133 a)) — the end
        // correction is IN the radiation load, so the peak ladder must
        // land on the corrected length, not the geometric one.
        let state = air();
        let (radius, length) = (0.0075, 0.5);
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let l_eff = length + 0.6133 * radius;
        let f_base = state.sound_speed / (4.0 * l_eff);
        let sweep = impedance_sweep(
            &duct,
            &state,
            2.0 * core::f64::consts::PI * 0.3 * f_base,
            2.0 * core::f64::consts::PI * 5.6 * f_base,
            6000,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("sweep");
        let peaks = impedance_peaks(&sweep);
        assert!(peaks.len() >= 3, "need three quarter-wave peaks");
        let mut rows = String::new();
        for (n, &idx) in peaks.iter().take(3).enumerate() {
            let f_peak = sweep[idx].omega / (2.0 * core::f64::consts::PI);
            let f_pred = (2.0 * n as f64 + 1.0) * f_base;
            let rel = (f_peak / f_pred - 1.0).abs();
            // Viscothermal dispersion flattens peaks slightly below the
            // lossless prediction; 0.7% absorbs it at this bore.
            assert!(
                rel < 7e-3,
                "peak {n}: {f_peak:.2} Hz vs corrected {f_pred:.2} Hz (rel {rel:.4})"
            );
            write!(
                rows,
                "{}{{\"n\":{n},\"f\":{f_peak:.2},\"pred\":{f_pred:.2}}}",
                if n == 0 { "" } else { "," }
            )
            .expect("write");
        }
        // Mutation-style contrast: against the UNcorrected length the
        // first peak misses by the full end-correction fraction.
        let f_uncorrected = state.sound_speed / (4.0 * length);
        let f1 = sweep[peaks[0]].omega / (2.0 * core::f64::consts::PI);
        assert!(
            (f1 / f_uncorrected - 1.0).abs() > 5e-3,
            "the corrected and uncorrected ladders must be distinguishable"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"quarter-wave-end-correction\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn cone_matrix_oracles() {
        // Four independent pins on the numerically-built cone 2-port:
        // (1) near-zero taper reproduces the cylinder matrix, (2) the
        // (p, U) determinant is exactly 1 (constant-Wronskian result:
        // r proportional to apex distance makes S/x^2 constant),
        // (3) two half-cones compose into the full cone, (4) a
        // near-complete closed-open cone resonates at ALL harmonics
        // n c / (2 L_apex) — the classic conical-bore result.
        let state = air();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let wave_cyl = segment_wave(&state, 0.01, omega, LossModel::Lossless).expect("wave");
        // (1) tiny taper vs cylinder.
        let cone = Segment::Cone {
            inlet_radius: 0.01,
            outlet_radius: 0.01 * (1.0 + 1e-6),
            length: 0.3,
        };
        let cyl = Segment::Cylinder {
            radius: 0.01,
            length: 0.3,
        };
        let m_cone = segment_matrix(&cone, &wave_cyl).expect("cone");
        let m_cyl = segment_matrix(&cyl, &wave_cyl).expect("cyl");
        for i in 0..4 {
            assert!(
                (m_cone[i] - m_cyl[i]).abs() < 1e-4 * m_cyl[i].abs().max(1e-6),
                "cylinder limit entry {i}: {:?} vs {:?}",
                m_cone[i],
                m_cyl[i]
            );
        }
        // (2) det = 1 for both segment kinds.
        for m in [&m_cone, &m_cyl] {
            let det_m = m[0] * m[3] - m[1] * m[2];
            assert!(
                (det_m - C64::ONE).abs() < 1e-10,
                "2-port determinant must be 1: {det_m:?}"
            );
        }
        // (3) composition exactness: half . half = whole.
        let whole = Segment::Cone {
            inlet_radius: 0.005,
            outlet_radius: 0.02,
            length: 0.4,
        };
        let first = Segment::Cone {
            inlet_radius: 0.005,
            outlet_radius: 0.0125,
            length: 0.2,
        };
        let second = Segment::Cone {
            inlet_radius: 0.0125,
            outlet_radius: 0.02,
            length: 0.2,
        };
        let wave_mean = segment_wave(&state, 0.0125, omega, LossModel::Lossless).expect("wave");
        let m_whole = segment_matrix(&whole, &wave_mean).expect("whole");
        let m_first = segment_matrix(&first, &wave_mean).expect("first");
        let m_second = segment_matrix(&second, &wave_mean).expect("second");
        let mut prod = [C64::ZERO; 4];
        for row in 0..2 {
            for col in 0..2 {
                prod[row * 2 + col] =
                    m_first[row * 2] * m_second[col] + m_first[row * 2 + 1] * m_second[2 + col];
            }
        }
        for i in 0..4 {
            assert!(
                (prod[i] - m_whole[i]).abs() < 1e-9 * m_whole[i].abs().max(1e-9),
                "composition entry {i}: {:?} vs {:?}",
                prod[i],
                m_whole[i]
            );
        }
        // (4) near-complete cone: harmonics at n c / (2 L_apex).
        let (r1, r2, length) = (0.001, 0.02, 0.5);
        let l_apex = length + r1 * length / (r2 - r1);
        let duct = Duct {
            segments: vec![Segment::Cone {
                inlet_radius: r1,
                outlet_radius: r2,
                length,
            }],
        };
        let f_base = state.sound_speed / (2.0 * l_apex);
        let sweep = impedance_sweep(
            &duct,
            &state,
            2.0 * core::f64::consts::PI * 0.5 * f_base,
            2.0 * core::f64::consts::PI * 3.6 * f_base,
            8000,
            LossModel::Lossless,
            Termination::IdealOpen,
        )
        .expect("sweep");
        let peaks = impedance_peaks(&sweep);
        assert!(peaks.len() >= 3, "need three cone harmonics");
        for (n, &idx) in peaks.iter().take(3).enumerate() {
            let f_peak = sweep[idx].omega / (2.0 * core::f64::consts::PI);
            let f_pred = (n as f64 + 1.0) * f_base;
            let rel = (f_peak / f_pred - 1.0).abs();
            assert!(
                rel < 0.03,
                "cone harmonic {n}: {f_peak:.1} vs {f_pred:.1} (rel {rel:.4})"
            );
        }
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"cone-oracles\",\"l_apex\":{l_apex:.4},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn resonance_q_matches_kirchhoff_wall_losses() {
        // INDEPENDENT loss oracle: for a closed cylinder the |Z_in|
        // resonance at k L = pi has quality factor Q = k/(2 alpha) with
        // alpha the Kirchhoff wall-loss coefficient. The test measures
        // Q from half-power bandwidth of the full TMM sweep and
        // compares against alpha computed from the CLOSED-FORM
        // nu-based expression — different arithmetic all the way down
        // from the implementation's shear-number form.
        let state = air();
        let (radius, length) = (0.0075, 0.5);
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let f0 = state.sound_speed / (2.0 * length);
        let lo = 2.0 * core::f64::consts::PI * 0.90 * f0;
        let hi = 2.0 * core::f64::consts::PI * 1.10 * f0;
        let sweep = impedance_sweep(
            &duct,
            &state,
            lo,
            hi,
            20_000,
            LossModel::WideTube,
            Termination::Closed,
        )
        .expect("sweep");
        let peaks = impedance_peaks(&sweep);
        assert_eq!(peaks.len(), 1, "one resonance in the window");
        let peak = peaks[0];
        let z_peak = sweep[peak].impedance.abs();
        let half = z_peak / core::f64::consts::SQRT_2;
        let mut i_lo = peak;
        while i_lo > 0 && sweep[i_lo].impedance.abs() > half {
            i_lo -= 1;
        }
        let mut i_hi = peak;
        while i_hi + 1 < sweep.len() && sweep[i_hi].impedance.abs() > half {
            i_hi += 1;
        }
        assert!(i_lo > 0 && i_hi + 1 < sweep.len(), "band inside window");
        let q_measured = sweep[peak].omega / (sweep[i_hi].omega - sweep[i_lo].omega);
        // Closed-form Kirchhoff alpha at the peak frequency.
        let omega = sweep[peak].omega;
        let nu = state.dynamic_viscosity / state.density;
        let alpha = (nu * omega / 2.0).sqrt() / (radius * state.sound_speed)
            * (1.0 + (state.gamma - 1.0) / state.prandtl.sqrt());
        let q_theory = (omega / state.sound_speed) / (2.0 * alpha);
        let rel = (q_measured / q_theory - 1.0).abs();
        assert!(
            rel < 0.03,
            "Q: measured {q_measured:.1} vs Kirchhoff {q_theory:.1} (rel {rel:.4})"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"kirchhoff-q\",\"q_measured\":{q_measured:.1},\"q_theory\":{q_theory:.1},\"rel\":{rel:.4},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn passivity_flattening_and_thermal_share() {
        // Sign pins: (a) Re Z_in >= 0 across the band for a radiating
        // viscothermal duct (a wrong eps sign pumps energy), (b) lossy
        // peaks sit BELOW their lossless counterparts (phase velocity
        // deficit), (c) the thermal loss share for air is the
        // documented ~47% on top of viscous — a dropped thermal term
        // fails loudly.
        let state = air();
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0075,
                length: 0.5,
            }],
        };
        let lo = 2.0 * core::f64::consts::PI * 60.0;
        let hi = 2.0 * core::f64::consts::PI * 1500.0;
        let lossy = impedance_sweep(
            &duct,
            &state,
            lo,
            hi,
            4000,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("lossy");
        for r in &lossy {
            assert!(
                r.impedance.re >= -1e-9,
                "passivity violated at omega {}: {:?}",
                r.omega,
                r.impedance
            );
        }
        let lossless = impedance_sweep(
            &duct,
            &state,
            lo,
            hi,
            4000,
            LossModel::Lossless,
            Termination::UnflangedOpen,
        )
        .expect("lossless");
        let p_lossy = impedance_peaks(&lossy);
        let p_lossless = impedance_peaks(&lossless);
        assert!(!p_lossy.is_empty() && p_lossy.len() == p_lossless.len());
        for (a, b) in p_lossy.iter().zip(p_lossless.iter()) {
            assert!(
                lossy[*a].omega < lossless[*b].omega,
                "viscothermal peaks must flatten below lossless"
            );
        }
        // Thermal share.
        let omega = 2.0 * core::f64::consts::PI * 500.0;
        let wave = segment_wave(&state, 0.0075, omega, LossModel::WideTube).expect("wave");
        let k0 = omega / state.sound_speed;
        let viscous_only = k0 / (core::f64::consts::SQRT_2 * wave.shear_number);
        let ratio = wave.wavenumber.im / viscous_only;
        assert!(
            (1.40..1.55).contains(&ratio),
            "thermal share for air must be ~1.47x viscous: {ratio:.3}"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"passivity-flattening-thermal\",\"thermal_ratio\":{ratio:.3},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn hot_duct_scales_with_the_independent_sqrt_t_law() {
        // Doctrine arm: the same duct at 700 K — every lossless
        // resonance scales by sqrt(700/293.15), pinned against the
        // independently computed constant (not against the gas module).
        let cold = air();
        let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 700.0, 101_325.0).expect("hot");
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0075,
                length: 0.5,
            }],
        };
        let peak_f = |state: &GasState| -> f64 {
            let f0 = state.sound_speed / (2.0 * 0.5);
            let sweep = impedance_sweep(
                &duct,
                state,
                2.0 * core::f64::consts::PI * 0.8 * f0,
                2.0 * core::f64::consts::PI * 1.2 * f0,
                8000,
                LossModel::Lossless,
                Termination::Closed,
            )
            .expect("sweep");
            let peaks = impedance_peaks(&sweep);
            sweep[peaks[0]].omega
        };
        let ratio = peak_f(&hot) / peak_f(&cold);
        let expected = (700.0f64 / 293.15).sqrt();
        assert!(
            (ratio - expected).abs() < 2e-4 * expected,
            "hot/cold resonance ratio {ratio:.5} vs independent sqrt(T) {expected:.5}"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"hot-duct-sqrt-t\",\"ratio\":{ratio:.5},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn named_refusals_fire_and_repeats_are_bitwise() {
        let state = air();
        // Narrow tube at low frequency: rv below the wide-tube floor.
        let narrow = segment_wave(
            &state,
            5e-4,
            2.0 * core::f64::consts::PI * 50.0,
            LossModel::WideTube,
        )
        .unwrap_err();
        assert!(narrow.to_string().contains("FS-DUCT-TOO-NARROW"));
        // Radiation fit ceiling.
        let fat = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.05,
                length: 0.3,
            }],
        };
        let err = input_impedance(
            &fat,
            &state,
            2.0 * core::f64::consts::PI * 2000.0,
            LossModel::Lossless,
            Termination::UnflangedOpen,
        )
        .unwrap_err();
        assert!(err.to_string().contains("FS-DUCT-RADIATION-KA"));
        assert!(matches!(
            input_impedance(
                &Duct { segments: vec![] },
                &state,
                1000.0,
                LossModel::Lossless,
                Termination::Closed
            ),
            Err(DuctError::EmptyDuct)
        ));
        assert!(matches!(
            input_impedance(
                &Duct {
                    segments: vec![Segment::Cylinder {
                        radius: -1.0,
                        length: 0.1
                    }]
                },
                &state,
                1000.0,
                LossModel::Lossless,
                Termination::Closed
            ),
            Err(DuctError::BadParameter { .. })
        ));
        // Bitwise determinism.
        let duct = Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: 0.0075,
                    length: 0.3,
                },
                Segment::Cone {
                    inlet_radius: 0.0075,
                    outlet_radius: 0.02,
                    length: 0.15,
                },
            ],
        };
        let omega = 2.0 * core::f64::consts::PI * 440.0;
        let a = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("a");
        let b = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("b");
        assert_eq!(a.impedance.re.to_bits(), b.impedance.re.to_bits());
        assert_eq!(a.impedance.im.to_bits(), b.impedance.im.to_bits());
    }
}
