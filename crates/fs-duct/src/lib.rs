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
//! follow-up). A locally reacting wall is an optional
//! [`fs_phs::WallPin`] on [`input_impedance_wall`].
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
//! `LossModel::Bessel` is the frequency-by-frequency Zwikker–Kosten
//! wall law (`fs_phs::zwikker_kosten_f`). A tone-hole chimney is a
//! short cylinder plus a flanged mouth (AllRegime when the bore
//! asked for WideTube, so a narrow neck does not refuse).
//! Lossy cones cascade spherical substations at local radius
//! (lossless stays the exact one-shot `e^{±ikx}/x` 2-port).
//! Deferred: fingering tables (slice 3 of the bead), multimodal
//! horn expansion (trigger: bell
//! mismatch beyond authored tolerance), and BEM-computed radiation
//! loads.

use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;
use fs_phs::{
    WallPin, side_hole_inner_length, side_hole_mutual_length, side_hole_neck_length,
    side_hole_series_length, wall_admittance_per_metre, zwikker_kosten_f,
};

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
    /// All-regime viscothermal law: wide-tube ZK above
    /// [`MIN_SHEAR_NUMBER`], Poiseuille + isothermal-tending thermal
    /// shunt below. Never silently drops losses.
    AllRegime,
    /// Frequency-by-frequency Zwikker–Kosten: `F(r_v) = 2 J₁(ζ)/(ζ J₀(ζ))`
    /// at this `ω`, valid at every shear number. A cone still uses the
    /// spherical `e^{\pm ikx}/x` basis. This is not a Foster fit.
    Bessel,
}

/// Open/closed state of a tone hole, or a vent fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HoleState {
    /// Hole open to the exterior through its chimney.
    Open,
    /// Hole sealed at the outer end (pad down): a closed cavity.
    Closed,
    /// Opening `σ ∈ (0, 1)`: `Y = σ Y_open + (1−σ) Y_closed`.
    Vent(f64),
}

fn hole_sigma(state: HoleState) -> f64 {
    match state {
        HoleState::Open => 1.0,
        HoleState::Closed => 0.0,
        HoleState::Vent(s) => s.clamp(0.0, 1.0),
    }
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
    /// 1D propagation). Lossy cones cascade substations at local
    /// radius; lossless is the one-shot `e^{±ikx}/x` 2-port.
    Cone {
        /// Inlet radius [m].
        inlet_radius: f64,
        /// Outlet radius [m].
        outlet_radius: f64,
        /// Axial length [m].
        length: f64,
    },
    /// A tone hole: a zero-axial-length side branch at this point in
    /// the chain — compact T-junction `series(Z_s/2) · shunt(Z_h) ·
    /// series(Z_s/2)` with Nederveen `t_s = −0.37 b²/a` on an open
    /// hole (closed pads stay a pure shunt).
    /// OPEN: chimney mass with Dalmont inner matching on `b/a`
    /// plus wall-flanged `0.8216 b` and flanged radiation resistance.
    /// CLOSED: the chimney cavity's compliance.
    ToneHole {
        /// Hole (chimney) radius b [m]; must be smaller than the bore
        /// radius.
        hole_radius: f64,
        /// Chimney height h (wall thickness plus pad lift geometry)
        /// [m].
        chimney_height: f64,
        /// Main-bore radius at the hole [m].
        bore_radius: f64,
        /// Open or closed.
        state: HoleState,
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
            Segment::ToneHole {
                hole_radius,
                chimney_height,
                bore_radius,
                ..
            } => {
                ok(hole_radius)
                    && ok(chimney_height)
                    && ok(bore_radius)
                    && hole_radius < bore_radius
            }
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
            Segment::ToneHole { bore_radius, .. } => bore_radius,
        }
    }

    /// Radius at the segment's outlet plane.
    #[must_use]
    pub fn outlet_radius(&self) -> f64 {
        match *self {
            Segment::Cylinder { radius, .. } => radius,
            Segment::Cone { outlet_radius, .. } => outlet_radius,
            Segment::ToneHole { bore_radius, .. } => bore_radius,
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
    /// Mouth volume-velocity over inlet volume-velocity from the ABCD
    /// chain and `Z_L`. Zero at a closed termination.
    pub u_mouth_over_u_in: C64,
    /// Mouth pressure over inlet volume-velocity (`Z_L U_mouth / U_in`).
    /// Zero at a closed termination.
    pub p_mouth_over_u_in: C64,
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
        LossModel::WideTube => wide_tube_wave(state, radius, omega, k0, z0, true),
        LossModel::AllRegime => {
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            if rv >= MIN_SHEAR_NUMBER {
                wide_tube_wave(state, radius, omega, k0, z0, false)
            } else {
                poiseuille_wave(state, radius, omega, rv)
            }
        }
        LossModel::Bessel => bessel_wave(state, radius, omega),
    }
}

/// Fold a locally reacting wall into a gas [`SegmentWave`].
///
/// Reconstructs the telegraph pair `Z' = i k Zc`, `Y' = i k / Zc`
/// (`Im Z' > 0` for inertance, same as Poiseuille), adds the wall
/// shunt with Im flipped from the `e^{-iωt}` [`WallPin`] into that
/// telegraph convention, and rebuilds `k` and `Zc`.
fn apply_wall_to_wave(
    wave: SegmentWave,
    radius: f64,
    omega: f64,
    wall: &WallPin,
) -> Result<SegmentWave, DuctError> {
    let y_phs =
        wall_admittance_per_metre(wall, radius, omega).map_err(|_| DuctError::BadParameter {
            what: "wall pin surface density, stiffness, and resistance",
        })?;
    let y_w = C64::new(y_phs.re, -y_phs.im);
    let j = C64::new(0.0, 1.0);
    let z_series = j * wave.wavenumber * wave.characteristic_impedance;
    let y_gas = j * wave.wavenumber * wave.characteristic_impedance.recip();
    let y_shunt = y_gas + y_w;
    let mut k = (z_series * y_shunt).sqrt();
    if k.im < 0.0 {
        k = k.scale(-1.0);
    }
    let mut zc = (z_series * y_shunt.recip()).sqrt();
    if zc.re < 0.0 {
        zc = zc.scale(-1.0);
    }
    let area = core::f64::consts::PI * radius * radius;
    Ok(SegmentWave {
        wavenumber: k,
        characteristic_impedance: zc,
        specific_impedance: zc.scale(area),
        shear_number: wave.shear_number,
    })
}

fn wide_tube_wave(
    state: &GasState,
    radius: f64,
    omega: f64,
    k0: f64,
    z0: f64,
    refuse_narrow: bool,
) -> Result<SegmentWave, DuctError> {
    let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
    if refuse_narrow && rv < MIN_SHEAR_NUMBER {
        return Err(DuctError::TooNarrow { shear_number: rv });
    }
    let thermal = (state.gamma - 1.0) / state.prandtl.sqrt();
    let scale = 1.0 / (core::f64::consts::SQRT_2 * rv);
    // (1 + i) scale (1 +- thermal): Im k > 0 decays under
    // e^{i(kx - omega t)}. The eps_k signs are pinned by the
    // passivity/flattening/Q oracles; eps_z is pinned by the
    // independent sqrt(Z_series/Y_shunt) route oracle (a
    // review found every eps_z mutation survived the physics
    // oracles — the impedance correction needs its own pin).
    let eps_k = C64::new(1.0, 1.0).scale(scale * (1.0 + thermal));
    let eps_z = C64::new(1.0, 1.0).scale(scale * (1.0 - thermal));
    Ok(SegmentWave {
        wavenumber: (C64::ONE + eps_k).scale(k0),
        characteristic_impedance: (C64::ONE + eps_z).scale(z0),
        specific_impedance: (C64::ONE + eps_z).scale(state.density * state.sound_speed),
        shear_number: rv,
    })
}

/// Narrow-tube (Poiseuille) series impedance plus an isothermal-tending
/// thermal shunt. Used only below [`MIN_SHEAR_NUMBER`].
fn poiseuille_wave(
    state: &GasState,
    radius: f64,
    omega: f64,
    rv: f64,
) -> Result<SegmentWave, DuctError> {
    let area = core::f64::consts::PI * radius * radius;
    let r_visc = 8.0 * state.dynamic_viscosity / (area * radius * radius);
    let l_ax = 4.0 * state.density / (3.0 * area);
    let c_ad = area / (state.density * state.sound_speed * state.sound_speed);
    let rt = rv * state.prandtl.sqrt();
    let iso_frac = 1.0 / (1.0 + 0.25 * rt * rt);
    let c_eff = c_ad * (1.0 + (state.gamma - 1.0) * iso_frac);
    let g_th = (state.gamma - 1.0) * c_ad * omega * (rt * rt / 16.0).min(0.5);
    let z_series = C64::new(r_visc, omega * l_ax);
    let y_shunt = C64::new(g_th.max(0.0), omega * c_eff);
    let mut k = (z_series * y_shunt).sqrt();
    if k.im < 0.0 {
        k = k.scale(-1.0);
    }
    let mut zc = (z_series * y_shunt.recip()).sqrt();
    if zc.re < 0.0 {
        zc = zc.scale(-1.0);
    }
    Ok(SegmentWave {
        wavenumber: k,
        characteristic_impedance: zc,
        specific_impedance: zc.scale(area),
        shear_number: rv,
    })
}

/// Full Zwikker–Kosten `k` and `Zc` from `F(r_v)` at this frequency.
fn bessel_wave(state: &GasState, radius: f64, omega: f64) -> Result<SegmentWave, DuctError> {
    let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
    if !(rv > 0.0 && rv.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "Bessel shear number",
        });
    }
    let f_v = zwikker_kosten_f(rv).map_err(|_| DuctError::BadParameter { what: "Bessel F_v" })?;
    let rt = rv * state.prandtl.sqrt();
    let mut f_t =
        zwikker_kosten_f(rt).map_err(|_| DuctError::BadParameter { what: "Bessel F_t" })?;
    let area = core::f64::consts::PI * radius * radius;
    let jw = C64::new(0.0, omega);
    let den = C64::ONE - f_v;
    if den.abs() < 1.0e-18 {
        return poiseuille_wave(state, radius, omega, rv);
    }
    let z_series = jw.scale(state.density / area) / den;
    let mut y_shunt = jw.scale(area / (state.density * state.sound_speed * state.sound_speed))
        * (C64::ONE + f_t.scale(state.gamma - 1.0));
    if y_shunt.re < 0.0 {
        f_t = C64::new(f_t.re, -f_t.im);
        y_shunt = jw.scale(area / (state.density * state.sound_speed * state.sound_speed))
            * (C64::ONE + f_t.scale(state.gamma - 1.0));
    }
    let mut k = (z_series * y_shunt).sqrt();
    if k.im < 0.0 {
        k = k.scale(-1.0);
    }
    let mut zc = (z_series * y_shunt.recip()).sqrt();
    if zc.re < 0.0 {
        zc = zc.scale(-1.0);
    }
    Ok(SegmentWave {
        wavenumber: k,
        characteristic_impedance: zc,
        specific_impedance: zc.scale(area),
        shear_number: rv,
    })
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
        Segment::ToneHole { .. } => {
            unreachable!("tone holes are handled by the chain loop, not the basis builder")
        }
    }
}

/// Series wall law on a short chimney of length `length`.
///
/// Returns `(R, L)` under `e^{-iωt}` (`Z = R - i ω L`). A compact
/// hole is too short for a 2-port wave: this is the same lumped
/// all-regime / Bessel pin the bore uses, never a WideTube
/// refusal (chimneys sit below `r_v = 10`).
fn chimney_series(
    state: &GasState,
    radius: f64,
    length: f64,
    omega: f64,
    loss: LossModel,
) -> Result<(f64, f64), DuctError> {
    let area = core::f64::consts::PI * radius * radius;
    let l0 = state.density * length / area;
    match loss {
        LossModel::Lossless => Ok((0.0, l0)),
        LossModel::WideTube => {
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            let r_s = if rv > 0.0 && rv.is_finite() {
                omega * l0 * core::f64::consts::SQRT_2 / rv
            } else {
                0.0
            };
            Ok((r_s.max(0.0), l0))
        }
        LossModel::AllRegime => {
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            if rv >= MIN_SHEAR_NUMBER {
                chimney_series(state, radius, length, omega, LossModel::WideTube)
            } else {
                let r_s = 8.0 * state.dynamic_viscosity / (area * radius * radius) * length;
                Ok((r_s.max(0.0), l0 * 4.0 / 3.0))
            }
        }
        LossModel::Bessel => {
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            if !(rv > 0.0 && rv.is_finite()) {
                return Ok((0.0, l0));
            }
            let f_v =
                zwikker_kosten_f(rv).map_err(|_| DuctError::BadParameter { what: "Bessel F_v" })?;
            let den = C64::ONE - f_v;
            if den.abs() < 1.0e-18 {
                return chimney_series(state, radius, length, omega, LossModel::AllRegime);
            }
            // z' = iω ρ / (S (1-F)) in the +iω bookkeeping of
            // [`bessel_wave`]; Re is series R per metre, Im/ω is L.
            let z_prime = C64::new(0.0, omega).scale(state.density / area) / den;
            Ok((
                (z_prime.re * length).max(0.0),
                (z_prime.im / omega * length).abs().max(l0),
            ))
        }
    }
}

fn chimney_thermal_g(
    state: &GasState,
    radius: f64,
    compliance: f64,
    omega: f64,
    loss: LossModel,
) -> f64 {
    if matches!(loss, LossModel::Lossless) || !(compliance > 0.0 && omega > 0.0) {
        return 0.0;
    }
    let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
    if !(rv > 0.0 && rv.is_finite()) {
        return 0.0;
    }
    if matches!(loss, LossModel::Bessel) {
        if let Ok(mut f_t) = zwikker_kosten_f(rv * state.prandtl.sqrt()) {
            // Y' = iω C (1 + (γ-1) F_t); Re Y is G per metre.
            // A sign flip of Im F keeps G ≥ 0 (same pin as bessel_wave).
            if f_t.im > 0.0 {
                f_t = C64::new(f_t.re, -f_t.im);
            }
            return (omega * compliance * (state.gamma - 1.0) * (-f_t.im)).max(0.0);
        }
    }
    if rv >= MIN_SHEAR_NUMBER {
        let eps = core::f64::consts::SQRT_2 / rv;
        return (omega * compliance * (state.gamma - 1.0) * eps / state.prandtl.sqrt()).max(0.0);
    }
    let rt = rv * state.prandtl.sqrt();
    ((state.gamma - 1.0) * compliance * omega * (rt * rt / 16.0).min(0.5)).max(0.0)
}

/// Shunt impedance of a tone hole at `omega`.
///
/// The chimney is a short cylinder: OPEN is that run plus a
/// flanged mouth (Dalmont inner matching on `b/a` is extra
/// length; the `0.8216 b` mass lives in the termination, or
/// the Rayleigh piston above `ka = 1`). CLOSED is the same
/// run with a rigid cap. A compact chimney reprints the
/// lumped `L` / `C` plus wall law; a long one carries its
/// own quarter-wave. WideTube on the chimney is AllRegime
/// so a narrow neck does not raise the bore's `r_v` refusal.
///
/// # Errors
/// [`DuctError`] from the chimney line or a bad radius.
pub fn tone_hole_shunt(
    state: &GasState,
    hole_radius: f64,
    chimney_height: f64,
    hole_state: HoleState,
    omega: f64,
    loss: LossModel,
    bore_radius: f64,
) -> Result<C64, DuctError> {
    let sigma = hole_sigma(hole_state);
    if sigma > 0.0 && sigma < 1.0 {
        let z_o = tone_hole_shunt(
            state,
            hole_radius,
            chimney_height,
            HoleState::Open,
            omega,
            loss,
            bore_radius,
        )?;
        let z_c = tone_hole_shunt(
            state,
            hole_radius,
            chimney_height,
            HoleState::Closed,
            omega,
            loss,
            bore_radius,
        )?;
        return Ok((z_o.recip().scale(sigma) + z_c.recip().scale(1.0 - sigma)).recip());
    }
    if sigma <= 0.0 {
        return chimney_line_impedance(
            state,
            hole_radius,
            chimney_height,
            omega,
            loss,
            Termination::Closed,
        );
    }
    let inner = side_hole_inner_length(hole_radius, bore_radius);
    chimney_line_impedance(
        state,
        hole_radius,
        chimney_height + inner,
        omega,
        loss,
        Termination::FlangedOpen,
    )
}

/// Input impedance of a one-segment chimney. WideTube is
/// AllRegime so a compact neck never trips the bore floor.
fn chimney_line_impedance(
    state: &GasState,
    radius: f64,
    length: f64,
    omega: f64,
    loss: LossModel,
    termination: Termination,
) -> Result<C64, DuctError> {
    if !(radius > 0.0 && length > 0.0 && radius.is_finite() && length.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "chimney radius and length must be positive and finite",
        });
    }
    let chimney_loss = match loss {
        LossModel::WideTube => LossModel::AllRegime,
        other => other,
    };
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    Ok(input_impedance(&duct, state, omega, chimney_loss, termination)?.impedance)
}

fn mul2(a: [C64; 4], b: [C64; 4]) -> [C64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

/// Compact T-junction: series `Z_s/2`, shunt `Z_h`, series `Z_s/2`.
///
/// `Z_s = −i ω ρ t_s / S` under `e^{-iωt}` with Nederveen
/// [`side_hole_series_length`] (open holes only; a closed pad is
/// still a pure shunt).
fn tone_hole_t_junction(
    state: &GasState,
    hole_radius: f64,
    chimney_height: f64,
    hole_state: HoleState,
    omega: f64,
    loss: LossModel,
    bore_radius: f64,
    extra_series_m: f64,
) -> Result<[C64; 4], DuctError> {
    let zh = tone_hole_shunt(
        state,
        hole_radius,
        chimney_height,
        hole_state,
        omega,
        loss,
        bore_radius,
    )?;
    let shunt = [C64::ONE, C64::ZERO, zh.recip(), C64::ONE];
    let sigma = hole_sigma(hole_state);
    if sigma <= 0.0 {
        return Ok(shunt);
    }
    let ts = side_hole_series_length(hole_radius, bore_radius) * sigma + extra_series_m;
    if ts == 0.0 {
        return Ok(shunt);
    }
    let area = core::f64::consts::PI * bore_radius * bore_radius;
    let zs = C64::new(0.0, -omega * state.density * ts / area);
    let half = zs.scale(0.5);
    let series = [C64::ONE, half, C64::ZERO, C64::ONE];
    Ok(mul2(mul2(series, shunt), series))
}

/// Enough spherical substations that each slice's `|Δr|/r` is
/// about 1/5, clamped to `[2, 12]`. Lossless cones stay one
/// transfer (the exact `e^{±ikx}/x` 2-port).
fn cone_loss_slices(r_in: f64, r_out: f64) -> usize {
    let rmin = r_in.min(r_out).max(1.0e-9);
    let n = ((r_out - r_in).abs() / (0.2 * rmin)).ceil() as usize;
    n.clamp(2, 12)
}

/// Lossy cone as cascaded spherical 2-ports, each with `k, Zc`
/// at the slice's own mid-radius. Mean-radius one-shot is the
/// lossless path only.
fn cone_lossy_matrix(
    inlet_radius: f64,
    outlet_radius: f64,
    length: f64,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    wall: Option<&WallPin>,
) -> Result<([C64; 4], f64), DuctError> {
    let n = cone_loss_slices(inlet_radius, outlet_radius);
    let dx = length / n as f64;
    let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let mut min_rv = f64::INFINITY;
    for i in 0..n {
        let t0 = i as f64 / n as f64;
        let t1 = (i + 1) as f64 / n as f64;
        let ra = inlet_radius + (outlet_radius - inlet_radius) * t0;
        let rb = inlet_radius + (outlet_radius - inlet_radius) * t1;
        let r_loc = f64::midpoint(ra, rb);
        let mut wave = match segment_wave(state, r_loc, omega, loss) {
            Ok(w) => w,
            Err(_) if matches!(loss, LossModel::Bessel) => {
                segment_wave(state, r_loc, omega, LossModel::AllRegime)?
            }
            Err(e) => return Err(e),
        };
        if let Some(w) = wall {
            wave = apply_wall_to_wave(wave, r_loc, omega, w)?;
        }
        min_rv = min_rv.min(wave.shear_number);
        let s = segment_matrix(
            &Segment::Cone {
                inlet_radius: ra,
                outlet_radius: rb,
                length: dx,
            },
            &wave,
        )?;
        m = mul2(m, s);
    }
    Ok((m, min_rv))
}

/// The segment 2-port `[p_in, U_in] = M [p_out, U_out]`, built
/// numerically from the analytic basis (no transcribed matrices).
fn segment_matrix(segment: &Segment, wave: &SegmentWave) -> Result<[C64; 4], DuctError> {
    let length = match *segment {
        Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => length,
        Segment::ToneHole { .. } => {
            unreachable!("tone holes are handled by the chain loop, not the basis builder")
        }
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

/// Compact-mouth radiation impedance of an open termination.
///
/// Public name for the same Levine–Schwinger / flanged fit the TMM
/// already uses. Ideal and closed terminations refuse — they are not
/// radiation loads.
///
/// # Errors
/// [`DuctError::RadiationKaTooLarge`] or [`DuctError::BadParameter`].
pub fn compact_radiation_impedance(
    termination: Termination,
    state: &GasState,
    radius: f64,
    omega: f64,
) -> Result<C64, DuctError> {
    match termination {
        Termination::UnflangedOpen | Termination::FlangedOpen => {
            termination_impedance(termination, state, radius, omega)?
                .0
                .ok_or(DuctError::BadParameter {
                    what: "open radiation load missing impedance",
                })
        }
        Termination::Closed | Termination::IdealOpen => Err(DuctError::BadParameter {
            what: "compact radiation impedance is defined for radiating mouths only",
        }),
    }
}

/// Free-field observer amplitude after compact radiation: spherical
/// spreading times ISO 9613-1 molecular absorption.
///
/// Humidity is explicit. This does **not** add Stokes–Kirchhoff on
/// top of ISO (ISO already includes a classical term).
///
/// # Errors
/// ISO window refusals or a non-positive range.
pub fn absorbed_spherical_pressure(
    mouth_pressure: f64,
    state: &GasState,
    omega: f64,
    range_m: f64,
    relative_humidity: f64,
) -> Result<f64, DuctError> {
    if !(range_m > 0.0 && range_m.is_finite() && mouth_pressure.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "observer range must be positive and pressures finite",
        });
    }
    let alpha = fs_material::iso9613::iso9613_absorption(state, relative_humidity, omega).map_err(
        |err| DuctError::BadParameter {
            what: match err {
                fs_material::MaterialError::Parameters { .. } => {
                    "ISO 9613 humidity/temperature/pressure window refused"
                }
                _ => "ISO 9613 evaluation refused",
            },
        },
    )?;
    Ok(mouth_pressure
        * fs_material::iso9613::range_factor(
            range_m,
            alpha,
            fs_material::iso9613::Spreading::Spherical,
        ))
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
                let (rr, xx) = fs_phs::baffled_piston_impedance(
                    state.density,
                    state.sound_speed,
                    radius,
                    omega,
                    8,
                )
                .map_err(|_| DuctError::BadParameter {
                    what: "baffled piston radiation",
                })?;
                return Ok((Some(C64::new(rr, xx)), ka));
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
    input_impedance_wall(duct, state, omega, loss, termination, None)
}

/// [`input_impedance`] with a locally reacting wall.
///
/// The same [`WallPin`] the ODE shunt uses: `Y'` of the gas plus
/// `2π a / Z'_w` with `Z'_w = r − iωσ + i K/ω` under `e^{-iωt}`.
/// `None` is a rigid wall. Chimneys stay rigid this path.
///
/// # Errors
/// As [`input_impedance`], plus a bad wall pin.
pub fn input_impedance_wall(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: Option<&WallPin>,
) -> Result<DuctResponse, DuctError> {
    if duct.segments.is_empty() {
        return Err(DuctError::EmptyDuct);
    }
    for segment in &duct.segments {
        segment.validate()?;
    }
    // Open-hole sites for Nederveen mutual series (near-field
    // beyond the plane-wave cells already in the chain).
    let mut hole_sites: Vec<(f64, f64, f64)> = Vec::new();
    let mut acc = 0.0;
    for segment in &duct.segments {
        if let Segment::ToneHole {
            hole_radius,
            bore_radius,
            state: hole_state,
            ..
        } = *segment
        {
            let sigma = hole_sigma(hole_state);
            if sigma > 0.0 {
                hole_sites.push((acc, hole_radius * sigma.sqrt(), bore_radius));
            }
        } else {
            acc += match *segment {
                Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => length,
                Segment::ToneHole { .. } => 0.0,
            };
        }
    }
    // Chain matrix input -> outlet.
    let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let mut min_rv = f64::INFINITY;
    let mut hole_i = 0usize;
    for segment in &duct.segments {
        let s = if let Segment::ToneHole {
            hole_radius,
            chimney_height,
            bore_radius,
            state: hole_state,
        } = *segment
        {
            // The chimney is often the narrowest element; its shear
            // number folds into the reported margin, and the same
            // LossModel now sits on the lumped neck (no WideTube
            // refusal — a compact hole is not a 2-port wave).
            if matches!(
                loss,
                LossModel::WideTube | LossModel::AllRegime | LossModel::Bessel
            ) {
                let rv = hole_radius * (state.density * omega / state.dynamic_viscosity).sqrt();
                min_rv = min_rv.min(rv);
            }
            let mut extra = 0.0;
            if hole_sigma(hole_state) > 0.0 {
                let (xi, bi, ai) = hole_sites[hole_i];
                for &(xj, bj, _) in &hole_sites {
                    extra += side_hole_mutual_length(bi, bj, ai, (xj - xi).abs());
                }
                hole_i += 1;
            }
            tone_hole_t_junction(
                state,
                hole_radius,
                chimney_height,
                hole_state,
                omega,
                loss,
                bore_radius,
                extra,
            )?
        } else if let Segment::Cone {
            inlet_radius,
            outlet_radius,
            length,
        } = *segment
            && !matches!(loss, LossModel::Lossless)
        {
            let (s, rv) = cone_lossy_matrix(
                inlet_radius,
                outlet_radius,
                length,
                state,
                omega,
                loss,
                wall,
            )?;
            min_rv = min_rv.min(rv);
            s
        } else {
            let mut wave = segment_wave(state, segment.mean_radius(), omega, loss)?;
            if let Some(w) = wall {
                wave = apply_wall_to_wave(wave, segment.mean_radius(), omega, w)?;
            }
            min_rv = min_rv.min(wave.shear_number);
            segment_matrix(segment, &wave)?
        };
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
    // Z_in = (A Z_L + B)/(C Z_L + D); U_mouth/U_in = 1/(C Z_L + D).
    // Closed (U_mouth = 0) is the C, D limit of Z_in with a zero mouth ratio.
    let (impedance, u_mouth_over_u_in, p_mouth_over_u_in) = match z_load {
        Some(zl) => {
            let denom = m[2] * zl + m[3];
            (
                (m[0] * zl + m[1]) * denom.recip(),
                denom.recip(),
                zl * denom.recip(),
            )
        }
        None => (m[0] * m[2].recip(), C64::ZERO, C64::ZERO),
    };
    Ok(DuctResponse {
        omega,
        impedance,
        min_shear_number: min_rv,
        mouth_ka,
        u_mouth_over_u_in,
        p_mouth_over_u_in,
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
    fn all_regime_keeps_losses_below_the_wide_tube_floor() {
        let state = air();
        // r = 1 mm, f = 40 Hz → rv ≈ 3.3 < 10.
        let wave = segment_wave(
            &state,
            0.001,
            2.0 * core::f64::consts::PI * 40.0,
            LossModel::AllRegime,
        )
        .expect("all-regime");
        assert!(wave.shear_number < MIN_SHEAR_NUMBER);
        assert!(
            wave.wavenumber.im > 0.0,
            "narrow-tube k must decay: {:?}",
            wave.wavenumber
        );
        assert!(wave.characteristic_impedance.re > 0.0);
        assert!(
            segment_wave(
                &state,
                0.001,
                2.0 * core::f64::consts::PI * 40.0,
                LossModel::WideTube
            )
            .is_err()
        );
    }

    #[test]
    fn bessel_wave_is_all_regime_and_not_the_wide_tube_pin() {
        let state = air();
        let omega_wide = 2.0 * core::f64::consts::PI * 800.0;
        let omega_narrow = 2.0 * core::f64::consts::PI * 40.0;
        let bessel_w = segment_wave(&state, 0.0075, omega_wide, LossModel::Bessel).expect("wide");
        let tube = segment_wave(&state, 0.0075, omega_wide, LossModel::WideTube).expect("zk1");
        assert!(bessel_w.shear_number > MIN_SHEAR_NUMBER);
        assert!(bessel_w.wavenumber.im > 0.0);
        assert!(bessel_w.characteristic_impedance.re > 0.0);
        let dk = (bessel_w.wavenumber.im - tube.wavenumber.im).abs();
        assert!(
            dk > 1.0e-6 * tube.wavenumber.im,
            "Bessel must not reprint first-order ZK (ΔIm k = {dk})"
        );
        let bessel_n =
            segment_wave(&state, 0.001, omega_narrow, LossModel::Bessel).expect("narrow");
        assert!(bessel_n.shear_number < MIN_SHEAR_NUMBER);
        assert!(bessel_n.wavenumber.im > 0.0);
        assert!(bessel_n.characteristic_impedance.re > 0.0);
        let poiseuille =
            segment_wave(&state, 0.001, omega_narrow, LossModel::AllRegime).expect("poiseuille");
        assert!(
            (bessel_n.wavenumber.im - poiseuille.wavenumber.im).abs()
                < 0.5 * poiseuille.wavenumber.im.max(1.0e-8),
            "narrow Bessel must sit near the Poiseuille floor"
        );
    }

    #[test]
    fn bessel_cone_is_spherical_not_a_mean_cylinder() {
        let state = air();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let cone = Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.006,
                outlet_radius: 0.018,
                length: 0.34,
            }],
        };
        let cyl = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.012,
                length: 0.34,
            }],
        };
        let z_cone = input_impedance(&cone, &state, omega, LossModel::Bessel, Termination::Closed)
            .expect("cone");
        let z_cyl = input_impedance(&cyl, &state, omega, LossModel::Bessel, Termination::Closed)
            .expect("cyl");
        assert!(z_cone.impedance.re >= -1.0e-9);
        assert!(
            (z_cone.impedance - z_cyl.impedance).abs() > 1.0e-3 * z_cyl.impedance.abs().max(1.0),
            "spherical-wave Bessel cone must not reprint the mean cylinder"
        );
    }

    #[test]
    fn cone_losses_use_local_radius_not_the_mean() {
        let state = air();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let cone = Segment::Cone {
            inlet_radius: 0.002,
            outlet_radius: 0.016,
            length: 0.30,
        };
        let duct = Duct {
            segments: vec![cone],
        };
        let z_local = input_impedance(&duct, &state, omega, LossModel::Bessel, Termination::Closed)
            .expect("local")
            .impedance;
        let wave = segment_wave(
            &state,
            f64::midpoint(0.002, 0.016),
            omega,
            LossModel::Bessel,
        )
        .expect("mean wave");
        let m = segment_matrix(&cone, &wave).expect("mean M");
        let z_mean = m[0] * m[2].recip();
        assert!(
            (z_local - z_mean).abs() > 1.0e-3 * z_mean.abs().max(1.0),
            "lossy cone must not reprint the mean-radius 2-port ({z_local:?} vs {z_mean:?})"
        );
        assert!(z_local.re >= -1.0e-9);
        let z_ll = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::Closed,
        )
        .expect("ll")
        .impedance;
        let wave_ll = segment_wave(
            &state,
            f64::midpoint(0.002, 0.016),
            omega,
            LossModel::Lossless,
        )
        .expect("ll wave");
        let m_ll = segment_matrix(&cone, &wave_ll).expect("ll M");
        let z_ll_mean = m_ll[0] * m_ll[2].recip();
        assert!(
            (z_ll - z_ll_mean).abs() < 1.0e-9 * z_ll_mean.abs().max(1.0),
            "lossless cone must stay the exact spherical 2-port"
        );
    }

    #[test]
    fn flanged_piston_lifts_the_ka_ceiling() {
        let state = air();
        let radius = 0.05;
        let omega = 2.0 * core::f64::consts::PI * 2_000.0;
        let ka = omega / state.sound_speed * radius;
        assert!(ka > 1.0, "fixture must sit above the compact-fit ceiling");
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.20,
            }],
        };
        assert!(
            input_impedance(
                &duct,
                &state,
                omega,
                LossModel::Lossless,
                Termination::UnflangedOpen
            )
            .is_err(),
            "unflanged must still refuse ka > 1"
        );
        let z = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::FlangedOpen,
        )
        .expect("flanged piston")
        .impedance;
        assert!(
            z.re > 0.0 && z.re.is_finite() && z.im.is_finite(),
            "Rayleigh piston must stay passive above ka = 1 ({z:?})"
        );
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
            assert_eq!(
                closed.u_mouth_over_u_in,
                C64::ZERO,
                "closed termination has zero mouth flow"
            );
            assert_eq!(closed.p_mouth_over_u_in, C64::ZERO);
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
            assert!(
                open.u_mouth_over_u_in.abs() > 0.0,
                "ideal-open mouth flow ratio must be live"
            );
            assert_eq!(
                open.p_mouth_over_u_in,
                C64::ZERO,
                "ideal-open Z_L = 0 so mouth pressure vanishes"
            );
        }
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"lossless-cylinder-closed-forms\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn quarter_wave_peaks_carry_the_end_correction() {
        use core::fmt::Write as _;
        // Two separated claims (first run conflated them): (a) the
        // LOSSLESS peak ladder of an unflanged open pipe lands on
        // f_n = (2n - 1) c / (4 (L + 0.6133 a)) — the end correction
        // lives in the radiation load; (b) the VISCOTHERMAL ladder
        // sits below it by the independently computed dispersion
        // deficit delta = (1 + (gamma-1)/sqrt(Pr))/(sqrt2 rv) (measured
        // 1.62% at this bore/frequency vs 1.66% predicted — the first
        // executed run flagged exactly this).
        let state = air();
        let (radius, length) = (0.0075, 0.5);
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let l_eff = length + 0.6133 * radius;
        let f_base = state.sound_speed / (4.0 * l_eff);
        let run = |loss: LossModel| -> Vec<f64> {
            let sweep = impedance_sweep(
                &duct,
                &state,
                2.0 * core::f64::consts::PI * 0.3 * f_base,
                2.0 * core::f64::consts::PI * 5.6 * f_base,
                12_000,
                loss,
                Termination::UnflangedOpen,
            )
            .expect("sweep");
            impedance_peaks(&sweep)
                .iter()
                .map(|&i| sweep[i].omega / (2.0 * core::f64::consts::PI))
                .collect()
        };
        let lossless = run(LossModel::Lossless);
        let lossy = run(LossModel::WideTube);
        assert!(lossless.len() >= 3 && lossy.len() >= 3);
        let mut rows = String::new();
        for n in 0..3 {
            let f_pred = (2.0 * n as f64 + 1.0) * f_base;
            let rel = (lossless[n] / f_pred - 1.0).abs();
            // The 0.6133a fit itself is weakly ka-dependent; 0.5% holds
            // through the third peak on this bore.
            assert!(
                rel < 5e-3,
                "lossless peak {n}: {:.2} Hz vs corrected {f_pred:.2} Hz (rel {rel:.4})",
                lossless[n]
            );
            // Viscothermal deficit vs the independent dispersion
            // formula at the peak frequency.
            let omega = 2.0 * core::f64::consts::PI * lossless[n];
            let rv = radius * (state.density * omega / state.dynamic_viscosity).sqrt();
            let delta = (1.0 + (state.gamma - 1.0) / state.prandtl.sqrt())
                / (core::f64::consts::SQRT_2 * rv);
            let measured_deficit = (lossless[n] - lossy[n]) / lossless[n];
            let ratio = measured_deficit / (delta / (1.0 + delta));
            assert!(
                (0.8..1.2).contains(&ratio),
                "peak {n} dispersion deficit {measured_deficit:.4} vs predicted {:.4}",
                delta / (1.0 + delta)
            );
            write!(
                rows,
                "{}{{\"n\":{n},\"lossless\":{:.2},\"lossy\":{:.2},\"deficit_ratio\":{ratio:.3}}}",
                if n == 0 { "" } else { "," },
                lossless[n],
                lossy[n]
            )
            .expect("write");
        }
        // The corrected and UNcorrected ladders stay distinguishable.
        let f_uncorrected = state.sound_speed / (4.0 * length);
        assert!(
            (lossless[0] / f_uncorrected - 1.0).abs() > 5e-3,
            "end correction must be observable against the geometric length"
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

#[cfg(test)]
mod review_regressions {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    #[test]
    fn zc_and_k_match_the_independent_sqrt_route() {
        // REVIEW FINDING (executed): every mutation of the eps_z
        // impedance correction survived the physics oracles — the
        // Zc correction was observationally decorative. This oracle
        // pins BOTH k and Zc through the independent transmission-line
        // route: series impedance with the viscous boundary-layer
        // fraction F_v = sqrt2 (1+i)/rv, shunt admittance with the
        // thermal fraction F_t = F_v / sqrt(Pr), then
        // k_ref = -i sqrt(Zs Yp), z_ref = sqrt(Zs / Yp) via COMPLEX
        // SQUARE ROOTS — different algebra and arithmetic from the
        // implementation's first-order eps form. Agreement must be
        // second-order small; the zeroed / conjugated / thermal-flipped
        // eps_z mutants sit far outside the band.
        let state = air();
        for &(radius, f) in &[(0.0075f64, 200.0f64), (0.0075, 1000.0), (0.02, 500.0)] {
            let omega = 2.0 * core::f64::consts::PI * f;
            let wave = segment_wave(&state, radius, omega, LossModel::WideTube).expect("wave");
            let area = core::f64::consts::PI * radius * radius;
            let rv = wave.shear_number;
            let f_v = C64::new(1.0, 1.0).scale(core::f64::consts::SQRT_2 / rv);
            let f_t = f_v.scale(1.0 / state.prandtl.sqrt());
            // Zs = (i omega rho / S)(1 + F_v);
            // Yp = (i omega S / rho c^2)(1 + (gamma - 1) F_t).
            let i_omega = C64::new(0.0, omega);
            let zs = i_omega.scale(state.density / area) * (C64::ONE + f_v);
            let yp = i_omega.scale(area / (state.density * state.sound_speed * state.sound_speed))
                * (C64::ONE + f_t.scale(state.gamma - 1.0));
            // k = -i sqrt(Zs Yp) with the PHYSICAL branch selected
            // explicitly (Zs Yp sits just below the negative real
            // axis, so the principal square root can land on either
            // side): pick Re k > 0 (forward propagation), and for the
            // impedance pick Re Z > 0 (passive).
            let mut k_ref = C64::new(0.0, -1.0) * (zs * yp).sqrt();
            if k_ref.re < 0.0 {
                k_ref = k_ref.scale(-1.0);
            }
            let mut z_ref = (zs * yp.recip()).sqrt();
            if z_ref.re < 0.0 {
                z_ref = z_ref.scale(-1.0);
            }
            let eps = 1.0 / (core::f64::consts::SQRT_2 * rv);
            // First-order forms differ from the sqrt route at O(eps^2).
            let tol = 6.0 * eps * eps;
            let k_rel = (wave.wavenumber - k_ref).abs() / k_ref.abs();
            let z_rel = (wave.characteristic_impedance - z_ref).abs() / z_ref.abs();
            assert!(
                k_rel < tol,
                "k vs sqrt route at r={radius}, f={f}: rel {k_rel:.3e} (tol {tol:.3e})"
            );
            assert!(
                z_rel < tol,
                "Zc vs sqrt route at r={radius}, f={f}: rel {z_rel:.3e} (tol {tol:.3e})"
            );
            // The review's surviving mutants must violate the band:
            // zeroed eps_z, conjugated eps_z, thermal-flipped eps_z.
            let z0 = state.density * state.sound_speed / area;
            let thermal = (state.gamma - 1.0) / state.prandtl.sqrt();
            let scale = 1.0 / (core::f64::consts::SQRT_2 * rv);
            let mutants = [
                C64::from_re(z0),
                (C64::ONE + C64::new(1.0, -1.0).scale(scale * (1.0 - thermal))).scale(z0),
                (C64::ONE + C64::new(1.0, 1.0).scale(scale * (1.0 + thermal))).scale(z0),
            ];
            for (idx, mutant) in mutants.iter().enumerate() {
                let rel = (*mutant - z_ref).abs() / z_ref.abs();
                assert!(
                    rel > tol,
                    "eps_z mutant {idx} must sit outside the sqrt-route band: {rel:.3e}"
                );
            }
        }
        println!("{{\"suite\":\"fs-duct\",\"case\":\"zc-sqrt-route-pin\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn contracting_cones_are_correct_not_just_claimed() {
        // REVIEW FINDING: contracting cones (negative taper, negative
        // apex coordinate) were handled by the algebra but never
        // tested. Two pins: (a) port reversal — the contracting cone's
        // 2-port equals the expanding cone's with A/D swapped
        // ([[D, B], [C, A]] for a reciprocal det-1 2-port), (b) a
        // CLOSED contracting cone at low frequency matches the lumped
        // cavity-compliance impedance Z = rho c^2 / (-i omega V) with
        // V the frustum volume.
        let state = air();
        let omega = 2.0 * core::f64::consts::PI * 300.0;
        let wave = segment_wave(&state, 0.015, omega, LossModel::Lossless).expect("wave");
        let expanding = Segment::Cone {
            inlet_radius: 0.01,
            outlet_radius: 0.02,
            length: 0.2,
        };
        let contracting = Segment::Cone {
            inlet_radius: 0.02,
            outlet_radius: 0.01,
            length: 0.2,
        };
        let me = segment_matrix(&expanding, &wave).expect("expanding");
        let mc = segment_matrix(&contracting, &wave).expect("contracting");
        // Port reversal of a det-1 2-port swaps the diagonal.
        let expected = [me[3], me[1], me[2], me[0]];
        for i in 0..4 {
            assert!(
                (mc[i] - expected[i]).abs() < 1e-10 * expected[i].abs().max(1e-9),
                "port-reversal entry {i}: {:?} vs {:?}",
                mc[i],
                expected[i]
            );
        }
        // Lumped compliance limit at 20 Hz.
        let (r1, r2, length) = (0.02, 0.01, 0.2);
        let duct = Duct {
            segments: vec![Segment::Cone {
                inlet_radius: r1,
                outlet_radius: r2,
                length,
            }],
        };
        let omega_low = 2.0 * core::f64::consts::PI * 20.0;
        let z = input_impedance(
            &duct,
            &state,
            omega_low,
            LossModel::Lossless,
            Termination::Closed,
        )
        .expect("closed cone")
        .impedance;
        let volume = core::f64::consts::PI * length / 3.0 * (r1 * r1 + r1 * r2 + r2 * r2);
        // Z = rho c^2/(-i omega V) = +i rho c^2/(omega V).
        let z_expected = C64::new(
            0.0,
            state.density * state.sound_speed * state.sound_speed / (omega_low * volume),
        );
        let rel = (z - z_expected).abs() / z_expected.abs();
        assert!(
            rel < 2e-3,
            "closed contracting cone must reduce to cavity compliance: {z:?} vs {z_expected:?} \
             (rel {rel:.2e})"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"contracting-cone-pins\",\"compliance_rel\":{rel:.2e},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn flanged_ladder_carries_the_0p8216_correction() {
        // REVIEW FINDING: the flanged termination had zero executable
        // backing. Same lossless quarter-wave pattern as the unflanged
        // test, on the 0.8216 a corrected length — and the two ladders
        // must be mutually distinguishable.
        let state = air();
        let (radius, length) = (0.0075, 0.5);
        let duct = Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        };
        let run = |termination: Termination| -> f64 {
            let l_eff = match termination {
                Termination::FlangedOpen => length + 0.8216 * radius,
                _ => length + 0.6133 * radius,
            };
            let f_base = state.sound_speed / (4.0 * l_eff);
            let sweep = impedance_sweep(
                &duct,
                &state,
                2.0 * core::f64::consts::PI * 0.5 * f_base,
                2.0 * core::f64::consts::PI * 1.5 * f_base,
                8000,
                LossModel::Lossless,
                termination,
            )
            .expect("sweep");
            let peaks = impedance_peaks(&sweep);
            sweep[peaks[0]].omega / (2.0 * core::f64::consts::PI)
        };
        let f_flanged = run(Termination::FlangedOpen);
        let f_pred = state.sound_speed / (4.0 * (length + 0.8216 * radius));
        let rel = (f_flanged / f_pred - 1.0).abs();
        assert!(
            rel < 5e-3,
            "flanged peak {f_flanged:.2} vs corrected {f_pred:.2} (rel {rel:.4})"
        );
        let f_unflanged = run(Termination::UnflangedOpen);
        assert!(
            f_flanged < f_unflanged,
            "the larger flanged correction must sit lower: {f_flanged:.2} vs {f_unflanged:.2}"
        );
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"flanged-ladder\",\"f\":{f_flanged:.2},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn chained_halves_equal_the_whole_through_input_impedance() {
        // REVIEW FINDING: no end-to-end chain oracle. Two half
        // cylinders of the same radius must reproduce the single long
        // cylinder exactly (identical mean radius, so even the lossy
        // arm is exact).
        let state = air();
        let whole = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0075,
                length: 0.5,
            }],
        };
        let halves = Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: 0.0075,
                    length: 0.25,
                },
                Segment::Cylinder {
                    radius: 0.0075,
                    length: 0.25,
                },
            ],
        };
        for loss in [LossModel::Lossless, LossModel::WideTube] {
            for &f in &[150.0, 440.0, 900.0] {
                let omega = 2.0 * core::f64::consts::PI * f;
                let a = input_impedance(&whole, &state, omega, loss, Termination::UnflangedOpen)
                    .expect("whole")
                    .impedance;
                let b = input_impedance(&halves, &state, omega, loss, Termination::UnflangedOpen)
                    .expect("halves")
                    .impedance;
                assert!(
                    (a - b).abs() < 1e-10 * a.abs(),
                    "chain exactness at {f} Hz ({loss:?}): {a:?} vs {b:?}"
                );
            }
        }
        println!("{{\"suite\":\"fs-duct\",\"case\":\"chain-exactness\",\"verdict\":\"pass\"}}");
    }
}

#[cfg(test)]
mod tone_hole_tests {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};

    fn air20() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    /// The Ernoult et al. 2021 four-hole cylinder (Acta Acustica 5:47,
    /// CC-BY 4.0, Table 1 caliper-measured geometry): main bore
    /// r = 2 mm, L = 287.5 mm; holes at 100/130/180/240 mm with radii
    /// 1.5/1.75/1.75/1.25 mm and chimney heights 1.7/1.3/1.5/1.4 mm.
    /// A FINGERING IS DATA: the instrument description composes the
    /// generic facility (program doctrine).
    fn ernoult_duct(states: [HoleState; 4]) -> Duct {
        let bore = 2.0e-3;
        let hole = |i: usize| -> Segment {
            let (r, h) = [
                (1.5e-3, 1.7e-3),
                (1.75e-3, 1.3e-3),
                (1.75e-3, 1.5e-3),
                (1.25e-3, 1.4e-3),
            ][i];
            Segment::ToneHole {
                hole_radius: r,
                chimney_height: h,
                bore_radius: bore,
                state: states[i],
            }
        };
        let cyl = |length: f64| Segment::Cylinder {
            radius: bore,
            length,
        };
        Duct {
            segments: vec![
                cyl(0.100),
                hole(0),
                cyl(0.030),
                hole(1),
                cyl(0.050),
                hole(2),
                cyl(0.060),
                hole(3),
                cyl(0.0475),
            ],
        }
    }

    fn first_peak_hz(duct: &Duct, lo_hz: f64, hi_hz: f64) -> f64 {
        let state = air20();
        let sweep = impedance_sweep(
            duct,
            &state,
            2.0 * core::f64::consts::PI * lo_hz,
            2.0 * core::f64::consts::PI * hi_hz,
            12_000,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("sweep");
        let peaks = impedance_peaks(&sweep);
        assert!(!peaks.is_empty(), "no peak in [{lo_hz}, {hi_hz}]");
        sweep[peaks[0]].omega / (2.0 * core::f64::consts::PI)
    }

    #[test]
    fn ernoult_2021_measured_fingering_ladder() {
        use core::fmt::Write as _;
        // VALIDATION AGAINST MEASURED DATA (the bead's final clause):
        // first impedance-peak frequencies of the five fingerings of
        // the Ernoult 2021 four-hole cylinder, measured by the paper's
        // two-microphone method with stated +-2 cent peak accuracy
        // (geometry from the CC-BY paper's Table 1; peak values
        // extracted from the openwind-published measured curves,
        // GPLv3, sessions Measure1-3 agree within ~3 Hz: 283/332/449/619/770 Hz).
        // Authored envelope: 30 cents. The old systematic -8..-20
        // cent flat bias was the isolated-disc inner length; Dalmont
        // matching on b/a removes it (xoxx sits at +0.3 cents, inside
        // the paper's ±2 cent peak accuracy). Do not require a flat
        // sign — that encoded the missing term.
        use HoleState::{Closed as X, Open as O};
        let cases: [([HoleState; 4], f64, &str); 5] = [
            ([X, X, X, X], 283.0, "xxxx"),
            ([X, X, X, O], 332.0, "xxxo"),
            ([X, X, O, X], 449.0, "xxox"),
            ([X, O, X, X], 619.0, "xoxx"),
            ([O, X, X, X], 770.0, "oxxx"),
        ];
        let mut rows = String::new();
        let mut previous = 0.0f64;
        for (i, (states, measured, name)) in cases.iter().enumerate() {
            let f = first_peak_hz(&ernoult_duct(*states), 150.0, 1000.0);
            let cents = 1200.0 * (f / measured).ln() / core::f64::consts::LN_2;
            assert!(
                cents.abs() < 30.0,
                "{name}: {f:.1} Hz vs measured {measured} Hz = {cents:+.0} cents"
            );
            assert!(
                cents.abs() < 10.0,
                "{name}: after Dalmont inner matching, peaks stay inside 10 cents ({cents:+.1})"
            );
            // Opening successive holes from the bell up must RAISE the
            // pitch monotonically — the fingering ladder is the
            // instrument-as-data doctrine test.
            assert!(
                f > previous,
                "fingering ladder must rise monotonically: {f:.1} after {previous:.1}"
            );
            previous = f;
            write!(
                rows,
                "{}{{\"fingering\":\"{name}\",\"model_hz\":{f:.1},\"measured_hz\":{measured},\"cents\":{cents:.0}}}",
                if i == 0 { "" } else { "," }
            )
            .expect("write");
        }
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"ernoult-2021-fingering-ladder\",\"citation\":\"Ernoult, Chabassier, Rodriguez, Humeau, Acta Acustica 5:47 (2021), CC-BY-4.0, Table 1; measured curves openwind (GPLv3)\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one coherent algebra + contrast + refusal battery
    fn tone_hole_algebra_and_contrasts() {
        // (a) Exact-cascade pin: input_impedance through a hole equals
        // the hand-composed shunt cascade to 1e-12 — the T-junction
        // wiring cannot drift. (b) Opening a hole RAISES the first
        // resonance (shortens the bore). (c) A closed small hole is a
        // tiny perturbation. (d) Refusals: hole >= bore; high-ka
        // open hole uses the flanged piston (no compact-kb refuse).
        // ceiling. (e) Bitwise determinism.
        let state = air20();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let bore = 7.5e-3;
        let up = Segment::Cylinder {
            radius: bore,
            length: 0.25,
        };
        let down = Segment::Cylinder {
            radius: bore,
            length: 0.25,
        };
        let hole = Segment::ToneHole {
            hole_radius: 3.0e-3,
            chimney_height: 2.0e-3,
            bore_radius: bore,
            state: HoleState::Open,
        };
        let with_hole = Duct {
            segments: vec![up, hole, down],
        };
        let z_full = input_impedance(
            &with_hole,
            &state,
            omega,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("with hole")
        .impedance;
        // Hand composition: Z_down through the downstream tube, then
        // parallel with the shunt, then back through the upstream tube.
        let z_down = input_impedance(
            &Duct {
                segments: vec![down],
            },
            &state,
            omega,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("down")
        .impedance;
        let shunt = tone_hole_shunt(
            &state,
            3.0e-3,
            2.0e-3,
            HoleState::Open,
            omega,
            LossModel::WideTube,
            bore,
        )
        .expect("shunt");
        let ts = side_hole_series_length(3.0e-3, bore);
        let zs = C64::new(
            0.0,
            -omega * state.density * ts / (core::f64::consts::PI * bore * bore),
        );
        let z_after = zs.scale(0.5) + z_down;
        let z_par = (z_after.recip() + shunt.recip()).recip();
        let z_at = zs.scale(0.5) + z_par;
        // Push the T-junction load through the upstream tube's 2-port.
        let wave = segment_wave(&state, bore, omega, LossModel::WideTube).expect("wave");
        let m = segment_matrix(&up, &wave).expect("matrix");
        let z_hand = (m[0] * z_at + m[1]) * (m[2] * z_at + m[3]).recip();
        assert!(
            (z_full - z_hand).abs() < 1e-12 * z_hand.abs(),
            "cascade algebra: {z_full:?} vs {z_hand:?}"
        );
        // (b) open vs closed contrast on the first peak.
        let peak = |hole_state: HoleState| -> f64 {
            let duct = Duct {
                segments: vec![
                    up,
                    Segment::ToneHole {
                        hole_radius: 3.0e-3,
                        chimney_height: 2.0e-3,
                        bore_radius: bore,
                        state: hole_state,
                    },
                    down,
                ],
            };
            first_peak_hz(&duct, 80.0, 400.0)
        };
        let f_open = peak(HoleState::Open);
        let f_closed = peak(HoleState::Closed);
        assert!(
            f_open > 1.15 * f_closed,
            "opening the hole must raise the resonance: {f_open:.1} vs {f_closed:.1}"
        );
        // (c) closed SMALL hole is a tiny perturbation vs no hole.
        let f_plain = first_peak_hz(
            &Duct {
                segments: vec![Segment::Cylinder {
                    radius: bore,
                    length: 0.5,
                }],
            },
            80.0,
            400.0,
        );
        let f_small_closed = {
            let duct = Duct {
                segments: vec![
                    up,
                    Segment::ToneHole {
                        hole_radius: 1.0e-3,
                        chimney_height: 1.0e-3,
                        bore_radius: bore,
                        state: HoleState::Closed,
                    },
                    down,
                ],
            };
            first_peak_hz(&duct, 80.0, 400.0)
        };
        assert!(
            (f_small_closed / f_plain - 1.0).abs() < 5e-3,
            "a closed pinhole must barely perturb: {f_small_closed:.2} vs {f_plain:.2}"
        );
        // (d) refusals.
        assert!(matches!(
            input_impedance(
                &Duct {
                    segments: vec![Segment::ToneHole {
                        hole_radius: 8.0e-3,
                        chimney_height: 1.0e-3,
                        bore_radius: bore,
                        state: HoleState::Open,
                    }],
                },
                &state,
                omega,
                LossModel::WideTube,
                Termination::Closed,
            ),
            Err(DuctError::BadParameter { .. })
        ));
        let z_hi = tone_hole_shunt(
            &state,
            6.0e-3,
            1.0e-3,
            HoleState::Open,
            2.0e5,
            LossModel::Lossless,
            0.01,
        )
        .expect("flanged piston on a hole above ka = 1");
        assert!(
            z_hi.re > 0.0 && z_hi.re.is_finite(),
            "a high-ka hole must stay passive via the Rayleigh piston"
        );
        // (e) determinism.
        let a = input_impedance(
            &with_hole,
            &state,
            omega,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("a")
        .impedance;
        assert_eq!(a.re.to_bits(), z_full.re.to_bits());
        assert_eq!(a.im.to_bits(), z_full.im.to_bits());
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"tone-hole-algebra\",\"f_open\":{f_open:.1},\"f_closed\":{f_closed:.1},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn chimney_wall_law_is_the_bore_pin() {
        let state = air20();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let inv = tone_hole_shunt(
            &state,
            1.5e-3,
            1.7e-3,
            HoleState::Open,
            omega,
            LossModel::Lossless,
            2.0e-3,
        )
        .expect("inviscid");
        let all = tone_hole_shunt(
            &state,
            1.5e-3,
            1.7e-3,
            HoleState::Open,
            omega,
            LossModel::AllRegime,
            2.0e-3,
        )
        .expect("all-regime");
        let bessel = tone_hole_shunt(
            &state,
            1.5e-3,
            1.7e-3,
            HoleState::Open,
            omega,
            LossModel::Bessel,
            2.0e-3,
        )
        .expect("bessel");
        assert!(
            all.re > inv.re * 1.05,
            "narrow chimney must add Poiseuille R ({:.4} vs {:.4})",
            all.re,
            inv.re
        );
        assert!(
            bessel.re > inv.re * 1.05,
            "Bessel F on the chimney must add series R ({:.4} vs {:.4})",
            bessel.re,
            inv.re
        );
        let closed_inv = tone_hole_shunt(
            &state,
            1.5e-3,
            1.7e-3,
            HoleState::Closed,
            omega,
            LossModel::Lossless,
            2.0e-3,
        )
        .expect("closed inviscid");
        let closed_all = tone_hole_shunt(
            &state,
            1.5e-3,
            1.7e-3,
            HoleState::Closed,
            omega,
            LossModel::AllRegime,
            2.0e-3,
        )
        .expect("closed all-regime");
        assert!(
            closed_all.re > 1.0e-6 && closed_inv.re.abs() < 1.0e-12,
            "thermal G on a closed chimney must give Re Z > 0 ({} vs {})",
            closed_all.re,
            closed_inv.re
        );
    }

    #[test]
    fn long_chimney_is_not_a_lumped_inductor() {
        let state = air20();
        let b = 4.0e-3;
        let h = 0.05;
        let bore = 8.0e-3;
        let omega = core::f64::consts::PI * state.sound_speed / (2.0 * h);
        let z_line = tone_hole_shunt(
            &state,
            b,
            h,
            HoleState::Open,
            omega,
            LossModel::Lossless,
            bore,
        )
        .expect("line");
        let t_eff = side_hole_neck_length(h, b, bore);
        let (_, l_h) =
            chimney_series(&state, b, t_eff, omega, LossModel::Lossless).expect("lumped");
        let z_lumped = C64::new(0.0, -omega * l_h);
        assert!(
            (z_line - z_lumped).abs() > 0.5 * z_lumped.abs().max(1.0),
            "a λ/4 chimney must not reprint lumped L ({z_line:?} vs {z_lumped:?})"
        );
        let z_closed = tone_hole_shunt(
            &state,
            b,
            h,
            HoleState::Closed,
            omega,
            LossModel::Lossless,
            bore,
        )
        .expect("closed line");
        let area = core::f64::consts::PI * b * b;
        let c_h = area * h / (state.density * state.sound_speed * state.sound_speed);
        let g = chimney_thermal_g(&state, b, c_h, omega, LossModel::Lossless);
        let z_c_lumped = C64::new(g, -omega * c_h).recip();
        assert!(
            (z_closed - z_c_lumped).abs() > 0.5 * z_c_lumped.abs().max(1.0),
            "a λ/4 closed chimney must not reprint lumped C ({z_closed:?} vs {z_c_lumped:?})"
        );
    }

    #[test]
    fn vent_fraction_is_the_admittance_mix() {
        let state = air20();
        let omega = 2.0 * core::f64::consts::PI * 400.0;
        let z_o = tone_hole_shunt(
            &state,
            3.0e-3,
            2.0e-3,
            HoleState::Open,
            omega,
            LossModel::Lossless,
            7.5e-3,
        )
        .expect("open");
        let z_c = tone_hole_shunt(
            &state,
            3.0e-3,
            2.0e-3,
            HoleState::Closed,
            omega,
            LossModel::Lossless,
            7.5e-3,
        )
        .expect("closed");
        let z_h = tone_hole_shunt(
            &state,
            3.0e-3,
            2.0e-3,
            HoleState::Vent(0.5),
            omega,
            LossModel::Lossless,
            7.5e-3,
        )
        .expect("half");
        let y_mix = z_o.recip().scale(0.5) + z_c.recip().scale(0.5);
        let err = (z_h - y_mix.recip()).abs();
        assert!(
            err < 1.0e-12 * z_h.abs(),
            "Vent(1/2) must be the admittance mix ({z_h:?} vs {:?})",
            y_mix.recip()
        );
        assert!(z_h.re.is_finite() && z_h.re >= 0.0);
    }

    #[test]
    fn locally_reacting_wall_is_not_a_rigid_tmm() {
        let state = air20();
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.012,
                length: 0.34,
            }],
        };
        let omega = 2.0 * core::f64::consts::PI * 220.0;
        let rigid = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::Closed,
        )
        .expect("rigid");
        let soft = WallPin {
            surface_density: 1.5,
            stiffness_per_area: 2.0e5,
            resistance: 0.0,
        };
        let yielding = input_impedance_wall(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::Closed,
            Some(&soft),
        )
        .expect("soft");
        let d_soft = (yielding.impedance - rigid.impedance).abs();
        assert!(
            d_soft > 0.05 * rigid.impedance.abs().max(1.0),
            "a soft wall must move TMM Z_in ({:?} vs {:?})",
            yielding.impedance,
            rigid.impedance
        );
        let stiff = WallPin {
            surface_density: 20.0,
            stiffness_per_area: 1.0e10,
            resistance: 0.0,
        };
        let heavy = input_impedance_wall(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::Closed,
            Some(&stiff),
        )
        .expect("stiff");
        let d_stiff = (heavy.impedance - rigid.impedance).abs();
        assert!(
            d_stiff < 0.05 * d_soft,
            "a stiff wall must sit nearer rigid than a soft wall ({d_stiff} vs {d_soft})"
        );
        let lossy = input_impedance_wall(
            &duct,
            &state,
            omega,
            LossModel::Bessel,
            Termination::Closed,
            Some(&soft),
        )
        .expect("bessel wall");
        assert!(lossy.impedance.re > 0.0);
        assert!(
            input_impedance_wall(
                &duct,
                &state,
                omega,
                LossModel::Lossless,
                Termination::Closed,
                Some(&WallPin {
                    surface_density: 0.0,
                    stiffness_per_area: 2.0e5,
                    resistance: 0.0,
                }),
            )
            .is_err()
        );
    }
}

#[cfg(test)]
mod radiation_and_free_field {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};

    #[test]
    fn compact_radiation_matches_the_tmm_load() {
        let state =
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
        let z =
            compact_radiation_impedance(Termination::FlangedOpen, &state, 0.02, 2.0e3).expect("z");
        assert!(z.re > 0.0 && z.im < 0.0);
        let zu = compact_radiation_impedance(Termination::UnflangedOpen, &state, 0.02, 2.0e3)
            .expect("zu");
        assert!(z.re > zu.re);
        assert!(compact_radiation_impedance(Termination::Closed, &state, 0.02, 2.0e3).is_err());
    }

    #[test]
    fn iso_absorption_on_a_spherical_observer_kills_high_frequency() {
        let state =
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
        let lo = absorbed_spherical_pressure(1.0, &state, 2.0e3, 200.0, 0.50).expect("lo");
        let hi = absorbed_spherical_pressure(1.0, &state, 2.0e4, 200.0, 0.50).expect("hi");
        assert!(hi < lo, "high frequency {hi} must fall below low {lo}");
        assert!(absorbed_spherical_pressure(1.0, &state, 2.0e3, 0.0, 0.50).is_err());
    }
}
