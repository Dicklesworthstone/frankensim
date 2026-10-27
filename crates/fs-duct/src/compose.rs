//! Stable pressure/volume-flow composition, exp(-i omega t). Incoming/outgoing
//! waves use the SAME positive real reference impedance at connected ports.
//! Exact local wave bases contain only exp(i k L), never exp(-i k L).
//! Reflection recursion and Redheffer composition are separate cross-checks.
//! Condition estimates describe feedback cancellation, not global ABCD growth:
//! a perfectly matched arbitrarily attenuating line must remain admissible.
use crate::{
    Duct, DuctError, DuctResponse, HoleState, LoadSpec, LossModel, OpenHoleFlow, Segment,
    SegmentWave, Termination,
};
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_phs::WallPin;
use crate::WallModel;

pub const MAX_COMPOSITION_BLOCKS: usize = 131_072;
pub const MAX_COMPOSITION_MODES: usize = 32;
pub const DEFAULT_MAX_CONDITION: f64 = 1e12;
/// Explicit legacy numerical envelope, independent of the physical pole guard.
pub const LEGACY_MAX_TRANSFER_CONDITION: f64 = 1e8;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Composition {
    ImpedanceRecursion,
    Scattering,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositionMode {
    LegacyAbcd,
    ImpedanceRecursion,
    Scattering,
}
#[derive(Debug, Clone, Copy)]
pub struct ComposeOptions {
    pub rule: Composition,
    pub maximum_condition: f64,
    pub maximum_blocks: usize,
    /// Preserve the old operation graph only inside its declared envelope.
    pub preserve_well_conditioned_bits: bool,
}
impl Default for ComposeOptions {
    fn default() -> Self {
        Self {
            rule: Composition::ImpedanceRecursion,
            maximum_condition: DEFAULT_MAX_CONDITION,
            maximum_blocks: MAX_COMPOSITION_BLOCKS,
            preserve_well_conditioned_bits: true,
        }
    }
}
fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn bound(value: f64, limit: f64) -> Result<(), DuctError> {
    if !limit.is_finite() || limit < 1.0 {
        return Err(bad("composition condition bound must be finite and >=1"));
    }
    if !value.is_finite() || value > limit {
        Err(DuctError::Singular)
    } else {
        Ok(())
    }
}
fn feedback(a: C64, b: C64) -> (C64, f64) {
    let product = a * b;
    let d = C64::ONE - product;
    (d, (1.0 + product.abs()) / d.abs())
}
fn poll(cancelled: &mut impl FnMut() -> bool, step: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled {
            completed_steps: step,
        })
    } else {
        Ok(())
    }
}
impl ComposeOptions {
    fn validate(self) -> Result<(), DuctError> {
        if !self.maximum_condition.is_finite()
            || self.maximum_condition < 1.0
            || self.maximum_blocks == 0
            || self.maximum_blocks > MAX_COMPOSITION_BLOCKS
        {
            return Err(bad(
                "composition needs finite condition >=1 and declared block cap",
            ));
        }
        Ok(())
    }
}

/// Power-wave scattering, row-major [S11,S12,S21,S22]. Physical p=sqrt(R)(a+b)
/// and inward U=(a-b)/sqrt(R); rightward flow has opposite sign at port 2.
/// A supplied feature owns its passivity/evanescent projection qualification.
#[derive(Debug, Clone, Copy)]
pub struct TwoPort {
    pub scattering: [C64; 4],
    pub reference_impedance: f64,
}
impl TwoPort {
    pub fn new(scattering: [C64; 4], reference_impedance: f64) -> Result<Self, DuctError> {
        if scattering.iter().any(|&v| !finite(v))
            || !reference_impedance.is_finite()
            || reference_impedance <= 0.0
        {
            return Err(bad(
                "finite scattering block and positive real reference impedance required",
            ));
        }
        Ok(Self {
            scattering,
            reference_impedance,
        })
    }
    pub fn through(reference_impedance: f64) -> Result<Self, DuctError> {
        Self::new(
            [C64::ZERO, C64::ONE, C64::ONE, C64::ZERO],
            reference_impedance,
        )
    }
    /// Convert an individually resolved finite local ABCD block. Long wave
    /// spans must use from_wave: no conversion can recover an overflowed ABCD.
    pub fn from_abcd(m: [C64; 4], reference: f64) -> Result<Self, DuctError> {
        Self::through(reference)?;
        if m.iter().any(|&v| !finite(v)) {
            return Err(DuctError::Singular);
        }
        let [a, b, c, d] = [
            m[0],
            m[1].scale(1.0 / reference),
            m[2].scale(reference),
            m[3],
        ];
        let denom = a + b + c + d;
        if !finite(denom) || denom.abs() == 0.0 {
            return Err(DuctError::Singular);
        }
        let det = a * d - b * c;
        Self::new(
            [
                (a + b - c - d) / denom,
                det.scale(2.0) / denom,
                C64::from_re(2.0) / denom,
                (-a + b - c + d) / denom,
            ],
            reference,
        )
    }
    /// Derive S from endpoint analytic wave bases. Backward waves are
    /// referenced at the RIGHT endpoint; no growing exponential is evaluated.
    pub fn from_wave(
        segment: &Segment,
        wave: &SegmentWave,
        reference: f64,
    ) -> Result<Self, DuctError> {
        Self::through(reference)?;
        let (rin, rout, length) = match *segment {
            Segment::Cylinder { radius, length } => (radius, radius, length),
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } => (inlet_radius, outlet_radius, length),
            _ => return Err(bad("wave block needs cylinder or spherical cone")),
        };
        segment.validate()?;
        if !finite(wave.wavenumber)
            || wave.wavenumber.im < 0.0
            || !finite(wave.specific_impedance)
            || wave.specific_impedance.abs() == 0.0
            || wave.wavenumber.abs() == 0.0
        {
            return Err(bad(
                "wave block requires finite forward-decaying k and impedance",
            ));
        }
        let e = crate::exp_i(wave.wavenumber.scale(length));
        let (ratio, ylp, ylm, yrp, yrm) = if (rout - rin).abs() < 1e-12 * rin {
            // Derivative of exp(+/-ikx), with the SAME constitutive series
            // relation used by basis_pair. Do not replace the wave owner.
            let area = core::f64::consts::PI * rin * rin;
            let ik = C64::new(0.0, 1.0) * wave.wavenumber;
            let series = (ik * wave.specific_impedance).recip();
            let y = ik.scale(area) * series;
            (1.0, y, -y, y, -y)
        } else {
            let slope = (rout - rin) / length;
            let xl = rin / slope;
            let xr = xl + length;
            let ik = C64::new(0.0, 1.0) * wave.wavenumber;
            let series = (ik * wave.specific_impedance).recip();
            let area_l = core::f64::consts::PI * rin * rin;
            let area_r = core::f64::consts::PI * rout * rout;
            (
                xl / xr,
                (ik - C64::from_re(1.0 / xl)).scale(area_l) * series,
                (-ik - C64::from_re(1.0 / xl)).scale(area_l) * series,
                (ik - C64::from_re(1.0 / xr)).scale(area_r) * series,
                (-ik - C64::from_re(1.0 / xr)).scale(area_r) * series,
            )
        };
        let yl = [ylp.scale(reference), ylm.scale(reference)];
        let yr = [yrp.scale(reference), yrm.scale(reference)];
        // Incoming/outgoing power waves from local forward-left/backward-right
        // pressure amplitudes. The common 1/(2 sqrt(R)) cancels in H G^-1.
        let el = e.scale(1.0 / ratio);
        let er = e.scale(ratio);
        let g = [
            C64::ONE + yl[0],
            el * (C64::ONE + yl[1]),
            er * (C64::ONE - yr[0]),
            C64::ONE - yr[1],
        ];
        let h = [
            C64::ONE - yl[0],
            el * (C64::ONE - yl[1]),
            er * (C64::ONE + yr[0]),
            C64::ONE + yr[1],
        ];
        let scale = g.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale == 0.0 {
            return Err(DuctError::Singular);
        }
        let gs = g.map(|v| v.scale(1.0 / scale));
        let hs = h.map(|v| v.scale(1.0 / scale));
        let lu = lu_complex(&gs, 2).map_err(|_| DuctError::Singular)?;
        let mut s = [C64::ZERO; 4];
        for col in 0..2 {
            let mut x = [C64::ZERO; 2];
            x[col] = C64::ONE;
            lu.solve(&mut x);
            for row in 0..2 {
                s[row * 2 + col] = hs[row * 2] * x[0] + hs[row * 2 + 1] * x[1];
            }
        }
        Self::new(s, reference)
    }
    /// Connect self's right port to right's left port. Solving their two
    /// internal wave equations gives denominator 1-self.S22*right.S11.
    pub fn star(self, right: Self, maximum_condition: f64) -> Result<(Self, f64), DuctError> {
        if self.reference_impedance != right.reference_impedance {
            return Err(bad("connected scattering references differ"));
        }
        let a = self.scattering;
        let b = right.scattering;
        let (d, condition) = feedback(a[3], b[0]);
        bound(condition, maximum_condition)?;
        let s = [
            a[0] + a[1] * b[0] * (a[2] / d),
            a[1] * (b[1] / d),
            b[2] * (a[2] / d),
            b[3] + b[2] * a[3] * (b[1] / d),
        ];
        Ok((Self::new(s, self.reference_impedance)?, condition))
    }
    /// ABCD export only where representable. Network assembly should use
    /// admittance(), which remains finite for an opaque/strongly damped line.
    pub fn abcd(self) -> Result<[C64; 4], DuctError> {
        let [a, b, c, d] = self.scattering;
        let denom = c.scale(2.0);
        let one = C64::ONE;
        let m = [
            ((one + a) * (one - d) + b * c) / denom,
            (((one + a) * (one + d) - b * c) / denom).scale(self.reference_impedance),
            (((one - a) * (one - d) - b * c) / denom).scale(1.0 / self.reference_impedance),
            ((one - a) * (one + d) + b * c) / denom,
        ];
        if m.iter().any(|&v| !finite(v)) {
            Err(DuctError::Singular)
        } else {
            Ok(m)
        }
    }
    pub fn admittance(self, maximum_condition: f64) -> Result<([C64; 4], f64), DuctError> {
        let s = self.scattering;
        let plus = [C64::ONE + s[0], s[1], s[2], C64::ONE + s[3]];
        let (inv, condition) = inverse(&plus, 2, maximum_condition, &mut || false)?;
        let minus = [C64::ONE - s[0], -s[1], -s[2], C64::ONE - s[3]];
        let y = crate::mul2(minus, [inv[0], inv[1], inv[2], inv[3]])
            .map(|v| v.scale(1.0 / self.reference_impedance));
        Ok((y, condition))
    }
}

#[derive(Debug, Clone)]
pub struct ChainResult {
    pub z_in: C64,
    pub u_out_over_u_in: C64,
    pub p_out_over_u_in: C64,
    pub condition_estimate: f64,
    pub scattering: TwoPort,
    /// None when the ABCD-growth indicator cannot be represented. The
    /// terminated physical problem may still be perfectly well-conditioned.
    pub transfer_condition_estimate: Option<f64>,
}
// Log magnitude avoids underflow followed by a later resonant gain. Phases
// are fixed-order products of unit complex numbers; only final output rounds
// to zero below the representable range. No exponential of +alpha*L occurs.
#[derive(Clone, Copy)]
struct Product {
    log_abs: f64,
    phase: C64,
    zero: bool,
}
impl Product {
    fn new() -> Self {
        Self {
            log_abs: 0.0,
            phase: C64::ONE,
            zero: false,
        }
    }
    fn multiply(&mut self, value: C64) -> Result<(), DuctError> {
        let magnitude = value.abs();
        if !magnitude.is_finite() {
            return Err(DuctError::Singular);
        }
        if magnitude == 0.0 {
            self.zero = true;
            return Ok(());
        }
        self.log_abs += det::ln(magnitude);
        self.phase = self.phase * value.scale(1.0 / magnitude);
        self.phase = self.phase.scale(1.0 / self.phase.abs());
        Ok(())
    }
    fn value(self) -> Result<C64, DuctError> {
        let value = if self.zero {
            C64::ZERO
        } else {
            self.phase.scale(det::exp(self.log_abs))
        };
        if finite(value) {
            Ok(value)
        } else {
            Err(DuctError::Singular)
        }
    }
}
/// Load is p/U in Pa*s/m^3; None is a rigid cap. Returning the pressure ratio
/// of a rigid cap is physically meaningful here (the legacy DuctResponse uses
/// zero as its historical closed-mouth radiation sentinel).
pub fn compose_chain(
    blocks: &[TwoPort],
    load: Option<C64>,
    options: ComposeOptions,
) -> Result<ChainResult, DuctError> {
    compose_chain_controlled(blocks, load, options, &mut || false)
}
pub fn compose_chain_controlled(
    blocks: &[TwoPort],
    load: Option<C64>,
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ChainResult, DuctError> {
    let (total, condition) = compose_scattering(blocks, options, &mut *cancelled)?;
    let reference = total.reference_impedance;
    let r_load = match load_reflection(load, reference)? {
        LoadReflection::Finite(r) => r,
        LoadReflection::Pole {
            numerator,
            denominator,
        } => return finish_projective(total, condition, numerator, denominator, options),
    };
    finish_chain(
        blocks, options, cancelled, total, condition, r_load, reference,
    )
}

enum LoadReflection {
    Finite(C64),
    Pole { numerator: C64, denominator: C64 },
}
fn load_reflection(load: Option<C64>, reference: f64) -> Result<LoadReflection, DuctError> {
    let r_load = match load {
        None => C64::ONE,
        Some(z) if finite(z) => {
            let scale = z.re.abs().max(z.im.abs()).max(reference);
            let z = z.scale(1.0 / scale);
            let r = C64::from_re(reference / scale);
            let numerator = z - r;
            let denominator = z + r;
            if denominator == C64::ZERO {
                return Ok(LoadReflection::Pole {
                    numerator,
                    denominator,
                });
            }
            numerator / denominator
        }
        _ => return Err(bad("nonfinite composition load")),
    };
    if !finite(r_load) {
        return Err(DuctError::Singular);
    }
    Ok(LoadReflection::Finite(r_load))
}

/// Terminated input reflection of a chain, referenced to the chain's positive
/// real reference impedance. A lossless resonance drives `|1-r|` to zero, so
/// `Z = R(1+r)/(1-r)` is unrepresentable there while `r` stays finite and
/// exact (DEF-01). This result never forms Z; its feedback refusals are those
/// of [`compose_chain`], without the input-impedance conversion bound.
#[derive(Debug, Clone, Copy)]
pub struct InputReflection {
    pub reflection: C64,
    pub reference_impedance: f64,
    pub condition_estimate: f64,
}
impl InputReflection {
    /// The same inlet wave referenced to another positive real impedance,
    /// through homogeneous `[R(1+r), Zc(1-r)]` coordinates: Z is never formed.
    /// Refuses where that ratio cancels (an active `Z = -Zc`).
    pub fn referenced_to(&self, impedance: f64) -> Result<C64, DuctError> {
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(bad("reflection reference must be finite and positive"));
        }
        let a = (C64::ONE + self.reflection).scale(self.reference_impedance);
        let b = (C64::ONE - self.reflection).scale(impedance);
        let denominator = a + b;
        bound(
            (a.abs() + b.abs()) / denominator.abs(),
            DEFAULT_MAX_CONDITION,
        )?;
        let r = (a - b) / denominator;
        if finite(r) {
            Ok(r)
        } else {
            Err(DuctError::Singular)
        }
    }
}
/// [`compose_chain`]'s input reflection without the impedance conversion.
/// An active load at the reflection-coordinate pole `Z_L = -R` refuses.
pub fn compose_chain_reflection(
    blocks: &[TwoPort],
    load: Option<C64>,
    options: ComposeOptions,
) -> Result<InputReflection, DuctError> {
    let (total, condition) = compose_scattering(blocks, options, &mut || false)?;
    let reference = total.reference_impedance;
    let LoadReflection::Finite(r_load) = load_reflection(load, reference)? else {
        return Err(DuctError::Singular);
    };
    let (reflection, _, condition) =
        terminated_reflection(blocks, options, &mut || false, total, condition, r_load)?;
    if !finite(reflection) {
        return Err(DuctError::Singular);
    }
    Ok(InputReflection {
        reflection,
        reference_impedance: reference,
        condition_estimate: condition,
    })
}
/// The chain and analytic load of [`crate::input_impedance_wall`], returning
/// the input reflection instead of Z. Use it where Z itself refuses at a
/// lossless pole; elsewhere both describe the same wave.
pub fn input_reflection_wall(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: Option<&WallPin>,
) -> Result<InputReflection, DuctError> {
    input_reflection_with_wall(duct, state, omega, loss, termination, WallModel::from_pin(wall))
}

/// Direct stable reflection with the complete local material-wall response.
pub fn input_reflection_with_wall(
    duct: &Duct, state: &GasState, omega: f64, loss: LossModel,
    termination: Termination, wall: WallModel<'_>,
) -> Result<InputReflection, DuctError> {
    let options = ComposeOptions::default();
    let chain = duct_chain_with_wall(duct, state, omega, loss, wall, options, true, || false)?;
    let radius = duct
        .segments
        .last()
        .ok_or(DuctError::EmptyDuct)?
        .outlet_radius();
    let (z_load, _) = crate::termination_impedance(termination, state, radius, omega)?;
    compose_chain_reflection(&chain.blocks, z_load, options)
}

/// Stable reflection for an explicitly supplied finite mouth impedance.
/// Uses the same local propagation and conditioning as the analytic-load path;
/// a pressure/flow impedance pole need not become a reflection-coordinate pole.
pub fn input_reflection_load_with_wall(
    duct: &Duct, state: &GasState, omega: f64, loss: LossModel,
    load: C64, wall: WallModel<'_>,
) -> Result<InputReflection, DuctError> {
    if !finite(load) { return Err(bad("finite mouth impedance required")); }
    let options = ComposeOptions::default();
    let chain = duct_chain_with_wall(duct, state, omega, loss, wall, options, true, || false)?;
    compose_chain_reflection(&chain.blocks, Some(load), options)
}

// A finite load Z=-R is a reflection-coordinate pole, not necessarily a
// physical pole. Contract its homogeneous [Z-R,Z+R] data directly instead of
// manufacturing an infinity or refusing a valid active algebraic probe.
fn finish_projective(
    total: TwoPort,
    prior: f64,
    n: C64,
    d: C64,
    options: ComposeOptions,
) -> Result<ChainResult, DuctError> {
    let [a, b, c, f] = total.scattering;
    let one = C64::ONE;
    let terms = [(one - a) * d, -(one - a) * f * n, -b * c * n];
    let denominator = terms[0] + terms[1] + terms[2];
    let condition = prior.max(terms.iter().map(|v| v.abs()).sum::<f64>() / denominator.abs());
    bound(condition, options.maximum_condition)?;
    let z = (((one + a) * d + (b * c - (one + a) * f) * n) / denominator)
        .scale(total.reference_impedance);
    let u = c * ((d - n) / denominator);
    let p = (c * ((d + n) / denominator)).scale(total.reference_impedance);
    if [z, u, p].iter().any(|&v| !finite(v)) {
        return Err(DuctError::Singular);
    }
    let indicator = 4.0 / c.norm_sq();
    Ok(ChainResult {
        z_in: z,
        u_out_over_u_in: u,
        p_out_over_u_in: p,
        condition_estimate: condition,
        scattering: total,
        transfer_condition_estimate: indicator.is_finite().then_some(indicator),
    })
}

/// Unloaded two-port composition, suitable for graph/feature stamping without
/// reconstructing ABCD from opposite signed impedance probes.
pub fn compose_scattering(
    blocks: &[TwoPort],
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(TwoPort, f64), DuctError> {
    options.validate()?;
    if blocks.is_empty() {
        return Err(DuctError::EmptyDuct);
    }
    if blocks.len() > options.maximum_blocks {
        return Err(bad("composition block work cap exceeded"));
    }
    let reference = blocks[0].reference_impedance;
    let mut total = TwoPort::through(reference)?;
    let mut condition = 1.0_f64;
    for (i, &b) in blocks.iter().enumerate() {
        poll(cancelled, i)?;
        TwoPort::new(b.scattering, b.reference_impedance)?;
        let (s, c) = total.star(b, options.maximum_condition)?;
        total = s;
        condition = condition.max(c);
    }
    Ok((total, condition))
}
fn finish_chain(
    blocks: &[TwoPort],
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
    total: TwoPort,
    condition: f64,
    r_load: C64,
    reference: f64,
) -> Result<ChainResult, DuctError> {
    let (r_in, transfer, mut condition) =
        terminated_reflection(blocks, options, cancelled, total, condition, r_load)?;
    let denom = C64::ONE - r_in;
    condition = condition.max((1.0 + r_in.abs()) / denom.abs());
    bound(condition, options.maximum_condition)?;
    let z_in = ((C64::ONE + r_in) / denom).scale(reference);
    let u_out_over_u_in = transfer * ((C64::ONE - r_load) / denom);
    let p_out_over_u_in = (transfer * ((C64::ONE + r_load) / denom)).scale(reference);
    if [z_in, u_out_over_u_in, p_out_over_u_in]
        .iter()
        .any(|&v| !finite(v))
    {
        return Err(DuctError::Singular);
    }
    let transmission = total.scattering[2].abs();
    let indicator = 4.0 / (transmission * transmission);
    let transfer_condition_estimate = indicator.is_finite().then_some(indicator);
    Ok(ChainResult {
        z_in,
        u_out_over_u_in,
        p_out_over_u_in,
        condition_estimate: condition,
        scattering: total,
        transfer_condition_estimate,
    })
}
/// Inlet reflection and load-to-inlet wave transfer with every feedback
/// condition bounded; no impedance has been formed yet.
fn terminated_reflection(
    blocks: &[TwoPort],
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
    total: TwoPort,
    mut condition: f64,
    r_load: C64,
) -> Result<(C64, C64, f64), DuctError> {
    let (r_in, transfer) = match options.rule {
        Composition::Scattering => {
            let s = total.scattering;
            let (d, c) = feedback(s[3], r_load);
            condition = condition.max(c);
            bound(condition, options.maximum_condition)?;
            (s[0] + s[1] * r_load * (s[2] / d), s[2] / d)
        }
        Composition::ImpedanceRecursion => {
            let mut r = r_load;
            let mut product = Product::new();
            for (i, block) in blocks.iter().enumerate().rev() {
                poll(cancelled, i)?;
                let s = block.scattering;
                let (d, c) = feedback(s[3], r);
                condition = condition.max(c);
                bound(condition, options.maximum_condition)?;
                let t = s[2] / d;
                product.multiply(t)?;
                r = s[0] + s[1] * r * t;
            }
            (r, product.value()?)
        }
    };
    Ok((r_in, transfer, condition))
}

/// Shared modal propagation: the existing exact multiplication order is kept.
pub fn propagate_reflection(r: &mut [C64], propagators: &[C64]) -> Result<(), DuctError> {
    let n = propagators.len();
    if n == 0
        || n > MAX_COMPOSITION_MODES
        || r.len() != n * n
        || propagators
            .iter()
            .any(|&v| !finite(v) || v.abs() > 1.0 + 8.0 * f64::EPSILON)
        || r.iter().any(|&v| !finite(v))
    {
        return Err(bad("bounded modal propagation shape/decay admission"));
    }
    for i in 0..n {
        for j in 0..n {
            r[i * n + j] = propagators[i] * r[i * n + j] * propagators[j];
        }
    }
    Ok(())
}

fn mul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut out = vec![C64::ZERO; n * n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                out[i * n + j] = out[i * n + j] + a[i * n + k] * b[k * n + j];
            }
        }
    }
    out
}
fn inverse(
    a: &[C64],
    n: usize,
    limit: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    let norm = |m: &[C64]| {
        m.chunks_exact(n)
            .map(|r| r.iter().map(|v| v.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max)
    };
    let lu =
        fs_la::eigen_complex::lu_complex_with_cancel(a, n, &mut *cancelled).map_err(
            |e| match e {
                fs_la::eigen_complex::LuComplexError::Cancelled => {
                    DuctError::Cancelled { completed_steps: 0 }
                }
                _ => DuctError::Singular,
            },
        )?;
    let mut out = vec![C64::ZERO; n * n];
    for col in 0..n {
        poll(cancelled, col)?;
        let mut b = vec![C64::ZERO; n];
        b[col] = C64::ONE;
        lu.solve_with_cancel(&mut b, &mut *cancelled)
            .map_err(|e| match e {
                fs_la::eigen_complex::LuComplexError::Cancelled => DuctError::Cancelled {
                    completed_steps: col,
                },
                _ => DuctError::Singular,
            })?;
        for row in 0..n {
            out[row * n + col] = b[row];
        }
    }
    let c = norm(a) * norm(&out);
    bound(c, limit)?;
    Ok((out, c))
}
/// Generic equal-mode Redheffer product for feature/multimodal owners. Each
/// input is row-major 2n x 2n, all connected power-wave references must match.
/// No passivity or modal-truncation certificate is inferred from supplied S.
pub fn redheffer_star(
    left: &[C64],
    right: &[C64],
    n: usize,
    maximum_condition: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    if n == 0
        || n > MAX_COMPOSITION_MODES
        || left.len() != 4 * n * n
        || right.len() != 4 * n * n
        || left.iter().chain(right).any(|&v| !finite(v))
        || !maximum_condition.is_finite()
        || maximum_condition < 1.0
    {
        return Err(bad("Redheffer product shape/mode/condition cap"));
    }
    let block = |m: &[C64], row: usize, col: usize| -> Vec<C64> {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| m[(row * n + i) * 2 * n + col * n + j]))
            .collect()
    };
    let a = core::array::from_fn::<_, 4, _>(|i| block(left, i / 2, i % 2));
    let b = core::array::from_fn::<_, 4, _>(|i| block(right, i / 2, i % 2));
    let mut feedback = mul(&b[0], &a[3], n);
    let cancellation_scale = 1.0
        + feedback
            .chunks_exact(n)
            .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max);
    for v in &mut feedback {
        *v = -*v;
    }
    for i in 0..n {
        feedback[i * n + i] = feedback[i * n + i] + C64::ONE;
    }
    let (inv, matrix_condition) = inverse(&feedback, n, maximum_condition, &mut cancelled)?;
    let condition = matrix_condition.max(
        cancellation_scale
            * inv
                .chunks_exact(n)
                .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
                .fold(0.0_f64, f64::max),
    );
    bound(condition, maximum_condition)?;
    let x = mul(&inv, &mul(&b[0], &a[2], n), n);
    let y = mul(&inv, &b[1], n);
    let mut s11 = mul(&a[1], &x, n);
    for i in 0..n * n {
        s11[i] = s11[i] + a[0][i];
    }
    let s12 = mul(&a[1], &y, n);
    let mut mid = mul(&a[3], &x, n);
    for i in 0..n * n {
        mid[i] = mid[i] + a[2][i];
    }
    let s21 = mul(&b[2], &mid, n);
    let mut s22 = mul(&b[2], &mul(&a[3], &y, n), n);
    for i in 0..n * n {
        s22[i] = s22[i] + b[3][i];
    }
    let blocks = [s11, s12, s21, s22];
    let mut out = vec![C64::ZERO; 4 * n * n];
    for i in 0..2 * n {
        poll(&mut cancelled, i)?;
        for j in 0..2 * n {
            out[i * 2 * n + j] = blocks[(i / n) * 2 + j / n][(i % n) * n + j % n];
        }
    }
    if out.iter().any(|&v| !finite(v)) {
        return Err(DuctError::Singular);
    }
    Ok((out, condition))
}

/// Exact owner-derived blocks, retaining the supplied cone/profile resolution.
pub struct PreparedChain {
    pub blocks: Vec<TwoPort>,
    pub reference_impedance: f64,
    pub min_shear_number: f64,
}
#[derive(Debug, Clone, Copy)]
pub(crate) struct PreparedOpenHole {
    pub block_index: usize,
    pub axial_position_m: f64,
    pub effective_radius_m: f64,
    pub chimney_length_m: f64,
    pub sigma: f64,
    pub half_series_impedance: C64,
    pub open_chimney_impedance: C64,
    pub chimney_u_mouth_over_u_in: C64,
    pub vortex_resistance_pa_s_m3: f64,
    pub vortex_critical_flow_m3_s: f64,
    pub vortex_quadratic_coeff_pa_s2_m6: f64,
}
fn push_block(blocks: &mut Vec<TwoPort>, limit: usize, block: TwoPort) -> Result<(), DuctError> {
    if blocks.len() >= limit {
        return Err(bad("composition block work cap exceeded"));
    }
    blocks.push(block);
    Ok(())
}
pub fn duct_chain(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    wall: Option<&WallPin>,
    options: ComposeOptions,
    compact_corrections: bool,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, DuctError> {
    duct_chain_with_wall(duct, state, omega, loss, WallModel::from_pin(wall),
        options, compact_corrections, cancelled)
}

/// The same stable cascade, evaluating a local wall law at every gas/geometry sample.
#[allow(clippy::too_many_arguments)]
pub fn duct_chain_with_wall(
    duct: &Duct, state: &GasState, omega: f64, loss: LossModel,
    wall: WallModel<'_>, options: ComposeOptions, compact_corrections: bool,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, DuctError> {
    duct_chain_with_wall_and_holes(
        duct,
        state,
        omega,
        loss,
        wall,
        options,
        compact_corrections,
        cancelled,
    )
    .map(|(chain, _)| chain)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn duct_chain_with_wall_and_holes(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    wall: WallModel<'_>,
    options: ComposeOptions,
    compact_corrections: bool,
    cancelled: impl FnMut() -> bool,
) -> Result<(PreparedChain, Vec<PreparedOpenHole>), DuctError> {
    duct_chain_in_exterior(duct, state, state, omega, loss, wall, options,
        compact_corrections, cancelled)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn duct_chain_in_exterior(
    duct: &Duct,
    state: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    wall: WallModel<'_>,
    options: ComposeOptions,
    compact_corrections: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PreparedChain, Vec<PreparedOpenHole>), DuctError> {
    options.validate()?;
    poll(&mut cancelled, 0)?;
    let first = duct.segments.first().ok_or(DuctError::EmptyDuct)?;
    crate::validate_gas(state)?;
    crate::exterior::admit(state, exterior)?;
    if !omega.is_finite() || omega <= 0.0 {
        return Err(bad(
            "composition requires positive finite gas and angular frequency",
        ));
    }
    for s in &duct.segments {
        s.validate()?;
    }
    let inlet = first.inlet_radius();
    let reference = state.density * state.sound_speed / (core::f64::consts::PI * inlet * inlet);
    let mut blocks = Vec::new();
    let mut open_holes = Vec::new();
    let mut min_shear = f64::INFINITY;
    let (sites, total_acoustic_len_m, outlet_radius_m) = crate::collect_open_hole_sites(duct);
    let mut hole_index = 0;
    for (index, s) in duct.segments.iter().enumerate() {
        poll(&mut cancelled, index)?;
        match *s {
            Segment::Profile { ref profile, steps } => {
                if steps.0 > options.maximum_blocks - blocks.len() {
                    return Err(bad("composition block work cap exceeded"));
                }
                let (local, rv) = crate::profile::profile_blocks_controlled(
                    profile,
                    state,
                    omega,
                    loss,
                    steps,
                    wall,
                    &mut cancelled,
                )?;
                min_shear = min_shear.min(rv);
                for m in local {
                    push_block(
                        &mut blocks,
                        options.maximum_blocks,
                        TwoPort::from_abcd(m, reference)?,
                    )?;
                }
            }
            Segment::ToneHole {
                hole_radius,
                chimney_height,
                bore_radius,
                state: hs,
            } => {
                if matches!(loss, LossModel::WideTube | LossModel::Bessel) {
                    min_shear = min_shear.min(
                        hole_radius * (state.density * omega / state.dynamic_viscosity).sqrt(),
                    );
                }
                let sigma = crate::hole_sigma(hs);
                let block_index = blocks.len();
                let m = if compact_corrections {
                    let mut extra = 0.0;
                    let mut r_vortex = 0.0;
                    if sigma > 0.0 {
                        let site = sites[hole_index];
                        let (xi, bi, ai) = (
                            site.axial_position_m,
                            site.effective_radius_m,
                            site.bore_radius_m,
                        );
                        for site_j in &sites {
                            extra += fs_phs::side_hole_mutual_length(
                                bi,
                                site_j.effective_radius_m,
                                ai,
                                (site_j.axial_position_m - xi).abs(),
                            );
                        }
                        let (mut open_resp, chimney_length_m) =
                            crate::chimney_open_response_with_wall(
                                state,
                                exterior,
                                hole_radius,
                                chimney_height,
                                bore_radius,
                                omega,
                                loss,
                                true,
                                wall,
                            )?;
                        let (vortex_r, vortex_u_crit, vortex_k_quad) =
                            crate::open_hole_vortex_spec_in_exterior(
                                state,
                                exterior,
                                loss,
                                hole_radius,
                                bore_radius,
                                chimney_height,
                                &sites,
                                hole_index,
                                total_acoustic_len_m,
                                outlet_radius_m,
                            );
                        r_vortex = vortex_r;
                        open_resp.impedance = open_resp.impedance + C64::from_re(r_vortex);
                        let ts =
                            fs_phs::side_hole_series_length(hole_radius, bore_radius) * sigma
                                + extra;
                        let half = if ts == 0.0 {
                            C64::ZERO
                        } else {
                            let area = core::f64::consts::PI * bore_radius * bore_radius;
                            C64::new(0.0, -omega * state.density * ts / area).scale(0.5)
                        };
                        open_holes.push(PreparedOpenHole {
                            block_index,
                            axial_position_m: xi,
                            effective_radius_m: bi,
                            chimney_length_m,
                            sigma,
                            half_series_impedance: half,
                            open_chimney_impedance: open_resp.impedance,
                            chimney_u_mouth_over_u_in: open_resp.u_mouth_over_u_in,
                            vortex_resistance_pa_s_m3: vortex_r,
                            vortex_critical_flow_m3_s: vortex_u_crit,
                            vortex_quadratic_coeff_pa_s2_m6: vortex_k_quad,
                        });
                        hole_index += 1;
                    }
                    crate::tone_hole_t_junction_with_wall_vortex(
                        state,
                        exterior,
                        hole_radius,
                        chimney_height,
                        hs,
                        omega,
                        loss,
                        bore_radius,
                        extra,
                        r_vortex,
                        wall,
                    )?
                } else {
                    if sigma > 0.0 {
                        let site = sites[hole_index];
                        let (xi, bi) = (site.axial_position_m, site.effective_radius_m);
                        let (open_resp, chimney_length_m) = crate::chimney_open_response_with_wall(
                            state,
                            exterior,
                            hole_radius,
                            chimney_height,
                            bore_radius,
                            omega,
                            loss,
                            false,
                            wall,
                        )?;
                        open_holes.push(PreparedOpenHole {
                            block_index,
                            axial_position_m: xi,
                            effective_radius_m: bi,
                            chimney_length_m,
                            sigma,
                            half_series_impedance: C64::ZERO,
                            open_chimney_impedance: open_resp.impedance,
                            chimney_u_mouth_over_u_in: open_resp.u_mouth_over_u_in,
                            vortex_resistance_pa_s_m3: 0.0,
                            vortex_critical_flow_m3_s: 0.0,
                            vortex_quadratic_coeff_pa_s2_m6: 0.0,
                        });
                        hole_index += 1;
                    }
                    let loss = if loss == LossModel::WideTube {
                        LossModel::Bessel
                    } else {
                        loss
                    };
                    let impedance = |termination| {
                        crate::chimney_line_impedance_with_wall(
                            state,
                            exterior,
                            hole_radius,
                            chimney_height,
                            omega,
                            loss,
                            termination,
                            wall,
                        )
                    };
                    let y = match hs {
                        HoleState::Open => impedance(Termination::FlangedOpen)?.recip(),
                        HoleState::Closed => impedance(Termination::Closed)?.recip(),
                        HoleState::Vent(sigma) => {
                            impedance(Termination::FlangedOpen)?.recip().scale(sigma)
                                + impedance(Termination::Closed)?.recip().scale(1.0 - sigma)
                        }
                    };
                    [C64::ONE, C64::ZERO, y, C64::ONE]
                };
                push_block(
                    &mut blocks,
                    options.maximum_blocks,
                    TwoPort::from_abcd(m, reference)?,
                )?;
            }
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } if !wall.is_rigid() || loss != LossModel::Lossless => {
                let count = crate::cone_loss_slices(inlet_radius, outlet_radius)?;
                if count > options.maximum_blocks - blocks.len() {
                    return Err(bad("composition block work cap exceeded"));
                }
                let dx = length / count as f64;
                for i in 0..count {
                    poll(&mut cancelled, i)?;
                    let t0 = i as f64 / count as f64;
                    let t1 = (i + 1) as f64 / count as f64;
                    let ra = inlet_radius + (outlet_radius - inlet_radius) * t0;
                    let rb = inlet_radius + (outlet_radius - inlet_radius) * t1;
                    let mut wave = crate::cone_midpoint_wave(state, ra, rb, omega, loss)?;
                    if !wall.is_rigid() {
                        let slope = (rb - ra) / dx;
                        wave = wall.apply(wave, f64::midpoint(ra, rb), omega,
                            (1.0 + slope * slope).sqrt())?;
                    }
                    min_shear = min_shear.min(wave.shear_number);
                    push_block(
                        &mut blocks,
                        options.maximum_blocks,
                        TwoPort::from_wave(
                            &Segment::Cone {
                                inlet_radius: ra,
                                outlet_radius: rb,
                                length: dx,
                            },
                            &wave,
                            reference,
                        )?,
                    )?;
                }
            }
            _ => {
                let mut wave = crate::segment_wave(state, s.mean_radius(), omega, loss)?;
                wave = wall.apply(wave, s.mean_radius(), omega, 1.0)?;
                min_shear = min_shear.min(wave.shear_number);
                push_block(
                    &mut blocks,
                    options.maximum_blocks,
                    TwoPort::from_wave(s, &wave, reference)?,
                )?;
            }
        }
    }
    Ok((
        PreparedChain {
            blocks,
            reference_impedance: reference,
            min_shear_number: min_shear,
        },
        open_holes,
    ))
}

fn open_hole_flows(
    chain: &PreparedChain,
    open_holes: &[PreparedOpenHole],
    z_load: Option<C64>,
    options: ComposeOptions,
) -> Result<Vec<OpenHoleFlow>, DuctError> {
    if open_holes.is_empty() {
        return Ok(Vec::new());
    }
    let reference = chain.reference_impedance;
    let LoadReflection::Finite(r_load) = load_reflection(z_load, reference)? else {
        return Err(DuctError::Singular);
    };
    let n = chain.blocks.len();
    let mut r_after = vec![C64::ZERO; n];
    let mut t_step = vec![C64::ZERO; n];
    let mut r = r_load;
    let mut condition = 1.0_f64;
    for (i, block) in chain.blocks.iter().enumerate().rev() {
        r_after[i] = r;
        let s = block.scattering;
        let (d, c) = feedback(s[3], r);
        condition = condition.max(c);
        bound(condition, options.maximum_condition)?;
        let t = s[2] / d;
        t_step[i] = t;
        r = s[0] + s[1] * r * t;
    }
    let denom = C64::ONE - r;
    condition = condition.max((1.0 + r.abs()) / denom.abs());
    bound(condition, options.maximum_condition)?;
    let mut prefix = Product::new();
    let mut next_block = 0usize;
    let mut left_block = 0usize;
    let mut r_left = C64::ZERO;
    let mut flows = Vec::with_capacity(open_holes.len());
    for hole in open_holes {
        while left_block < hole.block_index {
            let s = chain.blocks[left_block].scattering;
            let (d, c) = feedback(s[0], r_left);
            condition = condition.max(c);
            bound(condition, options.maximum_condition)?;
            r_left = s[3] + s[1] * r_left * (s[2] / d);
            left_block += 1;
        }
        while next_block <= hole.block_index {
            prefix.multiply(t_step[next_block])?;
            next_block += 1;
        }
        let a_right_over_a_in = prefix.value()?;
        let r_right = r_after[hole.block_index];
        let p_mid_over_r0_u_in = a_right_over_a_in
            * ((C64::ONE + r_right)
                + hole
                    .half_series_impedance
                    .scale(1.0 / reference)
                    * (C64::ONE - r_right))
            / denom;
        let u_chimney_in_over_u_in = p_mid_over_r0_u_in
            .scale(hole.sigma * reference)
            / hole.open_chimney_impedance;
        let u_hole_over_u_in = u_chimney_in_over_u_in * hole.chimney_u_mouth_over_u_in;
        if !finite(u_hole_over_u_in) || !finite(u_chimney_in_over_u_in) {
            return Err(DuctError::Singular);
        }
        let shunt_port_reflection = if hole.vortex_resistance_pa_s_m3 > 0.0 {
            let z_left = (C64::ONE + r_left).scale(reference) / (C64::ONE - r_left)
                + hole.half_series_impedance;
            let z_right = (C64::ONE + r_right).scale(reference) / (C64::ONE - r_right)
                + hole.half_series_impedance;
            let z_bore_mid = ((z_left * z_right) / (z_left + z_right)).scale(hole.sigma);
            let r2 = C64::from_re(hole.vortex_resistance_pa_s_m3);
            let z_port2_minus_r2 = (hole.open_chimney_impedance - r2.scale(2.0)) + z_bore_mid;
            let z_port2_plus_r2 = hole.open_chimney_impedance + z_bore_mid;
            let s22 = z_port2_minus_r2 / z_port2_plus_r2;
            if !finite(s22) {
                return Err(DuctError::Singular);
            }
            s22
        } else {
            C64::ZERO
        };
        flows.push(OpenHoleFlow {
            axial_position_m: hole.axial_position_m,
            effective_radius_m: hole.effective_radius_m,
            chimney_length_m: hole.chimney_length_m,
            u_hole_over_u_in,
            u_chimney_in_over_u_in,
            vortex_resistance_pa_s_m3: hole.vortex_resistance_pa_s_m3,
            vortex_critical_flow_m3_s: hole.vortex_critical_flow_m3_s,
            vortex_quadratic_coeff_pa_s2_m6: hole.vortex_quadratic_coeff_pa_s2_m6,
            shunt_port_reflection,
        });
    }
    Ok(flows)
}

pub(crate) fn evaluate(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: LoadSpec,
    wall: Option<&WallPin>,
) -> Result<DuctResponse, DuctError> {
    evaluate_with_wall(duct, state, omega, loss, load, WallModel::from_pin(wall))
}

pub(crate) fn evaluate_with_wall(
    duct: &Duct, state: &GasState, omega: f64, loss: LossModel,
    load: LoadSpec, wall: WallModel<'_>,
) -> Result<DuctResponse, DuctError> {
    evaluate_in_exterior(duct, state, state, omega, loss, load, wall)
}

// The legacy ABCD path has no separate radiation gas. Reusing it after the
// new chain was built would silently discard every changed side-hole load.
fn exterior_options(duct: &Duct, state: &GasState, exterior: &GasState) -> ComposeOptions {
    ComposeOptions { preserve_well_conditioned_bits: state == exterior
        || !duct.segments.iter().any(|segment| matches!(segment,
            Segment::ToneHole { state, .. } if crate::hole_sigma(*state) > 0.0)),
        ..ComposeOptions::default() }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate_in_exterior(
    duct: &Duct, state: &GasState, exterior: &GasState, omega: f64, loss: LossModel,
    load: LoadSpec, wall: WallModel<'_>,
) -> Result<DuctResponse, DuctError> {
    let options = exterior_options(duct, state, exterior);
    let (chain, _) = duct_chain_in_exterior(duct, state, exterior, omega, loss,
        wall, options, true, || false)?;
    evaluate_prepared_chain(duct, state, omega, loss, load, wall, options, &chain)
}

pub(crate) fn evaluate_with_wall_and_holes(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: LoadSpec,
    wall: WallModel<'_>,
) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
    evaluate_in_exterior_and_holes(duct, state, state, omega, loss, load, wall)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate_in_exterior_and_holes(
    duct: &Duct, state: &GasState, exterior: &GasState, omega: f64, loss: LossModel,
    load: LoadSpec, wall: WallModel<'_>,
) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
    let options = exterior_options(duct, state, exterior);
    let (chain, open_holes) =
        duct_chain_in_exterior(duct, state, exterior, omega, loss, wall, options, true, || false)?;
    let radius = duct
        .segments
        .last()
        .ok_or(DuctError::EmptyDuct)?
        .outlet_radius();
    let (z_load, _) = match load {
        LoadSpec::Analytic(t) => crate::termination_impedance(t, state, radius, omega)?,
        LoadSpec::Fixed(z) => (Some(z), omega * radius / state.sound_speed),
    };
    let response = evaluate_prepared_chain(duct, state, omega, loss, load, wall, options, &chain)?;
    let flows = open_hole_flows(&chain, &open_holes, z_load, options)?;
    Ok((response, flows))
}

#[allow(clippy::too_many_arguments)]
fn evaluate_prepared_chain(
    duct: &Duct,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: LoadSpec,
    wall: WallModel<'_>,
    options: ComposeOptions,
    chain: &PreparedChain,
) -> Result<DuctResponse, DuctError> {
    let radius = duct
        .segments
        .last()
        .ok_or(DuctError::EmptyDuct)?
        .outlet_radius();
    let (z_load, mouth_ka) = match load {
        LoadSpec::Analytic(t) => crate::termination_impedance(t, state, radius, omega)?,
        LoadSpec::Fixed(z) => (Some(z), omega * radius / state.sound_speed),
    };
    let solved = compose_chain(&chain.blocks, z_load, options)?;
    if let Some(pin) = wall.legacy_pin()
        && options.preserve_well_conditioned_bits
        && solved.condition_estimate <= LEGACY_MAX_TRANSFER_CONDITION
        && solved
            .transfer_condition_estimate
            .is_some_and(|v| v <= LEGACY_MAX_TRANSFER_CONDITION)
    {
        if let Ok(mut old) = crate::input_impedance_legacy(duct, state, omega, loss, load, pin) {
            if [old.impedance, old.u_mouth_over_u_in, old.p_mouth_over_u_in]
                .iter()
                .all(|&v| finite(v))
            {
                old.condition_estimate = solved.condition_estimate;
                old.composition = CompositionMode::LegacyAbcd;
                return Ok(old);
            }
        }
    }
    Ok(DuctResponse {
        omega,
        impedance: solved.z_in,
        min_shear_number: chain.min_shear_number,
        mouth_ka,
        u_mouth_over_u_in: solved.u_out_over_u_in,
        p_mouth_over_u_in: if z_load.is_none() {
            C64::ZERO
        } else {
            solved.p_out_over_u_in
        },
        condition_estimate: solved.condition_estimate,
        composition: CompositionMode::ImpedanceRecursion,
    })
}

#[cfg(test)]
mod tests;
