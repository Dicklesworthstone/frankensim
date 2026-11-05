//! Continuous resonance analysis in the original power-conjugate coordinates.
//!
//! For peak phasors with `exp(-i omega t)`, solve
//! `[-i omega E - (J-R) W] x = G u`, where `E = diag(I_diff, 0)` and
//! `W = diag(Q_diff, I_alg)`. Algebraic multipliers are EFFORTS, not stored
//! coordinates. They remain in the solve, including at component resonances;
//! neither E nor an isolated component impedance is inverted.
//!
//! The full declared quadratic storage, mutual coupling and relaxation states
//! are retained. No time-step warping, modal truncation, linewidth fitting or
//! numerical damping is introduced. Ordinary constructors require exact linear
//! storage; the explicit equilibrium adapter admits local nonlinear tangents.
//! Its optional resisting-flow derivative T adds +T to the spectral pencil.
//! An AC solution does not imply that free transients decay to it.

use crate::{DescriptorPortHamiltonian, PortHamiltonian, Storage};
use crate::interconnection::accumulation::Sum;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;

/// Small-signal analysis at an admitted equilibrium of the original nonlinear storage.
pub mod equilibrium;
mod frequency;
mod parameter;
mod poles;
/// Bracketed driven-output peaks and half-power bandwidth from the original system.
pub mod resonance;
/// Continuous wave feedback, closed before factorization with physical power accounts.
pub mod scattering;
pub use frequency::HarmonicFrequencyResponse;
pub use parameter::{HarmonicParameterDirection, HarmonicParameterResponse};
/// Prepared physical port participation of complete descriptor modes.
pub use poles::response;
pub use poles::{BeynOptions, HarmonicSpectrum, HarmonicSpectrumError, ResonanceRegion};
pub use poles::{HarmonicPole, HarmonicPoleError, PoleOptions};

/// Hard ceiling of this explicitly dense frequency-domain path.
pub const MAX_HARMONIC_DIMENSION: usize = 1024;

// Fixed, bounded corrections reuse the same physical factorization. This is
// not a tolerance relaxation or a replacement of nearly cancelling states.
const MAX_RESPONSE_REFINEMENTS: usize = 4;

/// Caller-specified capacity and algebraic accuracy. No default drops states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicBudget {
    /// Full state count, including every algebraic multiplier.
    pub max_dimension: usize,
    /// Conservative peak internal numerical payload (not process RSS).
    pub max_dense_bytes: usize,
    /// Bound `n^2 * (max(n, 16) + ports + 8)` for preparation, factorization or
    /// one response, including at most four same-factor residual corrections.
    pub max_dense_work: u64,
    /// Original-equation residual ceiling: componentwise for forced AC;
    /// homogeneous poles use row coefficient norms times the full state norm.
    pub relative_residual_tolerance: f64,
}

/// A response is never published after failed admission, solve or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum HarmonicError {
    /// Invalid data, unsupported storage, or an exhausted explicit capacity.
    Invalid { what: &'static str },
    /// No unique response at the requested frequency (including a true pole).
    Singular { index: usize },
    /// Original equations did not meet the caller's algebraic tolerance.
    Residual { relative: f64, tolerance: f64 },
    /// A cooperative cancellation request was observed.
    Cancelled,
}

impl core::fmt::Display for HarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "harmonic response refused: {what}"),
            Self::Singular { index } => write!(f, "harmonic operator singular at pivot {index}"),
            Self::Residual {
                relative,
                tolerance,
            } => write!(
                f,
                "harmonic backward residual {relative:e} exceeds {tolerance:e}"
            ),
            Self::Cancelled => f.write_str("harmonic response cancelled"),
        }
    }
}
impl std::error::Error for HarmonicError {}
impl From<LuComplexError> for HarmonicError {
    fn from(error: LuComplexError) -> Self {
        match error {
            LuComplexError::Cancelled => Self::Cancelled,
            LuComplexError::Singular { index } => Self::Singular { index },
        }
    }
}

fn invalid(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HarmonicError> {
    if cancelled() {
        Err(HarmonicError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn resource_bounds(n: usize, m: usize) -> Result<(usize, u64), HarmonicError> {
    let square = n
        .checked_mul(n)
        .ok_or_else(|| invalid("dimension overflow"))?;
    let ports = n
        .checked_mul(m)
        .ok_or_else(|| invalid("port dimension overflow"))?;
    let bytes = square
        .checked_mul(88)
        .and_then(|v| ports.checked_mul(16).and_then(|p| v.checked_add(p)))
        .and_then(|v| n.checked_mul(256).and_then(|p| v.checked_add(p)))
        .and_then(|v| m.checked_mul(96).and_then(|p| v.checked_add(p)))
        .ok_or_else(|| invalid("payload bound overflow"))?;
    let work = (square as u64)
        .checked_mul(
            (n as u64)
                .max(4 * MAX_RESPONSE_REFINEMENTS as u64)
                .checked_add(m as u64)
                .and_then(|v| v.checked_add(8))
                .ok_or_else(|| invalid("work bound overflow"))?,
        )
        .ok_or_else(|| invalid("work bound overflow"))?;
    Ok((bytes, work))
}

/// Immutable complete linear image; assembly is performed once per system.
/// Owned coefficients prevent a later source edit from silently changing a fit.
pub struct HarmonicModel {
    n: usize,
    n_diff: usize,
    m: usize,
    q: Vec<f64>,
    generator: Vec<f64>,
    // Original J-R, retained for constitutive chain rules without inverting Q.
    structure: Vec<f64>,
    // Optional total D_x+D_e W; independent of the stored original R and Q.
    // Allocated/admitted only by the explicit nonlinear equilibrium path.
    resisting_tangent: Option<Vec<f64>>,
    r: Vec<f64>,
    g: Vec<f64>,
    budget: HarmonicBudget,
}

/// Reusable factor at one positive angular frequency. Multiple port excitations
/// share the same LU; no response or accepted simulation state is held here.
pub struct HarmonicFactor<'a> {
    model: &'a HarmonicModel,
    omega: f64,
    row_scales: Vec<f64>,
    column_scales: Vec<f64>,
    lu: LuComplex,
}

/// Cycle-mean and oscillating power at one original conjugate input/output pair.
///
/// For the owner's peak `exp(-i omega t)` phasors, the complex exchange is
/// `conj(u) y / 2`. Its real part transfers energy; its imaginary part retains
/// phase information that mean work alone loses. Quadrature is positive for
/// an effort-driven inertia/inductance and negative for compliance/capacitance.
/// Swapping input and output reverses ONLY that quadrature sign. It must not
/// be relabelled as electrical reactive power without the voltage/current chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPortPower {
    /// Mean power, positive into the modeled system [W]. May be negative at
    /// a source absorbing energy from other ports.
    pub active_w: f64,
    /// `Im(conj(u) y)/2`, with the input/output convention above [W].
    /// This is reversible exchange, never an additional heat source.
    pub quadrature_w: f64,
    /// `|u| |y|/2` [W], equal to the product of the two RMS amplitudes.
    /// Also the amplitude of the instantaneous power's twice-frequency ripple.
    /// Apparent powers of distinct ports do not cancel like signed mean work.
    pub apparent_w: f64,
}

impl HarmonicPortPower {
    /// Observe one scalar conjugate pair without solving or changing its state.
    /// Both amplitudes must be PEAK phasors in original power-conjugate units.
    /// Nonfinite amplitudes and unrepresentable physical power refuse.
    pub fn from_peak(input: C64, output: C64) -> Result<Self, HarmonicError> {
        if !finite(input) || !finite(output) {
            return Err(invalid("finite peak conjugate port amplitudes required"));
        }
        // Retain product remainders and apply the peak factor to each scalar
        // pair. Scaling an entire phasor can erase its tiny other quadrature.
        let mut active = Sum::default();
        let mut quadrature = Sum::default();
        active.half_product(input.re, output.re)
            .and_then(|()| active.half_product(input.im, output.im))
            .and_then(|()| quadrature.half_product(input.re, output.im))
            .and_then(|()| quadrature.half_product(-input.im, output.re))
            .ok_or_else(|| invalid("conjugate port power is not representable"))?;
        let exchange = C64::new(
            active.value().ok_or_else(|| invalid("conjugate port power is not representable"))?,
            quadrature.value().ok_or_else(|| invalid("conjugate port power is not representable"))?,
        );
        let apparent_w = exchange.abs();
        if !finite(exchange) || !apparent_w.is_finite() {
            return Err(invalid("conjugate port power is not representable"));
        }
        Ok(Self {
            active_w: exchange.re,
            quadrature_w: exchange.im,
            apparent_w,
        })
    }

    /// Signed mean/apparent power. Zero-amplitude ports have no defined power
    /// factor; `None` does not silently assign them unity or zero efficiency.
    /// For one sinusoid this is the phase cosine, not a conversion efficiency.
    #[must_use]
    pub fn power_factor(&self) -> Option<f64> {
        (self.apparent_w > 0.0).then(|| self.active_w / self.apparent_w)
    }
}

/// Peak complex amplitudes and cycle averages in the source's original units.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicResponse {
    /// Angular frequency [rad/s]; convention is exp(-i omega t).
    pub omega_rad_s: f64,
    /// Original differential coordinates followed by algebraic multipliers.
    pub state: Vec<C64>,
    /// `[Q_diff x_diff, x_alg]`, including constraint reactions.
    pub effort: Vec<C64>,
    /// Original power-conjugate outputs `G^T effort`, in original port order.
    pub output: Vec<C64>,
    /// Mean stored energy `Re(x* Q x)/4` [J]; multipliers store no energy.
    pub mean_stored_energy_j: f64,
    /// Per-port mean supplied power `Re(conj(u_p) y_p)/2` [W].
    pub port_power_w: Vec<f64>,
    /// Complete scalar power components in the same original port order.
    /// These use actual boundary inputs, including physical wave feedback.
    pub port_power_components: Vec<HarmonicPortPower>,
    /// Aggregate original port exchange, signed into the system [W]. Carries
    /// the original G*u contraction through the effort before final rounding;
    /// work smaller than individual output/power rounding remains observable.
    pub supplied_power_w: f64,
    /// Independently evaluated Re(e* R e + e* T x)/2 [W]. T is present only
    /// for explicit nonlinear equilibrium models. Then this is signed
    /// INCREMENTAL work, not the total physical heat at the biased source.
    pub dissipated_power_w: f64,
    /// Dissipation minus supply; never used to manufacture physical loss [W].
    pub power_residual_w: f64,
    /// Maximum componentwise original-equation backward error.
    pub relative_residual: f64,
}

impl HarmonicModel {
    /// Prepare an ordinary linear pHS, retaining every declared internal state.
    /// Nonlinear or undeclared storage refuses instead of being probed or fitted.
    pub fn from_ode(
        system: &PortHamiltonian,
        budget: HarmonicBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HarmonicError> {
        Self::prepare(
            system.n,
            system.n,
            system.m,
            &system.j,
            &system.r,
            &system.g,
            system.storage.as_ref(),
            budget,
            &mut cancelled,
        )
    }

    /// Prepare the original constrained system, including Kirchhoff equations,
    /// ideal transformers and all coupled inductive/capacitive storage.
    /// A singular descriptor mass is expected; it is never inverted.
    pub fn from_descriptor(
        system: &DescriptorPortHamiltonian,
        budget: HarmonicBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HarmonicError> {
        Self::prepare(
            system.n,
            system.n_diff,
            system.m,
            &system.j,
            &system.r,
            &system.g,
            system.storage.as_ref(),
            budget,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        n: usize,
        n_diff: usize,
        m: usize,
        j: &[f64],
        r: &[f64],
        g: &[f64],
        storage: &dyn Storage,
        budget: HarmonicBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, HarmonicError> {
        poll(cancelled)?;
        let square = n
            .checked_mul(n)
            .ok_or_else(|| invalid("dimension overflow"))?;
        let ports = n
            .checked_mul(m)
            .ok_or_else(|| invalid("port dimension overflow"))?;
        let (bytes, work) = resource_bounds(n, m)?;
        if n == 0
            || n_diff > n
            || n > MAX_HARMONIC_DIMENSION
            || n > budget.max_dimension
            || bytes > budget.max_dense_bytes
            || work > budget.max_dense_work
            || !budget.relative_residual_tolerance.is_finite()
            || budget.relative_residual_tolerance <= 0.0
            || budget.relative_residual_tolerance >= 1.0
        {
            return Err(invalid(
                "invalid accuracy or insufficient dimension, byte or work budget",
            ));
        }
        if j.len() != square
            || r.len() != square
            || g.len() != ports
            || j.iter().chain(r).chain(g).any(|v| !v.is_finite())
        {
            return Err(invalid("invalid or nonfinite original structure"));
        }
        let part = storage
            .quadratic_part()
            .ok_or_else(|| invalid("storage has no declared quadratic law"))?;
        if !part.exact || part.matrix.len() != square || part.matrix.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "storage must declare its complete finite quadratic law",
            ));
        }
        let q = part.matrix;
        // Neither algebraic storage nor a differential/algebraic energy cross
        // term belongs to the descriptor convention of this crate.
        for i in 0..n {
            poll(cancelled)?;
            for k in 0..n {
                if (i >= n_diff || k >= n_diff) && q[i * n + k] != 0.0 {
                    return Err(invalid("algebraic multipliers cannot carry storage"));
                }
            }
        }
        let mut structure = vec![0.0; square];
        let mut generator = vec![0.0; square];
        for i in 0..n {
            poll(cancelled)?;
            for k in 0..n {
                structure[i * n + k] = j[i * n + k] - r[i * n + k];
            }
            if structure[i * n..(i + 1) * n].iter().any(|v| !v.is_finite()) {
                return Err(invalid("original structure is not representable"));
            }
            for k in 0..n {
                generator[i * n + k] = if k < n_diff {
                    (0..n_diff)
                        .map(|l| structure[i * n + l] * q[l * n + k])
                        .sum()
                } else {
                    structure[i * n + k]
                };
            }
        }
        if generator.iter().any(|v| !v.is_finite()) {
            return Err(invalid("linear generator is not representable"));
        }
        poll(cancelled)?;
        Ok(Self {
            n,
            n_diff,
            m,
            q,
            generator,
            structure,
            resisting_tangent: None,
            r: r.to_vec(),
            g: g.to_vec(),
            budget,
        })
    }

    /// Factor at one positive frequency; independent row/column scaling handles
    /// mixed SI units without changing constraints, damping or state ownership.
    /// A genuine pole or a redundant unconstrained coordinate still refuses.
    pub fn factor(
        &self,
        omega_rad_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFactor<'_>, HarmonicError> {
        poll(&mut cancelled)?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("AC frequency must be positive and finite"));
        }
        let n = self.n;
        let mut matrix: Vec<_> = self.generator.iter().map(|v| C64::from_re(-v)).collect();
        for i in 0..self.n_diff {
            matrix[i * n + i].im = -omega_rad_s;
        }
        let mut row_scales = vec![0.0; n];
        for (i, row) in matrix.chunks_exact_mut(n).enumerate() {
            poll(&mut cancelled)?;
            let scale = row.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return Err(HarmonicError::Singular { index: i });
            }
            if !scale.is_finite() {
                return Err(invalid("operator row scale overflow"));
            }
            row_scales[i] = scale;
            for value in row {
                *value = C64::new(value.re / scale, value.im / scale);
            }
        }
        let mut column_scales = vec![0.0; n];
        for k in 0..n {
            poll(&mut cancelled)?;
            let scale = (0..n)
                .map(|i| matrix[i * n + k].abs())
                .fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return Err(HarmonicError::Singular { index: k });
            }
            column_scales[k] = scale;
            for i in 0..n {
                let z = &mut matrix[i * n + k];
                *z = C64::new(z.re / scale, z.im / scale);
            }
        }
        let lu = lu_complex_with_cancel(&matrix, n, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(HarmonicFactor {
            model: self,
            omega: omega_rad_s,
            row_scales,
            column_scales,
            lu,
        })
    }
}

impl HarmonicFactor<'_> {
    /// Evaluate a new excitation with the already factored original operator.
    /// Phasors are peak amplitudes, not RMS. All work and loss are cycle means.
    /// For explicit nonlinear local models these are perturbation quantities;
    /// negative incremental loss does not mean the original resistor creates heat.
    pub fn response(
        &self,
        input: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = self.model;
        let (n, m) = (model.n, model.m);
        if input.len() != m || input.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "one finite peak phasor per source port is required",
            ));
        }
        let mut rhs = vec![C64::ZERO; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            rhs[i] = input_value(input_row_sums(&model.g[i*m..(i+1)*m], |p| input[p])?, C64::ZERO)?;
        }
        let (state, relative) = self.solve_rhs(&rhs, &mut cancelled)?;
        let mut effort = state.clone();
        let mut energy = 0.0;
        for i in 0..model.n_diff {
            poll(&mut cancelled)?;
            effort[i] = C64::ZERO;
            for k in 0..model.n_diff {
                effort[i] = effort[i] + state[k].scale(model.q[i * n + k]);
            }
            energy += (state[i].conj() * effort[i]).re * 0.25;
        }
        let mut loss = 0.0;
        let mut output = vec![C64::ZERO; m];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut resistance = C64::ZERO;
            for k in 0..n {
                resistance = resistance + effort[k].scale(model.r[i * n + k]);
            }
            if model.resisting_tangent.is_some() {
                resistance = resistance + model.nonlinear_loss_row(i, &state);
            }
            loss += (effort[i].conj() * resistance).re * 0.5;
        }
        for (p, value) in output.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            *value = port_output(&model.g, m, p, &effort)?;
        }
        let mut port_power_components = Vec::with_capacity(m);
        let mut port_power = Vec::with_capacity(m);
        for (&u, &y) in input.iter().zip(&output) {
            poll(&mut cancelled)?;
            let power = HarmonicPortPower::from_peak(u, y)?;
            port_power.push(power.active_w);
            port_power_components.push(power);
        }
        let supply = physical_power(model, input, &effort, &mut cancelled)?;
        let balance = loss - supply;
        if effort.iter().chain(&output).any(|&v| !finite(v))
            || port_power.iter().any(|v| !v.is_finite())
            || [energy, loss, supply, balance]
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(invalid("physical effort, power or storage overflow"));
        }
        poll(&mut cancelled)?;
        Ok(HarmonicResponse {
            omega_rad_s: self.omega,
            state,
            effort,
            output,
            mean_stored_energy_j: energy,
            port_power_w: port_power,
            port_power_components,
            supplied_power_w: supply,
            dissipated_power_w: loss,
            power_residual_w: balance,
            relative_residual: relative,
        })
    }

    // Original-coordinate solve shared by the drive and its physical tangents.
    // The factor, equilibration and original-equation refinement are identical.
    fn solve_rhs(
        &self,
        rhs: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, f64), HarmonicError> {
        let model = self.model;
        let n = model.n;
        let mut state: Vec<_> = rhs
            .iter()
            .zip(&self.row_scales)
            .map(|(z, s)| C64::new(z.re / s, z.im / s))
            .collect();
        if state.iter().any(|&v| !finite(v)) {
            return Err(invalid("scaled port forcing overflow"));
        }
        self.lu.solve_with_cancel(&mut state, &mut *cancelled)?;
        for (z, s) in state.iter_mut().zip(&self.column_scales) {
            *z = C64::new(z.re / s, z.im / s);
        }
        if state.iter().any(|&v| !finite(v)) {
            return Err(invalid("response is not representable"));
        }
        // Componentwise backward error preserves small algebraic constraints
        // alongside much larger dimensional differential rows. At an exact
        // antiresonance some states vanish while coupled states remain finite;
        // one LU solve can leave roundoff large relative to these small rows.
        // Correct the ORIGINAL residual using the existing factor, preserving
        // the requested gate, physical coordinates and genuine singularities.
        let mut correction = vec![C64::ZERO; n];
        let mut relative: f64 = 0.0;
        for pass in 0..=MAX_RESPONSE_REFINEMENTS {
            relative = 0.0;
            for i in 0..n {
                poll(cancelled)?;
                let mut residual = -rhs[i];
                let mut denominator = rhs[i].abs();
                for k in 0..n {
                    let a = C64::new(
                        -model.generator[i * n + k],
                        if i == k && i < model.n_diff {
                            -self.omega
                        } else {
                            0.0
                        },
                    );
                    residual = residual + a * state[k];
                    denominator += a.abs() * state[k].abs();
                }
                if !finite(residual) || !denominator.is_finite() {
                    return Err(invalid("original-equation residual overflow"));
                }
                let row_error = if denominator == 0.0 {
                    0.0
                } else {
                    residual.abs() / denominator
                };
                relative = relative.max(row_error);
                correction[i] = C64::new(
                    -residual.re / self.row_scales[i],
                    -residual.im / self.row_scales[i],
                );
            }
            if relative <= model.budget.relative_residual_tolerance
                || pass == MAX_RESPONSE_REFINEMENTS
            {
                break;
            }
            if correction.iter().any(|&v| !finite(v)) {
                return Err(invalid("scaled residual correction overflow"));
            }
            self.lu
                .solve_with_cancel(&mut correction, &mut *cancelled)?;
            for i in 0..n {
                state[i] = state[i]
                    + C64::new(
                        correction[i].re / self.column_scales[i],
                        correction[i].im / self.column_scales[i],
                    );
            }
            if state.iter().any(|&v| !finite(v)) {
                return Err(invalid("refined response is not representable"));
            }
        }
        if relative > model.budget.relative_residual_tolerance {
            return Err(HarmonicError::Residual {
                relative,
                tolerance: model.budget.relative_residual_tolerance,
            });
        }
        Ok((state, relative))
    }
}

fn port_output(g: &[f64], ports: usize, port: usize, effort: &[C64]) -> Result<C64, HarmonicError> {
    let re = crate::port_output_component(g, ports, port, |i| effort[i].re);
    let im = crate::port_output_component(g, ports, port, |i| effort[i].im);
    Ok(C64::new(
        re.ok_or_else(|| invalid("physical port output is not representable"))?,
        im.ok_or_else(|| invalid("physical port output is not representable"))?,
    ))
}

// Keep BOTH terms of d(G^T e) in the same accumulator. Rounding G^T de and
// dG^T e separately can erase a weak tangent when those large terms cancel.
fn port_output_tangent(
    g: &[f64], ports: usize, port: usize, effort_derivative: &[C64],
    port_derivative: Option<&[f64]>, effort: &[C64],
) -> Result<C64, HarmonicError> {
    let contract = || -> Option<C64> {
        let mut re = Sum::default();
        let mut im = Sum::default();
        for (i, row) in g.chunks_exact(ports).enumerate() {
            let derivative = effort_derivative[i];
            if !finite(derivative) { return None; }
            if row[port] != 0.0 {
                re.product(row[port], derivative.re)?;
                im.product(row[port], derivative.im)?;
            }
            if let Some(dg) = port_derivative {
                let original = effort[i];
                if !finite(original) { return None; }
                let coefficient = dg[i * ports + port];
                if coefficient != 0.0 {
                    re.product(coefficient, original.re)?;
                    im.product(coefficient, original.im)?;
                }
            }
        }
        Some(C64::new(re.value()?, im.value()?))
    };
    contract().ok_or_else(|| invalid("physical port output tangent is not representable"))
}

// Retain each physical input row's two words until its final consumer. In
// particular, power and tangent rows must not consume an already-rounded G*u.
fn input_row_sums(row: &[f64], input: impl Fn(usize) -> C64) -> Result<[Sum; 2], HarmonicError> {
    Ok([
        crate::port_input_row(row, |p| input(p).re)
            .ok_or_else(|| invalid("physical input contraction is not representable"))?,
        crate::port_input_row(row, |p| input(p).im)
            .ok_or_else(|| invalid("physical input contraction is not representable"))?,
    ])
}

fn input_tangent_sums(
    row: &[f64], direction: Option<&[f64]>, input: &[C64], derivative: Option<&[C64]>,
    keep: impl Fn(usize) -> bool,
) -> Result<[Sum; 2], HarmonicError> {
    let mut sums = if let Some(du) = derivative {
        input_row_sums(row, |p| if keep(p) { du[p] } else { C64::ZERO })?
    } else {
        [Sum::default(); 2]
    };
    if let Some(dg) = direction {
        let terms = input_row_sums(dg, |p| if keep(p) { input[p] } else { C64::ZERO })?;
        for k in 0..2 {
            sums[k].add_scaled(terms[k], 1.0)
                .ok_or_else(|| invalid("physical input tangent is not representable"))?;
        }
    }
    Ok(sums)
}

// The existing generator action or i*E*x joins the unrounded input terms.
// Adding it after input_value would lose weak motion when these terms cancel.
fn input_value([mut re, mut im]: [Sum; 2], initial: C64) -> Result<C64, HarmonicError> {
    re.product(initial.re, 1.0)
        .and_then(|()| im.product(initial.im, 1.0))
        .ok_or_else(|| invalid("physical input row is not representable"))?;
    Ok(C64::new(
        re.value().ok_or_else(|| invalid("physical input row is not representable"))?,
        im.value().ok_or_else(|| invalid("physical input row is not representable"))?,
    ))
}

fn physical_power(
    model: &HarmonicModel, input: &[C64], effort: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, HarmonicError> {
    physical_power_sum(model, input, effort, |_| true, cancelled)?.value()
        .ok_or_else(|| invalid("physical supplied power is not representable"))
}

// A wave account selects only its ordinary ports, then appends its independent
// incident/outgoing fluxes to this same retained sum before final rounding.
fn physical_power_sum(
    model: &HarmonicModel, input: &[C64], effort: &[C64], keep: impl Fn(usize) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Sum, HarmonicError> {
    let mut power = Sum::default();
    for (i, e) in effort.iter().enumerate() {
        poll(cancelled)?;
        let [re, im] = input_row_sums(&model.g[i*model.m..(i+1)*model.m],
            |p| if keep(p) { input[p] } else { C64::ZERO })?;
        power.add_half_scaled(re, e.re)
            .and_then(|()| power.add_half_scaled(im, e.im))
            .ok_or_else(|| invalid("physical supplied power is not representable"))?;
    }
    Ok(power)
}

fn physical_power_tangent(
    model: &HarmonicModel, input: &[C64], input_derivative: Option<&[C64]>,
    effort: &[C64], effort_derivative: &[C64], port_derivative: Option<&[f64]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, HarmonicError> {
    physical_power_tangent_sum(model, input, input_derivative, effort,
        effort_derivative, port_derivative, |_| true, cancelled)?.value()
        .ok_or_else(|| invalid("physical supplied-power tangent is not representable"))
}

#[allow(clippy::too_many_arguments)] // retain original efforts, source directions and port selection
fn physical_power_tangent_sum(
    model: &HarmonicModel, input: &[C64], input_derivative: Option<&[C64]>,
    effort: &[C64], effort_derivative: &[C64], port_derivative: Option<&[f64]>,
    keep: impl Fn(usize) -> bool, cancelled: &mut impl FnMut() -> bool,
) -> Result<Sum, HarmonicError> {
    let mut power = Sum::default();
    for (i, (e, de)) in effort.iter().zip(effort_derivative).enumerate() {
        poll(cancelled)?;
        let range = i*model.m..(i+1)*model.m;
        let original = input_row_sums(&model.g[range.clone()],
            |p| if keep(p) { input[p] } else { C64::ZERO })?;
        let derivative = input_tangent_sums(&model.g[range.clone()],
            port_derivative.map(|dg| &dg[range]), input, input_derivative, &keep)?;
        for (k, (e, de)) in [(e.re, de.re), (e.im, de.im)].into_iter().enumerate() {
            power.add_half_scaled(original[k], de)
                .and_then(|()| power.add_half_scaled(derivative[k], e))
                .ok_or_else(|| invalid("physical supplied-power tangent is not representable"))?;
        }
    }
    Ok(power)
}

fn peak_power_tangent(
    input: C64, output: C64, input_derivative: Option<C64>, output_derivative: C64,
) -> Result<f64, HarmonicError> {
    let mut power = Sum::default();
    if let Some(du) = input_derivative {
        power.half_product(du.re, output.re)
            .and_then(|()| power.half_product(du.im, output.im))
            .ok_or_else(|| invalid("conjugate port-power tangent is not representable"))?;
    }
    power.half_product(input.re, output_derivative.re)
        .and_then(|()| power.half_product(input.im, output_derivative.im))
        .ok_or_else(|| invalid("conjugate port-power tangent is not representable"))?;
    power.value().ok_or_else(|| invalid("conjugate port-power tangent is not representable"))
}
