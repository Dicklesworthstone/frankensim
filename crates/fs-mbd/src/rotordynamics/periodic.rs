//! Period maps of declared piecewise-constant linear rotor-bearing cycles.
//!
//! Every segment retains the original equation
//! M q'' + (C + Omega G) q' + K q = 0. Signed spin, full mass and skew G
//! stay fixed; full C/K may change, including cross-coupled or active entries.
//! Coordinates and velocities are continuous at each coefficient switch.
//! Changing stiffness can exchange energy with its modulator; no change is
//! silently assigned heat.
//!
//! Native bounded matrix exponentials advance each constant segment. Their
//! chronological product E_last ... E_first is the complete monodromy map.
//! Its numerical multipliers can expose parametric resonance even when every
//! frozen operating point has stable poles. No modal truncation, artificial
//! damping or eigenvalue clipping is applied.
//!
//! This is the declared segmented model. Approximating smooth modulation by
//! such segments additionally requires coefficient/time refinement. Neither
//! eigenvalue error bounds nor a stability certificate are inferred. In
//! particular, rounding around unit spectral radius is not classified.

use super::{
    MAX_ROTOR_DOFS, ROTOR_QR_SWEEPS_PER_WINDOW, RotorBearingModel, RotorError, finite, poll,
};
use fs_la::eigen_complex::{EigError, eig_with_cancel};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::linear_exact::{DenseExponentialBudget, exponential_with_budget};

/// Hard bound on coefficient switches in a single declared physical cycle.
pub const MAX_ROTOR_PERIODIC_SEGMENTS: usize = 256;

/// One positive-duration constant-coefficient portion of a periodic rotor cycle.
///
/// Models must share the same coordinate order, units, frame, full mass and G.
/// Their C/K are actual coefficients for this portion of the cycle. The final
/// segment is followed by the first, without an impulse or state reset.
#[derive(Clone, Copy, Debug)]
pub struct RotorPeriodicSegment<'a> {
    /// Original full physical matrices; borrowed unchanged.
    pub model: &'a RotorBearingModel,
    /// Physical duration [s], finite and strictly positive.
    pub duration_seconds: f64,
}

/// Explicit admission for the complete period map and its numerical spectrum.
///
/// The existing 32-coordinate rotor bound still applies. Dense bytes include
/// all owner workspace and returned arrays, excluding the borrowed models.
/// Work conservatively reserves 80 dense products per segment exponential,
/// and additionally mass solves, chronological products and bounded complex QR.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RotorFloquetBudget {
    /// Supplied segment count ceiling, in 1..=MAX_ROTOR_PERIODIC_SEGMENTS.
    pub max_segments: usize,
    /// Native QR sweeps per active window, in 1..=ROTOR_QR_SWEEPS_PER_WINDOW.
    pub max_qr_sweeps_per_window: usize,
    /// Aggregate live dense payload ceiling in bytes.
    pub max_dense_bytes: usize,
    /// Conservative scalar multiply/add term ceiling for the whole query.
    pub max_multiply_terms: u64,
}

impl Default for RotorFloquetBudget {
    fn default() -> Self {
        Self {
            max_segments: 64,
            max_qr_sweeps_per_window: ROTOR_QR_SWEEPS_PER_WINDOW,
            max_dense_bytes: 1024 * 1024,
            max_multiply_terms: 4_000_000_000,
        }
    }
}

impl RotorFloquetBudget {
    fn admit(self, n: usize, segments: usize) -> Result<DenseExponentialBudget, RotorError> {
        if self.max_segments == 0
            || self.max_segments > MAX_ROTOR_PERIODIC_SEGMENTS
            || segments > self.max_segments
            || self.max_qr_sweeps_per_window == 0
            || self.max_qr_sweeps_per_window > ROTOR_QR_SWEEPS_PER_WINDOW
        {
            return Err(RotorError::PeriodicBudget("segment or QR sweep admission"));
        }
        // n <= 32 and segments <= 256 were checked before these products.
        // Sixteen real d-by-d matrices plus vector storage cover exponential,
        // Cholesky, complex QR, current/candidate maps and returned payloads.
        let d = 2 * n;
        let bytes = (16 * d * d + 32 * d) * core::mem::size_of::<f64>();
        let cube = (d as u64).pow(3);
        let work = cube * (84 * segments as u64 + 64 * (self.max_qr_sweeps_per_window as u64 + 1));
        if bytes > self.max_dense_bytes || work > self.max_multiply_terms {
            return Err(RotorError::PeriodicBudget(
                "dense byte or whole-cycle work admission",
            ));
        }
        Ok(DenseExponentialBudget {
            max_dimension: d,
            max_dense_bytes: self.max_dense_bytes,
            max_multiply_terms: self.max_multiply_terms,
        })
    }
}

/// Complete numerical Floquet analysis of the original segmented model.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorFloquetResponse {
    /// Signed constant shaft speed used in every C + Omega G [rad/s].
    pub spin_rad_s: f64,
    /// Sum of the declared physical segment durations [s].
    pub period_seconds: f64,
    /// Number of original generalized displacement coordinates.
    pub dofs: usize,
    /// Original-coordinate row-major 2n-by-2n map, in [q, q_dot] state order.
    /// Acting on an initial state gives its homogeneous one-period endpoint.
    pub monodromy: Vec<f64>,
    /// All 2n multipliers, in native canonical (real, imaginary) order.
    /// These are numerical estimates without eigenvectors or error enclosures.
    pub multipliers: Vec<C64>,
    /// Maximum computed multiplier magnitude. Unit-radius rounding is retained.
    pub spectral_radius: f64,
    /// ln(spectral_radius) / period [1/s], without a stability classification.
    pub growth_rate_per_s: f64,
    /// Largest independently reconstructed, mass-scaled original column
    /// residual in forming the generator. Not an eigenvalue or motion bound.
    pub maximum_mass_residual: f64,
}

/// Analyze a homogeneous cycle of original constant-coefficient rotor segments.
///
/// M and G must agree coefficient by coefficient across all segments. Spin
/// remains constant, so no spin-acceleration, changing-inertia or impulse law
/// is inferred. All C/K entries survive, including nonconservative ones.
///
/// An invertible mass/time coordinate scaling conditions the full state.
/// Native Cholesky solves form the generator; original separate C and Omega G
/// enter the independent column-residual checks. The mass residual tolerance
/// must be finite and in (0,1). No inverse, mode reduction or damping floor is
/// substituted for the full physical matrices.
///
/// Cancellation is polled during assembly, products and native QR, and before
/// and after each shared exponential. One exponential is an indivisible call
/// bounded to 64 states, 32 Taylor products and 48 squarings; this is not an
/// operating-system latency guarantee. Any refusal publishes no partial map.
pub fn periodic_floquet(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RotorFloquetResponse, RotorError> {
    poll(&mut cancelled)?;
    if segments.is_empty()
        || segments.len() > MAX_ROTOR_PERIODIC_SEGMENTS
        || !spin_rad_s.is_finite()
        || !relative_mass_residual_tolerance.is_finite()
        || relative_mass_residual_tolerance <= 0.0
        || relative_mass_residual_tolerance >= 1.0
    {
        return Err(RotorError::Invalid(
            "finite spin, nonempty bounded cycle and mass tolerance in (0,1) required",
        ));
    }
    let reference = segments[0].model;
    let n = reference.dofs();
    let dense_budget = budget.admit(n, segments.len())?;
    let mut period = Dd::ZERO;
    for segment in segments {
        poll(&mut cancelled)?;
        if segment.model.dofs() != n
            || segment.model.mass() != reference.mass()
            || segment.model.gyro_per_rad_s() != reference.gyro_per_rad_s()
            || !segment.duration_seconds.is_finite()
            || segment.duration_seconds <= 0.0
        {
            return Err(RotorError::Invalid(
                "periodic segments require identical M/G and positive finite durations",
            ));
        }
        period = period + Dd::from_f64(segment.duration_seconds);
        finite(period.to_f64())?;
    }
    let period_seconds = finite(period.to_f64())?;
    let mass = MassSolve::new(reference, &mut cancelled)?;
    // z = [S q, T S q_dot], S_ii=sqrt(M_ii), T=physical period.
    // This retains free coordinates and makes the time scale common to all
    // segments; no per-segment basis changes contaminate the matrix product.
    let d = 2 * n;
    let mut state_scale = mass.scale.clone();
    for &scale in &mass.scale {
        let velocity_scale = finite(period_seconds * scale)?;
        if velocity_scale == 0.0 {
            return Err(RotorError::Invalid("periodic coordinate scale underflow"));
        }
        state_scale.push(velocity_scale);
    }
    let mut monodromy = vec![0.0; d * d];
    for i in 0..d {
        monodromy[i * d + i] = 1.0;
    }
    let mut maximum_mass_residual = 0.0_f64;
    for segment in segments {
        poll(&mut cancelled)?;
        let (generator, residual) = mass.segment_generator(
            segment,
            spin_rad_s,
            period_seconds,
            relative_mass_residual_tolerance,
            &mut cancelled,
        )?;
        maximum_mass_residual = maximum_mass_residual.max(residual);
        poll(&mut cancelled)?;
        let flow = exponential_with_budget(&generator, d, dense_budget)
            .map_err(RotorError::PeriodicExponential)?;
        poll(&mut cancelled)?;
        // Left multiplication is physical chronological order.
        monodromy = multiply(&flow, &monodromy, d, &mut cancelled)?;
    }
    let complex: Vec<_> = monodromy
        .iter()
        .map(|&value| C64::new(value, 0.0))
        .collect();
    let multipliers = eig_with_cancel(&complex, d, budget.max_qr_sweeps_per_window, &mut cancelled)
        .map_err(|error| match error {
            EigError::Cancelled => RotorError::Cancelled,
            other => RotorError::PeriodicSpectrum(other),
        })?;
    let mut spectral_radius = 0.0_f64;
    for value in &multipliers {
        finite(value.re)?;
        finite(value.im)?;
        let magnitude = finite(det::hypot(value.re, value.im))?;
        if magnitude == 0.0 {
            return Err(RotorError::Invalid(
                "periodic multiplier is unresolved or underflowed",
            ));
        }
        spectral_radius = spectral_radius.max(magnitude);
    }
    // A finite-time fundamental map is invertible. No vanished direction
    // can justify a zero multiplier or an infinite decay claim.
    if spectral_radius == 0.0 {
        return Err(RotorError::Invalid("periodic multiplier radius underflow"));
    }
    let growth_rate_per_s = finite(det::ln(spectral_radius) / period_seconds)?;
    for i in 0..d {
        poll(&mut cancelled)?;
        for j in 0..d {
            monodromy[i * d + j] =
                restore_entry(monodromy[i * d + j], state_scale[i], state_scale[j])?;
        }
    }
    poll(&mut cancelled)?;
    Ok(RotorFloquetResponse {
        spin_rad_s,
        period_seconds,
        dofs: n,
        monodromy,
        multipliers,
        spectral_radius,
        growth_rate_per_s,
        maximum_mass_residual,
    })
}

struct MassSolve<'a> {
    model: &'a RotorBearingModel,
    scale: Vec<f64>,
    factor: Cholesky,
}

impl<'a> MassSolve<'a> {
    fn new(
        model: &'a RotorBearingModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RotorError> {
        let n = model.dofs();
        let scale: Vec<_> = (0..n).map(|i| det::sqrt(model.mass()[i * n + i])).collect();
        let mut scaled = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                scaled[i * n + j] = finite((model.mass()[i * n + j] / scale[i]) / scale[j])?;
            }
        }
        let factor = cholesky(&scaled, n).map_err(RotorError::TransientLinear)?;
        poll(cancelled)?;
        Ok(Self {
            model,
            scale,
            factor,
        })
    }

    fn segment_generator(
        &self,
        segment: &RotorPeriodicSegment<'_>,
        spin: f64,
        period: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64), RotorError> {
        let n = self.model.dofs();
        let d = 2 * n;
        let h = segment.duration_seconds;
        let fraction = finite(h / period)?;
        if fraction == 0.0 {
            return Err(RotorError::Invalid("periodic segment duration underflow"));
        }
        let mut generator = vec![0.0; d * d];
        let mut maximum = 0.0_f64;
        for j in 0..n {
            poll(cancelled)?;
            generator[j * d + n + j] = fraction;
            for stiffness in [true, false] {
                let (column, residual) =
                    self.solve_column(segment.model, spin, stiffness, j, tolerance, cancelled)?;
                maximum = maximum.max(residual);
                for i in 0..n {
                    if stiffness {
                        generator[(n + i) * d + j] = product_three(-column[i], h, period)?;
                    } else {
                        generator[(n + i) * d + n + j] = product_three(-column[i], h, 1.0)?;
                    }
                }
            }
        }
        Ok((generator, maximum))
    }

    #[allow(clippy::too_many_arguments)] // explicit source coefficient and original column audit
    fn solve_column(
        &self,
        model: &RotorBearingModel,
        spin: f64,
        stiffness: bool,
        column: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; MAX_ROTOR_DOFS], f64), RotorError> {
        let n = self.model.dofs();
        let mut solution = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in solution.iter_mut().enumerate().take(n) {
            let index = i * n + column;
            let coefficient = if stiffness {
                model.stiffness()[index]
            } else {
                finite(spin.mul_add(model.gyro_per_rad_s()[index], model.damping()[index]))?
            };
            *value = finite((coefficient / self.scale[i]) / self.scale[column])?;
            if coefficient != 0.0 && *value == 0.0 {
                return Err(RotorError::Invalid(
                    "periodic coefficient scaling underflow",
                ));
            }
        }
        poll(cancelled)?;
        self.factor.solve(&mut solution[..n]);
        poll(cancelled)?;
        let mut residual_norm = 0.0_f64;
        let mut action_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let row_scale = Dd::from_f64(self.scale[i]);
            let mut action = Dd::ZERO;
            for (j, &value) in solution.iter().enumerate().take(n) {
                finite(value)?;
                // Reconstruct directly from the original mass, independently
                // of the rounded Cholesky matrix and its triangular factors.
                let coefficient = Dd::from_f64(self.model.mass()[i * n + j])
                    / row_scale
                    / Dd::from_f64(self.scale[j]);
                action = action + coefficient * Dd::from_f64(value);
            }
            let index = i * n + column;
            let original_rhs = if stiffness {
                Dd::from_f64(model.stiffness()[index])
            } else {
                Dd::from_f64(model.damping()[index])
                    + Dd::from_f64(spin) * Dd::from_f64(model.gyro_per_rad_s()[index])
            };
            let rhs = original_rhs / row_scale / Dd::from_f64(self.scale[column]);
            residual_norm = residual_norm.max(finite((action - rhs).to_f64().abs())?);
            action_norm = action_norm.max(finite(action.to_f64().abs())?);
            rhs_norm = rhs_norm.max(finite(rhs.to_f64().abs())?);
        }
        let scale = action_norm.max(rhs_norm);
        let relative = if scale == 0.0 {
            residual_norm
        } else {
            finite((residual_norm / scale) / (action_norm / scale + rhs_norm / scale))?
        };
        if relative > tolerance {
            return Err(RotorError::Residual {
                relative,
                limit: tolerance,
            });
        }
        Ok((solution, relative))
    }
}

fn multiply(
    left: &[f64],
    right: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, RotorError> {
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            let mut value = 0.0_f64;
            for k in 0..n {
                value = left[i * n + k].mul_add(right[k * n + j], value);
            }
            result[i * n + j] = finite(value)?;
        }
    }
    Ok(result)
}

// Avoid avoidable intermediate overflow/underflow in h*T*coefficient.
// This changes only arithmetic association, never the original coefficient.
fn product_three(a: f64, b: f64, c: f64) -> Result<f64, RotorError> {
    if a == 0.0 {
        return Ok(a);
    }
    for (x, y, z) in [(a, b, c), (a, c, b), (b, c, a)] {
        let pair = x * y;
        if pair.is_finite() && pair != 0.0 {
            let value = pair * z;
            if value.is_finite() && value != 0.0 {
                return Ok(value);
            }
        }
    }
    Err(RotorError::Invalid(
        "periodic generator scaling is unrepresentable",
    ))
}

// S^-1 P S for the original [q,v] state; zero entries stay exactly zero.
fn restore_entry(value: f64, row_scale: f64, column_scale: f64) -> Result<f64, RotorError> {
    if value == 0.0 {
        return Ok(value);
    }
    let ratio = column_scale / row_scale;
    if ratio.is_finite() && ratio > 0.0 {
        let restored = value * ratio;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    let product = value * column_scale;
    if product.is_finite() && product != 0.0 {
        let restored = product / row_scale;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    let quotient = value / row_scale;
    if quotient.is_finite() && quotient != 0.0 {
        let restored = quotient * column_scale;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    Err(RotorError::Invalid(
        "original periodic state map is unrepresentable",
    ))
}
