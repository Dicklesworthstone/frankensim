//! Exact transient derivatives of the ORIGINAL coupled dissipation laws.
//!
//! The physical state remains `[q, M*v, sqrt(k)*viscous_strain]`. At fixed
//! mass, stiffness, branch projections and force map, `x_theta' = A*x_theta
//! + A_theta*x + G*u_theta`. The upper-right block of
//! `exp(h * [[E, E_theta], [0, E]])`, `E=[[A,G],[0,0]]`, supplies the exact
//! held-input parameter derivative. The native block exponential is prepared
//! once; no mode differentiation, trajectory perturbation or new integrator
//! is involved. This read-only observation advances no physical clock.

use super::LinearMechanicalSystem;
use crate::interconnection::accumulation::Sum;
use crate::linear_exact::{
    DenseExponentialBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, UpperBlocks,
    upper_block_exponential_with_budget,
};
use crate::{PhsError, PreparedStepError};
use fs_math::det;

/// One authored physical direction per unit of a scalar parameter `theta`.
///
/// Mass, relaxed stiffness, branch stiffness/projections and the force map are
/// fixed. The direction can be signed and need not itself be passive; its BASE
/// `LinearMechanicalSystem` must already be admitted. This is a local derivative,
/// not permission to take a finite parameter step outside that passive domain.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalDampingDirection<'a> {
    /// Row-major n-by-n `dC/dtheta` in the original force/velocity coordinates.
    /// All off-diagonal and skew gyroscopic terms remain. This is dimensional
    /// C, NOT `M^-1 C`, a fitted modal ratio, or a derivative of its factor.
    pub damping: &'a [f64],
    /// `d(1/tau_j)/dtheta`, one entry per ORIGINAL Maxwell branch, in 1/s per
    /// unit parameter. Equal-time and zero-projection branches stay distinct.
    /// For a time-constant parameter use `-tau_theta/tau^2`; for a logarithmic
    /// time parameter use `-1/tau`. Stiffness and canonical branch units stay fixed.
    pub relaxation_rates: &'a [f64],
}

/// Immutable exact held-input tangent and pullback maps for one physical law.
///
/// Cache preparation binds the original source, parameter direction and step
/// duration. A changed source/direction/clock requires a new preparation. Warm
/// calls allocate nothing and perform compensated matrix-vector contractions
/// only. They neither retain nor publish physical state. Continue stepping the
/// source with its existing native owner, using the same held input and clock.
pub struct PreparedMechanicalDampingSensitivity {
    n: usize,
    m: usize,
    duration: f64,
    transition: Vec<f64>,
    direction: Vec<f64>,
}

impl LinearMechanicalSystem {
    /// Cold-prepare an exact local derivative of full C and Maxwell decay rates.
    ///
    /// Both copies of every state AND held force are included in the explicit
    /// dimension: `2*(state_dim+port_dim) <= 256`. Before allocating, this owner
    /// reserves `12*N^2*sizeof(f64)` bytes and `96*N^3` multiply-add terms, with
    /// N the doubled dimension, covering source assembly, the native block
    /// exponential and retained maps. Supplied limits are never increased.
    /// Cancellation is checked during assembly and on either side of the one
    /// bounded native exponential. This is not a latency or conditioning proof.
    pub fn prepare_damping_sensitivity(
        &self,
        direction: MechanicalDampingDirection<'_>,
        duration: f64,
        budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalDampingSensitivity, PreparedStepError> {
        poll(&mut cancelled)?;
        let n = self.system.state_dim();
        let m = self.system.port_dim();
        let d = n.checked_add(m).ok_or_else(capacity)?;
        let doubled = d.checked_mul(2).ok_or_else(capacity)?;
        let bytes = doubled
            .checked_mul(doubled)
            .and_then(|v| v.checked_mul(12 * core::mem::size_of::<f64>()));
        let work = (doubled as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(96));
        if doubled > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || doubled > budget.max_dimension
            || bytes.is_none_or(|v| v > budget.max_dense_bytes)
            || work.is_none_or(|v| v > budget.max_multiply_terms)
        {
            return Err(capacity().into());
        }
        let coordinates = self.coordinates;
        if direction.damping.len() != coordinates * coordinates
            || direction.relaxation_rates.len() != self.branches.len()
            || !finite(direction.damping)
            || !finite(direction.relaxation_rates)
            || !duration.is_finite()
            || duration <= 0.0
        {
            return Err(
                invalid("complete finite physical damping direction and positive clock").into(),
            );
        }
        let storage = self
            .system
            .storage()
            .quadratic_part()
            .filter(|part| part.exact)
            .ok_or_else(|| invalid("original exact mechanical storage"))?;
        let q = &storage.matrix;
        let (j, r, g) = self.system.structure();
        // Same diagonal storage equilibration as the ordinary exact-ZOH owner.
        // A free/gauge coordinate is retained with unit numerical scale.
        let mut scales = vec![1.0; d];
        for i in 0..n {
            if q[i * n + i] > 0.0 {
                scales[i] = det::sqrt(q[i * n + i]);
            }
        }
        let mut generator = vec![0.0; d * d];
        let mut derivative = vec![0.0; d * d];
        for row in 0..n {
            poll(&mut cancelled)?;
            for column in 0..n {
                let mut sum = Sum::default();
                for k in 0..n {
                    sum.product(j[row * n + k] - r[row * n + k], q[k * n + column])
                        .ok_or_else(arithmetic)?;
                }
                generator[row * d + column] = ratio(
                    sum.value().ok_or_else(arithmetic)?,
                    scales[row],
                    scales[column],
                )?;
            }
        }
        // Momentum derivative is -C_theta M^-1 p. No inverse-mass, factor or
        // eigenvector derivative belongs here: those physical sources are fixed.
        for row in 0..coordinates {
            poll(&mut cancelled)?;
            for column in 0..coordinates {
                let mut sum = Sum::default();
                for k in 0..coordinates {
                    sum.product(
                        -direction.damping[row * coordinates + k],
                        self.inverse_mass[k * coordinates + column],
                    )
                    .ok_or_else(arithmetic)?;
                }
                let i = coordinates + row;
                let k = coordinates + column;
                derivative[i * d + k] =
                    ratio(sum.value().ok_or_else(arithmetic)?, scales[i], scales[k])?;
            }
        }
        for (branch, &rate) in direction.relaxation_rates.iter().enumerate() {
            poll(&mut cancelled)?;
            let row = 2 * coordinates + branch;
            // z_j' = rho_j*(sqrt(k_j)*b_j^T q-z_j), in the public native
            // canonical history. It is not the unscaled stress/strain state.
            for column in 0..coordinates {
                // Q already retains the original admitted canonical coupling;
                // do not rebuild it from a second rounded square-root product.
                derivative[row * d + column] = ratio(
                    checked_product(-rate, q[row * n + column])?,
                    scales[row],
                    scales[column],
                )?;
            }
            derivative[row * d + row] = -rate;
        }
        // Constant held-force coordinates include their effect on the tangent
        // throughout the interval, not merely a terminal-force correction.
        for port in 0..m {
            poll(&mut cancelled)?;
            let mut scale = 0.0_f64;
            for row in 0..n {
                scale = scale.max(checked_product(scales[row], g[row * m + port])?.abs());
            }
            scales[n + port] = if scale == 0.0 { 1.0 } else { scale };
            for row in 0..n {
                generator[row * d + n + port] =
                    ratio(g[row * m + port], scales[row], scales[n + port])?;
            }
        }
        // Parameter units must not change the exponential's scaling burden.
        // Normalize only this numerical direction, then restore its exact units
        // in the cached original-coordinate derivative. No source is retuned.
        let direction_scale = derivative.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        for value in &mut generator {
            *value = checked_product(*value, duration)?;
        }
        if direction_scale != 0.0 {
            for value in &mut derivative {
                *value = checked_product(ratio(*value, 1.0, direction_scale)?, duration)?;
            }
        }
        let operator = UpperBlocks {
            first: generator.clone(),
            cross: derivative,
            last: generator,
        };
        poll(&mut cancelled)?;
        let flow = upper_block_exponential_with_budget(&operator, d, budget)?;
        poll(&mut cancelled)?;
        let mut transition = vec![0.0; n * d];
        let mut parameter = vec![0.0; n * d];
        for row in 0..n {
            poll(&mut cancelled)?;
            for column in 0..d {
                transition[row * d + column] =
                    ratio(flow.first[row * d + column], scales[column], scales[row])?;
                parameter[row * d + column] = checked_product(
                    ratio(flow.cross[row * d + column], scales[column], scales[row])?,
                    direction_scale,
                )?;
            }
        }
        Ok(PreparedMechanicalDampingSensitivity {
            n,
            m,
            duration,
            transition,
            direction: parameter,
        })
    }
}

impl PreparedMechanicalDampingSensitivity {
    pub const fn state_dim(&self) -> usize {
        self.n
    }
    pub const fn port_dim(&self) -> usize {
        self.m
    }
    pub const fn duration(&self) -> f64 {
        self.duration
    }

    /// Propagate `dx/dtheta`, including authored initial-state and held-force
    /// derivatives. All state arrays use the ORIGINAL canonical state order.
    /// For fixed initial conditions/forces their derivative arrays are zero.
    /// A physical velocity derivative is the original `M^-1` observation of
    /// this momentum derivative; `LinearMechanicalSystem::velocities_into`
    /// already supplies it because M is fixed. Outputs are candidate scratch
    /// and can change on arithmetic refusal; accepted physical state is untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_into(
        &self,
        initial: &[f64],
        initial_derivative: &[f64],
        input: &[f64],
        input_derivative: &[f64],
        derivative: &mut [f64],
    ) -> Result<(), PhsError> {
        self.admit(initial, input)?;
        self.admit(initial_derivative, input_derivative)?;
        if derivative.len() != self.n {
            return Err(invalid("mechanical tangent output extent"));
        }
        for (row, value) in derivative.iter_mut().enumerate() {
            let mut sum = self.parameter_row(row, initial, input)?;
            let offset = row * (self.n + self.m);
            for (column, &entry) in initial_derivative
                .iter()
                .chain(input_derivative)
                .enumerate()
            {
                sum.product(self.transition[offset + column], entry)
                    .ok_or_else(arithmetic)?;
            }
            *value = sum.value().ok_or_else(arithmetic)?;
        }
        Ok(())
    }

    /// Reverse one native exact interval. Return the local scalar `dJ/dtheta`
    /// at fixed initial state/input, and write their endpoint pullbacks.
    /// Sum the scalar over a reverse trajectory sweep; propagate the returned
    /// state cotangent to the preceding interval. Authored dependence of initial
    /// conditions or forces is an additional contraction with these cotangents.
    /// The parameter contribution is contracted before rounding each tangent
    /// row, retaining weak physical contributions beside cancelling components.
    /// Outputs are candidate scratch and may change on arithmetic refusal.
    pub fn pullback_into(
        &self,
        initial: &[f64],
        input: &[f64],
        endpoint_cotangent: &[f64],
        initial_cotangent: &mut [f64],
        input_cotangent: &mut [f64],
    ) -> Result<f64, PhsError> {
        self.admit(initial, input)?;
        if endpoint_cotangent.len() != self.n
            || !finite(endpoint_cotangent)
            || initial_cotangent.len() != self.n
            || input_cotangent.len() != self.m
        {
            return Err(invalid("complete finite mechanical pullback layouts"));
        }
        let mut parameter = Sum::default();
        for (row, &cotangent) in endpoint_cotangent.iter().enumerate() {
            parameter
                .add_scaled(self.parameter_row(row, initial, input)?, cotangent)
                .ok_or_else(arithmetic)?;
        }
        for (column, value) in initial_cotangent
            .iter_mut()
            .chain(input_cotangent)
            .enumerate()
        {
            let mut sum = Sum::default();
            for (row, &cotangent) in endpoint_cotangent.iter().enumerate() {
                sum.product(self.transition[row * (self.n + self.m) + column], cotangent)
                    .ok_or_else(arithmetic)?;
            }
            *value = sum.value().ok_or_else(arithmetic)?;
        }
        parameter.value().ok_or_else(arithmetic)
    }

    fn admit(&self, state: &[f64], input: &[f64]) -> Result<(), PhsError> {
        if state.len() != self.n || input.len() != self.m || !finite(state) || !finite(input) {
            return Err(invalid(
                "complete finite original mechanical state and held inputs",
            ));
        }
        Ok(())
    }

    fn parameter_row(&self, row: usize, state: &[f64], input: &[f64]) -> Result<Sum, PhsError> {
        let offset = row * (self.n + self.m);
        let mut sum = Sum::default();
        for (column, &value) in state.iter().chain(input).enumerate() {
            sum.product(self.direction[offset + column], value)
                .ok_or_else(arithmetic)?;
        }
        Ok(sum)
    }
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}
fn invalid(what: &'static str) -> PhsError {
    PhsError::Unadmitted { what }
}
fn arithmetic() -> PhsError {
    invalid("representable original mechanical damping derivative")
}
fn capacity() -> PhsError {
    PhsError::LinearSplitBudget {
        what: "complete doubled mechanical tangent dimension/bytes/work",
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}
fn checked_product(a: f64, b: f64) -> Result<f64, PhsError> {
    let value = a * b;
    if !value.is_finite() || (value == 0.0 && a != 0.0 && b != 0.0) {
        Err(arithmetic())
    } else {
        Ok(value)
    }
}

// Try both equivalent orderings before refusing an intermediate range loss.
// A nonzero original coefficient is never silently replaced with zero.
fn ratio(value: f64, numerator: f64, denominator: f64) -> Result<f64, PhsError> {
    if value == 0.0 {
        return Ok(value);
    }
    for candidate in [
        (value * numerator) / denominator,
        (value / denominator) * numerator,
    ] {
        if candidate.is_finite() && candidate != 0.0 {
            return Ok(candidate);
        }
    }
    Err(arithmetic())
}
