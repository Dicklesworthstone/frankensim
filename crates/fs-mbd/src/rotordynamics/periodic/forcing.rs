//! Held physical loads and periodic forced orbits through the shared cycle owner.
use super::{
    PreparedRotorCycle, RotorFloquetBudget, RotorFloquetResponse, RotorPeriodicSegment, finite,
    poll, prepare_cycle, restore_entry,
};
use crate::rotordynamics::{RotorError, transient::RotorTransientState};
use fs_la::factor::lu;
use fs_math::dd::Dd;

/// Actual original-coordinate states at every coefficient/force boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorPeriodicTrajectory {
    /// Elapsed times [s], including zero and the complete physical period.
    pub boundary_times_seconds: Vec<f64>,
    /// Displacement/angle and physical velocity at the corresponding times.
    /// Force switches do not reset either state.
    pub boundary_states: Vec<RotorTransientState>,
    /// Maximum row-relative defect of the original affine increment equation,
    /// reconstructed in double-double arithmetic on the actual returned states.
    /// This is not a continuum integration or state-error bound.
    pub relative_affine_residual: f64,
}

/// A unique numerical periodic particular orbit of the declared forced cycle.
///
/// This is not a claim of attraction. The accompanying prepared Floquet report
/// retains any growing perturbations about this particular orbit.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorPeriodicOrbit {
    /// One physical period, including the independently propagated endpoint.
    pub trajectory: RotorPeriodicTrajectory,
    /// Original affine fixed-point residual at the actual original initial state.
    /// Its row scale uses (P-I)x and the forcing offset, not large cancelling x/Px.
    pub relative_fixed_point_residual: f64,
    /// Actual endpoint/start mismatch in the common mass/time scaled maximum
    /// norm, relative to the full-state peak over returned boundaries. An
    /// exactly zero trajectory has an explicit zero-scale case.
    pub relative_closure_residual: f64,
    /// Estimated ||(I-P)^-1||_1 (1+||P||_1) in the COMMON mass/time scaled
    /// state coordinates. The native estimate is not a norm bound or a
    /// coordinate-invariant physical condition number.
    pub resolvent_gain_estimate: f64,
}

/// Prepare exact affine maps for a declared piecewise-constant physical load.
///
/// Each force row is the original generalized force/torque held over the whole
/// corresponding segment, in the same coordinate order and conjugate units as
/// the model. Coefficients and loads repeat after the declared total period.
/// This does not interpret midpoint samples as a continuous sinusoid, and a
/// force jump is not an impulse. Smooth forcing needs independent refinement.
/// Coefficient switches can exchange parametric work. This API returns motion
/// and numerical residuals, not a full actuator-work or energy receipt.
///
/// The SAME full mass/time scaling, original-equation checks, segment
/// exponentials and Floquet owner serve the homogeneous and forced APIs.
/// A normalized constant source coordinate raises the largest indivisible
/// exponential from 64 to 65 coordinates; it is not a physical state or mode.
/// Explicit budgets include retained affine maps and one trajectory/solve.
/// Failure or cancellation publishes no partially prepared plan.
pub fn prepare_forced_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    held_forces: &[&[f64]],
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorCycle, RotorError> {
    prepare_cycle(
        spin_rad_s,
        segments,
        Some(held_forces),
        relative_mass_residual_tolerance,
        budget,
        cancelled,
    )
}

impl PreparedRotorCycle {
    /// Full homogeneous spectrum and original-coordinate period matrix.
    pub fn floquet_report(&self) -> &RotorFloquetResponse {
        &self.report
    }

    /// Original [q, q_dot] affine offset b in x_next = P x + b.
    pub fn forcing_offset(&self) -> &[f64] {
        &self.physical_offset
    }

    /// Fixed boundary times [s], including zero and the total period.
    pub fn boundary_times_seconds(&self) -> &[f64] {
        &self.boundary_times
    }

    /// Advance one cycle from an actual original physical state.
    ///
    /// Cached maps cost O(segments * (2n)^2); no new factorization,
    /// exponential or steady-state assumption is introduced. Each subsequent
    /// segment starts from its actually returned original q/v boundary.
    /// The final check uses the ORIGINAL affine increment equation row by row:
    /// (x_end-x_start) - (P-I)x_start - b = 0. A small physical load is not
    /// hidden by the magnitude of a large unchanged state.
    pub fn advance(
        &self,
        initial: &RotorTransientState,
        relative_residual_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorPeriodicTrajectory, RotorError> {
        poll(&mut cancelled)?;
        admit_tolerance(relative_residual_tolerance)?;
        let initial_values = state_values(initial, self.report.dofs)?;
        let d = self.state_scale.len();
        let mut balanced = encode(&initial_values, &self.state_scale, &mut cancelled)?;
        let mut states = Vec::with_capacity(self.segment_maps.len() + 1);
        states.push(initial.clone());
        for segment in &self.segment_maps {
            poll(&mut cancelled)?;
            let candidate = apply_affine(
                &segment.transition,
                &segment.offset,
                &balanced,
                d,
                &mut cancelled,
            )?;
            let state = decode(&candidate, &self.state_scale, &mut cancelled)?;
            // Round-trip the public physical endpoint rather than silently
            // continuing a more precise, hidden balanced state.
            balanced = encode(
                &state_values(&state, self.report.dofs)?,
                &self.state_scale,
                &mut cancelled,
            )?;
            states.push(state);
        }
        let final_values = state_values(&states[states.len() - 1], self.report.dofs)?;
        let relative_affine_residual = affine_residual(
            &self.report.monodromy,
            &self.physical_offset,
            &initial_values,
            &final_values,
            relative_residual_tolerance,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(RotorPeriodicTrajectory {
            boundary_times_seconds: self.boundary_times.clone(),
            boundary_states: states,
            relative_affine_residual,
        })
    }

    /// Solve (I-P)x = b for a UNIQUE periodic particular orbit.
    ///
    /// Native partial-pivot LU and its bounded condition estimator operate in
    /// the common mass/time coordinates. The caller's finite gain limit (>=1)
    /// applies to ||(I-P)^-1||_1 (1+||P||_1), so a nearly identity period map
    /// cannot appear well resolved merely because cond(I-P) is small.
    /// This is a numerical admission estimate, not a certified error bound.
    ///
    /// Singular maps refuse even for compatible forcing with a nonunique
    /// periodic family. An unstable but nonsingular cycle may have a unique
    /// particular orbit; its homogeneous growth remains in floquet_report().
    /// The actual original initial state must pass a rowwise force-scale DD
    /// residual, and the separately propagated physical endpoint must close.
    pub fn periodic_orbit(
        &self,
        relative_residual_tolerance: f64,
        maximum_resolvent_gain: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorPeriodicOrbit, RotorError> {
        poll(&mut cancelled)?;
        admit_tolerance(relative_residual_tolerance)?;
        if !maximum_resolvent_gain.is_finite() || maximum_resolvent_gain < 1.0 {
            return Err(RotorError::Invalid(
                "finite periodic resolvent limit >=1 required",
            ));
        }
        let d = self.state_scale.len();
        let mut equation = vec![0.0; d * d];
        for i in 0..d {
            poll(&mut cancelled)?;
            for j in 0..d {
                let identity = if i == j { 1.0 } else { 0.0 };
                equation[i * d + j] = finite(identity - self.monodromy[i * d + j])?;
            }
        }
        poll(&mut cancelled)?;
        let factor = lu(&equation, d).map_err(RotorError::PeriodicLinear)?;
        finite(factor.growth)?;
        poll(&mut cancelled)?;
        let condition = finite(factor.condition_1(&equation))?;
        poll(&mut cancelled)?;
        let norm_equation = norm_one(&equation, d)?;
        let norm_period = norm_one(&self.monodromy, d)?;
        let resolvent_gain_estimate =
            finite((condition / norm_equation) * finite(1.0 + norm_period)?)?;
        if condition <= 0.0 || resolvent_gain_estimate <= 0.0 {
            return Err(RotorError::Invalid(
                "periodic resolvent estimate is unrepresentable",
            ));
        }
        if resolvent_gain_estimate > maximum_resolvent_gain {
            return Err(RotorError::PeriodicOrbitUnresolved {
                estimated_gain: resolvent_gain_estimate,
                limit: maximum_resolvent_gain,
            });
        }
        let mut balanced = self.offset.clone();
        factor.solve(&mut balanced);
        poll(&mut cancelled)?;
        let initial = decode(&balanced, &self.state_scale, &mut cancelled)?;
        let original = state_values(&initial, self.report.dofs)?;
        let relative_fixed_point_residual = affine_residual(
            &self.report.monodromy,
            &self.physical_offset,
            &original,
            &original,
            relative_residual_tolerance,
            &mut cancelled,
        )?;
        let trajectory = self.advance(&initial, relative_residual_tolerance, &mut cancelled)?;
        // A nodal or static velocity can be exactly zero physically while
        // its rounded endpoint is O(eps). Use the declared common full-state
        // norm here; the separate rowwise affine gates still resolve forces.
        let mut peak = 0.0_f64;
        for state in &trajectory.boundary_states {
            poll(&mut cancelled)?;
            for value in encode(
                &state_values(state, self.report.dofs)?,
                &self.state_scale,
                &mut cancelled,
            )? {
                peak = peak.max(value.abs());
            }
        }
        let last = state_values(
            &trajectory.boundary_states[trajectory.boundary_states.len() - 1],
            self.report.dofs,
        )?;
        let mut closure = 0.0_f64;
        for i in 0..d {
            let difference = (Dd::from_f64(last[i]) - Dd::from_f64(original[i]))
                * Dd::from_f64(self.state_scale[i]);
            closure = closure.max(finite(difference.to_f64().abs())?);
        }
        let relative_closure_residual = if peak == 0.0 {
            closure
        } else {
            finite(closure / peak)?
        };
        admit_residual(relative_closure_residual, relative_residual_tolerance)?;
        poll(&mut cancelled)?;
        Ok(RotorPeriodicOrbit {
            trajectory,
            relative_fixed_point_residual,
            relative_closure_residual,
            resolvent_gain_estimate,
        })
    }
}

fn admit_tolerance(tolerance: f64) -> Result<(), RotorError> {
    if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
        Err(RotorError::Invalid(
            "periodic residual tolerance in (0,1) required",
        ))
    } else {
        Ok(())
    }
}

fn admit_residual(relative: f64, limit: f64) -> Result<(), RotorError> {
    if relative > limit {
        Err(RotorError::Residual { relative, limit })
    } else {
        Ok(())
    }
}

fn state_values(state: &RotorTransientState, n: usize) -> Result<Vec<f64>, RotorError> {
    if state.displacement.len() != n
        || state.velocity.len() != n
        || state
            .displacement
            .iter()
            .chain(&state.velocity)
            .any(|value| !value.is_finite())
    {
        return Err(RotorError::Invalid(
            "matching finite original rotor state required",
        ));
    }
    Ok(state
        .displacement
        .iter()
        .chain(&state.velocity)
        .copied()
        .collect())
}

fn encode(
    values: &[f64],
    scale: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, RotorError> {
    let mut result = Vec::with_capacity(values.len());
    for (&value, &factor) in values.iter().zip(scale) {
        poll(cancelled)?;
        result.push(restore_entry(value, 1.0, factor)?);
    }
    Ok(result)
}

fn decode(
    values: &[f64],
    scale: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RotorTransientState, RotorError> {
    let n = values.len() / 2;
    let mut displacement = Vec::with_capacity(n);
    let mut velocity = Vec::with_capacity(n);
    for i in 0..n {
        poll(cancelled)?;
        displacement.push(restore_entry(finite(values[i])?, scale[i], 1.0)?);
        velocity.push(restore_entry(finite(values[n + i])?, scale[n + i], 1.0)?);
    }
    Ok(RotorTransientState {
        displacement,
        velocity,
    })
}

pub(super) fn apply_affine(
    matrix: &[f64],
    offset: &[f64],
    state: &[f64],
    d: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, RotorError> {
    let mut result = Vec::with_capacity(d);
    for i in 0..d {
        poll(cancelled)?;
        let mut value = Dd::from_f64(offset[i]);
        for j in 0..d {
            value = value + Dd::from_f64(matrix[i * d + j]) * Dd::from_f64(state[j]);
        }
        result.push(finite(value.to_f64())?);
    }
    Ok(result)
}

// Original, rowwise increment equation. Keeping P-I and endpoint differences
// in DD prevents the large unchanged state from setting a false residual scale.
fn affine_residual(
    matrix: &[f64],
    offset: &[f64],
    before: &[f64],
    after: &[f64],
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, RotorError> {
    let d = before.len();
    let mut maximum = 0.0_f64;
    let mut magnitudes = vec![0.0; d + 2];
    for i in 0..d {
        poll(cancelled)?;
        let increment = Dd::from_f64(after[i]) - Dd::from_f64(before[i]);
        let mut residual = increment - Dd::from_f64(offset[i]);
        magnitudes[0] = finite(increment.to_f64().abs())?;
        magnitudes[1] = offset[i].abs();
        for j in 0..d {
            let coefficient =
                Dd::from_f64(matrix[i * d + j]) - Dd::from_f64(if i == j { 1.0 } else { 0.0 });
            let term = coefficient * Dd::from_f64(before[j]);
            residual = residual - term;
            magnitudes[j + 2] = finite(term.to_f64().abs())?;
        }
        let scale = magnitudes.iter().copied().fold(0.0_f64, f64::max);
        let numerator = finite(residual.to_f64().abs())?;
        let relative = if scale == 0.0 {
            numerator
        } else {
            let denominator: f64 = magnitudes.iter().map(|value| value / scale).sum();
            finite((numerator / scale) / denominator)?
        };
        maximum = maximum.max(relative);
    }
    admit_residual(maximum, tolerance)?;
    Ok(maximum)
}

fn norm_one(matrix: &[f64], d: usize) -> Result<f64, RotorError> {
    let mut maximum = 0.0_f64;
    for j in 0..d {
        let mut sum = Dd::ZERO;
        for i in 0..d {
            sum = sum + Dd::from_f64(matrix[i * d + j].abs());
        }
        maximum = maximum.max(finite(sum.to_f64())?);
    }
    Ok(maximum)
}
