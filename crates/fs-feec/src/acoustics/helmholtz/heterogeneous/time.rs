//! Broadband acoustic waves in the SAME heterogeneous spatial operator.
//!
//! The state is pressure potential `psi` [Pa s] and actual pressure `p` [Pa]:
//! `p = psi_dot`, `v = -grad(psi)/rho`. Hence the original pressure mass and
//! stiffness give `M_beta p_dot + C_Y p + K_rho psi = b_in` and the physical
//! acoustic energy `(p^T M_beta p + psi^T K_rho psi)/2`. Treating pressure and
//! pressure rate as mechanical displacement/velocity would give a different,
//! nonphysical storage and boundary-work pair.
//!
//! Implicit midpoint eliminates to one sparse positive-definite operator
//! `M_beta + (h/2) C_Y + (h^2/4) K_rho`. The existing sparse LDLT owner factors
//! it once, including its native numeric-memory admission. Every original
//! pressure degree of freedom remains live; no modal reduction, dense mass
//! inverse or material resampling enters subsequent pulse steps.
//!
//! Sources are inward velocity traces evaluated at the time-step midpoint.
//! Boundary admittance is REAL, passive and memoryless. A frequency-dependent
//! or reactive complex impedance needs its own causal internal state and is
//! deliberately not admitted by this interface. The physical source is
//! frozen; there is no mean flow, bulk heat, vorticity mode or thermal evolution.
//!
//! Midpoint is second order and passive for these forms, but is not an exact
//! wave propagator: phase/dispersion error requires a time-step convergence
//! axis independently of geometry, h/p and material quadrature. Calls poll
//! between sparse rows and around original trace/factor/solve operations;
//! those existing owners have no internal cancellation callback. Returned
//! states are candidates: a refusal never mutates the caller's accepted state.

use super::{
    AcousticFeError, BoundaryAdmittance, HeterogeneousError, HeterogeneousHelmholtzFe, invalid,
    poll,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::direct::{DirectOrdering, LdltFactor, LdltOptions, SymbolicLdlt};
use fs_sparse::{Coo, Csr};

/// Instantaneous outward velocity law `v_out = Y p`, with `Y >= 0`.
#[derive(Clone, Debug)]
pub struct ResistiveBoundary {
    /// Original boundary triangles, unique across all groups.
    pub faces: Vec<[u32; 3]>,
    /// Specific admittance [m/(Pa s)], finite and nonnegative.
    pub admittance_m_pa_s: f64,
}

/// Boundary face and inward normal velocity [m/s] at its three vertices.
/// Values follow the supplied vertex order and represent the step MIDPOINT.
pub type RealVelocityTrace = ([u32; 3], [f64; 3]);

/// Admission for sparse preparation and physical step validation.
#[derive(Clone, Copy, Debug)]
pub struct HeterogeneousTimeLimits {
    /// Maximum original pressure degrees of freedom.
    pub max_dofs: usize,
    /// Native LDLT factor/front/scratch payload cap. Caller CSR, symbolic
    /// analysis, O(ndof) vectors and allocator/runtime overhead are additional.
    pub max_factor_bytes: usize,
    /// Componentwise residual of the original first-order physical equations.
    pub residual_tolerance: f64,
    /// Relative physical work-balance allowance, plus evaluated roundoff.
    pub energy_balance_tolerance: f64,
}

impl Default for HeterogeneousTimeLimits {
    fn default() -> Self {
        Self {
            max_dofs: 4096,
            max_factor_bytes: 64 * 1024 * 1024,
            residual_tolerance: 5e-12,
            energy_balance_tolerance: 5e-11,
        }
    }
}

/// Actual instantaneous acoustic storage [J], not a pressure-rate surrogate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeterogeneousTimeEnergy {
    /// `integral(rho |v|^2)/2 = psi^T K_rho psi/2`.
    pub kinetic_j: f64,
    /// `integral(p^2/(rho c^2))/2 = p^T M_beta p/2`.
    pub compression_j: f64,
}

impl HeterogeneousTimeEnergy {
    #[must_use]
    pub fn total_j(self) -> f64 {
        self.kinetic_j + self.compression_j
    }
}

/// Real fields bound to their ORIGINAL geometry and material operator.
/// Different time-step sizes may reuse a state on that same frozen operator.
#[derive(Clone)]
pub struct HeterogeneousTimeState<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    potential: Vec<f64>,
    pressure: Vec<f64>,
}

impl core::fmt::Debug for HeterogeneousTimeState<'_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HeterogeneousTimeState")
            .field("pressure_dofs", &self.pressure.len())
            .finish_non_exhaustive()
    }
}

impl HeterogeneousTimeState<'_, '_, '_> {
    /// Original hierarchical potential coefficients [Pa s]. A spatially
    /// constant, time-independent potential is an acoustic gauge.
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        &self.potential
    }

    /// Original hierarchical pressure coefficients [Pa]. Only vertex
    /// coordinates are point samples; higher-order coordinates are bubbles.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure
    }
}

/// Physical work over this exact discrete time window [J].
#[derive(Clone, Copy, Debug)]
pub struct HeterogeneousTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    /// Symmetric bilinear energy difference, avoiding endpoint subtraction.
    pub delta_energy_j: f64,
    /// `h * integral(p_mid v_in_mid)` with the original velocity trace.
    pub input_work_j: f64,
    /// `h * integral(Y p_mid^2)`, independently of stored-energy change.
    pub boundary_dissipation_j: f64,
    /// Maximum normalized residual of dynamics and pressure-potential kinematics.
    pub backward_error: f64,
}

impl HeterogeneousTimeWork {
    /// Signed defect; never set to zero or assigned to a heat channel.
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.boundary_dissipation_j - self.input_work_j
    }
}

/// A fully admitted candidate; accepted state and clock remain caller-owned.
#[derive(Clone, Debug)]
pub struct HeterogeneousTimeStep<'m, 'f, 'c> {
    pub state: HeterogeneousTimeState<'m, 'f, 'c>,
    pub work: HeterogeneousTimeWork,
}

/// Fixed-step sparse factor and unchanged original physical forms.
pub struct PreparedHeterogeneousTime<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    duration_s: f64,
    conductance: Csr,
    factor: LdltFactor,
    scales: Vec<f64>,
    limits: HeterogeneousTimeLimits,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare full-space midpoint evolution with the original sparse factor
    /// and trace owners. Real admittance can overlap an inward source face.
    pub fn prepare_time<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        limits: HeterogeneousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHeterogeneousTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let half_h = 0.5 * duration_s;
        let quarter_h2 = half_h * half_h;
        if !duration_s.is_finite()
            || half_h <= 0.0
            || !quarter_h2.is_finite()
            || quarter_h2 == 0.0
            || limits.max_dofs == 0
            || limits.max_factor_bytes == 0
            || !limits.residual_tolerance.is_finite()
            || limits.residual_tolerance <= 0.0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
        {
            return Err(invalid(
                "finite positive time step, solve tolerances and sparse limits required",
            )
            .into());
        }
        let n = self.ndof();
        if n > limits.max_dofs {
            return Err(HeterogeneousError::Limit {
                what: "transient pressure degrees of freedom",
                requested: n,
                maximum: limits.max_dofs,
            });
        }
        let mut admitted = Vec::with_capacity(boundaries.len());
        for boundary in boundaries {
            poll(cancelled)?;
            admitted.push(BoundaryAdmittance {
                faces: boundary.faces.clone(),
                admittance_m_pa_s: C64::from_re(boundary.admittance_m_pa_s),
            });
        }
        let (conductance, _) = self.admittance_matrices(&admitted, cancelled)?;
        let mut effective = Coo::new(n, n);
        for row in 0..n {
            poll(cancelled)?;
            for (matrix, multiplier) in [
                (&self.mass, 1.0),
                (&conductance, half_h),
                (&self.stiffness, quarter_h2),
            ] {
                let (cols, values) = matrix.row(row);
                for (&col, &value) in cols.iter().zip(values) {
                    let entry = multiplier * value;
                    if !entry.is_finite() {
                        return Err(invalid("transient pressure operator overflow").into());
                    }
                    effective.push(row, col, entry);
                }
            }
        }
        let effective = effective.assemble();
        let mut scales = vec![0.0; n];
        for (row, scale) in scales.iter_mut().enumerate() {
            poll(cancelled)?;
            let (cols, values) = effective.row(row);
            let diagonal = cols
                .binary_search(&row)
                .ok()
                .map(|i| values[i])
                .ok_or_else(|| invalid("transient pressure mass has an unsupported coordinate"))?;
            if !diagonal.is_finite() || diagonal <= 0.0 {
                return Err(
                    invalid("transient pressure diagonal must be positive and finite").into(),
                );
            }
            *scale = 1.0 / det::sqrt(diagonal);
            if !scale.is_finite() || *scale <= 0.0 {
                return Err(invalid("transient pressure equilibration is unrepresentable").into());
            }
        }
        let mut scaled = Coo::new(n, n);
        for row in 0..n {
            poll(cancelled)?;
            let (cols, values) = effective.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                if col < row {
                    continue;
                }
                let entry = (value * scales[row]) * scales[col];
                if !entry.is_finite() {
                    return Err(invalid("transient pressure equilibration overflow").into());
                }
                scaled.push(row, col, entry);
                if row != col {
                    scaled.push(col, row, entry);
                }
            }
        }
        let scaled = scaled.assemble();
        poll(cancelled)?;
        let symbolic =
            SymbolicLdlt::analyze(&scaled, DirectOrdering::Amd).map_err(AcousticFeError::from)?;
        poll(cancelled)?;
        let factor = symbolic
            .factor(
                &scaled,
                &LdltOptions {
                    max_factor_bytes: Some(limits.max_factor_bytes),
                    ..LdltOptions::default()
                },
            )
            .map_err(AcousticFeError::from)?;
        poll(cancelled)?;
        if factor.inertia().positive != n || factor.inertia().negative != 0 {
            return Err(invalid("transient pressure operator is not positive definite").into());
        }
        Ok(PreparedHeterogeneousTime {
            model: self,
            duration_s,
            conductance,
            factor,
            scales,
            limits,
        })
    }
}

impl<'m, 'f, 'c> PreparedHeterogeneousTime<'m, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Zero potential and pressure on the admitted original operator.
    #[must_use]
    pub fn zero_state(&self) -> HeterogeneousTimeState<'m, 'f, 'c> {
        HeterogeneousTimeState {
            model: self.model,
            potential: vec![0.0; self.model.ndof()],
            pressure: vec![0.0; self.model.ndof()],
        }
    }

    /// Import finite coefficients in this exact original space. No interpolation,
    /// source rebinding or energy normalization occurs.
    pub fn state(
        &self,
        pressure_potential_pa_s: &[f64],
        pressure_pa: &[f64],
    ) -> Result<HeterogeneousTimeState<'m, 'f, 'c>, HeterogeneousError> {
        let n = self.model.ndof();
        if pressure_potential_pa_s.len() != n
            || pressure_pa.len() != n
            || pressure_potential_pa_s
                .iter()
                .chain(pressure_pa)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "transient state requires finite original potential/pressure coefficients",
            )
            .into());
        }
        Ok(HeterogeneousTimeState {
            model: self.model,
            potential: pressure_potential_pa_s.to_vec(),
            pressure: pressure_pa.to_vec(),
        })
    }

    fn validate_state(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, state.model) {
            return Err(
                invalid("transient state belongs to another spatial/material operator").into(),
            );
        }
        Ok(())
    }

    /// Physical energy of the supplied state on the original operator.
    pub fn energy(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
    ) -> Result<HeterogeneousTimeEnergy, HeterogeneousError> {
        self.validate_state(state)?;
        self.energy_with_roundoff(state, &mut || false)
            .map(|pair| pair.0)
    }

    fn energy_with_roundoff(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(HeterogeneousTimeEnergy, f64), HeterogeneousError> {
        let (kinetic, rk) = positive_form(&self.model.stiffness, &state.potential, cancelled)?;
        let (compression, rm) = positive_form(&self.model.mass, &state.pressure, cancelled)?;
        Ok((
            HeterogeneousTimeEnergy {
                kinetic_j: 0.5 * kinetic,
                compression_j: 0.5 * compression,
            },
            0.5 * (rk + rm),
        ))
    }

    /// Advance one fixed time step from MIDPOINT inward velocity traces. The
    /// accepted input is immutable on cancellation, invalid input or refusal.
    pub fn step(
        &self,
        initial: &HeterogeneousTimeState<'_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousTimeStep<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate_state(initial)?;
        // The SAME original real trace projection supplies both dynamics and work.
        let load = self.model.geometry.boundary_linear_load(inward_velocity)?;
        poll(cancelled)?;
        let (initial_energy, initial_roundoff) = self.energy_with_roundoff(initial, cancelled)?;
        let mut next = self.state(&initial.potential, &initial.pressure)?;
        let mut midpoint_p = initial.pressure.clone();
        let mut midpoint_psi = initial.potential.clone();
        // One solve plus, only when needed, at most two residual corrections.
        // The initial candidate is the no-acceleration drift. Its residual is
        // the physical RHS; no copied evolution matrix supplies the check.
        for correction in 0..=3 {
            for i in 0..self.model.ndof() {
                midpoint_p[i] = f64::midpoint(initial.pressure[i], next.pressure[i]);
                next.potential[i] = self.duration_s.mul_add(midpoint_p[i], initial.potential[i]);
                midpoint_psi[i] = f64::midpoint(initial.potential[i], next.potential[i]);
            }
            if next
                .potential
                .iter()
                .chain(&next.pressure)
                .any(|x| !x.is_finite())
            {
                return Err(invalid("transient pressure candidate overflow").into());
            }
            let (residual, backward_error) = self.physical_residual(
                initial,
                &next,
                &midpoint_p,
                &midpoint_psi,
                &load,
                cancelled,
            )?;
            if backward_error <= self.limits.residual_tolerance {
                let work = self.work(
                    initial,
                    &next,
                    &midpoint_p,
                    &load,
                    initial_energy,
                    initial_roundoff,
                    backward_error,
                    cancelled,
                );
                match work {
                    Ok(work) => return Ok(HeterogeneousTimeStep { state: next, work }),
                    Err(HeterogeneousError::EnergyBalance { .. }) if correction < 3 => {}
                    Err(error) => return Err(error),
                }
            } else if correction == 3 {
                return Err(HeterogeneousError::Residual {
                    observed: backward_error,
                    limit: self.limits.residual_tolerance,
                });
            }
            let rhs = residual
                .iter()
                .zip(&self.scales)
                .map(|(&r, &s)| -r * s)
                .collect::<Vec<_>>();
            if rhs.iter().any(|x| !x.is_finite()) {
                return Err(invalid("transient pressure correction overflow").into());
            }
            poll(cancelled)?;
            let correction = self.factor.solve(&rhs);
            poll(cancelled)?;
            for (i, &delta) in correction.iter().enumerate() {
                next.pressure[i] = delta.mul_add(self.scales[i], next.pressure[i]);
            }
        }
        unreachable!("last physical residual or work check returns")
    }

    #[allow(clippy::too_many_arguments)]
    fn physical_residual(
        &self,
        initial: &HeterogeneousTimeState<'_, '_, '_>,
        next: &HeterogeneousTimeState<'_, '_, '_>,
        midpoint_p: &[f64],
        midpoint_psi: &[f64],
        load: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64), HeterogeneousError> {
        let h = Dd::from_f64(self.duration_s);
        let mut residual = vec![0.0; self.model.ndof()];
        let mut worst: f64 = 0.0;
        for row in 0..self.model.ndof() {
            poll(cancelled)?;
            let mut value = -h * Dd::from_f64(load[row]);
            let mut scale = self.duration_s * load[row].abs();
            let (cols, entries) = self.model.mass.row(row);
            for (&col, &entry) in cols.iter().zip(entries) {
                let delta = Dd::from_f64(next.pressure[col]) - Dd::from_f64(initial.pressure[col]);
                let term = Dd::from_f64(entry) * delta;
                value = value + term;
                scale += term.to_f64().abs();
            }
            for (matrix, field) in [
                (&self.conductance, midpoint_p),
                (&self.model.stiffness, midpoint_psi),
            ] {
                let (cols, entries) = matrix.row(row);
                for (&col, &entry) in cols.iter().zip(entries) {
                    let term = h * Dd::from_f64(entry) * Dd::from_f64(field[col]);
                    value = value + term;
                    scale += term.to_f64().abs();
                }
            }
            let value = value.to_f64();
            let kinematics = (Dd::from_f64(next.potential[row])
                - Dd::from_f64(initial.potential[row])
                - h * Dd::from_f64(midpoint_p[row]))
            .to_f64();
            let kinematic_scale = next.potential[row].abs()
                + initial.potential[row].abs()
                + self.duration_s * midpoint_p[row].abs();
            if !value.is_finite()
                || !scale.is_finite()
                || !kinematics.is_finite()
                || !kinematic_scale.is_finite()
            {
                return Err(invalid("unrepresentable transient physical residual").into());
            }
            worst = worst
                .max(relative(value.abs(), scale))
                .max(relative(kinematics.abs(), kinematic_scale));
            residual[row] = value;
        }
        Ok((residual, worst))
    }

    #[allow(clippy::too_many_arguments)]
    fn work(
        &self,
        initial: &HeterogeneousTimeState<'_, '_, '_>,
        next: &HeterogeneousTimeState<'_, '_, '_>,
        midpoint_p: &[f64],
        load: &[f64],
        initial_energy: HeterogeneousTimeEnergy,
        initial_roundoff: f64,
        backward_error: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousTimeWork, HeterogeneousError> {
        let (final_energy, final_roundoff) = self.energy_with_roundoff(next, cancelled)?;
        let (loss, loss_roundoff) = positive_form(&self.conductance, midpoint_p, cancelled)?;
        let mut input = Dd::from_f64(0.0);
        let mut input_absolute = 0.0;
        for (&p, &b) in midpoint_p.iter().zip(load) {
            let term = Dd::from_f64(p) * Dd::from_f64(b) * Dd::from_f64(self.duration_s);
            input = input + term;
            input_absolute += term.to_f64().abs();
        }
        let change = energy_change(
            &self.model.mass,
            &initial.pressure,
            &next.pressure,
            cancelled,
        )? + energy_change(
            &self.model.stiffness,
            &initial.potential,
            &next.potential,
            cancelled,
        )?;
        let work = HeterogeneousTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: change,
            input_work_j: input.to_f64(),
            boundary_dissipation_j: self.duration_s * loss,
            backward_error,
        };
        let physical_scale = change.abs() + work.input_work_j.abs() + work.boundary_dissipation_j;
        let roundoff = initial_roundoff
            + final_roundoff
            + self.duration_s * loss_roundoff
            + 64.0 * f64::EPSILON * input_absolute;
        let limit_j = self.limits.energy_balance_tolerance * physical_scale + roundoff;
        let residual = work.balance_residual_j();
        if !work.initial_energy_j.is_finite()
            || !work.final_energy_j.is_finite()
            || !residual.is_finite()
            || !limit_j.is_finite()
        {
            return Err(invalid("unrepresentable transient acoustic work").into());
        }
        if residual.abs() > limit_j {
            return Err(HeterogeneousError::EnergyBalance {
                observed_j: residual,
                limit_j,
            });
        }
        Ok(work)
    }
}

fn relative(defect: f64, scale: f64) -> f64 {
    if defect == 0.0 {
        0.0
    } else if scale == 0.0 {
        f64::INFINITY
    } else {
        defect / scale
    }
}

// Evaluate the original real quadratic form with a roundoff bound. Only a
// signed negative residue within that bound is zeroed, never physical heat.
fn positive_form(
    matrix: &Csr,
    field: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), HeterogeneousError> {
    let mut sum = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (row, &value) in field.iter().enumerate() {
        poll(cancelled)?;
        let (cols, entries) = matrix.row(row);
        for (&col, &entry) in cols.iter().zip(entries) {
            let term = Dd::from_f64(value) * Dd::from_f64(entry) * Dd::from_f64(field[col]);
            sum = sum + term;
            absolute += term.to_f64().abs();
        }
    }
    let value = sum.to_f64();
    let roundoff = 64.0 * f64::EPSILON * absolute;
    if !value.is_finite() || !absolute.is_finite() || value < -roundoff {
        return Err(invalid("negative or unrepresentable transient acoustic energy").into());
    }
    Ok((value.max(0.0), roundoff))
}

fn energy_change(
    matrix: &Csr,
    initial: &[f64],
    next: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, HeterogeneousError> {
    let mut change = Dd::from_f64(0.0);
    for row in 0..initial.len() {
        poll(cancelled)?;
        let delta = Dd::from_f64(next[row]) - Dd::from_f64(initial[row]);
        let (cols, entries) = matrix.row(row);
        for (&col, &entry) in cols.iter().zip(entries) {
            let sum = Dd::from_f64(next[col]) + Dd::from_f64(initial[col]);
            change = change + Dd::from_f64(0.5) * delta * Dd::from_f64(entry) * sum;
        }
    }
    let change = change.to_f64();
    if !change.is_finite() {
        return Err(invalid("transient acoustic energy change overflow").into());
    }
    Ok(change)
}

#[cfg(test)]
mod tests;
