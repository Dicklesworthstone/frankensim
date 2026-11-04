//! Full-space exponential evolution of heterogeneous pressure acoustics.
//!
//! Canonical coordinates are pressure potential `psi` and its momentum
//! `pi = M_beta p`. Their storage is `diag(K_rho, M_beta^-1)`, interconnection
//! is canonical, and the momentum resistance is the ORIGINAL real boundary
//! admittance `C_Y`. A diagonal change of units conditions these matrices;
//! neither the constant-potential gauge nor any other FE coordinate is removed.
//!
//! The existing pHS exponential and Gramian owner evolves this complete linear
//! system under held inward-velocity profiles. Source work uses the integrated
//! conjugate pressure, and absorption integrates `p^T C_Y p` on that SAME
//! trajectory. It is not the energy deficit or an endpoint quadrature rule.
//! Preparation is dense and explicitly bounded. Warm steps use cached matrix
//! contractions, with no factorization, material query or exponential.
//!
//! Exact refers to the semidiscrete, held-input solution up to numerical
//! exponential error. Spatial/geometry/material-quadrature convergence, input
//! bandwidth and output antialiasing remain separate obligations. The medium
//! is frozen, lossless in the bulk, without mean flow. Reactive walls require
//! their own causal states and are not replaced by a real admittance here.

use super::super::VolumeVelocitySource;
use super::{
    BoundaryAdmittance, C64, Dd, HeterogeneousError, HeterogeneousHelmholtzFe,
    HeterogeneousTimeEnergy, HeterogeneousTimeState, RealVelocityTrace, ResistiveBoundary, det,
    energy_change, invalid, poll, positive_form,
};
use fs_la::factor::cholesky;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage};

/// Fixed spatial inward-velocity shape, multiplied by a held amplitude [m/s].
#[derive(Clone, Debug)]
pub struct ExactVelocityPort {
    /// Original boundary triangles and DIMENSIONLESS linear shape values at
    /// their vertices. Each port has unique faces; different ports may overlap
    /// and interfere coherently. The original area quadrature enters once.
    pub profile: Vec<RealVelocityTrace>,
}

/// Admission of the complete spatial state and its cached exact trajectory.
#[derive(Clone, Copy, Debug)]
pub struct ExactTimeLimits {
    /// Number of original pressure coordinates; every one retains potential
    /// and momentum. The existing pHS augmented-state cap applies additionally.
    pub max_dofs: usize,
    /// Existing exponential/Gramian cache, byte and multiply-work limits.
    /// Coordinate maps additionally occupy O(max_dofs^2) storage.
    pub preparation: DenseExponentialBudget,
    /// Relative physical energy/work defect, plus evaluated quadratic roundoff.
    pub energy_balance_tolerance: f64,
}

impl Default for ExactTimeLimits {
    fn default() -> Self {
        Self {
            max_dofs: 63,
            preparation: DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: 64 * 1024 * 1024,
                max_multiply_terms: 20_000_000_000,
            },
            energy_balance_tolerance: 5e-11,
        }
    }
}

/// Physical same-trajectory energy and work [J], with independently integrated loss.
#[derive(Clone, Copy, Debug)]
pub struct ExactTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    /// Boundary and volume-source work, each integrated along this trajectory.
    pub input_work_j: f64,
    /// Signed volume-source subtotal [J], already included in input work.
    pub volume_source_work_j: f64,
    pub boundary_dissipation_j: f64,
}

impl ExactTimeWork {
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.boundary_dissipation_j - self.input_work_j
    }
}

/// Checked candidate, without mutating the caller's accepted field or clock.
#[derive(Clone, Debug)]
pub struct ExactTimeStep<'m, 'f, 'c> {
    pub state: HeterogeneousTimeState<'m, 'f, 'c>,
    pub work: ExactTimeWork,
    /// `integral_dt integral_boundary(shape * p) dA` [Pa m² s], in port order.
    /// Multiplication by each held velocity amplitude [m/s] gives its work.
    pub integrated_port_pressure_pa_m2_s: Vec<f64>,
    /// `integral_dt integral_volume(shape * p) dV` [Pa s], in source order.
    /// Multiply by each held volume-flow amplitude [m^3/s] for its work [J].
    pub integrated_volume_source_pressure_pa_s: Vec<f64>,
}

/// Immutable complete-field dynamics and fixed-clock source/absorption cache.
pub struct PreparedExactTime<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    duration_s: f64,
    potential_scale: Vec<f64>,
    mass_scale: Vec<f64>,
    normalized_mass: Vec<f64>,
    inverse_mass: Vec<f64>,
    port_count: usize,
    volume_source_count: usize,
    flow: PreparedQuadraticFlow,
    limits: ExactTimeLimits,
}

fn propagation(error: PreparedStepError) -> HeterogeneousError {
    match error {
        PreparedStepError::Cancelled => HeterogeneousError::Cancelled,
        PreparedStepError::Solver(error) => HeterogeneousError::Propagation(error),
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare the original FULL pressure field for exact held-source steps.
    /// Sources and boundary loss share the original trace quadrature. This
    /// never silently falls back to midpoint or a truncated modal basis.
    pub fn prepare_exact_time<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        ports: &[ExactVelocityPort],
        limits: ExactTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_volume_sources(
            duration_s,
            boundaries,
            ports,
            &[],
            limits,
            cancelled,
        )
    }

    /// Prepare the SAME full-state exponential with additional original-volume
    /// injection columns. Every source is bound to this exact material/geometry
    /// owner. Source shapes are reused without resampling, and their amplitudes
    /// [m^3/s] remain distinct from boundary-velocity amplitudes [m/s].
    pub fn prepare_exact_time_with_volume_sources<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        ports: &[ExactVelocityPort],
        sources: &[&VolumeVelocitySource<'_, '_, '_>],
        limits: ExactTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        super::super::source::admit_count(sources.len())?;
        let n = self.ndof();
        let port_count = ports
            .len()
            .checked_add(sources.len())
            .ok_or_else(|| invalid("exact acoustic source count overflow"))?;
        let states = n
            .checked_mul(2)
            .ok_or_else(|| invalid("exact acoustic state overflow"))?;
        let augmented = states
            .checked_add(port_count)
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(2))
            .ok_or_else(|| invalid("exact acoustic source extent overflow"))?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
            || limits.energy_balance_tolerance > 1e-6
            || limits.max_dofs == 0
        {
            return Err(invalid(
                "positive exact clock, bounded energy tolerance and full-field limits required",
            )
            .into());
        }
        if n > limits.max_dofs {
            return Err(HeterogeneousError::Limit {
                what: "exact transient pressure degrees of freedom",
                requested: n,
                maximum: limits.max_dofs,
            });
        }
        if augmented > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(HeterogeneousError::Limit {
                what: "exact transient augmented states",
                requested: augmented,
                maximum: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
            });
        }
        for source in sources {
            poll(cancelled)?;
            source.validate(self)?;
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
        let mut loads = Vec::with_capacity(port_count);
        for port in ports {
            poll(cancelled)?;
            if port.profile.is_empty() {
                return Err(
                    invalid("exact acoustic source requires a nonempty trace profile").into(),
                );
            }
            loads.push(self.geometry.boundary_linear_load(&port.profile)?);
        }
        for source in sources {
            poll(cancelled)?;
            loads.push(source.load.clone());
        }

        // pi=M p, z=D^-1 pi with D_ii=sqrt(M_ii). Thus H_z=Mhat^-1
        // where Mhat=D^-1 M D^-1. Cholesky never sees the pressure gauge.
        let mut mass_scale = Vec::with_capacity(n);
        let mut potential_scale = Vec::with_capacity(n);
        for row in 0..n {
            poll(cancelled)?;
            for (matrix, scales) in [
                (&self.mass, &mut mass_scale),
                (&self.stiffness, &mut potential_scale),
            ] {
                let (cols, values) = matrix.row(row);
                let diagonal = cols
                    .binary_search(&row)
                    .ok()
                    .map(|i| values[i])
                    .ok_or_else(|| {
                        invalid("original acoustic coordinate has no storage diagonal")
                    })?;
                let value = det::sqrt(diagonal);
                if !value.is_finite() || value <= 0.0 {
                    return Err(invalid(
                        "positive representable acoustic coordinate scales required",
                    )
                    .into());
                }
                scales.push(value);
            }
        }
        let mut mass = vec![0.0; n * n];
        for row in 0..n {
            poll(cancelled)?;
            let (cols, values) = self.mass.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                mass[row * n + col] = value / mass_scale[row] / mass_scale[col];
            }
        }
        let factor = cholesky(&mass, n).map_err(HeterogeneousError::MassFactor)?;
        let mut inverse = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for col in 0..n {
            poll(cancelled)?;
            rhs.fill(0.0);
            rhs[col] = 1.0;
            factor.solve(&mut rhs);
            for row in 0..n {
                inverse[row * n + col] = rhs[row];
            }
        }
        // Test ALL original mass entries; no lower-triangle-only substitution
        // or post-hoc symmetric/positive projection can change the source form.
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                let mut value = Dd::from_f64(0.0);
                let mut bound = 0.0;
                for k in 0..n {
                    let term = Dd::from_f64(mass[row * n + k]) * Dd::from_f64(inverse[k * n + col]);
                    value = value + term;
                    bound += term.to_f64().abs();
                }
                let expected: f64 = if row == col { 1.0 } else { 0.0 };
                if !bound.is_finite()
                    || !value.to_f64().is_finite()
                    || (value.to_f64() - expected).abs()
                        > 512.0 * f64::EPSILON * n as f64 * bound.max(expected)
                {
                    return Err(
                        invalid("original acoustic mass inverse failed its residual").into(),
                    );
                }
            }
        }
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut r = vec![0.0; states * states];
        let mut g = vec![0.0; states * port_count];
        for row in 0..n {
            poll(cancelled)?;
            let (cols, values) = self.stiffness.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                q[row * states + col] = value / potential_scale[row] / potential_scale[col];
            }
            for col in 0..n {
                q[(n + row) * states + n + col] = inverse[row * n + col];
            }
            let coupling = potential_scale[row] / mass_scale[row];
            j[row * states + n + row] = coupling;
            j[(n + row) * states + row] = -coupling;
            let (cols, values) = conductance.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                r[(n + row) * states + n + col] = value / mass_scale[row] / mass_scale[col];
            }
            for (port, load) in loads.iter().enumerate() {
                g[(n + row) * port_count + port] = load[row] / mass_scale[row];
            }
        }
        // Original loss C, pulled back through p=D^-1 Mhat^-1 z. This is
        // independent of endpoint energy and retains cross-coordinate terms.
        let mut right = vec![0.0; n * n];
        for row in 0..n {
            poll(cancelled)?;
            for k in 0..n {
                let value = r[(n + row) * states + n + k];
                if value != 0.0 {
                    for col in 0..n {
                        right[row * n + col] += value * inverse[k * n + col];
                    }
                }
            }
        }
        let mut loss = vec![0.0; states * states];
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                loss[(n + row) * states + n + col] = (0..n)
                    .map(|k| inverse[k * n + row] * right[k * n + col])
                    .sum();
            }
        }
        let observations: Vec<Vec<f64>> = loads
            .iter()
            .map(|load| {
                let mut observation = vec![0.0; states];
                for col in 0..n {
                    observation[n + col] = (0..n)
                        .map(|row| load[row] / mass_scale[row] * inverse[row * n + col])
                        .sum();
                }
                observation
            })
            .collect();
        let storage = QuadraticStorage::new(q, states).map_err(HeterogeneousError::Propagation)?;
        let system = PortHamiltonian::new(states, port_count, j, r, g, Box::new(storage))
            .map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: states,
                max_omega_dt: 1.0,
            },
        )
        .map_err(HeterogeneousError::Propagation)?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &[loss],
            &observations,
            limits.preparation,
            &mut *cancelled,
        )
        .map_err(propagation)?;
        poll(cancelled)?;
        Ok(PreparedExactTime {
            model: self,
            duration_s,
            potential_scale,
            mass_scale,
            normalized_mass: mass,
            inverse_mass: inverse,
            port_count: ports.len(),
            volume_source_count: sources.len(),
            flow,
            limits,
        })
    }
}

impl<'m, 'f, 'c> PreparedExactTime<'m, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    #[must_use]
    pub const fn port_count(&self) -> usize {
        self.port_count
    }

    #[must_use]
    pub const fn volume_source_count(&self) -> usize {
        self.volume_source_count
    }

    #[must_use]
    pub fn zero_state(&self) -> HeterogeneousTimeState<'m, 'f, 'c> {
        HeterogeneousTimeState {
            model: self.model,
            potential: vec![0.0; self.model.ndof()],
            pressure: vec![0.0; self.model.ndof()],
        }
    }

    /// Original hierarchical coefficients; no resampling or energy normalization.
    pub fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<HeterogeneousTimeState<'m, 'f, 'c>, HeterogeneousError> {
        if potential.len() != self.model.ndof()
            || pressure.len() != self.model.ndof()
            || potential.iter().chain(pressure).any(|x| !x.is_finite())
        {
            return Err(invalid(
                "exact transient state requires finite original potential/pressure coefficients",
            )
            .into());
        }
        Ok(HeterogeneousTimeState {
            model: self.model,
            potential: potential.to_vec(),
            pressure: pressure.to_vec(),
        })
    }

    fn validate(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, state.model) {
            return Err(invalid(
                "exact transient state belongs to another spatial/material operator",
            )
            .into());
        }
        Ok(())
    }

    /// Physical energy evaluated from the original sparse spatial matrices.
    pub fn energy(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
    ) -> Result<HeterogeneousTimeEnergy, HeterogeneousError> {
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }

    fn energy_with_roundoff(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(HeterogeneousTimeEnergy, f64), HeterogeneousError> {
        self.validate(state)?;
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

    /// Advance the complete field with source amplitudes held over this interval.
    /// Returned physical pressure/potential is directly usable by the existing
    /// heterogeneous probes. Foreign states, invalid inputs and cancellation
    /// refuse without publishing a partial endpoint or advancing any clock.
    pub fn step(
        &self,
        initial: &HeterogeneousTimeState<'_, '_, '_>,
        amplitude_m_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactTimeStep<'m, 'f, 'c>, HeterogeneousError> {
        self.step_with_volume_sources(initial, amplitude_m_s, &[], cancelled)
    }

    /// Evolve held boundary velocities and volume injection on one original
    /// trajectory. Source work uses its actual integrated conjugate pressure,
    /// including mutual interference and negative work under reversed driving.
    /// Every prepared source requires an explicit amplitude, including zero.
    pub fn step_with_volume_sources(
        &self,
        initial: &HeterogeneousTimeState<'_, '_, '_>,
        amplitude_m_s: &[f64],
        volume_flow_m3_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactTimeStep<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate(initial)?;
        if amplitude_m_s.len() != self.port_count
            || volume_flow_m3_s.len() != self.volume_source_count
            || amplitude_m_s
                .iter()
                .chain(volume_flow_m3_s)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite held boundary and volume-source amplitudes required",
            )
            .into());
        }
        // Preserve the zero-volume legacy path without an additional allocation.
        let combined;
        let amplitudes = if self.volume_source_count == 0 {
            amplitude_m_s
        } else {
            combined = amplitude_m_s
                .iter()
                .chain(volume_flow_m3_s)
                .copied()
                .collect::<Vec<_>>();
            &combined
        };
        let (initial_energy, initial_roundoff) = self.energy_with_roundoff(initial, cancelled)?;
        let n = self.model.ndof();
        let mut x = vec![0.0; 2 * n];
        for row in 0..n {
            poll(cancelled)?;
            x[row] = self.potential_scale[row] * initial.potential[row];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n {
                value = value
                    + Dd::from_f64(self.normalized_mass[row * n + col])
                        * Dd::from_f64(self.mass_scale[col])
                        * Dd::from_f64(initial.pressure[col]);
            }
            x[n + row] = value.to_f64();
        }
        let result = self
            .flow
            .evaluate(&x, amplitudes)
            .map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let mut next = self.zero_state();
        for row in 0..n {
            poll(cancelled)?;
            next.potential[row] = result.state[row] / self.potential_scale[row];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n {
                value = value
                    + Dd::from_f64(self.inverse_mass[row * n + col])
                        * Dd::from_f64(result.state[n + col]);
            }
            next.pressure[row] = value.to_f64() / self.mass_scale[row];
        }
        if next
            .potential
            .iter()
            .chain(&next.pressure)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("exact acoustic endpoint is unrepresentable").into());
        }
        let (final_energy, final_roundoff) = self.energy_with_roundoff(&next, cancelled)?;
        let mut input = Dd::from_f64(0.0);
        let mut volume_work = Dd::from_f64(0.0);
        let mut absolute_work = 0.0;
        for (index, (&amplitude, &integral)) in
            amplitudes.iter().zip(&result.integrals[1..]).enumerate()
        {
            let work = Dd::from_f64(amplitude) * Dd::from_f64(integral);
            input = input + work;
            if index >= self.port_count {
                volume_work = volume_work + work;
            }
            absolute_work += work.to_f64().abs();
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
        let work = ExactTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: change,
            input_work_j: input.to_f64(),
            volume_source_work_j: volume_work.to_f64(),
            boundary_dissipation_j: result.integrals[0],
        };
        let physical_scale = work
            .initial_energy_j
            .abs()
            .max(work.final_energy_j.abs())
            .max(change.abs())
            .max(absolute_work)
            .max(work.boundary_dissipation_j.abs());
        let limit_j = self.limits.energy_balance_tolerance * physical_scale
            + initial_roundoff
            + final_roundoff
            + 64.0 * f64::EPSILON * absolute_work;
        if !limit_j.is_finite()
            || !work.volume_source_work_j.is_finite()
            || !work.balance_residual_j().is_finite()
            || result.integrals.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("exact acoustic work is unrepresentable").into());
        }
        if work.balance_residual_j().abs() > limit_j || work.boundary_dissipation_j < -limit_j {
            return Err(HeterogeneousError::EnergyBalance {
                observed_j: work.balance_residual_j(),
                limit_j,
            });
        }
        poll(cancelled)?;
        Ok(ExactTimeStep {
            state: next,
            work,
            integrated_port_pressure_pa_m2_s: result.integrals[1..1 + self.port_count].to_vec(),
            integrated_volume_source_pressure_pa_s: result.integrals[1 + self.port_count..]
                .to_vec(),
        })
    }
}

#[cfg(test)]
mod tests;
