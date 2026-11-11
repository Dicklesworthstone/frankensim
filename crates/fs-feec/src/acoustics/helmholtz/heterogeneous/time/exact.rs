//! Complete-field exponential evolution of heterogeneous pressure acoustics.
//!
//! Canonical coordinates are pressure potential `psi` and its momentum
//! `pi = M_beta p`. Their storage is `diag(K_rho, M_beta^-1)`, interconnection
//! is canonical, and the momentum resistance is the ORIGINAL real boundary
//! admittance `C_Y`. Optional ideal pressure-release traces restrict these
//! original forms BEFORE the mass inverse: `pi_f = M_ff p_f`. Every remaining
//! coordinate is evolved, and the complete physical state is recovered by zero
//! lifting. With no release faces the constant-potential gauge remains present.
//! A diagonal change of units conditions these same physical matrices.
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

use super::super::{BoundaryVelocitySource, VolumeVelocitySource};
use super::{
    Dd, HeterogeneousError, HeterogeneousHelmholtzFe, HeterogeneousTimeEnergy,
    HeterogeneousTimeState, PressureReleaseConstraints, RealVelocityTrace, ResistiveBoundary, det,
    energy_change, invalid, poll, positive_form,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
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
    /// Number of unconstrained original pressure coordinates; every one retains
    /// potential and momentum. The existing pHS augmented-state cap applies
    /// additionally. Original sparse forms and O(original_ndof) constraint maps
    /// and lifted state/work vectors are additional to the dense cache.
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
    /// `integral_dt integral_boundary(shape * p) dA` [Pa m² s], with linear
    /// velocity ports first, followed by prepared boundary sources in input order.
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
    pressure_release: PressureReleaseConstraints,
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

fn canonical_extent(
    dofs: usize,
    source_count: usize,
    maximum_dofs: usize,
) -> Result<usize, HeterogeneousError> {
    let states = dofs
        .checked_mul(2)
        .ok_or_else(|| invalid("exact acoustic state overflow"))?;
    let augmented = states
        .checked_add(source_count)
        .and_then(|v| v.checked_add(1))
        .and_then(|v| v.checked_mul(2))
        .ok_or_else(|| invalid("exact acoustic source extent overflow"))?;
    if dofs > maximum_dofs {
        return Err(HeterogeneousError::Limit {
            what: "exact transient pressure degrees of freedom",
            requested: dofs,
            maximum: maximum_dofs,
        });
    }
    if augmented > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
        return Err(HeterogeneousError::Limit {
            what: "exact transient augmented states",
            requested: augmented,
            maximum: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        });
    }
    Ok(states)
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
        self.prepare_exact_time_with_pressure_release(
            duration_s,
            boundaries,
            ports,
            sources,
            &[],
            limits,
            cancelled,
        )
    }

    /// Prepare exact held-input evolution with ideal zero-pressure faces.
    /// Every original free pressure/potential coordinate is retained. Original
    /// mass, stiffness, real boundary loss and both source families are restricted
    /// before forming the canonical momentum and inverting the consistent mass;
    /// slicing a full-space inverse would define different fluid inertia.
    ///
    /// This path selects `psi = p = 0` on every hierarchical release trace.
    /// Imported states, including those from another prepared owner on this
    /// model, must satisfy both exactly; no stationary flow component is erased.
    /// The boundary law belongs to this preparation, not to a state vector.
    /// Returned states use complete original coordinates and the existing probes.
    ///
    /// Release faces cannot also own a real admittance or a velocity port.
    /// Adjacent faces may share constrained edge/vertex coordinates. Original
    /// volume sources are reused without resampling. Source work and passive
    /// loss use the existing same-trajectory observations and Gramian; the zero
    /// pressure boundary reactions do no work. Empty faces retain the original
    /// arithmetic and cheap full-space admission. All dense limits count the
    /// free canonical state, with no modal truncation or propagator fallback.
    /// This ideal boundary supplies no radiation or opening end correction.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_pressure_release<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        ports: &[ExactVelocityPort],
        sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_boundary_sources(
            duration_s,
            boundaries,
            ports,
            &[],
            sources,
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Prepare the same exact trajectory with cached quadratic boundary motion.
    /// The prepared source's original full-space load is reused without face
    /// resampling, corner-only reconstruction or area normalization. It must
    /// belong to this spatial/material operator. Pressure-release conflicts are
    /// refused even for an identically zero shape; separate driving columns may
    /// overlap each other and the real-admittance faces.
    ///
    /// Linear velocity ports come first, followed by `boundary_sources` in the
    /// same order. Both use held amplitudes [m/s] in [`PreparedExactTime::step`]
    /// and its volume-source variant. Their observations have units [Pa m² s].
    /// Volume-flow amplitudes [m³/s] and observations [Pa s] stay separate.
    /// Every column enters the existing augmented-state budget, forcing matrix
    /// and conjugate-work observation exactly once, without extra dynamic states.
    /// The pressure-release, strict state-admission and dense-cache semantics of
    /// [`Self::prepare_exact_time_with_pressure_release`] are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_boundary_sources<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactTime<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_exact_time_with_spatial_boundaries(
            duration_s,
            boundaries,
            &[],
            ports,
            boundary_sources,
            sources,
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Exact complete-field propagation with spatial memoryless wall losses.
    ///
    /// The same pointwise conductance form enters the dynamics and the
    /// trajectory-integrated dissipation. Uniform, spatial and release faces
    /// must be disjoint, and every spatial field must belong to this exact
    /// model. Held boundary/volume sources, original trace coordinates and
    /// pressure-release constraints retain their existing semantics. No
    /// frequency sample is converted into a broadband material law.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_exact_time_with_spatial_boundaries<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        spatial_boundaries: &[&super::SpatialResistiveBoundary<'_, '_, '_>],
        ports: &[ExactVelocityPort],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: ExactTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        super::super::boundary_source::admit_count(boundary_sources.len())?;
        super::super::source::admit_count(sources.len())?;
        let boundary_port_count = ports
            .len()
            .checked_add(boundary_sources.len())
            .ok_or_else(|| invalid("exact acoustic boundary source count overflow"))?;
        let port_count = boundary_port_count
            .checked_add(sources.len())
            .ok_or_else(|| invalid("exact acoustic source count overflow"))?;
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
        // Refuse an oversized legacy full field before constructing even its
        // identity map. With release faces, only the true free extent is capped.
        if pressure_release_faces.is_empty() {
            canonical_extent(self.ndof(), port_count, limits.max_dofs)?;
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let n = pressure_release.free_dofs.len();
        let states = canonical_extent(n, port_count, limits.max_dofs)?;
        for source in boundary_sources {
            poll(cancelled)?;
            source.validate(self)?;
            if !pressure_release.faces.is_empty() {
                for face in &source.faces {
                    poll(cancelled)?;
                    if pressure_release.faces.contains(face) {
                        return Err(invalid(
                            "pressure-release and prepared boundary-source faces must not overlap",
                        )
                        .into());
                    }
                }
            }
        }
        for source in sources {
            poll(cancelled)?;
            source.validate(self)?;
        }
        let conductance =
            self.resistive_matrix(boundaries, spatial_boundaries, &pressure_release, cancelled)?;
        let mut loads = Vec::with_capacity(port_count);
        for port in ports {
            poll(cancelled)?;
            if port.profile.is_empty() {
                return Err(
                    invalid("exact acoustic source requires a nonempty trace profile").into(),
                );
            }
            if !pressure_release.faces.is_empty() {
                for &(face, _) in &port.profile {
                    poll(cancelled)?;
                    let mut key = face;
                    key.sort_unstable();
                    if pressure_release.faces.contains(&key) {
                        return Err(invalid(
                            "pressure-release and exact velocity-port faces must not overlap",
                        )
                        .into());
                    }
                }
            }
            loads.push(self.geometry.boundary_linear_load(&port.profile)?);
        }
        for source in boundary_sources {
            poll(cancelled)?;
            loads.push(source.load.clone());
        }
        for source in sources {
            poll(cancelled)?;
            loads.push(source.load.clone());
        }

        // pi_f=M_ff p_f, z=D^-1 pi_f with D_ii=sqrt(M_ff,ii).
        // Thus H_z=Mhat_ff^-1, where Mhat_ff=D^-1 M_ff D^-1.
        // Assemble the restricted original entries directly; never invert the
        // full mass and then slice its inverse. Cholesky sees no pressure gauge.
        let mut mass_scale = Vec::with_capacity(n);
        let mut potential_scale = Vec::with_capacity(n);
        for &original_row in &pressure_release.free_dofs {
            poll(cancelled)?;
            for (matrix, scales) in [
                (&self.mass, &mut mass_scale),
                (&self.stiffness, &mut potential_scale),
            ] {
                let (cols, values) = matrix.row(original_row);
                let diagonal = cols
                    .binary_search(&original_row)
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
        for (row, &original_row) in pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (cols, values) = self.mass.row(original_row);
            for (&col, &value) in cols.iter().zip(values) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    mass[row * n + free_col] = value / mass_scale[row] / mass_scale[free_col];
                }
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
        // Test ALL retained original mass entries; no lower-triangle-only
        // substitution or post-hoc symmetric/positive projection can change them.
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
        for (row, &original_row) in pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (cols, values) = self.stiffness.row(original_row);
            for (&col, &value) in cols.iter().zip(values) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    q[row * states + free_col] =
                        value / potential_scale[row] / potential_scale[free_col];
                }
            }
            for col in 0..n {
                q[(n + row) * states + n + col] = inverse[row * n + col];
            }
            let coupling = potential_scale[row] / mass_scale[row];
            j[row * states + n + row] = coupling;
            j[(n + row) * states + row] = -coupling;
            let (cols, values) = conductance.row(original_row);
            for (&col, &value) in cols.iter().zip(values) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    r[(n + row) * states + n + free_col] =
                        value / mass_scale[row] / mass_scale[free_col];
                }
            }
            for (port, load) in loads.iter().enumerate() {
                g[(n + row) * port_count + port] = load[original_row] / mass_scale[row];
            }
        }
        // Original free loss C_ff pulled back through p_f=D^-1 Mhat_ff^-1 z,
        // independently of endpoint energy and retaining cross-coordinate terms.
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
                        .map(|row| {
                            load[pressure_release.free_dofs[row]] / mass_scale[row]
                                * inverse[row * n + col]
                        })
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
            pressure_release,
            port_count: boundary_port_count,
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

    /// Total boundary columns: original linear ports followed by prepared sources.
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
    /// Both potential and pressure must vanish on this owner's release trace.
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
        self.admit_pressure_release_trace(potential, pressure)?;
        Ok(HeterogeneousTimeState {
            model: self.model,
            potential: potential.to_vec(),
            pressure: pressure.to_vec(),
        })
    }

    fn admit_pressure_release_trace(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<(), HeterogeneousError> {
        if self.pressure_release.faces.is_empty() {
            return Ok(());
        }
        for (dof, &fixed) in self.pressure_release.constrained.iter().enumerate() {
            if fixed && (potential[dof] != 0.0 || pressure[dof] != 0.0) {
                return Err(invalid(
                    "exact pressure-release state requires zero potential and pressure on every trace DOF",
                )
                .into());
            }
        }
        Ok(())
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
        // Private state constructors already checked full sizes and finite
        // coefficients. Only compatibility with this boundary law is new.
        self.admit_pressure_release_trace(&state.potential, &state.pressure)
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
    /// Supply all linear-port amplitudes followed by all prepared boundary-source
    /// amplitudes, in their preparation order; every amplitude has units m/s.
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
        let n = self.pressure_release.free_dofs.len();
        let mut x = vec![0.0; 2 * n];
        for (row, &original_row) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            x[row] = self.potential_scale[row] * initial.potential[original_row];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n {
                value = value
                    + Dd::from_f64(self.normalized_mass[row * n + col])
                        * Dd::from_f64(self.mass_scale[col])
                        * Dd::from_f64(initial.pressure[self.pressure_release.free_dofs[col]]);
            }
            x[n + row] = value.to_f64();
        }
        let result = self
            .flow
            .evaluate(&x, amplitudes)
            .map_err(HeterogeneousError::Propagation)?;
        poll(cancelled)?;
        let mut next = self.zero_state();
        for (row, &original_row) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            next.potential[original_row] = result.state[row] / self.potential_scale[row];
            let mut value = Dd::from_f64(0.0);
            for col in 0..n {
                value = value
                    + Dd::from_f64(self.inverse_mass[row * n + col])
                        * Dd::from_f64(result.state[n + col]);
            }
            next.pressure[original_row] = value.to_f64() / self.mass_scale[row];
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
mod boundary_source_tests;
#[cfg(test)]
mod pressure_release_tests;
#[cfg(test)]
mod tests;
