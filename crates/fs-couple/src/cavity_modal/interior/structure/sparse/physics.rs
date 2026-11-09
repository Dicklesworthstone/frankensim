//! Sparse midpoint dynamics in original work-conjugate coordinates.
use fs_math::{dd::Dd, det};
use fs_sparse::direct::{
    DirectLimits, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt,
};
use fs_sparse::{Coo, Csr};

mod relaxation;
pub(super) use relaxation::BranchSource;
pub use relaxation::{SparseStructuralMemoryWork, SparseStructuralRelaxationLimits};

/// Limits for original sparse forms and the existing native LDLT owner.
#[derive(Clone, Copy, Debug)]
pub struct SparseStructuralTimeLimits {
    /// Structural plus pressure coordinates; each retains position and velocity.
    pub maximum_coordinates: usize,
    /// Sum of retained M/C/K/B entries, and each assembly staging extent.
    /// Source geometry, caller matrices and allocator metadata are additional.
    pub maximum_matrix_entries: usize,
    /// Number of authored normal-trace face terms and resistive faces.
    pub maximum_trace_terms: usize,
    /// Numeric factor/front/scratch payload, applied to each native factor.
    pub maximum_factor_bytes: usize,
    /// Native symbolic/fill/frontal/scalar-work ceilings, per factorization.
    pub factor: DirectLimits,
    /// Original componentwise momentum and kinematic equation tolerance.
    pub maximum_relative_residual: f64,
    /// Relative tolerance of EACH physical work account, plus roundoff.
    pub energy_balance_tolerance: f64,
}
impl Default for SparseStructuralTimeLimits {
    fn default() -> Self {
        Self {
            maximum_coordinates: 16_384,
            maximum_matrix_entries: 4_000_000,
            maximum_trace_terms: 1_000_000,
            maximum_factor_bytes: 128 * 1024 * 1024,
            factor: DirectLimits {
                max_symbolic_entries: 8_000_000,
                max_analysis_work: 100_000_000,
                max_numeric_work: 1_000_000_000,
                max_front_dim: 2048,
            },
            maximum_relative_residual: 5e-11,
            energy_balance_tolerance: 5e-10,
        }
    }
}

/// No partial factor or candidate state escapes a refusal.
#[derive(Debug)]
pub enum SparseStructuralTimeError {
    Invalid(&'static str),
    Limit {
        what: &'static str,
        requested: usize,
        maximum: usize,
    },
    Factor(LdltError),
    Cancelled,
    Residual {
        observed: f64,
        maximum: f64,
    },
    EnergyBalance {
        account: &'static str,
        residual_j: f64,
        maximum_j: f64,
    },
}
impl core::fmt::Display for SparseStructuralTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "sparse structural interior: {s}"),
            Self::Limit {
                what,
                requested,
                maximum,
            } => write!(
                f,
                "sparse structural interior {what}: {requested} exceeds {maximum}"
            ),
            Self::Factor(e) => e.fmt(f),
            Self::Cancelled => f.write_str("sparse structural interior cancelled"),
            Self::Residual { observed, maximum } => write!(
                f,
                "sparse structural interior residual {observed} exceeds {maximum}"
            ),
            Self::EnergyBalance {
                account,
                residual_j,
                maximum_j,
            } => write!(
                f,
                "sparse structural interior {account} defect {residual_j} J exceeds {maximum_j} J"
            ),
        }
    }
}
impl std::error::Error for SparseStructuralTimeError {}
impl From<LdltError> for SparseStructuralTimeError {
    fn from(e: LdltError) -> Self {
        if matches!(e, LdltError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Factor(e)
        }
    }
}
type Result<T> = core::result::Result<T, SparseStructuralTimeError>;
const ROUND: f64 = 4096.0 * f64::EPSILON;

impl SparseStructuralTimeLimits {
    pub(super) fn admit(&self, ns: usize, na: usize) -> Result<usize> {
        let n = ns
            .checked_add(na)
            .ok_or(SparseStructuralTimeError::Invalid(
                "coordinate extent overflow",
            ))?;
        if ns == 0
            || na == 0
            || self.maximum_matrix_entries == 0
            || self.maximum_trace_terms == 0
            || self.maximum_factor_bytes == 0
            || !self.maximum_relative_residual.is_finite()
            || self.maximum_relative_residual <= 0.0
            || self.maximum_relative_residual > 1e-6
            || !self.energy_balance_tolerance.is_finite()
            || self.energy_balance_tolerance <= 0.0
            || self.energy_balance_tolerance > 1e-6
        {
            return Err(SparseStructuralTimeError::Invalid(
                "positive bounded dimensions, payload and physical tolerances required",
            ));
        }
        bound("coordinates", n, self.maximum_coordinates)?;
        Ok(n)
    }
}

pub(super) fn bound(what: &'static str, requested: usize, maximum: usize) -> Result<()> {
    if requested > maximum {
        Err(SparseStructuralTimeError::Limit {
            what,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}
pub(super) fn poll(cancelled: &impl Fn() -> bool) -> Result<()> {
    if cancelled() {
        Err(SparseStructuralTimeError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original [structural, acoustic] forms; no lumping or modal coordinates.
pub(super) struct SparseForms {
    pub mass: [Csr; 2],
    pub stiffness: [Csr; 2],
    pub damping: [Csr; 2],
    /// Acoustic rows, structural columns: integral(phi_a * n.u_s).
    pub coupling: Csr,
}

/// Original sparse forms and a cached symmetrically signed midpoint factor.
pub(super) struct SparseMidpoint {
    forms: SparseForms,
    factor: LdltFactor,
    scales: Vec<f64>,
    duration_s: f64,
    limits: SparseStructuralTimeLimits,
    branches: relaxation::Branches,
}

/// Every original coefficient and a clock, bound to its exact preparation.
#[derive(Clone)]
pub struct SparseStructuralTimeState<'p> {
    owner: &'p SparseMidpoint,
    time_s: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
    stress: Vec<f64>,
}
impl core::fmt::Debug for SparseStructuralTimeState<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SparseStructuralTimeState")
            .field("time_s", &self.time_s)
            .field("structural_dofs", &self.owner.structural_dofs())
            .field("pressure_dofs", &self.owner.pressure_dofs())
            .finish_non_exhaustive()
    }
}
impl SparseStructuralTimeState<'_> {
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub fn structural_displacement(&self) -> &[f64] {
        &self.position[..self.owner.structural_dofs()]
    }
    #[must_use]
    pub fn structural_velocity(&self) -> &[f64] {
        &self.velocity[..self.owner.structural_dofs()]
    }
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        &self.position[self.owner.structural_dofs()..]
    }
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.velocity[self.owner.structural_dofs()..]
    }
    /// Actual original Maxwell stresses, in authored branch order.
    #[must_use]
    pub fn branch_stresses(&self) -> &[f64] {
        &self.stress
    }
}

/// Physical storage; pressure is the velocity of pressure potential.
#[derive(Clone, Copy, Debug)]
pub struct SparseStructuralTimeEnergy {
    pub structural_kinetic_j: f64,
    pub structural_elastic_j: f64,
    pub acoustic_kinetic_j: f64,
    pub acoustic_compression_j: f64,
    /// Actual sigma^2/(2k); relaxed spatial stiffness excludes these branches.
    pub relaxation_j: f64,
}
impl SparseStructuralTimeEnergy {
    #[must_use]
    pub fn structural_j(self) -> f64 {
        self.structural_kinetic_j + self.structural_elastic_j
    }
    #[must_use]
    pub fn acoustic_j(self) -> f64 {
        self.acoustic_kinetic_j + self.acoustic_compression_j
    }
    #[must_use]
    pub fn total_j(self) -> f64 {
        self.structural_j() + self.acoustic_j() + self.relaxation_j
    }
}

/// Independently contracted work on the same midpoint trajectory.
#[derive(Clone, Copy, Debug)]
pub struct SparseStructuralTimeWork {
    pub initial_energy: SparseStructuralTimeEnergy,
    pub final_energy: SparseStructuralTimeEnergy,
    /// Evaluated by polarized original forms, avoiding subtraction of large totals.
    pub structural_energy_change_j: f64,
    pub acoustic_energy_change_j: f64,
    pub structural_input_work_j: f64,
    pub acoustic_input_work_j: f64,
    /// Positive from fluid to structure; independently contracted from B^T p.
    pub fluid_to_structure_work_j: f64,
    /// The reciprocal contraction from p and B v. It must match the above.
    pub acoustic_interface_work_j: f64,
    pub structural_dissipation_j: f64,
    pub wall_dissipation_j: f64,
    pub relaxation_energy_change_j: f64,
    pub structural_to_relaxation_work_j: f64,
    pub acoustic_to_relaxation_work_j: f64,
    pub relaxation_dissipation_j: f64,
    pub maximum_relative_residual: f64,
}
impl SparseStructuralTimeWork {
    #[must_use]
    pub fn structural_balance_residual_j(self) -> f64 {
        (Dd::from_f64(self.structural_energy_change_j)
            - Dd::from_f64(self.structural_input_work_j)
            - Dd::from_f64(self.fluid_to_structure_work_j)
            + Dd::from_f64(self.structural_dissipation_j)
            + Dd::from_f64(self.structural_to_relaxation_work_j))
        .to_f64()
    }
    #[must_use]
    pub fn acoustic_balance_residual_j(self) -> f64 {
        (Dd::from_f64(self.acoustic_energy_change_j) - Dd::from_f64(self.acoustic_input_work_j)
            + Dd::from_f64(self.acoustic_interface_work_j)
            + Dd::from_f64(self.wall_dissipation_j)
            + Dd::from_f64(self.acoustic_to_relaxation_work_j))
        .to_f64()
    }
    #[must_use]
    pub fn interface_balance_residual_j(self) -> f64 {
        self.fluid_to_structure_work_j - self.acoustic_interface_work_j
    }
    #[must_use]
    pub fn relaxation_balance_residual_j(self) -> f64 {
        (Dd::from_f64(self.relaxation_energy_change_j)
            + Dd::from_f64(self.relaxation_dissipation_j)
            - Dd::from_f64(self.structural_to_relaxation_work_j)
            - Dd::from_f64(self.acoustic_to_relaxation_work_j))
        .to_f64()
    }
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        (Dd::from_f64(self.structural_energy_change_j)
            + Dd::from_f64(self.acoustic_energy_change_j)
            - Dd::from_f64(self.structural_input_work_j)
            - Dd::from_f64(self.acoustic_input_work_j)
            + Dd::from_f64(self.structural_dissipation_j)
            + Dd::from_f64(self.wall_dissipation_j)
            + Dd::from_f64(self.relaxation_energy_change_j)
            + Dd::from_f64(self.relaxation_dissipation_j))
        .to_f64()
    }
}

/// A checked candidate; the input state remains unchanged even on refusal.
#[derive(Debug)]
pub struct SparseStructuralTimeStep<'p> {
    pub state: SparseStructuralTimeState<'p>,
    pub work: SparseStructuralTimeWork,
    /// Midpoint original generalized fluid forces B^T p [force units].
    pub fluid_force: Vec<f64>,
    /// Midpoint original pressure-basis outward volume-flow load B v [m³/s].
    pub displaced_flow_m3_s: Vec<f64>,
    /// Separate original branch storage, work and heat, in source order.
    pub memory: Vec<SparseStructuralMemoryWork>,
}

impl SparseMidpoint {
    pub(super) fn new(
        forms: SparseForms,
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self> {
        Self::new_relaxing(
            forms,
            duration_s,
            limits,
            &[],
            SparseStructuralRelaxationLimits::default(),
            cancelled,
        )
    }

    pub(super) fn new_relaxing(
        forms: SparseForms,
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        sources: &[BranchSource<'_>],
        memory_limits: SparseStructuralRelaxationLimits,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let ns = forms.mass[0].nrows();
        let na = forms.mass[1].nrows();
        let n = limits.admit(ns, na)?;
        let t = 0.5 * duration_s;
        let tt = t * t;
        if !duration_s.is_finite() || t <= 0.0 || !tt.is_finite() || tt == 0.0 {
            return Err(SparseStructuralTimeError::Invalid(
                "finite representable positive midpoint duration required",
            ));
        }
        let mut entries = forms.coupling.nnz();
        for (d, size) in [ns, na].into_iter().enumerate() {
            for a in [&forms.mass[d], &forms.stiffness[d], &forms.damping[d]] {
                if a.nrows() != size || a.ncols() != size {
                    return Err(SparseStructuralTimeError::Invalid(
                        "one complete original space per matrix block required",
                    ));
                }
                entries =
                    entries
                        .checked_add(a.nnz())
                        .ok_or(SparseStructuralTimeError::Invalid(
                            "matrix entry extent overflow",
                        ))?;
            }
        }
        bound(
            "retained matrix entries",
            entries,
            limits.maximum_matrix_entries,
        )?;
        if forms.coupling.nrows() != na || forms.coupling.ncols() != ns {
            return Err(SparseStructuralTimeError::Invalid(
                "coupling must retain acoustic rows and structural columns",
            ));
        }
        for row in 0..na {
            poll(cancelled)?;
            if forms.coupling.row(row).1.iter().any(|x| !x.is_finite()) {
                return Err(SparseStructuralTimeError::Invalid(
                    "finite original interface coefficients required",
                ));
            }
        }
        for d in 0..2 {
            admit_form(&forms.mass[d], true, &limits, cancelled)?;
            admit_form(&forms.stiffness[d], false, &limits, cancelled)?;
            admit_form(&forms.damping[d], false, &limits, cancelled)?;
        }
        let branches = relaxation::Branches::prepare(
            sources,
            n,
            ns,
            t,
            entries,
            &limits,
            memory_limits,
            cancelled,
        )?;
        let dimension = n
            .checked_add(branches.len())
            .ok_or(SparseStructuralTimeError::Invalid(
                "augmented material coordinate extent overflow",
            ))?;
        // Midpoint velocities [v_mid,p_mid] are simultaneous unknowns:
        // [ A_s  -t B^T ][v] = [ M_s v0 + t F - t K_s u0 ]
        // [ -t B -A_a   ][p]   [-M_a p0 - t s + t K_a psi0].
        // Multiplying ONLY the acoustic equation by -1 makes this symmetric.
        // Neither K nor the acoustic static operator is inverted.
        let mut effective = Coo::new(dimension, dimension);
        for (d, size) in [ns, na].into_iter().enumerate() {
            let offset = if d == 0 { 0 } else { ns };
            let sign = if d == 0 { 1.0 } else { -1.0 };
            for (matrix, weight) in [
                (&forms.mass[d], sign),
                (&forms.damping[d], sign * t),
                (&forms.stiffness[d], sign * tt),
            ] {
                for row in 0..size {
                    poll(cancelled)?;
                    let (cols, vals) = matrix.row(row);
                    for (&col, &value) in cols.iter().zip(vals) {
                        push(
                            &mut effective,
                            offset + row,
                            offset + col,
                            weight * value,
                            &limits,
                        )?;
                    }
                }
            }
        }
        for row in 0..na {
            poll(cancelled)?;
            let (cols, vals) = forms.coupling.row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                push(&mut effective, ns + row, col, -t * value, &limits)?;
                push(&mut effective, col, ns + row, -t * value, &limits)?;
            }
        }
        // Scalar Maxwell histories keep a sparse augmented column, never b*b^T.
        branches.assemble(&mut effective, n, &limits, cancelled)?;
        let mut effective = effective.assemble();
        let mut scales = Vec::with_capacity(dimension);
        for row in 0..dimension {
            poll(cancelled)?;
            let diagonal = effective.get(row, row);
            let sign = if row < ns {
                1.0
            } else if row < n {
                -1.0
            } else {
                branches.sign(row - n)
            };
            if !diagonal.is_finite() || sign * diagonal <= 0.0 {
                return Err(SparseStructuralTimeError::Invalid(
                    "original midpoint block diagonal is not definite",
                ));
            }
            let scale = 1.0 / det::sqrt(diagonal.abs());
            if !scale.is_finite() || scale <= 0.0 {
                return Err(SparseStructuralTimeError::Invalid(
                    "midpoint congruence scale is unrepresentable",
                ));
            }
            scales.push(scale);
        }
        scale_form(&mut effective, &scales, cancelled)?;
        let factor = factor(&effective, &limits, cancelled)?;
        if factor.inertia().positive != ns + branches.acoustic
            || factor.inertia().negative != na + branches.len() - branches.acoustic
        {
            return Err(SparseStructuralTimeError::Invalid(
                "coupled midpoint factor has unexpected inertia",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            forms,
            factor,
            scales,
            duration_s,
            limits,
            branches,
        })
    }
    pub(super) fn branch_count(&self) -> usize {
        self.branches.len()
    }
    pub(super) fn structural_dofs(&self) -> usize {
        self.forms.mass[0].nrows()
    }
    pub(super) fn pressure_dofs(&self) -> usize {
        self.forms.mass[1].nrows()
    }
    pub(super) const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub(super) fn coupling(&self) -> &Csr {
        &self.forms.coupling
    }

    pub(super) fn zero_state(&self) -> SparseStructuralTimeState<'_> {
        let n = self.structural_dofs() + self.pressure_dofs();
        SparseStructuralTimeState {
            owner: self,
            time_s: 0.0,
            position: vec![0.0; n],
            velocity: vec![0.0; n],
            stress: vec![0.0; self.branches.len()],
        }
    }
    pub(super) fn state(
        &self,
        time_s: f64,
        u: &[f64],
        v: &[f64],
        psi: &[f64],
        p: &[f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<SparseStructuralTimeState<'_>> {
        self.state_with_stresses(
            time_s,
            u,
            v,
            psi,
            p,
            &vec![0.0; self.branches.len()],
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn state_with_stresses(
        &self,
        time_s: f64,
        u: &[f64],
        v: &[f64],
        psi: &[f64],
        p: &[f64],
        stresses: &[f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<SparseStructuralTimeState<'_>> {
        poll(cancelled)?;
        let ns = self.structural_dofs();
        let na = self.pressure_dofs();
        if !time_s.is_finite()
            || time_s < 0.0
            || u.len() != ns
            || v.len() != ns
            || psi.len() != na
            || p.len() != na
            || u.iter()
                .chain(v)
                .chain(psi)
                .chain(p)
                .any(|x| !x.is_finite())
        {
            return Err(SparseStructuralTimeError::Invalid(
                "complete finite original state and nonnegative clock required",
            ));
        }
        self.branches.energy(stresses, cancelled)?;
        let mut position = Vec::with_capacity(ns + na);
        position.extend_from_slice(u);
        position.extend_from_slice(psi);
        let mut velocity = Vec::with_capacity(ns + na);
        velocity.extend_from_slice(v);
        velocity.extend_from_slice(p);
        let state = SparseStructuralTimeState {
            owner: self,
            time_s,
            position,
            velocity,
            stress: stresses.to_vec(),
        };
        self.energy(&state, cancelled)?;
        Ok(state)
    }
    fn check(&self, state: &SparseStructuralTimeState<'_>) -> Result<()> {
        if !core::ptr::eq(self, state.owner) {
            Err(SparseStructuralTimeError::Invalid(
                "state must belong to this exact prepared spatial/material/time owner",
            ))
        } else {
            Ok(())
        }
    }
    pub(super) fn energy(
        &self,
        state: &SparseStructuralTimeState<'_>,
        cancelled: &impl Fn() -> bool,
    ) -> Result<SparseStructuralTimeEnergy> {
        Ok(self.energy_with_roundoff(state, cancelled)?.0)
    }
    fn energy_with_roundoff(
        &self,
        state: &SparseStructuralTimeState<'_>,
        cancelled: &impl Fn() -> bool,
    ) -> Result<(SparseStructuralTimeEnergy, [f64; 3])> {
        self.check(state)?;
        let (ks, kscale) = positive(&self.forms.mass[0], state.structural_velocity(), cancelled)?;
        let (es, escale) = positive(
            &self.forms.stiffness[0],
            state.structural_displacement(),
            cancelled,
        )?;
        let (ka, ascale) = positive(
            &self.forms.stiffness[1],
            state.pressure_potential_pa_s(),
            cancelled,
        )?;
        let (ca, cscale) = positive(&self.forms.mass[1], state.pressure_pa(), cancelled)?;
        let kinetic_s = 0.5 * ks;
        let elastic_s = 0.5 * es;
        let kinetic_a = 0.5 * ka;
        let compression = 0.5 * ca;
        let result = SparseStructuralTimeEnergy {
            structural_kinetic_j: kinetic_s,
            structural_elastic_j: elastic_s,
            acoustic_kinetic_j: kinetic_a,
            acoustic_compression_j: compression,
            relaxation_j: self.branches.energy(&state.stress, cancelled)?,
        };
        if !result.total_j().is_finite() {
            return Err(SparseStructuralTimeError::Invalid(
                "physical storage overflow",
            ));
        }
        Ok((
            result,
            [
                0.5 * (kscale + escale),
                0.5 * (ascale + cscale),
                result.relaxation_j,
            ],
        ))
    }
    pub(super) fn step<'p>(
        &'p self,
        initial: &SparseStructuralTimeState<'p>,
        force: &[f64],
        source: &[f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<SparseStructuralTimeStep<'p>> {
        poll(cancelled)?;
        self.check(initial)?;
        let ns = self.structural_dofs();
        let na = self.pressure_dofs();
        let n = ns + na;
        if force.len() != ns
            || source.len() != na
            || force.iter().chain(source).any(|x| !x.is_finite())
        {
            return Err(SparseStructuralTimeError::Invalid(
                "finite complete original force and volume-flow loads required",
            ));
        }
        let h = self.duration_s;
        let t = 0.5 * h;
        let next_time = initial.time_s + h;
        if !next_time.is_finite() || next_time <= initial.time_s {
            return Err(SparseStructuralTimeError::Invalid(
                "time increment is unrepresentable",
            ));
        }
        let mut midpoint = vec![0.0; n + self.branches.len()];
        for (d, load) in [force, source].into_iter().enumerate() {
            let (offset, size) = if d == 0 { (0, ns) } else { (ns, na) };
            for row in 0..size {
                poll(cancelled)?;
                let mass = row_product(
                    &self.forms.mass[d],
                    row,
                    &initial.velocity[offset..offset + size],
                )
                .0;
                let spring = row_product(
                    &self.forms.stiffness[d],
                    row,
                    &initial.position[offset..offset + size],
                )
                .0;
                let rhs = (mass + Dd::from_f64(t) * (Dd::from_f64(load[row]) - spring)).to_f64();
                midpoint[offset + row] =
                    rhs * self.scales[offset + row] * if d == 0 { 1.0 } else { -1.0 };
            }
        }
        self.branches
            .rhs(&initial.stress, &mut midpoint[n..], cancelled)?;
        finite(&midpoint)?;
        poll(cancelled)?;
        self.factor.solve_in_place(&mut midpoint);
        poll(cancelled)?;
        for (x, &scale) in midpoint.iter_mut().zip(&self.scales) {
            *x *= scale;
        }
        finite(&midpoint)?;
        let mut next = SparseStructuralTimeState {
            owner: self,
            time_s: next_time,
            position: vec![0.0; n],
            velocity: vec![0.0; n],
            stress: self
                .branches
                .next_stresses(&initial.stress, &midpoint[n..], cancelled)?,
        };
        for i in 0..n {
            next.position[i] = h.mul_add(midpoint[i], initial.position[i]);
            next.velocity[i] = (Dd::from_f64(2.0) * Dd::from_f64(midpoint[i])
                - Dd::from_f64(initial.velocity[i]))
            .to_f64();
        }
        finite(&next.position)?;
        finite(&next.velocity)?;
        // Recover physical midpoint fields from the candidate endpoints. The
        // audits therefore check the published states, not just the solve vector.
        let middle_q = (0..n)
            .map(|i| {
                (Dd::from_f64(initial.position[i]) * Dd::from_f64(0.5)
                    + Dd::from_f64(next.position[i]) * Dd::from_f64(0.5))
                .to_f64()
            })
            .collect::<Vec<_>>();
        for i in 0..n {
            midpoint[i] = (Dd::from_f64(initial.velocity[i]) * Dd::from_f64(0.5)
                + Dd::from_f64(next.velocity[i]) * Dd::from_f64(0.5))
            .to_f64();
        }
        let memory = self.branches.audit(
            &initial.stress,
            &next.stress,
            &midpoint[..n],
            h,
            &self.limits,
            cancelled,
        )?;
        let (fluid_force, displaced_flow) = interface(
            &self.forms.coupling,
            &midpoint[..ns],
            &midpoint[ns..n],
            cancelled,
        )?;
        let residual = self
            .residual(
                initial,
                &next,
                &midpoint,
                &middle_q,
                force,
                source,
                &fluid_force,
                &displaced_flow,
                &memory.reaction,
                cancelled,
            )?
            .max(memory.relative_residual);
        if residual > self.limits.maximum_relative_residual {
            return Err(SparseStructuralTimeError::Residual {
                observed: residual,
                maximum: self.limits.maximum_relative_residual,
            });
        }
        let (initial_energy, initial_roundoff) = self.energy_with_roundoff(initial, cancelled)?;
        let (final_energy, final_roundoff) = self.energy_with_roundoff(&next, cancelled)?;
        let mut changes = [0.0; 2];
        let mut change_scales = [0.0; 2];
        let mut losses = [0.0; 2];
        let mut loss_scales = [0.0; 2];
        for d in 0..2 {
            let (offset, size) = if d == 0 { (0, ns) } else { (ns, na) };
            let mut change = Dd::from_f64(0.0);
            for (a, x0, x1) in [
                (
                    &self.forms.mass[d],
                    &initial.velocity[offset..offset + size],
                    &next.velocity[offset..offset + size],
                ),
                (
                    &self.forms.stiffness[d],
                    &initial.position[offset..offset + size],
                    &next.position[offset..offset + size],
                ),
            ] {
                let (value, scale) = energy_change(a, x0, x1, cancelled)?;
                change = change + Dd::from_f64(value);
                change_scales[d] += scale;
            }
            changes[d] = change.to_f64();
            let (loss, scale) = positive(
                &self.forms.damping[d],
                &midpoint[offset..offset + size],
                cancelled,
            )?;
            losses[d] = h * loss;
            loss_scales[d] = h * scale;
        }
        let (input_s, input_s_scale) = pair(h, force, &midpoint[..ns])?;
        let (input_a, input_a_scale) = pair(h, source, &midpoint[ns..n])?;
        let (exchange_s, exchange_s_scale) = pair(h, &fluid_force, &midpoint[..ns])?;
        let (exchange_a, exchange_a_scale) = pair(h, &displaced_flow, &midpoint[ns..n])?;
        let work = SparseStructuralTimeWork {
            initial_energy,
            final_energy,
            structural_energy_change_j: changes[0],
            acoustic_energy_change_j: changes[1],
            structural_input_work_j: input_s,
            acoustic_input_work_j: input_a,
            fluid_to_structure_work_j: exchange_s,
            acoustic_interface_work_j: exchange_a,
            structural_dissipation_j: losses[0],
            wall_dissipation_j: losses[1],
            relaxation_energy_change_j: memory.energy_change,
            structural_to_relaxation_work_j: memory.work[0],
            acoustic_to_relaxation_work_j: memory.work[1],
            relaxation_dissipation_j: memory.dissipation,
            maximum_relative_residual: residual,
        };
        let scales = [
            change_scales[0]
                + input_s_scale
                + exchange_s_scale
                + loss_scales[0]
                + initial_roundoff[0]
                + final_roundoff[0]
                + memory.work_scale[0],
            change_scales[1]
                + input_a_scale
                + exchange_a_scale
                + loss_scales[1]
                + initial_roundoff[1]
                + final_roundoff[1]
                + memory.work_scale[1],
        ];
        for (account, defect, physical, roundoff) in [
            (
                "structural",
                work.structural_balance_residual_j(),
                changes[0].abs()
                    + input_s.abs()
                    + exchange_s.abs()
                    + losses[0].abs()
                    + memory.work[0].abs(),
                scales[0],
            ),
            (
                "acoustic",
                work.acoustic_balance_residual_j(),
                changes[1].abs()
                    + input_a.abs()
                    + exchange_a.abs()
                    + losses[1].abs()
                    + memory.work[1].abs(),
                scales[1],
            ),
            (
                "material memory",
                work.relaxation_balance_residual_j(),
                memory.energy_change.abs()
                    + memory.dissipation
                    + memory.work[0].abs()
                    + memory.work[1].abs(),
                memory.energy_scale + initial_roundoff[2] + final_roundoff[2],
            ),
            (
                "interface",
                work.interface_balance_residual_j(),
                exchange_s.abs() + exchange_a.abs(),
                exchange_s_scale + exchange_a_scale,
            ),
            (
                "total",
                work.balance_residual_j(),
                changes[0].abs()
                    + changes[1].abs()
                    + input_s.abs()
                    + input_a.abs()
                    + losses[0].abs()
                    + losses[1].abs()
                    + memory.energy_change.abs()
                    + memory.dissipation,
                scales[0]
                    + scales[1]
                    + memory.energy_scale
                    + initial_roundoff[2]
                    + final_roundoff[2],
            ),
        ] {
            let maximum = self.limits.energy_balance_tolerance * physical + ROUND * roundoff;
            if !maximum.is_finite() || !defect.is_finite() {
                return Err(SparseStructuralTimeError::Invalid(
                    "physical work account overflow",
                ));
            }
            if defect.abs() > maximum {
                return Err(SparseStructuralTimeError::EnergyBalance {
                    account,
                    residual_j: defect,
                    maximum_j: maximum,
                });
            }
        }
        poll(cancelled)?;
        Ok(SparseStructuralTimeStep {
            state: next,
            work,
            fluid_force,
            displaced_flow_m3_s: displaced_flow,
            memory: memory.memory,
        })
    }
    #[allow(clippy::too_many_arguments)]
    fn residual(
        &self,
        initial: &SparseStructuralTimeState<'_>,
        next: &SparseStructuralTimeState<'_>,
        middle_v: &[f64],
        middle_q: &[f64],
        force: &[f64],
        source: &[f64],
        fluid_force: &[f64],
        displaced_flow: &[f64],
        memory_reaction: &[f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<f64> {
        let ns = self.structural_dofs();
        let na = self.pressure_dofs();
        let h = Dd::from_f64(self.duration_s);
        let mut maximum = 0.0_f64;
        for (d, load) in [force, source].into_iter().enumerate() {
            let (offset, size) = if d == 0 { (0, ns) } else { (ns, na) };
            for row in 0..size {
                poll(cancelled)?;
                let mut momentum = Dd::from_f64(0.0);
                let mut momentum_scale = 0.0;
                let (cols, vals) = self.forms.mass[d].row(row);
                for (&col, &value) in cols.iter().zip(vals) {
                    let term = Dd::from_f64(value)
                        * (Dd::from_f64(next.velocity[offset + col])
                            - Dd::from_f64(initial.velocity[offset + col]));
                    momentum = momentum + term;
                    momentum_scale += value.abs()
                        * (next.velocity[offset + col].abs()
                            + initial.velocity[offset + col].abs());
                }
                let (drag, ds) = row_product(
                    &self.forms.damping[d],
                    row,
                    &middle_v[offset..offset + size],
                );
                let (spring, ks) = row_product(
                    &self.forms.stiffness[d],
                    row,
                    &middle_q[offset..offset + size],
                );
                let exchange = if d == 0 {
                    -fluid_force[row]
                } else {
                    displaced_flow[row]
                };
                let memory = memory_reaction.get(offset + row).copied().unwrap_or(0.0);
                let defect = (momentum
                    + h * (drag + spring + Dd::from_f64(exchange) - Dd::from_f64(load[row])
                        + Dd::from_f64(memory)))
                .to_f64();
                let scale = momentum_scale
                    + self.duration_s * (ds + ks + exchange.abs() + load[row].abs() + memory.abs());
                maximum = maximum.max(relative(defect, scale)?);
                let i = offset + row;
                let defect = (Dd::from_f64(next.position[i])
                    - Dd::from_f64(initial.position[i])
                    - h * Dd::from_f64(middle_v[i]))
                .to_f64();
                let scale = next.position[i].abs()
                    + initial.position[i].abs()
                    + self.duration_s * middle_v[i].abs();
                maximum = maximum.max(relative(defect, scale)?);
            }
        }
        Ok(maximum)
    }
}

fn push(
    a: &mut Coo,
    row: usize,
    col: usize,
    value: f64,
    limits: &SparseStructuralTimeLimits,
) -> Result<()> {
    if !value.is_finite() {
        return Err(SparseStructuralTimeError::Invalid(
            "sparse operator assembly overflow",
        ));
    }
    if value != 0.0 {
        bound(
            "staged matrix entries",
            a.staged()
                .checked_add(1)
                .ok_or(SparseStructuralTimeError::Invalid(
                    "sparse staging extent overflow",
                ))?,
            limits.maximum_matrix_entries,
        )?;
        a.push(row, col, value);
    }
    Ok(())
}
fn finite(x: &[f64]) -> Result<()> {
    if x.iter().any(|v| !v.is_finite()) {
        Err(SparseStructuralTimeError::Invalid(
            "nonfinite sparse state or action",
        ))
    } else {
        Ok(())
    }
}
fn relative(defect: f64, scale: f64) -> Result<f64> {
    if !defect.is_finite() || !scale.is_finite() {
        Err(SparseStructuralTimeError::Invalid(
            "physical equation audit overflow",
        ))
    } else if scale == 0.0 {
        if defect == 0.0 {
            Ok(0.0)
        } else {
            Err(SparseStructuralTimeError::Invalid(
                "unscaled physical equation defect",
            ))
        }
    } else {
        Ok(defect.abs() / scale)
    }
}
fn factor(
    a: &Csr,
    limits: &SparseStructuralTimeLimits,
    cancelled: &impl Fn() -> bool,
) -> Result<LdltFactor> {
    let symbolic =
        SymbolicLdlt::analyze_bounded(a, DirectOrdering::Amd, &limits.factor, cancelled)?;
    let result = symbolic.factor_bounded(
        a,
        &LdltOptions {
            max_factor_bytes: Some(limits.maximum_factor_bytes),
            ..LdltOptions::default()
        },
        &limits.factor,
        cancelled,
    )?;
    poll(cancelled)?;
    Ok(result)
}
fn scale_form(a: &mut Csr, scales: &[f64], cancelled: &impl Fn() -> bool) -> Result<()> {
    // Each triangular pair uses the SAME multiplication order, retaining exact
    // numerical symmetry after congruence even across disparate physical units.
    for row in 0..a.nrows() {
        poll(cancelled)?;
        let cols = a.row(row).0.to_vec();
        for (&col, value) in cols.iter().zip(a.row_values_mut(row)) {
            let (i, j) = if row <= col { (row, col) } else { (col, row) };
            *value = (*value * scales[i]) * scales[j];
            if !value.is_finite() {
                return Err(SparseStructuralTimeError::Invalid(
                    "sparse congruence overflow",
                ));
            }
        }
    }
    Ok(())
}
fn admit_form(
    a: &Csr,
    strict: bool,
    limits: &SparseStructuralTimeLimits,
    cancelled: &impl Fn() -> bool,
) -> Result<()> {
    let n = a.nrows();
    let mut scales = Vec::with_capacity(n);
    let mut any = false;
    for row in 0..n {
        poll(cancelled)?;
        let diagonal = a.get(row, row);
        if !diagonal.is_finite() || diagonal < 0.0 || (strict && diagonal == 0.0) {
            return Err(SparseStructuralTimeError::Invalid(
                "original mass must be positive and stiffness/damping nonnegative on every diagonal",
            ));
        }
        let (cols, vals) = a.row(row);
        for (&col, &value) in cols.iter().zip(vals) {
            // No symmetrization: gyroscopic/nonreciprocal structures use the
            // existing exact owner. Canonical native FE forms are symmetric.
            if !value.is_finite() || value != a.get(col, row) || (diagonal == 0.0 && value != 0.0) {
                return Err(SparseStructuralTimeError::Invalid(
                    "symmetric finite passive original forms required; no coupled zero diagonal",
                ));
            }
            any |= value != 0.0;
        }
        let scale = if diagonal > 0.0 {
            1.0 / det::sqrt(diagonal)
        } else {
            1.0
        };
        if !scale.is_finite() || scale <= 0.0 {
            return Err(SparseStructuralTimeError::Invalid(
                "source form congruence is unrepresentable",
            ));
        }
        scales.push(scale);
    }
    if !any {
        return Ok(());
    }
    let mut admitted = a.clone();
    scale_form(&mut admitted, &scales, cancelled)?;
    if !strict {
        // Inertia of A_scaled + ROUND I admits numerical PSD, including rigid
        // and constant-potential null spaces. This is ONLY an admission probe:
        // original A, every step, every energy and every heat remain unchanged.
        let mut shifted = Coo::new(n, n);
        for row in 0..n {
            poll(cancelled)?;
            let (cols, vals) = admitted.row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                push(&mut shifted, row, col, value, limits)?;
            }
            push(&mut shifted, row, row, ROUND, limits)?;
        }
        admitted = shifted.assemble();
    }
    let f = factor(&admitted, limits, cancelled)?;
    if f.inertia().positive != n || f.inertia().negative != 0 {
        return Err(SparseStructuralTimeError::Invalid(
            "original storage/dissipation has a resolved negative direction",
        ));
    }
    Ok(())
}
fn row_product(a: &Csr, row: usize, x: &[f64]) -> (Dd, f64) {
    let mut sum = Dd::from_f64(0.0);
    let mut scale = 0.0;
    let (cols, vals) = a.row(row);
    for (&col, &value) in cols.iter().zip(vals) {
        let term = Dd::from_f64(value) * Dd::from_f64(x[col]);
        sum = sum + term;
        scale += term.to_f64().abs();
    }
    (sum, scale)
}
fn positive(a: &Csr, x: &[f64], cancelled: &impl Fn() -> bool) -> Result<(f64, f64)> {
    let mut sum = Dd::from_f64(0.0);
    let mut scale = 0.0;
    for row in 0..a.nrows() {
        poll(cancelled)?;
        let (cols, vals) = a.row(row);
        for (&col, &value) in cols.iter().zip(vals) {
            let term = Dd::from_f64(x[row]) * Dd::from_f64(value) * Dd::from_f64(x[col]);
            sum = sum + term;
            scale += term.to_f64().abs();
        }
    }
    let value = sum.to_f64();
    if !value.is_finite() || !scale.is_finite() || value < -ROUND * scale {
        return Err(SparseStructuralTimeError::Invalid(
            "negative or nonfinite original physical storage/dissipation",
        ));
    }
    Ok((value, scale))
}
fn energy_change(
    a: &Csr,
    x0: &[f64],
    x1: &[f64],
    cancelled: &impl Fn() -> bool,
) -> Result<(f64, f64)> {
    let mut sum = Dd::from_f64(0.0);
    let mut scale = 0.0;
    for row in 0..a.nrows() {
        poll(cancelled)?;
        let difference = Dd::from_f64(x1[row]) - Dd::from_f64(x0[row]);
        let (cols, vals) = a.row(row);
        for (&col, &value) in cols.iter().zip(vals) {
            let term = Dd::from_f64(0.5)
                * difference
                * Dd::from_f64(value)
                * (Dd::from_f64(x1[col]) + Dd::from_f64(x0[col]));
            sum = sum + term;
            scale += term.to_f64().abs();
        }
    }
    let value = sum.to_f64();
    if !value.is_finite() || !scale.is_finite() {
        return Err(SparseStructuralTimeError::Invalid(
            "storage change overflow",
        ));
    }
    Ok((value, scale))
}
fn pair(h: f64, x: &[f64], y: &[f64]) -> Result<(f64, f64)> {
    let mut sum = Dd::from_f64(0.0);
    let mut scale = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let term = Dd::from_f64(h) * Dd::from_f64(a) * Dd::from_f64(b);
        sum = sum + term;
        scale += term.to_f64().abs();
    }
    let value = sum.to_f64();
    if !value.is_finite() || !scale.is_finite() {
        return Err(SparseStructuralTimeError::Invalid("physical work overflow"));
    }
    Ok((value, scale))
}
fn interface(
    b: &Csr,
    v: &[f64],
    p: &[f64],
    cancelled: &impl Fn() -> bool,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut force = vec![Dd::from_f64(0.0); b.ncols()];
    let mut flow = Vec::with_capacity(b.nrows());
    for row in 0..b.nrows() {
        poll(cancelled)?;
        let (cols, vals) = b.row(row);
        let mut sum = Dd::from_f64(0.0);
        for (&col, &value) in cols.iter().zip(vals) {
            sum = sum + Dd::from_f64(value) * Dd::from_f64(v[col]);
            force[col] = force[col] + Dd::from_f64(value) * Dd::from_f64(p[row]);
        }
        flow.push(sum.to_f64());
    }
    let force = force.into_iter().map(Dd::to_f64).collect::<Vec<_>>();
    finite(&force)?;
    finite(&flow)?;
    Ok((force, flow))
}

#[cfg(test)]
mod tests;
