//! Sparse complete structural/interior waves beyond the dense exponential lane.
//!
//! The original structural CSR mass, relaxed elastic stiffness and symmetric
//! memoryless damping couple reciprocally to the original heterogeneous
//! pressure-potential/pressure FE forms. The native normal trace gives B once;
//! BOTH signs of its work pairing remain in one simultaneous sparse midpoint
//! solve. No pressure mode, structural rigid mode, lateral elastic wave,
//! off-diagonal inertia or constant cavity compression is removed.
//!
//! This is a second-order passive time discretization of the declared frozen,
//! linear model. Temporal phase needs refinement independently of spatial and
//! material quadrature. The existing exact owner remains the appropriate small
//! system route for exact held-input evolution and gyroscopic structures.
//! Original Maxwell channels optionally retain structural material or acoustic
//! wall memory through their own stress states and positive heat. Relaxed K
//! excludes these branch stiffnesses. Frequency-local damping samples do not
//! authorize a time law. Mixed structural/acoustic channels, gyroscopic forms,
//! changing material laws, thermal diffusion, exterior radiation and pressure
//! release require their respective original dynamic owners.
//!
//! Preparation copies only sparse physical forms and caches one native bounded
//! LDLT factor. Steps allocate O(n) candidate/audit vectors and perform sparse
//! triangular solves and original-form contractions; no factorization, dense
//! mass inverse, material query or coordinate reduction occurs during a step.
//! The native solve is synchronous; cancellation is checked around it and
//! between sparse assembly/audit rows. A refusal leaves accepted state intact.

use super::{
    ConformingTetStructure, NormalDisplacementTrace, StructuralInteriorError, StructuralMatrices,
    conforming_tet_traces,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    HeterogeneousHelmholtzFe,
    time::{RealVelocityTrace, ResistiveBoundary},
};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

mod physics;
/// Original source columns on the existing complete sparse midpoint timeline.
pub mod sources;
pub use fs_modal::relaxation_modes::MaxwellDampingBranch;
use physics::{BranchSource, SparseForms, SparseMidpoint, bound, poll};
pub use physics::{
    SparseStructuralMemoryWork, SparseStructuralRelaxationLimits, SparseStructuralTimeEnergy,
    SparseStructuralTimeError, SparseStructuralTimeLimits, SparseStructuralTimeState,
    SparseStructuralTimeStep, SparseStructuralTimeWork,
};

/// Geometry/trace errors preserve the original owner; time refusals retain
/// native sparse pivot, memory, work and physical-account details.
#[derive(Debug)]
pub enum SparseStructuralInteriorError {
    Interior(StructuralInteriorError),
    Time(SparseStructuralTimeError),
}
impl core::fmt::Display for SparseStructuralInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Interior(e) => e.fmt(f),
            Self::Time(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for SparseStructuralInteriorError {}
impl From<StructuralInteriorError> for SparseStructuralInteriorError {
    fn from(e: StructuralInteriorError) -> Self {
        Self::Interior(e)
    }
}
impl From<SparseStructuralTimeError> for SparseStructuralInteriorError {
    fn from(e: SparseStructuralTimeError) -> Self {
        Self::Time(e)
    }
}

/// One original geometry/material source and its complete reciprocal time factor.
pub struct PreparedSparseStructuralInterior<'a, 'f, 'c> {
    interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
    evolution: SparseMidpoint,
}
impl<'a, 'f, 'c> PreparedSparseStructuralInterior<'a, 'f, 'c> {
    /// Prepare the ORIGINAL complete structural and acoustic coordinate spaces.
    ///
    /// Matrices must be exactly symmetric; no gyroscopic part is discarded.
    /// Native scaled inertia admits SPD mass and PSD stiffness/damping within
    /// 4096*epsilon of the diagonal-congruence coordinates. The PSD shift is
    /// used ONLY in admission; every stored coefficient, evolution, storage
    /// and heat contraction uses the unchanged original form. Resolved negative
    /// directions and a zero diagonal with nonzero coupling refuse.
    ///
    /// Normal displacement is positive OUT of the fluid. Each structural
    /// coordinate supplies its original trace; fixed/uncoupled coordinates may
    /// use an empty collection. Real wall admittance is relative to structural
    /// motion and contributes only its own independently integrated heat.
    /// Physical basis identity, units and constitutive validity remain authored.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        structure: StructuralMatrices<'_>,
        normal_traces: &[Vec<NormalDisplacementTrace>],
        walls: &[ResistiveBoundary],
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, SparseStructuralInteriorError> {
        Self::new_relaxing(
            interior,
            structure,
            normal_traces,
            walls,
            &[],
            duration_s,
            limits,
            SparseStructuralRelaxationLimits::default(),
            cx,
        )
    }

    /// Keep original Maxwell material/wall histories in the SAME full sparse
    /// solve. K must be RELAXED: branch stiffness is not also included in K.
    ///
    /// Each original strain follows complete [u,psi] coordinate order. A
    /// structural-only map supplies hereditary solid stress. A potential-only
    /// map gives the causal wall admittance k*tau/(1+s*tau)*b*b^T. Empty maps
    /// retain their original unforced stress relaxation. Mixed cross-domain
    /// maps refuse intact; the original general exact owner supports them.
    ///
    /// The cached sparse factor adds one scalar per branch with its original
    /// local sparse direction; it never forms dense b*b^T updates. Original
    /// storage, reaction work and sigma^2/(k*tau) heat remain separate. Implicit
    /// midpoint needs temporal refinement to resolve short relaxation times;
    /// it is passive but is neither exact nor an L-stable history integrator.
    #[allow(clippy::too_many_arguments)]
    pub fn new_relaxing(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        structure: StructuralMatrices<'_>,
        normal_traces: &[Vec<NormalDisplacementTrace>],
        walls: &[ResistiveBoundary],
        branches: &[MaxwellDampingBranch],
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        memory_limits: SparseStructuralRelaxationLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, SparseStructuralInteriorError> {
        let cancelled = || cx.checkpoint().is_err();
        poll(&cancelled)?;
        let ns = structure.mass.nrows();
        let na = interior.ndof();
        let n = limits.admit(ns, na)?;
        bound(
            "material branches",
            branches.len(),
            memory_limits.maximum_branches,
        )?;
        bound(
            "material source coefficients",
            n.checked_mul(branches.len())
                .ok_or(SparseStructuralTimeError::Invalid(
                    "material source extent overflow",
                ))?,
            memory_limits.maximum_source_entries,
        )?;
        if normal_traces.len() != ns {
            return Err(SparseStructuralTimeError::Invalid(
                "one trace collection per complete structural coordinate required",
            )
            .into());
        }
        let mut trace_terms = 0usize;
        for count in normal_traces
            .iter()
            .map(Vec::len)
            .chain(walls.iter().map(|w| w.faces.len()))
        {
            trace_terms = trace_terms
                .checked_add(count)
                .ok_or(SparseStructuralTimeError::Invalid("trace extent overflow"))?;
            bound("trace terms", trace_terms, limits.maximum_trace_terms)?;
        }
        let mut entries = 0usize;
        for (matrix, size) in [
            (structure.mass, ns),
            (structure.stiffness, ns),
            (interior.mass(), na),
            (interior.stiffness(), na),
        ]
        .into_iter()
        .chain(structure.damping.into_iter().map(|a| (a, ns)))
        {
            if matrix.nrows() != size || matrix.ncols() != size {
                return Err(SparseStructuralTimeError::Invalid(
                    "complete original matrix shapes required",
                )
                .into());
            }
            entries =
                entries
                    .checked_add(matrix.nnz())
                    .ok_or(SparseStructuralTimeError::Invalid(
                        "source matrix extent overflow",
                    ))?;
            bound(
                "source matrix entries",
                entries,
                limits.maximum_matrix_entries,
            )?;
        }
        let mut coupling = Coo::new(na, ns);
        for (s, traces) in normal_traces.iter().enumerate() {
            poll(&cancelled)?;
            if traces.is_empty() {
                continue;
            }
            let load = interior
                .geometry()
                .boundary_linear_load(traces)
                .map_err(StructuralInteriorError::from)?;
            for (a, value) in load.into_iter().enumerate() {
                if !value.is_finite() {
                    return Err(SparseStructuralTimeError::Invalid(
                        "native normal-trace load overflow",
                    )
                    .into());
                }
                if value != 0.0 {
                    bound(
                        "coupling entries",
                        coupling.staged() + 1,
                        limits.maximum_matrix_entries,
                    )?;
                    coupling.push(a, s, value);
                }
            }
        }
        let mut conductance = Coo::new(na, na);
        let mut seen = BTreeSet::new();
        for wall in walls {
            poll(&cancelled)?;
            if !wall.admittance_m_pa_s.is_finite() || wall.admittance_m_pa_s < 0.0 {
                return Err(SparseStructuralTimeError::Invalid(
                    "finite nonnegative relative-wall admittance required",
                )
                .into());
            }
            for &face in &wall.faces {
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(SparseStructuralTimeError::Invalid(
                        "relative-wall faces must be unique",
                    )
                    .into());
                }
            }
            let surface = interior
                .geometry()
                .boundary_mass(&wall.faces)
                .map_err(StructuralInteriorError::from)?;
            for row in 0..na {
                poll(&cancelled)?;
                let (cols, vals) = surface.row(row);
                for (&col, &value) in cols.iter().zip(vals) {
                    let entry = value * wall.admittance_m_pa_s;
                    if !entry.is_finite() {
                        return Err(SparseStructuralTimeError::Invalid(
                            "relative-wall conductance overflow",
                        )
                        .into());
                    }
                    if entry != 0.0 {
                        bound(
                            "wall matrix entries",
                            conductance.staged() + 1,
                            limits.maximum_matrix_entries,
                        )?;
                        conductance.push(row, col, entry);
                    }
                }
            }
        }
        let coupling = coupling.assemble();
        let conductance = conductance.assemble();
        entries = entries
            .checked_add(coupling.nnz())
            .and_then(|v| v.checked_add(conductance.nnz()))
            .ok_or(SparseStructuralTimeError::Invalid(
                "complete sparse form extent overflow",
            ))?;
        bound(
            "retained matrix entries",
            entries,
            limits.maximum_matrix_entries,
        )?;
        let branch_sources = branches
            .iter()
            .map(|b| BranchSource {
                strain: &b.strain,
                stiffness: b.stiffness,
                relaxation_time_s: b.relaxation_time_s,
            })
            .collect::<Vec<_>>();
        let forms = SparseForms {
            mass: [structure.mass.clone(), interior.mass().clone()],
            stiffness: [structure.stiffness.clone(), interior.stiffness().clone()],
            damping: [
                structure
                    .damping
                    .cloned()
                    .unwrap_or_else(|| Coo::new(ns, ns).assemble()),
                conductance,
            ],
            coupling,
        };
        let evolution = if branch_sources.is_empty() {
            SparseMidpoint::new(forms, duration_s, limits, &cancelled)?
        } else {
            SparseMidpoint::new_relaxing(
                forms,
                duration_s,
                limits,
                &branch_sources,
                memory_limits,
                &cancelled,
            )?
        };
        Ok(Self {
            interior,
            evolution,
        })
    }

    /// Use the SAME conforming tetrahedral interface adapter as the complete
    /// harmonic/exact owner: original free DOFs, exact vertex coincidence,
    /// owning-volume outward normal and opposite solid/fluid sidedness.
    /// No dense coupled model is constructed on the way to this sparse factor.
    #[allow(clippy::too_many_arguments)]
    pub fn from_tet_elastic(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        solid: ConformingTetStructure<'_>,
        damping: Option<&Csr>,
        walls: &[ResistiveBoundary],
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, SparseStructuralInteriorError> {
        Self::from_tet_elastic_relaxing(
            interior,
            solid,
            damping,
            walls,
            &[],
            duration_s,
            limits,
            SparseStructuralRelaxationLimits::default(),
            cx,
        )
    }

    /// Original conforming free-tetrahedron trace plus authored Maxwell maps
    /// on its complete [free structural DOFs, pressure potential] space.
    #[allow(clippy::too_many_arguments)]
    pub fn from_tet_elastic_relaxing(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        solid: ConformingTetStructure<'_>,
        damping: Option<&Csr>,
        walls: &[ResistiveBoundary],
        branches: &[MaxwellDampingBranch],
        duration_s: f64,
        limits: SparseStructuralTimeLimits,
        memory_limits: SparseStructuralRelaxationLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, SparseStructuralInteriorError> {
        limits.admit(solid.assembly.free_dofs.len(), interior.ndof())?;
        let maximum_terms =
            solid
                .interface
                .len()
                .checked_mul(9)
                .ok_or(SparseStructuralTimeError::Invalid(
                    "conforming trace extent overflow",
                ))?;
        bound(
            "conforming trace terms",
            maximum_terms,
            limits.maximum_trace_terms,
        )?;
        let traces = conforming_tet_traces(interior, &solid, cx)?;
        Self::new_relaxing(
            interior,
            StructuralMatrices {
                mass: &solid.assembly.mass,
                stiffness: &solid.assembly.stiffness,
                damping,
            },
            &traces,
            walls,
            branches,
            duration_s,
            limits,
            memory_limits,
            cx,
        )
    }
    /// Every supplied original scalar history, including dark branches.
    #[must_use]
    pub fn branch_count(&self) -> usize {
        self.evolution.branch_count()
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.evolution.duration_s()
    }
    #[must_use]
    pub fn structural_dofs(&self) -> usize {
        self.evolution.structural_dofs()
    }
    #[must_use]
    pub fn pressure_dofs(&self) -> usize {
        self.evolution.pressure_dofs()
    }
    /// Original sparse reciprocal map: acoustic rows, structural columns.
    #[must_use]
    pub fn coupling(&self) -> &Csr {
        self.evolution.coupling()
    }

    /// Original native full pressure-trace projection, with no extra area,
    /// density or angular frequency. Supply its physical midpoint values.
    pub fn project_inward_velocity(
        &self,
        traces: &[RealVelocityTrace],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SparseStructuralInteriorError> {
        poll(&|| cx.checkpoint().is_err())?;
        let load = self
            .interior
            .geometry()
            .boundary_linear_load(traces)
            .map_err(StructuralInteriorError::from)?;
        poll(&|| cx.checkpoint().is_err())?;
        Ok(load)
    }
    #[must_use]
    pub fn zero_state(&self) -> SparseStructuralTimeState<'_> {
        self.evolution.zero_state()
    }
    /// Import complete original coordinates with explicitly zero initial stress
    /// in every admitted Maxwell branch. Nothing is projected or removed.
    #[allow(clippy::too_many_arguments)]
    pub fn state(
        &self,
        time_s: f64,
        u: &[f64],
        v: &[f64],
        psi: &[f64],
        p: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SparseStructuralTimeState<'_>, SparseStructuralTimeError> {
        self.evolution
            .state(time_s, u, v, psi, p, &|| cx.checkpoint().is_err())
    }
    /// Import original nonzero stresses in authored order. They are independent
    /// initial histories, not inferred from displacement or assumed relaxed.
    #[allow(clippy::too_many_arguments)]
    pub fn state_with_stresses(
        &self,
        time_s: f64,
        u: &[f64],
        v: &[f64],
        psi: &[f64],
        p: &[f64],
        stresses: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SparseStructuralTimeState<'_>, SparseStructuralTimeError> {
        self.evolution
            .state_with_stresses(time_s, u, v, psi, p, stresses, &|| cx.checkpoint().is_err())
    }
    pub fn energy(
        &self,
        state: &SparseStructuralTimeState<'_>,
        cx: &Cx<'_>,
    ) -> Result<SparseStructuralTimeEnergy, SparseStructuralTimeError> {
        self.evolution.energy(state, &|| cx.checkpoint().is_err())
    }
    /// Midpoint structural forces and original projected inward volume flow.
    /// Each return includes independent structural, acoustic, interface,
    /// material-memory and total work residuals; every original branch ODE
    /// and branch energy balance is checked before publishing the candidate.
    pub fn step<'p>(
        &'p self,
        initial: &SparseStructuralTimeState<'p>,
        structural_force: &[f64],
        inward_flow_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SparseStructuralTimeStep<'p>, SparseStructuralTimeError> {
        self.evolution
            .step(initial, structural_force, inward_flow_load, &|| {
                cx.checkpoint().is_err()
            })
    }
}

#[cfg(test)]
mod tests;
