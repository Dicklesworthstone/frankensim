//! Source-bound mean and twice-frequency induction forces on a linear solid.
//!
//! The original P1 Lorentz loads feed the existing complete-coordinate
//! `ExactHarmonicSolid` at [0, 2*omega]. Both the static component and the
//! positive-imaginary sine component retain their original EM phase. Native
//! consistent mass, stiffness, optional physical viscosity, initial transients,
//! resonant growth and independent force-work/heat integration are unchanged.
//!
//! This is prescribed, one-way forcing on fixed geometry. Structural velocity
//! is not fed back into the electrical solve, and mechanical source work is
//! not claimed as a new closed EM debit. No magnetostriction, impressed-winding
//! force, magnetic interface traction, geometry update or acoustic radiation is
//! inferred. Conduction-current force must be the explicitly chosen load scope.

use fs_em::spatial::vector_potential::induction::{
    harmonic::{SpatialHarmonicInductionResponse, loads::SpatialHarmonicLorentzLoads},
    mechanical::{SpatialLorentzLoadError, SpatialLorentzLoadLimits},
};
use fs_exec::Cx;
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_solid::{
    dynamics::{
        ExactHarmonicSolid, ExactLinearSolidControls, ExactLinearSolidError, ExactLinearSolidStep,
    },
    linear3::{TetElasticError, TetLinearElasticProblem},
};
use fs_sparse::Csr;

/// Exact conforming correspondence; no nearest-node interpolation is inferred.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInductionSolidMap<'a> {
    /// Original EM vertex for each original mechanical node, in node order.
    pub field_vertices: &'a [usize],
    /// Original EM tetrahedron for each original mechanical tetrahedron.
    /// Body cells may include nonconductors; only conducting cells exert this load.
    pub field_cells: &'a [usize],
}

/// Native owners retain their own geometry, flow, work and heat allowances.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInductionVibrationControls {
    pub loads: SpatialLorentzLoadLimits,
    pub structure: ExactLinearSolidControls,
    /// Extra load/scatter, correspondence, phase and initial-state storage,
    /// including cold scratch, in real scalar slots. No implicit allocation cap.
    pub maximum_transfer_scalars: usize,
}

#[derive(Debug)]
pub enum SpatialInductionVibrationError {
    Invalid(&'static str),
    Load(SpatialLorentzLoadError),
    Assembly(TetElasticError),
    Structure(ExactLinearSolidError),
    Cancelled,
}
impl core::fmt::Display for SpatialInductionVibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial induction vibration refused: {self:?}")
    }
}
impl std::error::Error for SpatialInductionVibrationError {}
impl From<SpatialLorentzLoadError> for SpatialInductionVibrationError {
    fn from(error: SpatialLorentzLoadError) -> Self {
        match error {
            SpatialLorentzLoadError::Cancelled => Self::Cancelled,
            other => Self::Load(other),
        }
    }
}
impl From<ExactLinearSolidError> for SpatialInductionVibrationError {
    fn from(error: ExactLinearSolidError) -> Self {
        match error {
            ExactLinearSolidError::Cancelled => Self::Cancelled,
            other => Self::Structure(other),
        }
    }
}
impl From<TetElasticError> for SpatialInductionVibrationError {
    fn from(error: TetElasticError) -> Self {
        match error {
            TetElasticError::Cancelled => Self::Cancelled,
            other => Self::Assembly(other),
        }
    }
}
type Result<T> = core::result::Result<T, SpatialInductionVibrationError>;
fn invalid(message: &'static str) -> SpatialInductionVibrationError {
    SpatialInductionVibrationError::Invalid(message)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| SpatialInductionVibrationError::Cancelled)
}

/// Native complete-solid motion and the unchanged field that supplies its load.
///
/// The existing two-tone exponential/independent-work owner admits AT MOST 32
/// free mechanical coordinates. This is a bounded small-body response, not a
/// sparse industrial-mesh or reduced-modal scalability claim. No modes are
/// removed to meet that limit. Supports are homogeneous and stationary; loads
/// at supported nodes remain visible in `loads()` but do zero mechanical work.
pub struct SpatialInductionVibration<
    'a,
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: &'a SpatialHarmonicInductionResponse<'s, 'm, C>,
    problem: &'a TetLinearElasticProblem<'a>,
    field_vertices: Vec<usize>,
    field_cells: Vec<usize>,
    loads: SpatialHarmonicLorentzLoads,
    solid: ExactHarmonicSolid,
    cosine_force_n: Vec<f64>,
    sine_force_n: Vec<f64>,
    next_cosine_force_n: Vec<f64>,
    next_sine_force_n: Vec<f64>,
}

impl<'a, 's, 'm, C: SpatialConductorMaterial> SpatialInductionVibration<'a, 's, 'm, C> {
    /// Bind exact geometry and original material/constraint inputs before any
    /// force transfer. Every mechanical node and tetrahedron names a distinct
    /// original EM entity, with identical global coordinates and connectivity;
    /// local tetrahedron/node permutations are allowed. Mechanical cells can
    /// include a nonconducting support, but only mapped original conducting
    /// cells contribute Lorentz force. No force from another body is included.
    ///
    /// Initial displacement/velocity use the complete mechanical node order.
    /// Fixed coordinates must be exactly zero. Optional `damping_ns_m` follows
    /// the original assembly's increasing free-DOF order; None explicitly means
    /// zero viscosity. The native owner admits its full symmetric PSD law.
    /// Material identities are retained by the borrowed original problem.
    ///
    /// Source phase starts at the EM phasor's t=0. Under exp(-i omega t), the
    /// 2*omega force uses Re(F2)*cos + Im(F2)*sin, with no RMS rescaling.
    /// Preparations exceeding the native 32-free-coordinate two-tone limit
    /// refuse before assembly; no hidden modal truncation is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        field: &'a SpatialHarmonicInductionResponse<'s, 'm, C>,
        problem: &'a TetLinearElasticProblem<'a>,
        correspondence: SpatialInductionSolidMap<'_>,
        damping_ns_m: Option<&Csr>,
        dt_s: f64,
        initial_displacement_m: &[[f64; 3]],
        initial_velocity_m_s: &[[f64; 3]],
        controls: SpatialInductionVibrationControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let original = field.problem();
        let nodes = problem.nodes_m.len();
        let cells = problem.tetrahedra.len();
        let total_dofs = nodes
            .checked_mul(3)
            .ok_or(invalid("solid coordinate count overflow"))?;
        let payload = original
            .positions_m()
            .len()
            .checked_mul(12)
            .and_then(|n| cells.checked_mul(32).and_then(|v| n.checked_add(v)))
            .and_then(|n| nodes.checked_mul(64).and_then(|v| n.checked_add(v)));
        if nodes == 0
            || cells == 0
            || correspondence.field_vertices.len() != nodes
            || correspondence.field_cells.len() != cells
            || initial_displacement_m.len() != nodes
            || initial_velocity_m_s.len() != nodes
            || nodes > problem.budget.maximum_nodes
            || cells > problem.budget.maximum_tetrahedra
            || original.positions_m().len() > controls.loads.max_vertices
            || cells > controls.loads.max_cells
            || payload.is_none_or(|n| n > controls.maximum_transfer_scalars)
            || problem.fixed_dofs.len() > total_dofs
        {
            return Err(invalid(
                "complete bounded original solid/field correspondence required",
            ));
        }
        let mut fixed = problem.fixed_dofs.to_vec();
        fixed.sort_unstable();
        if fixed.last().is_some_and(|&dof| dof >= total_dofs)
            || fixed.windows(2).any(|pair| pair[0] == pair[1])
        {
            return Err(invalid(
                "solid supports must name distinct original coordinates",
            ));
        }
        let free = total_dofs - fixed.len();
        if free == 0 || free > 32 || free > controls.structure.maximum_free_dofs {
            return Err(invalid(
                "two-tone native solid requires 1 through 32 admitted free coordinates",
            ));
        }
        let mut vertices = correspondence.field_vertices.to_vec();
        vertices.sort_unstable();
        if vertices.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid(
                "mechanical nodes must name distinct field vertices",
            ));
        }
        for (node, &vertex) in correspondence.field_vertices.iter().enumerate() {
            poll(cx)?;
            if original.positions_m().get(vertex) != Some(&problem.nodes_m[node])
                || initial_displacement_m[node]
                    .iter()
                    .chain(&initial_velocity_m_s[node])
                    .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "solid node coordinates or initial state differ from admitted field chart",
                ));
            }
        }
        for &dof in &fixed {
            if initial_displacement_m[dof / 3][dof % 3] != 0.0
                || initial_velocity_m_s[dof / 3][dof % 3] != 0.0
            {
                return Err(invalid(
                    "fixed solid coordinates require zero initial motion",
                ));
            }
        }
        let mut field_cells = correspondence.field_cells.to_vec();
        field_cells.sort_unstable();
        if field_cells.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid(
                "mechanical tetrahedra must name distinct field cells",
            ));
        }
        let mut conducting = Vec::new();
        for (element, &cell) in correspondence.field_cells.iter().enumerate() {
            poll(cx)?;
            let source = original
                .mesh()
                .tets
                .get(cell)
                .ok_or(invalid("mechanical cell is outside the original field"))?;
            let mut expected = source.map(|node| node as usize);
            let mut actual = [0; 4];
            for local in 0..4 {
                actual[local] = *correspondence
                    .field_vertices
                    .get(problem.tetrahedra[element][local])
                    .ok_or(invalid("mechanical tetrahedron has an absent node"))?;
            }
            expected.sort_unstable();
            actual.sort_unstable();
            if expected != actual {
                return Err(invalid(
                    "mechanical tetrahedron differs from its original field cell",
                ));
            }
            if field
                .cells()
                .binary_search_by_key(&cell, |entry| entry.cell)
                .is_ok()
            {
                conducting.push(cell);
            }
        }
        let loads = field.lorentz_loads(&conducting, controls.loads, &mut || {
            cx.checkpoint().is_err()
        })?;
        let assembly = problem.assemble(cx)?;
        let mut displacement = Vec::with_capacity(free);
        let mut velocity = Vec::with_capacity(free);
        let mut cosine_force_n = Vec::with_capacity(2 * free);
        let mut sine_force_n = Vec::with_capacity(2 * free);
        for &dof in &assembly.free_dofs {
            poll(cx)?;
            let (node, axis) = (dof / 3, dof % 3);
            let vertex = correspondence.field_vertices[node];
            displacement.push(initial_displacement_m[node][axis]);
            velocity.push(initial_velocity_m_s[node][axis]);
            cosine_force_n.extend([
                loads.mean().nodal_forces_n()[vertex][axis],
                loads.second_harmonic_real().nodal_forces_n()[vertex][axis],
            ]);
            sine_force_n.extend([
                0.0,
                loads.second_harmonic_imaginary().nodal_forces_n()[vertex][axis],
            ]);
        }
        let frequencies = [0.0, loads.second_harmonic_omega_rad_s()];
        let solid = match damping_ns_m {
            Some(damping) => ExactHarmonicSolid::prepare_viscous(
                &assembly,
                damping,
                &frequencies,
                dt_s,
                &displacement,
                &velocity,
                controls.structure,
                cx,
            )?,
            None => ExactHarmonicSolid::prepare(
                &assembly,
                &frequencies,
                dt_s,
                &displacement,
                &velocity,
                controls.structure,
                cx,
            )?,
        };
        poll(cx)?;
        Ok(Self {
            field,
            problem,
            field_vertices: correspondence.field_vertices.to_vec(),
            field_cells: correspondence.field_cells.to_vec(),
            loads,
            solid,
            next_cosine_force_n: vec![0.0; cosine_force_n.len()],
            next_sine_force_n: vec![0.0; sine_force_n.len()],
            cosine_force_n,
            sine_force_n,
        })
    }

    pub fn field(&self) -> &'a SpatialHarmonicInductionResponse<'s, 'm, C> {
        self.field
    }
    pub fn problem(&self) -> &'a TetLinearElasticProblem<'a> {
        self.problem
    }
    pub fn field_vertices(&self) -> &[usize] {
        &self.field_vertices
    }
    pub fn field_cells(&self) -> &[usize] {
        &self.field_cells
    }
    /// All original nodal forces, including supported coordinates.
    pub fn loads(&self) -> &SpatialHarmonicLorentzLoads {
        &self.loads
    }
    /// Read-only original native motion, energies, clock and free-DOF labels.
    pub fn solid(&self) -> &ExactHarmonicSolid {
        &self.solid
    }
    /// Current interval's [DC, 2*omega] coefficients per original free DOF [N].
    pub fn cosine_force_n(&self) -> &[f64] {
        &self.cosine_force_n
    }
    pub fn sine_force_n(&self) -> &[f64] {
        &self.sine_force_n
    }

    /// Advance actual motion with all DC/AC/transient cross terms in source
    /// work and physical viscosity. Native rejection preserves state, clock,
    /// accounts AND source phase. No Joule heat is added to mechanical heat.
    pub fn advance(&mut self, cx: &Cx<'_>) -> Result<ExactLinearSolidStep> {
        poll(cx)?;
        // Compute the next source in scratch before the native transaction.
        // Following accepted motion, the infallible swaps publish its phase.
        for index in 0..self.cosine_force_n.len() {
            poll(cx)?;
            let (cosine, sine) = self.solid.advance_force_phase(
                index % 2,
                self.cosine_force_n[index],
                self.sine_force_n[index],
            )?;
            self.next_cosine_force_n[index] = cosine;
            self.next_sine_force_n[index] = sine;
        }
        let step = self
            .solid
            .advance(&self.cosine_force_n, &self.sine_force_n, cx)?;
        core::mem::swap(&mut self.cosine_force_n, &mut self.next_cosine_force_n);
        core::mem::swap(&mut self.sine_force_n, &mut self.next_sine_force_n);
        Ok(step)
    }
}
