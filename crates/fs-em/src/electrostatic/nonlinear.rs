//! Reversible nonlinear electrostatics on the original planar P1 complex.
//!
//! The voltage functional uses integral D.dE; physical field storage uses
//! integral E.dD. They coincide only for linear response. Floating electrodes
//! retain their total charge and one shared potential, solved with the interior
//! field; no unit-voltage capacitance superposition substitutes for this state.
//!
//! Geometry is prepared once. The analytic material tangent is applied cell by
//! cell through the original checked-CG owner, with diagonal equilibration and
//! bounded coenergy backtracking. No global tangent matrix is assembled.
//! These are fixed-source electrostatics and discrete balance checks, not
//! dielectric memory, thermal evolution, electrostriction or spatial error bounds.

mod electrodes;
pub mod force;
mod solve;
pub mod tangent;

use super::{ElectrostaticError, finite, geometry, poll, positive};
use crate::{planar::gradients, scalar_potential};
use fs_feec::{
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::CellWeight2,
};
use fs_material::{
    field_storage::{
        dielectric::reversible::ReversibleDielectricMaterial,
        planar::{PlanarDielectricLaw, PlanarPermittivity},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, SolveReport};
use std::sync::Arc;

pub use electrodes::{ElectrodeControl, NonlinearElectrodes};

/// Original source laws may coexist in one field. Linear anisotropy retains
/// its explicit material/mesh frame; nonlinear response here is isotropic.
#[derive(Clone, Debug, PartialEq)]
pub enum DielectricMaterial {
    Linear {
        material: ResolvedMaterialStatePoint,
        law: PlanarDielectricLaw,
    },
    Reversible(ReversibleDielectricMaterial),
}
impl DielectricMaterial {
    pub fn isotropic(material: ResolvedMaterialStatePoint) -> Self {
        Self::Linear {
            material,
            law: PlanarDielectricLaw::Isotropic,
        }
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Linear { material, .. } => material,
            Self::Reversible(law) => law.material(),
        }
    }
}

/// Independent nonlinear accuracy, physical admission and finite work limits.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearElectrostaticControls {
    pub absolute_residual_tolerance_c: f64,
    pub relative_residual_tolerance: f64,
    pub charge_tolerance_c: f64,
    /// Residual of integral(E.D) against actual volume/electrode pairing [J].
    pub energy_tolerance_j: f64,
    pub max_newton_iterations: usize,
    pub max_backtracks: usize,
    /// Across all original checked-CG calls, including residual recomputation.
    pub max_operator_applications: usize,
    pub linear: CheckedCgConfig,
}

#[derive(Clone, Debug)]
pub struct NonlinearElectrostaticReport {
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub backtracks: usize,
    pub operator_applications: usize,
    pub residual_norm_c: f64,
    pub residual_target_c: f64,
}

#[derive(Debug)]
pub enum NonlinearElectrostaticError {
    Electric(ElectrostaticError),
    FieldRange {
        cell: usize,
        observed_v_m: f64,
        maximum_v_m: f64,
    },
    Nonconvergence(NonlinearElectrostaticReport),
    OperatorBudget,
}
impl From<ElectrostaticError> for NonlinearElectrostaticError {
    fn from(error: ElectrostaticError) -> Self {
        Self::Electric(error)
    }
}
impl From<scalar_potential::PotentialError> for NonlinearElectrostaticError {
    fn from(error: scalar_potential::PotentialError) -> Self {
        Self::Electric(error.into())
    }
}
impl core::fmt::Display for NonlinearElectrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear electrostatics: {self:?}")
    }
}
impl std::error::Error for NonlinearElectrostaticError {}
type Error = NonlinearElectrostaticError;
fn invalid(reason: &'static str) -> Error {
    ElectrostaticError::Invalid(reason).into()
}
fn weighted(error: WeightedError) -> Error {
    match error {
        WeightedError::Cancelled => ElectrostaticError::Cancelled.into(),
        other => ElectrostaticError::Assembly(other).into(),
    }
}
fn add(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), Error> {
    let term = finite(term)?;
    let next = finite(*sum + term)?;
    *correction = finite(
        *correction
            + if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
fn norm(values: &[f64]) -> f64 {
    values.iter().copied().fold(0.0, f64::hypot)
}

#[derive(Clone, Copy)]
struct Element {
    vertices: [usize; 3],
    /// Gradients of local shape functions 1 and 2. Differences from vertex0
    /// retain the exact constant-potential nullspace in every kernel.
    gradient: [[f64; 2]; 2],
    volume_m3: f64,
}
impl Element {
    fn gradient_of(&self, values: [f64; 3]) -> [f64; 2] {
        let a = values[1] - values[0];
        let b = values[2] - values[0];
        [
            self.gradient[0][0].mul_add(a, self.gradient[1][0] * b),
            self.gradient[0][1].mul_add(a, self.gradient[1][1] * b),
        ]
    }
    fn field(&self, potential: &[f64]) -> Result<[f64; 2], Error> {
        let g = self.gradient_of(self.vertices.map(|v| potential[v]));
        Ok([-finite(g[0])?, -finite(g[1])?])
    }
    fn reactions(&self, displacement: [f64; 2]) -> Result<[f64; 3], Error> {
        let a = finite(
            -self.volume_m3
                * self.gradient[0][0]
                    .mul_add(displacement[0], self.gradient[0][1] * displacement[1]),
        )?;
        let b = finite(
            -self.volume_m3
                * self.gradient[1][0]
                    .mul_add(displacement[0], self.gradient[1][1] * displacement[1]),
        )?;
        Ok([-finite(a + b)?, a, b])
    }
}

/// Prepared original mesh and source laws; cloned material curves share tables.
pub struct PlanarNonlinearElectrostatics {
    mesh: TriComplex2,
    materials: Vec<DielectricMaterial>,
    linear: Vec<Option<PlanarPermittivity>>,
    elements: Vec<Element>,
}

struct Evaluation {
    electric: Vec<[f64; 2]>,
    displacement: Vec<[f64; 2]>,
    tangent: Vec<[[f64; 2]; 2]>,
    cell_energy: Vec<f64>,
    cell_coenergy: Vec<f64>,
    nodal: Vec<f64>,
    correction: Vec<f64>,
    energy: f64,
    coenergy: f64,
}
impl Evaluation {
    fn new(problem: &PlanarNonlinearElectrostatics) -> Self {
        let c = problem.elements.len();
        let n = problem.mesh.vertices().len();
        Self {
            electric: vec![[0.0; 2]; c],
            displacement: vec![[0.0; 2]; c],
            tangent: vec![[[0.0; 2]; 2]; c],
            cell_energy: vec![0.0; c],
            cell_coenergy: vec![0.0; c],
            nodal: vec![0.0; n],
            correction: vec![0.0; n],
            energy: 0.0,
            coenergy: 0.0,
        }
    }
}

impl PlanarNonlinearElectrostatics {
    /// `max_triplets` bounds the nine local tangent contributions per triangle;
    /// no global triplet or CSR tangent storage is allocated. Existing native
    /// geometry and interval coefficient admission remain the shared owners.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<DielectricMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. })
            || mesh.faces().is_empty()
            || materials.len() != mesh.faces().len()
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|v| v > limits.max_triplets)
        {
            return Err(invalid(
                "bounded planar mesh and one original material law per triangle required",
            ));
        }
        let mut linear = Vec::with_capacity(materials.len());
        let mut elements = Vec::with_capacity(materials.len());
        let mut frame = None;
        for (cell, (material, face)) in materials.iter().zip(mesh.faces()).enumerate() {
            poll(cancelled)?;
            let coefficient = match material {
                DielectricMaterial::Linear { material, law } => {
                    if let Some(axes) = law.frame() {
                        if frame.is_some_and(|id| id != axes.geometry_frame()) {
                            return Err(invalid(
                                "oriented dielectric cells must name one mesh frame",
                            ));
                        }
                        frame = Some(axes.geometry_frame());
                    }
                    Some(
                        law.resolve(material)
                            .map_err(ElectrostaticError::Material)?,
                    )
                }
                DielectricMaterial::Reversible(law) => {
                    let epsilon = law
                        .evaluate(0.0)
                        .map_err(ElectrostaticError::Material)?
                        .differential_permittivity_f_m;
                    CellWeight2::Scalar(epsilon)
                        .validate(cell, 1)
                        .map_err(weighted)?;
                    None
                }
            };
            if let Some(epsilon) = coefficient {
                match epsilon {
                    PlanarPermittivity::Scalar(v) => CellWeight2::Scalar(v),
                    PlanarPermittivity::Tensor(v) => CellWeight2::Tensor(v),
                }
                .validate(cell, 1)
                .map_err(weighted)?;
            }
            let g = gradients(face.map(|v| mesh.vertices()[v as usize])).map_err(geometry)?;
            let volume = mesh
                .face_measure(cell)
                .ok_or_else(|| invalid("missing dielectric cell volume"))?;
            if !positive(volume) {
                return Err(invalid("positive finite dielectric cell volume required"));
            }
            elements.push(Element {
                vertices: face.map(|v| v as usize),
                gradient: [g[1], g[2]],
                volume_m3: volume,
            });
            linear.push(coefficient);
        }
        poll(cancelled)?;
        Ok(Self {
            mesh,
            materials,
            linear,
            elements,
        })
    }
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    pub fn materials(&self) -> &[DielectricMaterial] {
        &self.materials
    }

    fn evaluate(
        &self,
        potential: &[f64],
        out: &mut Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        out.nodal.fill(0.0);
        out.correction.fill(0.0);
        let (mut energy, mut ce, mut coenergy, mut cc) = (0.0, 0.0, 0.0, 0.0);
        for (cell, element) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let e = element.field(potential)?;
            let magnitude = finite(e[0].hypot(e[1]))?;
            let (d, tangent, u, w) = if let Some(linear) = self.linear[cell] {
                let d = linear.displacement_c_m2(e);
                let d = [finite(d[0])?, finite(d[1])?];
                let u = finite(0.5 * e[0].mul_add(d[0], e[1] * d[1]))?;
                (d, linear.matrix_f_m(), u, u)
            } else {
                let DielectricMaterial::Reversible(law) = &self.materials[cell] else {
                    unreachable!()
                };
                if magnitude > law.maximum_field_strength_v_m() {
                    return Err(Error::FieldRange {
                        cell,
                        observed_v_m: magnitude,
                        maximum_v_m: law.maximum_field_strength_v_m(),
                    });
                }
                let response = law
                    .evaluate(magnitude)
                    .map_err(ElectrostaticError::Material)?;
                let s = response.secant_permittivity_f_m;
                let p = response.differential_permittivity_f_m;
                let (d, t) = if magnitude == 0.0 {
                    ([0.0; 2], [[p, 0.0], [0.0, p]])
                } else {
                    let n = [e[0] / magnitude, e[1] / magnitude];
                    let xy = (p - s) * n[0] * n[1];
                    (
                        [
                            finite(response.displacement_c_m2 * n[0])?,
                            finite(response.displacement_c_m2 * n[1])?,
                        ],
                        [
                            [p * n[0] * n[0] + s * n[1] * n[1], xy],
                            [xy, p * n[1] * n[1] + s * n[0] * n[0]],
                        ],
                    )
                };
                (
                    d,
                    t,
                    response.energy_density_j_m3,
                    response.coenergy_density_j_m3,
                )
            };
            let u = finite(element.volume_m3 * u)?;
            let w = finite(element.volume_m3 * w)?;
            if u < 0.0
                || w < 0.0
                || (magnitude > 0.0 && (u == 0.0 || w == 0.0))
                || tangent.iter().flatten().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "dielectric storage or tangent is not representable",
                ));
            }
            for (vertex, q) in element.vertices.into_iter().zip(element.reactions(d)?) {
                add(&mut out.nodal[vertex], &mut out.correction[vertex], q)?;
            }
            out.electric[cell] = e;
            out.displacement[cell] = d;
            out.tangent[cell] = tangent;
            out.cell_energy[cell] = u;
            out.cell_coenergy[cell] = w;
            add(&mut energy, &mut ce, u)?;
            add(&mut coenergy, &mut cc, w)?;
        }
        for (v, c) in out.nodal.iter_mut().zip(&out.correction) {
            *v = finite(*v + c)?;
        }
        out.energy = finite(energy + ce)?;
        out.coenergy = finite(coenergy + cc)?;
        Ok(())
    }

    fn coenergy_change(
        &self,
        before: &Evaluation,
        after: &Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let (mut sum, mut correction) = (0.0, 0.0);
        for (cell, element) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let a = before.electric[cell];
            let b = after.electric[cell];
            let change = if let Some(linear) = self.linear[cell] {
                let delta = [b[0] - a[0], b[1] - a[1]];
                let midpoint = [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])];
                let d = linear.displacement_c_m2(midpoint);
                let mut value = delta[0].mul_add(d[0], delta[1] * d[1]);
                // The represented midpoint may differ from a+delta/2.
                let offset = [
                    0.5 * delta[0] - (midpoint[0] - a[0]),
                    0.5 * delta[1] - (midpoint[1] - a[1]),
                ];
                let rounding = linear.displacement_c_m2(offset);
                value += delta[0].mul_add(rounding[0], delta[1] * rounding[1]);
                finite(value)?
            } else {
                let DielectricMaterial::Reversible(law) = &self.materials[cell] else {
                    unreachable!()
                };
                law.coenergy_change_vector_density_j_m3(a, b)
                    .map_err(ElectrostaticError::Material)?
            };
            add(&mut sum, &mut correction, element.volume_m3 * change)?;
        }
        Ok(finite(sum + correction)?)
    }
}

/// Actual nonlinear field at its operating point. Relative potentials preserve
/// small fields under a large common electrode bias. The stored tangent belongs
/// to this accepted state and cannot be treated as a global Q=C*V law.
pub struct NonlinearElectricField<'m> {
    problem: &'m PlanarNonlinearElectrostatics,
    electrodes: Arc<electrodes::Layout>,
    controls: Vec<ElectrodeControl>,
    dofs: electrodes::Dofs,
    reference_terminal: usize,
    reference_voltage_v: f64,
    relative_voltage_v: Vec<f64>,
    terminal_voltage_v: Vec<f64>,
    terminal_charge_c: Vec<f64>,
    nodal_boundary_charge_c: Vec<Option<f64>>,
    free_charge_density_c_m3: Vec<f64>,
    evaluation: Evaluation,
    volume_charge_c: f64,
    charge_balance_residual_c: f64,
    maximum_free_charge_residual_c: f64,
    maximum_floating_charge_residual_c: f64,
    energy_identity_residual_j: f64,
    report: NonlinearElectrostaticReport,
}
impl<'m> NonlinearElectricField<'m> {
    pub fn problem(&self) -> &'m PlanarNonlinearElectrostatics {
        self.problem
    }
    pub fn reference_terminal(&self) -> usize {
        self.reference_terminal
    }
    pub fn reference_voltage_v(&self) -> f64 {
        self.reference_voltage_v
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn terminal_voltage_v(&self) -> &[f64] {
        &self.terminal_voltage_v
    }
    pub fn terminal_charge_c(&self) -> &[f64] {
        &self.terminal_charge_c
    }
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        &self.nodal_boundary_charge_c
    }
    /// Original cell sources, retained for force admission and explicit
    /// transport under a virtual geometry change.
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    /// Actual voltage and total-charge controls of this accepted field.
    pub fn electrode_controls(&self) -> &[ElectrodeControl] {
        &self.controls
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.evaluation.electric
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        &self.evaluation.displacement
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.evaluation.cell_energy
    }
    pub fn cell_coenergy_j(&self) -> &[f64] {
        &self.evaluation.cell_coenergy
    }
    pub fn energy_j(&self) -> f64 {
        self.evaluation.energy
    }
    pub fn coenergy_j(&self) -> f64 {
        self.evaluation.coenergy
    }
    pub fn volume_charge_c(&self) -> f64 {
        self.volume_charge_c
    }
    pub fn charge_balance_residual_c(&self) -> f64 {
        self.charge_balance_residual_c
    }
    pub fn maximum_free_charge_residual_c(&self) -> f64 {
        self.maximum_free_charge_residual_c
    }
    pub fn maximum_floating_charge_residual_c(&self) -> f64 {
        self.maximum_floating_charge_residual_c
    }
    pub fn energy_identity_residual_j(&self) -> f64 {
        self.energy_identity_residual_j
    }
    pub fn report(&self) -> &NonlinearElectrostaticReport {
        &self.report
    }
}
