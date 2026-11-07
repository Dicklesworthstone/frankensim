//! Nonlinear Gauss fields on the original three-dimensional tetrahedral mesh.
//!
//! Fixed-voltage and total-charge electrodes share the original nonlinear
//! iteration and electrical-component gauges. Each metal is a union of complete
//! exterior faces; a floating electrode retains one unknown voltage and its
//! requested total free charge. Unlisted faces have D.n=0.
//!
//! Native tetrahedral geometry and source-resolved scalar/oriented linear laws
//! coexist with the original isotropic reversible D(|E|) curves. Physical
//! storage is integral E.dD, coenergy is integral D.dE, and the analytic tangent
//! is applied without a global matrix. Source bounds are never extrapolated.
//! This fixed-geometry snapshot has no polarization history, charge transport,
//! imposed polarization, sheet charge or infinite-domain boundary inference.

mod constitutive;
pub mod force;
mod solve;
pub mod tangent;

use super::electrodes::{Dofs, Layout, VoltageReferences};
use super::{
    ElectrodeControl, ElectrostaticError, Error, NonlinearElectrostaticControls,
    NonlinearElectrostaticReport, add, finite, gauss, invalid, poll, positive, weighted,
};
use crate::spatial::{SpatialGeometry, SpatialMagneticError, prepare_geometry};
use fs_feec::weighted::{CellWeight, WeightedAssemblyLimits};
use fs_material::{
    field_storage::{
        dielectric::reversible::ReversibleDielectricMaterial,
        spatial::{SpatialDielectricLaw, SpatialPermittivity},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TetComplex;
use fs_solver::{CheckedCgConfig, SolveReport};
use std::sync::Arc;

/// Original fixed-condition dielectric law for each original tetrahedron.
/// Nonlinear curves are isotropic; linear tensors retain all three principal
/// source properties and the caller's admitted material frame.
#[derive(Clone, Debug, PartialEq)]
pub enum SpatialDielectricMaterial {
    Linear {
        material: ResolvedMaterialStatePoint,
        law: SpatialDielectricLaw,
    },
    Reversible(ReversibleDielectricMaterial),
}
impl SpatialDielectricMaterial {
    pub fn isotropic(material: ResolvedMaterialStatePoint) -> Self {
        Self::Linear {
            material,
            law: SpatialDielectricLaw::Isotropic,
        }
    }

    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Linear { material, .. } => material,
            Self::Reversible(law) => law.material(),
        }
    }
}

fn native(error: SpatialMagneticError) -> Error {
    match error {
        SpatialMagneticError::Cancelled => ElectrostaticError::Cancelled.into(),
        SpatialMagneticError::Invalid(reason) => invalid(reason),
        SpatialMagneticError::Material(error) | SpatialMagneticError::Recoil { error, .. } => {
            ElectrostaticError::Material(error).into()
        }
        SpatialMagneticError::Assembly(error) => weighted(error),
        SpatialMagneticError::Nonconvergence(report) => {
            ElectrostaticError::Nonconvergence(report).into()
        }
        SpatialMagneticError::Geometry(_) => invalid("invalid original tetrahedral geometry"),
    }
}

fn trace_error(error: crate::spatial::electrostatic::SpatialElectrostaticError) -> Error {
    use crate::spatial::electrostatic::SpatialElectrostaticError;
    match error {
        SpatialElectrostaticError::Native(error) => native(error),
        SpatialElectrostaticError::Material(error) => ElectrostaticError::Material(error).into(),
        SpatialElectrostaticError::Assembly(error) => weighted(error),
        SpatialElectrostaticError::Invalid(reason) => invalid(reason),
        SpatialElectrostaticError::Cancelled => ElectrostaticError::Cancelled.into(),
        SpatialElectrostaticError::Balance { .. } => invalid("invalid native electrode trace"),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

/// A borrowed view of the native P1 geometry, with no duplicate gradient cache.
#[derive(Clone, Copy)]
struct Element<'a> {
    vertices: &'a [u32; 4],
    gradients: &'a [[f64; 3]; 4],
    volume_m3: f64,
}
impl Element<'_> {
    fn gradient_of(&self, values: [f64; 4]) -> [f64; 3] {
        let differences = [
            values[1] - values[0],
            values[2] - values[0],
            values[3] - values[0],
        ];
        std::array::from_fn(|axis| {
            differences[0].mul_add(
                self.gradients[1][axis],
                differences[1].mul_add(
                    self.gradients[2][axis],
                    differences[2] * self.gradients[3][axis],
                ),
            )
        })
    }

    fn field(&self, potential: &[f64]) -> Result<[f64; 3], Error> {
        let gradient = self.gradient_of(self.vertices.map(|v| potential[v as usize]));
        Ok([
            -finite(gradient[0])?,
            -finite(gradient[1])?,
            -finite(gradient[2])?,
        ])
    }

    fn reactions(&self, displacement: [f64; 3]) -> Result<[f64; 4], Error> {
        let mut reaction = [0.0; 4];
        for local in 1..4 {
            reaction[local] = finite(-self.volume_m3 * dot(self.gradients[local], displacement))?;
        }
        // Construct the constant-potential nullspace from the same represented
        // contributions, rather than relying on a rounded sum of gradients.
        reaction[0] = -finite(reaction[1] + reaction[2] + reaction[3])?;
        Ok(reaction)
    }
}

/// Prepared original tetrahedral geometry and source laws, shared across drives.
/// Native gradients are retained once; cloning shares the geometry and curves.
#[derive(Clone)]
pub struct SpatialNonlinearElectrostatics {
    geometry: Arc<SpatialGeometry>,
    materials: Vec<SpatialDielectricMaterial>,
    linear: Arc<[Option<SpatialPermittivity>]>,
}
impl SpatialNonlinearElectrostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<SpatialDielectricMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if materials.len() != tetrahedra.len() {
            return Err(invalid(
                "one original dielectric law per tetrahedron required",
            ));
        }
        let geometry =
            prepare_geometry(positions_m, tetrahedra, limits, cancelled).map_err(native)?;
        let mut linear = Vec::with_capacity(materials.len());
        for (cell, material) in materials.iter().enumerate() {
            poll(cancelled)?;
            let coefficient = match material {
                SpatialDielectricMaterial::Linear { material, law } => {
                    let coefficient = law
                        .resolve(material)
                        .map_err(ElectrostaticError::Material)?;
                    match coefficient {
                        SpatialPermittivity::Scalar(value) => CellWeight::Scalar(value),
                        SpatialPermittivity::Tensor(tensor) => CellWeight::Tensor(tensor),
                    }
                    .validate(cell, 1)
                    .map_err(weighted)?;
                    Some(coefficient)
                }
                SpatialDielectricMaterial::Reversible(law) => {
                    let response = law.evaluate(0.0).map_err(ElectrostaticError::Material)?;
                    CellWeight::Scalar(response.differential_permittivity_f_m)
                        .validate(cell, 1)
                        .map_err(weighted)?;
                    None
                }
            };
            if !positive(geometry.elements().vol_signed[cell].abs()) {
                return Err(invalid(
                    "positive finite tetrahedral dielectric volume required",
                ));
            }
            linear.push(coefficient);
        }
        poll(cancelled)?;
        Ok(Self {
            geometry: Arc::new(geometry),
            materials,
            linear: linear.into(),
        })
    }

    pub fn mesh(&self) -> &TetComplex {
        self.geometry.mesh()
    }

    pub fn positions_m(&self) -> &[[f64; 3]] {
        self.geometry.positions_m()
    }

    pub fn materials(&self) -> &[SpatialDielectricMaterial] {
        &self.materials
    }

    pub fn exterior_faces(&self) -> Vec<usize> {
        self.geometry.exterior_faces()
    }

    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.geometry.exterior_normal(face)
    }

    fn element(&self, cell: usize) -> Element<'_> {
        Element {
            vertices: &self.mesh().tets[cell],
            gradients: &self.geometry.elements().grads[cell],
            volume_m3: self.geometry.elements().vol_signed[cell].abs(),
        }
    }

    /// Each input is a nonempty set of original exterior FACE indices.
    /// Distinct metals cannot share a vertex. One metal can join disconnected
    /// dielectric volumes; every resulting electrical component needs a
    /// voltage-controlled reference when solving.
    pub fn prepare_electrodes<'m>(
        &'m self,
        terminals: &[&[usize]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearElectrodes<'m>, Error> {
        poll(cancelled)?;
        if terminals.is_empty() || terminals.len() > self.mesh().vertex_count {
            return Err(invalid("at least one bounded physical electrode required"));
        }
        let (faces, terminals) = crate::spatial::electrostatic::electrode_faces_on_geometry(
            &self.geometry,
            terminals,
            cancelled,
        )
        .map_err(trace_error)?;
        let mut components = self.geometry.components().to_vec();
        fn root(parent: &mut [usize], mut v: usize) -> usize {
            while parent[v] != v {
                parent[v] = parent[parent[v]];
                v = parent[v];
            }
            v
        }
        let mut node_terminal = vec![None; self.mesh().vertex_count];
        for (terminal, trace) in terminals.iter().enumerate() {
            for &v in trace {
                poll(cancelled)?;
                let a = root(&mut components, trace[0]);
                let b = root(&mut components, v);
                components[a.max(b)] = a.min(b);
                node_terminal[v] = Some(terminal);
            }
        }
        for v in 0..components.len() {
            poll(cancelled)?;
            components[v] = root(&mut components, v);
        }
        poll(cancelled)?;
        Ok(SpatialNonlinearElectrodes {
            problem: self,
            faces: Arc::new(faces),
            layout: Arc::new(Layout {
                terminals,
                node_terminal,
                components,
            }),
        })
    }
}

/// Original exterior electrodes, reused across voltage/charge operating points.
pub struct SpatialNonlinearElectrodes<'m> {
    problem: &'m SpatialNonlinearElectrostatics,
    faces: Arc<Vec<Vec<usize>>>,
    layout: Arc<Layout>,
}
impl<'m> SpatialNonlinearElectrodes<'m> {
    pub fn problem(&self) -> &'m SpatialNonlinearElectrostatics {
        self.problem
    }

    pub fn terminal_faces(&self) -> &[Vec<usize>] {
        &self.faces
    }

    pub fn terminal_vertices(&self) -> &[Vec<usize>] {
        &self.layout.terminals
    }

    /// Solve div(D)=rho with exact P1 moments of one constant density per cell.
    /// Initial voltages are relative to EACH electrical component's reference;
    /// floating trace vertices must share one initial voltage. The initial
    /// constitutive field must already lie inside every original source range.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        controls: &[ElectrodeControl],
        initial_relative_voltage_v: &[f64],
        settings: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearElectricField<'m>, Error> {
        solve::solve(
            self,
            free_charge_density_c_m3,
            controls,
            initial_relative_voltage_v,
            settings,
            cancelled,
        )
    }
}

pub(super) struct Evaluation {
    electric: Vec<[f64; 3]>,
    displacement: Vec<[f64; 3]>,
    tangent: Vec<[[f64; 3]; 3]>,
    cell_energy: Vec<f64>,
    cell_coenergy: Vec<f64>,
    nodal: Vec<f64>,
    correction: Vec<f64>,
    energy: f64,
    coenergy: f64,
}

/// Accepted original 3D field, physical free charge and nonlinear storage.
/// The retained tangent describes this operating point, never a global Q=C V.
pub struct SpatialNonlinearElectricField<'m> {
    problem: &'m SpatialNonlinearElectrostatics,
    faces: Arc<Vec<Vec<usize>>>,
    electrodes: Arc<Layout>,
    controls: Vec<ElectrodeControl>,
    dofs: Dofs,
    references: VoltageReferences,
    relative_voltage_v: Vec<f64>,
    terminal_voltage_v: Vec<f64>,
    terminal_charge_c: Vec<f64>,
    nodal_boundary_charge_c: Vec<Option<f64>>,
    free_charge_density_c_m3: Vec<f64>,
    evaluation: Evaluation,
    volume_charge_c: f64,
    charge_balance_residual_c: f64,
    maximum_component_charge_residual_c: f64,
    maximum_free_charge_residual_c: f64,
    maximum_floating_charge_residual_c: f64,
    energy_identity_residual_j: f64,
    report: NonlinearElectrostaticReport,
}
impl<'m> SpatialNonlinearElectricField<'m> {
    pub fn problem(&self) -> &'m SpatialNonlinearElectrostatics {
        self.problem
    }
    pub fn terminal_faces(&self) -> &[Vec<usize>] {
        &self.faces
    }
    pub fn terminal_vertices(&self) -> &[Vec<usize>] {
        &self.electrodes.terminals
    }
    pub fn electrode_controls(&self) -> &[ElectrodeControl] {
        &self.controls
    }
    /// First physical reference; other components retain their own references.
    pub fn reference_terminal(&self) -> usize {
        self.references.first_terminal
    }
    pub fn reference_voltage_v(&self) -> f64 {
        self.references.first_voltage_v
    }
    pub fn reference_terminal_by_vertex(&self) -> &[usize] {
        &self.references.terminal_by_vertex
    }
    pub fn reference_voltage_by_vertex_v(&self) -> &[f64] {
        &self.references.voltage_by_vertex_v
    }
    pub fn reference_terminal_for_terminal(&self, terminal: usize) -> Option<usize> {
        self.electrodes
            .terminals
            .get(terminal)
            .map(|trace| self.references.terminal_by_vertex[trace[0]])
    }
    pub fn reference_voltage_for_terminal_v(&self, terminal: usize) -> Option<f64> {
        self.electrodes
            .terminals
            .get(terminal)
            .map(|trace| self.references.voltage_by_vertex_v[trace[0]])
    }
    /// Voltages relative to each vertex's independent physical reference.
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
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.evaluation.electric
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 3]] {
        &self.evaluation.displacement
    }
    /// Analytic dD/dE at the accepted bias. Source knots retain the source's
    /// one-sided convention; terminal differentiation separately admits kinks.
    pub fn differential_permittivity_f_m(&self) -> &[[[f64; 3]; 3]] {
        &self.evaluation.tangent
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
    pub fn maximum_component_charge_residual_c(&self) -> f64 {
        self.maximum_component_charge_residual_c
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
