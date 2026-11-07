//! Three-dimensional multipole dielectric memory and Ohmic charge relaxation.
//!
//! The original tetrahedral Whitney Hodge supplies the scalar potential
//! operator. Cell-local vector polarization, total-current continuity and
//! energy recovery use the same evolution kernel as planar/axisymmetric fields.
//! There are no additional global polarization unknowns. Fixed-size steps
//! reuse their prepared endpoint operator and the original material sources.
//!
//! Every connected component needs a positive-area exterior voltage trace.
//! Other nodal free charges are initialized from the supplied potential and
//! polarization, then evolve by total-current continuity. No additional
//! impressed charge or current is silently inserted. Prescribed traces exchange
//! both stored electrode charge and transported Ohmic charge with their sources.
//! Cell heat separates polarization loss from Joule loss; field storage is
//! never reclassified as heat. Geometry and material conditions remain fixed.
//!
//! The mesh admission and finite-domain limitations are those of
//! [`super::SpatialElectrostatics`]. This isotropic electroquasistatic model
//! does not include induction, carrier transport or breakdown. Midpoint is
//! second order and passive, but large steps can oscillate; relaxation times
//! require independent temporal resolution.
//! Model: <https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>.

use std::sync::Arc;

/// Inertial dielectric memory on the same native tetrahedral geometry.
pub mod lorentz;

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::field_storage::dielectric::bulk::ResolvedDebyeMaterial;
#[cfg(feature = "resolved-leaky-dielectric")]
use fs_material::field_storage::dielectric::leaky::ResolvedLeakyDielectric;
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

use super::super::{SpatialGeometry, SpatialMagneticError, prepare_geometry};
use crate::electrostatic::{
    ElectrostaticControls, ElectrostaticError,
    dispersion::{
        DielectricFieldCandidate, DielectricFieldState,
        kernel::{DebyeGeometry, Field, Prepared},
        sealed,
    },
};

type Result<T> = core::result::Result<T, ElectrostaticError>;

fn invalid(reason: &'static str) -> ElectrostaticError {
    ElectrostaticError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "tetrahedral dielectric arithmetic is not representable",
        ))
    }
}
fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(ElectrostaticError::Cancelled)
    } else {
        Ok(())
    }
}
fn accumulate(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
fn geometry_error(error: SpatialMagneticError) -> ElectrostaticError {
    match error {
        SpatialMagneticError::Cancelled => ElectrostaticError::Cancelled,
        SpatialMagneticError::Invalid(reason) => invalid(reason),
        SpatialMagneticError::Assembly(error) => ElectrostaticError::Assembly(error),
        _ => invalid("invalid or unrepresentable tetrahedral dielectric geometry"),
    }
}

struct TetrahedralGeometry {
    native: SpatialGeometry,
    volumes_m3: Vec<f64>,
}
impl TetrahedralGeometry {
    fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Arc<Self>> {
        let native =
            prepare_geometry(positions_m, tetrahedra, limits, cancelled).map_err(geometry_error)?;
        let mut volumes_m3 = Vec::with_capacity(native.mesh.tets.len());
        for &volume in &native.geometry.vol_signed {
            poll(cancelled)?;
            volumes_m3.push(finite(volume.abs())?);
        }
        Ok(Arc::new(Self { native, volumes_m3 }))
    }

    fn boundary(
        &self,
        voltages: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>> {
        poll(cancelled)?;
        let n = self.native.mesh.vertex_count;
        if voltages.len() > n {
            return Err(invalid(
                "voltage trace exceeds the tetrahedral vertex count",
            ));
        }
        let mut prescribed = vec![None; n];
        for &(vertex, voltage) in voltages {
            poll(cancelled)?;
            if vertex >= n || !voltage.is_finite() || prescribed[vertex].is_some() {
                return Err(invalid(
                    "distinct valid voltage vertices and finite values required",
                ));
            }
            prescribed[vertex] = Some(voltage);
        }
        // A physical voltage trace is a union of complete exterior triangles.
        // A point, edge, interior node or numerical gauge is not an electrode.
        let mut covered = vec![false; n];
        let mut anchored = vec![false; n];
        for (vertices, face) in self.native.mesh.faces.iter().zip(&self.native.faces) {
            poll(cancelled)?;
            if face.cells.len() == 1 && vertices.iter().all(|&v| prescribed[v as usize].is_some()) {
                for &vertex in vertices {
                    covered[vertex as usize] = true;
                    anchored[self.native.components[vertex as usize]] = true;
                }
            }
        }
        for vertex in 0..n {
            poll(cancelled)?;
            if prescribed[vertex].is_some() && !covered[vertex] {
                return Err(invalid(
                    "voltage vertices must cover complete exterior triangular faces",
                ));
            }
            if !anchored[self.native.components[vertex]] {
                return Err(invalid(
                    "every dielectric component needs a physical voltage trace",
                ));
            }
        }
        Ok(prescribed)
    }

    fn initial_boundary(
        &self,
        voltage_v: &[f64],
        voltage_vertices: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>> {
        if voltage_v.len() != self.native.mesh.vertex_count
            || voltage_vertices.len() > voltage_v.len()
            || voltage_vertices.iter().any(|&v| v >= voltage_v.len())
        {
            return Err(invalid(
                "complete initial voltage state and valid trace vertices required",
            ));
        }
        let boundary: Vec<_> = voltage_vertices
            .iter()
            .map(|&v| (v, voltage_v[v]))
            .collect();
        self.boundary(&boundary, cancelled)
    }
}

impl DebyeGeometry<3> for TetrahedralGeometry {
    fn vertex_count(&self) -> usize {
        self.native.mesh.vertex_count
    }
    fn cell_count(&self) -> usize {
        self.native.mesh.tets.len()
    }
    fn volumes_m3(&self) -> &[f64] {
        &self.volumes_m3
    }
    fn assemble(
        &self,
        coefficients: &[f64],
        dimensions: Dims,
        limits: WeightedAssemblyLimits,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Csr> {
        let weights: Vec<_> = coefficients
            .iter()
            .copied()
            .map(CellWeight::Scalar)
            .collect();
        let mass = weighted_mass_matrix(
            &self.native.mesh,
            &self.native.geometry,
            1,
            &weights,
            dimensions,
            limits,
            &mut || cancelled(),
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => ElectrostaticError::Cancelled,
            other => ElectrostaticError::Assembly(other),
        })?;
        poll(cancelled)?;
        let matrix = stiffness(&incidence_to_csr(&self.native.mesh.d0()), mass.matrix());
        for row in 0..matrix.nrows() {
            poll(cancelled)?;
            if matrix.row(row).1.iter().any(|value| !value.is_finite()) {
                return Err(invalid(
                    "tetrahedral dielectric operator is not representable",
                ));
            }
        }
        Ok(matrix)
    }
    fn fields(
        &self,
        voltage: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>> {
        let mut fields = Vec::with_capacity(self.cell_count());
        for (vertices, gradients) in self
            .native
            .mesh
            .tets
            .iter()
            .zip(&self.native.geometry.grads)
        {
            poll(cancelled)?;
            let mut field = [0.0; 3];
            for local in 1..4 {
                let difference =
                    finite(voltage[vertices[local] as usize] - voltage[vertices[0] as usize])?;
                for axis in 0..3 {
                    field[axis] =
                        finite((-difference).mul_add(gradients[local][axis], field[axis]))?;
                }
            }
            fields.push(field);
        }
        Ok(fields)
    }
    fn polarization_load(
        &self,
        polarization: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let mut rhs = vec![0.0; self.vertex_count()];
        let mut corrections = vec![0.0; rhs.len()];
        for (cell, vertices) in self.native.mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            let p = polarization[cell];
            let gradients = self.native.geometry.grads[cell];
            let volume = self.volumes_m3[cell];
            let mut load = [0.0; 4];
            for local in 1..4 {
                let g = gradients[local];
                load[local] =
                    finite(volume * finite(g[0].mul_add(p[0], g[1].mul_add(p[1], g[2] * p[2])))?)?;
            }
            // Constant cell polarization has an exact P1 weak moment. Close
            // partition of unity locally before reducing into shared vertices.
            load[0] = -finite(finite(load[1] + load[2])? + load[3])?;
            for (&vertex, load) in vertices.iter().zip(load) {
                accumulate(
                    &mut rhs[vertex as usize],
                    &mut corrections[vertex as usize],
                    load,
                )?;
            }
        }
        for (value, correction) in rhs.iter_mut().zip(corrections) {
            *value = finite(*value + correction)?;
        }
        Ok(rhs)
    }
}

/// Accepted physical 3D dielectric field with retained vector pole history.
/// Clones share immutable geometry/materials; trial steps never mutate history.
#[derive(Clone)]
pub struct SpatialDebyeDielectric {
    geometry: Arc<TetrahedralGeometry>,
    core: Field<3>,
}

/// Reusable endpoint operator for one fixed tetrahedral material/geometry owner.
pub struct PreparedSpatialDebyeStep {
    core: Prepared<3>,
}
impl PreparedSpatialDebyeStep {
    pub fn time_step_s(&self) -> f64 {
        self.core.time_step_s()
    }
}
/// Transactional candidate sharing the same physical receipts as other metrics.
pub type SpatialDebyeFieldCandidate = DielectricFieldCandidate<SpatialDebyeDielectric, 3>;

impl SpatialDebyeDielectric {
    /// Initialize an explicitly insulating dielectric. The supplied voltages
    /// and cell/pole vectors define the initial weak free charge, including any
    /// precharge. Each component needs a complete exterior voltage face.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedDebyeMaterial>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; 3]>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let geometry = TetrahedralGeometry::new(positions_m, tetrahedra, limits, cancelled)?;
        let prescribed = geometry.initial_boundary(&voltage_v, voltage_vertices, cancelled)?;
        let core = Field::new(
            Arc::clone(&geometry) as Arc<dyn DebyeGeometry<3>>,
            materials,
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )?;
        Ok(Self { geometry, core })
    }

    /// Initialize independent DC leakage and Debye storage from one complete
    /// source bundle per tetrahedron. Sourced zero conductivity is admitted
    /// exactly, without a conductivity floor or a separate singular solve.
    #[cfg(feature = "resolved-leaky-dielectric")]
    #[allow(clippy::too_many_arguments)]
    pub fn new_leaky(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedLeakyDielectric>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; 3]>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let geometry = TetrahedralGeometry::new(positions_m, tetrahedra, limits, cancelled)?;
        let prescribed = geometry.initial_boundary(&voltage_v, voltage_vertices, cancelled)?;
        let core = Field::new_leaky(
            Arc::clone(&geometry) as Arc<dyn DebyeGeometry<3>>,
            materials,
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )?;
        Ok(Self { geometry, core })
    }
    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.native.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.native.positions_m
    }
    pub fn materials(&self) -> &[ResolvedDebyeMaterial] {
        self.core.materials()
    }
    pub fn voltage_v(&self) -> &[f64] {
        self.core.voltage_v()
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.core.relative_voltage_v()
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        self.core.electric_field_v_m()
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 3]] {
        self.core.displacement_c_m2()
    }
    pub fn polarization_c_m2(&self) -> &[Vec<[f64; 3]>] {
        self.core.states()
    }
    pub fn nodal_charge_c(&self) -> &[f64] {
        self.core.nodal_charge_c()
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        self.core.cell_energy_j()
    }
    pub fn cell_volumes_m3(&self) -> &[f64] {
        self.core.cell_volumes_m3()
    }
    pub fn conductivity_s_m(&self) -> &[f64] {
        self.core.conductivity_s_m()
    }
    pub fn conductivity_properties(&self) -> &[Option<String>] {
        self.core.conductivity_properties()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.core.stored_energy_j()
    }
    pub fn time_s(&self) -> f64 {
        self.core.time_s()
    }
    pub fn same_problem(&self, other: &Self) -> bool {
        self.core.same_problem(&other.core)
    }

    /// Prepare once and reuse at the same dt without repeating Hodge assembly.
    pub fn prepare_step(
        &self,
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedSpatialDebyeStep> {
        Ok(PreparedSpatialDebyeStep {
            core: self.core.prepare_step(dt_s, cancelled)?,
        })
    }

    /// Advance with linearly ramped endpoint voltage traces. The original
    /// voltage-vertex set is fixed; cancellation/refusal leaves history intact.
    pub fn step_candidate(
        &self,
        prepared: &PreparedSpatialDebyeStep,
        endpoint_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialDebyeFieldCandidate> {
        let prescribed = self.geometry.boundary(endpoint_voltage_v, cancelled)?;
        let (core, step) =
            self.core
                .step_candidate(&prepared.core, prescribed, controls, cancelled)?;
        Ok(DielectricFieldCandidate {
            field: Self {
                geometry: Arc::clone(&self.geometry),
                core,
            },
            step,
        })
    }
}

impl sealed::Sealed for SpatialDebyeDielectric {}
impl DielectricFieldState for SpatialDebyeDielectric {
    type Material = ResolvedDebyeMaterial;
    fn materials(&self) -> &[ResolvedDebyeMaterial] {
        self.core.materials()
    }
    fn cell_volumes_m3(&self) -> &[f64] {
        self.core.cell_volumes_m3()
    }
    fn conductivity_s_m(&self) -> &[f64] {
        self.core.conductivity_s_m()
    }
    fn conductivity_properties(&self) -> &[Option<String>] {
        self.core.conductivity_properties()
    }
    fn stored_energy_j(&self) -> f64 {
        self.core.stored_energy_j()
    }
    fn time_s(&self) -> f64 {
        self.core.time_s()
    }
    fn same_problem(&self, other: &Self) -> bool {
        self.core.same_problem(&other.core)
    }
}
