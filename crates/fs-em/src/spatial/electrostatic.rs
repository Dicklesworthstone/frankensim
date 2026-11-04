//! Three-dimensional source-bound electrostatics on the original P1 complex.
//!
//! Solve div(D)=rho, with E=-grad(V) and D=epsilon E+P_fixed.
//! The tetrahedral Whitney Hodge, geometry, trace integration and checked CG
//! are the SAME owners used by the other spatial fields. Permittivity is an
//! absolute SI dielectric source, never a conductivity or a relative constant.
//!
//! Prescribed D.n is positive OUT of the dielectric; its weak load is negative.
//! Dirichlet reactions are free electrode charge. A pure-Neumann component
//! must balance its volume charge against its outward displacement before a
//! deterministic numerical gauge is introduced. That gauge supplies no charge.
//!
//! This is a fixed-geometry, linear, nondispersive snapshot. Unlisted exterior
//! faces have D.n=0; a finite enclosing box is not an infinite-domain model.
//! Cell fields retain material jumps. Reported weak residuals and raw interface
//! jumps are not continuum error bounds or an equilibrated H(div) certificate.
//! Impressed polarization is additional to the induced response in epsilon.
//! Its maintaining source owns any changing-polarization work or storage law.
//! Formulation reference: <https://mfem.org/electromagnetics/#electrostatics>.

mod force;
mod terminals;
mod virtual_work;
pub mod network;
pub use network::{SpatialCapacitanceControls, SpatialCapacitanceMatrix, SpatialCapacitanceResponse};
pub use terminals::SpatialCapacitance;

use super::{
    ScalarBoundary, SpatialGeometry, SpatialMagneticError, dot, prepare_geometry,
    prepare_scalar_boundary_traces_with_sources, solve_scalar_system,
};
use crate::electrostatic::ElectrostaticControls;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError,
        spatial::{SpatialDielectricLaw, SpatialPermittivity},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TetComplex;
use fs_solver::SolveReport;
use fs_sparse::Csr;
use std::sync::Arc;

/// P1 voltage trace on an original exterior face, in sorted face-vertex order.
/// Shared vertices must receive the same potential from every incident trace.
#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticVoltageFace {
    pub face: usize,
    pub voltage_v: [f64; 3],
}

/// Uniform OUTWARD electric displacement D.n [C/m2] on an exterior face.
/// This is a distributed Neumann condition, not a floating equipotential metal.
#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticDisplacementFace {
    pub face: usize,
    pub outward_displacement_c_m2: f64,
}

/// Uniform free sheet charge [C/m2] on one original INTERNAL face.
/// Positive charge is independent of face orientation. This is an impressed
/// distribution, not an equipotential conductor or an inferred bound charge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticChargeSheet {
    pub face: usize,
    pub free_charge_density_c_m2: f64,
}

/// Explicit fixed sources for the original three-dimensional dielectric law.
/// Volume charge has one entry per cell. When present, impressed polarization
/// also has one mesh-frame vector per cell and is ADDED to epsilon E; it must
/// not duplicate the induced polarization already included in permittivity.
/// Sheets use distinct internal faces. Their weak P1 face moments and those of
/// constant cell sources are exact for the represented finite-element data.
#[derive(Clone, Copy, Debug)]
pub struct SpatialElectrostaticSources<'a> {
    pub free_charge_density_c_m3: &'a [f64],
    pub impressed_polarization_c_m2: Option<&'a [[f64; 3]]>,
    pub sheets: &'a [ElectrostaticChargeSheet],
}

/// Unlisted faces have zero normal displacement. Voltage and displacement
/// traces cannot occupy the same face. Disconnected components retain their
/// own reference potential and charge-compatibility check.
#[derive(Clone, Debug, Default)]
pub struct SpatialElectrostaticBoundary {
    pub voltage: Vec<ElectrostaticVoltageFace>,
    pub displacement: Vec<ElectrostaticDisplacementFace>,
}

/// No partially admitted field returns on a material, solve or balance error.
#[derive(Debug)]
pub enum SpatialElectrostaticError {
    /// Shared volume geometry, boundary, compatibility or checked-solver error.
    Native(SpatialMagneticError),
    Material(FieldStorageError),
    Assembly(WeightedError),
    Invalid(&'static str),
    Cancelled,
    Balance {
        charge_residual_c: f64,
        maximum_component_charge_residual_c: f64,
        maximum_free_charge_residual_c: f64,
        energy_residual_j: f64,
        hodge_energy_residual_j: f64,
    },
}
impl From<SpatialMagneticError> for SpatialElectrostaticError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SpatialElectrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral electrostatics: {self:?}")
    }
}
impl std::error::Error for SpatialElectrostaticError {}
type Result<T> = core::result::Result<T, SpatialElectrostaticError>;

fn invalid(reason: &'static str) -> SpatialElectrostaticError {
    SpatialElectrostaticError::Invalid(reason)
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialElectrostaticError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable electrostatic arithmetic"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

/// Immutable complete material states and their original tetrahedral operator.
/// Geometry and element gradients are prepared once and reused across sources
/// and electrode excitations. Scalar and oriented tensor cells may coexist.
pub struct SpatialElectrostatics {
    geometry: Arc<SpatialGeometry>,
    materials: Vec<ResolvedMaterialStatePoint>,
    laws: Vec<SpatialDielectricLaw>,
    permittivity: Vec<SpatialPermittivity>,
    stiffness: Csr,
}

struct PreparedSources {
    rho: Vec<f64>,
    polarization: Option<Vec<[f64; 3]>>,
    sheets: Vec<ElectrostaticChargeSheet>,
    volume_load: Vec<f64>,
    sheet_load: Vec<f64>,
    polarization_load: Vec<f64>,
    total_load: Vec<f64>,
    sheet_charge_by_face: Vec<f64>,
    charged_sheet_vertices: Vec<bool>,
}

impl PreparedSources {
    fn weak_load_c(&self) -> &[f64] {
        if self.total_load.is_empty() {
            &self.volume_load
        } else {
            &self.total_load
        }
    }
}

impl SpatialElectrostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if tetrahedra.len() != materials.len() || tetrahedra.len() > limits.max_cells {
            return Err(invalid(
                "one material per tetrahedron within the explicit cell budget required",
            ));
        }
        let laws = vec![SpatialDielectricLaw::Isotropic; materials.len()];
        Self::new_with_laws(positions_m, tetrahedra, materials, laws, limits, cancelled)
    }

    pub fn new_with_laws(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedMaterialStatePoint>,
        laws: Vec<SpatialDielectricLaw>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if tetrahedra.len() != materials.len() || laws.len() != materials.len() {
            return Err(invalid(
                "one complete dielectric source and law per tetrahedron required",
            ));
        }
        let geometry = Arc::new(prepare_geometry(
            positions_m,
            tetrahedra,
            limits,
            cancelled,
        )?);
        let mut permittivity = Vec::with_capacity(materials.len());
        let mut weights = Vec::with_capacity(materials.len());
        for (material, law) in materials.iter().zip(&laws) {
            check(cancelled)?;
            let epsilon = law
                .resolve(material)
                .map_err(SpatialElectrostaticError::Material)?;
            weights.push(match epsilon {
                SpatialPermittivity::Scalar(value) => CellWeight::Scalar(value),
                SpatialPermittivity::Tensor(matrix) => CellWeight::Tensor(matrix),
            });
            permittivity.push(epsilon);
        }
        let mass = weighted_mass_matrix(
            &geometry.mesh,
            &geometry.geometry,
            1,
            &weights,
            ELECTRIC_PERMITTIVITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => SpatialElectrostaticError::Cancelled,
            other => SpatialElectrostaticError::Assembly(other),
        })?;
        check(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&geometry.mesh.d0()), mass.matrix());
        for row in 0..stiffness.nrows() {
            check(cancelled)?;
            for &value in stiffness.row(row).1 {
                finite(value)?;
            }
        }
        Ok(Self {
            geometry,
            materials,
            laws,
            permittivity,
            stiffness,
        })
    }

    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        &self.materials
    }
    pub fn permittivity_laws(&self) -> &[SpatialDielectricLaw] {
        &self.laws
    }
    pub fn permittivity_coefficients(&self) -> &[SpatialPermittivity] {
        &self.permittivity
    }
    pub fn cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.geometry.geometry.vol_signed.get(cell).map(|v| v.abs())
    }
    pub fn exterior_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.geometry
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, face)| (face.cells.len() == 1).then_some(i))
    }
    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.geometry
            .faces
            .get(face)
            .filter(|f| f.cells.len() == 1)
            .map(|f| f.outward_normal)
    }
    pub fn face_area_m2(&self, face: usize) -> Option<f64> {
        self.geometry.faces.get(face).map(|f| f.area_m2)
    }

    /// Solve fixed per-cell free volume charge and the explicitly supplied
    /// face traces. P1 volume and face loads use exact constant-source moments.
    /// Field recovery keeps component-relative potentials throughout; adding
    /// a large common voltage does not pollute charge reactions or stored energy.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'_>> {
        self.solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3,
                impressed_polarization_c_m2: None,
                sheets: &[],
            },
            boundary,
            controls,
            cancelled,
        )
    }

    /// Solve volume/sheet free charge and fixed impressed polarization using
    /// the original immutable operator. The weak equation is
    /// `K V = q_volume + q_sheet + integral grad(N).P_fixed - integral N D.n`.
    /// Interior sheets impose the normal-D jump without smoothing either cell
    /// field. Polarization discontinuities supply their own bound charge and
    /// must not also be authored as free sheets.
    ///
    /// Reported storage remains `integral E.epsilon.E/2`, the reversible linear
    /// dielectric store at fixed P. It is not `integral E.D/2` when P is nonzero.
    /// No constitutive evolution, hysteresis, polarization-source energy or
    /// changing-material/mechanical work is inferred from these fixed inputs.
    pub fn solve_with_sources(
        &self,
        sources: SpatialElectrostaticSources<'_>,
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'_>> {
        let (sources, traces) = self.prepare_sources(sources, boundary, controls, cancelled)?;
        let rhs = sources
            .weak_load_c()
            .iter()
            .zip(&traces.flux_load_wb)
            .map(|(&volume, &surface)| finite(volume + surface))
            .collect::<Result<Vec<_>>>()?;
        let (relative, report) = solve_scalar_system(
            &self.stiffness,
            &rhs,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_field(
            relative, sources, traces, boundary, controls, report, cancelled,
        )
    }

    fn prepare_boundary(
        &self,
        rho: &[f64],
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PreparedSources, ScalarBoundary)> {
        self.prepare_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: rho,
                impressed_polarization_c_m2: None,
                sheets: &[],
            },
            boundary,
            controls,
            cancelled,
        )
    }

    fn prepare_sources(
        &self,
        sources: SpatialElectrostaticSources<'_>,
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PreparedSources, ScalarBoundary)> {
        check(cancelled)?;
        let rho = sources.free_charge_density_c_m3;
        if rho.len() != self.mesh().tets.len()
            || rho.iter().any(|v| !v.is_finite())
            || sources.impressed_polarization_c_m2.is_some_and(|p| {
                p.len() != self.mesh().tets.len() || p.iter().flatten().any(|v| !v.is_finite())
            })
            || sources.sheets.len() > self.geometry.faces.len()
            || boundary.voltage.len() > self.geometry.faces.len()
            || boundary.displacement.len() > self.geometry.faces.len()
            || !positive(controls.linear.tolerance)
            || controls.linear.tolerance >= 1.0
            || !positive(controls.charge_tolerance_c)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(invalid(
                "complete finite cell sources, bounded sheets/traces and positive tolerances required",
            ));
        }
        let g = &self.geometry;
        let mut volume_load = vec![0.0; g.mesh.vertex_count];
        let mut corrections = vec![0.0; g.mesh.vertex_count];
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let load = finite(rho[cell] * g.geometry.vol_signed[cell].abs() / 4.0)?;
            for &v in tet {
                add(
                    &mut volume_load[v as usize],
                    &mut corrections[v as usize],
                    load,
                )?;
            }
        }
        let mut sheets = sources.sheets.to_vec();
        sheets.sort_unstable_by_key(|s| s.face);
        let mut sheet_load = Vec::new();
        let mut sheet_charge_by_face = Vec::new();
        let mut charged_sheet_vertices = Vec::new();
        if !sheets.is_empty() {
            sheet_load.resize(g.mesh.vertex_count, 0.0);
            sheet_charge_by_face.resize(g.faces.len(), 0.0);
            charged_sheet_vertices.resize(g.mesh.vertex_count, false);
            corrections.fill(0.0);
        }
        for (i, sheet) in sheets.iter().enumerate() {
            check(cancelled)?;
            if g.faces.get(sheet.face).is_none_or(|f| f.cells.len() != 2)
                || i > 0 && sheets[i - 1].face == sheet.face
                || !sheet.free_charge_density_c_m2.is_finite()
            {
                return Err(invalid(
                    "finite free sheets require distinct original internal faces",
                ));
            }
            let charge = finite(sheet.free_charge_density_c_m2 * g.faces[sheet.face].area_m2)?;
            sheet_charge_by_face[sheet.face] = charge;
            for v in g.mesh.faces[sheet.face] {
                let v = v as usize;
                add(&mut sheet_load[v], &mut corrections[v], charge / 3.0)?;
                charged_sheet_vertices[v] |= sheet.free_charge_density_c_m2 != 0.0;
            }
        }
        // No extra nodal/face arrays are needed for the original volume-only
        // path. When sheets are present, reuse their complete free-charge
        // buffer for the eventual polarized weak load after gauge admission.
        let mut total_load = if sheets.is_empty() {
            Vec::new()
        } else {
            volume_load
                .iter()
                .zip(&sheet_load)
                .map(|(&volume, &sheet)| finite(volume + sheet))
                .collect::<Result<Vec<_>>>()?
        };
        let traces = prepare_scalar_boundary_traces_with_sources(
            &g.mesh,
            &g.faces,
            &g.components,
            boundary.voltage.iter().map(|b| (b.face, b.voltage_v)),
            boundary
                .displacement
                .iter()
                .map(|b| (b.face, -b.outward_displacement_c_m2)),
            // Polarization creates no net free charge. Check the physical
            // compatibility of volume + sheet charge against outward D.n;
            // never let a numerical gauge hide uncompensated free charge.
            Some(if total_load.is_empty() {
                &volume_load
            } else {
                &total_load
            }),
            controls.charge_tolerance_c,
            cancelled,
        )?;
        let mut polarization_load = Vec::new();
        if let Some(polarization) = sources.impressed_polarization_c_m2 {
            polarization_load.resize(g.mesh.vertex_count, 0.0);
            corrections.fill(0.0);
            for (cell, tet) in g.mesh.tets.iter().enumerate() {
                check(cancelled)?;
                let (mut first, mut first_correction) = (0.0, 0.0);
                // The constant test function has zero gradient. Form its
                // zeroth contribution from the other three, matching the
                // voltage-difference gradient used by field recovery.
                for local in 1..4 {
                    let term = finite(
                        g.geometry.vol_signed[cell].abs()
                            * dot(g.geometry.grads[cell][local], polarization[cell]),
                    )?;
                    let v = tet[local] as usize;
                    add(&mut polarization_load[v], &mut corrections[v], term)?;
                    add(&mut first, &mut first_correction, -term)?;
                }
                let v = tet[0] as usize;
                add(&mut polarization_load[v], &mut corrections[v], first)?;
            }
            if total_load.is_empty() {
                total_load.clone_from(&volume_load);
            }
            for (load, &polarized) in total_load.iter_mut().zip(&polarization_load) {
                check(cancelled)?;
                *load = finite(*load + polarized)?;
            }
        }
        check(cancelled)?;
        Ok((
            PreparedSources {
                rho: rho.to_vec(),
                polarization: sources.impressed_polarization_c_m2.map(<[_]>::to_vec),
                sheets,
                volume_load,
                sheet_load,
                polarization_load,
                total_load,
                sheet_charge_by_face,
                charged_sheet_vertices,
            },
            traces,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn recover_field(
        &self,
        relative: Vec<f64>,
        sources: PreparedSources,
        traces: ScalarBoundary,
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'_>> {
        let ScalarBoundary {
            physical_prescribed,
            flux_load_wb: surface_load,
            gauges,
            references,
            ..
        } = traces;
        let g = &self.geometry;
        let n = g.mesh.vertex_count;
        let mut nodal_boundary_charge_c = vec![None; n];
        let mut component_charge = vec![0.0; n];
        let mut component_correction = vec![0.0; n];
        let (mut volume_charge_c, mut volume_correction) = (0.0, 0.0);
        let (mut sheet_charge_c, mut sheet_correction) = (0.0, 0.0);
        let (mut boundary_charge_c, mut boundary_correction) = (0.0, 0.0);
        let (mut pairing_j, mut pairing_correction) = (0.0, 0.0);
        let (mut hodge_j, mut hodge_correction) = (0.0, 0.0);
        let mut maximum_free_charge_residual_c = 0.0_f64;
        for v in 0..n {
            check(cancelled)?;
            let (columns, values) = self.stiffness.row(v);
            let (mut applied, mut correction) = (0.0, 0.0);
            for (&w, &coefficient) in columns.iter().zip(values) {
                add(
                    &mut applied,
                    &mut correction,
                    finite(coefficient * finite(relative[w] - relative[v])?)?,
                )?;
            }
            add(
                &mut hodge_j,
                &mut hodge_correction,
                finite(0.5 * relative[v] * applied)?,
            )?;
            let sheet_load = sources.sheet_load.get(v).copied().unwrap_or(0.0);
            let polarization_load = sources.polarization_load.get(v).copied().unwrap_or(0.0);
            let free_load = finite(sources.volume_load[v] + sheet_load)?;
            let free_applied = finite(applied - polarization_load)?;
            let charge = if physical_prescribed[v] {
                // Includes prescribed displacement on adjacent faces ONCE.
                Some(finite(free_applied - free_load)?)
            } else {
                let residual =
                    finite(applied - finite(sources.weak_load_c()[v] + surface_load[v])?)?;
                maximum_free_charge_residual_c = maximum_free_charge_residual_c.max(residual.abs());
                (surface_load[v] != 0.0).then_some(surface_load[v])
            };
            let boundary_charge = charge.unwrap_or(0.0);
            let supplied = if physical_prescribed[v] {
                // Keep the exact operator reaction instead of subtracting then
                // re-adding potentially much larger prescribed source loads.
                applied
            } else {
                finite(sources.weak_load_c()[v] + boundary_charge)?
            };
            add(
                &mut pairing_j,
                &mut pairing_correction,
                finite(0.5 * relative[v] * supplied)?,
            )?;
            add(
                &mut volume_charge_c,
                &mut volume_correction,
                sources.volume_load[v],
            )?;
            add(&mut sheet_charge_c, &mut sheet_correction, sheet_load)?;
            add(
                &mut boundary_charge_c,
                &mut boundary_correction,
                boundary_charge,
            )?;
            let component = g.components[v];
            let free_supplied = if physical_prescribed[v] {
                free_applied
            } else {
                finite(free_load + boundary_charge)?
            };
            add(
                &mut component_charge[component],
                &mut component_correction[component],
                free_supplied,
            )?;
            nodal_boundary_charge_c[v] = charge;
        }
        let mut electric_field_v_m = Vec::with_capacity(self.materials.len());
        let mut displacement_c_m2 = Vec::with_capacity(self.materials.len());
        let mut cell_energy_j = Vec::with_capacity(self.materials.len());
        let (mut energy_j, mut energy_correction) = (0.0, 0.0);
        let (mut polarization_pairing_j, mut polarization_pairing_correction) = (0.0, 0.0);
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let mut e = [0.0; 3];
            for local in 1..4 {
                let difference = finite(relative[tet[local] as usize] - relative[tet[0] as usize])?;
                for (axis, out) in e.iter_mut().enumerate() {
                    *out =
                        finite((-difference).mul_add(g.geometry.grads[cell][local][axis], *out))?;
                }
            }
            let epsilon_e = self.permittivity[cell].displacement_c_m2(e);
            let d = if let Some(p) = &sources.polarization {
                add(
                    &mut polarization_pairing_j,
                    &mut polarization_pairing_correction,
                    finite(dot(e, p[cell]) * g.geometry.vol_signed[cell].abs())?,
                )?;
                std::array::from_fn(|axis| epsilon_e[axis] + p[cell][axis])
            } else {
                epsilon_e
            };
            for value in d {
                finite(value)?;
            }
            // Fixed impressed polarization is a source, not an additional
            // linear storage coefficient. In an open-circuit polarized body
            // D may vanish while the depolarizing field stores positive energy.
            let energy = finite(0.5 * dot(e, epsilon_e) * g.geometry.vol_signed[cell].abs())?;
            if energy < 0.0 || (e != [0.0; 3] && energy == 0.0) {
                return Err(invalid(
                    "nonzero dielectric field lost representably positive storage",
                ));
            }
            add(&mut energy_j, &mut energy_correction, energy)?;
            electric_field_v_m.push(e);
            displacement_c_m2.push(d);
            cell_energy_j.push(energy);
        }
        let charge_balance_residual_c =
            finite(finite(volume_charge_c + sheet_charge_c)? + boundary_charge_c)?;
        let maximum_component_charge_residual_c =
            component_charge.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let energy_balance_residual_j = finite(energy_j - pairing_j)?;
        let hodge_energy_residual_j = finite(energy_j - hodge_j)?;
        if charge_balance_residual_c.abs() > controls.charge_tolerance_c
            || maximum_component_charge_residual_c > controls.charge_tolerance_c
            || maximum_free_charge_residual_c > controls.charge_tolerance_c
            || energy_balance_residual_j.abs() > controls.energy_tolerance_j
            || hodge_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(SpatialElectrostaticError::Balance {
                charge_residual_c: charge_balance_residual_c,
                maximum_component_charge_residual_c,
                maximum_free_charge_residual_c,
                energy_residual_j: energy_balance_residual_j,
                hodge_energy_residual_j,
            });
        }
        let mut interface_displacement_jumps_c = Vec::new();
        let mut interface_displacement_residuals_c = Vec::new();
        for (face, data) in g.faces.iter().enumerate() {
            check(cancelled)?;
            if let [left, right] = data.cells.as_slice() {
                let jump = std::array::from_fn(|axis| {
                    displacement_c_m2[*left][axis] - displacement_c_m2[*right][axis]
                });
                let jump = finite(dot(jump, data.outward_normal) * data.area_m2)?;
                interface_displacement_jumps_c.push((face, jump));
                // n points from the first cell to the second. Gauss' pillbox
                // gives (D_second-D_first).n = sigma_free on a sheet.
                if !sources.sheets.is_empty() {
                    interface_displacement_residuals_c
                        .push((face, finite(jump + sources.sheet_charge_by_face[face])?));
                }
            }
        }
        let mut voltage_v = relative
            .iter()
            .enumerate()
            .map(|(v, &value)| finite(value + references[g.components[v]]))
            .collect::<Result<Vec<_>>>()?;
        // Preserve caller-supplied physical traces bit-for-bit for reporting;
        // recovered fields continue to use the unrounded relative solution.
        for trace in &boundary.voltage {
            for (&v, &value) in g.mesh.faces[trace.face].iter().zip(&trace.voltage_v) {
                voltage_v[v as usize] = value;
            }
        }
        check(cancelled)?;
        Ok(SpatialElectricField {
            problem: self,
            free_charge_density_c_m3: sources.rho,
            impressed_polarization_c_m2: sources.polarization,
            sheets: sources.sheets,
            charged_sheet_vertices: sources.charged_sheet_vertices,
            has_imposed_displacement: boundary.displacement.iter()
                .any(|trace| trace.outward_displacement_c_m2 != 0.0),
            voltage_v,
            relative_voltage_v: relative,
            electric_field_v_m,
            displacement_c_m2,
            cell_energy_j,
            energy_j,
            nodal_boundary_charge_c,
            nodal_volume_charge_c: sources.volume_load,
            nodal_sheet_charge_c: sources.sheet_load,
            nodal_polarization_load_c: sources.polarization_load,
            nodal_applied_surface_charge_c: surface_load,
            volume_charge_c,
            sheet_charge_c,
            boundary_charge_c,
            polarization_pairing_j,
            charge_balance_residual_c,
            maximum_component_charge_residual_c,
            maximum_free_charge_residual_c,
            energy_balance_residual_j,
            hodge_energy_residual_j,
            interface_displacement_jumps_c,
            interface_displacement_residuals_c,
            gauges,
            report,
        })
    }
}

/// Accepted field retaining its exact source/geometry owner. Relative voltages
/// preserve a small physical signal even when an absolute output rounds under
/// a very large common bias. Field, charge and energy use the relative solution.
pub struct SpatialElectricField<'a> {
    problem: &'a SpatialElectrostatics,
    free_charge_density_c_m3: Vec<f64>,
    impressed_polarization_c_m2: Option<Vec<[f64; 3]>>,
    sheets: Vec<ElectrostaticChargeSheet>,
    charged_sheet_vertices: Vec<bool>,
    // Preserve boundary intent even when distinct face traces cancel nodally.
    has_imposed_displacement: bool,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 3]>,
    displacement_c_m2: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    nodal_boundary_charge_c: Vec<Option<f64>>,
    nodal_volume_charge_c: Vec<f64>,
    nodal_sheet_charge_c: Vec<f64>,
    nodal_polarization_load_c: Vec<f64>,
    nodal_applied_surface_charge_c: Vec<f64>,
    volume_charge_c: f64,
    sheet_charge_c: f64,
    boundary_charge_c: f64,
    polarization_pairing_j: f64,
    charge_balance_residual_c: f64,
    maximum_component_charge_residual_c: f64,
    maximum_free_charge_residual_c: f64,
    energy_balance_residual_j: f64,
    hodge_energy_residual_j: f64,
    interface_displacement_jumps_c: Vec<(usize, f64)>,
    interface_displacement_residuals_c: Vec<(usize, f64)>,
    gauges: Vec<usize>,
    report: Option<SolveReport>,
}
impl SpatialElectricField<'_> {
    pub fn problem(&self) -> &SpatialElectrostatics {
        self.problem
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    /// Original additional fixed polarization in the mesh frame. None means
    /// exactly zero; induced response remains represented by epsilon.
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 3]]> {
        self.impressed_polarization_c_m2.as_deref()
    }
    /// Original impressed free sheets, in canonical face order.
    pub fn charge_sheets(&self) -> &[ElectrostaticChargeSheet] {
        &self.sheets
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.electric_field_v_m
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 3]] {
        &self.displacement_c_m2
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    /// Reversible linear dielectric storage, integral E.epsilon.E/2, at fixed
    /// impressed polarization. No polarization-source energy is included.
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Free electrode charge, including any supplied face-displacement charge
    /// at shared trace vertices. Free unforced nodes carry None, including gauges.
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        &self.nodal_boundary_charge_c
    }
    pub fn nodal_volume_charge_c(&self) -> &[f64] {
        &self.nodal_volume_charge_c
    }
    /// Exact P1 moments of internal free sheet charge [C]. Empty means no
    /// sheets were supplied, hence exactly zero at every node.
    pub fn nodal_sheet_charge_c(&self) -> &[f64] {
        &self.nodal_sheet_charge_c
    }
    /// Nodal weak load integral grad(N).P_fixed [C]. This is not free charge.
    /// Empty means no impressed polarization was supplied, hence exactly zero.
    pub fn nodal_polarization_load_c(&self) -> &[f64] {
        &self.nodal_polarization_load_c
    }
    /// Signed free surface load = -integral(N_i D.n dA) [C].
    pub fn nodal_applied_surface_charge_c(&self) -> &[f64] {
        &self.nodal_applied_surface_charge_c
    }
    pub fn volume_charge_c(&self) -> f64 {
        self.volume_charge_c
    }
    pub fn sheet_charge_c(&self) -> f64 {
        self.sheet_charge_c
    }
    pub fn boundary_charge_c(&self) -> f64 {
        self.boundary_charge_c
    }
    /// Integral E.P_fixed [J], evaluated on the original cells. With the
    /// complete free-charge load, the weak identity is
    /// `2 energy_j = phi.q_free - polarization_pairing_j`.
    /// This signed pairing is not a separately inferred energy of the
    /// apparatus maintaining P or time-integrated changing-polarization work.
    pub fn polarization_pairing_j(&self) -> f64 {
        self.polarization_pairing_j
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
    pub fn energy_balance_residual_j(&self) -> f64 {
        self.energy_balance_residual_j
    }
    pub fn hodge_energy_residual_j(&self) -> f64 {
        self.hodge_energy_residual_j
    }
    /// Raw P1 normal-D discontinuity integrated on each interior face [C].
    /// The normal points from the first adjacent cell to the second; an exact
    /// free sheet contributes minus its charge to this first-minus-second jump.
    pub fn interface_displacement_jumps_c(&self) -> &[(usize, f64)] {
        &self.interface_displacement_jumps_c
    }
    /// Raw interface jump PLUS prescribed free sheet charge [C]. Zero is the
    /// resolved jump condition; P1 discretization residuals remain unsmoothed.
    /// This is a diagnostic, not an equilibrated flux or continuum error bound.
    pub fn interface_displacement_residuals_c(&self) -> &[(usize, f64)] {
        if self.interface_displacement_residuals_c.is_empty() {
            &self.interface_displacement_jumps_c
        } else {
            &self.interface_displacement_residuals_c
        }
    }
    pub fn gauge_vertices(&self) -> &[usize] {
        &self.gauges
    }
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}
