//! Source-resolved, current-carrying 3D Whitney edge magnetostatics.
//!
//! K=d1^T M_(nu,2) d1, B=curl(A), H=nu(B-Br). Impressed current is complete
//! piecewise-constant J [A/m2] on the original tets; its Whitney load is exact.
//! Exterior tangential-A traces are edge integrals [Wb]. Unlisted exterior
//! faces have the natural n cross H=0 condition, not an infinite-domain model.
//!
//! Gauge uses the original cohomology owner's coexact metric construction.
//! BOTH auxiliary face projection and physical curl-curl are semidefinite in
//! their full spaces. CG is admitted only for compatible invariant subspaces,
//! starting at zero; every returned solve has a recomputed true residual.
//! The source must be orthogonal to the ENTIRE curl kernel, including handles.
//! No incompatible source is silently projected away. The final coexact A is
//! M-orthogonal to all curl-free variations and must preserve every face flux.
//!
//! This is a static, lowest-order, isotropic linear/fixed-recoil rung. There is
//! no Coulomb saddle/HX/higher-order/eddy-current/saturation or TEAM claim.
//! Weak-form reference: https://mfem.org/tesla-notes/ .

pub mod nonlinear;

use super::{
    SpatialGeometry, SpatialMagneticError, SpatialMagneticMaterial, cross, dot, finite, invalid,
    poll, prepare_geometry,
};
use crate::planar::permanent::MagneticEnergyReference;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    cohomology::coexact::{CoexactProjectionError, PreparedCoexact1},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::field_storage::{linear_magnetic_permeability, permanent::RecoilOperatingPoint};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::Precond};
use std::collections::BTreeSet;
type Result<T> = core::result::Result<T, SpatialMagneticError>;

/// Tangential A on one exterior face in sorted canonical edge order
/// [(vertex0,vertex1),(vertex0,vertex2),(vertex1,vertex2)]. Values are Wb.
#[derive(Clone, Copy, Debug)]
pub struct VectorPotentialFace {
    pub face: usize,
    pub edge_integrals_wb: [f64; 3],
}
#[derive(Clone, Debug, Default)]
pub struct VectorMagneticBoundary {
    pub potential: Vec<VectorPotentialFace>,
}
#[derive(Clone, Copy, Debug)]
pub struct VectorMagneticControls {
    pub field: CheckedCgConfig,
    pub gauge: CheckedCgConfig,
    /// Sum of squared active face incidences per free edge, before SpGEMM.
    pub maximum_gauge_product_terms: usize,
    /// Physical edge residual, current compatibility and interface current [A].
    pub current_tolerance_a: f64,
    /// Gauge-preserved face flux, reconstructed flux and div(B) allowance [Wb].
    pub flux_tolerance_wb: f64,
    pub energy_tolerance_j: f64,
}
struct CellVectors {
    edges: [usize; 6],
    curls: [[f64; 3]; 6],
    integrals: [[f64; 3]; 6],
}

/// Fixed geometry/material curl-curl preparation, reused across current loads.
pub struct SpatialVectorMagnetostatics {
    pub(super) geometry: SpatialGeometry,
    materials: Vec<SpatialMagneticMaterial>,
    permeability: Vec<f64>,
    remanence: Vec<[f64; 3]>,
    cells: Vec<CellVectors>,
    curl: Csr,
    stiffness: Csr,
    remanence_load: Vec<f64>,
}
fn gauge_error(e: CoexactProjectionError) -> SpatialMagneticError {
    match e {
        CoexactProjectionError::Cancelled => SpatialMagneticError::Cancelled,
        CoexactProjectionError::Invalid(s) => invalid(s),
    }
}
fn apply(matrix: &Csr, x: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
    let mut y = Vec::with_capacity(matrix.nrows());
    for r in 0..matrix.nrows() {
        poll(cancelled)?;
        let (c, v) = matrix.row(r);
        y.push(finite(c.iter().zip(v).map(|(&c, &v)| v * x[c]).sum())?);
    }
    Ok(y)
}
struct BorrowedMatrix<'a>(&'a Csr);
impl LinearOp for BorrowedMatrix<'_> {
    fn n(&self) -> usize {
        self.0.nrows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.spmv(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.spmv(x, y);
    }
}
struct Jacobi(Vec<f64>);
impl Precond for Jacobi {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((z, r), d) in z.iter_mut().zip(r).zip(&self.0) {
            *z = r * d;
        }
    }
}
/// The caller establishes RHS compatibility with this PSD operator's range.
/// Zero-start PCG is positive definite on its generated compatible Krylov
/// subspace; no full-space SPD or unique auxiliary face solution is asserted.
fn compatible_solve(
    matrix: &Csr,
    rhs: &[f64],
    controls: CheckedCgConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Option<SolveReport>)> {
    poll(cancelled)?;
    if rhs.is_empty() {
        return Ok((vec![], None));
    }
    let mut inverse = Vec::with_capacity(rhs.len());
    for i in 0..rhs.len() {
        poll(cancelled)?;
        let (c, v) = matrix.row(i);
        let diagonal = c.binary_search(&i).ok().map(|j| v[j]).unwrap_or(0.0);
        if !diagonal.is_finite() || diagonal <= 0.0 || !diagonal.recip().is_finite() {
            return Err(invalid(
                "positive finite diagonal required on active curl/projection rows",
            ));
        }
        inverse.push(diagonal.recip());
    }
    let solution = checked_cg(
        &BorrowedMatrix(matrix),
        &Jacobi(inverse),
        rhs,
        controls,
        |_| poll(cancelled),
    )
    .map_err(|e| match e {
        CheckedCgError::InvalidInput(s) => invalid(s),
        CheckedCgError::Interrupted(e) => e,
    })?;
    if !solution.report.converged_euclidean() {
        return Err(SpatialMagneticError::Nonconvergence(solution.report));
    }
    Ok((solution.x, Some(solution.report)))
}
fn project(
    projector: &PreparedCoexact1,
    x: &[f64],
    controls: CheckedCgConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Vec<f64>, Option<SolveReport>)> {
    let rhs = projector.rhs(x, cancelled).map_err(gauge_error)?;
    let (w, report) = compatible_solve(projector.matrix(), &rhs, controls, cancelled)?;
    let primal = projector.recover(&w, cancelled).map_err(gauge_error)?;
    let dual = projector.recover_dual(&w, cancelled).map_err(gauge_error)?;
    Ok((primal, dual, report))
}
fn max_difference(a: &[f64], b: &[f64]) -> Result<f64> {
    a.iter().zip(b).try_fold(0.0_f64, |maximum, (a, b)| {
        Ok(maximum.max(finite(a - b)?.abs()))
    })
}
impl SpatialVectorMagnetostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<SpatialMagneticMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if tetrahedra.len() != materials.len() {
            return Err(invalid("complete original magnetic material map required"));
        }
        let geometry = prepare_geometry(positions_m, tetrahedra, limits, cancelled)?;
        let mut permeability = Vec::with_capacity(materials.len());
        let mut remanence = Vec::with_capacity(materials.len());
        let mut weights = Vec::with_capacity(materials.len());
        for material in &materials {
            poll(cancelled)?;
            let (mu, br) = match material {
                SpatialMagneticMaterial::Linear(m) => (linear_magnetic_permeability(m)?, [0.0; 3]),
                SpatialMagneticMaterial::Recoil(m) => {
                    (m.material().permeability_h_m(), m.remanence_t())
                }
            };
            let nu = finite(mu.recip())?;
            if nu <= 0.0 {
                return Err(invalid(
                    "positive representable magnetic reluctivity required",
                ));
            }
            permeability.push(mu);
            remanence.push(br);
            weights.push(CellWeight::Scalar(nu));
        }
        let mass = weighted_mass_matrix(
            &geometry.mesh,
            &geometry.geometry,
            2,
            &weights,
            Dims([-1, -1, 2, 0, 2, 0]),
            limits,
            cancelled,
        )
        .map_err(|e| {
            if e == WeightedError::Cancelled {
                SpatialMagneticError::Cancelled
            } else {
                SpatialMagneticError::Assembly(e)
            }
        })?;
        let curl = incidence_to_csr(&geometry.mesh.d1());
        let stiffness = stiffness(&curl, mass.matrix());
        let cells = prepare_cells(&geometry, cancelled)?;
        let mut remanence_load = vec![0.0; geometry.mesh.edges.len()];
        for (cell, vectors) in cells.iter().enumerate() {
            poll(cancelled)?;
            let volume = geometry.geometry.vol_signed[cell].abs();
            for local in 0..6 {
                let edge = vectors.edges[local];
                let load = finite(
                    volume * dot(vectors.curls[local], remanence[cell]) / permeability[cell],
                )?;
                remanence_load[edge] = finite(remanence_load[edge] + load)?;
            }
        }
        poll(cancelled)?;
        Ok(Self {
            geometry,
            materials,
            permeability,
            remanence,
            cells,
            curl,
            stiffness,
            remanence_load,
        })
    }
    pub const fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn materials(&self) -> &[SpatialMagneticMaterial] {
        &self.materials
    }
    pub fn exterior_faces(&self) -> Vec<usize> {
        self.geometry
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| (f.cells.len() == 1).then_some(i))
            .collect()
    }
    pub fn face_edges(&self, face: usize) -> Option<[usize; 3]> {
        let [a, b, c] = *self.mesh().faces.get(face)?;
        let mut result = [0; 3];
        for (i, key) in [[a, b], [a, c], [b, c]].iter().enumerate() {
            result[i] = self.mesh().edges.binary_search(key).ok()?;
        }
        Some(result)
    }
    /// Solve the authored static current and boundary problem. Natural faces
    /// cannot carry impressed current through the exterior; potential faces
    /// may, with compatibility checked against the actual relative curl kernel.
    /// Internal normal-current continuity is required for the P0 current field.
    pub fn solve(
        &self,
        current_a_m2: &[[f64; 3]],
        boundary: &VectorMagneticBoundary,
        controls: VectorMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialVectorField<'_>> {
        let VectorSources {
            prescribed,
            free,
            projector,
            current_load,
            maximum_current_flux_defect_a,
            source_compatibility_defect_a,
            source_report,
        } = VectorTopology {
            geometry: &self.geometry,
            cells: &self.cells,
        }
        .prepare_sources(current_a_m2, boundary, controls, cancelled)?;
        let ne = self.mesh().edges.len();
        let mut index = vec![usize::MAX; ne];
        for (i, &e) in free.iter().enumerate() {
            index[e] = i;
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        let mut rhs = Vec::with_capacity(free.len());
        for &edge in &free {
            poll(cancelled)?;
            let mut load = finite(current_load[edge] + self.remanence_load[edge])?;
            let (cols, vals) = self.stiffness.row(edge);
            for (&col, &val) in cols.iter().zip(vals) {
                if let Some(a) = prescribed[col] {
                    load = finite((-val).mul_add(a, load))?;
                } else {
                    columns.push(index[col]);
                    values.push(val);
                }
            }
            row_ptr.push(columns.len());
            rhs.push(load);
        }
        let reduced = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
        let (raw, field_report) = compatible_solve(&reduced, &rhs, controls.field, cancelled)?;
        let (gauged, _, gauge_report) = project(&projector, &raw, controls.gauge, cancelled)?;
        let raw_flux = projector.rhs(&raw, cancelled).map_err(gauge_error)?;
        let gauge_flux = projector.rhs(&gauged, cancelled).map_err(gauge_error)?;
        let gauge_flux_defect_wb = max_difference(&raw_flux, &gauge_flux)?;
        if gauge_flux_defect_wb > controls.flux_tolerance_wb {
            return Err(invalid("coexact gauge changed physical magnetic face flux"));
        }
        let mut potential_wb: Vec<_> = prescribed.iter().map(|v| v.unwrap_or(0.0)).collect();
        for (&edge, &value) in free.iter().zip(&gauged) {
            potential_wb[edge] = value;
        }
        let applied = apply(&self.stiffness, &potential_wb, cancelled)?;
        let mut residual_a = Vec::with_capacity(ne);
        let mut maximum_residual_a = 0.0_f64;
        for e in 0..ne {
            poll(cancelled)?;
            let residual = finite(applied[e] - current_load[e] - self.remanence_load[e])?;
            if prescribed[e].is_none() {
                maximum_residual_a = maximum_residual_a.max(residual.abs());
            }
            residual_a.push(residual);
        }
        if maximum_residual_a > controls.current_tolerance_a {
            return Err(invalid(
                "full physical curl-curl residual exceeds the current allowance",
            ));
        }
        let face_flux_wb = apply(&self.curl, &potential_wb, cancelled)?;
        let mut maximum_divergence_flux_wb = 0.0_f64;
        for row in self.mesh().d2().rows {
            poll(cancelled)?;
            let flux = finite(
                row.iter()
                    .map(|&(f, s)| f64::from(s) * face_flux_wb[f])
                    .sum(),
            )?;
            maximum_divergence_flux_wb = maximum_divergence_flux_wb.max(flux.abs());
        }
        if maximum_divergence_flux_wb > controls.flux_tolerance_wb {
            return Err(invalid(
                "discrete magnetic flux fails the original div-curl identity allowance",
            ));
        }
        let mut b = Vec::with_capacity(self.cells.len());
        let mut h = Vec::with_capacity(self.cells.len());
        let mut cell_energy_j = Vec::with_capacity(self.cells.len());
        let mut energy_j = 0.0;
        let mut curl_energy_j = 0.0;
        let mut recoil_operating_points = Vec::with_capacity(self.cells.len());
        for (cell, vectors) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let mut flux = [0.0; 3];
            for local in 0..6 {
                for (j, value) in flux.iter_mut().enumerate() {
                    *value = finite(
                        potential_wb[vectors.edges[local]].mul_add(vectors.curls[local][j], *value),
                    )?;
                }
            }
            let field = std::array::from_fn(|j| {
                (flux[j] - self.remanence[cell][j]) / self.permeability[cell]
            });
            for value in field {
                finite(value)?;
            }
            let operating = if let SpatialMagneticMaterial::Recoil(source) = &self.materials[cell] {
                let parallel = finite(dot(field, source.direction()))?;
                let transverse = super::norm(std::array::from_fn(|j| {
                    field[j] - parallel * source.direction()[j]
                }));
                Some(
                    source
                        .material()
                        .admit_operating_point(parallel, transverse)
                        .map_err(|error| SpatialMagneticError::Recoil { cell, error })?,
                )
            } else {
                None
            };
            let volume = self.geometry.geometry.vol_signed[cell].abs();
            let energy = finite(0.5 * self.permeability[cell] * dot(field, field) * volume)?;
            curl_energy_j =
                finite(curl_energy_j + 0.5 * dot(flux, flux) * volume / self.permeability[cell])?;
            energy_j = finite(energy_j + energy)?;
            cell_energy_j.push(energy);
            b.push(flux);
            h.push(field);
            recoil_operating_points.push(operating);
        }
        let hodge_energy_j = finite(
            0.5 * potential_wb
                .iter()
                .zip(&applied)
                .map(|(a, b)| a * b)
                .sum::<f64>(),
        )?;
        if (curl_energy_j - hodge_energy_j).abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "recovered magnetic energy disagrees with the original curl-curl Hodge",
            ));
        }
        let mut maximum_recovered_flux_defect_wb = 0.0_f64;
        for (face_id, face) in self.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            let [a, bp, c] = self.mesh().faces[face_id].map(|v| self.positions_m()[v as usize]);
            let raw = cross(super::sub(bp, a), super::sub(c, a));
            let sign = if dot(raw, face.outward_normal) > 0.0 {
                1.0
            } else {
                -1.0
            };
            for &cell in &face.cells {
                let flux = finite(face.area_m2 * dot(b[cell], face.outward_normal))?;
                maximum_recovered_flux_defect_wb = maximum_recovered_flux_defect_wb
                    .max(finite(flux - sign * face_flux_wb[face_id])?.abs());
            }
        }
        if maximum_recovered_flux_defect_wb > controls.flux_tolerance_wb {
            return Err(invalid(
                "recovered cell B does not reproduce its conforming magnetic face flux",
            ));
        }
        let impressed_pairing_j = finite(
            current_load
                .iter()
                .zip(&potential_wb)
                .map(|(j, a)| j * a)
                .sum(),
        )?;
        let boundary_pairing_j = finite(
            residual_a
                .iter()
                .zip(&potential_wb)
                .enumerate()
                .filter(|(e, _)| prescribed[*e].is_some())
                .map(|(_, (r, a))| r * a)
                .sum(),
        )?;
        let constitutive_pairing_j = finite(
            h.iter()
                .zip(&b)
                .enumerate()
                .map(|(cell, (h, b))| self.geometry.geometry.vol_signed[cell].abs() * dot(*h, *b))
                .sum(),
        )?;
        let work_residual_j =
            finite(constitutive_pairing_j - impressed_pairing_j - boundary_pairing_j)?;
        if work_residual_j.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "static current/boundary and constitutive magnetic pairings disagree",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialVectorField {
            problem: self,
            potential_wb,
            face_flux_wb,
            flux_density_t: b,
            field_strength_a_m: h,
            current_density_a_m2: current_a_m2.to_vec(),
            cell_energy_j,
            energy_j,
            recoil_operating_points,
            residual_a,
            prescribed_edges: prescribed.iter().map(Option::is_some).collect(),
            impressed_pairing_j,
            boundary_pairing_j,
            work_residual_j,
            maximum_current_flux_defect_a,
            source_compatibility_defect_a,
            gauge_flux_defect_wb,
            maximum_residual_a,
            maximum_divergence_flux_wb,
            maximum_recovered_flux_defect_wb,
            field_report,
            source_report,
            gauge_report,
        })
    }
}
/// Accepted actual 3D edge potential and H(div)-conforming curl field.
/// Static source/trace pairings are not time-integrated work or Joule heat.
pub struct SpatialVectorField<'m> {
    problem: &'m SpatialVectorMagnetostatics,
    potential_wb: Vec<f64>,
    face_flux_wb: Vec<f64>,
    flux_density_t: Vec<[f64; 3]>,
    field_strength_a_m: Vec<[f64; 3]>,
    current_density_a_m2: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    residual_a: Vec<f64>,
    prescribed_edges: Vec<bool>,
    pub impressed_pairing_j: f64,
    pub boundary_pairing_j: f64,
    pub work_residual_j: f64,
    pub maximum_current_flux_defect_a: f64,
    pub source_compatibility_defect_a: f64,
    pub gauge_flux_defect_wb: f64,
    pub maximum_residual_a: f64,
    pub maximum_divergence_flux_wb: f64,
    pub maximum_recovered_flux_defect_wb: f64,
    pub field_report: Option<SolveReport>,
    pub source_report: Option<SolveReport>,
    pub gauge_report: Option<SolveReport>,
}
impl<'m> SpatialVectorField<'m> {
    pub const fn problem(&self) -> &'m SpatialVectorMagnetostatics {
        self.problem
    }
    pub fn potential_wb(&self) -> &[f64] {
        &self.potential_wb
    }
    pub fn face_flux_wb(&self) -> &[f64] {
        &self.face_flux_wb
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    pub fn field_strength_a_m(&self) -> &[[f64; 3]] {
        &self.field_strength_a_m
    }
    pub fn impressed_current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.current_density_a_m2
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.energy_j
    }
    pub fn recoil_operating_points(&self) -> &[Option<RecoilOperatingPoint>] {
        &self.recoil_operating_points
    }
    /// Original full edge equations: fixed-trace entries are boundary reactions,
    /// every other entry is independently checked against the physical limit.
    pub fn edge_equation_residual_a(&self) -> &[f64] {
        &self.residual_a
    }
    pub fn prescribed_edges(&self) -> &[bool] {
        &self.prescribed_edges
    }
    /// Closed exterior Maxwell traction through the shared spatial integrator.
    /// Both incident cells must be current-free, induced-linear, and bound to
    /// the same complete material snapshot; recoil and current cells require
    /// an exterior collar before force can be sampled.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        placement: &fs_ga::Se3,
        moment_origin_m: [f64; 3],
        limits: super::force::SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<fs_ga::Wrench> {
        let mut reference: Option<&fs_material::state_point::ResolvedMaterialStatePoint> = None;
        super::force::integrate(
            super::force::SurfaceGeometry {
                mesh: self.problem.mesh(),
                positions_m: self.problem.positions_m(),
                faces: &self.problem.geometry.faces,
            },
            enclosed_cells,
            placement,
            moment_origin_m,
            limits,
            |inner, outer| {
                for cell in [inner, outer] {
                    if self.current_density_a_m2[cell].iter().any(|&j| j != 0.0) {
                        return Err(invalid(
                            "Maxwell surface requires a current-free exterior collar",
                        ));
                    }
                    let SpatialMagneticMaterial::Linear(material) = &self.problem.materials[cell]
                    else {
                        return Err(invalid(
                            "Maxwell surface requires induced-linear exterior material",
                        ));
                    };
                    if reference.is_some_and(|previous| previous != material) {
                        return Err(invalid(
                            "Maxwell surface requires one homogeneous original material snapshot",
                        ));
                    }
                    reference = Some(material);
                }
                Ok((self.flux_density_t[outer], self.field_strength_a_m[outer]))
            },
            cancelled,
        )
    }
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if self
            .problem
            .materials
            .iter()
            .any(|m| matches!(m, SpatialMagneticMaterial::Recoil(_)))
        {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
}

fn prepare_cells(
    geometry: &SpatialGeometry,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<CellVectors>> {
    let mut cells = Vec::with_capacity(geometry.mesh.tets.len());
    for (cell, tet) in geometry.mesh.tets.iter().enumerate() {
        poll(cancelled)?;
        let grads = &geometry.geometry.grads[cell];
        let volume = geometry.geometry.vol_signed[cell].abs();
        let mut edges = [0; 6];
        let mut curls = [[0.0; 3]; 6];
        let mut integrals = [[0.0; 3]; 6];
        let mut e = 0;
        for a in 0..4 {
            for b in a + 1..4 {
                let (a, b) = if tet[a] < tet[b] { (a, b) } else { (b, a) };
                let edge = geometry
                    .mesh
                    .edges
                    .binary_search(&[tet[a], tet[b]])
                    .map_err(|_| invalid("missing canonical magnetic edge"))?;
                let curl = cross(grads[a], grads[b]).map(|v| 2.0 * v);
                let integral = std::array::from_fn(|j| volume * (grads[b][j] - grads[a][j]) / 4.0);
                for value in curl.into_iter().chain(integral) {
                    finite(value)?;
                }
                edges[e] = edge;
                curls[e] = curl;
                integrals[e] = integral;
                e += 1;
            }
        }
        cells.push(CellVectors {
            edges,
            curls,
            integrals,
        });
    }
    Ok(cells)
}

/// Borrowed original geometry shared by the linear and nonlinear current owners.
struct VectorTopology<'a> {
    geometry: &'a SpatialGeometry,
    cells: &'a [CellVectors],
}
struct VectorSources {
    prescribed: Vec<Option<f64>>,
    free: Vec<usize>,
    projector: PreparedCoexact1,
    current_load: Vec<f64>,
    maximum_current_flux_defect_a: f64,
    source_compatibility_defect_a: f64,
    source_report: Option<SolveReport>,
}
impl VectorTopology<'_> {
    fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    fn face_edges(&self, face: usize) -> Option<[usize; 3]> {
        let [a, b, c] = *self.mesh().faces.get(face)?;
        let mut result = [0; 3];
        for (i, key) in [[a, b], [a, c], [b, c]].iter().enumerate() {
            result[i] = self.mesh().edges.binary_search(key).ok()?;
        }
        Some(result)
    }

    fn prepare_sources(
        &self,
        current_a_m2: &[[f64; 3]],
        boundary: &VectorMagneticBoundary,
        controls: VectorMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorSources> {
        poll(cancelled)?;
        if current_a_m2.len() != self.cells.len()
            || current_a_m2.iter().flatten().any(|v| !v.is_finite())
            || [
                controls.current_tolerance_a,
                controls.flux_tolerance_wb,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || [controls.field.tolerance, controls.gauge.tolerance]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0 || *v >= 1.0)
        {
            return Err(invalid(
                "complete finite cell current and positive physical/solver allowances required",
            ));
        }
        let ne = self.mesh().edges.len();
        let mut prescribed = vec![None; ne];
        let mut fixed_faces = BTreeSet::new();
        for trace in &boundary.potential {
            poll(cancelled)?;
            if self
                .geometry
                .faces
                .get(trace.face)
                .is_none_or(|f| f.cells.len() != 1)
                || !fixed_faces.insert(trace.face)
                || trace.edge_integrals_wb.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "unique exterior tangential-potential faces with finite edge integrals required",
                ));
            }
            let edges = self
                .face_edges(trace.face)
                .ok_or(invalid("missing potential trace edge"))?;
            for (edge, value) in edges.into_iter().zip(trace.edge_integrals_wb) {
                if prescribed[edge].is_some_and(|old| old != value) {
                    return Err(invalid(
                        "conflicting tangential potential on a shared boundary edge",
                    ));
                }
                prescribed[edge] = Some(value);
            }
        }
        let mut maximum_current_flux_defect_a = 0.0_f64;
        for (i, face) in self.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            let current = if face.cells.len() == 2 {
                std::array::from_fn(|j| {
                    current_a_m2[face.cells[0]][j] - current_a_m2[face.cells[1]][j]
                })
            } else if fixed_faces.contains(&i) {
                continue;
            } else {
                current_a_m2[face.cells[0]]
            };
            let defect = finite(face.area_m2 * dot(current, face.outward_normal))?.abs();
            maximum_current_flux_defect_a = maximum_current_flux_defect_a.max(defect);
        }
        if maximum_current_flux_defect_a > controls.current_tolerance_a {
            return Err(invalid(
                "impressed current violates interface or natural-exterior charge continuity",
            ));
        }
        let mut current_load = vec![0.0; ne];
        for (cell, vectors) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            for local in 0..6 {
                let e = vectors.edges[local];
                current_load[e] =
                    finite(current_load[e] + dot(current_a_m2[cell], vectors.integrals[local]))?;
            }
        }
        let free: Vec<_> = (0..ne).filter(|&e| prescribed[e].is_none()).collect();
        let projector = PreparedCoexact1::new(
            self.mesh(),
            self.positions_m(),
            &self.geometry.geometry,
            &free,
            controls.maximum_gauge_product_terms,
            cancelled,
        )
        .map_err(gauge_error)?;
        let free_current: Vec<_> = free.iter().map(|&e| current_load[e]).collect();
        let riesz: Vec<_> = free_current
            .iter()
            .zip(projector.metric())
            .map(|(j, m)| finite(j / m))
            .collect::<Result<_>>()?;
        let (_, projected_source, source_report) =
            project(&projector, &riesz, controls.gauge, cancelled)?;
        let source_compatibility_defect_a = max_difference(&free_current, &projected_source)?;
        if source_compatibility_defect_a > controls.current_tolerance_a {
            return Err(invalid(
                "impressed current contains a gradient or harmonic component incompatible with the magnetic boundary",
            ));
        }
        Ok(VectorSources {
            prescribed,
            free,
            projector,
            current_load,
            maximum_current_flux_defect_a,
            source_compatibility_defect_a,
            source_report,
        })
    }
}
