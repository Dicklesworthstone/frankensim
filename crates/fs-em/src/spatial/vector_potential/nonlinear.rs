//! Saturating current-carrying 3D Whitney magnetostatics.
//!
//! The original source energy psi(B) gives H and its exact symmetric radial
//! tangent. Newton minimizes integral psi(curl A) - J.A at fixed tangential A.
//! Geometry, exact current loading, complete kernel compatibility and coexact
//! gauge come from the original vector owner. Source ranges are never extended.
//! This is reversible static magnetism, not hysteresis or transient induction.

use super::super::{
    SpatialGeometry, SpatialMagneticError, cross, dot, finite, invalid,
    nonlinear::SpatialNonlinearMaterial, norm, poll, prepare_geometry, sub,
};
use super::{
    CellVectors, VectorMagneticBoundary, VectorMagneticControls, VectorSources, VectorTopology,
    apply, compatible_solve, max_difference, prepare_cells, project,
};
use crate::planar::permanent::MagneticEnergyReference;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::field_storage::{linear_magnetic_permeability, permanent::RecoilOperatingPoint};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

#[derive(Clone, Copy, Debug)]
pub struct NonlinearVectorControls {
    pub field: VectorMagneticControls,
    pub relative_residual_tolerance: f64,
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
}
#[derive(Clone, Debug)]
pub struct NonlinearVectorReport {
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub backtracks: usize,
    pub residual_norm_a: f64,
    pub residual_target_a: f64,
    pub maximum_residual_a: f64,
    pub source_compatibility_defect_a: f64,
    pub maximum_gauge_flux_defect_wb: f64,
    pub maximum_newton_compatibility_defect_a: f64,
}
#[derive(Debug)]
pub enum NonlinearVectorError {
    Spatial(SpatialMagneticError),
    FluxRange {
        cell: usize,
        observed_t: f64,
        maximum_t: f64,
    },
    Nonconvergence(NonlinearVectorReport),
}
impl From<SpatialMagneticError> for NonlinearVectorError {
    fn from(e: SpatialMagneticError) -> Self {
        Self::Spatial(e)
    }
}
impl core::fmt::Display for NonlinearVectorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear vector magnetic field: {self:?}")
    }
}
impl std::error::Error for NonlinearVectorError {}
type Result<T> = core::result::Result<T, NonlinearVectorError>;

pub struct SpatialNonlinearVectorMagnetostatics {
    pub(in crate::spatial) geometry: SpatialGeometry,
    materials: Vec<SpatialNonlinearMaterial>,
    linear_permeability: Vec<Option<f64>>,
    pub(super) cells: Vec<CellVectors>,
    pub(super) curl: Csr,
    pub(super) limits: WeightedAssemblyLimits,
}
pub(super) struct Response {
    pub(super) h: [f64; 3],
    pub(super) energy: f64,
    coenergy: f64,
    pub(super) tangent: [[f64; 3]; 3],
    recoil: Option<RecoilOperatingPoint>,
}
struct Evaluation {
    b: Vec<[f64; 3]>,
    h: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    coenergy_j: f64,
    objective_j: f64,
    residual_a: Vec<f64>,
    tangent: Vec<CellWeight>,
    recoil: Vec<Option<RecoilOperatingPoint>>,
    legendre_defect_j: f64,
}
impl SpatialNonlinearVectorMagnetostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<SpatialNonlinearMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if materials.len() != tetrahedra.len() {
            return Err(invalid("one original magnetic source per tetrahedron required").into());
        }
        let geometry = prepare_geometry(positions_m, tetrahedra, limits, cancelled)?;
        let cells = prepare_cells(&geometry, cancelled)?;
        let mut linear_permeability = Vec::with_capacity(materials.len());
        for material in &materials {
            poll(cancelled)?;
            linear_permeability.push(match material {
                SpatialNonlinearMaterial::Linear(m) => {
                    Some(linear_magnetic_permeability(m).map_err(SpatialMagneticError::Material)?)
                }
                SpatialNonlinearMaterial::Recoil(m) => Some(m.material().permeability_h_m()),
                SpatialNonlinearMaterial::Reversible(m) => {
                    m.evaluate(0.0).map_err(SpatialMagneticError::Material)?;
                    None
                }
            });
        }
        let curl = incidence_to_csr(&geometry.mesh.d1());
        poll(cancelled)?;
        Ok(Self {
            geometry,
            materials,
            linear_permeability,
            cells,
            curl,
            limits,
        })
    }
    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn materials(&self) -> &[SpatialNonlinearMaterial] {
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
        VectorTopology {
            geometry: &self.geometry,
            cells: &self.cells,
        }
        .face_edges(face)
    }
    pub(super) fn response(&self, cell: usize, b: [f64; 3]) -> Result<Response> {
        let magnitude = finite(norm(b))?;
        let (h, energy, coenergy, secant, differential, recoil) = match &self.materials[cell] {
            SpatialNonlinearMaterial::Reversible(law) => {
                if magnitude > law.maximum_flux_density_t() {
                    return Err(NonlinearVectorError::FluxRange {
                        cell,
                        observed_t: magnitude,
                        maximum_t: law.maximum_flux_density_t(),
                    });
                }
                let value = law
                    .evaluate(magnitude)
                    .map_err(SpatialMagneticError::Material)?;
                let differential = value.differential_reluctivity_m_h;
                let secant = if magnitude == 0.0 {
                    differential
                } else {
                    finite(value.field_strength_a_m / magnitude)?
                };
                (
                    b.map(|v| v * secant),
                    value.energy_density_j_m3,
                    value.coenergy_density_j_m3,
                    secant,
                    differential,
                    None,
                )
            }
            SpatialNonlinearMaterial::Linear(_) | SpatialNonlinearMaterial::Recoil(_) => {
                let mu = self.linear_permeability[cell].expect("prepared linear source");
                let br = if let SpatialNonlinearMaterial::Recoil(m) = &self.materials[cell] {
                    m.remanence_t()
                } else {
                    [0.0; 3]
                };
                let h = std::array::from_fn(|j| (b[j] - br[j]) / mu);
                for v in h {
                    finite(v)?;
                }
                let energy = finite(0.5 * mu * dot(h, h))?;
                let coenergy = finite(dot(b, h) - energy)?;
                let operating = if let SpatialNonlinearMaterial::Recoil(m) = &self.materials[cell] {
                    let parallel = finite(dot(h, m.direction()))?;
                    let transverse =
                        norm(std::array::from_fn(|j| h[j] - parallel * m.direction()[j]));
                    Some(
                        m.material()
                            .admit_operating_point(parallel, transverse)
                            .map_err(|error| SpatialMagneticError::Recoil { cell, error })?,
                    )
                } else {
                    None
                };
                (
                    h,
                    energy,
                    coenergy,
                    finite(mu.recip())?,
                    finite(mu.recip())?,
                    operating,
                )
            }
        };
        if !secant.is_finite() || !differential.is_finite() || secant <= 0.0 || differential <= 0.0
        {
            return Err(invalid("positive representable differential reluctivity required").into());
        }
        for v in h {
            finite(v)?;
        }
        let u = if magnitude == 0.0 {
            [0.0; 3]
        } else {
            b.map(|v| v / magnitude)
        };
        let mut tangent = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                let v = finite(
                    ((differential - secant) * u[i])
                        .mul_add(u[j], if i == j { secant } else { 0.0 }),
                )?;
                tangent[i][j] = v;
                tangent[j][i] = v;
            }
        }
        Ok(Response {
            h,
            energy,
            coenergy,
            tangent,
            recoil,
        })
    }
    fn evaluate(
        &self,
        a: &[f64],
        source: &VectorSources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        let mut out = Evaluation {
            b: Vec::with_capacity(self.cells.len()),
            h: Vec::with_capacity(self.cells.len()),
            cell_energy_j: Vec::with_capacity(self.cells.len()),
            energy_j: 0.0,
            coenergy_j: 0.0,
            objective_j: 0.0,
            residual_a: source.current_load.iter().map(|v| -v).collect(),
            tangent: Vec::with_capacity(self.cells.len()),
            recoil: Vec::with_capacity(self.cells.len()),
            legendre_defect_j: 0.0,
        };
        for (cell, vectors) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let mut b = [0.0; 3];
            for local in 0..6 {
                for j in 0..3 {
                    b[j] = finite(a[vectors.edges[local]].mul_add(vectors.curls[local][j], b[j]))?;
                }
            }
            let response = self.response(cell, b)?;
            let volume = self.geometry.geometry.vol_signed[cell].abs();
            let energy = finite(volume * response.energy)?;
            out.energy_j = finite(out.energy_j + energy)?;
            out.coenergy_j = finite(out.coenergy_j + volume * response.coenergy)?;
            out.legendre_defect_j = finite(
                out.legendre_defect_j
                    + volume * (response.energy + response.coenergy - dot(b, response.h)).abs(),
            )?;
            for local in 0..6 {
                let e = vectors.edges[local];
                out.residual_a[e] =
                    finite((volume * dot(vectors.curls[local], response.h)) + out.residual_a[e])?;
            }
            out.b.push(b);
            out.h.push(response.h);
            out.cell_energy_j.push(energy);
            out.tangent.push(CellWeight::Tensor(response.tangent));
            out.recoil.push(response.recoil);
        }
        let current_pairing = finite(a.iter().zip(&source.current_load).map(|(a, j)| a * j).sum())?;
        out.objective_j = finite(out.energy_j - current_pairing)?;
        poll(cancelled)?;
        Ok(out)
    }

    /// Solve at held current and tangential A from a complete source-admitted
    /// initial edge cochain. Initial and Newton gauges preserve physical curl.
    /// Every original free-edge equation must pass, including all gauge modes.
    pub fn solve(
        &self,
        current_a_m2: &[[f64; 3]],
        boundary: &VectorMagneticBoundary,
        initial_potential_wb: &[f64],
        controls: NonlinearVectorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearVectorField<'_>> {
        poll(cancelled)?;
        if initial_potential_wb.len() != self.mesh().edges.len()
            || initial_potential_wb.iter().any(|v| !v.is_finite())
            || !controls.relative_residual_tolerance.is_finite()
            || controls.relative_residual_tolerance <= 0.0
            || controls.relative_residual_tolerance >= 1.0
        {
            return Err(invalid(
                "complete finite initial A and positive bounded Newton controls required",
            )
            .into());
        }
        let sources = VectorTopology {
            geometry: &self.geometry,
            cells: &self.cells,
        }
        .prepare_sources(current_a_m2, boundary, controls.field, cancelled)?;
        let mut a = initial_potential_wb.to_vec();
        for (e, value) in sources.prescribed.iter().enumerate() {
            if let Some(value) = value {
                a[e] = *value;
            }
        }
        let gauge = |values: &[f64], cancelled: &mut _| -> Result<(Vec<f64>, f64, usize)> {
            let (projected, _, report) =
                project(&sources.projector, values, controls.field.gauge, cancelled)?;
            let before = sources
                .projector
                .rhs(values, cancelled)
                .map_err(super::gauge_error)?;
            let after = sources
                .projector
                .rhs(&projected, cancelled)
                .map_err(super::gauge_error)?;
            let defect = max_difference(&before, &after)?;
            if defect > controls.field.flux_tolerance_wb {
                return Err(invalid("nonlinear coexact gauge changed physical flux").into());
            }
            Ok((projected, defect, report.map_or(0, |r| r.iters)))
        };
        let raw: Vec<_> = sources.free.iter().map(|&e| a[e]).collect();
        let (initial_gauge, initial_defect, gauge_iterations) = gauge(&raw, cancelled)?;
        for (&e, &v) in sources.free.iter().zip(&initial_gauge) {
            a[e] = v;
        }
        let measure = |r: &[f64]| -> Result<(f64, f64)> {
            let norm = finite(sources.free.iter().fold(0.0_f64, |n, &e| n.hypot(r[e])))?;
            let max = sources.free.iter().fold(0.0_f64, |n, &e| n.max(r[e].abs()));
            Ok((norm, max))
        };
        let mut value = self.evaluate(&a, &sources, cancelled)?;
        let (norm, max) = measure(&value.residual_a)?;
        let target = controls
            .field
            .current_tolerance_a
            .max(finite(controls.relative_residual_tolerance * norm)?);
        let mut report = NonlinearVectorReport {
            newton_iterations: 0,
            linear_iterations: gauge_iterations,
            backtracks: 0,
            residual_norm_a: norm,
            residual_target_a: target,
            maximum_residual_a: max,
            source_compatibility_defect_a: sources.source_compatibility_defect_a,
            maximum_gauge_flux_defect_wb: initial_defect,
            maximum_newton_compatibility_defect_a: 0.0,
        };
        let mut index = vec![usize::MAX; a.len()];
        for (i, &e) in sources.free.iter().enumerate() {
            index[e] = i;
        }
        loop {
            poll(cancelled)?;
            if report.residual_norm_a <= target
                && report.maximum_residual_a <= controls.field.current_tolerance_a
            {
                break;
            }
            if report.newton_iterations >= controls.maximum_newton_iterations {
                return Err(NonlinearVectorError::Nonconvergence(report));
            }
            let mass = weighted_mass_matrix(
                self.mesh(),
                &self.geometry.geometry,
                2,
                &value.tangent,
                Dims([-1, -1, 2, 0, 2, 0]),
                self.limits,
                cancelled,
            )
            .map_err(|e| {
                if e == WeightedError::Cancelled {
                    SpatialMagneticError::Cancelled
                } else {
                    SpatialMagneticError::Assembly(e)
                }
            })?;
            let full = stiffness(&self.curl, mass.matrix());
            let mut rows = vec![0];
            let mut cols = Vec::new();
            let mut vals = Vec::new();
            for &e in &sources.free {
                poll(cancelled)?;
                let (c, v) = full.row(e);
                for (&c, &v) in c.iter().zip(v) {
                    if index[c] != usize::MAX {
                        cols.push(index[c]);
                        vals.push(v);
                    }
                }
                rows.push(cols.len());
            }
            let matrix = Csr::from_parts(sources.free.len(), sources.free.len(), rows, cols, vals);
            let rhs: Vec<_> = sources.free.iter().map(|&e| -value.residual_a[e]).collect();
            // Newton lives in the admitted coexact tangent space. Project its
            // assembled dual residual before a relative Krylov solve: floating
            // curl cancellation must not become an incompatible kernel RHS as
            // the physical residual becomes small. The original current and
            // full nonlinear residual remain unchanged and independently gated.
            let riesz: Vec<_> = rhs
                .iter()
                .zip(sources.projector.metric())
                .map(|(r, m)| finite(r / m))
                .collect::<core::result::Result<_, _>>()?;
            let (_, projected_rhs, compatibility_report) =
                project(&sources.projector, &riesz, controls.field.gauge, cancelled)?;
            let compatibility_defect = max_difference(&rhs, &projected_rhs)?;
            if compatibility_defect > controls.field.current_tolerance_a {
                return Err(
                    invalid("Newton residual leaves the admitted coexact tangent space").into(),
                );
            }
            report.maximum_newton_compatibility_defect_a = report
                .maximum_newton_compatibility_defect_a
                .max(compatibility_defect);
            report.linear_iterations += compatibility_report.map_or(0, |r| r.iters);
            let (raw, linear) =
                compatible_solve(&matrix, &projected_rhs, controls.field.field, cancelled)?;
            report.linear_iterations += linear.map_or(0, |r| r.iters);
            let (direction, defect, iters) = gauge(&raw, cancelled)?;
            report.linear_iterations += iters;
            report.maximum_gauge_flux_defect_wb = report.maximum_gauge_flux_defect_wb.max(defect);
            let slope = finite(
                sources
                    .free
                    .iter()
                    .zip(&direction)
                    .map(|(&e, d)| value.residual_a[e] * d)
                    .sum(),
            )?;
            if slope >= 0.0 {
                return Err(NonlinearVectorError::Nonconvergence(report));
            }
            let mut factor: f64 = 1.0;
            let mut accepted = None;
            for backtrack in 0..=controls.maximum_backtracks {
                poll(cancelled)?;
                let mut trial = a.clone();
                for (&e, &d) in sources.free.iter().zip(&direction) {
                    trial[e] = finite(factor.mul_add(d, a[e]))?;
                }
                match self.evaluate(&trial, &sources, cancelled) {
                    Ok(next)
                        if next.objective_j
                            <= finite(value.objective_j + 1e-4 * factor * slope)? =>
                    {
                        accepted = Some((trial, next));
                        break;
                    }
                    Ok(_)
                    | Err(NonlinearVectorError::FluxRange { .. })
                    | Err(NonlinearVectorError::Spatial(SpatialMagneticError::Recoil { .. })) => {}
                    Err(e) => return Err(e),
                }
                if backtrack < controls.maximum_backtracks {
                    factor *= 0.5;
                    report.backtracks += 1;
                }
            }
            let Some((next, next_value)) = accepted else {
                return Err(NonlinearVectorError::Nonconvergence(report));
            };
            a = next;
            value = next_value;
            report.newton_iterations += 1;
            (report.residual_norm_a, report.maximum_residual_a) = measure(&value.residual_a)?;
        }
        let face_flux_wb = apply(&self.curl, &a, cancelled)?;
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
        let mut maximum_recovered_flux_defect_wb = 0.0_f64;
        for (i, face) in self.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            let [a, b, c] = self.mesh().faces[i].map(|v| self.positions_m()[v as usize]);
            let sign = if dot(cross(sub(b, a), sub(c, a)), face.outward_normal) > 0.0 {
                1.0
            } else {
                -1.0
            };
            for &cell in &face.cells {
                let flux = finite(face.area_m2 * dot(value.b[cell], face.outward_normal))?;
                maximum_recovered_flux_defect_wb = maximum_recovered_flux_defect_wb
                    .max(finite(flux - sign * face_flux_wb[i])?.abs());
            }
        }
        if maximum_divergence_flux_wb > controls.field.flux_tolerance_wb
            || maximum_recovered_flux_defect_wb > controls.field.flux_tolerance_wb
        {
            return Err(invalid("nonlinear B violates original conforming div-curl flux").into());
        }
        let impressed_pairing_j = finite(
            a.iter()
                .zip(&sources.current_load)
                .map(|(a, j)| a * j)
                .sum(),
        )?;
        // Only actual prescribed-edge reactions are physical boundary supply.
        let boundary_pairing_j = finite(
            a.iter()
                .zip(&value.residual_a)
                .enumerate()
                .filter(|(e, _)| sources.prescribed[*e].is_some())
                .map(|(_, (a, r))| a * r)
                .sum(),
        )?;
        let work_residual_j =
            finite(value.energy_j + value.coenergy_j - impressed_pairing_j - boundary_pairing_j)?;
        if value.legendre_defect_j > controls.field.energy_tolerance_j
            || work_residual_j.abs() > controls.field.energy_tolerance_j
        {
            return Err(invalid(
                "independent nonlinear energy/coenergy/current/trace pairing failed",
            )
            .into());
        }
        poll(cancelled)?;
        Ok(SpatialNonlinearVectorField {
            problem: self,
            potential_wb: a,
            face_flux_wb,
            flux_density_t: value.b,
            field_strength_a_m: value.h,
            current_density_a_m2: current_a_m2.to_vec(),
            prescribed_edges: sources.prescribed.iter().map(Option::is_some).collect(),
            cell_energy_j: value.cell_energy_j,
            energy_j: value.energy_j,
            coenergy_j: value.coenergy_j,
            recoil_operating_points: value.recoil,
            edge_residual_a: value.residual_a,
            impressed_pairing_j,
            boundary_pairing_j,
            work_residual_j,
            legendre_defect_j: value.legendre_defect_j,
            maximum_divergence_flux_wb,
            maximum_recovered_flux_defect_wb,
            report,
        })
    }
}
pub struct SpatialNonlinearVectorField<'m> {
    problem: &'m SpatialNonlinearVectorMagnetostatics,
    potential_wb: Vec<f64>,
    face_flux_wb: Vec<f64>,
    flux_density_t: Vec<[f64; 3]>,
    field_strength_a_m: Vec<[f64; 3]>,
    current_density_a_m2: Vec<[f64; 3]>,
    prescribed_edges: Vec<bool>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    coenergy_j: f64,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    edge_residual_a: Vec<f64>,
    pub impressed_pairing_j: f64,
    pub boundary_pairing_j: f64,
    pub work_residual_j: f64,
    pub legendre_defect_j: f64,
    pub maximum_divergence_flux_wb: f64,
    pub maximum_recovered_flux_defect_wb: f64,
    pub report: NonlinearVectorReport,
}
impl<'m> SpatialNonlinearVectorField<'m> {
    pub fn problem(&self) -> &'m SpatialNonlinearVectorMagnetostatics {
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
    /// Actual tangential-A trace edges, distinct from numerical gauge choices.
    pub fn prescribed_edges(&self) -> &[bool] {
        &self.prescribed_edges
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.energy_j
    }
    pub fn coenergy_j(&self) -> f64 {
        self.coenergy_j
    }
    pub fn recoil_operating_points(&self) -> &[Option<RecoilOperatingPoint>] {
        &self.recoil_operating_points
    }
    pub fn edge_equation_residual_a(&self) -> &[f64] {
        &self.edge_residual_a
    }
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if self
            .problem
            .materials
            .iter()
            .any(|m| matches!(m, SpatialNonlinearMaterial::Recoil(_)))
        {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
    /// Static Maxwell wrench through the original source-identical current-free
    /// linear collar. Saturating/recoil material may be enclosed, never sampled.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        placement: &fs_ga::Se3,
        moment_origin_m: [f64; 3],
        limits: super::super::force::SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<fs_ga::Wrench, SpatialMagneticError> {
        let mut reference = None;
        super::super::force::integrate(
            super::super::force::SurfaceGeometry {
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
                    if self.current_density_a_m2[cell].iter().any(|&v| v != 0.0) {
                        return Err(invalid(
                            "Maxwell surface requires a current-free exterior collar",
                        ));
                    }
                    let SpatialNonlinearMaterial::Linear(material) = &self.problem.materials[cell]
                    else {
                        return Err(invalid(
                            "Maxwell surface requires original linear exterior material",
                        ));
                    };
                    let source = *reference.get_or_insert(material);
                    if source != material {
                        return Err(invalid("Maxwell collar source snapshots differ"));
                    }
                }
                Ok((self.flux_density_t[outer], self.field_strength_a_m[outer]))
            },
            cancelled,
        )
    }
}
