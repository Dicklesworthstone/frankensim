//! Saturating current-free 3D magnetostatics with the original sourced H(B).
//!
//! The scalar potential gives H=-grad(phi). The strictly increasing source law
//! is inverted segment by segment to obtain B(H), without a second fit. Its
//! integrated coenergy is the convex Newton objective; its exact differential
//! permeability supplies the original tensor Whitney Hodge. Energy-decreasing
//! backtracking never extends a source range. Linear regions and fixed oriented
//! recoil magnets share this weak flux equation and its boundary/gauge owner.
//!
//! This is a fixed-state, single-valued scalar-potential problem: no free
//! currents, hysteresis, induction, parameter evolution or mechanical work is
//! inferred. Finite exterior boundaries and raw P1 interface flux jumps retain
//! the original spatial formulation's meaning. Magnetic energy and coenergy
//! differ in saturation; neither is replaced by B.H/2.

use super::{
    ScalarBoundary, SpatialGeometry, SpatialMagneticBoundary, SpatialMagneticControls,
    SpatialMagneticError, SpatialRecoilMaterial, dot, finite, invalid, norm, poll,
    prepare_geometry, prepare_scalar_boundary, solve_scalar_system, sub,
};
use crate::planar::permanent::MagneticEnergyReference;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, linear_magnetic_permeability,
        nonlinear::ReversibleMagneticMaterial, permanent::RecoilOperatingPoint,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

/// One original source at fixed external conditions per original tetrahedron.
#[derive(Clone, Debug, PartialEq)]
pub enum SpatialNonlinearMaterial {
    Linear(ResolvedMaterialStatePoint),
    Reversible(ReversibleMagneticMaterial),
    Recoil(SpatialRecoilMaterial),
}
impl SpatialNonlinearMaterial {
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Linear(material) => material,
            Self::Reversible(law) => law.material(),
            Self::Recoil(law) => law.material().material(),
        }
    }
}

/// Bounded Newton/CG work with dimensional flux and energy acceptance.
#[derive(Clone, Copy, Debug)]
pub struct SpatialNonlinearControls {
    /// Original CG, maximum weak flux residual and energy allowances.
    pub field: SpatialMagneticControls,
    /// Euclidean residual target is max(flux_tolerance, relative * initial norm).
    /// The original maximum physical weak residual must ALSO pass, including
    /// every artificial gauge vertex.
    pub relative_residual_tolerance: f64,
    /// Zero checks an initial field without admitting a Newton update.
    pub maximum_newton_iterations: usize,
    /// Maximum halvings after each full Newton trial.
    pub maximum_backtracks: usize,
}
#[derive(Clone, Debug)]
pub struct SpatialNonlinearReport {
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub backtracks: usize,
    pub residual_norm_wb: f64,
    pub residual_target_wb: f64,
    pub maximum_weak_residual_wb: f64,
}
#[derive(Debug)]
pub enum SpatialNonlinearError {
    Spatial(SpatialMagneticError),
    FieldRange {
        cell: usize,
        observed_field_a_m: f64,
        maximum_field_a_m: f64,
    },
    Nonconvergence(SpatialNonlinearReport),
}
impl From<SpatialMagneticError> for SpatialNonlinearError {
    fn from(e: SpatialMagneticError) -> Self {
        Self::Spatial(e)
    }
}
impl core::fmt::Display for SpatialNonlinearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear spatial magnetic field: {self:?}")
    }
}
impl std::error::Error for SpatialNonlinearError {}
type Result<T> = core::result::Result<T, SpatialNonlinearError>;

/// Immutable original geometry and physical laws; Newton rebuilds only their
/// changing material tangent Hodge, with one cached topological incidence.
pub struct SpatialNonlinearMagnetostatics {
    pub(super) geometry: SpatialGeometry,
    pub(super) materials: Vec<SpatialNonlinearMaterial>,
    linear_permeability_h_m: Vec<Option<f64>>,
    incidence: Csr,
    limits: WeightedAssemblyLimits,
}
struct Response {
    b: [f64; 3],
    energy: f64,
    coenergy: f64,
    tangent: [[f64; 3]; 3],
    recoil: Option<RecoilOperatingPoint>,
}
struct Evaluation {
    h: Vec<[f64; 3]>,
    b: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    cell_coenergy_j: Vec<f64>,
    energy_j: f64,
    coenergy_j: f64,
    objective_j: f64,
    legendre_residual_j: f64,
    residual_wb: Vec<f64>,
    tangent: Vec<CellWeight>,
    recoil: Vec<Option<RecoilOperatingPoint>>,
}

impl SpatialNonlinearMagnetostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<SpatialNonlinearMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if materials.len() != tetrahedra.len() {
            return Err(
                invalid("one nonlinear magnetic source per original tetrahedron required").into(),
            );
        }
        let geometry = prepare_geometry(positions_m, tetrahedra, limits, cancelled)?;
        let mut linear_permeability_h_m = Vec::with_capacity(materials.len());
        for material in &materials {
            poll(cancelled)?;
            let mu = match material {
                SpatialNonlinearMaterial::Linear(source) => Some(
                    linear_magnetic_permeability(source).map_err(SpatialMagneticError::Material)?,
                ),
                SpatialNonlinearMaterial::Recoil(law) => Some(law.material().permeability_h_m()),
                SpatialNonlinearMaterial::Reversible(law) => {
                    law.evaluate(0.0).map_err(SpatialMagneticError::Material)?;
                    None
                }
            };
            linear_permeability_h_m.push(mu);
        }
        let incidence = incidence_to_csr(&geometry.mesh.d0());
        poll(cancelled)?;
        Ok(Self {
            geometry,
            materials,
            linear_permeability_h_m,
            incidence,
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
    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.geometry
            .faces
            .get(face)
            .filter(|f| f.cells.len() == 1)
            .map(|f| f.outward_normal)
    }

    fn response(&self, cell: usize, h: [f64; 3]) -> Result<Response> {
        let magnitude = finite(norm(h))?;
        let (b, energy, coenergy, secant, differential, recoil) = match &self.materials[cell] {
            SpatialNonlinearMaterial::Reversible(law) => {
                if magnitude > law.maximum_field_strength_a_m() {
                    return Err(SpatialNonlinearError::FieldRange {
                        cell,
                        observed_field_a_m: magnitude,
                        maximum_field_a_m: law.maximum_field_strength_a_m(),
                    });
                }
                let flux = law
                    .flux_density_at_field_strength(magnitude)
                    .map_err(SpatialMagneticError::Material)?;
                let response = law.evaluate(flux).map_err(SpatialMagneticError::Material)?;
                let differential = finite(response.differential_reluctivity_m_h.recip())?;
                let secant = if magnitude == 0.0 {
                    differential
                } else {
                    finite(flux / magnitude)?
                };
                let b = h.map(|v| v * secant);
                (
                    b,
                    response.energy_density_j_m3,
                    response.coenergy_density_j_m3,
                    secant,
                    differential,
                    None,
                )
            }
            SpatialNonlinearMaterial::Linear(_) | SpatialNonlinearMaterial::Recoil(_) => {
                let mu = self.linear_permeability_h_m[cell].expect("prepared linear/recoil source");
                let energy = finite((0.5 * mu * magnitude) * magnitude)?;
                let (remanence, recoil) =
                    if let SpatialNonlinearMaterial::Recoil(law) = &self.materials[cell] {
                        let parallel = finite(dot(h, law.direction()))?;
                        let transverse = norm(std::array::from_fn(|j| {
                            h[j] - parallel * law.direction()[j]
                        }));
                        let operating = law
                            .material()
                            .admit_operating_point(parallel, transverse)
                            .map_err(|error| SpatialMagneticError::Recoil { cell, error })?;
                        (law.remanence_t(), Some(operating))
                    } else {
                        ([0.0; 3], None)
                    };
                let b = std::array::from_fn(|j| mu.mul_add(h[j], remanence[j]));
                let coenergy = finite(energy + dot(remanence, h))?;
                (b, energy, coenergy, mu, mu, recoil)
            }
        };
        for value in b {
            finite(value)?;
        }
        if secant <= 0.0 || differential <= 0.0 || magnitude > 0.0 && energy <= 0.0 {
            return Err(invalid(
                "positive differential magnetic response and representable energy required",
            )
            .into());
        }
        let direction = if magnitude == 0.0 {
            [0.0; 3]
        } else {
            h.map(|v| v / magnitude)
        };
        let mut tangent = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                let value = finite(
                    ((differential - secant) * direction[i])
                        .mul_add(direction[j], if i == j { secant } else { 0.0 }),
                )?;
                tangent[i][j] = value;
                tangent[j][i] = value;
            }
        }
        Ok(Response {
            b,
            energy,
            coenergy,
            tangent,
            recoil,
        })
    }

    fn evaluate(
        &self,
        potential: &[f64],
        boundary: &ScalarBoundary,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        let mut out = Evaluation {
            h: Vec::with_capacity(self.materials.len()),
            b: Vec::with_capacity(self.materials.len()),
            cell_energy_j: Vec::with_capacity(self.materials.len()),
            cell_coenergy_j: Vec::with_capacity(self.materials.len()),
            energy_j: 0.0,
            coenergy_j: 0.0,
            objective_j: 0.0,
            legendre_residual_j: 0.0,
            residual_wb: boundary.flux_load_wb.clone(),
            tangent: Vec::with_capacity(self.materials.len()),
            recoil: Vec::with_capacity(self.materials.len()),
        };
        for (cell, tet) in self.geometry.mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            let gradients = &self.geometry.geometry.grads[cell];
            let volume = self.geometry.geometry.vol_signed[cell].abs();
            let mut h = [0.0; 3];
            for local in 1..4 {
                let change = finite(potential[tet[local] as usize] - potential[tet[0] as usize])?;
                for j in 0..3 {
                    h[j] = finite((-change).mul_add(gradients[local][j], h[j]))?;
                }
            }
            let response = self.response(cell, h)?;
            let energy = finite(volume * response.energy)?;
            let coenergy = finite(volume * response.coenergy)?;
            let pairing = finite(volume * dot(h, response.b))?;
            let defect = finite(pairing - energy - coenergy)?;
            out.legendre_residual_j = finite(out.legendre_residual_j + defect.abs())?;
            out.energy_j = finite(out.energy_j + energy)?;
            out.coenergy_j = finite(out.coenergy_j + coenergy)?;
            out.cell_energy_j.push(energy);
            out.cell_coenergy_j.push(coenergy);
            out.h.push(h);
            out.b.push(response.b);
            out.tangent.push(CellWeight::Tensor(response.tangent));
            out.recoil.push(response.recoil);
            let mut first = 0.0;
            for local in 1..4 {
                let contribution = finite(-volume * dot(gradients[local], response.b))?;
                let v = tet[local] as usize;
                out.residual_wb[v] = finite(out.residual_wb[v] + contribution)?;
                first = finite(first - contribution)?;
            }
            let v = tet[0] as usize;
            out.residual_wb[v] = finite(out.residual_wb[v] + first)?;
        }
        let mut boundary_work = 0.0;
        for (&p, &flux) in potential.iter().zip(&boundary.flux_load_wb) {
            boundary_work = finite(p.mul_add(flux, boundary_work))?;
        }
        out.objective_j = finite(out.coenergy_j + boundary_work)?;
        poll(cancelled)?;
        Ok(out)
    }

    /// Solve from a complete finite nodal initial potential [A]. Physical
    /// Dirichlet traces take precedence; an automatic gauge is shifted as a
    /// constant on its component, preserving the supplied initial H.
    /// The initial field must already lie inside every source envelope.
    /// No accepted field is returned on source, line-search or work failure.
    pub fn solve(
        &self,
        boundary: &SpatialMagneticBoundary,
        initial_potential_a: &[f64],
        controls: SpatialNonlinearControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearField<'_>> {
        poll(cancelled)?;
        let physical = controls.field;
        let n = self.geometry.positions_m.len();
        if initial_potential_a.len() != n
            || initial_potential_a.iter().any(|v| !v.is_finite())
            || !physical.flux_tolerance_wb.is_finite()
            || physical.flux_tolerance_wb <= 0.0
            || !physical.energy_tolerance_j.is_finite()
            || physical.energy_tolerance_j <= 0.0
            || !physical.cg.tolerance.is_finite()
            || !(0.0..1.0).contains(&physical.cg.tolerance)
            || physical.cg.tolerance == 0.0
            || !controls.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_residual_tolerance)
            || controls.relative_residual_tolerance == 0.0
        {
            return Err(invalid(
                "complete finite initial potential and positive bounded solver controls required",
            )
            .into());
        }
        let boundary = prepare_scalar_boundary(
            &self.geometry.mesh,
            &self.geometry.faces,
            &self.geometry.components,
            boundary,
            physical.flux_tolerance_wb,
            cancelled,
        )?;
        let mut shifts = boundary.references.clone();
        for &gauge in &boundary.gauges {
            shifts[gauge] = initial_potential_a[gauge];
        }
        let mut potential = Vec::with_capacity(n);
        for (v, &initial) in initial_potential_a.iter().enumerate() {
            potential.push(match boundary.prescribed[v] {
                Some(value) => value,
                None => finite(initial - shifts[self.geometry.components[v]])?,
            });
        }
        let free: Vec<_> = (0..n)
            .filter(|&v| boundary.prescribed[v].is_none())
            .collect();
        let zero_trace: Vec<_> = boundary.prescribed.iter().map(|v| v.map(|_| 0.0)).collect();
        let measure = |values: &[f64]| -> Result<(f64, f64)> {
            let norm = finite(
                free.iter()
                    .fold(0.0_f64, |length, &v| length.hypot(values[v])),
            )?;
            let maximum = (0..n)
                .filter(|&v| !boundary.physical_prescribed[v])
                .fold(0.0_f64, |largest, v| largest.max(values[v].abs()));
            Ok((norm, maximum))
        };
        let mut evaluation = self.evaluate(&potential, &boundary, cancelled)?;
        let (initial_norm, initial_max) = measure(&evaluation.residual_wb)?;
        let target = physical
            .flux_tolerance_wb
            .max(finite(controls.relative_residual_tolerance * initial_norm)?);
        let mut report = SpatialNonlinearReport {
            newton_iterations: 0,
            linear_iterations: 0,
            backtracks: 0,
            residual_norm_wb: initial_norm,
            residual_target_wb: target,
            maximum_weak_residual_wb: initial_max,
        };
        loop {
            poll(cancelled)?;
            if report.residual_norm_wb <= target
                && report.maximum_weak_residual_wb <= physical.flux_tolerance_wb
            {
                break;
            }
            if report.newton_iterations >= controls.maximum_newton_iterations {
                return Err(SpatialNonlinearError::Nonconvergence(report));
            }
            let mass = weighted_mass_matrix(
                &self.geometry.mesh,
                &self.geometry.geometry,
                1,
                &evaluation.tangent,
                MAGNETIC_PERMEABILITY_DIMS,
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
            let tangent = stiffness(&self.incidence, mass.matrix());
            let rhs: Vec<_> = evaluation.residual_wb.iter().map(|r| -*r).collect();
            let (direction, linear_report) =
                solve_scalar_system(&tangent, &rhs, &zero_trace, physical.cg, cancelled)?;
            if let Some(linear) = linear_report {
                report.linear_iterations += linear.iters;
            }
            let slope = finite(
                free.iter()
                    .map(|&v| evaluation.residual_wb[v] * direction[v])
                    .sum(),
            )?;
            if slope >= 0.0 {
                return Err(SpatialNonlinearError::Nonconvergence(report));
            }
            let mut factor: f64 = 1.0;
            let mut accepted = None;
            for backtrack in 0..=controls.maximum_backtracks {
                poll(cancelled)?;
                let mut trial = potential.clone();
                for &v in &free {
                    trial[v] = finite(factor.mul_add(direction[v], potential[v]))?;
                }
                let candidate = self.evaluate(&trial, &boundary, cancelled);
                match candidate {
                    Ok(value)
                        if value.objective_j
                            <= finite((1e-4 * factor * slope) + evaluation.objective_j)? =>
                    {
                        accepted = Some((trial, value));
                        break;
                    }
                    Ok(_)
                    | Err(SpatialNonlinearError::FieldRange { .. })
                    | Err(SpatialNonlinearError::Spatial(SpatialMagneticError::Recoil {
                        ..
                    })) => {}
                    Err(e) => return Err(e),
                }
                if backtrack < controls.maximum_backtracks {
                    factor *= 0.5;
                    report.backtracks += 1;
                }
            }
            let Some((next, next_evaluation)) = accepted else {
                return Err(SpatialNonlinearError::Nonconvergence(report));
            };
            potential = next;
            evaluation = next_evaluation;
            report.newton_iterations += 1;
            (report.residual_norm_wb, report.maximum_weak_residual_wb) =
                measure(&evaluation.residual_wb)?;
        }
        let nodal_boundary_flux_wb: Vec<_> = boundary
            .flux_load_wb
            .iter()
            .zip(&evaluation.residual_wb)
            .map(|(&g, &r)| finite(g - r))
            .collect::<core::result::Result<_, _>>()?;
        let mut trace_work = 0.0;
        // Only physical Dirichlet nodes supply unknown reactions. Elsewhere
        // use the authored flux: an interior residual cannot supply boundary work.
        for (v, &p) in potential.iter().enumerate() {
            let flux = if boundary.physical_prescribed[v] {
                nodal_boundary_flux_wb[v]
            } else {
                boundary.flux_load_wb[v]
            };
            trace_work = finite(p.mul_add(flux, trace_work))?;
        }
        let energy_pairing_residual_j =
            finite(evaluation.energy_j + evaluation.coenergy_j + trace_work)?;
        if evaluation.legendre_residual_j > physical.energy_tolerance_j
            || energy_pairing_residual_j.abs() > physical.energy_tolerance_j
        {
            return Err(invalid(
                "nonlinear source energy/coenergy and independent weak field work disagree",
            )
            .into());
        }
        let mut interface_flux_jumps_wb = Vec::new();
        for (face_id, face) in self.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            if face.cells.len() == 2 {
                interface_flux_jumps_wb.push((
                    face_id,
                    finite(
                        face.area_m2
                            * dot(
                                sub(evaluation.b[face.cells[0]], evaluation.b[face.cells[1]]),
                                face.outward_normal,
                            ),
                    )?,
                ));
            }
        }
        let potential_a = potential
            .iter()
            .enumerate()
            .map(|(v, &p)| finite(p + boundary.references[self.geometry.components[v]]))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        poll(cancelled)?;
        Ok(SpatialNonlinearField {
            problem: self,
            potential_a,
            relative_potential_a: potential,
            field_strength_a_m: evaluation.h,
            flux_density_t: evaluation.b,
            cell_energy_j: evaluation.cell_energy_j,
            cell_coenergy_j: evaluation.cell_coenergy_j,
            energy_j: evaluation.energy_j,
            coenergy_j: evaluation.coenergy_j,
            recoil_operating_points: evaluation.recoil,
            nodal_boundary_flux_wb,
            interface_flux_jumps_wb,
            gauge_vertices: boundary.gauges,
            legendre_residual_j: evaluation.legendre_residual_j,
            energy_pairing_residual_j,
            report,
        })
    }
}

/// Accepted original-cell fields with distinct integrated source energy and
/// coenergy. Recoil coenergy uses the fixed-remanence Legendre reference and
/// can be signed; its stored energy remains nonnegative.
pub struct SpatialNonlinearField<'m> {
    problem: &'m SpatialNonlinearMagnetostatics,
    potential_a: Vec<f64>,
    relative_potential_a: Vec<f64>,
    field_strength_a_m: Vec<[f64; 3]>,
    flux_density_t: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    cell_coenergy_j: Vec<f64>,
    energy_j: f64,
    coenergy_j: f64,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    nodal_boundary_flux_wb: Vec<f64>,
    interface_flux_jumps_wb: Vec<(usize, f64)>,
    gauge_vertices: Vec<usize>,
    legendre_residual_j: f64,
    energy_pairing_residual_j: f64,
    report: SpatialNonlinearReport,
}
impl<'m> SpatialNonlinearField<'m> {
    pub fn problem(&self) -> &'m SpatialNonlinearMagnetostatics {
        self.problem
    }
    pub fn potential_a(&self) -> &[f64] {
        &self.potential_a
    }
    pub fn relative_potential_a(&self) -> &[f64] {
        &self.relative_potential_a
    }
    pub fn field_strength_a_m(&self) -> &[[f64; 3]] {
        &self.field_strength_a_m
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn cell_coenergy_j(&self) -> &[f64] {
        &self.cell_coenergy_j
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
    pub fn nodal_boundary_flux_wb(&self) -> &[f64] {
        &self.nodal_boundary_flux_wb
    }
    pub fn interface_flux_jumps_wb(&self) -> &[(usize, f64)] {
        &self.interface_flux_jumps_wb
    }
    pub fn gauge_vertices(&self) -> &[usize] {
        &self.gauge_vertices
    }
    pub fn legendre_residual_j(&self) -> f64 {
        self.legendre_residual_j
    }
    pub fn energy_pairing_residual_j(&self) -> f64 {
        self.energy_pairing_residual_j
    }
    pub fn report(&self) -> &SpatialNonlinearReport {
        &self.report
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
}
