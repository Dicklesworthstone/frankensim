//! Source-bound reversible saturation and recoil in axisymmetric magnetic fields.
//!
//! The unknown and cylindrical curl are the original `u=A_phi/r` formulation.
//! Original `IsotropicMagneticMaterial` laws supply H, differential reluctivity,
//! energy and coenergy together. Recoil directions are `[radial, axial]` in the
//! meridional mesh frame; a radial recoil direction cannot extend to the axis.
//!
//! Positive tensor Gauss quadrature on each triangle defines one consistent
//! discrete potential, residual and Hessian. Order three integrates the linear
//! limit exactly. Nonlinear quadrature is an estimate: increase its order and
//! refine the physical mesh, especially where a saturation segment changes.
//! Vertex admission bounds the entire affine B field against every original
//! source range; quadrature-point checks alone cannot replace this condition.
//!
//! The native Newton--Krylov owner provides bounded globalization and actual
//! residual checks. Sparse tangents and Jacobi diagonals are cached once per
//! Newton point. Trials leaving a material domain are rejected without clamping
//! the field or changing remanence. This is fixed-state magnetostatics, without
//! hysteresis, irreversible demagnetization, heat, induction or force authority.

use super::{
    AxisymmetricMagneticSample, Cell, MagneticError, assembly_error, current_load, finite, poll,
    prescribed_boundary, sample_components,
};
use crate::planar::{
    nonlinear::IsotropicMagneticMaterial,
    permanent::{MagneticEnergyReference, PlanarRecoilMaterial},
};
use fs_feec::{
    gauss_legendre,
    weighted::{WeightedAssemblyLimits, WeightedError, assemble_triplets},
};
use fs_material::field_storage::permanent::RecoilOperatingPoint;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
    SolverRunProgress, norm2,
};
use fs_sparse::Csr;
use std::cell::{Cell as Flag, RefCell};

/// Geometry and quadrature work admitted before material/field preparation.
#[derive(Clone, Copy, Debug)]
pub struct AxisymmetricNonlinearAssembly {
    /// Original cell, nodal-coordinate and sparse-contribution limits.
    pub limits: WeightedAssemblyLimits,
    /// Native Gauss order on each Duffy coordinate, in 3..=16. There are
    /// order squared positive points per triangle; order three is exact for
    /// the original linear curl pairing and cellwise-constant recoil loads.
    pub quadrature_order: usize,
    /// Explicit upper bound on all retained triangle quadrature points.
    pub max_quadrature_points: usize,
}

/// Native nonlinear work limits, with residuals measured in A m2 for u [T].
#[derive(Clone, Copy, Debug)]
pub struct AxisymmetricNonlinearControls {
    /// Native globalization, true residual tolerances, and Krylov work limits.
    /// Absolute residual tolerance has units A m2; relative tolerance uses the
    /// original initial free-coordinate residual, without a gauge subtraction.
    pub solver: NewtonKrylovConfig,
    /// Maximum accepted or rejected native Newton attempts. Zero checks only
    /// the initial field; no source amplitude continuation is invented.
    pub max_newton_iterations: usize,
    /// Maximum scalar slots in the restarted native V/Z/H Krylov basis.
    /// The effective restart is at most the number of free coordinates.
    pub max_krylov_basis_entries: usize,
}

/// A source refusal or a bounded solve that cannot publish a magnetic field.
#[derive(Debug)]
pub enum AxisymmetricNonlinearError {
    /// Original magnetic source, metric, arithmetic, assembly or cancellation error.
    Magnetic(MagneticError),
    /// The affine cell field exceeds the original reversible source's B range.
    FieldRange {
        /// Original triangle index.
        triangle: usize,
        /// Observed magnitude [T].
        observed_t: f64,
        /// Original source endpoint [T].
        maximum_t: f64,
    },
    /// Actual H exceeds the original source-backed recoil operating envelope.
    RecoilRange {
        /// Original triangle index.
        triangle: usize,
        /// Component along the original specimen direction [A/m].
        parallel_field_a_m: f64,
        /// Magnitude perpendicular to that direction [A/m].
        transverse_field_a_m: f64,
    },
    /// Native solver configuration or initial numerical-state refusal.
    Solver(NewtonError),
    /// Native work exhausted or globalization stalled; no field is published.
    Nonconvergence(NewtonReport),
}
impl From<MagneticError> for AxisymmetricNonlinearError {
    fn from(error: MagneticError) -> Self {
        Self::Magnetic(error)
    }
}
impl core::fmt::Display for AxisymmetricNonlinearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "axisymmetric nonlinear magnetostatics: {self:?}")
    }
}
impl std::error::Error for AxisymmetricNonlinearError {}
type Result<T> = core::result::Result<T, AxisymmetricNonlinearError>;

struct QuadraturePoint {
    shape: [f64; 3],
    curl: [[f64; 2]; 3],
    volume_m3: f64,
}

struct Response {
    h: [f64; 2],
    direction: [f64; 2],
    transverse_reluctivity: f64,
    longitudinal_reluctivity: f64,
    energy_density: f64,
    coenergy_density: f64,
    recoil: Option<RecoilOperatingPoint>,
}

/// Original source laws and geometry prepared once for repeated static solves.
/// Models are the existing linear, reversible H(B), and oriented recoil owners.
/// Recoil vectors refer to the meridional [r,z] frame, with azimuthal symmetry.
pub struct AxisymmetricNonlinearMagnetostatics {
    mesh: TriComplex2,
    materials: Vec<IsotropicMagneticMaterial>,
    cells: Vec<Cell>,
    quadrature: Vec<Vec<QuadraturePoint>>,
    quadrature_order: usize,
    angular_span: f64,
}

struct Evaluation {
    flux: Vec<[[f64; 2]; 3]>,
    recoil: Vec<Option<[RecoilOperatingPoint; 3]>>,
    energies: Vec<f64>,
    coenergies: Vec<f64>,
    tangent: Vec<[[f64; 3]; 3]>,
    residual: Vec<f64>,
    energy: f64,
    coenergy: f64,
}

impl AxisymmetricNonlinearMagnetostatics {
    /// Bind original material models to the original meridional cells. Positive
    /// quadrature, all work extents, and axis/source compatibility are admitted
    /// before a prepared owner is returned. No constitutive fit is generated.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<IsotropicMagneticMaterial>,
        assembly: AxisymmetricNonlinearAssembly,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let Metric2::Axisymmetric { angular_span } = mesh.metric() else {
            return Err(
                MagneticError::Invalid("nonlinear A_phi requires axisymmetric metric").into(),
            );
        };
        let count = mesh.faces().len();
        let order = assembly.quadrature_order;
        if count == 0
            || materials.len() != count
            || !(3..=16).contains(&order)
            || count > assembly.limits.max_cells
            || mesh.vertices().len() > assembly.limits.max_dofs
            || count
                .checked_mul(9)
                .is_none_or(|n| n > assembly.limits.max_triplets)
            || count
                .checked_mul(order)
                .and_then(|n| n.checked_mul(order))
                .is_none_or(|n| n > assembly.max_quadrature_points)
        {
            return Err(MagneticError::Invalid(
                "complete nonlinear materials and bounded positive triangle quadrature required",
            )
            .into());
        }
        let (nodes, weights) = gauss_legendre(order);
        let mut cells = Vec::with_capacity(count);
        let mut quadrature = Vec::with_capacity(count);
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            materials[cell].response(0.0)?;
            let points = face.map(|v| mesh.vertices()[v as usize]);
            if let IsotropicMagneticMaterial::Recoil(recoil) = &materials[cell]
                && recoil.remanence_t()[0] != 0.0
                && points.iter().any(|p| p[0] == 0.0)
            {
                return Err(MagneticError::Invalid(
                    "radial remanence cannot extend to the symmetry axis",
                )
                .into());
            }
            let geometry = Cell::new(points, angular_span)?;
            let volume = mesh
                .face_measure(cell)
                .ok_or(MagneticError::Invalid("missing meridional cell volume"))?;
            let mean_radius = finite(points.iter().map(|p| p[0] / 3.0).sum())?;
            if volume <= 0.0 || mean_radius <= 0.0 {
                return Err(MagneticError::Invalid("positive swept cell volume required").into());
            }
            let mut samples = Vec::with_capacity(order * order);
            for (&x, &wx) in nodes.iter().zip(&weights) {
                for (&y, &wy) in nodes.iter().zip(&weights) {
                    let xi = 0.5 * (x + 1.0);
                    let eta = 0.5 * (y + 1.0);
                    let shape = [1.0 - xi, xi * (1.0 - eta), xi * eta];
                    let radius = finite((0..3).map(|i| shape[i] * points[i][0]).sum())?;
                    // Duffy coordinates: dA = 2*area*xi*d_xi*d_eta. Native
                    // Gauss weights live on [-1,1], hence the factor 1/2.
                    let volume_m3 = finite(volume * (radius / mean_radius) * (0.5 * xi * wx * wy))?;
                    if volume_m3 <= 0.0 {
                        return Err(MagneticError::Invalid(
                            "positive representable quadrature weight required",
                        )
                        .into());
                    }
                    let curl = core::array::from_fn(|i| {
                        [
                            -radius * geometry.gradients[i][1],
                            radius.mul_add(geometry.gradients[i][0], 2.0 * shape[i]),
                        ]
                    });
                    if curl.iter().flatten().any(|v| !v.is_finite()) {
                        return Err(
                            MagneticError::Invalid("unrepresentable cylindrical basis").into()
                        );
                    }
                    samples.push(QuadraturePoint {
                        shape,
                        curl,
                        volume_m3,
                    });
                }
            }
            cells.push(geometry);
            quadrature.push(samples);
        }
        poll(cancelled)?;
        Ok(Self {
            mesh,
            materials,
            cells,
            quadrature,
            quadrature_order: order,
            angular_span,
        })
    }

    /// Original mesh and durable feature identities.
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete original linear, saturation or recoil source laws.
    pub fn materials(&self) -> &[IsotropicMagneticMaterial] {
        &self.materials
    }
    /// Retained native Gauss order on each Duffy coordinate.
    pub fn quadrature_order(&self) -> usize {
        self.quadrature_order
    }
    /// Declared rotationally invariant sweep angle [rad].
    pub fn angular_span_rad(&self) -> f64 {
        self.angular_span
    }

    fn response(&self, triangle: usize, b: [f64; 2]) -> Result<Response> {
        let material = &self.materials[triangle];
        let br = match material {
            IsotropicMagneticMaterial::Recoil(recoil) => recoil.remanence_t(),
            _ => [0.0; 2],
        };
        let induced = [finite(b[0] - br[0])?, finite(b[1] - br[1])?];
        let magnitude = finite(induced[0].hypot(induced[1]))?;
        if let IsotropicMagneticMaterial::Reversible(law) = material
            && magnitude > law.maximum_flux_density_t()
        {
            return Err(AxisymmetricNonlinearError::FieldRange {
                triangle,
                observed_t: magnitude,
                maximum_t: law.maximum_flux_density_t(),
            });
        }
        let response = material.response(magnitude)?;
        let nu = response.secant_reluctivity_m_h;
        let longitudinal = response.differential_reluctivity_m_h;
        if !nu.is_finite() || nu <= 0.0 || !longitudinal.is_finite() || longitudinal <= 0.0 {
            return Err(MagneticError::Invalid(
                "positive source magnetic differential and secant reluctivities required",
            )
            .into());
        }
        let h = [finite(nu * induced[0])?, finite(nu * induced[1])?];
        let recoil = match material {
            IsotropicMagneticMaterial::Recoil(law) => Some(admit_recoil(triangle, law, h)?),
            _ => None,
        };
        let coenergy_density =
            finite(br[0].mul_add(h[0], br[1].mul_add(h[1], response.coenergy_density_j_m3)))?;
        Ok(Response {
            h,
            direction: if magnitude == 0.0 {
                [1.0, 0.0]
            } else {
                induced.map(|v| v / magnitude)
            },
            transverse_reluctivity: nu,
            longitudinal_reluctivity: longitudinal,
            energy_density: response.energy_density_j_m3,
            coenergy_density,
            recoil,
        })
    }

    fn evaluate(
        &self,
        potential: &[f64],
        rhs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        poll(cancelled)?;
        if potential.len() != self.mesh.vertices().len() || potential.iter().any(|v| !v.is_finite())
        {
            return Err(
                MagneticError::Invalid("complete finite meridional potential required").into(),
            );
        }
        let count = self.cells.len();
        let mut evaluation = Evaluation {
            flux: Vec::with_capacity(count),
            recoil: Vec::with_capacity(count),
            energies: Vec::with_capacity(count),
            coenergies: Vec::with_capacity(count),
            tangent: Vec::with_capacity(count),
            residual: rhs.iter().map(|v| -v).collect(),
            energy: 0.0,
            coenergy: 0.0,
        };
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let points = face.map(|v| self.mesh.vertices()[v as usize]);
            let values = face.map(|v| potential[v as usize]);
            let flux = self.cells[cell].flux_vertices(values, points)?;
            // Norms of affine fields are convex. Vertex B bounds therefore
            // cover every interior point. Recoil H is also affine; parallel
            // interval and transverse-norm bounds have the same property.
            let mut operating = [None; 3];
            for vertex in 0..3 {
                operating[vertex] = self.response(cell, flux[vertex])?.recoil;
            }
            let mut tangent = [[0.0; 3]; 3];
            let mut energy = 0.0;
            let mut coenergy = 0.0;
            for q in &self.quadrature[cell] {
                poll(cancelled)?;
                let mut b = [0.0; 2];
                for (vertex, &shape) in q.shape.iter().enumerate() {
                    for (axis, component) in b.iter_mut().enumerate() {
                        *component = finite(shape.mul_add(flux[vertex][axis], *component))?;
                    }
                }
                let response = self.response(cell, b)?;
                energy = finite(q.volume_m3.mul_add(response.energy_density, energy))?;
                coenergy = finite(q.volume_m3.mul_add(response.coenergy_density, coenergy))?;
                let direction = response.direction;
                let along: [f64; 3] = core::array::from_fn(|i| {
                    direction[0].mul_add(q.curl[i][0], direction[1] * q.curl[i][1])
                });
                let across: [f64; 3] = core::array::from_fn(|i| {
                    direction[0].mul_add(q.curl[i][1], -direction[1] * q.curl[i][0])
                });
                for i in 0..3 {
                    let pairing =
                        finite(response.h[0].mul_add(q.curl[i][0], response.h[1] * q.curl[i][1]))?;
                    let v = face[i] as usize;
                    evaluation.residual[v] =
                        finite(q.volume_m3.mul_add(pairing, evaluation.residual[v]))?;
                    for j in i..3 {
                        // Positive directional Gram form avoids cancellation
                        // between large secant and small differential slopes.
                        let gram = finite(
                            if response.longitudinal_reluctivity == response.transverse_reluctivity
                            {
                                response.transverse_reluctivity
                                    * q.curl[i][0]
                                        .mul_add(q.curl[j][0], q.curl[i][1] * q.curl[j][1])
                            } else {
                                (response.longitudinal_reluctivity * along[i]).mul_add(
                                    along[j],
                                    response.transverse_reluctivity * across[i] * across[j],
                                )
                            },
                        )?;
                        tangent[i][j] = finite(q.volume_m3.mul_add(gram, tangent[i][j]))?;
                        tangent[j][i] = tangent[i][j];
                    }
                }
            }
            evaluation.energy = finite(evaluation.energy + energy)?;
            evaluation.coenergy = finite(evaluation.coenergy + coenergy)?;
            evaluation.energies.push(energy);
            evaluation.coenergies.push(coenergy);
            evaluation.flux.push(flux);
            evaluation.recoil.push(match operating {
                [Some(a), Some(b), Some(c)] => Some([a, b, c]),
                _ => None,
            });
            evaluation.tangent.push(tangent);
        }
        Ok(evaluation)
    }

    /// Solve one static excitation, starting from a complete admissible u [T].
    /// Currents are original cellwise J_phi [A/m2]. Boundary and axis semantics
    /// are those of the original linear axisymmetric owner. Initial values must
    /// match prescribed traces exactly and satisfy every complete source range.
    /// No accepted state is mutated by failed trials, cancellation or refusal.
    pub fn solve(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_over_radius_t: &[(usize, f64)],
        initial_potential_over_radius_t: &[f64],
        controls: AxisymmetricNonlinearControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricNonlinearField<'_>> {
        poll(cancelled)?;
        controls
            .solver
            .validate()
            .map_err(AxisymmetricNonlinearError::Solver)?;
        let prescribed =
            prescribed_boundary(&self.mesh, boundary_potential_over_radius_t, cancelled)?;
        let n = prescribed.len();
        if initial_potential_over_radius_t.len() != n
            || initial_potential_over_radius_t
                .iter()
                .any(|v| !v.is_finite())
            || prescribed
                .iter()
                .zip(initial_potential_over_radius_t)
                .any(|(p, &value)| p.is_some_and(|p| p != value))
        {
            return Err(MagneticError::Invalid(
                "initial u must be finite, complete, and match its traces",
            )
            .into());
        }
        let rhs = current_load(&self.mesh, &self.cells, current_density_a_m2, cancelled)?;
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; n];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut config = controls.solver;
        config.linear_restart = config.linear_restart.min(free.len().max(1));
        let width = config.linear_restart;
        let basis_entries = if free.is_empty() || controls.max_newton_iterations == 0 {
            Some(0)
        } else {
            width
                .checked_mul(2)
                .and_then(|v| v.checked_add(1))
                .and_then(|v| v.checked_mul(free.len()))
                .and_then(|v| {
                    width
                        .checked_add(1)
                        .and_then(|w| width.checked_mul(w))
                        .and_then(|w| v.checked_add(w))
                })
        };
        if basis_entries.is_none_or(|v| v > controls.max_krylov_basis_entries) {
            return Err(MagneticError::Invalid(
                "native Krylov basis exceeds declared work storage",
            )
            .into());
        }
        let initial = self.evaluate(initial_potential_over_radius_t, &rhs, cancelled)?;
        let x: Vec<_> = free
            .iter()
            .map(|&v| initial_potential_over_radius_t[v])
            .collect();
        let initial_norm = norm2(
            &free
                .iter()
                .map(|&v| initial.residual[v])
                .collect::<Vec<_>>(),
        );
        let target = config
            .absolute_tolerance
            .max(config.relative_tolerance * initial_norm);
        let cancellation = RefCell::new(&mut *cancelled);
        let interrupted = Flag::new(false);
        let check = || {
            if interrupted.get() || (*cancellation.borrow_mut())() {
                interrupted.set(true);
                true
            } else {
                false
            }
        };
        let adapter = NativeProblem {
            owner: self,
            rhs: &rhs,
            prescribed: &prescribed,
            free: &free,
            index: &index,
            cancelled: &check,
            failure: RefCell::new(None),
            cache: RefCell::new(Some(CachedEvaluation {
                x: x.clone(),
                value: initial,
                matrix: None,
            })),
        };
        let initial_solver = NewtonKrylovState::new(&adapter, x, config);
        if let Some(error) = adapter.failure.borrow_mut().take() {
            return Err(error);
        }
        if interrupted.get() {
            return Err(MagneticError::Cancelled.into());
        }
        let mut solver = initial_solver.map_err(AxisymmetricNonlinearError::Solver)?;
        let completed =
            solver.run_controlled(&adapter, controls.max_newton_iterations, &mut || check());
        if let Some(error) = adapter.failure.borrow_mut().take() {
            return Err(error);
        }
        if interrupted.get() || matches!(completed.progress, SolverRunProgress::Paused) {
            return Err(MagneticError::Cancelled.into());
        }
        if completed.progress != SolverRunProgress::Complete {
            return Err(
                MagneticError::Invalid("native nonlinear callback did not complete").into(),
            );
        }
        if !completed.report.converged {
            return Err(AxisymmetricNonlinearError::Nonconvergence(completed.report));
        }
        let mut potential = initial_potential_over_radius_t.to_vec();
        for (&v, &value) in free.iter().zip(&solver.x) {
            potential[v] = value;
        }
        // Recover from the actual accepted physical field, independently of
        // native cached residuals and preconditioned Krylov histories.
        let evaluated = self.evaluate(&potential, &rhs, &mut || check())?;
        let final_norm = norm2(
            &free
                .iter()
                .map(|&v| evaluated.residual[v])
                .collect::<Vec<_>>(),
        );
        let mut report = completed.report;
        report.residual_norm = final_norm;
        report.converged = final_norm.is_finite() && final_norm <= target;
        if !report.converged {
            return Err(AxisymmetricNonlinearError::Nonconvergence(report));
        }
        let mut current_pairing = 0.0;
        let mut boundary_pairing = 0.0;
        for (v, &u) in potential.iter().enumerate() {
            poll(&mut || check())?;
            current_pairing = finite(u.mul_add(rhs[v], current_pairing))?;
            if prescribed[v].is_some() {
                boundary_pairing = finite(u.mul_add(evaluated.residual[v], boundary_pairing))?;
            }
        }
        let work_balance_error =
            finite(evaluated.energy + evaluated.coenergy - current_pairing - boundary_pairing)?;
        poll(&mut || check())?;
        Ok(AxisymmetricNonlinearField {
            problem: self,
            potential_over_radius_t: potential,
            impressed_current_density_a_m2: current_density_a_m2.to_vec(),
            cell_vertex_flux_density_t: evaluated.flux,
            recoil_operating_points: evaluated.recoil,
            cell_energy_j: evaluated.energies,
            cell_coenergy_j: evaluated.coenergies,
            energy_j: evaluated.energy,
            coenergy_j: evaluated.coenergy,
            impressed_current_pairing_j: current_pairing,
            boundary_pairing_j: boundary_pairing,
            work_balance_error_j: work_balance_error,
            report,
        })
    }
}

fn admit_recoil(
    triangle: usize,
    recoil: &PlanarRecoilMaterial,
    h: [f64; 2],
) -> Result<RecoilOperatingPoint> {
    let direction = recoil.direction();
    let parallel = finite(direction[0].mul_add(h[0], direction[1] * h[1]))?;
    let transverse = finite(direction[0].mul_add(h[1], -direction[1] * h[0]))?.abs();
    let law = recoil.material();
    let limits = law.parallel_range_a_m();
    if parallel < limits[0] || parallel > limits[1] || transverse > law.maximum_transverse_a_m() {
        return Err(AxisymmetricNonlinearError::RecoilRange {
            triangle,
            parallel_field_a_m: parallel,
            transverse_field_a_m: transverse,
        });
    }
    law.admit_operating_point(parallel, transverse)
        .map_err(|e| MagneticError::Material(e).into())
}

struct CachedEvaluation {
    x: Vec<f64>,
    value: Evaluation,
    matrix: Option<(Csr, Vec<f64>)>,
}

struct NativeProblem<'a> {
    owner: &'a AxisymmetricNonlinearMagnetostatics,
    rhs: &'a [f64],
    prescribed: &'a [Option<f64>],
    free: &'a [usize],
    index: &'a [usize],
    cancelled: &'a dyn Fn() -> bool,
    failure: RefCell<Option<AxisymmetricNonlinearError>>,
    cache: RefCell<Option<CachedEvaluation>>,
}
impl NativeProblem<'_> {
    fn prepare(&self, x: &[f64], tangent: bool) -> Result<()> {
        poll(&mut || (self.cancelled)())?;
        let matches = self
            .cache
            .borrow()
            .as_ref()
            .is_some_and(|cache| cache.x == x);
        if !matches {
            let mut potential: Vec<_> = self.prescribed.iter().map(|v| v.unwrap_or(0.0)).collect();
            for (&v, &value) in self.free.iter().zip(x) {
                potential[v] = value;
            }
            let value = self
                .owner
                .evaluate(&potential, self.rhs, &mut || (self.cancelled)())?;
            *self.cache.borrow_mut() = Some(CachedEvaluation {
                x: x.to_vec(),
                value,
                matrix: None,
            });
        }
        if tangent
            && self
                .cache
                .borrow()
                .as_ref()
                .is_some_and(|cache| cache.matrix.is_none())
        {
            let mut cache = self.cache.borrow_mut();
            let cache = cache.as_mut().expect("prepared native field");
            let mut entries = Vec::new();
            entries
                .try_reserve_exact(self.owner.mesh.faces().len() * 9)
                .map_err(|_| MagneticError::Assembly(WeightedError::Allocation))?;
            for (cell, face) in self.owner.mesh.faces().iter().enumerate() {
                poll(&mut || (self.cancelled)())?;
                for i in 0..3 {
                    for j in 0..3 {
                        let a = self.index[face[i] as usize];
                        let b = self.index[face[j] as usize];
                        if a != usize::MAX && b != usize::MAX {
                            entries.push((a, b, cache.value.tangent[cell][i][j]));
                        }
                    }
                }
            }
            let matrix = assemble_triplets(self.free.len(), entries, &mut || (self.cancelled)())
                .map_err(assembly_error)?;
            let mut diagonal = Vec::with_capacity(self.free.len());
            for v in 0..self.free.len() {
                poll(&mut || (self.cancelled)())?;
                let value = matrix.get(v, v);
                if value <= 0.0 || !value.recip().is_finite() {
                    return Err(MagneticError::Invalid(
                        "positive representable nonlinear magnetic diagonal required",
                    )
                    .into());
                }
                diagonal.push(value.recip());
            }
            cache.matrix = Some((matrix, diagonal));
        }
        Ok(())
    }
    fn admit_callback(&self, x: &[f64], tangent: bool, output: &mut [f64]) -> bool {
        if self.failure.borrow().is_some() {
            output.fill(f64::NAN);
            return false;
        }
        match self.prepare(x, tangent) {
            Ok(()) => true,
            Err(error) => {
                output.fill(f64::NAN);
                // The native line search contracts an inadmissible trial.
                // A cancellation, arithmetic or assembly failure is terminal.
                if !matches!(
                    error,
                    AxisymmetricNonlinearError::FieldRange { .. }
                        | AxisymmetricNonlinearError::RecoilRange { .. }
                ) {
                    *self.failure.borrow_mut() = Some(error);
                }
                false
            }
        }
    }
}
impl NonlinearProblem for NativeProblem<'_> {
    fn dimension(&self) -> usize {
        self.free.len()
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        if !self.admit_callback(x, false, output) {
            return;
        }
        let cache = self.cache.borrow();
        let value = &cache.as_ref().expect("prepared native field").value;
        for (out, &v) in output.iter_mut().zip(self.free) {
            *out = value.residual[v];
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        if !self.admit_callback(x, true, output) {
            return;
        }
        let cache = self.cache.borrow();
        cache
            .as_ref()
            .expect("prepared native field")
            .matrix
            .as_ref()
            .expect("prepared tangent")
            .0
            .spmv(direction, output);
    }
    fn preconditioner_apply(
        &self,
        x: &[f64],
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        if !self.admit_callback(x, true, output) {
            return;
        }
        let cache = self.cache.borrow();
        let diagonal = &cache
            .as_ref()
            .expect("prepared native field")
            .matrix
            .as_ref()
            .expect("prepared tangent")
            .1;
        for ((out, &value), &inverse) in output.iter_mut().zip(residual).zip(diagonal) {
            *out = value * inverse;
        }
    }
}

/// Accepted static estimate tied to its original material and geometry owners.
/// A_phi is quadratic and B affine per cell. H follows the original nonlinear
/// law at each query point and is never interpolated across a material interface.
pub struct AxisymmetricNonlinearField<'a> {
    problem: &'a AxisymmetricNonlinearMagnetostatics,
    potential_over_radius_t: Vec<f64>,
    impressed_current_density_a_m2: Vec<f64>,
    cell_vertex_flux_density_t: Vec<[[f64; 2]; 3]>,
    recoil_operating_points: Vec<Option<[RecoilOperatingPoint; 3]>>,
    cell_energy_j: Vec<f64>,
    cell_coenergy_j: Vec<f64>,
    energy_j: f64,
    coenergy_j: f64,
    impressed_current_pairing_j: f64,
    boundary_pairing_j: f64,
    work_balance_error_j: f64,
    report: NewtonReport,
}
impl<'a> AxisymmetricNonlinearField<'a> {
    /// Original immutable mesh, material sources and quadrature preparation.
    pub fn problem(&self) -> &'a AxisymmetricNonlinearMagnetostatics {
        self.problem
    }
    /// Original nodal u=A_phi/r [T], including finite axis values.
    pub fn potential_over_radius_t(&self) -> &[f64] {
        &self.potential_over_radius_t
    }
    /// Original cellwise azimuthal impressed current [A/m2].
    pub fn impressed_current_density_a_m2(&self) -> &[f64] {
        &self.impressed_current_density_a_m2
    }
    /// One-sided affine B values at each original triangle's three vertices [T].
    pub fn cell_vertex_flux_density_t(&self) -> &[[[f64; 2]; 3]] {
        &self.cell_vertex_flux_density_t
    }
    /// Whole-cell recoil admission is established from these three affine-H
    /// vertex states. Margins concern this numerical field, not an error bound.
    pub fn recoil_operating_points(&self) -> &[Option<[RecoilOperatingPoint; 3]>] {
        &self.recoil_operating_points
    }
    /// Per-cell source-integrated energy, using retained positive quadrature [J].
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    /// Per-cell Legendre-dual coenergy, including Br.H for recoil cells [J].
    pub fn cell_coenergy_j(&self) -> &[f64] {
        &self.cell_coenergy_j
    }
    /// Total source-integrated energy [J]; B.H/2 is generally not this quantity.
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Total Legendre-dual coenergy [J], which may be negative with recoil.
    pub fn coenergy_j(&self) -> f64 {
        self.coenergy_j
    }
    /// Storage is zero at B=0 in induced cells, and at B=Br in recoil cells.
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if self.recoil_operating_points.iter().any(Option::is_some) {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
    /// Static original-current pairing integral J_phi*A_phi dV [J].
    pub fn impressed_current_pairing_j(&self) -> f64 {
        self.impressed_current_pairing_j
    }
    /// Static prescribed-u reaction pairing [J].
    pub fn boundary_pairing_j(&self) -> f64 {
        self.boundary_pairing_j
    }
    /// W+W* minus current and boundary pairings [J]. This is an actual static
    /// residual diagnostic, without a transient-work or continuum-bound claim.
    pub fn work_balance_error_j(&self) -> f64 {
        self.work_balance_error_j
    }
    /// Original native Newton/Krylov report; residual norms have units A m2.
    pub fn report(&self) -> &NewtonReport {
        &self.report
    }
    /// Sample A_phi, B and actual constitutive H inside one original cell.
    pub fn sample(&self, cell: usize, barycentric: [f64; 3]) -> Result<AxisymmetricMagneticSample> {
        let (a, b) = sample_components(
            &self.problem.mesh,
            &self.potential_over_radius_t,
            &self.cell_vertex_flux_density_t,
            cell,
            barycentric,
        )?;
        let response = self.problem.response(cell, b)?;
        Ok(AxisymmetricMagneticSample {
            azimuthal_potential_wb_m: a,
            flux_density_t: b,
            magnetic_field_a_m: response.h,
        })
    }
}
