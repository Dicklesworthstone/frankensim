//! History-dependent small-slope bending on the original DKT triangles.
//!
//! Each original area point owns a resolved plane-stress history through the
//! thickness. The signed DKT curvature gives `epsilon(z) = z * kappa`; its
//! twisting entry is engineering shear. Moments and the algorithmic bending
//! tangent integrate the SAME returned stresses and condensed material tangent.
//! Composite three-point Gauss–Lobatto thickness quadrature includes both faces
//! and integrates elastic `z^2` exactly. Plastic fronts require refinement.
//!
//! Midsurface extension is constrained to zero. Membrane resultants are reported
//! as reactions, not silently relaxed. This is a homogeneous, flat, thin-plate
//! bending model, not a finite-rotation shell, a membrane/forming-limit model,
//! or a prediction of necking, fracture, thickness change or plastic heating.

use crate::{PlateError, PlateMesh, dkt_curvature_quadrature, edge_coeffs, h_vectors};
use fs_material::SmallStrainLaw;
use fs_material::plastic::reduction::{
    PlaneStressReduction, PlaneStressWorkspace, ReducedState, ReductionConfig, ReductionError,
};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};
use std::sync::Arc;

/// Physical domain, thickness resolution and retained-work bounds.
#[derive(Clone, Copy, Debug)]
pub struct BendingConfig {
    /// Uniform reference thickness [m].
    pub thickness_m: f64,
    /// Equal thickness cells, with `2 * cells + 1` distinct Lobatto points.
    pub thickness_cells: usize,
    /// Native plane-stress return and full tensor-strain controls.
    pub reduction: ReductionConfig,
    /// Maximum norm of the quadratic DKT rotation field, in `(0, 0.05]`.
    pub maximum_slope: f64,
    /// Maximum full plate displacement/slope coordinates.
    pub maximum_dofs: usize,
    /// Maximum retained three-dimensional material-point histories.
    pub maximum_material_points: usize,
}

/// Bounded quasistatic equilibrium with independently dimensioned tolerances.
#[derive(Clone, Copy, Debug)]
pub struct BendingSolveControls {
    /// Maximum Newton corrections in this physical increment.
    pub maximum_iterations: usize,
    /// Maximum halvings of one Newton correction.
    pub maximum_backtracks: usize,
    /// Maximum free transverse nodal-force imbalance [N].
    pub force_tolerance_n: f64,
    /// Maximum free nodal-moment imbalance [N m].
    pub moment_tolerance_nm: f64,
    /// Shared sparse LDLT numeric factor and scratch memory limit [bytes].
    pub maximum_factor_bytes: usize,
}

/// Refusal leaves the caller's accepted displacement and histories unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum BendingError {
    /// Invalid, unsupported, nonfinite or out-of-budget input.
    Invalid(&'static str),
    /// Original geometry or shared sparse-solver refusal.
    Plate(PlateError),
    /// Original plane-stress material refusal at one integration point.
    Material {
        /// Triangle index.
        element: usize,
        /// Original DKT area point, in `0..3`.
        point: usize,
        /// Through-thickness point, from bottom to top.
        layer: usize,
        /// Native local-solve failure.
        source: ReductionError,
    },
    /// No admissible equilibrium in the requested Newton/line-search budget.
    NonConvergence {
        /// Newton corrections attempted.
        iterations: usize,
        /// Largest free force imbalance [N].
        force_residual_n: f64,
        /// Largest free moment imbalance [N m].
        moment_residual_nm: f64,
    },
}
impl core::fmt::Display for BendingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-PLASTIC-BENDING: {self:?}")
    }
}
impl std::error::Error for BendingError {}
impl From<PlateError> for BendingError {
    fn from(error: PlateError) -> Self {
        Self::Plate(error)
    }
}

/// Immutable accepted history, sealed to the original plate/material owner.
#[derive(Clone, Debug, PartialEq)]
pub struct BendingState<S> {
    owner: Arc<()>,
    displacement: Vec<f64>,
    points: Vec<ReducedState<S>>,
}
impl<S> BendingState<S> {
    /// Full `(w [m], wx, wy)` per node, with slopes dimensionless.
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    /// Original three-dimensional strain/history, element-major, then the
    /// three DKT area points, then bottom-to-top thickness points.
    pub fn material_points(&self) -> &[ReducedState<S>] {
        &self.points
    }
}

/// Candidate history and its original force/tangent response.
#[derive(Debug)]
pub struct BendingResponse<S> {
    /// Publish only on acceptance of this physical increment.
    pub candidate: BendingState<S>,
    /// Full transverse forces [N] and slope-conjugate moments [N m].
    pub internal_force: Vec<f64>,
    /// Full derivative of `internal_force` with respect to nodal coordinates.
    pub tangent: Csr,
    /// `[Mxx, Myy, Mxy]` [N m / m], one per original DKT area point.
    pub moments: Vec<[f64; 3]>,
    /// `[Nxx, Nyy, Nxy]` [N/m], reactions to zero midsurface extension.
    pub membrane_resultants: Vec<[f64; 3]>,
}

/// Accepted equilibrium for one caller-chosen forming or unloading increment.
#[derive(Debug)]
pub struct BendingEquilibrium<S> {
    /// Candidate state and actual converged material response.
    pub response: BendingResponse<S>,
    /// Newton corrections, excluding local material return iterations.
    pub iterations: usize,
    /// Largest free force imbalance [N].
    pub force_residual_n: f64,
    /// Largest free moment imbalance [N m].
    pub moment_residual_nm: f64,
}

#[derive(Debug)]
struct Element {
    dofs: [usize; 9],
    curvature: [[[f64; 9]; 3]; 3],
    rotations: [[[f64; 9]; 2]; 6],
    area_weight: f64,
}

/// Prepared DKT geometry and one native homogeneous small-strain law.
///
/// The law's axes are the mesh x/y axes. Any material-frame transformation
/// belongs to that law. Geometry and constitutive coefficients cannot change
/// behind an accepted history. Preparation caches the exact original DKT
/// operators; one small plane-stress workspace serves all local return maps.
#[derive(Debug)]
pub struct PlasticBendingPlate<L> {
    owner: Arc<()>,
    config: BendingConfig,
    reduction: PlaneStressReduction<L>,
    elements: Vec<Element>,
    thickness: Vec<(f64, f64)>,
    dofs: usize,
    point_count: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BendingError> {
    if cancelled() {
        Err(PlateError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> BendingError {
    BendingError::Invalid(what)
}
fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

impl<L: SmallStrainLaw> PlasticBendingPlate<L> {
    /// Prepare original geometry and a checked material/thickness discretization.
    /// No material history is inferred from an elastic D matrix or yield scalar.
    pub fn new(
        mesh: &PlateMesh,
        law: L,
        config: BendingConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BendingError> {
        poll(&mut cancelled)?;
        let layers = config
            .thickness_cells
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("thickness-point count overflow"))?;
        let dofs = mesh
            .nodes
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("DOF count overflow"))?;
        let point_count = mesh
            .tris
            .len()
            .checked_mul(3)
            .and_then(|n| n.checked_mul(layers))
            .ok_or_else(|| invalid("material-point count overflow"))?;
        if config.thickness_cells == 0
            || config.thickness_cells > 4096
            || !config.thickness_m.is_finite()
            || config.thickness_m <= 0.0
            || !config.maximum_slope.is_finite()
            || config.maximum_slope <= 0.0
            || config.maximum_slope > 0.05
            || dofs > config.maximum_dofs
            || point_count > config.maximum_material_points
            || !law.admissibility().tangent_symmetric
        {
            return Err(invalid(
                "invalid thickness, small-slope domain, symmetric law or work budget",
            ));
        }
        mesh.validate()?;
        let reduction = PlaneStressReduction::new(law, config.reduction)
            .map_err(|_| invalid("invalid plane-stress controls"))?;
        let dz = config.thickness_m / (layers - 1) as f64;
        let thickness: Vec<_> = (0..layers)
            .map(|i| {
                let z = (i as f64 - config.thickness_cells as f64) * dz;
                let weight = dz / 3.0
                    * if i == 0 || i + 1 == layers {
                        1.0
                    } else if i % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                (z, weight)
            })
            .collect();
        if thickness
            .iter()
            .any(|(z, w)| !z.is_finite() || !w.is_finite() || *w <= 0.0)
        {
            return Err(invalid("unrepresentable thickness quadrature"));
        }
        let mut elements = Vec::with_capacity(mesh.tris.len());
        for (e, tri) in mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let (x0, y0) = mesh.nodes[tri[0]];
            let x = tri.map(|n| mesh.nodes[n].0 - x0);
            let y = tri.map(|n| mesh.nodes[n].1 - y0);
            let (curvature, twice_area) = dkt_curvature_quadrature(&x, &y, e)?;
            let cf = edge_coeffs(&x, &y);
            let mut rotations = [[[0.0; 9]; 2]; 6];
            // Quadratic Bernstein controls of the SAME DKT rotation field.
            // Convexity bounds its norm throughout the triangle, not just at nodes.
            for i in 0..6 {
                let mut shape = [0.0; 6];
                shape[i] = 1.0;
                let (rx, ry) = h_vectors(&cf, &shape);
                rotations[i] = [rx, ry];
            }
            for (edge, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                for c in 0..2 {
                    for d in 0..9 {
                        rotations[3 + edge][c][d] = 2.0 * rotations[3 + edge][c][d]
                            - 0.5 * (rotations[a][c][d] + rotations[b][c][d]);
                    }
                }
            }
            if curvature
                .iter()
                .flatten()
                .flatten()
                .chain(rotations.iter().flatten().flatten())
                .any(|v| !v.is_finite())
                || !(twice_area / 6.0 > 0.0)
            {
                return Err(invalid("unrepresentable DKT operators"));
            }
            elements.push(Element {
                dofs: std::array::from_fn(|i| 3 * tri[i / 3] + i % 3),
                curvature,
                rotations,
                area_weight: twice_area / 6.0,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            owner: Arc::new(()),
            config,
            reduction,
            elements,
            thickness,
            dofs,
            point_count,
        })
    }

    /// Original quadrature `(signed depth [m], integration weight [m])`.
    pub fn thickness_points(&self) -> &[(f64, f64)] {
        &self.thickness
    }

    /// Virgin flat history. A law with eigenstrain may initially carry reactions.
    pub fn initial_state(&self) -> BendingState<L::State> {
        BendingState {
            owner: self.owner.clone(),
            displacement: vec![0.0; self.dofs],
            points: vec![self.reduction.initial_state(); self.point_count],
        }
    }

    /// Evaluate imposed bending from immutable accepted history. Newton and
    /// line-search trials must all pass the SAME old state until convergence.
    pub fn evaluate(
        &self,
        old: &BendingState<L::State>,
        displacement: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BendingResponse<L::State>, BendingError> {
        let mut workspace = PlaneStressWorkspace::new()
            .map_err(|_| invalid("plane-stress workspace allocation failed"))?;
        self.evaluate_inner(old, displacement, &mut workspace, &mut cancelled)
    }

    fn evaluate_inner(
        &self,
        old: &BendingState<L::State>,
        displacement: &[f64],
        workspace: &mut PlaneStressWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BendingResponse<L::State>, BendingError> {
        poll(cancelled)?;
        if !Arc::ptr_eq(&self.owner, &old.owner)
            || old.points.len() != self.point_count
            || displacement.len() != self.dofs
            || displacement.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("foreign history or invalid displacement field"));
        }
        let mut states = Vec::with_capacity(self.point_count);
        let mut force = vec![0.0; self.dofs];
        let mut tangent = Coo::new(self.dofs, self.dofs);
        let mut moments = Vec::with_capacity(3 * self.elements.len());
        let mut membrane_resultants = Vec::with_capacity(3 * self.elements.len());
        for (e, element) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let u = element.dofs.map(|d| displacement[d]);
            if element.rotations.iter().any(|r| {
                let magnitude = dot(&r[0], &u).hypot(dot(&r[1], &u));
                !magnitude.is_finite() || magnitude > self.config.maximum_slope
            }) {
                return Err(invalid("DKT rotation field exceeds the small-slope domain"));
            }
            let mut ke = [[0.0; 9]; 9];
            for (q, b) in element.curvature.iter().enumerate() {
                let curvature = b.map(|row| dot(&row, &u));
                let mut m = [0.0; 3];
                let mut membrane = [0.0; 3];
                let mut d = [[0.0; 3]; 3];
                for (layer, &(z, weight)) in self.thickness.iter().enumerate() {
                    poll(cancelled)?;
                    let strain = [z * curvature[0], z * curvature[1], 0.5 * z * curvature[2]];
                    let response = self
                        .reduction
                        .evaluate_with_workspace(strain, &old.points[states.len()], workspace)
                        .map_err(|source| BendingError::Material {
                            element: e,
                            point: q,
                            layer,
                            source,
                        })?;
                    for i in 0..3 {
                        m[i] += weight * z * response.stress_pa[i];
                        membrane[i] += weight * response.stress_pa[i];
                        for j in 0..3 {
                            // Native tangent differentiates tensor shear; DKT uses gamma_xy.
                            d[i][j] += weight
                                * z
                                * z
                                * response.tangent_pa[i][j]
                                * if j == 2 { 0.5 } else { 1.0 };
                        }
                    }
                    states.push(response.state);
                }
                let scale = d.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs()));
                if !scale.is_finite()
                    || d.iter().flatten().any(|v| !v.is_finite())
                    || m.iter().chain(&membrane).any(|v| !v.is_finite())
                {
                    return Err(invalid("unrepresentable section response"));
                }
                for i in 0..3 {
                    for j in i + 1..3 {
                        if (d[i][j] - d[j][i]).abs() > 1e-10 * scale {
                            return Err(invalid("material bending tangent violates work symmetry"));
                        }
                        let symmetric = 0.5 * d[i][j] + 0.5 * d[j][i];
                        d[i][j] = symmetric;
                        d[j][i] = symmetric;
                    }
                }
                for i in 0..9 {
                    force[element.dofs[i]] +=
                        element.area_weight * (0..3).map(|a| b[a][i] * m[a]).sum::<f64>();
                    for j in i..9 {
                        let mut value = 0.0;
                        for a in 0..3 {
                            for c in 0..3 {
                                value += b[a][i] * d[a][c] * b[c][j];
                            }
                        }
                        ke[i][j] += element.area_weight * value;
                    }
                }
                moments.push(m);
                membrane_resultants.push(membrane);
            }
            for i in 0..9 {
                for j in i..9 {
                    if !ke[i][j].is_finite() {
                        return Err(invalid("unrepresentable bending stiffness"));
                    }
                    tangent.push(element.dofs[i], element.dofs[j], ke[i][j]);
                    if i != j {
                        tangent.push(element.dofs[j], element.dofs[i], ke[i][j]);
                    }
                }
            }
        }
        if force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable bending force"));
        }
        let tangent = tangent.assemble();
        if (0..self.dofs).any(|r| tangent.row(r).1.iter().any(|v| !v.is_finite())) {
            return Err(invalid("unrepresentable assembled tangent"));
        }
        poll(cancelled)?;
        Ok(BendingResponse {
            candidate: BendingState {
                owner: self.owner.clone(),
                displacement: displacement.to_vec(),
                points: states,
            },
            internal_force: force,
            tangent,
            moments,
            membrane_resultants,
        })
    }

    /// One quasistatic increment with retained plastic history and springback.
    /// `load` is the full dead force/moment vector; `prescribed` gives absolute
    /// `(DOF, value)` pairs in the original `(w [m], wx, wy)` coordinates.
    /// The caller controls the physical load path by accepting returned states.
    /// A required Newton correction with a singular or indefinite free tangent
    /// refuses. A balanced endpoint alone does not prove uniqueness or stability.
    pub fn equilibrate(
        &self,
        old: &BendingState<L::State>,
        load: &[f64],
        prescribed: &[(usize, f64)],
        controls: BendingSolveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BendingEquilibrium<L::State>, BendingError> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.owner, &old.owner)
            || load.len() != self.dofs
            || load.iter().any(|v| !v.is_finite())
            || controls.maximum_iterations == 0
            || controls.maximum_backtracks > 64
            || controls.maximum_factor_bytes == 0
            || !controls.force_tolerance_n.is_finite()
            || controls.force_tolerance_n <= 0.0
            || !controls.moment_tolerance_nm.is_finite()
            || controls.moment_tolerance_nm <= 0.0
        {
            return Err(invalid("foreign history or invalid load/solve controls"));
        }
        let mut fixed = vec![false; self.dofs];
        let mut u = old.displacement.clone();
        for &(d, value) in prescribed {
            if d >= self.dofs || fixed[d] || !value.is_finite() {
                return Err(invalid("invalid or repeated prescribed coordinate"));
            }
            fixed[d] = true;
            u[d] = value;
        }
        let free: Vec<_> = (0..self.dofs).filter(|&d| !fixed[d]).collect();
        let mut symbolic = None;
        let mut workspace = PlaneStressWorkspace::new()
            .map_err(|_| invalid("plane-stress workspace allocation failed"))?;
        let mut response = self.evaluate_inner(old, &u, &mut workspace, &mut cancelled)?;
        for iteration in 0..=controls.maximum_iterations {
            let (force, moment, score) = residual(&response.internal_force, load, &free, controls);
            if !score.is_finite() {
                return Err(invalid("nonfinite equilibrium residual"));
            }
            if score <= 1.0 {
                poll(&mut cancelled)?;
                return Ok(BendingEquilibrium {
                    response,
                    iterations: iteration,
                    force_residual_n: force,
                    moment_residual_nm: moment,
                });
            }
            let stalled = || BendingError::NonConvergence {
                iterations: iteration,
                force_residual_n: force,
                moment_residual_nm: moment,
            };
            if iteration == controls.maximum_iterations {
                return Err(stalled());
            }
            let rhs: Vec<_> = free
                .iter()
                .map(|&d| load[d] - response.internal_force[d])
                .collect();
            let correction = solve_correction(
                &response.tangent,
                &rhs,
                &free,
                &mut symbolic,
                controls.maximum_factor_bytes,
                &mut cancelled,
            )?;
            let mut fraction = 1.0;
            let mut next = None;
            for _ in 0..=controls.maximum_backtracks {
                poll(&mut cancelled)?;
                let mut trial = u.clone();
                for (i, &d) in free.iter().enumerate() {
                    trial[d] += fraction * correction[i];
                }
                match self.evaluate_inner(old, &trial, &mut workspace, &mut cancelled) {
                    Ok(candidate) => {
                        let (_, _, trial_score) =
                            residual(&candidate.internal_force, load, &free, controls);
                        if trial_score <= 1.0 || trial_score < (1.0 - 1e-4 * fraction) * score {
                            next = Some((trial, candidate));
                            break;
                        }
                    }
                    Err(BendingError::Plate(PlateError::Cancelled)) => {
                        return Err(PlateError::Cancelled.into());
                    }
                    Err(_) => {} // A smaller admissible material/domain step may succeed.
                }
                fraction *= 0.5;
            }
            let Some((next_u, next_response)) = next else {
                return Err(stalled());
            };
            u = next_u;
            response = next_response;
        }
        unreachable!("final bounded Newton iteration returns")
    }
}

fn residual(
    force: &[f64],
    load: &[f64],
    free: &[usize],
    c: BendingSolveControls,
) -> (f64, f64, f64) {
    let (mut f, mut m) = (0.0_f64, 0.0_f64);
    for &d in free {
        let r = (force[d] - load[d]).abs();
        if !r.is_finite() {
            return (r, r, f64::INFINITY);
        }
        if d % 3 == 0 {
            f = f.max(r);
        } else {
            m = m.max(r);
        }
    }
    (
        f,
        m,
        (f / c.force_tolerance_n).max(m / c.moment_tolerance_nm),
    )
}

// Domain-specific Newton correction using the existing sparse LDLT owner.
// Diagonal congruence removes mixed displacement/slope units; the symbolic
// pattern is reused throughout an increment. No constitutive trial is committed.
fn solve_correction(
    tangent: &Csr,
    rhs: &[f64],
    free: &[usize],
    symbolic: &mut Option<SymbolicLdlt>,
    maximum_factor_bytes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, BendingError> {
    poll(cancelled)?;
    let mut map = vec![None; tangent.nrows()];
    let mut scale = Vec::with_capacity(free.len());
    for (i, &d) in free.iter().enumerate() {
        map[d] = Some(i);
        let (cols, values) = tangent.row(d);
        let diagonal = cols.binary_search(&d).map(|j| values[j]).unwrap_or(0.0);
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(invalid("nonpositive free bending tangent"));
        }
        scale.push(1.0 / diagonal.sqrt());
    }
    let mut coo = Coo::new(free.len(), free.len());
    for (i, &d) in free.iter().enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        let (cols, values) = tangent.row(d);
        for (&col, &value) in cols.iter().zip(values) {
            if let Some(j) = map[col] {
                let value = (value * scale[i]) * scale[j];
                if !value.is_finite() {
                    return Err(invalid("unrepresentable equilibrated tangent"));
                }
                coo.push(i, j, value);
            }
        }
    }
    let matrix = coo.assemble();
    poll(cancelled)?;
    if symbolic.is_none() {
        *symbolic = Some(
            SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd)
                .map_err(PlateError::StaticFactor)?,
        );
    }
    poll(cancelled)?;
    let factor = symbolic
        .as_ref()
        .expect("initialized symbolic factor")
        .factor(
            &matrix,
            &LdltOptions {
                max_factor_bytes: Some(maximum_factor_bytes),
                ..LdltOptions::default()
            },
        )
        .map_err(PlateError::StaticFactor)?;
    poll(cancelled)?;
    if factor.inertia().positive != free.len() || factor.inertia().negative != 0 {
        return Err(invalid("free bending tangent is not positive definite"));
    }
    let b: Vec<_> = rhs.iter().zip(&scale).map(|(r, s)| r * s).collect();
    if b.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable scaled force"));
    }
    let mut x = factor.solve(&b);
    let mut remaining = vec![0.0; b.len()];
    let rhs_norm = b.iter().fold(0.0_f64, |n, v| n.hypot(*v));
    if !rhs_norm.is_finite() {
        return Err(invalid("unrepresentable scaled force norm"));
    }
    for attempt in 0..3 {
        poll(cancelled)?;
        matrix.spmv(&x, &mut remaining);
        for (r, b) in remaining.iter_mut().zip(&b) {
            *r = *b - *r;
        }
        let error = remaining.iter().fold(0.0_f64, |n, v| n.hypot(*v));
        if error.is_finite() && error <= 1e-10 * rhs_norm {
            for (v, s) in x.iter_mut().zip(&scale) {
                *v *= s;
            }
            if x.iter().any(|v| !v.is_finite()) {
                return Err(invalid("unrepresentable Newton correction"));
            }
            return Ok(x);
        }
        if attempt < 2 {
            let delta = factor.solve(&remaining);
            for (v, d) in x.iter_mut().zip(delta) {
                *v += d;
            }
        }
    }
    Err(invalid(
        "shared sparse bending correction missed its true-residual gate",
    ))
}
