//! 3D exterior potential flow: constant SOURCE panels with collocation
//! Neumann conditions. The influence matrix row is
//! `A_ij = n_i · ∇_x G(x_i, y_j) · area_j` with the flat-panel
//! centroid approximation for well-separated pairs and the outside-limit
//! jump `-σ/2` on the diagonal — a documented screening-grade
//! discretization whose measured convergence on the sphere IS the
//! gate. The GMRES matvec runs three fs-fmm passes (the gradient's
//! components are each a smooth kernel) dotted with target normals;
//! the dense assembly is the oracle it must match.

use crate::BemError;
use fs_fmm::{Fmm, Kernel};
use fs_geom::Point3;
use fs_rep_mesh::shapes::icosphere;
use fs_solver::krylov::{GmresState, SolveReport};
use fs_solver::op::LinearOp;
use std::cell::RefCell;

/// Largest generated sphere panelization admitted by the v1 constructors.
pub const MAX_SURFACE_PANELS: usize = 81_920;
/// Largest dense influence matrix admitted by the v1 oracle path.
pub const MAX_DENSE_PANELS: usize = 2_048;
/// Largest subdivision count whose 20 * 4^s triangle count is admitted.
pub const MAX_ICOSPHERE_SUBDIVISIONS: u32 = 6;
const FMM_LEAF_CAPACITY: usize = 40;
const MIN_SPHERE_RADIUS: f64 = 1.0e-12;
const MAX_SPHERE_RADIUS: f64 = 1.0e50;
const MAX_ABS_SURFACE_COORDINATE: f64 = 1.0e100;
const MAX_PANEL_AREA: f64 = 1.0e200;

/// A converged exterior solve and its full Krylov evidence.
#[derive(Debug, Clone)]
pub struct ExteriorSolution {
    /// Source density per surface panel.
    pub sigma: Vec<f64>,
    /// Full convergence report, including residual history.
    pub report: SolveReport,
}

/// Exterior solve refusal. Unconverged iterates remain inspectable evidence,
/// but cannot be mistaken for an ordinary successful solution.
#[derive(Debug, Clone)]
pub enum ExteriorSolveError {
    /// Invalid input or an inadmissible BEM/FMM work request.
    Input(BemError),
    /// GMRES exhausted or broke down before satisfying the tolerance.
    NotConverged {
        /// Last iterate, exposed for diagnosis only.
        sigma: Vec<f64>,
        /// Full solver report and stall diagnosis.
        report: SolveReport,
    },
    /// The fallible FMM-backed operator refused an application during GMRES.
    OperatorRefused {
        /// First operator error.
        source: BemError,
        /// Last iterate, exposed for diagnosis only.
        sigma: Vec<f64>,
        /// Solver state at refusal.
        report: SolveReport,
    },
}

impl core::fmt::Display for ExteriorSolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input(source) => write!(f, "exterior solve input rejected: {source}"),
            Self::NotConverged { report, .. } => write!(
                f,
                "exterior GMRES did not converge after {} iterations (relative residual {}, diagnosis {:?})",
                report.iters, report.rel_residual, report.diagnosis
            ),
            Self::OperatorRefused { source, report, .. } => write!(
                f,
                "exterior FMM operator refused after {} iterations: {source}",
                report.iters
            ),
        }
    }
}

impl std::error::Error for ExteriorSolveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input(source) => Some(source),
            Self::NotConverged { .. } => None,
            Self::OperatorRefused { source, .. } => Some(source),
        }
    }
}

impl From<BemError> for ExteriorSolveError {
    fn from(source: BemError) -> Self {
        Self::Input(source)
    }
}

/// One TRUE gradient component of the Laplace kernel:
/// `∂G/∂x_c = −(x_c − y_c)/(4π|x−y|³)` — sign conventions verified by
/// the uniform-sphere Gauss identity in the battery (row action on
/// ones ≈ −1).
struct GradKernel {
    c: usize,
}

impl Kernel for GradKernel {
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> f64 {
        let d = [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 < 1e-300 {
            return 0.0;
        }
        let r = r2.sqrt();
        -d[self.c] / (4.0 * std::f64::consts::PI * r2 * r)
    }
}

/// A panelized closed surface: centroids, outward normals, areas.
#[derive(Debug, Clone)]
pub struct SpherePanels {
    /// Panel centroids.
    centroids: Vec<[f64; 3]>,
    /// Outward unit normals.
    normals: Vec<[f64; 3]>,
    /// Panel areas.
    areas: Vec<f64>,
}

impl SpherePanels {
    /// Construct a panel surface from parallel centroid/normal/area vectors.
    pub fn new(
        centroids: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        areas: Vec<f64>,
    ) -> Result<Self, BemError> {
        let panels = Self {
            centroids,
            normals,
            areas,
        };
        panels.validate()?;
        Ok(panels)
    }

    /// Panel centroids.
    #[must_use]
    pub fn centroids(&self) -> &[[f64; 3]] {
        &self.centroids
    }

    /// Outward unit normals.
    #[must_use]
    pub fn normals(&self) -> &[[f64; 3]] {
        &self.normals
    }

    /// Panel areas.
    #[must_use]
    pub fn areas(&self) -> &[f64] {
        &self.areas
    }

    /// Panelize an icosphere (fs-rep-mesh) of given radius/subdivisions.
    pub fn icosphere(radius: f64, subdivisions: u32) -> Result<SpherePanels, BemError> {
        if !radius.is_finite() || !(MIN_SPHERE_RADIUS..=MAX_SPHERE_RADIUS).contains(&radius) {
            return Err(BemError::InvalidScalar {
                name: "sphere radius",
                value: radius,
                requirement: "finite and in [1e-12, 1e50]",
            });
        }
        if subdivisions > MAX_ICOSPHERE_SUBDIVISIONS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "icosphere subdivisions",
                requested: subdivisions as usize,
                max: MAX_ICOSPHERE_SUBDIVISIONS as usize,
            });
        }
        let panel_count = 20usize.checked_mul(4usize.pow(subdivisions)).ok_or(
            BemError::WorkEnvelopeExceeded {
                operation: "icosphere panels",
                requested: usize::MAX,
                max: MAX_SURFACE_PANELS,
            },
        )?;
        if panel_count > MAX_SURFACE_PANELS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "icosphere panels",
                requested: panel_count,
                max: MAX_SURFACE_PANELS,
            });
        }
        let soup = icosphere(Point3::new(0.0, 0.0, 0.0), radius, subdivisions);
        let mut centroids = Vec::new();
        let mut normals = Vec::new();
        let mut areas = Vec::new();
        for (values, operation) in [
            (&mut centroids, "sphere centroids"),
            (&mut normals, "sphere normals"),
        ] {
            values
                .try_reserve_exact(panel_count)
                .map_err(|_| BemError::AllocationFailed { operation })?;
        }
        areas
            .try_reserve_exact(panel_count)
            .map_err(|_| BemError::AllocationFailed {
                operation: "sphere areas",
            })?;
        for t in 0..soup.triangles.len() {
            let [a, b, c] = soup.tri(t);
            let cx = [
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            ];
            let u = [b.x - a.x, b.y - a.y, b.z - a.z];
            let v = [c.x - a.x, c.y - a.y, c.z - a.z];
            let mut n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let area = 0.5 * len;
            for x in &mut n {
                *x /= len;
            }
            // Outward (sphere centered at origin): flip if pointing in.
            if n[0] * cx[0] + n[1] * cx[1] + n[2] * cx[2] < 0.0 {
                for x in &mut n {
                    *x = -*x;
                }
            }
            centroids.push(cx);
            normals.push(n);
            areas.push(area);
        }
        SpherePanels::new(centroids, normals, areas)
    }

    /// Validate panel shape, finiteness, normal length, and work bounds.
    pub fn validate(&self) -> Result<(), BemError> {
        let n = self.centroids.len();
        if n == 0 || n > MAX_SURFACE_PANELS {
            return Err(BemError::InvalidPanelCount {
                count: n,
                min: 1,
                max: MAX_SURFACE_PANELS,
                even_required: false,
            });
        }
        if self.normals.len() != n {
            return Err(BemError::PanelDataLength {
                field: "panel normals",
                expected: n,
                actual: self.normals.len(),
            });
        }
        if self.areas.len() != n {
            return Err(BemError::PanelDataLength {
                field: "panel areas",
                expected: n,
                actual: self.areas.len(),
            });
        }
        for (i, ((centroid, normal), area)) in self
            .centroids
            .iter()
            .zip(&self.normals)
            .zip(&self.areas)
            .enumerate()
        {
            if centroid
                .iter()
                .any(|x| !x.is_finite() || x.abs() > MAX_ABS_SURFACE_COORDINATE)
            {
                return Err(BemError::InvalidSurfacePanel {
                    index: i,
                    field: "centroid",
                });
            }
            if normal.iter().any(|x| !x.is_finite()) {
                return Err(BemError::InvalidSurfacePanel {
                    index: i,
                    field: "normal",
                });
            }
            let norm2 = normal[0].mul_add(
                normal[0],
                normal[1].mul_add(normal[1], normal[2] * normal[2]),
            );
            if (norm2 - 1.0).abs() > 1e-8 {
                return Err(BemError::InvalidSurfacePanel {
                    index: i,
                    field: "unit normal",
                });
            }
            if !area.is_finite() || *area <= 0.0 || *area > MAX_PANEL_AREA {
                return Err(BemError::InvalidSurfacePanel {
                    index: i,
                    field: "area",
                });
            }
        }
        Ok(())
    }

    /// Dense influence matrix (row-major n×n): normal velocity at
    /// centroid i induced by unit source density on panel j.
    pub fn dense_matrix(&self) -> Result<Vec<f64>, BemError> {
        self.validate()?;
        let n = self.centroids.len();
        if n > MAX_DENSE_PANELS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "dense 3D BEM panels",
                requested: n,
                max: MAX_DENSE_PANELS,
            });
        }
        let entries = n.checked_mul(n).ok_or(BemError::WorkEnvelopeExceeded {
            operation: "dense 3D BEM matrix",
            requested: usize::MAX,
            max: MAX_DENSE_PANELS * MAX_DENSE_PANELS,
        })?;
        let gk = [
            GradKernel { c: 0 },
            GradKernel { c: 1 },
            GradKernel { c: 2 },
        ];
        let mut a = zeroed_f64(entries, "dense 3D BEM matrix")?;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    a[i * n + j] = -0.5; // outside-limit jump −σ/2
                    continue;
                }
                let mut v = 0.0;
                for (c, k) in gk.iter().enumerate() {
                    v += self.normals[i][c] * k.eval(self.centroids[i], self.centroids[j]);
                }
                a[i * n + j] = v * self.areas[j];
            }
        }
        if a.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "dense 3D BEM matrix",
            });
        }
        Ok(a)
    }

    /// FMM-accelerated matvec of the SAME operator (three gradient
    /// passes dotted with target normals; the P2P near field inside
    /// fs-fmm keeps close pairs exact at centroid resolution).
    pub fn fmm_matvec(&self, sigma: &[f64], order: usize) -> Result<Vec<f64>, BemError> {
        self.validate()?;
        let n = self.centroids.len();
        validate_vector("source densities", sigma, n)?;
        Fmm::validate_request(&self.centroids, order, FMM_LEAF_CAPACITY)?;
        let weighted = weighted_charges(sigma, &self.areas)?;
        let mut out = zeroed_f64(n, "3D FMM panel output")?;
        for c in 0..3 {
            let k = GradKernel { c };
            let fmm = Fmm::new(
                &k,
                copied_points(&self.centroids, "3D FMM panel points")?,
                order,
                FMM_LEAF_CAPACITY,
            )?;
            let comp = fmm.potentials(&weighted)?;
            for i in 0..n {
                out[i] += self.normals[i][c] * comp[i];
            }
        }
        for i in 0..n {
            out[i] += -0.5 * sigma[i];
        }
        if out.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "3D FMM panel matvec",
            });
        }
        Ok(out)
    }

    /// FMM-accelerated transpose matvec for the same nonsymmetric
    /// collocation operator. The panel area belongs to the column
    /// index, so the transpose uses normal-weighted source charges and
    /// applies the area at the target after swapping kernel arguments.
    pub fn fmm_transpose_matvec(&self, x: &[f64], order: usize) -> Result<Vec<f64>, BemError> {
        self.validate()?;
        let n = self.centroids.len();
        validate_vector("transpose values", x, n)?;
        Fmm::validate_request(&self.centroids, order, FMM_LEAF_CAPACITY)?;
        let mut out = zeroed_f64(n, "3D FMM transpose output")?;
        for c in 0..3 {
            let mut charges = Vec::new();
            charges
                .try_reserve_exact(n)
                .map_err(|_| BemError::AllocationFailed {
                    operation: "3D transpose charges",
                })?;
            charges.extend(
                x.iter()
                    .zip(&self.normals)
                    .map(|(xi, normal)| xi * normal[c]),
            );
            let k = GradKernel { c };
            let fmm = Fmm::new(
                &k,
                copied_points(&self.centroids, "3D FMM transpose points")?,
                order,
                FMM_LEAF_CAPACITY,
            )?;
            let comp = fmm.potentials(&charges)?;
            for (oi, (ci, area)) in out.iter_mut().zip(comp.iter().zip(&self.areas)) {
                *oi -= area * ci;
            }
        }
        for (oi, xi) in out.iter_mut().zip(x) {
            *oi += -0.5 * xi;
        }
        if out.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "3D FMM transpose matvec",
            });
        }
        Ok(out)
    }
}

fn validate_vector(field: &'static str, values: &[f64], expected: usize) -> Result<(), BemError> {
    if values.len() != expected {
        return Err(BemError::VectorLength {
            field,
            expected,
            actual: values.len(),
        });
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult { operation: field });
    }
    Ok(())
}

fn weighted_charges(values: &[f64], areas: &[f64]) -> Result<Vec<f64>, BemError> {
    let mut weighted = Vec::new();
    weighted
        .try_reserve_exact(values.len())
        .map_err(|_| BemError::AllocationFailed {
            operation: "3D weighted source densities",
        })?;
    weighted.extend(values.iter().zip(areas).map(|(value, area)| value * area));
    if weighted.iter().any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "3D weighted source densities",
        });
    }
    Ok(weighted)
}

fn zeroed_f64(len: usize, operation: &'static str) -> Result<Vec<f64>, BemError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| BemError::AllocationFailed { operation })?;
    values.resize(len, 0.0);
    Ok(values)
}

fn copied_points(points: &[[f64; 3]], operation: &'static str) -> Result<Vec<[f64; 3]>, BemError> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(points.len())
        .map_err(|_| BemError::AllocationFailed { operation })?;
    copy.extend_from_slice(points);
    Ok(copy)
}

/// The GMRES operator wrapping the FMM matvec.
struct FmmOp<'a> {
    panels: &'a SpherePanels,
    order: usize,
    first_error: RefCell<Option<BemError>>,
}

impl FmmOp<'_> {
    fn apply_fallible(
        &self,
        y: &mut [f64],
        operation: impl FnOnce() -> Result<Vec<f64>, BemError>,
    ) {
        if self.first_error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        match operation() {
            Ok(values) => y.copy_from_slice(&values),
            Err(error) => {
                *self.first_error.borrow_mut() = Some(error);
                y.fill(f64::NAN);
            }
        }
    }
}

impl LinearOp for FmmOp<'_> {
    fn n(&self) -> usize {
        self.panels.centroids.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.apply_fallible(y, || self.panels.fmm_matvec(x, self.order));
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply_fallible(y, || self.panels.fmm_transpose_matvec(x, self.order));
    }
}

/// Solve the exterior Neumann problem for uniform onset flow `u_inf`:
/// source densities σ with `A·σ = −u_inf·n`.
///
/// A successful value always has `report.converged == true`. Budget
/// exhaustion and breakdown return [`ExteriorSolveError::NotConverged`] with
/// both the last iterate and complete [`SolveReport`] preserved.
pub fn solve_exterior(
    panels: &SpherePanels,
    u_inf: [f64; 3],
    order: usize,
    tol: f64,
) -> Result<ExteriorSolution, ExteriorSolveError> {
    panels.validate()?;
    if u_inf.iter().any(|value| !value.is_finite()) {
        return Err(BemError::InvalidScalar {
            name: "exterior freestream component",
            value: u_inf
                .iter()
                .copied()
                .find(|value| !value.is_finite())
                .unwrap_or(f64::NAN),
            requirement: "finite",
        }
        .into());
    }
    if !tol.is_finite() || tol <= 0.0 || tol >= 1.0 {
        return Err(BemError::InvalidScalar {
            name: "GMRES relative tolerance",
            value: tol,
            requirement: "finite and in (0, 1)",
        }
        .into());
    }
    Fmm::validate_request(&panels.centroids, order, FMM_LEAF_CAPACITY).map_err(BemError::from)?;
    let n = panels.centroids.len();
    let mut rhs = Vec::new();
    rhs.try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "3D exterior RHS",
        })?;
    rhs.extend((0..n).map(|i| {
        -(u_inf[0] * panels.normals[i][0]
            + u_inf[1] * panels.normals[i][1]
            + u_inf[2] * panels.normals[i][2])
    }));
    if rhs.iter().all(|value| *value == 0.0) {
        return Ok(ExteriorSolution {
            sigma: zeroed_f64(n, "zero-flow exterior solution")?,
            report: SolveReport {
                iters: 0,
                rel_residual: 0.0,
                converged: true,
                history: Vec::new(),
                diagnosis: None,
            },
        });
    }
    let op = FmmOp {
        panels,
        order,
        first_error: RefCell::new(None),
    };
    let mut st = GmresState::new(&rhs, 60.min(n));
    let report = st.run(&op, &rhs, tol, 8, false);
    let sigma = st.x;
    if let Some(source) = op.first_error.into_inner() {
        return Err(ExteriorSolveError::OperatorRefused {
            source,
            sigma,
            report,
        });
    }
    if report.converged && sigma.iter().all(|value| value.is_finite()) {
        Ok(ExteriorSolution { sigma, report })
    } else {
        Err(ExteriorSolveError::NotConverged { sigma, report })
    }
}

/// Surface velocity at panel centroids for a solved σ (onset +
/// induced; the tangential projection is the physical speed on the
/// body since the normal component vanishes by construction).
pub fn surface_velocity(
    panels: &SpherePanels,
    sigma: &[f64],
    u_inf: [f64; 3],
    order: usize,
) -> Result<Vec<[f64; 3]>, BemError> {
    panels.validate()?;
    if u_inf.iter().any(|value| !value.is_finite()) {
        return Err(BemError::InvalidScalar {
            name: "exterior freestream component",
            value: u_inf
                .iter()
                .copied()
                .find(|value| !value.is_finite())
                .unwrap_or(f64::NAN),
            requirement: "finite",
        });
    }
    let n = panels.centroids.len();
    validate_vector("source densities", sigma, n)?;
    Fmm::validate_request(&panels.centroids, order, FMM_LEAF_CAPACITY)?;
    let weighted = weighted_charges(sigma, &panels.areas)?;
    let mut out = Vec::new();
    out.try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "3D surface velocities",
        })?;
    out.resize(n, u_inf);
    for c in 0..3 {
        let k = GradKernel { c };
        let fmm = Fmm::new(
            &k,
            copied_points(&panels.centroids, "3D surface-velocity points")?,
            order,
            FMM_LEAF_CAPACITY,
        )?;
        let comp = fmm.potentials(&weighted)?;
        for (o, ci) in out.iter_mut().zip(&comp) {
            o[c] += ci;
        }
    }
    // Project out the (small residual) normal component.
    for (o, nrm) in out.iter_mut().zip(&panels.normals) {
        let vn = o[0] * nrm[0] + o[1] * nrm[1] + o[2] * nrm[2];
        for (oc, nc) in o.iter_mut().zip(nrm) {
            *oc -= vn * nc;
        }
    }
    if out.iter().flatten().any(|component| !component.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "3D surface velocity",
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::{FmmOp, SpherePanels};
    use fs_solver::op::LinearOp;
    use std::cell::RefCell;

    #[test]
    fn fmm_operator_transpose_path_matches_dense_oracle() {
        let panels = SpherePanels::icosphere(1.0, 2).expect("valid fixture");
        let n = panels.centroids.len();
        let op = FmmOp {
            panels: &panels,
            order: 6,
            first_error: RefCell::new(None),
        };
        #[allow(clippy::cast_precision_loss)]
        let x: Vec<f64> = (0..n)
            .map(|i| fs_math::det::cos((i as f64) * 0.19))
            .collect();
        let mut got = vec![0.0; n];
        op.apply_transpose(&x, &mut got);

        let dense = panels.dense_matrix().expect("admitted dense fixture");
        let mut want = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                want[j] += dense[i * n + j] * x[i];
            }
        }
        let scale = want.iter().map(|v| v * v).sum::<f64>().sqrt();
        let rel = want
            .iter()
            .zip(&got)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
            / scale;
        assert!(
            rel < 1e-4,
            "LinearOp::apply_transpose must match dense transpose; rel={rel:.3e}"
        );
    }

    #[test]
    fn invalid_panel_vectors_are_refused_before_fmm_math() {
        let panels = SpherePanels::icosphere(1.0, 1).expect("valid fixture");
        let mut areas = panels.areas.clone();
        areas.pop();

        assert!(SpherePanels::new(panels.centroids, panels.normals, areas).is_err());
    }

    #[test]
    fn fmm_operator_retains_the_first_fallible_error() {
        let panels = SpherePanels::icosphere(1.0, 1).expect("valid fixture");
        let op = FmmOp {
            panels: &panels,
            order: 4,
            first_error: RefCell::new(None),
        };
        let mut input = vec![0.0; panels.centroids.len()];
        input[0] = f64::NAN;
        let mut output = vec![0.0; input.len()];

        op.apply(&input, &mut output);

        assert!(output.iter().all(|value| value.is_nan()));
        assert!(op.first_error.into_inner().is_some());
    }
}
