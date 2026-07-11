//! 2D Hess–Smith airfoil panel method: constant source density per
//! panel plus ONE bound vortex density shared by all panels, closed by
//! the KUTTA condition (equal tangential speeds leaving the two
//! trailing-edge panels — smooth flow off the edge; circulation is
//! DETERMINED, not assumed, which is the §8.3 tie to the harmonic-
//! circulation story). Dense solve at section scale; the ADJOINT of
//! the lift coefficient is one transposed solve, FD-gated in the
//! battery. Inviscid screening honesty label applies.

use crate::BemError;
use fs_la::factor::lu;
use fs_math::det;

/// Largest admitted dense Hess-Smith section.
pub const MAX_AIRFOIL_PANELS: usize = 2_048;
const MIN_AIRFOIL_PANELS: usize = 3;
const MAX_ABS_NODE_COORDINATE: f64 = 1.0e12;
const MIN_SECTION_SCALE: f64 = 1.0e-9;
const MAX_SECTION_SCALE: f64 = 1.0e9;

fn zeroed_vec(len: usize, operation: &'static str) -> Result<Vec<f64>, BemError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| BemError::AllocationFailed { operation })?;
    values.resize(len, 0.0);
    Ok(values)
}

fn copied_vec(values: &[f64], operation: &'static str) -> Result<Vec<f64>, BemError> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| BemError::AllocationFailed { operation })?;
    copy.extend_from_slice(values);
    Ok(copy)
}

/// Cross-ISA determinism (bead 6ure): every transcendental in this
/// kernel routes through fs-math's strict `det::` implementations —
/// platform libm (macOS vs glibc) differs by ulps on sin/cos/atan2/ln/
/// sqrt, which surfaced as a bit-divergent ornith ROA between the M4
/// Pro and the 5995WX. `hypot` is composed as `det::sqrt(x·x + y·y)`
/// (panel coordinates are O(1); no overflow regime).
/// PANEL-KERNEL BIT-SEMANTICS VERSION (bead 6ure, per the y4pt golden
/// discipline): bump on ANY change that can move this kernel's output
/// bits — transcendental routing, quadrature, assembly order.
/// Downstream goldens pin this in golden-couplings.json.
pub const PANEL_BIT_SEMANTICS_VERSION: u32 = 1;

fn det_hypot(x: f64, y: f64) -> f64 {
    det::sqrt(x.mul_add(x, y * y))
}

/// A closed airfoil section: nodes ordered clockwise from the trailing
/// edge along the lower surface and back over the upper surface.
#[derive(Debug, Clone)]
pub struct Airfoil2d {
    /// Panel endpoints (closed: node N = node 0 implicitly).
    nodes: Vec<[f64; 2]>,
}

impl Airfoil2d {
    /// Construct and validate a clockwise closed section.
    pub fn new(nodes: Vec<[f64; 2]>) -> Result<Self, BemError> {
        let foil = Self { nodes };
        foil.validate()?;
        Ok(foil)
    }

    /// Read-only section nodes. Mutation must pass through [`Self::new`] so
    /// solve, adjoint, and wake paths share one geometry contract.
    #[must_use]
    pub fn nodes(&self) -> &[[f64; 2]] {
        &self.nodes
    }

    /// Validate bounds, finiteness, scale, orientation, and simple-polygon
    /// topology before any dense allocation or panel integral.
    pub fn validate(&self) -> Result<(), BemError> {
        let n = self.nodes.len();
        if !(MIN_AIRFOIL_PANELS..=MAX_AIRFOIL_PANELS).contains(&n) {
            return Err(BemError::InvalidPanelCount {
                count: n,
                min: MIN_AIRFOIL_PANELS,
                max: MAX_AIRFOIL_PANELS,
                even_required: false,
            });
        }
        for (index, node) in self.nodes.iter().enumerate() {
            for (axis, &value) in node.iter().enumerate() {
                if !value.is_finite() || value.abs() > MAX_ABS_NODE_COORDINATE {
                    return Err(BemError::InvalidNode { index, axis, value });
                }
            }
        }

        let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
        for &[x, y] in &self.nodes {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
        let scale = (max_x - min_x).max(max_y - min_y);
        if !scale.is_finite() || !(MIN_SECTION_SCALE..=MAX_SECTION_SCALE).contains(&scale) {
            return Err(BemError::InvalidScalar {
                name: "airfoil section scale",
                value: scale,
                requirement: "finite and in [1e-9, 1e9]",
            });
        }
        let min_length = scale * 1.0e-12;
        for i in 0..n {
            let a = self.nodes[i];
            let b = self.nodes[(i + 1) % n];
            if det_hypot(b[0] - a[0], b[1] - a[1]) <= min_length {
                return Err(BemError::DegeneratePanel { index: i });
            }
        }

        let mut twice_area = 0.0;
        for i in 0..n {
            let a = self.nodes[i];
            let b = self.nodes[(i + 1) % n];
            twice_area += a[0].mul_add(b[1], -(b[0] * a[1]));
        }
        let linear_tol = scale * 1.0e-14;
        let area_tol = scale * scale * 1.0e-14;
        if !twice_area.is_finite() || twice_area.abs() <= area_tol {
            return Err(BemError::DegenerateAirfoil);
        }
        if twice_area > 0.0 {
            return Err(BemError::WrongAirfoilOrientation);
        }

        for first in 0..n {
            let a = self.nodes[first];
            let b = self.nodes[(first + 1) % n];
            for second in first + 1..n {
                if second == first + 1 || (first == 0 && second == n - 1) {
                    continue;
                }
                let c = self.nodes[second];
                let d = self.nodes[(second + 1) % n];
                if segments_intersect(a, b, c, d, linear_tol, area_tol) {
                    return Err(BemError::SelfIntersectingAirfoil { first, second });
                }
            }
        }
        Ok(())
    }
}

/// The solved panel state.
#[derive(Debug, Clone)]
pub struct PanelSolution2d {
    /// Source densities per panel.
    pub sources: Vec<f64>,
    /// The shared vortex density.
    pub gamma: f64,
    /// Lift coefficient from pressure integration of the solved surface field.
    pub cl: f64,
    /// Tangential speeds at panel midpoints.
    pub vt: Vec<f64>,
}

/// Symmetric NACA 4-digit thickness form (e.g. t = 0.12 for 0012),
/// `n` panels, closed trailing edge, unit chord.
#[allow(clippy::cast_precision_loss)]
pub fn naca4_symmetric(t: f64, n: usize) -> Result<Airfoil2d, BemError> {
    if !t.is_finite() || !(0.0..=0.4).contains(&t) || t == 0.0 {
        return Err(BemError::InvalidScalar {
            name: "thickness ratio",
            value: t,
            requirement: "finite and in (0, 0.4]",
        });
    }
    if !(8..=MAX_AIRFOIL_PANELS).contains(&n) || !n.is_multiple_of(2) {
        return Err(BemError::InvalidPanelCount {
            count: n,
            min: 8,
            max: MAX_AIRFOIL_PANELS,
            even_required: true,
        });
    }
    let half = n / 2;
    let thick = |x: f64| {
        5.0 * t
            * (0.2969 * det::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
                - 0.1036 * x * x * x * x)
    };
    let mut nodes = Vec::new();
    nodes
        .try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "NACA section nodes",
        })?;
    // Cosine-clustered x from TE (1) along the LOWER surface to LE (0).
    for k in 0..half {
        let x = f64::midpoint(1.0, det::cos(std::f64::consts::PI * k as f64 / half as f64));
        nodes.push([x, -thick(x)]);
    }
    // LE to TE along the UPPER surface.
    for k in 0..half {
        let x = 0.5 * (1.0 - det::cos(std::f64::consts::PI * k as f64 / half as f64));
        nodes.push([x, thick(x)]);
    }
    Airfoil2d::new(nodes)
}

fn orientation(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]).mul_add(c[1] - a[1], -((b[1] - a[1]) * (c[0] - a[0])))
}

fn on_segment(a: [f64; 2], b: [f64; 2], p: [f64; 2], tol: f64) -> bool {
    p[0] >= a[0].min(b[0]) - tol
        && p[0] <= a[0].max(b[0]) + tol
        && p[1] >= a[1].min(b[1]) - tol
        && p[1] <= a[1].max(b[1]) + tol
}

fn segments_intersect(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
    linear_tol: f64,
    area_tol: f64,
) -> bool {
    let (o1, o2, o3, o4) = (
        orientation(a, b, c),
        orientation(a, b, d),
        orientation(c, d, a),
        orientation(c, d, b),
    );
    if ((o1 > area_tol && o2 < -area_tol) || (o1 < -area_tol && o2 > area_tol))
        && ((o3 > area_tol && o4 < -area_tol) || (o3 < -area_tol && o4 > area_tol))
    {
        return true;
    }
    (o1.abs() <= area_tol && on_segment(a, b, c, linear_tol))
        || (o2.abs() <= area_tol && on_segment(a, b, d, linear_tol))
        || (o3.abs() <= area_tol && on_segment(c, d, a, linear_tol))
        || (o4.abs() <= area_tol && on_segment(c, d, b, linear_tol))
}

struct Geometry {
    mid: Vec<[f64; 2]>,
    tangent: Vec<[f64; 2]>,
    normal: Vec<[f64; 2]>,
    len: Vec<f64>,
}

fn geometry(foil: &Airfoil2d) -> Result<Geometry, BemError> {
    foil.validate()?;
    let n = foil.nodes.len();
    let mut g = Geometry {
        mid: Vec::new(),
        tangent: Vec::new(),
        normal: Vec::new(),
        len: Vec::new(),
    };
    g.mid
        .try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "2D panel midpoints",
        })?;
    g.tangent
        .try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "2D panel tangents",
        })?;
    g.normal
        .try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "2D panel normals",
        })?;
    g.len
        .try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "2D panel lengths",
        })?;
    for i in 0..n {
        let a = foil.nodes[i];
        let b = foil.nodes[(i + 1) % n];
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let l = det_hypot(dx, dy);
        g.mid
            .push([f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])]);
        g.tangent.push([dx / l, dy / l]);
        g.normal.push([-dy / l, dx / l]);
        g.len.push(l);
    }
    Ok(g)
}

/// Velocity at `p` induced by unit source density on panel `j`
/// (analytic constant-source panel integral in panel coordinates).
fn source_velocity(foil: &Airfoil2d, g: &Geometry, j: usize, p: [f64; 2]) -> [f64; 2] {
    let a = foil.nodes[j];
    let t = g.tangent[j];
    let nrm = g.normal[j];
    // Panel-local coordinates.
    let dx = p[0] - a[0];
    let dy = p[1] - a[1];
    let xl = dx * t[0] + dy * t[1];
    let yl = dx * nrm[0] + dy * nrm[1];
    let l = g.len[j];
    let r1 = det_hypot(xl, yl).max(1e-12);
    let r2 = det_hypot(xl - l, yl).max(1e-12);
    let theta1 = det::atan2(yl, xl);
    let theta2 = det::atan2(yl, xl - l);
    let two_pi = std::f64::consts::TAU;
    // v_n = (θ₂ − θ₁)/2π — the panel-probe battery pinned the
    // reversed order (a source panel must push AWAY on both sides).
    let ul = det::ln(r1 / r2) / two_pi;
    let vl = (theta2 - theta1) / two_pi;
    [ul * t[0] + vl * nrm[0], ul * t[1] + vl * nrm[1]]
}

/// Velocity at `p` induced by unit vortex density on panel `j`
/// (the 90°-rotated source field).
fn vortex_velocity(foil: &Airfoil2d, g: &Geometry, j: usize, p: [f64; 2]) -> [f64; 2] {
    let s = source_velocity(foil, g, j, p);
    // (u, v) → (v, −u) rotates the source sheet into a vortex sheet.
    [s[1], -s[0]]
}

/// Solve at angle of attack `alpha` (radians), unit freestream.
///
/// # Errors
/// Returns [`BemError`] for invalid geometry/angle, an inadmissible dense
/// request, a singular panel system, or a non-finite result.
#[allow(clippy::too_many_lines)] // assembly + Kutta + postprocess, one narrative
pub fn solve(foil: &Airfoil2d, alpha: f64) -> Result<PanelSolution2d, BemError> {
    validate_alpha(alpha)?;
    let g = geometry(foil)?;
    let n = foil.nodes.len();
    let (dim, _) = dense_dimension(n)?;
    let u_inf = [det::cos(alpha), det::sin(alpha)];
    let (a, rhs) = assemble(foil, &g, alpha)?;
    if a.iter().chain(&rhs).any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "2D panel assembly",
        });
    }
    let f = lu(&a, dim).map_err(|source| BemError::LinearSolve {
        stage: "primal panel",
        source,
    })?;
    let mut x = rhs;
    f.solve(&mut x);
    if x.iter().any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "2D panel solution",
        });
    }
    let gamma = x[n];
    // Tangential speeds for Cp.
    let mut vt = Vec::new();
    vt.try_reserve_exact(n)
        .map_err(|_| BemError::AllocationFailed {
            operation: "2D tangential velocities",
        })?;
    for i in 0..n {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        let mut v = [u_inf[0], u_inf[1]];
        for (j, &xj) in x.iter().take(n).enumerate() {
            if i == j {
                v[0] += 0.5 * gamma * tan[0];
                v[1] += 0.5 * gamma * tan[1];
                continue;
            }
            let sv = source_velocity(foil, &g, j, mid);
            let vv = vortex_velocity(foil, &g, j, mid);
            v[0] += sv[0] * xj + vv[0] * gamma;
            v[1] += sv[1] * xj + vv[1] * gamma;
        }
        vt.push(v[0] * tan[0] + v[1] * tan[1]);
    }
    // Lift by PRESSURE INTEGRATION (v_n = 0 on the body, so
    // |v|² = v_t²): Cl = −∮ Cp·n_y ds, then rotate the force normal to
    // the freestream. Trustworthy bookkeeping — the surface field is
    // what the collocation system actually enforces.
    let mut fx = 0.0;
    let mut fy = 0.0;
    for (i, &vti) in vt.iter().enumerate() {
        let cp = 1.0 - vti * vti;
        fx += -cp * g.normal[i][0] * g.len[i];
        fy += -cp * g.normal[i][1] * g.len[i];
    }
    let cl = fy * det::cos(alpha) - fx * det::sin(alpha);
    if !cl.is_finite() || vt.iter().any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "2D pressure integration",
        });
    }
    let sources = copied_vec(&x[..n], "2D source densities")?;
    Ok(PanelSolution2d {
        sources,
        gamma,
        cl,
        vt,
    })
}

/// dCl/dα by the ADJOINT of the panel system: one transposed solve
/// for the solution sensitivity (the committed gradient path); the
/// output function's own partials (`∂f/∂x`, `∂f/∂α` at FIXED
/// solution) are cheap solve-free differences. The gradient is therefore an
/// adjoint-assisted finite-difference derivative, FD-gated end to end; it is
/// not an exact symbolic derivative.
pub fn dcl_dalpha_adjoint(foil: &Airfoil2d, alpha: f64) -> Result<f64, BemError> {
    validate_alpha(alpha)?;
    let g = geometry(foil)?;
    let n = foil.nodes.len();
    let (dim, entries) = dense_dimension(n)?;
    let (a, rhs) = assemble(foil, &g, alpha)?;
    if a.iter().chain(&rhs).any(|value| !value.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "2D adjoint assembly",
        });
    }
    let f = lu(&a, dim).map_err(|source| BemError::LinearSolve {
        stage: "primal panel",
        source,
    })?;
    let mut x = rhs;
    f.solve(&mut x);
    // Output partials by directed differences (no solves involved).
    let eps = 1e-7;
    let cl_of = |xv: &[f64], al: f64| -> f64 { output_cl(foil, &g, xv, al) };
    let mut gvec = zeroed_vec(dim, "2D lift output gradient")?;
    let mut xp = copied_vec(&x, "2D lift positive perturbation")?;
    let mut xm = copied_vec(&x, "2D lift negative perturbation")?;
    for k in 0..dim {
        xp[k] += eps;
        xm[k] -= eps;
        gvec[k] = (cl_of(&xp, alpha) - cl_of(&xm, alpha)) / (2.0 * eps);
        xp[k] = x[k];
        xm[k] = x[k];
    }
    let df_dalpha = (cl_of(&x, alpha + eps) - cl_of(&x, alpha - eps)) / (2.0 * eps);
    // Adjoint solve Aᵀλ = g.
    let mut at = zeroed_vec(entries, "2D adjoint transpose matrix")?;
    for r in 0..dim {
        for c in 0..dim {
            at[c * dim + r] = a[r * dim + c];
        }
    }
    let ft = lu(&at, dim).map_err(|source| BemError::LinearSolve {
        stage: "adjoint panel",
        source,
    })?;
    let mut lambda = gvec;
    ft.solve(&mut lambda);
    // dB/dα.
    let u_inf_d = [-det::sin(alpha), det::cos(alpha)];
    let mut db = zeroed_vec(dim, "2D panel RHS derivative")?;
    for (i, slot) in db.iter_mut().take(n).enumerate() {
        let nrm = g.normal[i];
        *slot = -(u_inf_d[0] * nrm[0] + u_inf_d[1] * nrm[1]);
    }
    for &i in &[0usize, n - 1] {
        let tan = g.tangent[i];
        db[n] -= u_inf_d[0] * tan[0] + u_inf_d[1] * tan[1];
    }
    let derivative = df_dalpha + lambda.iter().zip(&db).map(|(l, d)| l * d).sum::<f64>();
    if !derivative.is_finite() {
        return Err(BemError::NonFiniteResult {
            operation: "2D lift derivative",
        });
    }
    Ok(derivative)
}

fn validate_alpha(alpha: f64) -> Result<(), BemError> {
    if !alpha.is_finite() {
        return Err(BemError::InvalidScalar {
            name: "angle of attack",
            value: alpha,
            requirement: "finite radians",
        });
    }
    Ok(())
}

fn dense_dimension(panel_count: usize) -> Result<(usize, usize), BemError> {
    let dim = panel_count
        .checked_add(1)
        .ok_or(BemError::WorkEnvelopeExceeded {
            operation: "2D panel matrix",
            requested: usize::MAX,
            max: (MAX_AIRFOIL_PANELS + 1) * (MAX_AIRFOIL_PANELS + 1),
        })?;
    let entries = dim.checked_mul(dim).ok_or(BemError::WorkEnvelopeExceeded {
        operation: "2D panel matrix",
        requested: usize::MAX,
        max: (MAX_AIRFOIL_PANELS + 1) * (MAX_AIRFOIL_PANELS + 1),
    })?;
    let max = (MAX_AIRFOIL_PANELS + 1) * (MAX_AIRFOIL_PANELS + 1);
    if entries > max {
        return Err(BemError::WorkEnvelopeExceeded {
            operation: "2D panel matrix",
            requested: entries,
            max,
        });
    }
    Ok((dim, entries))
}

/// Assemble the Hess–Smith system (shared by solve and adjoint).
#[allow(clippy::too_many_lines)]
fn assemble(foil: &Airfoil2d, g: &Geometry, alpha: f64) -> Result<(Vec<f64>, Vec<f64>), BemError> {
    let n = foil.nodes.len();
    let u_inf = [det::cos(alpha), det::sin(alpha)];
    let dim = n + 1;
    let entries = dim.checked_mul(dim).ok_or(BemError::WorkEnvelopeExceeded {
        operation: "2D panel matrix",
        requested: usize::MAX,
        max: (MAX_AIRFOIL_PANELS + 1) * (MAX_AIRFOIL_PANELS + 1),
    })?;
    let mut a = zeroed_vec(entries, "2D panel matrix")?;
    let mut rhs = zeroed_vec(dim, "2D panel RHS")?;
    for i in 0..n {
        let (mid, nrm) = (g.mid[i], g.normal[i]);
        for j in 0..n {
            let sv = if i == j {
                [0.5 * nrm[0], 0.5 * nrm[1]]
            } else {
                source_velocity(foil, g, j, mid)
            };
            a[i * dim + j] = sv[0] * nrm[0] + sv[1] * nrm[1];
            let vv = if i == j {
                [0.0, 0.0]
            } else {
                vortex_velocity(foil, g, j, mid)
            };
            a[i * dim + n] += vv[0] * nrm[0] + vv[1] * nrm[1];
        }
        rhs[i] = -(u_inf[0] * nrm[0] + u_inf[1] * nrm[1]);
    }
    let te = [0usize, n - 1];
    for &i in &te {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        for j in 0..n {
            let sv = if i == j {
                [0.0, 0.0]
            } else {
                source_velocity(foil, g, j, mid)
            };
            a[n * dim + j] += sv[0] * tan[0] + sv[1] * tan[1];
            let vv = if i == j {
                [0.5 * tan[0], 0.5 * tan[1]]
            } else {
                vortex_velocity(foil, g, j, mid)
            };
            a[n * dim + n] += vv[0] * tan[0] + vv[1] * tan[1];
        }
        rhs[n] -= u_inf[0] * tan[0] + u_inf[1] * tan[1];
    }
    Ok((a, rhs))
}

/// Pressure-integrated Cl as a pure function of (solution, α).
fn output_cl(foil: &Airfoil2d, g: &Geometry, x: &[f64], alpha: f64) -> f64 {
    let n = foil.nodes.len();
    let u_inf = [det::cos(alpha), det::sin(alpha)];
    let gamma = x[n];
    let mut fx = 0.0;
    let mut fy = 0.0;
    for i in 0..n {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        let mut v = [u_inf[0], u_inf[1]];
        for (j, &xj) in x.iter().take(n).enumerate() {
            if i == j {
                v[0] += 0.5 * gamma * tan[0];
                v[1] += 0.5 * gamma * tan[1];
                continue;
            }
            let sv = source_velocity(foil, g, j, mid);
            let vv = vortex_velocity(foil, g, j, mid);
            v[0] += sv[0] * xj + vv[0] * gamma;
            v[1] += sv[1] * xj + vv[1] * gamma;
        }
        let vt = v[0] * tan[0] + v[1] * tan[1];
        let cp = 1.0 - vt * vt;
        fx += -cp * g.normal[i][0] * g.len[i];
        fy += -cp * g.normal[i][1] * g.len[i];
    }
    fy * det::cos(alpha) - fx * det::sin(alpha)
}
