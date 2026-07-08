//! 2D Hess–Smith airfoil panel method: constant source density per
//! panel plus ONE bound vortex density shared by all panels, closed by
//! the KUTTA condition (equal tangential speeds leaving the two
//! trailing-edge panels — smooth flow off the edge; circulation is
//! DETERMINED, not assumed, which is the §8.3 tie to the harmonic-
//! circulation story). Dense solve at section scale; the ADJOINT of
//! the lift coefficient is one transposed solve, FD-gated in the
//! battery. Inviscid screening honesty label applies.

use fs_la::factor::lu;

/// A closed airfoil section: nodes ordered clockwise from the trailing
/// edge along the lower surface and back over the upper surface.
pub struct Airfoil2d {
    /// Panel endpoints (closed: node N = node 0 implicitly).
    pub nodes: Vec<[f64; 2]>,
}

/// The solved panel state.
pub struct PanelSolution2d {
    /// Source densities per panel.
    pub sources: Vec<f64>,
    /// The shared vortex density.
    pub gamma: f64,
    /// Lift coefficient (Kutta–Joukowski from total circulation).
    pub cl: f64,
    /// Tangential speeds at panel midpoints.
    pub vt: Vec<f64>,
}

/// Symmetric NACA 4-digit thickness form (e.g. t = 0.12 for 0012),
/// `n` panels, closed trailing edge, unit chord.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn naca4_symmetric(t: f64, n: usize) -> Airfoil2d {
    assert!(n >= 8 && n % 2 == 0, "even panel count >= 8");
    let half = n / 2;
    let thick = |x: f64| {
        5.0 * t
            * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
                - 0.1036 * x * x * x * x)
    };
    let mut nodes = Vec::with_capacity(n);
    // Cosine-clustered x from TE (1) along the LOWER surface to LE (0).
    for k in 0..half {
        let x = 0.5 * (1.0 + (std::f64::consts::PI * k as f64 / half as f64).cos());
        nodes.push([x, -thick(x)]);
    }
    // LE to TE along the UPPER surface.
    for k in 0..half {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / half as f64).cos());
        nodes.push([x, thick(x)]);
    }
    Airfoil2d { nodes }
}

struct Geometry {
    mid: Vec<[f64; 2]>,
    tangent: Vec<[f64; 2]>,
    normal: Vec<[f64; 2]>,
    len: Vec<f64>,
}

fn geometry(foil: &Airfoil2d) -> Geometry {
    let n = foil.nodes.len();
    let mut g = Geometry {
        mid: Vec::with_capacity(n),
        tangent: Vec::with_capacity(n),
        normal: Vec::with_capacity(n),
        len: Vec::with_capacity(n),
    };
    for i in 0..n {
        let a = foil.nodes[i];
        let b = foil.nodes[(i + 1) % n];
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let l = dx.hypot(dy);
        g.mid.push([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]);
        g.tangent.push([dx / l, dy / l]);
        g.normal.push([-dy / l, dx / l]);
        g.len.push(l);
    }
    g
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
    let r1 = xl.hypot(yl).max(1e-12);
    let r2 = (xl - l).hypot(yl).max(1e-12);
    let theta1 = yl.atan2(xl);
    let theta2 = yl.atan2(xl - l);
    let two_pi = std::f64::consts::TAU;
    // v_n = (θ₂ − θ₁)/2π — the panel-probe battery pinned the
    // reversed order (a source panel must push AWAY on both sides).
    let ul = (r1 / r2).ln() / two_pi;
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
/// # Panics
/// On singular geometry (degenerate panels).
#[must_use]
#[allow(clippy::too_many_lines)] // assembly + Kutta + postprocess, one narrative
pub fn solve(foil: &Airfoil2d, alpha: f64) -> PanelSolution2d {
    let g = geometry(foil);
    let n = foil.nodes.len();
    let u_inf = [alpha.cos(), alpha.sin()];
    // Unknowns: n sources + 1 vortex density. Equations: n
    // no-penetration + 1 Kutta.
    let dim = n + 1;
    let mut a = vec![0.0f64; dim * dim];
    let mut rhs = vec![0.0f64; dim];
    for i in 0..n {
        let (mid, nrm) = (g.mid[i], g.normal[i]);
        for j in 0..n {
            let sv = if i == j {
                // Self-induction of a constant source panel: σ/2 along
                // the normal.
                [0.5 * nrm[0], 0.5 * nrm[1]]
            } else {
                source_velocity(foil, &g, j, mid)
            };
            a[i * dim + j] = sv[0] * nrm[0] + sv[1] * nrm[1];
            let vv = if i == j {
                // Vortex self-term: ±γ/2 along the tangent (no normal
                // component at the panel's own midpoint).
                [0.0, 0.0]
            } else {
                vortex_velocity(foil, &g, j, mid)
            };
            a[i * dim + n] += vv[0] * nrm[0] + vv[1] * nrm[1];
        }
        rhs[i] = -(u_inf[0] * nrm[0] + u_inf[1] * nrm[1]);
    }
    // Kutta: tangential speeds off the two TE panels match:
    // V_t(panel 0) + V_t(panel n−1) = 0 with tangents oriented along
    // the marching direction (lower TE→LE, upper LE→TE).
    let te = [0usize, n - 1];
    for &i in &te {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        for j in 0..n {
            let sv = if i == j {
                [0.0, 0.0]
            } else {
                source_velocity(foil, &g, j, mid)
            };
            a[n * dim + j] += sv[0] * tan[0] + sv[1] * tan[1];
            let vv = if i == j {
                // Vortex sheet's own tangential jump: γ/2.
                [0.5 * tan[0], 0.5 * tan[1]]
            } else {
                vortex_velocity(foil, &g, j, mid)
            };
            a[n * dim + n] += vv[0] * tan[0] + vv[1] * tan[1];
        }
        rhs[n] -= u_inf[0] * tan[0] + u_inf[1] * tan[1];
    }
    let f = lu(&a, dim).expect("panel system is nonsingular");
    let mut x = rhs.clone();
    f.solve(&mut x);
    let gamma = x[n];
    // Tangential speeds for Cp.
    let mut vt = Vec::with_capacity(n);
    for i in 0..n {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        let mut v = [u_inf[0], u_inf[1]];
        for j in 0..n {
            if i == j {
                v[0] += 0.5 * gamma * tan[0];
                v[1] += 0.5 * gamma * tan[1];
                continue;
            }
            let sv = source_velocity(foil, &g, j, mid);
            let vv = vortex_velocity(foil, &g, j, mid);
            v[0] += sv[0] * x[j] + vv[0] * gamma;
            v[1] += sv[1] * x[j] + vv[1] * gamma;
        }
        vt.push(v[0] * tan[0] + v[1] * tan[1]);
    }
    // Lift by PRESSURE INTEGRATION (v_n = 0 on the body, so
    // |v|² = v_t²): Cl = −∮ Cp·n_y ds, then rotate the force normal to
    // the freestream. Trustworthy bookkeeping — the surface field is
    // what the collocation system actually enforces.
    let mut fx = 0.0;
    let mut fy = 0.0;
    for i in 0..n {
        let cp = 1.0 - vt[i] * vt[i];
        fx += -cp * g.normal[i][0] * g.len[i];
        fy += -cp * g.normal[i][1] * g.len[i];
    }
    let cl = fy * alpha.cos() - fx * alpha.sin();
    PanelSolution2d {
        sources: x[..n].to_vec(),
        gamma,
        cl,
        vt,
    }
}

/// dCl/dα by the ADJOINT of the panel system: one transposed solve
/// for the solution sensitivity (the committed gradient path); the
/// output function's own partials (`∂f/∂x`, `∂f/∂α` at FIXED
/// solution) are cheap solve-free differences. FD-gated end to end.
#[must_use]
pub fn dcl_dalpha_adjoint(foil: &Airfoil2d, alpha: f64) -> f64 {
    let g = geometry(foil);
    let n = foil.nodes.len();
    let dim = n + 1;
    let (a, rhs) = assemble(foil, &g, alpha);
    let f = lu(&a, dim).expect("panel system nonsingular");
    let mut x = rhs.clone();
    f.solve(&mut x);
    // Output partials by directed differences (no solves involved).
    let eps = 1e-7;
    let cl_of = |xv: &[f64], al: f64| -> f64 { output_cl(foil, &g, xv, al) };
    let mut gvec = vec![0.0f64; dim];
    for k in 0..dim {
        let mut xp = x.clone();
        xp[k] += eps;
        let mut xm = x.clone();
        xm[k] -= eps;
        gvec[k] = (cl_of(&xp, alpha) - cl_of(&xm, alpha)) / (2.0 * eps);
    }
    let df_dalpha = (cl_of(&x, alpha + eps) - cl_of(&x, alpha - eps)) / (2.0 * eps);
    // Adjoint solve Aᵀλ = g.
    let mut at = vec![0.0f64; dim * dim];
    for r in 0..dim {
        for c in 0..dim {
            at[c * dim + r] = a[r * dim + c];
        }
    }
    let ft = lu(&at, dim).expect("adjoint system nonsingular");
    let mut lambda = gvec;
    ft.solve(&mut lambda);
    // dB/dα.
    let u_inf_d = [-alpha.sin(), alpha.cos()];
    let mut db = vec![0.0f64; dim];
    for i in 0..n {
        let nrm = g.normal[i];
        db[i] = -(u_inf_d[0] * nrm[0] + u_inf_d[1] * nrm[1]);
    }
    for &i in &[0usize, n - 1] {
        let tan = g.tangent[i];
        db[n] -= u_inf_d[0] * tan[0] + u_inf_d[1] * tan[1];
    }
    df_dalpha + lambda.iter().zip(&db).map(|(l, d)| l * d).sum::<f64>()
}

/// Assemble the Hess–Smith system (shared by solve and adjoint).
#[allow(clippy::too_many_lines)]
fn assemble(foil: &Airfoil2d, g: &Geometry, alpha: f64) -> (Vec<f64>, Vec<f64>) {
    let n = foil.nodes.len();
    let u_inf = [alpha.cos(), alpha.sin()];
    let dim = n + 1;
    let mut a = vec![0.0f64; dim * dim];
    let mut rhs = vec![0.0f64; dim];
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
    (a, rhs)
}

/// Pressure-integrated Cl as a pure function of (solution, α).
fn output_cl(foil: &Airfoil2d, g: &Geometry, x: &[f64], alpha: f64) -> f64 {
    let n = foil.nodes.len();
    let u_inf = [alpha.cos(), alpha.sin()];
    let gamma = x[n];
    let mut fx = 0.0;
    let mut fy = 0.0;
    for i in 0..n {
        let (mid, tan) = (g.mid[i], g.tangent[i]);
        let mut v = [u_inf[0], u_inf[1]];
        for j in 0..n {
            if i == j {
                v[0] += 0.5 * gamma * tan[0];
                v[1] += 0.5 * gamma * tan[1];
                continue;
            }
            let sv = source_velocity(foil, g, j, mid);
            let vv = vortex_velocity(foil, g, j, mid);
            v[0] += sv[0] * x[j] + vv[0] * gamma;
            v[1] += sv[1] * x[j] + vv[1] * gamma;
        }
        let vt = v[0] * tan[0] + v[1] * tan[1];
        let cp = 1.0 - vt * vt;
        fx += -cp * g.normal[i][0] * g.len[i];
        fy += -cp * g.normal[i][1] * g.len[i];
    }
    fy * alpha.cos() - fx * alpha.sin()
}

