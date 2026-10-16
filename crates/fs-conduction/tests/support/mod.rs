//! Shared harness for the fs-conduction batteries: a deterministic `Cx`,
//! an EXACT tetrahedral quadrature for error norms, and the manufactured
//! solutions the G1 ladders and analytic cases use.
//!
//! ## Why the quadrature has to be this accurate
//!
//! The L2 error of a P₁ solve is `O(h²)`, so `∫(T_h − T)² dV` is `O(h⁴)`.
//! A degree-3 rule integrates that with an `O(1)` RELATIVE error (a 1-D
//! Simpson check on the exact interpolation remainder gives 25% low), so
//! the "measured" order would be an artefact of the rule. This harness
//! uses a Duffy-transformed 6×6×6 Gauss–Legendre product rule. The Duffy
//! Jacobian `(1−u)²(1−v)` raises a degree-`d` tet polynomial to degree
//! `d+2` in `u`, and 6-point Gauss is exact to degree 11, so the rule is
//! EXACT for tet polynomials of total degree ≤ 9. The manufactured
//! solutions below top out at degree 4, so `(T_h − T)²` has degree ≤ 8
//! and the error norms are computed EXACTLY — no quadrature term enters
//! the fitted order at all.
//!
//! ## Why one of the solutions is QUARTIC
//!
//! P₁ Galerkin on these structured Kuhn meshes reproduces exact
//! solutions up to CUBIC at the nodes — measured, not assumed; see
//! `conformance.rs::polynomial_reproduction_at_the_nodes`. The 1-D
//! mechanism is the textbook one: the second difference
//! `−u_{i−1} + 2u_i − u_{i+1}` equals `−h²u'' − h⁴u''''/12`, which is
//! exact whenever `u'''' = 0` and the load is integrated exactly.
//!
//! A ladder run on a quadratic or cubic solution therefore measures the
//! INTERPOLATION error `‖T − I_h T‖`, not the scheme's own approximation
//! error — a number that does not depend on `K` at all. Such a ladder
//! still gates order 2 and still catches a mis-assembled operator (a
//! wrong `K` destroys the nodal-reproduction identity immediately), but
//! it is weaker evidence than it looks. [`Quartic`] has `u'''' ≠ 0`, so
//! the Dirichlet ladders fit a genuinely non-degenerate discretization
//! error, and the conformance suite pins the reproduction fact itself so
//! the degeneracy is recorded rather than mistaken for strength.

#![allow(dead_code)]

use fs_conduction::mesh::ConductionMesh;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

/// Run `f` under a deterministic `Cx`.
pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(|_, cx| f(cx))
}

/// Run `f` with access to the cancellation gate, so a test can request
/// cancellation mid-flight.
pub fn with_gate<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C0DE_0C71_0000,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

/// A pre-cancelled `Cx`.
pub fn with_cancelled_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(|gate, cx| {
        gate.request();
        f(cx)
    })
}

/// Duffy-transformed Gauss–Legendre points on the reference tet, as
/// `(barycentric, weight)` with weights summing to `1/6` (the reference
/// tet's volume). Exact for tet polynomials of total degree ≤ 9.
pub fn tet_quadrature() -> Vec<([f64; 4], f64)> {
    const N: usize = 6;
    let (nodes, weights) = fs_feec::gauss_legendre(N);
    // Map [-1, 1] → [0, 1].
    let unit: Vec<(f64, f64)> = nodes
        .iter()
        .zip(&weights)
        .map(|(&x, &w)| (f64::midpoint(x, 1.0), w / 2.0))
        .collect();
    let mut out = Vec::with_capacity(N * N * N);
    for &(u, wu) in &unit {
        for &(v, wv) in &unit {
            for &(w, ww) in &unit {
                let l1 = u;
                let l2 = v * (1.0 - u);
                let l3 = w * (1.0 - u) * (1.0 - v);
                let l0 = 1.0 - l1 - l2 - l3;
                let jac = (1.0 - u) * (1.0 - u) * (1.0 - v);
                out.push(([l0, l1, l2, l3], wu * wv * ww * jac));
            }
        }
    }
    out
}

/// `‖T_h − T‖_{L2(Ω)}`, computed with [`tet_quadrature`].
pub fn l2_error(mesh: &ConductionMesh, nodal: &[f64], exact: &dyn Fn([f64; 3]) -> f64) -> f64 {
    let rule = tet_quadrature();
    let mut acc = 0.0f64;
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        let volume = mesh.element_volume(e);
        let p: [[f64; 3]; 4] = [
            mesh.positions()[tet[0] as usize],
            mesh.positions()[tet[1] as usize],
            mesh.positions()[tet[2] as usize],
            mesh.positions()[tet[3] as usize],
        ];
        let t: [f64; 4] = [
            nodal[tet[0] as usize],
            nodal[tet[1] as usize],
            nodal[tet[2] as usize],
            nodal[tet[3] as usize],
        ];
        let mut local = 0.0f64;
        for (bary, w) in &rule {
            let mut x = [0.0f64; 3];
            let mut th = 0.0f64;
            for a in 0..4 {
                for k in 0..3 {
                    x[k] += bary[a] * p[a][k];
                }
                th += bary[a] * t[a];
            }
            let d = th - exact(x);
            local += w * d * d;
        }
        // The rule's weights sum to the REFERENCE volume 1/6, so scale
        // by 6·V to reach the physical element.
        acc += 6.0 * volume * local;
    }
    fs_math::det::sqrt(acc)
}

/// `|T_h − T|_{H1(Ω)}` — the gradient semi-norm. `∇T_h` is constant per
/// element, so only the exact gradient is sampled.
pub fn h1_error(
    mesh: &ConductionMesh,
    nodal: &[f64],
    exact_grad: &dyn Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    let rule = tet_quadrature();
    let mut acc = 0.0f64;
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        let volume = mesh.element_volume(e);
        let g = &mesh.geometry().grads[e];
        let p: [[f64; 3]; 4] = [
            mesh.positions()[tet[0] as usize],
            mesh.positions()[tet[1] as usize],
            mesh.positions()[tet[2] as usize],
            mesh.positions()[tet[3] as usize],
        ];
        let mut grad_h = [0.0f64; 3];
        for a in 0..4 {
            let ta = nodal[tet[a] as usize];
            for k in 0..3 {
                grad_h[k] += ta * g[a][k];
            }
        }
        let mut local = 0.0f64;
        for (bary, w) in &rule {
            let mut x = [0.0f64; 3];
            for a in 0..4 {
                for k in 0..3 {
                    x[k] += bary[a] * p[a][k];
                }
            }
            let ge = exact_grad(x);
            let mut d2 = 0.0f64;
            for k in 0..3 {
                let d = grad_h[k] - ge[k];
                d2 += d * d;
            }
            local += w * d2;
        }
        acc += 6.0 * volume * local;
    }
    fs_math::det::sqrt(acc)
}

/// Escape a detail string for the JSON-lines verdicts. Refusal messages
/// quote region names, so an unescaped detail produces a log line no
/// tool can parse — and CONVENTIONS requires every failure to be
/// reproducible from its log line alone.
pub fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push(' '),
            c => out.push(c),
        }
    }
    out
}

/// The max nodal deviation from an exact solution.
pub fn max_nodal_error(
    mesh: &ConductionMesh,
    nodal: &[f64],
    exact: &dyn Fn([f64; 3]) -> f64,
) -> f64 {
    mesh.positions()
        .iter()
        .zip(nodal)
        .map(|(&p, &t)| (t - exact(p)).abs())
        .fold(0.0f64, f64::max)
}

// ---------------------------------------------------------------------
// Manufactured solutions (all QUADRATIC, so the error quadrature above
// is exact and no transcendental enters any battery).
// ---------------------------------------------------------------------

/// The full quadratic with every cross term: exercises all six second
/// derivatives, so an anisotropic tensor's off-diagonal entries appear
/// in the manufactured source.
pub struct FullQuadratic;

impl FullQuadratic {
    pub fn value(p: [f64; 3]) -> f64 {
        let [x, y, z] = p;
        300.0
            + 2.0 * x
            + 3.0 * y
            + 4.0 * z
            + 5.0 * x * x
            + 3.0 * y * y
            + 2.0 * z * z
            + x * y
            + y * z
            + z * x
    }

    pub fn gradient(p: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = p;
        [
            2.0 + 10.0 * x + y + z,
            3.0 + 6.0 * y + x + z,
            4.0 + 4.0 * z + y + x,
        ]
    }

    /// The constant Hessian `[∂xx, ∂yy, ∂zz, ∂xy, ∂xz, ∂yz]`.
    pub const HESSIAN: [f64; 6] = [10.0, 6.0, 4.0, 1.0, 1.0, 1.0];

    /// `f = −∇·(K∇T)` for a constant tensor `K`.
    pub fn source(k: [[f64; 3]; 3]) -> f64 {
        let [hxx, hyy, hzz, hxy, hxz, hyz] = Self::HESSIAN;
        -(k[0][0] * hxx
            + k[1][1] * hyy
            + k[2][2] * hzz
            + 2.0 * (k[0][1] * hxy + k[0][2] * hxz + k[1][2] * hyz))
    }
}

/// A CUBIC manufactured solution. P₁ Galerkin does NOT reproduce it at
/// the nodes, so a ladder run on it fits the scheme's own discretization
/// error. Its Hessian is linear, so `f = −∇·(K∇T)` is linear and the
/// nodal `V(1+δ)/20` load rule stays EXACT — the ladder measures the
/// discretization and nothing else.
pub struct Cubic;

impl Cubic {
    pub fn value(p: [f64; 3]) -> f64 {
        let [x, y, z] = p;
        300.0 + 2.0 * x + 3.0 * y + 4.0 * z + x * x * y + y * y * z + z * z * x + x * x * x
    }

    pub fn gradient(p: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = p;
        [
            2.0 + 2.0 * x * y + z * z + 3.0 * x * x,
            3.0 + x * x + 2.0 * y * z,
            4.0 + y * y + 2.0 * z * x,
        ]
    }

    /// `[∂xx, ∂yy, ∂zz, ∂xy, ∂xz, ∂yz]` — linear in `p`.
    pub fn hessian(p: [f64; 3]) -> [f64; 6] {
        let [x, y, z] = p;
        [
            2.0f64.mul_add(y, 6.0 * x),
            2.0 * z,
            2.0 * x,
            2.0 * x,
            2.0 * z,
            2.0 * y,
        ]
    }

    /// `f = −∇·(K∇T)` for a constant tensor `K` — linear in `p`.
    pub fn source(k: [[f64; 3]; 3], p: [f64; 3]) -> f64 {
        let [hxx, hyy, hzz, hxy, hxz, hyz] = Self::hessian(p);
        -(k[0][0] * hxx
            + k[1][1] * hyy
            + k[2][2] * hzz
            + 2.0 * (k[0][1] * hxy + k[0][2] * hxz + k[1][2] * hyz))
    }
}

/// A QUARTIC manufactured solution: `u'''' ≠ 0`, so P₁ Galerkin does
/// NOT reproduce it at the nodes and a ladder run on it fits the
/// scheme's own discretization error rather than the interpolation
/// error. Its Hessian is quadratic, so `f = −∇·(K∇T)` is quadratic and
/// the nodal `V(1+δ)/20` load rule is second-order rather than exact —
/// an `O(h²)` consistency term, the same order as the discretization
/// error it sits beside, so the fitted order stays 2 while the constant
/// changes. That is stated rather than hidden.
pub struct Quartic;

impl Quartic {
    pub fn value(p: [f64; 3]) -> f64 {
        let [x, y, z] = p;
        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;
        300.0 + 2.0 * x + 3.0 * y + 4.0 * z + x2 * x2 + y2 * y2 + z2 * z2 + x2 * y2
    }

    pub fn gradient(p: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = p;
        [
            2.0 + 4.0 * x * x * x + 2.0 * x * y * y,
            3.0 + 4.0 * y * y * y + 2.0 * x * x * y,
            4.0 + 4.0 * z * z * z,
        ]
    }

    /// `[∂xx, ∂yy, ∂zz, ∂xy, ∂xz, ∂yz]` — quadratic in `p`.
    pub fn hessian(p: [f64; 3]) -> [f64; 6] {
        let [x, y, z] = p;
        [
            12.0f64.mul_add(x * x, 2.0 * y * y),
            12.0f64.mul_add(y * y, 2.0 * x * x),
            12.0 * z * z,
            4.0 * x * y,
            0.0,
            0.0,
        ]
    }

    /// `f = −∇·(K∇T)` for a constant tensor `K` — quadratic in `p`.
    pub fn source(k: [[f64; 3]; 3], p: [f64; 3]) -> f64 {
        let [hxx, hyy, hzz, hxy, hxz, hyz] = Self::hessian(p);
        -(k[0][0] * hxx
            + k[1][1] * hyy
            + k[2][2] * hzz
            + 2.0 * (k[0][1] * hxy + k[0][2] * hxz + k[1][2] * hyz))
    }
}

/// A quadratic that is LINEAR on the `x = 0` face, so a Robin reference
/// temperature `T_ref = T − q_n/h` is linear there and the boundary data
/// enters the P₁ trace space exactly — the Robin ladder then measures
/// the discretization, not the interpolation of its own data.
pub struct FaceLinearQuadratic;

impl FaceLinearQuadratic {
    pub fn value(p: [f64; 3]) -> f64 {
        let [x, y, z] = p;
        300.0 + 2.0 * x + 3.0 * y + 4.0 * z + 5.0 * x * x + x * y + x * z
    }

    pub fn gradient(p: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = p;
        [2.0 + 10.0 * x + y + z, 3.0 + x, 4.0 + x]
    }

    /// `Δ` of this solution (constant).
    pub const LAPLACIAN: f64 = 10.0;
}
