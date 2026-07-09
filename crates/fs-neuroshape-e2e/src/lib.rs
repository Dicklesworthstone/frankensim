//! fs-neuroshape-e2e — NeuroShapeCert: a PROVEN neural implicit shape.
//! Layer: L5 (LUMEN).
//!
//! # The campaign
//!
//! A learned neural SDF renders a shape, but gives no guarantees: how far can a
//! sphere-tracing ray step without tunneling through a thin feature? How many
//! components does the level set actually have? This proves those, composing
//! crates never designed to meet:
//!
//! - **The field** ([`fs_rep_neural`]): a small `tanh`-MLP SDF whose
//!   spectral-normalized effective form is `≈ 2.12·Σ tanh(3(±coord − 0.7)) + 6.5`
//!   — provably negative near the origin, provably positive on a surrounding ring.
//! - **A certified Lipschitz constant** — `L = Π σᵢ` (product of spectral norms).
//!   Then `safe_step_radius(f, L) = |f|/L` is a certified sphere-trace step that
//!   CANNOT tunnel through the surface.
//! - **A topology certificate by interval arithmetic** — the network's sound
//!   Interval Bound Propagation (`eval_interval`) proves a central box is
//!   strictly inside (`hi < 0`) and that the FOUR edge strips of a bounding box
//!   are strictly outside (`lo > 0`). Those strips tile the box boundary into a
//!   CLOSED frame (corners overlap), so `{f<0}` provably cannot cross it: the
//!   interior is proven NON-EMPTY and BOUNDED — a proof, not a mesh. (Discrete
//!   ring boxes would leave angular gaps and prove nothing about boundedness.)
//! - **A Morse cross-check** ([`fs_viz`]): the field has a single interior
//!   minimum (`classify_hessian → Minimum`) — evidence (not proof) that the
//!   bounded region is a single component; `isocontour_crossings` localizes the
//!   zero set, all inside the certified frame.
//! - **Honest colors** ([`fs_evidence`]): every certificate is `Verified`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::Color;
use fs_rep_neural::{Layer, MlpSdf, safe_step_radius};
use fs_viz::{CriticalKind, Grid2, Vec2, classify_hessian};

/// The blob SDF network. `MlpSdf::new` spectral-normalizes every layer to
/// exactly `bound`, so with `bound = √18` the effective hidden slope is
/// `√18/√2 = 3` (a wall at `|coord| = 0.7`) and the effective output weight is
/// `√18/2 ≈ 2.12`. With the biases below the effective field is
/// `f ≈ 2.12·Σ tanh(3(±coord − 0.7)) + 6.5`: negative near the origin, positive
/// on a surrounding ring. `L = bound² = 18`.
#[must_use]
pub fn blob_sdf_net() -> MlpSdf {
    // Hidden layer: one tanh wall per ±axis direction (bias −2.1 ⇒ wall at 0.7).
    let l1 = Layer::new(
        vec![
            vec![3.0, 0.0],
            vec![-3.0, 0.0],
            vec![0.0, 3.0],
            vec![0.0, -3.0],
        ],
        vec![-2.1, -2.1, -2.1, -2.1],
    );
    // Linear output: sum the walls, lift by +6.5 (bias is not normalized).
    let l2 = Layer::new(vec![vec![1.0, 1.0, 1.0, 1.0]], vec![6.5]);
    MlpSdf::new(vec![l1, l2], (18.0_f64).sqrt())
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct NeuroShapeReport {
    /// The certified global Lipschitz constant `L`.
    pub lipschitz: f64,
    /// The field value at the origin.
    pub origin_value: f64,
    /// A certified no-tunnel sphere-trace step at the origin (`|f|/L`).
    pub safe_radius: f64,
    /// IBP enclosure of `f` over the central box.
    pub inside_interval: (f64, f64),
    /// Is the central box certified strictly inside (`hi < 0`)?
    pub certified_inside: bool,
    /// How many of the box-boundary strips are certified strictly outside.
    pub boundary_certified: usize,
    /// Total boundary strips (4 — a CLOSED frame around the box).
    pub boundary_segments: usize,
    /// Is the surface bounded — the whole boundary frame certified outside, so
    /// `{f<0}` provably cannot cross it (a closed barrier, not spot checks)?
    pub bounded: bool,
    /// Morse: does the field have a single interior minimum?
    pub single_minimum: bool,
    /// Number of zero-set crossings found on the visualization grid.
    pub surface_crossings: usize,
    /// The largest radius at which a crossing was found (must be inside the ring).
    pub max_crossing_radius: f64,
    /// The smallest radius at which a crossing was found — the NEAREST surface
    /// point; the safe step radius must under-estimate it (no-tunnel soundness).
    pub nearest_surface_radius: f64,
    /// The topology certificate color (`Verified` — IBP is sound).
    pub topology_color: Color,
}

fn radius(p: Vec2) -> f64 {
    p[0].hypot(p[1])
}

/// Run the NeuroShapeCert campaign on `net` with a bounding box of half-width
/// `ring_r` (its four edge strips form the closed barrier) and a central
/// certified-inside box of half-width `inner`.
#[must_use]
pub fn run_campaign(net: &MlpSdf, ring_r: f64, inner: f64) -> NeuroShapeReport {
    let lipschitz = net.lipschitz();
    let origin_value = net.eval(&[0.0, 0.0]);
    let safe_radius = safe_step_radius(origin_value, lipschitz);

    // Interval topology certificate.
    let inside_interval = net.eval_interval(&[-inner, -inner], &[inner, inner]);
    let certified_inside = inside_interval.1 < 0.0;

    // A CLOSED barrier: the four edge strips of the box [−R, R]² tile the whole
    // boundary frame (corners overlap), so certifying every strip strictly
    // outside (lo > 0) RIGOROUSLY traps {f<0} inside the box. Eight discrete
    // boxes would leave angular gaps the surface could escape through.
    let r = ring_r;
    let w = 0.4;
    let strips = [
        ([-r, r - w], [r, r]),   // top
        ([-r, -r], [r, -r + w]), // bottom
        ([-r, -r], [-r + w, r]), // left
        ([r - w, -r], [r, r]),   // right
    ];
    let boundary_segments = strips.len();
    let mut boundary_certified = 0usize;
    for (lo_pt, hi_pt) in strips {
        let (lo, _hi) = net.eval_interval(&lo_pt, &hi_pt);
        if lo > 0.0 {
            boundary_certified += 1;
        }
    }
    let bounded = boundary_certified == boundary_segments;

    // Morse cross-check: a single interior minimum (Hessian by finite diff).
    let h = 1e-3;
    let f00 = origin_value;
    let fxx = (net.eval(&[h, 0.0]) - 2.0 * f00 + net.eval(&[-h, 0.0])) / (h * h);
    let fyy = (net.eval(&[0.0, h]) - 2.0 * f00 + net.eval(&[0.0, -h])) / (h * h);
    let fxy = (net.eval(&[h, h]) - net.eval(&[h, -h]) - net.eval(&[-h, h]) + net.eval(&[-h, -h]))
        / (4.0 * h * h);
    let crit = classify_hessian([[fxx, fxy], [fxy, fyy]], 1e-6);
    let single_minimum = crit.kind == CriticalKind::Minimum;

    // Localize the zero set on a visualization grid.
    let grid = Grid2::from_fn(
        81,
        81,
        [-ring_r - 0.5, -ring_r - 0.5],
        [ring_r + 0.5, ring_r + 0.5],
        |p| net.eval(&[p[0], p[1]]),
    );
    let crossings = grid.isocontour_crossings(0.0);
    let max_crossing_radius = crossings.iter().copied().map(radius).fold(0.0, f64::max);
    let nearest_surface_radius = crossings
        .iter()
        .copied()
        .map(radius)
        .fold(f64::INFINITY, f64::min);

    // The topology claim is Verified iff the interval certificate closed: a
    // non-empty interior enclosed by the certified-positive boundary frame. The
    // certified containment is the frame's INNER edge `ring_r − w` (max-norm).
    let topology_color = if certified_inside && bounded {
        Color::Verified {
            lo: 0.0,
            hi: ring_r - w,
        }
    } else {
        Color::Estimated {
            estimator: "ibp-open".to_string(),
            dispersion: f64::INFINITY,
        }
    };

    NeuroShapeReport {
        lipschitz,
        origin_value,
        safe_radius,
        inside_interval,
        certified_inside,
        boundary_certified,
        boundary_segments,
        bounded,
        single_minimum,
        surface_crossings: crossings.len(),
        max_crossing_radius,
        nearest_surface_radius,
        topology_color,
    }
}
