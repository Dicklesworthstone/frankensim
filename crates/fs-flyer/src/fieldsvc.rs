//! Field-sampling service (bead wf-root-guzez.8.1.1, E7.1-i, plan
//! §5.5). `sample_field(state, grid_spec, component_mask)` returns
//! the full §5.5 kinematic set: u, grad_u, omega, the
//! analytic-vs-finite-difference divergence DUAL, strain magnitude,
//! Q-criterion, lambda2, kinetic speed gradient, validity /
//! singularity-core / solid-exclusion masks, provenance, and meta.
//!
//! Doctrine:
//!   - NO volume pressure (Round-1 P0): a kinematically synthesized
//!     solenoidal field does not determine pressure — no `grad_p_hat`
//!     exists anywhere in this API.
//!   - Omitted-component HONESTY: `meta.omitted_components` names
//!     every unselected component, and [`claim_total_flow`] refuses
//!     whenever a force-coupled component the active model supports
//!     is absent from the sum.
//!   - Ambient components are analytic (fs-atmo exact gradients); the
//!     bound-circulation and ground-image components use the local
//!     horseshoe kernel with a declared singularity core. The
//!     physical-wake (fs-vpm hybrid) and propeller-induced branches
//!     are follow-up wiring — they are DECLARED unsupported by this
//!     state version, so their omission is named but does not forbid
//!     the total-flow claim.

use crate::{Refusal, refuse};
use fs_atmo::Atmosphere;
use fs_blake3::hash_domain;
use fs_math::det;

/// Component bit: mean atmosphere (log law).
pub const C_MEAN_ATMO: u32 = 1 << 0;
/// Component bit: turbulent atmosphere.
pub const C_TURB_ATMO: u32 = 1 << 1;
/// Component bit: gust event (unsupported in state v1).
pub const C_GUST_EVENT: u32 = 1 << 2;
/// Component bit: bound circulation (horseshoe).
pub const C_BOUND_CIRCULATION: u32 = 1 << 3;
/// Component bit: physical wake (fs-vpm hybrid; unsupported in v1).
pub const C_PHYSICAL_WAKE: u32 = 1 << 4;
/// Component bit: ground images of the bound system.
pub const C_GROUND_IMAGES: u32 = 1 << 5;
/// Component bit: propeller induced field (unsupported in v1).
pub const C_PROP_INDUCED: u32 = 1 << 6;
/// Component bit: visualization-only embellishments (never physics).
pub const C_VIS_ONLY: u32 = 1 << 7;

/// All defined bits.
pub const C_ALL_DEFINED: u32 = 0xff;

/// Names, bit-order aligned.
pub const COMPONENT_NAMES: [&str; 8] = [
    "mean-atmosphere",
    "turbulent-atmosphere",
    "gust-event",
    "bound-circulation",
    "physical-wake",
    "ground-images",
    "propeller-induced",
    "visualization-only",
];

/// Force-coupled components (plan: body displacement flow, propeller
/// slipstream, viscous wakes, and the bound system that carries lift).
pub const FORCE_COUPLED: u32 = C_BOUND_CIRCULATION | C_PHYSICAL_WAKE | C_PROP_INDUCED;

/// Point cap per request (browser-tier budget).
pub const MAX_POINTS: usize = 262_144;

/// FD probe step [m].
const FD_H: f64 = 5.0e-4;

/// Gradient-norm floor under which normalized divergence is masked.
pub const GRAD_NORM_FLOOR: f64 = 1.0e-6;

/// The sampling grid (uniform, axis-aligned; z is height above the
/// certified plane).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridSpec {
    /// Southwest-bottom corner [m].
    pub origin_m: [f64; 3],
    /// Uniform spacing [m].
    pub dx_m: f64,
    /// Points along x.
    pub nx: usize,
    /// Points along y.
    pub ny: usize,
    /// Points along z.
    pub nz: usize,
}

impl GridSpec {
    /// Total points.
    #[must_use]
    pub fn n_points(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Admit the grid.
    ///
    /// # Errors
    /// `grid-spec-invalid` (zero extent, non-finite, dx ≤ 0, or above
    /// the point cap — AT the cap admits, one more point refuses).
    pub fn admit(&self) -> Result<(), Refusal> {
        let finite = self.origin_m.iter().all(|v| v.is_finite()) && self.dx_m.is_finite();
        if finite
            && self.dx_m > 0.0
            && self.nx >= 1
            && self.ny >= 1
            && self.nz >= 1
            && self.n_points() <= MAX_POINTS
        {
            Ok(())
        } else {
            Err(refuse(
                "grid-spec-invalid",
                format!("{self:?} ({} points)", self.n_points()),
                "finite origin, dx > 0, 1..=262144 points",
            ))
        }
    }

    /// Point i (x-fastest ordering, deterministic).
    #[must_use]
    pub fn point(&self, i: usize) -> [f64; 3] {
        let ix = i % self.nx;
        let iy = (i / self.nx) % self.ny;
        let iz = i / (self.nx * self.ny);
        [
            self.origin_m[0] + self.dx_m * ix as f64,
            self.origin_m[1] + self.dx_m * iy as f64,
            self.origin_m[2] + self.dx_m * iz as f64,
        ]
    }
}

/// A simple bound-circulation source: one horseshoe (span line + two
/// trailing legs downstream to a declared cutoff) with core radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundSystem {
    /// Circulation [m²/s].
    pub gamma_m2ps: f64,
    /// Left wingtip [m].
    pub tip_left_m: [f64; 3],
    /// Right wingtip [m].
    pub tip_right_m: [f64; 3],
    /// Trailing-leg length [m] (downstream −x).
    pub trail_m: f64,
    /// Singularity core radius [m].
    pub core_m: f64,
}

/// The v1 field-source state (native side; the wasm lease binds a
/// FieldSourceSnapshotV1 onto this in E7.1-ii).
#[derive(Clone, Debug)]
pub struct FieldSourceStateV1 {
    /// Source tick.
    pub tick: u64,
    /// Digest of the physical state this field derives from.
    pub source_state_digest: String,
    /// Atmosphere (mean + turbulence + air constants).
    pub atmosphere: Atmosphere,
    /// Bound system, if the model published one.
    pub bound: Option<BoundSystem>,
    /// Image plane active (FlatnessCertificate flat)?
    pub images_active: bool,
}

impl FieldSourceStateV1 {
    /// Components this STATE can physically produce.
    #[must_use]
    pub fn supported_components(&self) -> u32 {
        let mut m = C_MEAN_ATMO | C_TURB_ATMO;
        if self.bound.is_some() {
            m |= C_BOUND_CIRCULATION;
            if self.images_active {
                m |= C_GROUND_IMAGES;
            }
        }
        m
    }

    /// The state-bound field-source snapshot id.
    #[must_use]
    pub fn field_source_snapshot_id(&self) -> String {
        let mut b = self.source_state_digest.as_bytes().to_vec();
        b.extend_from_slice(&self.tick.to_le_bytes());
        b.extend_from_slice(&self.supported_components().to_le_bytes());
        hash_domain("org.frankensim.wf.field-source-snapshot.v1", &b).to_hex()
    }
}

/// Request-level metadata (§5.5 meta block).
#[derive(Clone, Debug, PartialEq)]
pub struct FieldMeta {
    /// State-bound snapshot id.
    pub field_source_snapshot_id: String,
    /// Source tick.
    pub source_tick: u64,
    /// Source state digest.
    pub source_state_digest: String,
    /// Components the state supports.
    pub source_modes: u32,
    /// Force-coupled subset of the DEFINED components.
    pub force_coupled_components: u32,
    /// Visualization-only subset.
    pub visualization_only_components: u32,
    /// Names of every component NOT in this sample's sum.
    pub omitted_components: Vec<&'static str>,
    /// Singularity core radius [m] (0 when no singular component).
    pub core_radius_m: f64,
    /// Export precision (bits of the payload floats).
    pub export_precision_bits: u32,
}

/// The §5.5 return set (arrays are grid-ordered, length n_points).
#[derive(Clone, Debug)]
pub struct FieldSampleSet {
    /// Velocity [m/s].
    pub u: Vec<[f64; 3]>,
    /// Velocity gradient du_i/dx_j.
    pub grad_u: Vec<[[f64; 3]; 3]>,
    /// Vorticity ∇×u [1/s].
    pub omega: Vec<[f64; 3]>,
    /// Divergence from the analytic gradient [1/s].
    pub div_analytic: Vec<f64>,
    /// Divergence by central differences of u [1/s] (the dual).
    pub div_finite_difference: Vec<f64>,
    /// ‖S‖ Frobenius [1/s].
    pub strain_magnitude: Vec<f64>,
    /// Q = ½(‖Ω‖² − ‖S‖²) [1/s²].
    pub q_criterion: Vec<f64>,
    /// λ₂ of S² + Ω² [1/s²].
    pub lambda2: Vec<f64>,
    /// ∇|u| [1/s] (zero-masked below the speed floor).
    pub kinetic_speed_gradient: Vec<[f64; 3]>,
    /// Sample validity (domain + finite).
    pub validity_mask: Vec<bool>,
    /// Inside a singularity core (normalized quantities masked).
    pub singularity_core_mask: Vec<bool>,
    /// Inside solid exclusion (below the plane).
    pub solid_exclusion_mask: Vec<bool>,
    /// The mask this sample summed.
    pub component_mask: u32,
    /// Per-component provenance strings, selected order.
    pub provenance: Vec<&'static str>,
    /// Meta block.
    pub meta: FieldMeta,
}

fn seg_velocity_cored(p: [f64; 3], a: [f64; 3], b: [f64; 3], core: f64) -> [f64; 3] {
    // Biot–Savart finite segment with a smooth core (Vatistas n=2
    // class): factor = c²/(c⁴ + core⁴)^{1/2} — deterministic sqrt only.
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let c = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    let l1 = det::sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    let l2 = det::sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]);
    if l1 < 1e-12 || l2 < 1e-12 {
        return [0.0; 3];
    }
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let dot = ab[0] * (r1[0] / l1 - r2[0] / l2)
        + ab[1] * (r1[1] / l1 - r2[1] / l2)
        + ab[2] * (r1[2] / l1 - r2[2] / l2);
    let denom = det::sqrt(c2 * c2 + core * core * core * core);
    if denom < 1e-24 {
        return [0.0; 3];
    }
    let k = dot * c2 / denom / (4.0 * core::f64::consts::PI * c2.max(1e-24));
    [k * c[0], k * c[1], k * c[2]]
}

fn horseshoe_velocity(p: [f64; 3], b: &BoundSystem, mirror: bool) -> [f64; 3] {
    let (sgn, refl) = if mirror { (-1.0, -1.0) } else { (1.0, 1.0) };
    let m = |q: [f64; 3]| [q[0], q[1], refl * q[2]];
    let l = m(b.tip_left_m);
    let r = m(b.tip_right_m);
    let lt = [l[0] - b.trail_m, l[1], l[2]];
    let rt = [r[0] - b.trail_m, r[1], r[2]];
    let g = sgn * b.gamma_m2ps;
    let mut v = [0.0f64; 3];
    for (a, bb, gg) in [(lt, l, g), (l, r, g), (r, rt, g)] {
        let s = seg_velocity_cored(p, a, bb, b.core_m);
        v[0] += gg * s[0];
        v[1] += gg * s[1];
        v[2] += gg * s[2];
    }
    v
}

/// Distance from p to the bound span line (core-mask predicate).
fn dist_to_span(p: [f64; 3], b: &BoundSystem) -> f64 {
    let a = b.tip_left_m;
    let c = b.tip_right_m;
    let ab = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let l2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let t = if l2 > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let q = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    det::sqrt(
        (p[0] - q[0]) * (p[0] - q[0])
            + (p[1] - q[1]) * (p[1] - q[1])
            + (p[2] - q[2]) * (p[2] - q[2]),
    )
}

/// Velocity of the SELECTED sum at a point (used for both the direct
/// samples and the FD probes — one code path, no drift).
fn u_at(state: &FieldSourceStateV1, mask: u32, p: [f64; 3], tick: u64) -> [f64; 3] {
    let mut u = [0.0f64; 3];
    if p[2] >= 0.0 {
        if mask & C_TURB_ATMO != 0 {
            let ts = state.atmosphere.turbulence.sample(p[0], p[1], p[2], tick);
            u[0] += ts.u[0];
            u[1] += ts.u[1];
            u[2] += ts.u[2];
        }
        if mask & C_MEAN_ATMO != 0 {
            u[0] += state.atmosphere.mean.speed(p[2]);
        }
    }
    if let Some(b) = &state.bound {
        if mask & C_BOUND_CIRCULATION != 0 {
            let v = horseshoe_velocity(p, b, false);
            u[0] += v[0];
            u[1] += v[1];
            u[2] += v[2];
        }
        if mask & C_GROUND_IMAGES != 0 && state.images_active {
            let v = horseshoe_velocity(p, b, true);
            u[0] += v[0];
            u[1] += v[1];
            u[2] += v[2];
        }
    }
    u
}

/// Analytic gradient of the ANALYTIC components (atmosphere); the
/// singular components contribute through the FD dual only, which the
/// grad field completes by central differences for those bits.
fn analytic_grad(state: &FieldSourceStateV1, mask: u32, p: [f64; 3], tick: u64) -> [[f64; 3]; 3] {
    let mut g = [[0.0f64; 3]; 3];
    if p[2] >= 0.0 {
        if mask & C_TURB_ATMO != 0 {
            let ts = state.atmosphere.turbulence.sample(p[0], p[1], p[2], tick);
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += ts.grad[i][j];
                }
            }
        }
        if mask & C_MEAN_ATMO != 0 {
            g[0][2] += state.atmosphere.mean.dspeed_dh(p[2]);
        }
    }
    g
}

/// Sample the field (plan §5.5).
///
/// # Errors
/// Grid refusals; `component-mask-empty`; `component-mask-unknown`
/// (undefined bits); `component-unsupported` only when the caller
/// requests a component the STATE cannot produce (honest refusal, not
/// a silent zero).
pub fn sample_field(
    state: &FieldSourceStateV1,
    grid: &GridSpec,
    component_mask: u32,
) -> Result<FieldSampleSet, Refusal> {
    grid.admit()?;
    if component_mask == 0 {
        return Err(refuse(
            "component-mask-empty",
            "no components selected".into(),
            "select at least one component bit",
        ));
    }
    if component_mask & !C_ALL_DEFINED != 0 {
        return Err(refuse(
            "component-mask-unknown",
            format!("undefined bits in {component_mask:#x}"),
            "use the C_* component constants",
        ));
    }
    let supported = state.supported_components();
    if component_mask & !supported != 0 {
        let missing: Vec<&str> = (0..8)
            .filter(|b| component_mask & !supported & (1 << b) != 0)
            .map(|b| COMPONENT_NAMES[b as usize])
            .collect();
        return Err(refuse(
            "component-unsupported",
            format!("state v1 cannot produce: {missing:?}"),
            "drop the unsupported bits; their omission is named in meta",
        ));
    }
    let n = grid.n_points();
    let tick = state.tick;
    let has_singular =
        component_mask & (C_BOUND_CIRCULATION | C_GROUND_IMAGES) != 0 && state.bound.is_some();
    let core = state.bound.as_ref().map_or(0.0, |b| b.core_m);
    let mut out = FieldSampleSet {
        u: Vec::with_capacity(n),
        grad_u: Vec::with_capacity(n),
        omega: Vec::with_capacity(n),
        div_analytic: Vec::with_capacity(n),
        div_finite_difference: Vec::with_capacity(n),
        strain_magnitude: Vec::with_capacity(n),
        q_criterion: Vec::with_capacity(n),
        lambda2: Vec::with_capacity(n),
        kinetic_speed_gradient: Vec::with_capacity(n),
        validity_mask: Vec::with_capacity(n),
        singularity_core_mask: Vec::with_capacity(n),
        solid_exclusion_mask: Vec::with_capacity(n),
        component_mask,
        provenance: (0..8)
            .filter(|b| component_mask & (1 << b) != 0)
            .map(|b| COMPONENT_NAMES[b as usize])
            .collect(),
        meta: FieldMeta {
            field_source_snapshot_id: state.field_source_snapshot_id(),
            source_tick: tick,
            source_state_digest: state.source_state_digest.clone(),
            source_modes: supported,
            force_coupled_components: FORCE_COUPLED,
            visualization_only_components: C_VIS_ONLY,
            omitted_components: (0..8)
                .filter(|b| component_mask & (1 << b) == 0)
                .map(|b| COMPONENT_NAMES[b as usize])
                .collect(),
            core_radius_m: if has_singular { core } else { 0.0 },
            export_precision_bits: 64,
        },
    };
    let singular_mask_bits = component_mask & (C_BOUND_CIRCULATION | C_GROUND_IMAGES);
    for i in 0..n {
        let p = grid.point(i);
        let solid = p[2] < 0.0;
        let u = u_at(state, component_mask, p, tick);
        // Gradient: analytic for the atmosphere + central differences
        // for the singular components (one shared u_at path).
        let mut g = analytic_grad(state, component_mask, p, tick);
        if singular_mask_bits != 0 {
            for j in 0..3 {
                let mut pp = p;
                let mut pm = p;
                pp[j] += FD_H;
                pm[j] -= FD_H;
                let up = u_at(state, singular_mask_bits, pp, tick);
                let um = u_at(state, singular_mask_bits, pm, tick);
                for gi in 0..3 {
                    g[gi][j] += (up[gi] - um[gi]) / (2.0 * FD_H);
                }
            }
        }
        // FD divergence dual over the FULL selected sum.
        let mut div_fd = 0.0;
        for j in 0..3 {
            let mut pp = p;
            let mut pm = p;
            pp[j] += FD_H;
            pm[j] -= FD_H;
            div_fd += (u_at(state, component_mask, pp, tick)[j]
                - u_at(state, component_mask, pm, tick)[j])
                / (2.0 * FD_H);
        }
        let div_an = g[0][0] + g[1][1] + g[2][2];
        // Tensor invariants.
        let mut s = [[0.0f64; 3]; 3];
        let mut w = [[0.0f64; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] = 0.5 * (g[a][b] + g[b][a]);
                w[a][b] = 0.5 * (g[a][b] - g[b][a]);
            }
        }
        let frob = |m: &[[f64; 3]; 3]| -> f64 { m.iter().flatten().map(|v| v * v).sum::<f64>() };
        let s2 = frob(&s);
        let w2 = frob(&w);
        let q = 0.5 * (w2 - s2);
        // λ₂ of A = S² + Ω² (symmetric): closed-form eigenvalues.
        let mut a_m = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    a_m[r][c] += s[r][k] * s[k][c] + w[r][k] * w[k][c];
                }
            }
        }
        let l2v = sym3_middle_eigenvalue(&a_m);
        // Kinetic speed gradient: ∇|u| = (∇u)ᵀ û.
        let speed = det::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        let ksg = if speed > 1e-9 {
            [
                (g[0][0] * u[0] + g[1][0] * u[1] + g[2][0] * u[2]) / speed,
                (g[0][1] * u[0] + g[1][1] * u[1] + g[2][1] * u[2]) / speed,
                (g[0][2] * u[0] + g[1][2] * u[1] + g[2][2] * u[2]) / speed,
            ]
        } else {
            [0.0; 3]
        };
        let in_core = has_singular
            && state
                .bound
                .as_ref()
                .is_some_and(|b| dist_to_span(p, b) < 4.0 * b.core_m);
        let finite_ok =
            u.iter().all(|v| v.is_finite()) && g.iter().flatten().all(|v| v.is_finite());
        out.u.push(u);
        out.grad_u.push(g);
        out.omega
            .push([g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]]);
        out.div_analytic.push(div_an);
        out.div_finite_difference.push(div_fd);
        out.strain_magnitude.push(det::sqrt(s2));
        out.q_criterion.push(q);
        out.lambda2.push(l2v);
        out.kinetic_speed_gradient.push(ksg);
        out.validity_mask.push(finite_ok && !solid);
        out.singularity_core_mask.push(in_core);
        out.solid_exclusion_mask.push(solid);
    }
    Ok(out)
}

impl FieldSampleSet {
    /// Normalized divergence ε = |∇·u| / ‖∇u‖ at point i, masked
    /// (None) under the gradient floor or inside a singularity core —
    /// the divergence overlay shows BOTH this and the absolute value.
    #[must_use]
    pub fn normalized_divergence(&self, i: usize) -> Option<f64> {
        let g = &self.grad_u[i];
        let norm = det::sqrt(g.iter().flatten().map(|v| v * v).sum::<f64>());
        if norm < GRAD_NORM_FLOOR || self.singularity_core_mask[i] || !self.validity_mask[i] {
            None
        } else {
            Some(self.div_analytic[i].abs() / norm)
        }
    }

    /// Content digest over the payload (goldens).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = Vec::new();
        for (u, g) in self.u.iter().zip(self.grad_u.iter()) {
            for v in u.iter().chain(g.iter().flatten()) {
                b.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
        for v in self
            .div_analytic
            .iter()
            .chain(self.div_finite_difference.iter())
            .chain(self.q_criterion.iter())
            .chain(self.lambda2.iter())
        {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.wf.fieldsvc.v1", &b).to_hex()
    }
}

/// The total-flow CLAIM gate (plan §5.5): the sum may be labeled
/// "total flow" only when every force-coupled component the active
/// model supports is present.
///
/// # Errors
/// `forbidden-claim-total-flow` (names the absent components).
pub fn claim_total_flow(sample: &FieldSampleSet) -> Result<(), Refusal> {
    let owed = sample.meta.source_modes & FORCE_COUPLED;
    let missing = owed & !sample.component_mask;
    if missing == 0 {
        Ok(())
    } else {
        let names: Vec<&str> = (0..8)
            .filter(|b| missing & (1 << b) != 0)
            .map(|b| COMPONENT_NAMES[b as usize])
            .collect();
        Err(refuse(
            "forbidden-claim-total-flow",
            format!("force-coupled components absent from the sum: {names:?}"),
            "include them, or present the overlay under its component list",
        ))
    }
}

/// Middle eigenvalue of a symmetric 3×3 (deterministic closed form,
/// det:: trig only).
fn sym3_middle_eigenvalue(a: &[[f64; 3]; 3]) -> f64 {
    let q = (a[0][0] + a[1][1] + a[2][2]) / 3.0;
    let mut b = *a;
    for i in 0..3 {
        b[i][i] -= q;
    }
    let p2 = b.iter().flatten().map(|v| v * v).sum::<f64>() / 6.0;
    let p = det::sqrt(p2);
    if p < 1e-30 {
        return q;
    }
    // det(B)/2p³ clamped into [-1, 1].
    let detb = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    let r = (detb / (2.0 * p * p2)).clamp(-1.0, 1.0);
    let phi = det::atan2(det::sqrt((1.0 - r * r).max(0.0)), r) / 3.0;
    let e1 = q + 2.0 * p * det::cos(phi);
    let e3 = q + 2.0 * p * det::cos(phi + 2.0 * core::f64::consts::PI / 3.0);
    // λ₂ = trace − λ₁ − λ₃.
    3.0 * q - e1 - e3
}
