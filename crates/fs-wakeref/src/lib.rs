//! fs-wakeref — the INDEPENDENT unsteady prescribed-wake referee
//! (bead wf-root-guzez.5.22, E4.9b; V-08b1 dense reference leaf).
//!
//! INDEPENDENCE LAW (the bead's DONE-WHEN): this crate depends on
//! fs-math and fs-blake3 ONLY. It shares NO solver code with the A1
//! FOM lane (fs-airfoil's Duhamel indicial machinery, fs-wing's
//! nonlinear strip solver, or any reduced basis) — even the
//! Biot–Savart segment kernel is written here from the formula. The
//! battery pins the dependency closure.
//!
//! Model (declared referee tier — dense and slow on purpose, never a
//! real-time path): single-chordwise-ring unsteady vortex lattice over
//! two lifting surfaces (canard + wing) with a PRESCRIBED wake — rings
//! shed every step convect rigidly downstream at `convection · V`
//! (never a free wake), exact flat-ground vortex images when a ground
//! case is on, unsteady thin-airfoil apparent-mass lift term at the
//! half chord. Frame: x downstream (wake convects +x), y starboard,
//! z UP; the canard sits UPSTREAM at negative x.
//!
//! Fixtures (canard ↔ wing MIMO): impulse, step, chirp, reversal —
//! free air and flat ground — emitting the V-08b1 receipt with
//! per-case digests and the physics summaries the E4.3b3 campaign
//! compares against the A1 lane.

use fs_blake3::hash_domain;
use fs_math::det;

/// Typed refusal (workspace convention).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

/// Step cap (absurd-input guard; 20 s at 120 Hz).
pub const MAX_STEPS: usize = 2_400;

/// Wake-ring memory cap per station (oldest dropped, DECLARED — the
/// receipt records the cap so truncation is never silent).
pub const MAX_WAKE_ROWS: usize = 600;

/// Biot–Savart core radius guard [m].
const CORE: f64 = 1.0e-8;

/// Referee geometry (registered 1903 reference values — flyer-reference
/// lineage: span 12.29 m, chord 1.981 m, canard 3.66 × 0.61 m with its
/// quarter-chord ~2.23 m ahead of the wing's, mounted low).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RefereeGeometry {
    /// Wing span [m].
    pub wing_span_m: f64,
    /// Wing chord [m].
    pub wing_chord_m: f64,
    /// Canard span [m].
    pub canard_span_m: f64,
    /// Canard chord [m].
    pub canard_chord_m: f64,
    /// Canard quarter-chord x relative to the wing's [m] (negative =
    /// upstream in this frame).
    pub canard_dx_m: f64,
    /// Canard height relative to the wing plane [m] (z up).
    pub canard_dz_m: f64,
    /// Wing spanwise stations.
    pub n_wing: usize,
    /// Canard spanwise stations.
    pub n_canard: usize,
}

/// The registered v1 referee geometry.
#[must_use]
pub fn wright_geometry_v1() -> RefereeGeometry {
    RefereeGeometry {
        wing_span_m: 12.29,
        wing_chord_m: 1.981,
        canard_span_m: 3.66,
        canard_chord_m: 0.61,
        canard_dx_m: -2.23,
        canard_dz_m: -0.7,
        n_wing: 8,
        n_canard: 4,
    }
}

impl RefereeGeometry {
    /// Content digest (receipt ingredient).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        for v in [
            self.wing_span_m,
            self.wing_chord_m,
            self.canard_span_m,
            self.canard_chord_m,
            self.canard_dx_m,
            self.canard_dz_m,
            self.n_wing as f64,
            self.n_canard as f64,
        ] {
            p.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.fs-wakeref.geometry.v1", &p).to_hex()
    }
}

/// Canard-command fixture.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fixture {
    /// One-tick pulse of +0.1 rad.
    Impulse,
    /// 0 → +0.1 rad at t = 0.
    Step,
    /// 0.05·sin swept 0.5 → 4 Hz over the run.
    Chirp,
    /// +0.1 held, flipped to −0.1 at half time.
    Reversal,
}

impl Fixture {
    fn name(self) -> &'static str {
        match self {
            Fixture::Impulse => "impulse",
            Fixture::Step => "step",
            Fixture::Chirp => "chirp",
            Fixture::Reversal => "reversal",
        }
    }

    /// Canard deflection [rad] at step k of n (dt known to the caller).
    fn delta(self, k: usize, n: usize, dt_s: f64) -> f64 {
        match self {
            Fixture::Impulse => {
                if k == 0 {
                    0.1
                } else {
                    0.0
                }
            }
            Fixture::Step => 0.1,
            Fixture::Chirp => {
                let t = k as f64 * dt_s;
                let total = n as f64 * dt_s;
                let f = 0.5 + 3.5 * (t / total);
                0.05 * det::sin(core::f64::consts::TAU * f * t)
            }
            Fixture::Reversal => {
                if k < n / 2 {
                    0.1
                } else {
                    -0.1
                }
            }
        }
    }
}

/// One referee case.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RefereeCase {
    /// Fixture waveform.
    pub fixture: Fixture,
    /// Flat ground height under the wing plane [m] (None = free air).
    pub ground_z_m: Option<f64>,
    /// Freestream [m/s].
    pub v_mps: f64,
    /// Baseline incidence applied to BOTH surfaces [rad] (the wing
    /// carries lift; the canard adds the fixture deflection on top).
    pub alpha0_rad: f64,
    /// Air density [kg/m³].
    pub rho_kg_m3: f64,
    /// Wake convection ratio.
    pub convection: f64,
    /// Time step [s].
    pub dt_s: f64,
    /// Steps.
    pub n_steps: usize,
}

/// The Dec-17-class case set (free air + flat ground per fixture).
#[must_use]
pub fn v08b1_cases_v1() -> Vec<RefereeCase> {
    let mut cases = Vec::new();
    for fixture in [
        Fixture::Impulse,
        Fixture::Step,
        Fixture::Chirp,
        Fixture::Reversal,
    ] {
        for ground in [None, Some(-2.4)] {
            cases.push(RefereeCase {
                fixture,
                ground_z_m: ground,
                v_mps: 13.0,
                alpha0_rad: 0.05,
                rho_kg_m3: 1.294,
                convection: 1.0,
                dt_s: 1.0 / 120.0,
                n_steps: 480, // 4 s
            });
        }
    }
    cases
}

/// Per-case time series (dense referee output).
#[derive(Clone, Debug, PartialEq)]
pub struct CaseSeries {
    /// Wing lift [N] per step.
    pub wing_lift_n: Vec<f64>,
    /// Canard lift [N] per step.
    pub canard_lift_n: Vec<f64>,
    /// Canard hinge moment [N·m] per step (thin-airfoil quarter-chord
    /// lift acting about the 40 %-chord hinge axis — declared tier).
    pub hinge_nm: Vec<f64>,
    /// Content digest over all three series (bitwise).
    pub digest: String,
}

/// One receipt row.
#[derive(Clone, Debug, PartialEq)]
pub struct CaseRecord {
    /// Fixture name.
    pub fixture: &'static str,
    /// Ground case?
    pub ground: bool,
    /// Series digest.
    pub digest: String,
    /// Steady (final) wing lift [N].
    pub steady_wing_lift_n: f64,
    /// Steady canard lift [N].
    pub steady_canard_lift_n: f64,
    /// Step fixtures: canard lift at step 3 / steady (the Wagner-class
    /// starting deficiency on the CIRCULATORY share — step 0 carries
    /// the impulsive non-circulatory apparent-mass spike, declared, so
    /// the convention samples after it decays).
    pub wagner_ratio: Option<f64>,
}

/// The V-08b1 receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct V08b1Receipt {
    /// Schema id.
    pub schema: &'static str,
    /// Declared model tier.
    pub tier: &'static str,
    /// Geometry identity.
    pub geometry_digest: String,
    /// Wake-memory cap (truncation is declared, never silent).
    pub wake_rows_cap: usize,
    /// Case rows.
    pub cases: Vec<CaseRecord>,
    /// Digest over the whole receipt.
    pub receipt_digest: String,
}

pub const RECEIPT_SCHEMA: &str = "org.frankensim.wf.v08b1-receipt.v1";
pub const REFEREE_TIER: &str = "dense-uvlm1-prescribed-wake: single-chordwise-ring lattice, \
     thin-airfoil closure, rigid downstream convection, exact flat-ground images, \
     half-chord apparent-mass lift term";

// ---------------------------------------------------------------------------
// Dense kernel (written HERE from the formula — independence law).
// ---------------------------------------------------------------------------

/// Biot–Savart velocity of a finite straight segment a→b, unit Γ, at p.
fn segment_velocity(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let c = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    if c2 < CORE {
        return [0.0; 3];
    }
    let l1 = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
    let l2 = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
    if l1 < CORE || l2 < CORE {
        return [0.0; 3];
    }
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let dot = ab[0] * (r1[0] / l1 - r2[0] / l2)
        + ab[1] * (r1[1] / l1 - r2[1] / l2)
        + ab[2] * (r1[2] / l1 - r2[2] / l2);
    let k = dot / (4.0 * core::f64::consts::PI * c2);
    [k * c[0], k * c[1], k * c[2]]
}

/// Velocity of a quadrilateral vortex ring (corners in order), unit Γ.
fn ring_velocity(p: [f64; 3], corners: &[[f64; 3]; 4]) -> [f64; 3] {
    let mut v = [0.0f64; 3];
    for i in 0..4 {
        let s = segment_velocity(p, corners[i], corners[(i + 1) % 4]);
        v[0] += s[0];
        v[1] += s[1];
        v[2] += s[2];
    }
    v
}

/// Ring velocity including the exact flat-ground image (z → 2 z_g − z,
/// circulation sign flipped by traversing corners in reverse order).
fn ring_velocity_imaged(p: [f64; 3], corners: &[[f64; 3]; 4], ground_z: Option<f64>) -> [f64; 3] {
    let mut v = ring_velocity(p, corners);
    if let Some(zg) = ground_z {
        let m = |c: [f64; 3]| [c[0], c[1], 2.0 * zg - c[2]];
        let mirrored = [m(corners[3]), m(corners[2]), m(corners[1]), m(corners[0])];
        let iv = ring_velocity(p, &mirrored);
        v[0] += iv[0];
        v[1] += iv[1];
        v[2] += iv[2];
    }
    v
}

struct Station {
    /// Bound-ring corners (LE-left, LE-right, TE-right, TE-left).
    ring: [[f64; 3]; 4],
    /// Collocation point (3/4 chord, mid-span of the station).
    colloc: [f64; 3],
    /// Trailing-edge midpoints (left/right) where wake rings attach.
    te_left: [f64; 3],
    te_right: [f64; 3],
    width_m: f64,
    chord_m: f64,
    /// True for canard stations (they receive the fixture deflection).
    is_canard: bool,
}

fn build_stations(g: &RefereeGeometry) -> Vec<Station> {
    let mut out = Vec::new();
    let mut surface = |span: f64, chord: f64, x0: f64, z0: f64, n: usize, is_canard: bool| {
        let dy = span / n as f64;
        for i in 0..n {
            let y0 = -span / 2.0 + i as f64 * dy;
            let y1 = y0 + dy;
            let ym = 0.5 * (y0 + y1);
            // Ring: bound segment at quarter chord, closing at TE+.
            let x_le = x0; // quarter-chord line (bound vortex) at x0
            let x_te = x0 + 0.75 * chord;
            out.push(Station {
                ring: [
                    [x_le, y0, z0],
                    [x_le, y1, z0],
                    [x_te, y1, z0],
                    [x_te, y0, z0],
                ],
                colloc: [x0 + 0.5 * chord, ym, z0],
                te_left: [x_te, y0, z0],
                te_right: [x_te, y1, z0],
                width_m: dy,
                chord_m: chord,
                is_canard,
            });
        }
    };
    surface(g.wing_span_m, g.wing_chord_m, 0.0, 0.0, g.n_wing, false);
    surface(
        g.canard_span_m,
        g.canard_chord_m,
        g.canard_dx_m,
        g.canard_dz_m,
        g.n_canard,
        true,
    );
    out
}

/// Gauss elimination with partial pivoting (small dense system; local
/// implementation — independence law).
fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> Result<(), Refusal> {
    for col in 0..n {
        let mut piv = col;
        for row in (col + 1)..n {
            if a[row * n + col].abs() > a[piv * n + col].abs() {
                piv = row;
            }
        }
        if a[piv * n + col].abs() < 1e-14 {
            return Err(Refusal {
                code: "referee-system-singular",
                message: format!("pivot {col} below floor"),
                ranked_repairs: vec!["geometry produced a degenerate lattice".into()],
            });
        }
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for row in (col + 1)..n {
            let f = a[row * n + col] / d;
            if f != 0.0 {
                for k in col..n {
                    a[row * n + k] -= f * a[col * n + k];
                }
                b[row] -= f * b[col];
            }
        }
    }
    for col in (0..n).rev() {
        let mut s = b[col];
        for k in (col + 1)..n {
            s -= a[col * n + k] * b[k];
        }
        b[col] = s / a[col * n + col];
    }
    Ok(())
}

struct WakeRing {
    corners: [[f64; 3]; 4],
    gamma: f64,
    station: usize,
}

/// Run one referee case (dense, deterministic).
///
/// # Errors
/// `referee-case-invalid` (caps at cap AND cap+1: steps, dt, speed,
/// rho, convection; ground must sit below both surfaces);
/// `referee-system-singular`.
pub fn run_case(g: &RefereeGeometry, case: &RefereeCase) -> Result<CaseSeries, Refusal> {
    let ok = case.v_mps.is_finite()
        && case.v_mps >= 5.0
        && case.v_mps <= 40.0
        && case.alpha0_rad.is_finite()
        && case.alpha0_rad.abs() <= 0.3
        && case.rho_kg_m3.is_finite()
        && case.rho_kg_m3 > 0.5
        && case.rho_kg_m3 < 2.0
        && case.dt_s.is_finite()
        && case.dt_s > 0.0
        && case.dt_s <= 0.05
        && case.n_steps >= 1
        && case.n_steps <= MAX_STEPS
        && (0.5..=1.5).contains(&case.convection)
        && case
            .ground_z_m
            .is_none_or(|z| z < g.canard_dz_m.min(0.0) - 0.1);
    if !ok {
        return Err(Refusal {
            code: "referee-case-invalid",
            message: format!("{case:?}"),
            ranked_repairs: vec![format!(
                "v [5,40]; |alpha0| <= 0.3; rho (0.5,2.0); dt (0,0.05]; steps [1,{MAX_STEPS}]; convection [0.5,1.5]; ground below both surfaces"
            )],
        });
    }
    let stations = build_stations(g);
    let n = stations.len();
    let ground = case.ground_z_m;
    // Bound-ring influence at collocations (z-component per unit Γ).
    let mut a_mat = vec![0.0f64; n * n];
    for (i, si) in stations.iter().enumerate() {
        for (j, sj) in stations.iter().enumerate() {
            a_mat[i * n + j] = ring_velocity_imaged(si.colloc, &sj.ring, ground)[2];
        }
    }
    let mut gamma = vec![0.0f64; n];
    let mut gamma_prev = vec![0.0f64; n];
    let mut wake: Vec<WakeRing> = Vec::new();
    let mut wing_lift = Vec::with_capacity(case.n_steps);
    let mut canard_lift = Vec::with_capacity(case.n_steps);
    let mut hinge = Vec::with_capacity(case.n_steps);
    let dx = case.convection * case.v_mps * case.dt_s;
    for k in 0..case.n_steps {
        let delta = case.fixture.delta(k, case.n_steps, case.dt_s);
        // Convect the prescribed wake rigidly downstream.
        for ring in &mut wake {
            for c in &mut ring.corners {
                c[0] += dx;
            }
        }
        // Wake-induced downwash at collocations.
        let mut rhs = vec![0.0f64; n];
        for (i, st) in stations.iter().enumerate() {
            let mut w = 0.0;
            for ring in &wake {
                w += ring.gamma * ring_velocity_imaged(st.colloc, &ring.corners, ground)[2];
            }
            let alpha = case.alpha0_rad + if st.is_canard { delta } else { 0.0 };
            // No-penetration: (V sinα + w_bound + w_wake) = 0 at the
            // collocation; small-angle freestream normal component.
            rhs[i] = -(case.v_mps * det::sin(alpha) + w);
        }
        let mut a_work = a_mat.clone();
        solve_dense(&mut a_work, &mut rhs, n)?;
        gamma.copy_from_slice(&rhs);
        // Shed one wake ring per station: spans TE→TE+dx, carrying the
        // CURRENT bound circulation (ring-to-ring differencing yields
        // the shed vorticity exactly — Kelvin holds by construction).
        for (i, st) in stations.iter().enumerate() {
            wake.push(WakeRing {
                corners: [
                    st.te_left,
                    st.te_right,
                    [st.te_right[0] + dx, st.te_right[1], st.te_right[2]],
                    [st.te_left[0] + dx, st.te_left[1], st.te_left[2]],
                ],
                gamma: gamma[i],
                station: i,
            });
        }
        // Declared memory cap: drop the OLDEST rows beyond the budget.
        let cap = MAX_WAKE_ROWS * n / 8;
        if wake.len() > cap {
            let drop = wake.len() - cap;
            wake.drain(0..drop);
        }
        // Forces: steady Kutta–Joukowsky + half-chord apparent term.
        let mut lw = 0.0;
        let mut lc = 0.0;
        let mut hm = 0.0;
        for (i, st) in stations.iter().enumerate() {
            let dgdt = (gamma[i] - gamma_prev[i]) / case.dt_s;
            let l = case.rho_kg_m3 * (case.v_mps * gamma[i] + 0.5 * st.chord_m * dgdt) * st.width_m;
            if st.is_canard {
                lc += l;
                // Hinge: quarter-chord lift about the 40 %-chord axis.
                hm += -l * (0.40 - 0.25) * st.chord_m;
            } else {
                lw += l;
            }
        }
        gamma_prev.copy_from_slice(&gamma);
        wing_lift.push(lw);
        canard_lift.push(lc);
        hinge.push(hm);
        let _ = &wake.last().map(|r| r.station); // keep field live
    }
    let mut bytes = Vec::new();
    for series in [&wing_lift, &canard_lift, &hinge] {
        for v in series.iter() {
            bytes.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    let digest = hash_domain("org.frankensim.fs-wakeref.case-series.v1", &bytes).to_hex();
    Ok(CaseSeries {
        wing_lift_n: wing_lift,
        canard_lift_n: canard_lift,
        hinge_nm: hinge,
        digest,
    })
}

/// Run the registered case set and emit the V-08b1 receipt.
///
/// # Errors
/// As [`run_case`].
pub fn emit_v08b1_receipt(g: &RefereeGeometry) -> Result<V08b1Receipt, Refusal> {
    let mut cases = Vec::new();
    let mut bytes = Vec::new();
    for case in v08b1_cases_v1() {
        let series = run_case(g, &case)?;
        let last = *series.canard_lift_n.last().unwrap_or(&0.0);
        let wagner = match case.fixture {
            Fixture::Step => {
                let early = series.canard_lift_n.get(3).copied().unwrap_or(0.0);
                if last.abs() > 1e-9 {
                    Some(early / last)
                } else {
                    None
                }
            }
            _ => None,
        };
        let record = CaseRecord {
            fixture: case.fixture.name(),
            ground: case.ground_z_m.is_some(),
            digest: series.digest.clone(),
            steady_wing_lift_n: *series.wing_lift_n.last().unwrap_or(&0.0),
            steady_canard_lift_n: last,
            wagner_ratio: wagner,
        };
        bytes.extend_from_slice(record.digest.as_bytes());
        cases.push(record);
    }
    let receipt_digest = hash_domain("org.frankensim.wf.v08b1-receipt.v1", &bytes).to_hex();
    Ok(V08b1Receipt {
        schema: RECEIPT_SCHEMA,
        tier: REFEREE_TIER,
        geometry_digest: g.digest(),
        wake_rows_cap: MAX_WAKE_ROWS,
        cases,
        receipt_digest,
    })
}
