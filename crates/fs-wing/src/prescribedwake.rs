//! PrescribedWakeReference3d (bead wf-root-guzez.5.7, E4.3b1). Plan
//! Round-3/A1 lane: the REFERENCE wake operator — a prescribed (never
//! free) trailing-vortex system, IMAGE-AWARE FROM BIRTH: every wake row
//! carries its complete FlatPlaneVortexImageExact mirror below the
//! certified ground plane, so ground effect is in the operator from the
//! first segment, not patched on. A FROZEN operating/linearization grid
//! (h/b × pitch × roll × canard × warp × convection) is identity-
//! registered; the A1 FOM/ROM lane and the E4.9b referee consume the
//! operator AT these points only.
//!
//! Discretization: each span station sheds a trailing line broken into
//! `rows` straight segments along the prescribed convection path
//! (downstream = −x, pitched with the operating point), closed by a
//! far-cap semi-infinite leg. The operator maps station circulations to
//! induced velocity at probe points via the shared Biot–Savart segment
//! kernel + the exact mirror (z → 2·z_g − z with the sign flip carried
//! by traversal order).

use crate::images::CertifiedGround;
use crate::{Refusal, refuse};
use fs_math::det;

/// Wake-row cap (absurd-input guard; 120 rows ≈ 1 s of shed wake).
pub const MAX_WAKE_ROWS: usize = 512;

/// Far cap distance appended after the last row [m].
pub const FAR_CAP_M: f64 = 1.0e4;

/// Segment Biot–Savart core guard.
const CORE: f64 = 1.0e-8;

/// One frozen operating point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeOperatingPoint {
    /// Height over span.
    pub h_over_b: f64,
    /// Pitch attitude [rad].
    pub pitch_rad: f64,
    /// Roll attitude [rad].
    pub roll_rad: f64,
    /// Canard deflection [rad].
    pub canard_rad: f64,
    /// Warp command [rad].
    pub warp_rad: f64,
    /// Wake convection ratio (wake speed / freestream).
    pub convection: f64,
}

/// The frozen grid (identity-registered).
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenGrid {
    /// Points in canonical (declared) order.
    pub points: Vec<WakeOperatingPoint>,
    /// Registered identity digest.
    pub digest: String,
}

/// The declared v1 axes (the registered spans; the grid is their full
/// product in canonical nesting order).
#[must_use]
pub fn frozen_grid_v1() -> FrozenGrid {
    let h_over_b = [0.1, 0.2, 0.3, 10.0];
    let pitch = [-0.05, 0.05, 0.15];
    let roll = [0.0, 0.10];
    let canard = [-0.2, 0.0, 0.2];
    let warp = [0.0, 0.08];
    let convection = [0.85, 1.0];
    let mut points = Vec::new();
    for &h in &h_over_b {
        for &p in &pitch {
            for &r in &roll {
                for &c in &canard {
                    for &w in &warp {
                        for &cv in &convection {
                            points.push(WakeOperatingPoint {
                                h_over_b: h,
                                pitch_rad: p,
                                roll_rad: r,
                                canard_rad: c,
                                warp_rad: w,
                                convection: cv,
                            });
                        }
                    }
                }
            }
        }
    }
    let digest = grid_digest(&points);
    FrozenGrid { points, digest }
}

fn grid_digest(points: &[WakeOperatingPoint]) -> String {
    let mut p = Vec::new();
    for pt in points {
        for v in [
            pt.h_over_b,
            pt.pitch_rad,
            pt.roll_rad,
            pt.canard_rad,
            pt.warp_rad,
            pt.convection,
        ] {
            p.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    fs_blake3::hash_domain("org.frankensim.fs-wing.wake-grid.v1", &p).to_hex()
}

/// The assembled prescribed-wake operator at one operating point:
/// `w_normal[probe][station]` maps station circulation to the induced
/// z-velocity at each probe.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedWakeOperator {
    /// Probe count.
    pub n_probes: usize,
    /// Station count.
    pub n_stations: usize,
    /// Row-major influence matrix [m/s per unit Γ], z-component.
    pub w_normal: Vec<f64>,
    /// The operating point this operator is FOR.
    pub point: WakeOperatingPoint,
    /// Whether images were included (always true through the public
    /// constructor with ground; the image-free twin exists for the
    /// battery's liveness checks only).
    pub image_aware: bool,
}

/// Biot–Savart velocity of a finite straight segment a→b at p.
fn segment_velocity(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let n1 = det::sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    let n2 = det::sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]);
    let c = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    if c2 < CORE || n1 < CORE || n2 < CORE {
        return [0.0; 3];
    }
    let r0 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let k = ((r0[0] * r1[0] + r0[1] * r1[1] + r0[2] * r1[2]) / n1
        - (r0[0] * r2[0] + r0[1] * r2[1] + r0[2] * r2[2]) / n2)
        / (4.0 * core::f64::consts::PI * c2);
    [c[0] * k, c[1] * k, c[2] * k]
}

/// Trailing-line polyline for one station under an operating point:
/// starts at the shed point and convects downstream (−x), pitched with
/// the attitude; `span_b` scales the height.
fn wake_polyline(
    shed: [f64; 3],
    point: &WakeOperatingPoint,
    rows: usize,
    row_dx_m: f64,
) -> Vec<[f64; 3]> {
    let mut line = Vec::with_capacity(rows + 2);
    let mut x = shed;
    line.push(x);
    // The prescribed path: straight at the pitched convection direction.
    let dir = [
        -det::cos(point.pitch_rad),
        0.0,
        -det::sin(point.pitch_rad), // pitched wake plane
    ];
    for k in 0..rows {
        let step = row_dx_m * point.convection;
        x = [
            x[0] + dir[0] * step,
            x[1] + dir[1] * step,
            x[2] + dir[2] * step,
        ];
        line.push(x);
        let _ = k;
    }
    // Far cap.
    line.push([
        x[0] + dir[0] * FAR_CAP_M,
        x[1] + dir[1] * FAR_CAP_M,
        x[2] + dir[2] * FAR_CAP_M,
    ]);
    line
}

/// Assemble the operator: probes × stations, each station's trailing
/// line + (image_aware) its exact mirror below the ground plane.
///
/// # Errors
/// `wake-rows-invalid` (0 or above [`MAX_WAKE_ROWS`] — cap AND cap+1);
/// `wake-grid-empty`; ground-certificate refusals pass through;
/// `wake-geometry-invalid` (aircraft below the ground plane).
#[allow(clippy::too_many_arguments)]
pub fn assemble_operator(
    shed_points: &[[f64; 3]],
    probes: &[[f64; 3]],
    ground: &CertifiedGround,
    point: &WakeOperatingPoint,
    rows: usize,
    row_dx_m: f64,
    span_b_m: f64,
) -> Result<PrescribedWakeOperator, Refusal> {
    ground.admit()?;
    if rows == 0 || rows > MAX_WAKE_ROWS {
        return Err(refuse(
            "wake-rows-invalid",
            format!("{rows} rows outside [1, {MAX_WAKE_ROWS}]"),
            "120 rows is the 1-second reference class",
        ));
    }
    if shed_points.is_empty() || probes.is_empty() {
        return Err(refuse(
            "wake-grid-empty",
            "no stations or probes".into(),
            "supply the surface trailing-edge stations",
        ));
    }
    let zg = ground.z_m;
    // The operating point's h/b sets the aircraft height ABOVE ground:
    // shed z = zg − h (FRD +z down: above ground = smaller z).
    let dz = -(point.h_over_b * span_b_m);
    for s in shed_points {
        if s[2] + dz >= zg {
            return Err(refuse(
                "wake-geometry-invalid",
                format!(
                    "shed point below the ground plane at h/b {}",
                    point.h_over_b
                ),
                "the aircraft must sit above the certified plane",
            ));
        }
    }
    let n_probes = probes.len();
    let n_stations = shed_points.len();
    let mut w = vec![0.0f64; n_probes * n_stations];
    for (js, shed0) in shed_points.iter().enumerate() {
        let shed = [shed0[0], shed0[1], shed0[2] + dz];
        let line = wake_polyline(shed, point, rows, row_dx_m);
        // Mirror line: z → 2·zg − z, traversed in the SAME order — the
        // sign flip of the image circulation is carried by negating the
        // contribution (equivalent to reversing traversal).
        let mirror: Vec<[f64; 3]> = line.iter().map(|q| [q[0], q[1], 2.0 * zg - q[2]]).collect();
        for (ip, probe0) in probes.iter().enumerate() {
            let probe = [probe0[0], probe0[1], probe0[2] + dz];
            let mut wz = 0.0f64;
            for seg in line.windows(2) {
                wz += segment_velocity(probe, seg[0], seg[1])[2];
            }
            for seg in mirror.windows(2) {
                wz -= segment_velocity(probe, seg[0], seg[1])[2];
            }
            w[ip * n_stations + js] = wz;
        }
    }
    Ok(PrescribedWakeOperator {
        n_probes,
        n_stations,
        w_normal: w,
        point: *point,
        image_aware: true,
    })
}

/// The IMAGE-FREE twin (battery liveness only — the reference lane
/// never consumes it; the flag says so).
///
/// # Errors
/// Same validation as [`assemble_operator`].
#[allow(clippy::too_many_arguments)]
pub fn assemble_operator_no_images(
    shed_points: &[[f64; 3]],
    probes: &[[f64; 3]],
    ground: &CertifiedGround,
    point: &WakeOperatingPoint,
    rows: usize,
    row_dx_m: f64,
    span_b_m: f64,
) -> Result<PrescribedWakeOperator, Refusal> {
    let mut op = assemble_operator(shed_points, probes, ground, point, rows, row_dx_m, span_b_m)?;
    // Rebuild without mirrors.
    let zg = ground.z_m;
    let dz = -(point.h_over_b * span_b_m);
    let n_stations = shed_points.len();
    for (js, shed0) in shed_points.iter().enumerate() {
        let shed = [shed0[0], shed0[1], shed0[2] + dz];
        let line = wake_polyline(shed, point, rows, row_dx_m);
        for (ip, probe0) in probes.iter().enumerate() {
            let probe = [probe0[0], probe0[1], probe0[2] + dz];
            let mut wz = 0.0f64;
            for seg in line.windows(2) {
                wz += segment_velocity(probe, seg[0], seg[1])[2];
            }
            op.w_normal[ip * n_stations + js] = wz;
        }
    }
    op.image_aware = false;
    let _ = zg;
    Ok(op)
}

/// Full-velocity probe of one station's wake (real + image) at a point
/// in the SAME shifted frame the operator uses — the battery's
/// plane-symmetry oracle needs all three components.
#[must_use]
pub fn probe_velocity(
    shed_shifted: [f64; 3],
    probe_shifted: [f64; 3],
    ground_z: f64,
    point: &WakeOperatingPoint,
    rows: usize,
    row_dx_m: f64,
) -> [f64; 3] {
    let line = wake_polyline(shed_shifted, point, rows, row_dx_m);
    let mirror: Vec<[f64; 3]> = line
        .iter()
        .map(|q| [q[0], q[1], 2.0 * ground_z - q[2]])
        .collect();
    let mut v = [0.0f64; 3];
    for seg in line.windows(2) {
        let s = segment_velocity(probe_shifted, seg[0], seg[1]);
        v[0] += s[0];
        v[1] += s[1];
        v[2] += s[2];
    }
    for seg in mirror.windows(2) {
        let s = segment_velocity(probe_shifted, seg[0], seg[1]);
        v[0] -= s[0];
        v[1] -= s[1];
        v[2] -= s[2];
    }
    v
}
