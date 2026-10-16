//! Swept critical-feature event localization (bead wf-root-guzez.4.9,
//! E3.4b, plan §5.1.5 Round-3 S-17 + Round-4 Q3). The 240 Hz penalty
//! update is PRECEDED by deterministic swept point/segment/CAPSULE vs
//! heightfield event localization for skids, canard frame, wingtips,
//! and PHASE-RESOLVED propeller blade capsules — at ~330–350 rpm the
//! blade phase changes materially within one step, so a disk proxy
//! false-triggers where neither blade is.
//!
//! `BladeCollisionProxyV1`: 16 capsules/blade baseline (two spanwise
//! rails at 25 % and 75 % local chord × eight radial intervals),
//! deterministic refinement to ≤ 24 when the cover certificate fails,
//! radii GENERATED from the provenance-bound blade-solid samples plus
//! declared margins (hand-entered radii have no API), hub void
//! excluded (the hub/drivetrain is a separate proxy), radial segments
//! ≤ 0.125 R. > 24 needed → the blade-strike CLAIM refuses and only
//! the conservative disk-clearance WARNING remains.
//!
//! A blade or hub strike emits `TerminalEvent::DamageModelUnavailable`
//! semantics: the physical run closes at the LOCALIZED event time —
//! cinematic continuation is presentation, never physics. This module
//! is the fs-flyer primitive; the offline fs-contact pass stays the
//! certified referee.

use crate::prelaunch::TerrainGrid;
use crate::{Refusal, simloop::TerminalEvent};
use fs_airscrew::Rotor;
use fs_blake3::hash_domain;
use fs_math::det;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// Proxy construction algorithm version (enters the artifact id).
pub const PROXY_ALGO_VERSION: &str = "blade-collision-proxy-v1";

/// Spanwise rails at fractions of local chord.
pub const RAILS: [f64; 2] = [0.25, 0.75];

/// Baseline radial intervals per blade (16 capsules with 2 rails).
pub const BASE_INTERVALS: usize = 8;

/// Refined radial intervals per blade (24 capsules with 2 rails).
pub const REFINED_INTERVALS: usize = 12;

/// Declared construction margins [m]: digitization, station
/// interpolation, numerical (summed into every generated radius).
pub const MARGIN_M: f64 = 0.002 + 0.002 + 0.0005;

/// Cover-excess floor [m] (plan: excess ≤ max(5 mm, geometry
/// uncertainty)).
pub const EXCESS_FLOOR_M: f64 = 0.005;

/// Radial segment cap as a fraction of R.
pub const MAX_SEGMENT_FRAC: f64 = 0.125;

/// Hub/drivetrain proxy radius [m] (separate from the blade claims).
pub const HUB_PROXY_RADIUS_M: f64 = 0.15;

/// Disk-clearance warning threshold [m] (conservative, never terminal).
pub const DISK_WARN_CLEARANCE_M: f64 = 0.05;

/// Rotation-rate admission cap [rad/s].
pub const MAX_OMEGA_RAD_S: f64 = 1_000.0;

/// Bisection iterations (fixed count — deterministic).
const BISECT_ITERS: usize = 48;

/// Phase step per swept substep [rad] (drives the substep count).
const PHASE_STEP_RAD: f64 = 0.05;

/// Substep hard cap.
pub const MAX_SUBSTEPS: usize = 4_096;

/// Local chord and twist at a radial fraction (linear between the
/// registered stations, clamped at the ends — the E1.6 table is the
/// provenance-bound source; no other geometry input exists).
fn chord_beta_at(rotor: &Rotor, r_over_r: f64) -> (f64, f64) {
    let st = &rotor.stations;
    if r_over_r <= st[0].r_over_r {
        return (st[0].chord_m, st[0].beta_rad);
    }
    for w in st.windows(2) {
        if r_over_r <= w[1].r_over_r {
            let t = (r_over_r - w[0].r_over_r) / (w[1].r_over_r - w[0].r_over_r);
            return (
                w[0].chord_m + t * (w[1].chord_m - w[0].chord_m),
                w[0].beta_rad + t * (w[1].beta_rad - w[0].beta_rad),
            );
        }
    }
    let last = st.last().expect("admitted rotor has stations");
    (last.chord_m, last.beta_rad)
}

/// A point on the blade surface in ROTOR frame at azimuth `phi`:
/// axis = +x, radial direction e_r(phi) = (0, cos phi, sin phi),
/// tangential e_t = axis × e_r, section rotated by beta about e_r,
/// chord fraction `f` measured from the leading rail datum (the
/// quarter-chord sits on the 25 % rail).
fn blade_point(rotor: &Rotor, r_over_r: f64, f: f64, phi: f64) -> [f64; 3] {
    let (c, beta) = chord_beta_at(rotor, r_over_r);
    let r = r_over_r * rotor.radius_m;
    let (sp, cp) = (det::sin(phi), det::cos(phi));
    let e_r = [0.0, cp, sp];
    let e_t = [0.0, -sp, cp];
    let (sb, cb) = (det::sin(beta), det::cos(beta));
    let s = c * (f - 0.25);
    [
        r * e_r[0] + s * (cb * e_t[0] + sb),
        r * e_r[1] + s * cb * e_t[1],
        r * e_r[2] + s * cb * e_t[2],
    ]
}

fn dist_point_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let l2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let t = if l2 > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let q = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// One capsule in the blade TEMPLATE frame (azimuth 0): axis endpoints
/// as (r/R, chord fraction) pairs plus the generated radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BladeCapsule {
    /// Inner axis endpoint radial fraction.
    pub r0_over_r: f64,
    /// Outer axis endpoint radial fraction.
    pub r1_over_r: f64,
    /// Rail chord fraction.
    pub rail_f: f64,
    /// GENERATED radius [m] (max assigned-sample distance + margins).
    pub radius_m: f64,
}

/// The cover certificate (recomputed against an INDEPENDENT, finer
/// validation sample set — construction samples never certify
/// themselves).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoverCertificate {
    /// Capsules per blade.
    pub capsules_per_blade: usize,
    /// Every validation sample inside its capsule?
    pub full_cover: bool,
    /// Worst uncovered margin [m] (0 when fully covered).
    pub worst_uncovered_m: f64,
    /// Worst per-capsule excess: radius − max validation need [m].
    pub excess_worst_m: f64,
    /// Longest radial segment as a fraction of R.
    pub seg_len_max_frac: f64,
    /// Smallest radial coordinate any capsule reaches, /R (hub void
    /// exclusion witness — must clear the first-station radius).
    pub hub_min_r_frac: f64,
}

/// The registered blade-collision proxy (Round-4 Q3 contract).
#[derive(Clone, Debug, PartialEq)]
pub struct BladeCollisionProxyV1 {
    /// Template capsules (blade 0; other blades are phase offsets).
    pub capsules: Vec<BladeCapsule>,
    /// Blade count.
    pub n_blades: usize,
    /// Tip radius [m].
    pub radius_m: f64,
    /// The accepted cover certificate.
    pub certificate: CoverCertificate,
    /// Content hash of { source geometry, algorithm version, capsule
    /// count/rails/radii, cover certificate } — the id the manifests
    /// carry (plan §"BladeCollisionProxyArtifactId").
    pub artifact_id: String,
}

fn rotor_digest_bytes(rotor: &Rotor) -> Vec<u8> {
    let mut b = Vec::new();
    for v in [rotor.radius_m, rotor.camber_ratio, rotor.n_blades as f64] {
        b.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    for s in &rotor.stations {
        for v in [s.r_over_r, s.chord_m, s.beta_rad] {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    b
}

fn build_at_resolution(rotor: &Rotor, n_intervals: usize) -> (Vec<BladeCapsule>, CoverCertificate) {
    let r_hub = rotor.stations[0].r_over_r;
    let span = 1.0 - r_hub;
    let mut capsules = Vec::with_capacity(n_intervals * RAILS.len());
    for k in 0..n_intervals {
        let r0 = r_hub + span * k as f64 / n_intervals as f64;
        let r1 = r_hub + span * (k + 1) as f64 / n_intervals as f64;
        for rail in RAILS {
            // Construction samples: 9 radial × 9 chordwise in this
            // interval's rail half of the section.
            let a = blade_point(rotor, r0, rail, 0.0);
            let b = blade_point(rotor, r1, rail, 0.0);
            let mut need: f64 = 0.0;
            for ir in 0..=8 {
                let r = r0 + (r1 - r0) * ir as f64 / 8.0;
                for ic in 0..=8 {
                    let f = if rail < 0.5 {
                        0.5 * ic as f64 / 8.0
                    } else {
                        0.5 + 0.5 * ic as f64 / 8.0
                    };
                    let p = blade_point(rotor, r, f, 0.0);
                    need = need.max(dist_point_segment(p, a, b));
                }
            }
            capsules.push(BladeCapsule {
                r0_over_r: r0,
                r1_over_r: r1,
                rail_f: rail,
                radius_m: need + MARGIN_M,
            });
        }
    }
    // Certificate from an INDEPENDENT validation sweep: 13 radial x
    // 11 chordwise per interval, offset from the construction grid.
    // Excess is judged PER RADIAL STATION (radius minus the local
    // need there): a constant-radius capsule over a varying chord
    // sticks out at the thin stations, and THAT is the conservatism
    // the 5 mm bound caps — an aggregate max would hide it.
    let mut full_cover = true;
    let mut worst_uncovered: f64 = 0.0;
    let mut excess_worst: f64 = 0.0;
    for cap in &capsules {
        let a = blade_point(rotor, cap.r0_over_r, cap.rail_f, 0.0);
        let b = blade_point(rotor, cap.r1_over_r, cap.rail_f, 0.0);
        for ir in 0..13 {
            let r = cap.r0_over_r + (cap.r1_over_r - cap.r0_over_r) * (ir as f64 + 0.37) / 13.0;
            let mut local_need: f64 = 0.0;
            // Chordwise validation INCLUDES the extremes (f = 0 and
            // 0.5 half-chord bounds) — otherwise a sampling gap of
            // ~0.03·c pollutes the excess floor and hides the real
            // radial-taper conservatism the 5 mm bound is for.
            for ic in 0..=10 {
                let half = ic as f64 / 10.0 * 0.5;
                let f = if cap.rail_f < 0.5 { half } else { 0.5 + half };
                let p = blade_point(rotor, r, f, 0.0);
                local_need = local_need.max(dist_point_segment(p, a, b));
            }
            if local_need > cap.radius_m {
                full_cover = false;
                worst_uncovered = worst_uncovered.max(local_need - cap.radius_m);
            }
            excess_worst = excess_worst.max(cap.radius_m - local_need);
        }
    }
    let cert = CoverCertificate {
        capsules_per_blade: capsules.len(),
        full_cover,
        worst_uncovered_m: worst_uncovered,
        excess_worst_m: excess_worst,
        seg_len_max_frac: span / n_intervals as f64,
        hub_min_r_frac: r_hub,
    };
    (capsules, cert)
}

fn cert_accepts(cert: &CoverCertificate, geometry_uncertainty_m: f64) -> bool {
    cert.full_cover
        && cert.excess_worst_m <= EXCESS_FLOOR_M.max(geometry_uncertainty_m) + MARGIN_M
        && cert.seg_len_max_frac <= MAX_SEGMENT_FRAC
        && cert.hub_min_r_frac > 0.0
}

/// Battery/diagnostic surface: the certificate each resolution would
/// produce (never an acceptance path).
#[doc(hidden)]
#[must_use]
pub fn certs_at_resolutions(rotor: &Rotor) -> Vec<CoverCertificate> {
    [BASE_INTERVALS, REFINED_INTERVALS]
        .iter()
        .map(|n| build_at_resolution(rotor, *n).1)
        .collect()
}

/// Battery/diagnostic surface: per-capsule (r0, r1, rail, radius) at
/// a resolution (never an acceptance path).
#[doc(hidden)]
#[must_use]
pub fn capsules_at_resolution(rotor: &Rotor, n_intervals: usize) -> Vec<BladeCapsule> {
    build_at_resolution(rotor, n_intervals).0
}

/// Build the proxy: baseline 16/blade, deterministic refinement to 24
/// when the certificate fails, typed refusal beyond.
///
/// # Errors
/// `blade-proxy-invalid` (rotor refusals or non-finite uncertainty);
/// `blade-cover-uncertifiable` (> 24 capsules/blade would be needed —
/// the blade-strike claim refuses; only the disk warning remains).
pub fn build_blade_proxy(
    rotor: &Rotor,
    geometry_uncertainty_m: f64,
) -> Result<BladeCollisionProxyV1, Refusal> {
    rotor.admit().map_err(|e| {
        refuse(
            "blade-proxy-invalid",
            format!("rotor: {}", e.message),
            "admit the E1.6 rotor first",
        )
    })?;
    if !(geometry_uncertainty_m.is_finite() && geometry_uncertainty_m >= 0.0) {
        return Err(refuse(
            "blade-proxy-invalid",
            format!("geometry uncertainty {geometry_uncertainty_m}"),
            "finite non-negative uncertainty",
        ));
    }
    for n_intervals in [BASE_INTERVALS, REFINED_INTERVALS] {
        let (capsules, certificate) = build_at_resolution(rotor, n_intervals);
        if cert_accepts(&certificate, geometry_uncertainty_m) {
            let mut b = rotor_digest_bytes(rotor);
            b.extend_from_slice(PROXY_ALGO_VERSION.as_bytes());
            for c in &capsules {
                for v in [c.r0_over_r, c.r1_over_r, c.rail_f, c.radius_m] {
                    b.extend_from_slice(&v.to_bits().to_le_bytes());
                }
            }
            for v in [
                certificate.excess_worst_m,
                certificate.seg_len_max_frac,
                certificate.hub_min_r_frac,
            ] {
                b.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            let artifact_id =
                hash_domain("org.frankensim.wf.blade-collision-proxy.v1", &b).to_hex();
            return Ok(BladeCollisionProxyV1 {
                capsules,
                n_blades: rotor.n_blades as usize,
                radius_m: rotor.radius_m,
                certificate,
                artifact_id,
            });
        }
    }
    Err(refuse(
        "blade-cover-uncertifiable",
        "cover certificate fails at 24 capsules/blade".into(),
        "the blade-strike claim refuses; keep the conservative disk warning only",
    ))
}

/// Rigid swept motion of the prop over one contact step: hub
/// translates linearly p0 → p1 while the rotor spins at omega from
/// theta0 (rotor angle and rate are checkpointed states).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweptPropMotion {
    /// Hub position at t = 0 [m].
    pub hub0_m: [f64; 3],
    /// Hub position at t = dt [m].
    pub hub1_m: [f64; 3],
    /// Rotor angle at t = 0 [rad].
    pub theta0_rad: f64,
    /// Rotor rate [rad/s].
    pub omega_rad_s: f64,
    /// Step length [s].
    pub dt_s: f64,
}

impl SweptPropMotion {
    /// Admit the motion.
    ///
    /// # Errors
    /// `swept-motion-invalid` (non-finite fields, dt ≤ 0, |omega|
    /// beyond the cap — AT the cap admits, beyond refuses).
    pub fn admit(&self) -> Result<(), Refusal> {
        let finite = self
            .hub0_m
            .iter()
            .chain(self.hub1_m.iter())
            .all(|v| v.is_finite())
            && self.theta0_rad.is_finite()
            && self.omega_rad_s.is_finite()
            && self.dt_s.is_finite();
        if finite && self.dt_s > 0.0 && self.omega_rad_s.abs() <= MAX_OMEGA_RAD_S {
            Ok(())
        } else {
            Err(refuse(
                "swept-motion-invalid",
                format!("{self:?}"),
                "finite motion, dt > 0, |omega| <= 1000 rad/s",
            ))
        }
    }
}

/// The swept outcome for one step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SweptOutcome {
    /// No feature event this step.
    Clear,
    /// A blade capsule reached the terrain at the localized time.
    BladeStrike {
        /// Localized event time within the step [s].
        t_event_s: f64,
        /// Blade index.
        blade: usize,
        /// Capsule index within the blade template.
        capsule: usize,
    },
    /// The hub/drivetrain proxy reached the terrain.
    HubStrike {
        /// Localized event time within the step [s].
        t_event_s: f64,
    },
}

impl SweptOutcome {
    /// TerminalEvent semantics: a strike whose continuation would need
    /// an unavailable breakage model closes the physical run.
    #[must_use]
    pub fn terminal(&self) -> Option<TerminalEvent> {
        match self {
            SweptOutcome::Clear => None,
            SweptOutcome::BladeStrike { .. } | SweptOutcome::HubStrike { .. } => {
                Some(TerminalEvent::DamageModelUnavailable)
            }
        }
    }
}

/// The conservative disk-envelope WARNING (separate channel — never
/// terminal, never a blade claim).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiskWarning {
    /// Minimum disk-envelope clearance seen this step [m].
    pub min_clearance_m: f64,
    /// Substep time where it occurred [s].
    pub t_s: f64,
}

/// The swept-step report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweptReport {
    /// Terminal-capable outcome.
    pub outcome: SweptOutcome,
    /// Disk warning, if the envelope came within the threshold.
    pub disk_warning: Option<DiskWarning>,
}

fn hub_at(m: &SweptPropMotion, t: f64) -> [f64; 3] {
    let s = t / m.dt_s;
    [
        m.hub0_m[0] + s * (m.hub1_m[0] - m.hub0_m[0]),
        m.hub0_m[1] + s * (m.hub1_m[1] - m.hub0_m[1]),
        m.hub0_m[2] + s * (m.hub1_m[2] - m.hub0_m[2]),
    ]
}

/// Signed clearance of one capsule at time t (min over 5 axis samples
/// of point clearance minus radius).
fn capsule_clearance(
    rotor: &Rotor,
    cap: &BladeCapsule,
    m: &SweptPropMotion,
    blade: usize,
    n_blades: usize,
    terrain: &TerrainGrid,
    t: f64,
) -> Result<f64, Refusal> {
    let phi =
        m.theta0_rad + m.omega_rad_s * t + core::f64::consts::TAU * blade as f64 / n_blades as f64;
    let hub = hub_at(m, t);
    let a = blade_point(rotor, cap.r0_over_r, cap.rail_f, phi);
    let b = blade_point(rotor, cap.r1_over_r, cap.rail_f, phi);
    let mut min_c = f64::INFINITY;
    for i in 0..5 {
        let s = i as f64 / 4.0;
        let p = [
            hub[0] + a[0] + s * (b[0] - a[0]),
            hub[1] + a[1] + s * (b[1] - a[1]),
            hub[2] + a[2] + s * (b[2] - a[2]),
        ];
        let h = terrain.height_m(p[0], p[1])?;
        min_c = min_c.min(p[2] - h - cap.radius_m);
    }
    Ok(min_c)
}

fn hub_clearance(m: &SweptPropMotion, terrain: &TerrainGrid, t: f64) -> Result<f64, Refusal> {
    let hub = hub_at(m, t);
    let h = terrain.height_m(hub[0], hub[1])?;
    Ok(hub[2] - h - HUB_PROXY_RADIUS_M)
}

fn disk_clearance(
    radius_m: f64,
    m: &SweptPropMotion,
    terrain: &TerrainGrid,
    t: f64,
) -> Result<f64, Refusal> {
    let hub = hub_at(m, t);
    let mut min_c = f64::INFINITY;
    for i in 0..24 {
        let phi = core::f64::consts::TAU * i as f64 / 24.0;
        let p = [
            hub[0],
            hub[1] + radius_m * det::cos(phi),
            hub[2] + radius_m * det::sin(phi),
        ];
        let h = terrain.height_m(p[0], p[1])?;
        min_c = min_c.min(p[2] - h);
    }
    Ok(min_c)
}

/// Phase-resolved swept localization for one 1/240 s contact step.
/// Substep count follows the phase rate (never fewer than 8); the
/// earliest negative-clearance feature is bisected to the event time
/// with a fixed iteration count (deterministic).
///
/// # Errors
/// Motion refusals; terrain-domain refusals propagate (`swept` never
/// invents heights); `swept-substeps-exhausted` if the phase rate
/// demands more than the cap.
pub fn swept_prop_step(
    rotor: &Rotor,
    proxy: &BladeCollisionProxyV1,
    motion: &SweptPropMotion,
    terrain: &TerrainGrid,
) -> Result<SweptReport, Refusal> {
    motion.admit()?;
    let needed = ((motion.omega_rad_s.abs() * motion.dt_s / PHASE_STEP_RAD).ceil() as usize).max(8);
    if needed > MAX_SUBSTEPS {
        return Err(refuse(
            "swept-substeps-exhausted",
            format!("{needed} substeps > {MAX_SUBSTEPS}"),
            "shorten the step or the rate",
        ));
    }
    // The minimum clearance across every blade capsule + hub at t.
    let worst_at = |t: f64| -> Result<(f64, SweptOutcome), Refusal> {
        let mut worst = f64::INFINITY;
        let mut what = SweptOutcome::Clear;
        for blade in 0..proxy.n_blades {
            for (ci, cap) in proxy.capsules.iter().enumerate() {
                let c = capsule_clearance(rotor, cap, motion, blade, proxy.n_blades, terrain, t)?;
                if c < worst {
                    worst = c;
                    what = SweptOutcome::BladeStrike {
                        t_event_s: t,
                        blade,
                        capsule: ci,
                    };
                }
            }
        }
        let hc = hub_clearance(motion, terrain, t)?;
        if hc < worst {
            worst = hc;
            what = SweptOutcome::HubStrike { t_event_s: t };
        }
        Ok((worst, what))
    };
    // Disk warning sweep (separate channel).
    let mut disk_warning: Option<DiskWarning> = None;
    for i in 0..=needed {
        let t = motion.dt_s * i as f64 / needed as f64;
        let dc = disk_clearance(proxy.radius_m, motion, terrain, t)?;
        if dc < DISK_WARN_CLEARANCE_M && disk_warning.is_none_or(|w| dc < w.min_clearance_m) {
            disk_warning = Some(DiskWarning {
                min_clearance_m: dc,
                t_s: t,
            });
        }
    }
    // Event scan + fixed-count bisection.
    let (c0, _) = worst_at(0.0)?;
    if c0 <= 0.0 {
        let (_, what) = worst_at(0.0)?;
        return Ok(SweptReport {
            outcome: what,
            disk_warning,
        });
    }
    let mut prev_t = 0.0;
    for i in 1..=needed {
        let t = motion.dt_s * i as f64 / needed as f64;
        let (c, _) = worst_at(t)?;
        if c <= 0.0 {
            let (mut lo, mut hi) = (prev_t, t);
            for _ in 0..BISECT_ITERS {
                let mid = 0.5 * (lo + hi);
                let (cm, _) = worst_at(mid)?;
                if cm <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let (_, what) = worst_at(hi)?;
            let outcome = match what {
                SweptOutcome::BladeStrike { blade, capsule, .. } => SweptOutcome::BladeStrike {
                    t_event_s: hi,
                    blade,
                    capsule,
                },
                SweptOutcome::HubStrike { .. } => SweptOutcome::HubStrike { t_event_s: hi },
                SweptOutcome::Clear => SweptOutcome::Clear,
            };
            return Ok(SweptReport {
                outcome,
                disk_warning,
            });
        }
        prev_t = t;
    }
    Ok(SweptReport {
        outcome: SweptOutcome::Clear,
        disk_warning,
    })
}

/// Swept SEGMENT feature (skid, canard frame, wingtip-as-degenerate)
/// vs the heightfield over one step: endpoints move linearly, interior
/// sampled at 5 axis points, scan + fixed-count bisection. Returns the
/// localized event time, or None when clear.
///
/// # Errors
/// `swept-feature-invalid` (non-finite geometry, dt ≤ 0, negative
/// radius); terrain-domain refusals propagate.
pub fn swept_feature_event(
    seg0: ([f64; 3], [f64; 3]),
    seg1: ([f64; 3], [f64; 3]),
    radius_m: f64,
    dt_s: f64,
    terrain: &TerrainGrid,
) -> Result<Option<f64>, Refusal> {
    let finite = seg0
        .0
        .iter()
        .chain(seg0.1.iter())
        .chain(seg1.0.iter())
        .chain(seg1.1.iter())
        .all(|v| v.is_finite());
    if !finite || !(dt_s > 0.0 && dt_s.is_finite()) || !(radius_m >= 0.0) {
        return Err(refuse(
            "swept-feature-invalid",
            "non-finite segment, bad dt, or negative radius".into(),
            "finite endpoints, dt > 0, radius >= 0",
        ));
    }
    let clearance = |t: f64| -> Result<f64, Refusal> {
        let s = t / dt_s;
        let lerp = |a: [f64; 3], b: [f64; 3]| {
            [
                a[0] + s * (b[0] - a[0]),
                a[1] + s * (b[1] - a[1]),
                a[2] + s * (b[2] - a[2]),
            ]
        };
        let a = lerp(seg0.0, seg1.0);
        let b = lerp(seg0.1, seg1.1);
        let mut min_c = f64::INFINITY;
        for i in 0..5 {
            let u = i as f64 / 4.0;
            let p = [
                a[0] + u * (b[0] - a[0]),
                a[1] + u * (b[1] - a[1]),
                a[2] + u * (b[2] - a[2]),
            ];
            let h = terrain.height_m(p[0], p[1])?;
            min_c = min_c.min(p[2] - h - radius_m);
        }
        Ok(min_c)
    };
    if clearance(0.0)? <= 0.0 {
        return Ok(Some(0.0));
    }
    const SCAN: usize = 16;
    let mut prev_t = 0.0;
    for i in 1..=SCAN {
        let t = dt_s * i as f64 / SCAN as f64;
        if clearance(t)? <= 0.0 {
            let (mut lo, mut hi) = (prev_t, t);
            for _ in 0..BISECT_ITERS {
                let mid = 0.5 * (lo + hi);
                if clearance(mid)? <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return Ok(Some(hi));
        }
        prev_t = t;
    }
    Ok(None)
}
