//! E4.6b-ii battery (bead wf-root-guzez.5.14.2, V-02b qualitative
//! gates): equilibrium SET vs an analytic cubic (per-branch oracles),
//! set-valued stiction intervals + sliding directions, the OVERBALANCE
//! sign/tendency gate on the REAL coupled solution (mid-band hinge is
//! self-driving near center, matching Orville's diary; the forward-hinge
//! falsifier twin shows the tendency ABSENT), sample caps at cap AND cap+1,
//! determinism, golden.
//! Repro: cargo test -p fs-flyer --test freecontrol_battery

use fs_flyer::Refusal;
use fs_flyer::freecontrol::{
    BranchStability, MAX_SWEEP_SAMPLES, MIN_SWEEP_SAMPLES, SweepSpec, free_control_analysis,
};
use fs_wing::hinge::{HingeAxis, SectionCouple, hinge_load};
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec, solve_nonlinear};
use fs_wing::{SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v02b\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;

/// The coupled wing+canard hinge moment as a function of canard
/// deflection at a frozen flight state (prop wash omitted for the
/// released-control sweep — the hinge tendency is a canard-local
/// question; documented tier choice). Hinge axis at `pct`% chord.
fn coupled_hinge_moment(pct: f64) -> impl Fn(f64) -> Result<f64, Refusal> {
    move |delta: f64| -> Result<f64, Refusal> {
        let map_err = |e: fs_wing::Refusal| Refusal {
            code: e.code,
            message: e.message,
            ranked_repairs: e.ranked_repairs,
        };
        let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
        p.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
        let base_canard = p.len();
        p.extend(flat_surface(SurfaceId::CanardLower, 3.66, 0.61, 2.9, 1.05, 4, 1).unwrap());
        p.extend(flat_surface(SurfaceId::CanardUpper, 3.66, 0.61, 2.9, 0.35, 4, 1).unwrap());
        let mut strips = Vec::new();
        for plane in 0..2 {
            let base = plane * 16;
            for s in 0..8 {
                strips.push(StripSpec {
                    panel_indices: vec![base + s, base + 8 + s],
                    chord_m: 1.981,
                    twist_rad: 0.0,
                });
            }
        }
        for plane in 0..2 {
            let base = base_canard + plane * 4;
            for s in 0..4 {
                strips.push(StripSpec {
                    panel_indices: vec![base + s],
                    chord_m: 0.61,
                    twist_rad: delta,
                });
            }
        }
        let alpha = 0.06f64;
        let fs_v = [V * alpha.cos(), 0.0, V * alpha.sin()];
        let op = InfluenceOperator::build(&p, fs_v, RHO).map_err(map_err)?;
        let closure = |_s: usize, a: f64| -> (f64, StripRegime) {
            let attached = core::f64::consts::TAU * (a + 0.1);
            let abs = a.abs();
            if abs <= 0.30 {
                (attached, StripRegime::Attached)
            } else if abs < 0.45 {
                let t = (abs - 0.30) / 0.15;
                let s = t * t * (3.0 - 2.0 * t);
                let sep = 1.98 * a.sin() * a.cos();
                (attached * (1.0 - s) + sep * s, StripRegime::Blended)
            } else {
                (1.98 * a.sin() * a.cos(), StripRegime::Separated)
            }
        };
        let sol =
            solve_nonlinear(&op, &p, &strips, fs_v, RHO, &closure, None, None).map_err(map_err)?;
        let axis = HingeAxis {
            point_m: [2.9 - pct / 100.0 * 0.61, 0.0, 0.7],
            axis_unit: [0.0, 1.0, 0.0],
        };
        // Section cm0 couples for the 8 canard strips.
        let q = 0.5 * RHO * V * V;
        let cm0 = -core::f64::consts::PI * 0.05;
        let m_per_strip = q * 0.61 * 0.61 * 0.915 * cm0;
        let couples: Vec<SectionCouple> = (0..8)
            .map(|_| SectionCouple {
                moment_nm: m_per_strip,
                span_unit: [0.0, 1.0, 0.0],
            })
            .collect();
        let rep = hinge_load(
            &p,
            &sol.gamma,
            fs_v,
            RHO,
            &[SurfaceId::CanardLower, SurfaceId::CanardUpper],
            &axis,
            &couples,
        )
        .map_err(map_err)?;
        Ok(rep.total_nm)
    }
}

/// The classic overbalance-shaped cubic: roots at −r, 0, +r; center
/// unstable, outer roots stable.
fn cubic(delta: f64) -> Result<f64, Refusal> {
    let (a, r) = (400.0, 0.3);
    Ok(-a * delta * (delta - r) * (delta + r))
}

#[test]
fn equilibrium_set_matches_the_analytic_cubic_per_branch() {
    let spec = SweepSpec {
        delta_min_rad: -0.5,
        delta_max_rad: 0.5,
        samples: 101,
        stiction_nm: 0.0,
    };
    let rep = free_control_analysis(&cubic, &spec).unwrap();
    assert_eq!(rep.equilibria.len(), 3, "the SET has three branches");
    let expect = [
        (-0.3, BranchStability::Stable),
        (0.0, BranchStability::Unstable),
        (0.3, BranchStability::Stable),
    ];
    for (e, (d_ref, s_ref)) in rep.equilibria.iter().zip(expect) {
        assert!(
            (e.delta_rad - d_ref).abs() < 1e-9,
            "root {} vs {d_ref}",
            e.delta_rad
        );
        assert_eq!(e.stability, s_ref, "branch at {d_ref}");
        // Analytic slope: d/dδ[−a·δ(δ²−r²)] = −a(3δ²−r²).
        let slope_ref = -400.0 * (3.0 * d_ref * d_ref - 0.09);
        assert!(
            (e.slope_nm_per_rad - slope_ref).abs() < 1e-3 * slope_ref.abs().max(1.0),
            "slope {} vs {slope_ref}",
            e.slope_nm_per_rad
        );
    }
    // Branch ids are distinct.
    assert_ne!(rep.equilibria[0].branch_id, rep.equilibria[1].branch_id);
    assert_eq!(rep.self_driving_near_center, Some(true));
    jlog("cubic-set", &format!("\"n\":{}", rep.equilibria.len()));
}

#[test]
fn stiction_set_and_sliding_directions_are_set_valued() {
    let spec = SweepSpec {
        delta_min_rad: -0.5,
        delta_max_rad: 0.5,
        samples: 401,
        stiction_nm: 2.0,
    };
    let rep = free_control_analysis(&cubic, &spec).unwrap();
    // Three stiction intervals, one around each root.
    assert_eq!(rep.stiction.len(), 3, "stiction SET: {:?}", rep.stiction);
    for (iv, root) in rep.stiction.iter().zip([-0.3, 0.0, 0.3]) {
        assert!(
            iv.lo_rad < root && root < iv.hi_rad,
            "interval {iv:?} must bracket {root}"
        );
    }
    // Four sliding segments with alternating directions [+1,−1,+1,−1]:
    // each drives AWAY from the unstable center toward a stable branch.
    assert_eq!(rep.sliding.len(), 4, "sliding: {:?}", rep.sliding);
    let dirs: Vec<i8> = rep.sliding.iter().map(|s| s.direction).collect();
    assert_eq!(dirs, vec![1, -1, 1, -1]);
    jlog(
        "stiction-sliding",
        &format!("\"stiction_n\":{},\"dirs\":{dirs:?}", rep.stiction.len()),
    );
}

#[test]
fn overbalance_tendency_on_the_real_coupled_solution() {
    // V-02b sign/tendency gate (the ONLY licensed claim class): with the
    // hinge mid-band (40% chord — 'balanced too near the center'), the
    // released canard is SELF-DRIVING near center; at the forward prior
    // edge (25%, on the qc line) the tendency is ABSENT.
    let spec = SweepSpec {
        delta_min_rad: -0.35,
        delta_max_rad: 0.35,
        samples: 29,
        stiction_nm: 0.5,
    };
    let mid = free_control_analysis(&coupled_hinge_moment(40.0), &spec).unwrap();
    assert!(
        !mid.equilibria.is_empty(),
        "the released mid-band canard must have equilibria in travel"
    );
    assert_eq!(
        mid.self_driving_near_center,
        Some(true),
        "mid-band hinge must reproduce the diary tendency: {:?}",
        mid.equilibria
    );
    let fwd = free_control_analysis(&coupled_hinge_moment(25.0), &spec).unwrap();
    // Measured: at the forward prior edge the circulatory arm vanishes
    // and the (constant, nose-down) camber couple dominates one-signed —
    // NO equilibrium in travel; the released surface floats to a stop
    // instead of hunting about center. Either way the diary tendency is
    // ABSENT, which is what the twin must show.
    assert_ne!(
        fwd.self_driving_near_center,
        Some(true),
        "forward hinge must NOT be self-driving (falsifier twin): {:?}",
        fwd.equilibria
    );
    if fwd.equilibria.is_empty() {
        // One-signed moment: every sliding segment drifts the same way.
        assert!(!fwd.sliding.is_empty());
        let d0 = fwd.sliding[0].direction;
        assert!(
            fwd.sliding.iter().all(|s| s.direction == d0),
            "no-equilibrium case must slide one-directionally: {:?}",
            fwd.sliding
        );
    }
    jlog(
        "overbalance",
        &format!(
            "\"mid_eq\":{},\"mid_selfdrive\":true,\"fwd_eq\":{},\"fwd_tendency_absent\":true",
            mid.equilibria.len(),
            fwd.equilibria.len()
        ),
    );
}

#[test]
fn sweep_caps_at_cap_and_cap_plus_one() {
    let mk = |samples: usize, stiction: f64| SweepSpec {
        delta_min_rad: -0.5,
        delta_max_rad: 0.5,
        samples,
        stiction_nm: stiction,
    };
    assert!(free_control_analysis(&cubic, &mk(MIN_SWEEP_SAMPLES, 0.0)).is_ok());
    assert_eq!(
        free_control_analysis(&cubic, &mk(MIN_SWEEP_SAMPLES - 1, 0.0))
            .unwrap_err()
            .code,
        "sweep-spec-invalid"
    );
    assert!(free_control_analysis(&cubic, &mk(MAX_SWEEP_SAMPLES, 0.0)).is_ok());
    assert_eq!(
        free_control_analysis(&cubic, &mk(MAX_SWEEP_SAMPLES + 1, 0.0))
            .unwrap_err()
            .code,
        "sweep-spec-invalid"
    );
    assert_eq!(
        free_control_analysis(&cubic, &mk(64, -1e-300))
            .unwrap_err()
            .code,
        "sweep-spec-invalid"
    );
    let bad_bounds = SweepSpec {
        delta_min_rad: 0.5,
        delta_max_rad: -0.5,
        samples: 64,
        stiction_nm: 0.0,
    };
    assert_eq!(
        free_control_analysis(&cubic, &bad_bounds).unwrap_err().code,
        "sweep-spec-invalid"
    );
    let nonfinite = |_d: f64| -> Result<f64, Refusal> { Ok(f64::NAN) };
    assert_eq!(
        free_control_analysis(&nonfinite, &mk(64, 0.0))
            .unwrap_err()
            .code,
        "sweep-moment-nonfinite"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn determinism_and_golden() {
    let spec = SweepSpec {
        delta_min_rad: -0.35,
        delta_max_rad: 0.35,
        samples: 29,
        stiction_nm: 0.5,
    };
    let a = free_control_analysis(&coupled_hinge_moment(40.0), &spec).unwrap();
    let b = free_control_analysis(&coupled_hinge_moment(40.0), &spec).unwrap();
    assert_eq!(a, b, "bitwise repeat");
    let mut payload = Vec::new();
    for e in &a.equilibria {
        payload.extend_from_slice(&e.delta_rad.to_bits().to_le_bytes());
        payload.extend_from_slice(&e.slope_nm_per_rad.to_bits().to_le_bytes());
    }
    for &(d, m) in &a.sweep {
        payload.extend_from_slice(&d.to_bits().to_le_bytes());
        payload.extend_from_slice(&m.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v02b-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "4726ee130dc831c58b8c425ed279afc9cb818376f609e9b5f235c6ad9b49bd39",
        "free-control golden moved — determinism regression or an \
         intentional model change requiring the golden-bump protocol"
    );
}
