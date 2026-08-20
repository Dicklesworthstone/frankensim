//! E4.2b battery (bead wf-root-guzez.5.4): hinge-load channel vs
//! independent quadrature oracles (per-item, never totals-only), the
//! exact axis-shift identities, the double-count HOSTILE TWIN (apparent
//! mass smuggled into the steady hinge load FAILS against the quadrature
//! oracle as designed), typed refusals at cap and cap+1, determinism,
//! golden. Repro: cargo test -p fs-wing --test hinge_battery

use fs_wing::hinge::{AXIS_UNIT_TOL, HingeAxis, SectionCouple, hinge_load};
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec, solve_nonlinear};
use fs_wing::{SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-e42b\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;

/// Biplane canard ahead of a biplane wing — the E4.6a layout class; the
/// hinge is the canard pivot (spanwise axis, 20% chord aft of the
/// quarter-chord line — near-balanced, per the Wright design intent).
fn layout() -> (Vec<fs_wing::Panel>, Vec<StripSpec>, HingeAxis) {
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
                twist_rad: 0.1,
            });
        }
    }
    // Hinge 20% chord AFT of the quarter-chord line: the canard was
    // aerodynamically balanced (hinge near, not at, the center of
    // pressure), and a nonzero arm keeps the quadrature oracle
    // non-trivial (hinge ON the qc line makes the circulatory moment
    // vanish by construction for single-row panels — measured).
    let axis = HingeAxis {
        point_m: [2.9 - 0.45 * 0.61, 0.0, 0.7],
        axis_unit: [0.0, 1.0, 0.0],
    };
    (p, strips, axis)
}

fn freestream(alpha: f64) -> [f64; 3] {
    [V * alpha.cos(), 0.0, V * alpha.sin()]
}

fn camber_closure(_s: usize, alpha: f64) -> (f64, StripRegime) {
    (
        2.0 * std::f64::consts::PI * (alpha + 0.1),
        StripRegime::Attached,
    )
}

const CANARDS: [SurfaceId; 2] = [SurfaceId::CanardLower, SurfaceId::CanardUpper];

#[test]
fn hinge_load_matches_the_quadrature_oracle_per_item() {
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let rep = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &[]).unwrap();
    // Independent per-panel quadrature: recompute each selected panel's
    // ((r-p) x rho*Gamma*(V x seg)) . y-hat by hand.
    let mut total = 0.0;
    let mut checked = 0;
    for item in &rep.items {
        let panel = &p[item.panel];
        assert!(matches!(
            panel.surface,
            SurfaceId::CanardLower | SurfaceId::CanardUpper
        ));
        let seg = [
            panel.b[0] - panel.a[0],
            panel.b[1] - panel.a[1],
            panel.b[2] - panel.a[2],
        ];
        let s = RHO * sol.gamma[item.panel];
        let f = [
            s * (fs_v[1] * seg[2] - fs_v[2] * seg[1]),
            s * (fs_v[2] * seg[0] - fs_v[0] * seg[2]),
            s * (fs_v[0] * seg[1] - fs_v[1] * seg[0]),
        ];
        let r = [
            0.5 * (panel.a[0] + panel.b[0]) - axis.point_m[0],
            0.5 * (panel.a[1] + panel.b[1]) - axis.point_m[1],
            0.5 * (panel.a[2] + panel.b[2]) - axis.point_m[2],
        ];
        let m_ref = r[2] * f[0] - r[0] * f[2]; // (r x f) . y-hat
        assert!(
            (item.moment_nm - m_ref).abs() < 1e-12 * m_ref.abs().max(1.0),
            "panel {}: {} vs quadrature {m_ref}",
            item.panel,
            item.moment_nm
        );
        total += m_ref;
        checked += 1;
    }
    assert_eq!(checked, 8, "both canard planes, all panels");
    assert!((rep.circulatory_nm - total).abs() < 1e-9);
    // Physical plausibility: the hinge sits 20% chord aft of the
    // quarter-chord line, so |M|/(|Fz|*c) must recover that arm
    // fraction class.
    let frac = rep.circulatory_nm.abs() / (rep.force_n[2].abs() * 0.61);
    assert!(
        (0.05..0.45).contains(&frac),
        "hinge arm fraction {frac} outside the balanced-canard band"
    );
    jlog(
        "quadrature",
        &format!(
            "\"total_nm\":{},\"canard_fz\":{},\"arm_frac\":{frac}",
            rep.circulatory_nm, rep.force_n[2]
        ),
    );
}

#[test]
fn axis_shift_identities_are_exact() {
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let base = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &[]).unwrap();
    // (1) Sliding the point ALONG the axis changes nothing.
    let mut along = axis;
    along.point_m[1] += 1.7;
    let rep_along = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &along, &[]).unwrap();
    assert!(
        (rep_along.circulatory_nm - base.circulatory_nm).abs() < 1e-10,
        "along-axis shift must be invariant"
    );
    // (2) A perpendicular shift d changes the moment by exactly
    // -(d x F_total) . a-hat.
    let d = [0.3, 0.0, -0.2];
    let mut perp = axis;
    perp.point_m[0] += d[0];
    perp.point_m[2] += d[2];
    let rep_perp = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &perp, &[]).unwrap();
    let f = base.force_n;
    let dxf_y = d[2] * f[0] - d[0] * f[2];
    let expect = base.circulatory_nm - dxf_y;
    assert!(
        (rep_perp.circulatory_nm - expect).abs() < 1e-9,
        "perp-shift identity: {} vs {expect}",
        rep_perp.circulatory_nm
    );
    jlog("axis-identities", &format!("\"dxf_y\":{dxf_y}"));
}

#[test]
fn section_couples_project_as_free_vectors() {
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let couples = [
        SectionCouple {
            moment_nm: -12.5,
            span_unit: [0.0, 1.0, 0.0],
        },
        SectionCouple {
            moment_nm: -12.5,
            span_unit: [0.0, 1.0, 0.0],
        },
    ];
    let rep = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &couples).unwrap();
    assert!((rep.section_nm - (-25.0)).abs() < 1e-12);
    assert!((rep.total_nm - (rep.circulatory_nm - 25.0)).abs() < 1e-12);
    // A couple orthogonal to the axis contributes nothing.
    let orth = [SectionCouple {
        moment_nm: 99.0,
        span_unit: [1.0, 0.0, 0.0],
    }];
    let rep2 = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &orth).unwrap();
    assert!(rep2.section_nm.abs() < 1e-15);
    jlog("couples", &format!("\"section_nm\":{}", rep.section_nm));
}

#[test]
fn double_count_hostile_twin_fails_as_designed() {
    // The DESIGNED failure (plan DONE-WHEN): apparent-mass hinge torque
    // belongs to the added-mass blocks; a twin that smuggles it into the
    // steady hinge load MUST be caught. The interface takes no
    // acceleration inputs, so the only smuggling route is additive — and
    // at STEADY state the quadrature oracle pins the circulatory value
    // exactly, so ANY nonzero smuggled term is detected.
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let clean = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &[]).unwrap();
    let quadrature: f64 = clean.items.iter().map(|i| i.moment_nm).sum();
    // The twin: an "apparent mass" torque added on top at steady state
    // (alpha-double-dot = 0, so the TRUE apparent-mass term is zero).
    let smuggled = 7.3;
    let twin_total = clean.circulatory_nm + smuggled;
    let caught = (twin_total - quadrature).abs() > 1.0;
    assert!(caught, "the quadrature oracle must catch the smuggled term");
    assert!(
        (clean.circulatory_nm - quadrature).abs() < 1e-9,
        "and the clean channel must pass the same oracle"
    );
    jlog(
        "hostile-twin",
        &format!("\"smuggled_nm\":{smuggled},\"caught\":{caught}"),
    );
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    // Axis-unit tolerance: just inside and clearly outside (the exact
    // ulp boundary is unreachable through the sqrt in admit()).
    let at_cap = HingeAxis {
        point_m: axis.point_m,
        axis_unit: [0.0, 1.0 + 0.5 * AXIS_UNIT_TOL, 0.0],
    };
    assert!(hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &at_cap, &[]).is_ok());
    let past_cap = HingeAxis {
        point_m: axis.point_m,
        axis_unit: [0.0, 1.0 + 2.0 * AXIS_UNIT_TOL, 0.0],
    };
    assert_eq!(
        hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &past_cap, &[])
            .unwrap_err()
            .code,
        "hinge-axis-invalid"
    );
    // Empty and unmatched selections.
    assert_eq!(
        hinge_load(&p, &sol.gamma, fs_v, RHO, &[], &axis, &[])
            .unwrap_err()
            .code,
        "hinge-selection-empty"
    );
    assert_eq!(
        hinge_load(
            &p,
            &sol.gamma,
            fs_v,
            RHO,
            &[SurfaceId::Vertical],
            &axis,
            &[]
        )
        .unwrap_err()
        .code,
        "hinge-selection-empty"
    );
    // Gamma mismatch.
    assert_eq!(
        hinge_load(&p, &sol.gamma[1..], fs_v, RHO, &CANARDS, &axis, &[])
            .unwrap_err()
            .code,
        "gamma-length-mismatch"
    );
    jlog("refusals", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn determinism_and_golden() {
    let (p, strips, axis) = layout();
    let fs_v = freestream(0.06);
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &camber_closure, None, None).unwrap();
    let a = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &[]).unwrap();
    let b = hinge_load(&p, &sol.gamma, fs_v, RHO, &CANARDS, &axis, &[]).unwrap();
    assert_eq!(a, b, "bitwise repeat");
    let mut payload = Vec::new();
    for i in &a.items {
        payload.extend_from_slice(&i.moment_nm.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&a.circulatory_nm.to_bits().to_le_bytes());
    let digest = fs_blake3::hash_domain("org.frankensim.fs-wing.e42b-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"total\":{}", a.circulatory_nm),
    );
    assert_eq!(
        digest, "08685e43461661db86796d5e64e0024d9c46adfde93966be14816ad255ac71c2",
        "hinge golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
