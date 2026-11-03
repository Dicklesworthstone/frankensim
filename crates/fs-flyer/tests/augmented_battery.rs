//! E8.2-i battery (bead wf-root-guzez.9.2.1): eigensolver vs analytic
//! block fixtures AND vs the independent 4-state solver, the
//! block-decoupling union oracle, unstable-pole persistence in the
//! coupled 7-state engine (shift REPORTED), pilot-column activation
//! rule, the typed trim-refusal path naming the limiting subsystem,
//! structural mode attribution, determinism, golden.
//! Repro: cargo test -p fs-flyer --test augmented_battery

use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::augmented::{
    AugmentedEngine, ModeFamily, ROTOR_INERTIA_KG_M2, attribute_modes, build_engine, eig_dense,
    eigenvalues, wrap_trim_refusal,
};
use fs_flyer::canardmech::CANARD_MECH_V1;
use fs_flyer::longitudinal::{IYY_KG_M2, eig4, linearize};
use fs_flyer::pilot::PilotWrightModel;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e82i\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;

fn rigid_and_column() -> (
    fs_flyer::aircraft::OpenLoopDesign,
    fs_flyer::longitudinal::LongitudinalReport,
    [f64; 4],
) {
    rigid_and_column_for(wright_openloop_v1())
}

fn rigid_and_column_for(
    d: fs_flyer::aircraft::OpenLoopDesign,
) -> (
    fs_flyer::aircraft::OpenLoopDesign,
    fs_flyer::longitudinal::LongitudinalReport,
    [f64; 4],
) {
    let trim = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let rep = linearize(&d, &trim, RHO).unwrap();
    let h = 0.004;
    let bp = d
        .force_buildup(
            trim.v_mps,
            trim.alpha_rad,
            trim.delta_canard_rad + h,
            trim.omega_prop_rad_s,
            0.0,
            RHO,
        )
        .unwrap();
    let bm = d
        .force_buildup(
            trim.v_mps,
            trim.alpha_rad,
            trim.delta_canard_rad - h,
            trim.omega_prop_rad_s,
            0.0,
            RHO,
        )
        .unwrap();
    let m = d.gross_mass_kg;
    let b = [
        (bp.force_n[0] - bm.force_n[0]) / (2.0 * h) / m,
        (bp.force_n[2] - bm.force_n[2]) / (2.0 * h) / m,
        (bp.moment_y_nm - bm.moment_y_nm) / (2.0 * h) / IYY_KG_M2,
        0.0,
    ];
    (d, rep, b)
}

#[test]
fn shared_shaft_jacobian_matches_coupled_torque_in_each_physical_coordinate() {
    let mut design = wright_openloop_v1();
    design.disks[0].center_m[2] = -1.35;
    design.disks[1].center_m[2] = -0.25;
    let (design, rigid, control) = rigid_and_column_for(design);
    let engine = build_engine(&design, &CANARD_MECH_V1, None, &rigid, &control, RHO).unwrap();
    let trim = &rigid.trim;
    let state = [
        trim.v_mps * trim.alpha_rad.cos(),
        trim.v_mps * trim.alpha_rad.sin(),
        0.0,
        trim.delta_canard_rad,
        trim.omega_prop_rad_s,
    ];
    let shaft_acceleration = |x: [f64; 5]| {
        let speed = (x[0] * x[0] + x[1] * x[1]).sqrt();
        let alpha = x[1].atan2(x[0]);
        let loads = design
            .force_buildup(speed, alpha, x[3], x[4], x[2], RHO)
            .unwrap();
        // Independent total shaft angular-momentum balance: two engine
        // torque shares and both accepted rotor loads accelerate TWO inertias.
        let total_torque = 2.0 * fs_airscrew::engine_torque_at_prop_nm(x[4])
            - loads.coupled.torque_nm.iter().sum::<f64>();
        total_torque / (2.0 * ROTOR_INERTIA_KG_M2)
    };
    for (coordinate, column, h) in [
        (0, 0, 0.002),
        (1, 1, 0.002),
        (2, 2, 0.0004),
        (3, 4, 0.00016),
        (4, 6, 0.01),
    ] {
        let mut plus = state;
        let mut minus = state;
        plus[coordinate] += h;
        minus[coordinate] -= h;
        // Five-times smaller perturbations than the engine, sampled directly
        // off trim: this checks the nonlinear response, not the same FD code.
        let measured = (shaft_acceleration(plus) - shaft_acceleration(minus)) / (2.0 * h);
        let predicted = engine.a[6][column];
        assert!(
            (measured - predicted).abs() < 2e-4 * (1.0 + measured.abs()),
            "shaft derivative column {column}: {predicted} vs off-trim {measured}"
        );
        assert!(
            measured.abs() > 1e-6,
            "column {column} must expose real coupling"
        );
        jlog(
            "shaft-jacobian-coordinate",
            &format!("\"column\":{column},\"matrix\":{predicted},\"measured\":{measured}"),
        );
    }
    // The supplied rigid and actuator-control owners are retained exactly;
    // attitude is independent of air velocity and has no bare rotor torque.
    for row in 0..4 {
        assert_eq!(&engine.a[row][..4], &rigid.a[row]);
        assert_eq!(engine.a[row][4], control[row]);
    }
    assert_eq!(engine.a[6][3], 0.0);
    assert_eq!(engine.a[6][5], 0.0);
    jlog("shaft-jacobian", "\"all_five_physical_coordinates\":true");
}

#[test]
fn rotor_speed_column_tracks_washed_body_loads_and_both_unequal_hub_arms() {
    let mut design = wright_openloop_v1();
    design.disks[0].center_m[2] = -1.35;
    design.disks[1].center_m[2] = -0.25;
    let (mut design, rigid, control) = rigid_and_column_for(design);
    let engine = build_engine(&design, &CANARD_MECH_V1, None, &rigid, &control, RHO).unwrap();
    let trim = &rigid.trim;
    let h = 0.01;
    let sample = |omega| {
        design
            .force_buildup(
                trim.v_mps,
                trim.alpha_rad,
                trim.delta_canard_rad,
                omega,
                0.0,
                RHO,
            )
            .unwrap()
    };
    let plus = sample(trim.omega_prop_rad_s + h);
    let minus = sample(trim.omega_prop_rad_s - h);
    let measured = [
        (plus.force_n[0] - minus.force_n[0]) / (2.0 * h) / design.gross_mass_kg,
        (plus.force_n[2] - minus.force_n[2]) / (2.0 * h) / design.gross_mass_kg,
        (plus.moment_y_nm - minus.moment_y_nm) / (2.0 * h) / IYY_KG_M2,
    ];
    for row in 0..3 {
        assert!(
            (engine.a[row][6] - measured[row]).abs() < 2e-4 * (1.0 + measured[row].abs()),
            "body row {row}: {} vs off-trim {}",
            engine.a[row][6],
            measured[row]
        );
    }
    assert!(
        measured[1].abs() > 1e-4,
        "prop wash must affect vertical momentum"
    );
    let thrust_only_pitch: f64 = design
        .disks
        .iter()
        .enumerate()
        .map(|(k, disk)| {
            (plus.thrust_n[k] - minus.thrust_n[k]) / (2.0 * h) * (disk.center_m[2] - design.cg_m[2])
                / IYY_KG_M2
        })
        .sum();
    assert!(
        (measured[2] - thrust_only_pitch).abs() > 1e-5,
        "the pitch column must include aerodynamic wash moments, not only thrust arms"
    );
    // Disk labels cannot change physical response. The old single-left-arm
    // approximation violates this on the same unequal-hub aircraft.
    design.disks.swap(0, 1);
    let relabeled = build_engine(&design, &CANARD_MECH_V1, None, &rigid, &control, RHO).unwrap();
    for row in [0, 1, 2, 6] {
        for column in [0, 1, 2, 4, 6] {
            assert!(
                (engine.a[row][column] - relabeled.a[row][column]).abs()
                    < 1e-8 * (1.0 + engine.a[row][column].abs()),
                "disk relabeling changed A[{row}][{column}]"
            );
        }
    }
    jlog(
        "rotor-body-jacobian",
        "\"washed_loads_and_unequal_arms\":true",
    );
}

#[test]
fn eig_dense_matches_analytic_fixtures_and_the_independent_solver() {
    // Analytic 7x7: diagonal reals + one rotation pair.
    let mut a = vec![vec![0.0f64; 7]; 7];
    let (sg, om) = (0.4, 3.0);
    a[0][0] = sg;
    a[0][1] = om;
    a[1][0] = -om;
    a[1][1] = sg;
    for (i, v) in [(-0.5, 2), (-1.0, 3), (-2.0, 4), (-4.0, 5), (1.5, 6)] {
        a[v][v] = i;
    }
    let p = eig_dense(&a).unwrap();
    assert!((p[0].re - 1.5).abs() < 1e-8, "leading real: {:?}", p[0]);
    let pair: Vec<_> = p.iter().filter(|x| x.im.abs() > 1e-9).collect();
    assert_eq!(pair.len(), 2, "one complex pair");
    assert!((pair[0].re - sg).abs() < 1e-8 && (pair[0].im.abs() - om).abs() < 1e-8);
    // Cross-solver: the REAL rigid A through eig_dense vs eig4.
    let (_d, rep, _b) = rigid_and_column();
    let a4: Vec<Vec<f64>> = rep.a.iter().map(|r| r.to_vec()).collect();
    let dense = eig_dense(&a4).unwrap();
    let quart = eig4(&rep.a);
    for (x, y) in dense.iter().zip(quart.iter()) {
        assert!(
            (x.re - y.re).abs() < 1e-7 && (x.im.abs() - y.im.abs()).abs() < 1e-7,
            "cross-solver mismatch: {x:?} vs {y:?}"
        );
    }
    jlog("eig-oracles", "\"analytic_and_cross_solver\":true");
}

#[test]
fn decoupled_engine_is_the_union_of_block_spectra() {
    let (d, rep, b) = rigid_and_column();
    let engine = build_engine(&d, &CANARD_MECH_V1, None, &rep, &b, RHO).unwrap();
    // Freeze ALL couplings: zero off-diagonal blocks.
    let ranges: [std::ops::Range<usize>; 3] = [0..4, 4..6, 6..7];
    let mut frozen = engine.a.clone();
    for i in 0..engine.n {
        for j in 0..engine.n {
            let bi = ranges.iter().position(|r| r.contains(&i));
            let bj = ranges.iter().position(|r| r.contains(&j));
            if bi != bj {
                frozen[i][j] = 0.0;
            }
        }
    }
    let spectrum = eig_dense(&frozen).unwrap();
    // Union oracle: rigid eig4 + actuator {0, -c_eff/I} + rotor a66.
    let mut expected: Vec<(f64, f64)> = eig4(&rep.a).iter().map(|p| (p.re, p.im)).collect();
    let c_eff =
        CANARD_MECH_V1.viscous_nm_s + CANARD_MECH_V1.coulomb_nm / CANARD_MECH_V1.friction_reg_rad_s;
    expected.push((0.0, 0.0));
    expected.push((-c_eff / CANARD_MECH_V1.inertia_kg_m2, 0.0));
    expected.push((engine.a[6][6], 0.0));
    for (re, im) in expected {
        let d_min = spectrum
            .iter()
            .map(|p| ((p.re - re).powi(2) + (p.im.abs() - im.abs()).powi(2)).sqrt())
            .fold(f64::INFINITY, f64::min);
        assert!(
            d_min < 1e-6 * (1.0 + re.abs()),
            "expected block eigenvalue ({re},{im}) missing: nearest {d_min}"
        );
    }
    jlog(
        "block-union",
        &format!("\"actuator_pole\":{}", -c_eff / 6.0),
    );
}

#[test]
fn coupled_engine_keeps_the_unstable_pole_and_reports_the_shift() {
    let (d, rep, b) = rigid_and_column();
    let engine = build_engine(&d, &CANARD_MECH_V1, None, &rep, &b, RHO).unwrap();
    assert_eq!(engine.n, 7);
    let poles = eigenvalues(&engine).unwrap();
    let max_re = poles.iter().map(|p| p.re).fold(f64::NEG_INFINITY, f64::max);
    let rigid_max = rep.max_re;
    assert!(
        max_re > 0.0,
        "the instability must persist in the 7-state engine"
    );
    let shift = max_re / rigid_max - 1.0;
    assert!(
        shift.abs() < 0.3,
        "coupling shifted the unstable pole implausibly: {shift}"
    );
    jlog(
        "coupled",
        &format!(
            "\"rigid_pole\":{rigid_max},\"augmented_pole\":{max_re},\"shift\":{shift},\"declared\":{}",
            engine.declared_approximations.len()
        ),
    );
}

#[test]
fn pilot_columns_activate_only_with_a_pilot() {
    let (d, rep, b) = rigid_and_column();
    let without = build_engine(&d, &CANARD_MECH_V1, None, &rep, &b, RHO).unwrap();
    assert_eq!(without.n, 7);
    assert!(!without.pilot_active);
    assert_eq!(without.labels.len(), 7);
    let pilot = PilotWrightModel::new(2, 9).unwrap();
    let with = build_engine(&d, &CANARD_MECH_V1, Some(&pilot), &rep, &b, RHO).unwrap();
    assert_eq!(with.n, 11);
    assert!(with.pilot_active);
    assert!(with.labels.contains(&"pilot_pade"));
    assert!(
        with.declared_approximations
            .iter()
            .any(|s| s.contains("Pade")),
        "the Pade approximation must be DECLARED"
    );
    let poles = eigenvalues(&with).unwrap();
    assert_eq!(poles.len(), 11);
    // The closed pilot loop's dominant growth must be reported; H-02c
    // showed the latency loop is unstable — assert an RHP pole exists
    // and its magnitude class is the H-02c one (report, wide band).
    let max_re = poles.iter().map(|p| p.re).fold(f64::NEG_INFINITY, f64::max);
    assert!(
        max_re > 0.0 && max_re < 10.0,
        "pilot-loop growth class: {max_re}"
    );
    jlog(
        "pilot-columns",
        &format!("\"n_with\":{},\"pilot_loop_max_re\":{max_re}", with.n),
    );
}

#[test]
fn trim_refusal_names_the_limiting_subsystem() {
    let mut sliver = wright_openloop_v1();
    sliver.canard_span_m = 0.2;
    sliver.canard_chord_m = 0.05;
    let inner = sliver.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap_err();
    let wrapped = wrap_trim_refusal(inner);
    assert_eq!(wrapped.code, "augmented-trim-refused");
    assert!(
        wrapped.message.contains("canard-authority"),
        "must NAME the limiting subsystem: {}",
        wrapped.message
    );
    jlog("trim-refusal", "\"named_subsystem\":true");
}

#[test]
fn structural_attribution_labels_the_expected_families() {
    let (d, rep, b) = rigid_and_column();
    let engine = build_engine(&d, &CANARD_MECH_V1, None, &rep, &b, RHO).unwrap();
    let labeled = attribute_modes(&engine).unwrap();
    assert_eq!(labeled.len(), 7);
    // The heavily-damped actuator pole (~ -c_eff/I ~ -16.8) must be
    // attributed to the Actuator family.
    let act = labeled
        .iter()
        .min_by(|x, y| x.pole.re.partial_cmp(&y.pole.re).unwrap())
        .unwrap();
    // (most negative pole is the fast subsidence — find the one nearest
    // the actuator's analytic value instead)
    let c_eff =
        CANARD_MECH_V1.viscous_nm_s + CANARD_MECH_V1.coulomb_nm / CANARD_MECH_V1.friction_reg_rad_s;
    let target = -c_eff / CANARD_MECH_V1.inertia_kg_m2;
    let nearest = labeled
        .iter()
        .min_by(|x, y| {
            (x.pole.re - target)
                .abs()
                .partial_cmp(&(y.pole.re - target).abs())
                .unwrap()
        })
        .unwrap();
    assert_eq!(
        nearest.family,
        ModeFamily::Actuator,
        "pole near {target} must attribute to Actuator: {nearest:?}"
    );
    let _ = act;
    // Per-item: every label carries a positive attribution shift.
    for l in &labeled {
        assert!(
            l.attribution_shift >= 0.0 && l.attribution_shift.is_finite(),
            "attribution receipt: {l:?}"
        );
    }
    jlog(
        "attribution",
        &format!(
            "\"actuator_pole\":{},\"family\":\"Actuator\"",
            nearest.pole.re
        ),
    );
}

#[test]
fn determinism_and_golden() {
    let (d, rep, b) = rigid_and_column();
    let engine = build_engine(&d, &CANARD_MECH_V1, None, &rep, &b, RHO).unwrap();
    let p1 = eigenvalues(&engine).unwrap();
    let p2 = eigenvalues(&engine).unwrap();
    for (x, y) in p1.iter().zip(p2.iter()) {
        assert_eq!(x.re.to_bits(), y.re.to_bits());
        assert_eq!(x.im.to_bits(), y.im.to_bits());
    }
    let mut payload = Vec::new();
    for p in &p1 {
        payload.extend_from_slice(&p.re.to_bits().to_le_bytes());
        payload.extend_from_slice(&p.im.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e82i-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    // Golden bump 2026-08-21 (bead guzez.7.2.1): det::-routing of
    // Prandtl acos + contact tanh moved the trim point by ulps;
    // the spectrum golden follows the trim it linearizes about.
    // Pinned per-ISA (augmented::eig_dense uses platform f64::hypot at the
    // final 2x1 Givens step; bit-identical across debug and release on each ISA).
    // Golden bump 2026-10-28: the residual-certified native station solve
    // moves trim roundoff and its augmented spectrum. Only x86_64 was
    // remeasured here; the aarch64 pin retains its previous measurement
    // and requires a native aarch64 remeasurement for the station and
    // disk-induced-wash corrections. No aarch64 result is inferred from x86.
    // Golden bump 2026-10-31: the rigid matrix now holds pitch fixed while
    // varying air velocity, removing spurious gravity derivatives. Only the
    // x86_64 spectrum was measured; all physical acceptance bands are intact.
    // Golden bump 2026-11-03: the rigid block now differentiates full local
    // pitch velocity and consistently washed loads. The prior x86_64 digest
    // was reproduced before correction; other ISA pins await remeasurement.
    // Golden bump 2026-11-03: rotor/body derivatives now use accepted coupled
    // aircraft loads and the common-shaft torque balance. The old owner
    // reproduces the preceding x86_64 pin and fails both coupling regressions.
    // Golden bump 2026-11-03: transverse circulatory lift corrects induced
    // drag and the trim spectrum. Only the measured x86_64 pin changes.
    #[cfg(target_arch = "aarch64")]
    const EXPECTED_DIGEST: &str =
        "389e2a7c4f5a97e8c8507e011979571d10cb87d50c3f57d0b06cd38955d29321";
    #[cfg(target_arch = "x86_64")]
    const EXPECTED_DIGEST: &str =
        "8534327a3f9d497e407d6db6805247d80653e50f21f61980a25052ea6d2c766d";
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    const EXPECTED_DIGEST: &str =
        "48828f654100a2407311265ac37a06525e09ed8543aaca5737ddf86ed1756dcc";
    assert_eq!(
        digest, EXPECTED_DIGEST,
        "augmented-spectrum golden moved — determinism regression or an \
         intentional model change requiring the golden-bump protocol"
    );
}

// Silence a false-positive unused warning path for AugmentedEngine in
// some cfg combinations (the type is used above).
#[allow(dead_code)]
fn _type_anchor(e: &AugmentedEngine) -> usize {
    e.n
}
