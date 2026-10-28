//! Resolved turbulence loads must reach roll/yaw without a synthetic gust gain.
use fs_flyer::aircraft::{ForceBuildup, OpenLoopDesign, wright_openloop_v1};
use fs_flyer::lateral::{
    FLYER_ROLL_INERTIA_KG_M2, REDUCED_YAW_INERTIA_KG_M2,
    LateralModel, LateralMomentsNm, LateralState, RudderLinkage,
};

const RHO: f64 = 1.294;
const DC: f64 = 0.05;

fn fixture(shear: f64) -> (OpenLoopDesign, ForceBuildup) {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let velocity = [13.86 * fs_math::det::cos(0.06), 0.0, 13.86 * fs_math::det::sin(0.06)];
    let loads = design.force_buildup_with_incident(
        13.86, 0.06, DC, 40.0, 0.0, RHO, None,
        &|point| Ok([velocity[0], velocity[1], velocity[2] + shear * point[1]]),
    ).unwrap();
    (design, loads)
}

#[test]
fn distributed_moments_match_independent_panel_and_strip_integrals() {
    let (design, loads) = fixture(0.025);
    let actual = design.distributed_lateral_loads(&loads, DC, RHO).unwrap();
    let (panels, strips) = design.layout(DC).unwrap();
    let local = loads.coupled.strip_velocity_mps.as_ref().unwrap();
    let (mut roll, mut yaw) = (0.0, 0.0);
    let (mut profile_roll, mut profile_yaw) = (0.0, 0.0);
    for (strip, v) in strips.iter().zip(local) {
        for &j in &strip.panel_indices {
            let p = panels[j];
            let y = (p.a[1] + p.b[1]) / 2.0 - design.cg_m[1];
            let span = p.b[1] - p.a[1];
            // These flat bound vortices point along y, so Fx=-rho Gamma w b
            // and Fz=rho Gamma u b. Independently integrate y Fz and -y Fx.
            roll += y * RHO * loads.coupled.gamma[j] * v[0] * span;
            yaw += y * RHO * loads.coupled.gamma[j] * v[2] * span;
        }
        let p = panels[strip.panel_indices[0]];
        let y = (p.a[1] + p.b[1]) / 2.0 - design.cg_m[1];
        let span = (p.b[1] - p.a[1]).abs();
        let speed = (v[0] * v[0] + v[2] * v[2]).sqrt();
        let factor = 0.5 * RHO * speed * strip.chord_m * span * design.cd0_wing;
        profile_roll -= y * factor * v[2];
        profile_yaw += y * factor * v[0];
    }
    assert!((actual.circulation.roll_nm - roll).abs() < 1e-9);
    assert!((actual.circulation.yaw_nm - yaw).abs() < 1e-9);
    assert!((actual.profile_drag.roll_nm - profile_roll).abs() < 1e-10);
    assert!((actual.profile_drag.yaw_nm - profile_yaw).abs() < 1e-10);
    assert!(actual.net().roll_nm < -1.0, "right-side upgust must roll left");
}

#[test]
fn uniform_flow_has_no_spurious_torque_and_mirrored_gust_reverses_it() {
    let (design, uniform) = fixture(0.0);
    let symmetric = design.distributed_lateral_loads(&uniform, DC, RHO).unwrap().net();
    assert!(symmetric.roll_nm.abs() < 1e-8);
    assert!(symmetric.yaw_nm.abs() < 1e-8);
    let (_, right) = fixture(0.025);
    let (_, left) = fixture(-0.025);
    let a = design.distributed_lateral_loads(&right, DC, RHO).unwrap().net();
    let b = design.distributed_lateral_loads(&left, DC, RHO).unwrap().net();
    assert!(a.roll_nm.abs() > 1.0);
    assert!((a.roll_nm + b.roll_nm).abs() < 1e-7);
    assert!((a.yaw_nm + b.yaw_nm).abs() < 1e-7);
}

#[test]
fn differential_propeller_thrust_uses_the_actual_lateral_arm_once() {
    let (design, mut loads) = fixture(0.0);
    let before = design.distributed_lateral_loads(&loads, DC, RHO).unwrap();
    loads.thrust_n[0] += 10.0;
    loads.thrust_n[1] -= 10.0;
    loads.coupled.thrust_n = loads.thrust_n;
    let after = design.distributed_lateral_loads(&loads, DC, RHO).unwrap();
    let expected = 10.0 * (design.disks[1].center_m[1] - design.disks[0].center_m[1]);
    assert!((after.propellers.yaw_nm - before.propellers.yaw_nm - expected).abs() < 1e-12);
    assert_eq!(before.propellers.roll_nm, after.propellers.roll_nm);
    assert_eq!(before.circulation, after.circulation);
    assert_eq!(before.profile_drag, after.profile_drag);
}

#[test]
fn incident_torque_and_reduced_controls_share_one_rate_and_attitude_update() {
    let model = LateralModel::wright_v1(RudderLinkage::Linked { gain_nm_per_rad: 66.0 });
    let original = LateralState { phi_rad: 0.1, p_rad_s: 0.02, psi_rad: -0.1, r_rad_s: -0.03 };
    let dt = 1.0 / 120.0;
    let incident = LateralMomentsNm { roll_nm: -250.0, yaw_nm: 34.0 };
    let mut plain = original;
    let base = model.step(&mut plain, 0.01, 13.86, RHO, dt).unwrap();
    let mut forced = original;
    let out = model.step_with_moments(&mut forced, 0.01, 13.86, RHO, dt, incident).unwrap();
    let dp = dt * incident.roll_nm / FLYER_ROLL_INERTIA_KG_M2;
    let dr = dt * incident.yaw_nm / REDUCED_YAW_INERTIA_KG_M2;
    assert!((forced.p_rad_s - plain.p_rad_s - dp).abs() < 1e-15);
    assert!((forced.r_rad_s - plain.r_rad_s - dr).abs() < 1e-15);
    assert!((forced.phi_rad - plain.phi_rad - dt * dp).abs() < 1e-15);
    assert!((forced.psi_rad - plain.psi_rad - dt * dr).abs() < 1e-15);
    assert_eq!(out.reduced.yaw, base.yaw, "do not relabel gust yaw as rudder/drag calibration");
    assert_eq!(out.net_yaw_moment_nm, base.yaw.net() + incident.yaw_nm);
    assert!((out.net_roll_moment_nm / FLYER_ROLL_INERTIA_KG_M2 * dt
        - (forced.p_rad_s - original.p_rad_s)).abs() < 1e-15);
    let mut zero = original;
    let unchanged = model.step_with_moments(&mut zero, 0.01, 13.86, RHO, dt,
        LateralMomentsNm::default()).unwrap();
    assert_eq!(unchanged.reduced, base);
    assert_eq!(zero, plain);
}

#[test]
fn refused_incident_torque_is_atomic_including_yaw_failure_after_roll_trial() {
    let model = LateralModel::wright_v1(RudderLinkage::Decoupled);
    let original = LateralState { phi_rad: 0.1, p_rad_s: 0.02, psi_rad: 0.2, r_rad_s: 0.01 };
    for incident in [
        LateralMomentsNm { roll_nm: f64::NAN, yaw_nm: 0.0 },
        LateralMomentsNm { roll_nm: 1e308, yaw_nm: 0.0 },
        LateralMomentsNm { roll_nm: 10.0, yaw_nm: 1e308 },
    ] {
        let mut state = original;
        assert_eq!(model.step_with_moments(&mut state, 0.0, 13.86, RHO, 1.0 / 120.0,
            incident).unwrap_err().code, "lateral-envelope-exceeded");
        assert_eq!(state, original);
    }
}

#[test]
fn inconsistent_accepted_fields_refuse_without_uniform_flow_fallback() {
    let (design, loads) = fixture(0.0);
    let mut missing = loads.clone();
    missing.coupled.strip_velocity_mps = None;
    let mut short = loads.clone();
    short.coupled.gamma.pop();
    let mut nonfinite = loads.clone();
    nonfinite.coupled.strip_velocity_mps.as_mut().unwrap()[0][2] = f64::NAN;
    let mut short_local = loads.clone();
    short_local.coupled.strip_velocity_mps.as_mut().unwrap().pop();
    for bad in [missing, short, nonfinite, short_local] {
        assert_eq!(design.distributed_lateral_loads(&bad, DC, RHO).unwrap_err().code,
            "distributed-loads-invalid");
    }
    assert!(design.distributed_lateral_loads(&loads, DC, f64::NAN).is_err());
}
