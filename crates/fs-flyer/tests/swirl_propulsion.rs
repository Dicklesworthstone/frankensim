//! Real wing/swirl loads, signed shaft ports and complete body/rotor momentum.
use fs_airscrew::swirl::Controls;
use fs_airscrew::{BladeStation, Rotor};
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::propcoupling::{CANDIDATE_B, PropCouplingSolverSpec};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState};
use fs_flyer::spine::rotors::propulsion::swirl::SwirlPropulsion;
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};

const ZERO: [f64; 3] = [0.0; 3];
const AXES: [[f64; 3]; 2] = [[1.0, 0.0, 0.0]; 2];
fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn design() -> OpenLoopDesign {
    let mut design = wright_openloop_v1();
    design.coupling = PropCouplingSolverSpec {
        tol: 1e-10,
        cap: 64,
        ..CANDIDATE_B
    };
    design
}
fn assembly(mass: f64, axes: [[f64; 3]; 2], ratios: [f64; 2]) -> GearedAssembly<2> {
    GearedAssembly::new(
        mass,
        Mat3 {
            m: [1000.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0, 1200.0],
        },
        core::array::from_fn(|k| GearedRotor {
            axis_body: axes[k],
            ratio: ratios[k],
            axial_inertia_kgm2: [1.5, 2.5][k],
        }),
    )
    .unwrap()
}
fn mounts(axes: [[f64; 3]; 2], ratios: [f64; 2]) -> [PropellerMount; 2] {
    core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        thrust_direction_body: axes[k],
        handedness: if ratios[k] > 0.0 {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}
fn state(assembly: &GearedAssembly<2>, speed: f64, alpha: f64, omega: f64) -> GearedState {
    assembly
        .state_from_speed(
            SixDofState {
                pos_m: [2.0, -3.0, -10.0],
                vel_mps: [speed * alpha.cos(), 0.0, speed * alpha.sin()],
                quat: [1.0, 0.0, 0.0, 0.0],
                omega_body: ZERO,
            },
            0.0,
            omega,
        )
        .unwrap()
}
fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        motor_nm: 80.0,
        gravity_world_mps2: ZERO,
    }
}
fn orientation(state: &GearedState) -> So3 {
    let q = state.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn angular_momentum(a: &GearedAssembly<2>, state: &GearedState, origin: [f64; 3]) -> Vec3 {
    let mut spin = a.locked_inertia_kgm2().apply(vec(state.rigid.omega_body));
    let relative = a.speed(state).unwrap();
    for r in a.rotors() {
        spin = spin + vec(r.axis_body).scale(r.axial_inertia_kgm2 * r.ratio * relative);
    }
    orientation(state).rotate(spin).unwrap()
        + (vec(state.rigid.pos_m) - vec(origin)).cross(vec(state.rigid.vel_mps).scale(a.mass_kg()))
}

#[test]
fn actual_swirl_changes_two_way_wing_loads_and_reaches_signed_shaft_ports() {
    let d = design();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new(&d, &a, mounts(AXES, [1.0, -1.0])).unwrap();
    let model = SwirlPropulsion::new(&base, Controls::default());
    let state = state(&a, 12.0, 0.08, 35.0);
    let plain = base
        .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let result = model
        .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let (panels, _) = d.layout(inputs().canard_rad).unwrap();
    let coupled = &result.coupled.coupled;
    assert_ne!(coupled.spec_digest, plain.coupled.spec_digest);
    assert!((coupled.thrust_n[0] - plain.coupled.thrust_n[0]).abs() > 1e-4);
    assert!((coupled.torque_nm[0] - plain.coupled.torque_nm[0]).abs() > 1e-4);
    assert!(
        coupled
            .gamma
            .iter()
            .zip(&plain.coupled.gamma)
            .any(|(a, b)| (a - b).abs() > 1e-7)
    );
    for k in 0..2 {
        let rotor = &result.coupled.rotors[k];
        let sign = [1.0, -1.0][k];
        assert_eq!(result.loads.aerodynamic_nm[k], -sign * rotor.torque_nm);
        assert_eq!(coupled.thrust_n[k], rotor.thrust_n);
        assert_eq!(coupled.torque_nm[k], rotor.torque_nm);
        assert_eq!(coupled.power[k], rotor.power);
        assert_eq!(
            rotor.power.shaft_w,
            -result.loads.aerodynamic_nm[k] * result.propeller_absolute_speed_rad_s[k]
        );
        assert!(rotor.power.swirl_w > 0.0);
        let wing = fs_wing::induced_velocity_free(
            d.disks[k].center_m,
            &panels,
            &coupled.gamma,
            state.rigid.vel_mps,
        );
        let axial = state.rigid.vel_mps[0] + wing[0]; // separated disks have no mutual overlap
        let seeds: Vec<_> = rotor.stations.iter().map(|s| s.induction).collect();
        let independent = fs_airscrew::swirl::bemt_solve(
            &d.rotor,
            1.225,
            axial,
            35.0,
            Controls::default(),
            Some(&seeds),
        )
        .unwrap();
        near(independent.thrust_n, rotor.thrust_n, 1e-10);
        near(independent.torque_nm, rotor.torque_nm, 1e-10);
        let vi = coupled.w_slip[k];
        let area = std::f64::consts::PI * d.rotor.radius_m * d.rotor.radius_m;
        near(2.0 * 1.225 * area * (axial + vi) * vi, rotor.thrust_n, 1e-9);
        for station in &rotor.stations {
            near(
                station.thrust_n,
                station.axial_momentum_flux_n,
                3e-10 * station.thrust_n.abs(),
            );
            near(
                station.torque_nm,
                station.angular_momentum_flux_nm,
                3e-10 * station.torque_nm.abs(),
            );
        }
    }
    let continued = model
        .evaluate(
            &state,
            inputs(),
            &|_| Ok(ZERO),
            Some(&result.coupled.warm_start()),
        )
        .unwrap();
    near(
        continued.coupled.coupled.thrust_n[0],
        coupled.thrust_n[0],
        1e-7,
    );
    assert_eq!(continued.coupled.coupled.spec_digest, coupled.spec_digest);
}

#[test]
fn actual_midpoint_swirl_loads_conserve_global_momentum_without_double_shaft_recoil() {
    let d = design();
    // Co-rotation gives a nonzero net shaft moment. It must stay out of the
    // carrier-only wrench while remaining external to the whole assembly.
    let a = assembly(d.gross_mass_kg, AXES, [1.0, 1.0]);
    let base = GearedPropulsion::new(&d, &a, mounts(AXES, [1.0, 1.0])).unwrap();
    let model = SwirlPropulsion::new(&base, Controls::default());
    let start = state(&a, 12.0, 0.08, 35.0);
    let origin = start.rigid.pos_m;
    let control = PropulsionInputs {
        motor_nm: 0.0,
        ..inputs()
    };
    let result = model
        .step(
            &start,
            0.0,
            0.002,
            origin,
            None,
            |_, _| Ok(control),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert!(result.initial.loads.carrier.moment_nm[0].abs() < 1e-9);
    assert!(
        result
            .initial
            .loads
            .aerodynamic_nm
            .iter()
            .sum::<f64>()
            .abs()
            > 1.0
    );
    assert!(result.motion.external_angular_impulse_world[0].abs() > 1e-4);
    assert_ne!(
        result.initial.coupled.rotors[0].torque_nm,
        result.midpoint.coupled.rotors[0].torque_nm
    );
    let delta_l =
        angular_momentum(&a, &result.motion.state, origin) - angular_momentum(&a, &start, origin);
    assert!((delta_l - vec(result.motion.external_angular_impulse_world)).norm() < 2e-10);
    let delta_p =
        (vec(result.motion.state.rigid.vel_mps) - vec(start.rigid.vel_mps)).scale(a.mass_kg());
    assert!((delta_p - vec(result.motion.external_linear_impulse_world)).norm() < 2e-10);
    assert_eq!(result.motion.motor_work_j, 0.0);
}

/// An attached-flow windmilling rotor with prescribed w=-2 m/s, u=-0.15 m/s
/// at V=10 m/s and Omega=60 rad/s, constructed from both momentum laws.
fn windmill_rotor() -> Rotor {
    let stations = [0.3, 0.5, 0.7, 0.9]
        .into_iter()
        .map(|r| {
            let va = 8.0_f64;
            let vt = 60.0 * r + 0.15;
            let phi = va.atan2(vt);
            let ratio = (phi.cos() + 0.075 * phi.sin()) / (0.075 * phi.cos() - phi.sin());
            let b = std::f64::consts::TAU;
            let alpha = 0.024 * ratio / (b + (b * b - 4.0 * 0.8 * 0.012 * ratio * ratio).sqrt());
            let loss = |distance: f64| {
                2.0 / std::f64::consts::PI * (-1.5 * distance / (r * phi.sin())).exp().acos()
            };
            let f = (loss(1.0 - r) * loss(r - 0.24)).clamp(1e-3, 1.0);
            let target = 4.0 * std::f64::consts::PI * 1.225 * r * f * va * (-2.0);
            let unit_force = 1.5
                * 1.225
                * (va * va + vt * vt)
                * (b * alpha * phi.cos() - (0.012 + 0.8 * alpha * alpha) * phi.sin());
            BladeStation {
                r_over_r: r,
                chord_m: target / unit_force,
                beta_rad: phi + alpha,
            }
        })
        .collect();
    Rotor {
        radius_m: 1.0,
        n_blades: 3,
        camber_ratio: 0.0,
        stations,
    }
}

#[test]
fn signed_windmilling_transfers_air_momentum_into_real_geared_spin() {
    let mut d = design();
    d.camber_ratio = 0.0;
    d.rotor = windmill_rotor();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new(&d, &a, mounts(AXES, [1.0, -1.0])).unwrap();
    let model = SwirlPropulsion::new(&base, Controls::default());
    let start = state(&a, 10.0, 0.0, 60.0);
    let control = PropulsionInputs {
        canard_rad: 0.0,
        motor_nm: 0.0,
        ..inputs()
    };
    let result = model
        .step(
            &start,
            0.0,
            0.001,
            start.rigid.pos_m,
            None,
            |_, _| Ok(control),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    for stage in [&result.initial, &result.midpoint] {
        for k in 0..2 {
            let rotor = &stage.coupled.rotors[k];
            assert!(rotor.thrust_n < 0.0 && rotor.torque_nm < 0.0);
            assert!(rotor.power.shaft_w < 0.0 && rotor.power.swirl_w > 0.0);
            assert!(stage.loads.aerodynamic_nm[k] * stage.propeller_absolute_speed_rad_s[k] > 0.0);
        }
    }
    assert!(result.motion.speed_rad_s > a.speed(&start).unwrap());
    assert!(result.motion.state.rigid.vel_mps[0] < start.rigid.vel_mps[0]);
    assert!(vec(result.motion.angular_momentum_residual).norm() < 1e-10);
}

#[test]
fn spatial_swirl_preserves_body_loads_under_a_world_frame_rotation() {
    let unit = |x: [f64; 3]| array(vec(x).scale(1.0 / vec(x).norm()));
    let axes = [unit([1.0, 0.04, 0.06]), unit([1.0, -0.03, 0.02])];
    let d = design();
    let a = assembly(d.gross_mass_kg, axes, [1.0, -1.0]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(axes, [1.0, -1.0])).unwrap();
    let model = SwirlPropulsion::new(&base, Controls::default());
    let start = state(&a, 12.0, 0.08, 35.0);
    let first = model
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    let mut moved = start;
    moved.rigid.pos_m = array(rotation.rotate(vec(start.rigid.pos_m)).unwrap());
    moved.rigid.vel_mps = array(rotation.rotate(vec(start.rigid.vel_mps)).unwrap());
    let q = rotation.as_quat();
    moved.rigid.quat = [q.w, q.x, q.y, q.z];
    let second = model
        .evaluate(&moved, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert!(
        (vec(first.aerodynamic_force_body_n) - vec(second.aerodynamic_force_body_n)).norm() < 1e-8
    );
    assert!(
        (rotation.rotate(vec(first.loads.carrier.force_n)).unwrap()
            - vec(second.loads.carrier.force_n))
        .norm()
            < 1e-8
    );
    for k in 0..2 {
        near(
            first.loads.aerodynamic_nm[k],
            second.loads.aerodynamic_nm[k],
            1e-9,
        );
    }
}

#[test]
fn invalid_swirl_or_late_stage_refuses_without_publishing_partial_motion() {
    let d = design();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new(&d, &a, mounts(AXES, [1.0, -1.0])).unwrap();
    let start = state(&a, 12.0, 0.08, 35.0);
    let original = start;
    let bounded = SwirlPropulsion::new(
        &base,
        Controls {
            maximum_station_evaluations: 3,
            ..Controls::default()
        },
    );
    assert_eq!(
        bounded
            .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
            .unwrap_err()
            .code,
        "swirl-nonconvergence"
    );
    let model = SwirlPropulsion::new(&base, Controls::default());
    let mut calls = 0;
    let result = model.step(
        &start,
        0.0,
        0.001,
        ZERO,
        None,
        |_, _| {
            calls += 1;
            Ok(PropulsionInputs {
                motor_nm: if calls == 2 { f64::NAN } else { 0.0 },
                ..inputs()
            })
        },
        |_, _| Ok(ZERO),
    );
    assert_eq!(result.unwrap_err().code, "propulsion-state-invalid");
    assert_eq!(calls, 2);
    assert_eq!(start, original);
}
