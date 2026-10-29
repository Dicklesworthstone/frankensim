//! G0/G1/G3 evolving native blade/wake/gear/airframe momentum.
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState};
use fs_flyer::spine::rotors::propulsion::dynamic::{DynamicPropulsion, DynamicState};
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::{Mat3, Quat, So3, Vec3};

const RHO: f64 = 1.225;
const ZERO: [f64; 3] = [0.0; 3];

fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let x = v(a).scale(1.0 / v(a).norm());
    [x.x, x.y, x.z]
}

fn axes() -> [[f64; 3]; 2] {
    [unit([1.0, 0.08, -0.02]), unit([1.0, -0.03, 0.06])]
}

fn assembly(mass: f64, reverse_coordinates: bool) -> GearedAssembly<2> {
    let sign = if reverse_coordinates { -1.0 } else { 1.0 };
    GearedAssembly::new(
        mass,
        Mat3 {
            m: [1800.0, 10.0, 5.0, 10.0, 420.0, 2.0, 5.0, 2.0, 1850.0],
        },
        core::array::from_fn(|k| GearedRotor {
            axis_body: axes()[k].map(|x| sign * x),
            ratio: sign * [1.0, -1.0][k],
            axial_inertia_kgm2: [1.5, 2.5][k],
        }),
    )
    .unwrap()
}

fn mounts(reverse_coordinates: bool) -> [PropellerMount; 2] {
    core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        thrust_direction_body: axes()[k],
        handedness: if (k == 0) != reverse_coordinates {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}

fn state(assembly: &GearedAssembly<2>) -> DynamicState {
    DynamicState {
        motion: assembly
            .state_from_speed(
                SixDofState {
                    pos_m: [1.0, -2.0, -10.0],
                    vel_mps: [12.0, 0.1, 0.8],
                    quat: [1.0, 0.0, 0.0, 0.0],
                    omega_body: [0.04, -0.02, 0.03],
                },
                0.2,
                35.0,
            )
            .unwrap(),
        induced_mps: [0.4, 0.8],
    }
}

fn inputs(time: f64) -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: RHO,
        canard_rad: 0.03,
        motor_nm: 100.0 + 30.0 * time,
        gravity_world_mps2: [0.0, 0.0, 9.80665],
    }
}

fn momentum(assembly: &GearedAssembly<2>, state: &GearedState, origin: [f64; 3]) -> Vec3 {
    let q = state.rigid.quat;
    let rotation = So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap();
    let mut spin = assembly
        .locked_inertia_kgm2()
        .apply(v(state.rigid.omega_body));
    for rotor in assembly.rotors() {
        spin = spin
            + v(rotor.axis_body)
                .scale(rotor.axial_inertia_kgm2 * rotor.ratio * assembly.speed(state).unwrap());
    }
    rotation.rotate(spin).unwrap()
        + (v(state.rigid.pos_m) - v(origin)).cross(v(state.rigid.vel_mps).scale(assembly.mass_kg()))
}

#[test]
fn canted_transient_closes_aircraft_impulse_and_each_physical_wake_balance() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let propulsion = DynamicPropulsion::new(&base, RHO).unwrap();
    let mut state = state(&assembly);
    state.induced_mps = [0.0; 2]; // startup of both wakes, not a developed equilibrium
    let initial = state;
    let origin = [-5.0, 3.0, 2.0];
    let initial_angular = momentum(&assembly, &state.motion, origin);
    let initial_linear = v(state.motion.rigid.vel_mps).scale(assembly.mass_kg());
    let mut angular_impulse = Vec3::default();
    let mut linear_impulse = Vec3::default();
    let mut wake_impulse = [0.0; 2];
    let mut wake_work = [0.0; 2];
    let mut maximum_rigid_defect: f64 = 0.0;
    let mut maximum_wake_defect: f64 = 0.0;
    for tick in 0..32 {
        let step = propulsion
            .step(
                &state,
                f64::from(tick) * 0.004,
                0.004,
                origin,
                |time, _| Ok(inputs(time)),
                |_, point| Ok([0.02 * point[1], 0.0, 0.0]),
            )
            .unwrap();
        assert_eq!(step.initial.aerodynamics.induced_mps, state.induced_mps);
        assert!(step.midpoint.aerodynamics.wing_residual < 1e-10);
        assert!(v(step.motion.angular_momentum_residual).norm() < 2e-9);
        assert!(v(step.motion.linear_momentum_residual).norm() < 2e-10);
        maximum_rigid_defect =
            maximum_rigid_defect.max(v(step.motion.angular_momentum_residual).norm());
        angular_impulse = angular_impulse + v(step.motion.external_angular_impulse_world);
        linear_impulse = linear_impulse + v(step.motion.external_linear_impulse_world);
        for k in 0..2 {
            let flow = step.inflow[k];
            assert!(flow.momentum_residual_ns.abs() < 2e-13);
            assert!(flow.energy_residual_j.abs() < 2e-12);
            maximum_wake_defect = maximum_wake_defect.max(flow.momentum_residual_ns.abs());
            wake_impulse[k] += flow.thrust_impulse_ns - flow.transport_impulse_ns;
            wake_work[k] += flow.induced_work_j - flow.transport_work_j;
        }
        state = step.state;
    }
    let angular_defect =
        (momentum(&assembly, &state.motion, origin) - initial_angular - angular_impulse).norm();
    let linear_defect =
        (v(state.motion.rigid.vel_mps).scale(assembly.mass_kg()) - initial_linear - linear_impulse)
            .norm();
    assert!(angular_defect < 2e-8);
    assert!(linear_defect < 2e-9);
    let mass = propulsion.inflow_model().apparent_mass_kg();
    for k in 0..2 {
        assert!(
            (mass * (state.induced_mps[k] - initial.induced_mps[k]) - wake_impulse[k]).abs()
                < 2e-12
        );
        let energy = 0.5 * mass * (state.induced_mps[k].powi(2) - initial.induced_mps[k].powi(2));
        assert!((energy - wake_work[k]).abs() < 2e-11);
        assert!(state.induced_mps[k] > initial.induced_mps[k] + 0.1);
    }
    assert!((assembly.speed(&state.motion).unwrap() - 35.0).abs() > 0.01);
    println!(
        "32-tick dynamic propulsion: angular defect={angular_defect:e}, linear defect={linear_defect:e}, max rigid angular defect={maximum_rigid_defect:e}, max wake impulse defect={maximum_wake_defect:e}, induction={:?}",
        state.induced_mps
    );
}

#[test]
fn reversing_mechanical_axis_coordinates_preserves_the_entire_dynamic_trajectory() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, false);
    let b = assembly(design.gross_mass_kg, true);
    let base_a = GearedPropulsion::new_spatial(&design, &a, mounts(false)).unwrap();
    let base_b = GearedPropulsion::new_spatial(&design, &b, mounts(true)).unwrap();
    let dynamic_a = DynamicPropulsion::new(&base_a, RHO).unwrap();
    let dynamic_b = DynamicPropulsion::new(&base_b, RHO).unwrap();
    let mut state_a = state(&a);
    let mut state_b = state(&b);
    for tick in 0..8 {
        let run = |model: &DynamicPropulsion<'_, '_, 2>, state: &DynamicState| {
            model
                .step(
                    state,
                    f64::from(tick) * 0.004,
                    0.004,
                    ZERO,
                    |t, _| Ok(inputs(t)),
                    |_, _| Ok(ZERO),
                )
                .unwrap()
        };
        let next_a = run(&dynamic_a, &state_a);
        let next_b = run(&dynamic_b, &state_b);
        assert_eq!(next_a.state, next_b.state);
        assert_eq!(next_a.inflow, next_b.inflow);
        assert_eq!(next_a.midpoint.aerodynamics, next_b.midpoint.aerodynamics);
        assert_eq!(
            next_a.motion.external_angular_impulse_world,
            next_b.motion.external_angular_impulse_world
        );
        for k in 0..2 {
            assert_eq!(
                next_a.midpoint.loads.aerodynamic_nm[k],
                -next_b.midpoint.loads.aerodynamic_nm[k]
            );
        }
        state_a = next_a.state;
        state_b = next_b.state;
    }
}

#[test]
fn real_coupled_shaft_wake_and_airframe_have_second_order_timestep_refinement() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let propulsion = DynamicPropulsion::new(&base, RHO).unwrap();
    let mut endpoints = Vec::new();
    for count in [8, 16, 32, 64] {
        let dt = 0.12 / f64::from(count);
        let mut state = state(&assembly);
        for tick in 0..count {
            state = propulsion
                .step(
                    &state,
                    f64::from(tick) * dt,
                    dt,
                    ZERO,
                    |t, _| Ok(inputs(t)),
                    |_, _| Ok(ZERO),
                )
                .unwrap()
                .state;
        }
        endpoints.push([
            state.motion.rigid.vel_mps[0],
            state.motion.rigid.vel_mps[1],
            state.motion.rigid.vel_mps[2],
            state.motion.rigid.omega_body[0],
            state.motion.rigid.omega_body[1],
            state.motion.rigid.omega_body[2],
            assembly.speed(&state.motion).unwrap(),
            state.induced_mps[0],
            state.induced_mps[1],
        ]);
    }
    let differences: Vec<_> = endpoints
        .windows(2)
        .map(|pair| {
            pair[0]
                .iter()
                .zip(pair[1])
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    println!("dynamic coupled timestep differences: {differences:?}");
    for pair in differences.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "observed refinement ratio {ratio}"
        );
    }
}

#[test]
fn joint_midpoint_receives_predicted_wakes_and_late_density_refusal_is_atomic() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let propulsion = DynamicPropulsion::new(&base, RHO).unwrap();
    let before = state(&assembly);
    let mut samples = Vec::new();
    let error = propulsion
        .step(
            &before,
            0.0,
            0.004,
            ZERO,
            |time, stage| {
                samples.push((time, *stage));
                let mut values = inputs(time);
                if time > 0.0 {
                    values.density_kg_m3 *= 1.001;
                }
                Ok(values)
            },
            |_, _| Ok(ZERO),
        )
        .unwrap_err();
    assert_eq!(error.code, "dynamic-propulsion-invalid");
    assert_eq!(samples.len(), 2);
    assert_eq!(samples[0], (0.0, before));
    assert_eq!(samples[1].0, 0.002);
    assert_ne!(samples[1].1.induced_mps, before.induced_mps);
    assert_ne!(samples[1].1.motion.momentum, before.motion.momentum);
    assert_eq!(before, state(&assembly));
}
