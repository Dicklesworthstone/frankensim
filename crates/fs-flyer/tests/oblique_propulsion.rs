//! G0/G1/G3: real oblique blades, canted shafts and aircraft momentum.
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
fn a(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn unit(x: [f64; 3]) -> [f64; 3] {
    a(v(x).scale(1.0 / v(x).norm()))
}
fn axes() -> [[f64; 3]; 2] {
    [unit([1.0, 0.08, -0.02]), unit([1.0, -0.03, 0.06])]
}

fn assembly(mass: f64, reverse: bool) -> GearedAssembly<2> {
    let sign = if reverse { -1.0 } else { 1.0 };
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

fn mounts(reverse: bool) -> [PropellerMount; 2] {
    core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        thrust_direction_body: axes()[k],
        handedness: if (k == 0) != reverse {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}

fn initial(assembly: &GearedAssembly<2>) -> DynamicState {
    DynamicState {
        motion: assembly
            .state_from_speed(
                SixDofState {
                    pos_m: [1.0, -2.0, -10.0],
                    vel_mps: [12.0, 2.0, 0.8],
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

fn inputs(t: f64) -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: RHO,
        canard_rad: 0.03,
        motor_nm: 100.0 + 30.0 * t,
        gravity_world_mps2: [0.0, 0.0, 9.80665],
    }
}

fn angular_momentum(assembly: &GearedAssembly<2>, state: &GearedState, origin: [f64; 3]) -> Vec3 {
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
fn full_crossflow_force_and_transverse_hub_moment_reach_the_original_carrier_once() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let axial = DynamicPropulsion::new(&base, RHO).unwrap();
    let oblique = DynamicPropulsion::new_oblique(&base, RHO, 64).unwrap();
    let state = initial(&assembly);
    let old = axial.evaluate(&state, inputs(0.0), &|_| Ok(ZERO)).unwrap();
    let new = oblique
        .evaluate(&state, inputs(0.0), &|_| Ok(ZERO))
        .unwrap();
    assert_eq!(old.aerodynamics.gamma, new.aerodynamics.gamma);
    assert_eq!(
        old.aerodynamics.strip_velocity_mps,
        new.aerodynamics.strip_velocity_mps
    );
    let loads = new.aerodynamics.oblique.unwrap();
    let mut delta_force = Vec3::default();
    let mut delta_moment = Vec3::default();
    let mut transverse = 0.0;
    for k in 0..2 {
        let difference =
            v(loads[k].force_n) - v(axes()[k]).scale(old.aerodynamics.propellers[k].thrust_n);
        delta_force = delta_force + difference;
        delta_moment = delta_moment
            + (v(design.disks[k].center_m) - v(design.cg_m)).cross(difference)
            + v(loads[k].hub_moment_nm);
        let offaxis = v(loads[k].force_n) - v(axes()[k]).scale(loads[k].thrust_n);
        transverse += offaxis.norm() + v(loads[k].hub_moment_nm).norm();
        assert_eq!(
            new.loads.aerodynamic_nm[k],
            [-1.0, 1.0][k] * loads[k].torque_nm
        );
        assert!(loads[k].power_residual_w.abs() < 1e-9);
        assert!(loads[k].profile_dissipation_w > 0.0);
        assert!(
            (loads[k].kinematic_power_w + loads[k].induced_work_w + loads[k].profile_dissipation_w)
                .abs()
                < 1e-9
        );
    }
    assert!(transverse > 0.1);
    assert!(
        (v(new.aerodynamic_force_body_n) - v(old.aerodynamic_force_body_n) - delta_force).norm()
            < 2e-10
    );
    assert!(
        (v(new.loads.carrier.moment_nm) - v(old.loads.carrier.moment_nm) - delta_moment).norm()
            < 2e-10
    );
    // Axial propeller torque was excluded from delta_moment: adding it to
    // both the carrier and shaft ports would make this reconstruction fail.
}

#[test]
fn actual_oblique_body_shaft_wake_trajectory_closes_world_linear_and_angular_impulse() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let model = DynamicPropulsion::new_oblique(&base, RHO, 64).unwrap();
    let mut state = initial(&assembly);
    let origin = [-5.0, 3.0, 2.0];
    let before = state;
    let p0 = v(state.motion.rigid.vel_mps).scale(assembly.mass_kg());
    let l0 = angular_momentum(&assembly, &state.motion, origin);
    let mut dp = Vec3::default();
    let mut dl = Vec3::default();
    let mut max_power: f64 = 0.0;
    for tick in 0..24 {
        let result = model
            .step(
                &state,
                tick as f64 * 0.004,
                0.004,
                origin,
                |t, _| Ok(inputs(t)),
                |_, point| Ok([0.02 * point[1], 0.0, 0.0]),
            )
            .unwrap();
        dp = dp + v(result.motion.external_linear_impulse_world);
        dl = dl + v(result.motion.external_angular_impulse_world);
        for load in result.midpoint.aerodynamics.oblique.unwrap() {
            max_power = max_power.max(load.power_residual_w.abs());
        }
        for wake in result.inflow {
            assert!(wake.momentum_residual_ns.abs() < 3e-13);
            assert!(wake.energy_residual_j.abs() < 3e-12);
        }
        state = result.state;
    }
    let p_defect = (v(state.motion.rigid.vel_mps).scale(assembly.mass_kg()) - p0 - dp).norm();
    let l_defect = (angular_momentum(&assembly, &state.motion, origin) - l0 - dl).norm();
    println!(
        "24-tick oblique propulsion: linear defect={p_defect:e}, angular defect={l_defect:e}, max section power defect={max_power:e}"
    );
    assert!(p_defect < 2e-9);
    assert!(l_defect < 2e-8);
    assert!(max_power < 2e-9);
    assert!(
        (assembly.speed(&state.motion).unwrap() - assembly.speed(&before.motion).unwrap()).abs()
            > 0.01
    );
    assert!((state.induced_mps[0] - before.induced_mps[0]).abs() > 0.01);
}

#[test]
fn mechanical_axis_sign_conventions_preserve_the_full_oblique_trajectory() {
    let design = wright_openloop_v1();
    let aa = assembly(design.gross_mass_kg, false);
    let ab = assembly(design.gross_mass_kg, true);
    let ba = GearedPropulsion::new_spatial(&design, &aa, mounts(false)).unwrap();
    let bb = GearedPropulsion::new_spatial(&design, &ab, mounts(true)).unwrap();
    let ma = DynamicPropulsion::new_oblique(&ba, RHO, 32).unwrap();
    let mb = DynamicPropulsion::new_oblique(&bb, RHO, 32).unwrap();
    let mut sa = initial(&aa);
    let mut sb = initial(&ab);
    for tick in 0..8 {
        let run = |model: &DynamicPropulsion<'_, '_, 2>, state: &DynamicState| {
            model
                .step(
                    state,
                    tick as f64 * 0.004,
                    0.004,
                    ZERO,
                    |t, _| Ok(inputs(t)),
                    |_, _| Ok(ZERO),
                )
                .unwrap()
        };
        let a = run(&ma, &sa);
        let b = run(&mb, &sb);
        assert_eq!(a.state, b.state);
        assert_eq!(a.inflow, b.inflow);
        assert_eq!(a.midpoint.aerodynamics, b.midpoint.aerodynamics);
        assert_eq!(
            a.motion.external_angular_impulse_world,
            b.motion.external_angular_impulse_world
        );
        for k in 0..2 {
            assert_eq!(
                a.midpoint.loads.aerodynamic_nm[k],
                -b.midpoint.loads.aerodynamic_nm[k]
            );
        }
        sa = a.state;
        sb = b.state;
    }
}
