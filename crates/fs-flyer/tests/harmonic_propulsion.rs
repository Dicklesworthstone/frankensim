//! G0/G1/G3: asymmetric wake memory in the actual aircraft/shaft trajectory.
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState};
use fs_flyer::spine::rotors::propulsion::dynamic::first_harmonic::{
    HarmonicPropulsion, HarmonicState,
};
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
fn initial(assembly: &GearedAssembly<2>) -> HarmonicState {
    HarmonicState {
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
        induction: [
            HarmonicInduction {
                mean_mps: 0.4,
                first_mps: [0.2, -0.15],
            },
            HarmonicInduction {
                mean_mps: 0.8,
                first_mps: [-0.25, 0.12],
            },
        ],
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
fn asymmetric_wake_reaches_both_real_blades_and_wing_strip_forces() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let model = HarmonicPropulsion::new(&base, RHO, 64).unwrap();
    let old_model = DynamicPropulsion::new_oblique(&base, RHO, 64).unwrap();
    let mut state = initial(&assembly);
    let asymmetric = state;
    for u in &mut state.induction {
        u.first_mps = [0.0; 2];
    }
    let old = old_model
        .evaluate(
            &DynamicState {
                motion: state.motion,
                induced_mps: state.induction.map(|u| u.mean_mps),
            },
            inputs(0.0),
            &|_| Ok(ZERO),
        )
        .unwrap();
    let zero = model.evaluate(&state, inputs(0.0), &|_| Ok(ZERO)).unwrap();
    assert_eq!(old.loads, zero.loads);
    assert_eq!(
        old.aerodynamics.oblique.unwrap(),
        zero.aerodynamics.propellers
    );
    assert_eq!(old.aerodynamics.gamma, zero.aerodynamics.gamma);
    let new = model
        .evaluate(&asymmetric, inputs(0.0), &|_| Ok(ZERO))
        .unwrap();
    assert_ne!(zero.aerodynamics.gamma, new.aerodynamics.gamma);
    let wash_delta: f64 = new
        .aerodynamics
        .strip_velocity_mps
        .iter()
        .zip(&zero.aerodynamics.strip_velocity_mps)
        .map(|(new, old)| (v(*new) - v(*old)).norm())
        .sum();
    assert!(wash_delta > 0.01);
    let frames = model.disk_frames();
    for k in 0..2 {
        let load = new.aerodynamics.propellers[k];
        let q = new.aerodynamics.wake_force_n[k];
        let expected = v(frames[k].thrust_axis())
            .cross(v(load.hub_moment_nm))
            .scale(1.0 / design.rotor.radius_m);
        assert!((q[1] - expected.dot(v(frames[k].radial_reference()))).abs() < 1e-12);
        assert!((q[2] - expected.dot(v(frames[k].azimuth_reference()))).abs() < 1e-12);
        assert!((v(q).dot(v(asymmetric.induction[k].modes())) - load.induced_work_w).abs() < 1e-10);
        assert!(load.power_residual_w.abs() < 2e-9);
        assert_eq!(new.loads.aerodynamic_nm[k], [-1.0, 1.0][k] * load.torque_nm);
    }
}

#[test]
fn evolving_harmonics_preserve_actual_aircraft_and_shaft_momentum_accounts() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let model = HarmonicPropulsion::new(&base, RHO, 64).unwrap();
    let before = initial(&assembly);
    let mut state = before;
    let origin = [-5.0, 3.0, 2.0];
    let p0 = v(state.motion.rigid.vel_mps).scale(assembly.mass_kg());
    let l0 = angular_momentum(&assembly, &state.motion, origin);
    let mut dp = Vec3::default();
    let mut dl = Vec3::default();
    let mut max_modal: f64 = 0.0;
    let mut max_work: f64 = 0.0;
    for tick in 0..24 {
        let step = model
            .step(
                &state,
                tick as f64 * 0.004,
                0.004,
                origin,
                |t, _| Ok(inputs(t)),
                |_, p| Ok([0.02 * p[1], 0.0, 0.0]),
            )
            .unwrap();
        dp = dp + v(step.motion.external_linear_impulse_world);
        dl = dl + v(step.motion.external_angular_impulse_world);
        for k in 0..2 {
            let wake = step.inflow[k];
            max_modal = max_modal.max(v(wake.momentum_residual_ns).norm());
            max_work = max_work.max(wake.energy_residual_j.abs());
            assert!(wake.frame_work_j.abs() < 1e-15);
            assert!(wake.transport_work_j >= 0.0);
            assert!(
                step.midpoint.aerodynamics.propellers[k]
                    .power_residual_w
                    .abs()
                    < 2e-9
            );
            let expected = wake.induced_work_j
                - 0.004 * step.midpoint.aerodynamics.propellers[k].induced_work_w;
            assert_eq!(expected, step.coupling_work_residual_j[k]);
        }
        state = step.state;
    }
    let p_defect = (v(state.motion.rigid.vel_mps).scale(assembly.mass_kg()) - p0 - dp).norm();
    let l_defect = (angular_momentum(&assembly, &state.motion, origin) - l0 - dl).norm();
    println!(
        "24-tick harmonic propulsion: P={p_defect:e}, L={l_defect:e}, modal={max_modal:e}, wake work={max_work:e}"
    );
    assert!(p_defect < 2e-9);
    assert!(l_defect < 2e-8);
    assert!(max_modal < 2e-12);
    assert!(max_work < 2e-11);
    assert!((state.induction[0].first_mps[0] - before.induction[0].first_mps[0]).abs() > 0.01);
    assert!(
        (assembly.speed(&state.motion).unwrap() - assembly.speed(&before.motion).unwrap()).abs()
            > 0.01
    );
}

fn integrate(
    model: &HarmonicPropulsion<'_, '_, 2>,
    start: HarmonicState,
    n: usize,
) -> (HarmonicState, f64) {
    let mut state = start;
    let mut work = 0.0;
    let dt = 0.096 / n as f64;
    for tick in 0..n {
        let step = model
            .step(
                &state,
                tick as f64 * dt,
                dt,
                ZERO,
                |t, _| Ok(inputs(t)),
                |t, p| Ok([0.02 * p[1] + 0.1 * t, 0.0, 0.0]),
            )
            .unwrap();
        work += step
            .coupling_work_residual_j
            .into_iter()
            .map(f64::abs)
            .sum::<f64>();
        state = step.state;
    }
    (state, work)
}
fn distance(a: &HarmonicState, b: &HarmonicState) -> f64 {
    let modes: f64 = a
        .induction
        .iter()
        .zip(&b.induction)
        .flat_map(|(a, b)| a.modes().into_iter().zip(b.modes()))
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    (modes
        + (v(a.motion.rigid.pos_m) - v(b.motion.rigid.pos_m))
            .norm()
            .powi(2)
        + (v(a.motion.rigid.vel_mps) - v(b.motion.rigid.vel_mps))
            .norm()
            .powi(2)
        + (v(a.motion.rigid.omega_body) - v(b.motion.rigid.omega_body))
            .norm()
            .powi(2)
        + (a.motion.momentum - b.motion.momentum).powi(2)
        + (a.motion.phase_rad - b.motion.phase_rad).powi(2)
        + a.motion
            .rigid
            .quat
            .iter()
            .zip(b.motion.rigid.quat)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>())
    .sqrt()
}

#[test]
fn complete_motion_and_partitioned_induced_work_converge_at_second_order() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let model = HarmonicPropulsion::new(&base, RHO, 32).unwrap();
    let start = initial(&assembly);
    let reference = integrate(&model, start, 384).0;
    let runs = [
        integrate(&model, start, 24),
        integrate(&model, start, 48),
        integrate(&model, start, 96),
    ];
    let errors = runs.each_ref().map(|r| distance(&r.0, &reference));
    let work = runs.each_ref().map(|r| r.1);
    println!("coupled harmonic errors: {errors:?}; induced coupling work defects: {work:?}");
    for e in errors.windows(2) {
        assert!((3.5..4.6).contains(&(e[0] / e[1])));
    }
    for e in work.windows(2) {
        assert!((3.5..4.6).contains(&(e[0] / e[1])));
    }
}

#[test]
fn shaft_axis_conventions_and_late_refusal_preserve_complete_state() {
    let design = wright_openloop_v1();
    let aa = assembly(design.gross_mass_kg, false);
    let ab = assembly(design.gross_mass_kg, true);
    let ba = GearedPropulsion::new_spatial(&design, &aa, mounts(false)).unwrap();
    let bb = GearedPropulsion::new_spatial(&design, &ab, mounts(true)).unwrap();
    let ma = HarmonicPropulsion::new(&ba, RHO, 32).unwrap();
    let mb = HarmonicPropulsion::new(&bb, RHO, 32).unwrap();
    let sa = initial(&aa);
    let sb = initial(&ab);
    let a = integrate(&ma, sa, 24);
    let b = integrate(&mb, sb, 24);
    assert_eq!(a, b);
    let mut calls = 0;
    let error = ma
        .step(
            &sa,
            0.0,
            0.004,
            ZERO,
            |t, _| {
                calls += 1;
                let mut x = inputs(t);
                if calls == 2 {
                    x.density_kg_m3 = 1.0;
                }
                Ok(x)
            },
            |_, _| Ok(ZERO),
        )
        .unwrap_err();
    assert_eq!(calls, 2);
    assert_eq!(error.code, "harmonic-propulsion-invalid");
    assert_eq!(sa, initial(&aa));
    let run = || {
        ma.step(&sa, 0.0, 0.004, ZERO, |t, _| Ok(inputs(t)), |_, _| Ok(ZERO))
            .unwrap()
    };
    assert_eq!(run(), run());
}
