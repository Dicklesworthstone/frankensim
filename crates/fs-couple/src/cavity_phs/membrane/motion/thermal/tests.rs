//! Native physics regressions; shared synthetic source fixtures are not measurements.
use super::*;
use crate::cavity_phs::membrane::tests::motion::{motion_controls, real_gas};
use crate::cavity_phs::membrane::tests::{chart, near, pyramid, source_gas};
use fs_solid::mesh2::Mesh2;

fn controls() -> GasMembraneThermalControls {
    let old = motion_controls();
    GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        exterior_gas: None,
        time: old.time,
        maximum_evaluations: 2000,
        maximum_steps: old.maximum_steps,
        maximum_krylov_restart: old.maximum_krylov_restart,
        energy_tolerance_j: 1e-5,
        pressure_work_tolerance_j: 1e-8,
    }
}

#[test]
fn heat_changes_pressure_and_motion_without_restarting_an_accepted_trajectory() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.0),
        &reference,
        &fixed,
        100000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut u = reference.clone();
    u[0][2] = 1e-4;
    let mut v = [[0.0; 3]; 5];
    v[0][2] = 0.03;
    let loads = [[0.0; 3]; 5];
    let mut prefix =
        GasMembraneMotion::new(&model, &u, &v, &loads, motion_controls(), || false).unwrap();
    for _ in 0..2 {
        prefix
            .step(1e-5, &loads, motion_controls(), || false)
            .unwrap();
    }
    let old = prefix.state().clone();
    let old_gas = prefix.accepted().gas.clone();
    let mut hot = prefix.into_thermal();
    assert_eq!(hot.state(), &old);
    assert_eq!(hot.accepted().gas, old_gas);
    let mut unheated = hot.clone();
    let (dt, heat) = (2e-5, 0.4);
    let frame = hot
        .step(dt, &loads, heat, None, controls(), || false)
        .unwrap();
    unheated
        .step(dt, &loads, 0.0, None, controls(), || false)
        .unwrap();
    assert!(hot.accepted().pressure.pressure_pa > unheated.accepted().pressure.pressure_pa);
    assert!(hot.state().q[0] > unheated.state().q[0]);
    assert_eq!(hot.state().steps, old.steps + 1);
    assert_eq!(
        &hot.state().history[..old.history.len()],
        old.history.as_slice()
    );
    assert_eq!(hot.accepted().gas.mass_kg(), old_gas.mass_kg());
    let v0 = old_gas.volume_m3();
    let v1 = hot.accepted().gas.volume_m3();
    let p0 = old_gas.pressure_response().unwrap().pressure_pa;
    // Independent constant-Cv midpoint pressure for the 75:25 source fixture.
    let gamma_minus_one = 1.0 / 2.25;
    let pm = (p0 * v0 + 0.5 * gamma_minus_one * heat)
        / (0.5 * (v0 + v1) + 0.5 * gamma_minus_one * (v1 - v0));
    near(frame.midpoint_pressure_pa, pm, 1e-6, 1e-11);
    let mut midpoint = reference.clone();
    midpoint[0][2] = 0.5 * (old.q[0] + hot.state().q[0]);
    let elastic = patch
        .evaluate(&midpoint, || false)
        .unwrap()
        .internal_gradient_n[0][2];
    let b = patch.pressure_load_n(&midpoint, 1.0, || false).unwrap()[0][2];
    near(
        patch.nodal_mass_kg()[0] * (hot.state().v[0] - old.v[0]) / dt + elastic
            - (pm - 100000.0) * b,
        0.0,
        3e-7,
        0.0,
    );
    assert_eq!(frame.gas.supplied_heat_j, heat);
    assert!(frame.gas_trials < frame.evaluations);
}

#[test]
fn curved_heated_tangent_uses_native_ideal_and_pr_discrete_derivatives() {
    let (patch, mut reference, fixed) =
        chart(Mesh2::triangles(0.3, 0.2, 3, 2), &[15, 16, 17, 18, 19, 20]);
    reference[5][2] = 0.006;
    reference[6][2] = 0.009;
    let mut settings = controls();
    settings.gas.minimum_temperature_k = 390.0;
    settings.gas.maximum_temperature_k = 450.0;
    settings.exterior_gas = Some(GasEnergySolve {
        minimum_temperature_k: 280.0,
        maximum_temperature_k: 330.0,
        ..settings.gas
    });
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            gas: real_gas().into(),
            volume_interval_m3: [0.0007, 0.0015],
            entropy: GasEntropySolve {
                minimum_temperature_k: 280.0,
                maximum_temperature_k: 330.0,
                ..motion_controls().gas
            },
        },
        || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; reference.len()];
    let motion =
        GasMembraneMotion::new(&model, &reference, &zero, &zero, motion_controls(), || {
            false
        })
        .unwrap()
        .into_thermal();
    let shift = [0.0001, -0.00008, 0.00015, -0.00006, 0.00005, -0.00012];
    let stage: Vec<_> = motion
        .state()
        .q
        .iter()
        .zip(shift)
        .map(|(&q, d)| q + d)
        .collect();
    let endpoint = |stage: &[f64]| {
        motion
            .state()
            .q
            .iter()
            .zip(stage)
            .map(|(&old, &q)| old + 2.0 * (q - old))
            .collect::<Vec<_>>()
    };
    let trial_at = |stage: &[f64]| {
        gas::Trial::new(
            &model,
            motion.accepted(),
            stage,
            &endpoint(stage),
            0.2,
            Some(-0.1),
            settings,
            &mut || false,
        )
        .unwrap()
    };
    let trial = trial_at(&stage);
    let MembraneGas::Ideal(first) = &motion.accepted().gas else {
        panic!("ideal source")
    };
    let native = first
        .midpoint_heat_volume(trial.endpoint.volume, 0.2, settings.gas, || false)
        .unwrap();
    assert_eq!(
        trial.gas.endpoint,
        MembraneGas::Ideal(native.endpoint().clone())
    );
    assert_eq!(trial.gas.report, *native.step());
    let MembraneGas::PengRobinson(second) = motion.accepted().exterior_gas.as_ref().unwrap() else {
        panic!("original PR source")
    };
    let native = second
        .clone()
        .into_isochoric()
        .midpoint_heat_volume(
            trial.endpoint.exterior_volume.unwrap(),
            -0.1,
            settings.exterior_gas.unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(
        trial.exterior.as_ref().unwrap().endpoint,
        MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(
            native.endpoint().clone()
        ))
    );
    assert_eq!(trial.exterior.as_ref().unwrap().report, *native.step());
    let direction = [0.4, -0.3, 0.2, -0.2, 0.1, -0.5];
    let h = 2e-7;
    let lo: Vec<_> = stage
        .iter()
        .zip(direction)
        .map(|(&q, d)| q - h * d)
        .collect();
    let hi: Vec<_> = stage
        .iter()
        .zip(direction)
        .map(|(&q, d)| q + h * d)
        .collect();
    let a = trial_at(&lo);
    let b = trial_at(&hi);
    let tangent = trial.tangent(&model, &direction, &mut || false).unwrap();
    for (&dof, actual) in model.free_dofs().iter().zip(tangent) {
        let force = |t: &gas::Trial| {
            t.stage.internal[dof / 3][dof % 3]
                - t.net_pressure_pa * t.stage.gradient[dof / 3][dof % 3]
        };
        near(actual, (force(&b) - force(&a)) / (2.0 * h), 2e-2, 5e-6);
    }
    let full = model.expand(&direction, true);
    let geometric = patch
        .pressure_tangent_apply(
            &trial.stage.displacement,
            &full,
            trial.net_pressure_pa,
            || false,
        )
        .unwrap();
    assert!(
        model
            .free_dofs()
            .iter()
            .any(|&i| geometric[i / 3][i % 3].abs() > 1.0)
    );
}

fn curved_geometry() -> (
    MembranePatch,
    Vec<[f64; 3]>,
    Vec<Option<f64>>,
    Vec<[f64; 3]>,
) {
    // At least one triangle has THREE free vertices: its swept-volume chart is
    // genuinely cubic, so midpoint traction does not integrate its work exactly.
    let (patch, mut u, fixed) = chart(
        Mesh2::triangles(0.3, 0.3, 3, 3),
        &[15, 16, 17, 18, 19, 20, 27, 28, 29, 30, 31, 32],
    );
    let mut v = vec![[0.0; 3]; u.len()];
    for (node, height, velocity) in [
        (5, 0.004, [1.2, -0.7, 2.1]),
        (6, 0.007, [-0.8, 1.4, -1.1]),
        (9, 0.008, [0.6, 0.9, -1.7]),
        (10, 0.006, [-1.5, -0.4, 0.8]),
    ] {
        u[node][2] = height;
        v[node] = velocity;
    }
    (patch, u, fixed, v)
}
fn backing() -> SealedMembraneExterior {
    let mut gas = source_gas(0.75, 0.002);
    gas.set_volume_isentropic(0.0012, motion_controls().gas, || false)
        .unwrap();
    SealedMembraneExterior {
        gas: gas.into(),
        volume_interval_m3: [0.0007, 0.0015],
        entropy: motion_controls().gas,
    }
}
fn analytic_energy(old: &MembraneGas, new: &MembraneGas) -> f64 {
    let (MembraneGas::Ideal(a), MembraneGas::Ideal(b)) = (old, new) else {
        panic!("ideal fixture")
    };
    let (t0, t1) = (a.state().state().temperature, b.state().state().temperature);
    a.mass_kg() * a.state().state().specific_gas_constant * (t1 - t0) * (2.25 + 0.00075 * (t1 + t0))
}

#[test]
fn two_gas_heat_work_accounts_and_curved_pressure_error_refine_with_native_steps() {
    let (patch, reference, fixed, velocity) = curved_geometry();
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        backing(),
        || false,
    )
    .unwrap();
    let loads = model
        .evaluate(
            &reference,
            &vec![[0.0; 3]; reference.len()],
            motion_controls().gas,
            || false,
        )
        .unwrap()
        .residual_n;
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        ..controls()
    };
    let run = |steps: usize| {
        let mut motion = GasMembraneMotion::new(
            &model,
            &reference,
            &velocity,
            &loads,
            motion_controls(),
            || false,
        )
        .unwrap()
        .into_thermal();
        let dt = 8e-5 / (steps as f64);
        let mut maximum_pressure_defect = 0.0_f64;
        for _ in 0..steps {
            let before = motion.accepted().clone();
            let k0 = motion.kinetic_energy_j().unwrap();
            let frame = motion
                .step(dt, &loads, 2000.0 * dt, Some(-500.0 * dt), settings, || {
                    false
                })
                .unwrap();
            let after = motion.accepted();
            let du1 = analytic_energy(&before.gas, &after.gas);
            let du2 = analytic_energy(
                before.exterior_gas.as_ref().unwrap(),
                after.exterior_gas.as_ref().unwrap(),
            );
            near(frame.gas.stored_energy_change_j, du1, 1e-9, 1e-8);
            near(
                frame.exterior_gas.unwrap().stored_energy_change_j,
                du2,
                1e-9,
                1e-8,
            );
            let work: f64 = loads
                .iter()
                .zip(&before.displacement_m)
                .zip(&after.displacement_m)
                .flat_map(|((f, a), b)| (0..3).map(move |i| f[i] * (b[i] - a[i])))
                .sum();
            let enclosing = motion.kinetic_energy_j().unwrap() - k0 + after.strain_energy_j
                - before.strain_energy_j
                + du1
                + du2
                - work
                - 1500.0 * dt;
            near(frame.integration_defect_j, enclosing, 3e-9, 1e-7);
            assert_eq!(frame.exterior_pressure_work_j, 0.0);
            near(
                after.gas.volume_m3() + after.exterior_gas.as_ref().unwrap().volume_m3(),
                0.0022,
                1e-18,
                0.0,
            );
            maximum_pressure_defect =
                maximum_pressure_defect.max(frame.pressure_work_defect_j.abs());
        }
        (motion.state().q.clone(), maximum_pressure_defect)
    };
    let (coarse, ec) = run(4);
    let (fine, ef) = run(8);
    let (reference, er) = run(16);
    let distance = |q: &[f64]| {
        q.iter()
            .zip(&reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    assert!(distance(&coarse) > 2.5 * distance(&fine));
    assert!(
        ec > 1e-12 && ec > 3.0 * ef && ef > 3.0 * er,
        "pressure errors: {ec:e}, {ef:e}, {er:e}"
    );
}

#[test]
fn source_and_late_work_energy_cancel_refusals_preserve_both_gases_and_history() {
    let (patch, reference, fixed, velocity) = curved_geometry();
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        backing(),
        || false,
    )
    .unwrap();
    let loads = model
        .evaluate(
            &reference,
            &vec![[0.0; 3]; reference.len()],
            motion_controls().gas,
            || false,
        )
        .unwrap()
        .residual_n;
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        ..controls()
    };
    let mut motion = GasMembraneMotion::new(
        &model,
        &reference,
        &velocity,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    motion
        .step(2e-5, &loads, 0.04, Some(-0.01), settings, || false)
        .unwrap();
    let state = motion.state().clone();
    let accepted = motion.accepted().clone();
    let unchanged = |motion: &GasMembraneThermalMotion<'_, '_>| {
        assert_eq!(motion.state(), &state);
        assert_eq!(motion.accepted(), &accepted);
    };
    assert!(matches!(
        motion.step(
            2e-5,
            &loads,
            0.04,
            Some(-0.01),
            GasMembraneThermalControls {
                pressure_work_tolerance_j: 0.0,
                ..settings
            },
            || false
        ),
        Err(GasMembraneMotionError::PressureWork { .. })
    ));
    unchanged(&motion);
    assert!(matches!(
        motion.step(
            2e-5,
            &loads,
            0.04,
            Some(-0.01),
            GasMembraneThermalControls {
                energy_tolerance_j: 0.0,
                ..settings
            },
            || false
        ),
        Err(GasMembraneMotionError::Energy { .. })
    ));
    unchanged(&motion);
    let MembraneGas::Ideal(back) = accepted.exterior_gas.as_ref().unwrap() else {
        panic!("ideal exterior")
    };
    let t = back.state().state().temperature;
    let source_limit = GasMembraneThermalControls {
        exterior_gas: Some(GasEnergySolve {
            minimum_temperature_k: t - 1e-6,
            maximum_temperature_k: t + 1e-6,
            ..settings.gas
        }),
        ..settings
    };
    assert!(matches!(
        motion.step(2e-5, &loads, 0.04, Some(100.0), source_limit, || false),
        Err(GasMembraneMotionError::Heat(_))
    ));
    unchanged(&motion);
    let mut probe = motion.clone();
    let mut polls = 0usize;
    probe
        .step(2e-5, &loads, 0.04, Some(-0.01), settings, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0usize;
    assert!(matches!(
        motion.step(2e-5, &loads, 0.04, Some(-0.01), settings, || {
            seen += 1;
            seen >= polls
        }),
        Err(GasMembraneMotionError::Cancelled)
    ));
    unchanged(&motion);
    motion
        .step(2e-5, &loads, 0.04, Some(-0.01), settings, || false)
        .unwrap();
    assert_eq!(motion.state(), probe.state());
    assert_eq!(motion.accepted(), probe.accepted());
}
