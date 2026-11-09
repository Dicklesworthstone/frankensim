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
        gas::Trial::<GasHeatVolumeStep>::new(
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
    let tangent = trial
        .tangent(
            &model,
            &motion.accepted().displacement_m,
            &direction,
            &mut || false,
        )
        .unwrap();
    for (&dof, actual) in model.free_dofs().iter().zip(tangent) {
        let force = |t: &gas::Trial| {
            t.stage.internal[dof / 3][dof % 3]
                - t.net_pressure_pa * t.mean_volume_gradient[dof / 3][dof % 3]
        };
        near(actual, (force(&b) - force(&a)) / (2.0 * h), 2e-2, 5e-6);
    }
    let full = model.expand(&direction, true);
    let geometric = patch
        .pressure_segment_tangent_apply(
            &motion.accepted().displacement_m,
            &trial.endpoint.displacement,
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
fn two_gas_curved_motion_closes_pressure_work_while_mechanical_error_refines() {
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
        pressure_work_tolerance_j: 1e-11,
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
        let mut maximum_midpoint_defect = 0.0_f64;
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
            // Independent discarded-midpoint comparison: the original cubic
            // chart must actually exercise the geometric pressure correction.
            let midpoint: Vec<_> = before
                .displacement_m
                .iter()
                .zip(&after.displacement_m)
                .map(|(a, b)| std::array::from_fn::<_, 3, _>(|i| a[i] + 0.5 * (b[i] - a[i])))
                .collect();
            let gradient = patch.pressure_load_n(&midpoint, 1.0, || false).unwrap();
            let midpoint_volume: f64 = gradient
                .iter()
                .zip(&before.displacement_m)
                .zip(&after.displacement_m)
                .flat_map(|((g, a), b)| (0..3).map(move |i| g[i] * (b[i] - a[i])))
                .sum();
            let midpoint_defect = (frame.midpoint_pressure_pa
                - frame.exterior_midpoint_pressure_pa)
                * midpoint_volume
                - frame.gas.work_by_gas_j
                - frame.exterior_gas.unwrap().work_by_gas_j;
            maximum_midpoint_defect = maximum_midpoint_defect.max(midpoint_defect.abs());
        }
        (
            motion.state().q.clone(),
            maximum_pressure_defect,
            maximum_midpoint_defect,
        )
    };
    let (coarse, ec, midpoint_coarse) = run(4);
    let (fine, ef, _) = run(8);
    let (reference, er, _) = run(16);
    let distance = |q: &[f64]| {
        q.iter()
            .zip(&reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    assert!(distance(&coarse) > 2.5 * distance(&fine));
    assert!(
        ec < 1e-11 && ef < 1e-11 && er < 1e-11,
        "pressure transfer residuals: {ec:e}, {ef:e}, {er:e}"
    );
    assert!(
        midpoint_coarse > 1e-12 && ec < 0.1 * midpoint_coarse,
        "path mean {ec:e} must resolve the actual midpoint defect {midpoint_coarse:e}"
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
    let mut work_probe = motion.clone();
    let work = work_probe
        .step(2e-5, &loads, 0.04, Some(-0.01), settings, || false)
        .unwrap();
    assert!(work.pressure_work_defect_j.abs() < 1e-11);
    // The transfer gate still refuses actual roundoff when requested. Exact
    // zero closure is valid and must not be required to manufacture an error.
    if work.pressure_work_defect_j != 0.0 {
        assert!(matches!(
            motion.step(
                2e-5,
                &loads,
                0.04,
                Some(-0.01),
                GasMembraneThermalControls {
                    pressure_work_tolerance_j: 0.5 * work.pressure_work_defect_j.abs(),
                    ..settings
                },
                || false
            ),
            Err(GasMembraneMotionError::PressureWork { .. })
        ));
        unchanged(&motion);
    }
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

#[test]
fn g1_gibbs_curved_tangent_retains_original_ideal_and_pr_work_efforts() {
    let (patch, reference, fixed, _) = curved_geometry();
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
    let settings = GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 390.0,
            maximum_temperature_k: 450.0,
            ..controls().gas
        },
        exterior_gas: Some(GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 330.0,
            ..controls().gas
        }),
        ..controls()
    };
    let direction = [
        0.4, -0.3, 0.2, -0.2, 0.1, -0.5, 0.3, 0.2, -0.4, -0.1, 0.5, 0.3,
    ];
    let stage: Vec<_> = motion
        .state()
        .q
        .iter()
        .zip(direction)
        .map(|(&q, d)| q + 0.003 * d)
        .collect();
    let trial_at = |stage: &[f64]| {
        let endpoint: Vec<_> = motion
            .state()
            .q
            .iter()
            .zip(stage)
            .map(|(&q, &s)| q + 2.0 * (s - q))
            .collect();
        gas::Trial::<GasDiscreteGradientStep>::new(
            &model,
            motion.accepted(),
            stage,
            &endpoint,
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
        .discrete_gradient_heat_volume(trial.endpoint.volume, 0.2, settings.gas, || false)
        .unwrap();
    assert_eq!(
        trial.gas.endpoint,
        MembraneGas::Ideal(native.endpoint().clone())
    );
    assert_eq!(trial.gas.report, *native.step());
    assert_eq!(trial.gas.work_pressure_pa, native.step().work_pressure_pa);
    assert_eq!(
        trial.gas.sampling_temperature_k,
        native.thermodynamic_mean().state().temperature
    );
    assert_eq!(
        trial.gas.sampling_pressure_pa,
        native.thermodynamic_mean().state().pressure
    );
    let MembraneGas::PengRobinson(second) = motion.accepted().exterior_gas.as_ref().unwrap() else {
        panic!("original PR source")
    };
    let native = second
        .clone()
        .into_isochoric()
        .discrete_gradient_heat_volume(
            trial.endpoint.exterior_volume.unwrap(),
            -0.1,
            settings.exterior_gas.unwrap(),
            || false,
        )
        .unwrap();
    let exterior = trial.exterior.as_ref().unwrap();
    assert_eq!(
        exterior.endpoint,
        MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(
            native.endpoint().clone()
        ))
    );
    assert_eq!(exterior.report, *native.step());
    assert_eq!(exterior.work_pressure_pa, native.step().work_pressure_pa);
    assert_eq!(
        exterior.sampling_pressure_pa,
        native.thermodynamic_mean().pressure_pa()
    );
    assert_eq!(
        exterior.sampling_temperature_k,
        native.thermodynamic_mean().temperature_k()
    );
    assert!(
        (exterior.work_pressure_pa - exterior.sampling_pressure_pa).abs() > 1e-5,
        "the PR fixture must distinguish EOS sample and discrete work pressure"
    );
    assert_eq!(
        trial.net_pressure_pa,
        trial.gas.work_pressure_pa - exterior.work_pressure_pa
    );
    let h = 2e-7;
    let low: Vec<_> = stage
        .iter()
        .zip(direction)
        .map(|(&q, d)| q - h * d)
        .collect();
    let high: Vec<_> = stage
        .iter()
        .zip(direction)
        .map(|(&q, d)| q + h * d)
        .collect();
    let (a, b) = (trial_at(&low), trial_at(&high));
    let tangent = trial
        .tangent(
            &model,
            &motion.accepted().displacement_m,
            &direction,
            &mut || false,
        )
        .unwrap();
    for (&dof, action) in model.free_dofs().iter().zip(tangent) {
        let force = |t: &gas::Trial<GasDiscreteGradientStep>| {
            t.stage.internal[dof / 3][dof % 3]
                - t.net_pressure_pa * t.mean_volume_gradient[dof / 3][dof % 3]
        };
        near(action, (force(&b) - force(&a)) / (2.0 * h), 3e-2, 8e-6);
    }
}

#[test]
fn g1_gibbs_membrane_motion_exchanges_heat_and_preserves_adiabatic_entropy() {
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
    let zero = vec![[0.0; 3]; reference.len()];
    let loads = model
        .evaluate(&reference, &zero, motion_controls().gas, || false)
        .unwrap()
        .residual_n;
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        pressure_work_tolerance_j: 1e-11,
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
    // Continue a genuine accepted midpoint trajectory through the new law.
    motion
        .step(1e-5, &loads, 0.01, Some(-0.002), settings, || false)
        .unwrap();
    let prefix = motion.state().clone();
    for (heat, exterior_heat) in [(0.0, 0.0), (0.04, -0.01), (-0.02, 0.005), (0.0, 0.0)] {
        let before = motion.accepted().clone();
        let old_kinetic = motion.kinetic_energy_j().unwrap();
        let frame = motion
            .step_discrete_gradient(1e-5, &loads, heat, Some(exterior_heat), settings, || false)
            .unwrap();
        let after = motion.accepted();
        let native = frame.exterior_gas.unwrap();
        for (old, new, receipt, theta, supplied) in [
            (
                &before.gas,
                &after.gas,
                frame.gas,
                frame.gas_sampling_temperature_k,
                heat,
            ),
            (
                before.exterior_gas.as_ref().unwrap(),
                after.exterior_gas.as_ref().unwrap(),
                native,
                frame.exterior_sampling_temperature_k.unwrap(),
                exterior_heat,
            ),
        ] {
            let (MembraneGas::Ideal(a), MembraneGas::Ideal(b)) = (old, new) else {
                panic!("variable-Cv ideal source fixture")
            };
            let (t0, t1) = (a.state().state().temperature, b.state().state().temperature);
            let r = a.state().state().specific_gas_constant;
            let energy = analytic_energy(old, new);
            // Independent integral of Cv(T)/T for the original 75:25 NASA fixture.
            let entropy = a.mass_kg()
                * r
                * (2.25 * (t1 / t0).ln()
                    + 0.0015 * (t1 - t0)
                    + (b.volume_m3() / a.volume_m3()).ln());
            let allowance = a.mass_kg() * settings.gas.energy_tolerance_j_kg;
            near(receipt.stored_energy_change_j, energy, 1e-9, 1e-8);
            near(receipt.entropy_change_j_k, entropy, 1e-12, 1e-7);
            near(receipt.entropy_transfer_j_k, supplied / theta, 1e-15, 1e-13);
            assert!(receipt.entropy_residual_j_k.abs() <= allowance / theta);
            assert!(receipt.gibbs_residual_j.abs() <= allowance);
            assert!(receipt.energy_residual_j.abs() <= allowance);
            if supplied == 0.0 {
                assert_eq!(receipt.entropy_transfer_j_k, 0.0);
                assert!(entropy.abs() <= allowance / theta + 1e-12);
            }
            assert_eq!(old.mass_kg(), new.mass_kg());
        }
        assert_eq!(frame.work_pressure_pa, frame.gas.work_pressure_pa);
        assert_eq!(frame.exterior_work_pressure_pa, native.work_pressure_pa);
        near(
            after.gas.volume_m3() + after.exterior_gas.as_ref().unwrap().volume_m3(),
            0.0022,
            1e-18,
            0.0,
        );
        let work: f64 = loads
            .iter()
            .zip(&before.displacement_m)
            .zip(&after.displacement_m)
            .flat_map(|((f, a), b)| (0..3).map(move |i| f[i] * (b[i] - a[i])))
            .sum();
        let balance = motion.kinetic_energy_j().unwrap() - old_kinetic + after.strain_energy_j
            - before.strain_energy_j
            + analytic_energy(&before.gas, &after.gas)
            + analytic_energy(
                before.exterior_gas.as_ref().unwrap(),
                after.exterior_gas.as_ref().unwrap(),
            )
            - work
            - heat
            - exterior_heat;
        near(frame.integration_defect_j, balance, 3e-9, 1e-7);
        assert!(frame.pressure_work_defect_j.abs() < 1e-11);
    }
    assert_eq!(motion.state().steps, prefix.steps + 4);
    assert_eq!(
        &motion.state().history[..prefix.history.len()],
        prefix.history.as_slice()
    );
    // Scheme selection never resets storage or the original mechanical history.
    let last = motion.state().clone();
    motion
        .step(1e-5, &loads, 0.0, Some(0.0), settings, || false)
        .unwrap();
    assert_eq!(motion.state().steps, last.steps + 1);
    assert_eq!(
        &motion.state().history[..last.history.len()],
        last.history.as_slice()
    );
}

#[test]
fn g4_gibbs_membrane_source_energy_and_late_cancel_refusals_are_atomic() {
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
    let zero = vec![[0.0; 3]; reference.len()];
    let loads = model
        .evaluate(&reference, &zero, motion_controls().gas, || false)
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
        .step_discrete_gradient(1e-5, &loads, 0.02, Some(-0.005), settings, || false)
        .unwrap();
    let state = motion.state().clone();
    let accepted = motion.accepted().clone();
    let unchanged = |motion: &GasMembraneThermalMotion<'_, '_>| {
        assert_eq!(motion.state(), &state);
        assert_eq!(motion.accepted(), &accepted);
    };
    let mut reference_run = motion.clone();
    let mut polls = 0;
    let reference_step = reference_run
        .step_discrete_gradient(1e-5, &loads, 0.02, Some(-0.005), settings, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(
        reference_step.integration_defect_j != 0.0,
        "curved midpoint material force must exercise the independent energy gate"
    );
    assert!(matches!(
        motion.step_discrete_gradient(
            1e-5,
            &loads,
            0.02,
            Some(-0.005),
            GasMembraneThermalControls {
                energy_tolerance_j: 0.0,
                ..settings
            },
            || false,
        ),
        Err(GasMembraneMotionError::Energy { .. })
    ));
    unchanged(&motion);
    if reference_step.pressure_work_defect_j != 0.0 {
        assert!(matches!(
            motion.step_discrete_gradient(
                1e-5,
                &loads,
                0.02,
                Some(-0.005),
                GasMembraneThermalControls {
                    pressure_work_tolerance_j: 0.5 * reference_step.pressure_work_defect_j.abs(),
                    ..settings
                },
                || false,
            ),
            Err(GasMembraneMotionError::PressureWork { .. })
        ));
        unchanged(&motion);
    }
    let MembraneGas::Ideal(back) = accepted.exterior_gas.as_ref().unwrap() else {
        panic!("ideal exterior")
    };
    let t = back.state().state().temperature;
    assert!(matches!(
        motion.step_discrete_gradient(
            1e-5,
            &loads,
            0.02,
            Some(100.0),
            GasMembraneThermalControls {
                exterior_gas: Some(GasEnergySolve {
                    minimum_temperature_k: t - 1e-6,
                    maximum_temperature_k: t + 1e-6,
                    ..settings.gas
                }),
                ..settings
            },
            || false,
        ),
        Err(GasMembraneMotionError::Heat(_))
    ));
    unchanged(&motion);
    let mut count = 0;
    assert!(matches!(
        motion.step_discrete_gradient(1e-5, &loads, 0.02, Some(-0.005), settings, || {
            count += 1;
            count >= polls
        },),
        Err(GasMembraneMotionError::Cancelled)
    ));
    unchanged(&motion);
    let retry = motion
        .step_discrete_gradient(1e-5, &loads, 0.02, Some(-0.005), settings, || false)
        .unwrap();
    assert_eq!(motion.state(), reference_run.state());
    assert_eq!(motion.accepted(), reference_run.accepted());
    assert_eq!(retry.gas, reference_step.gas);
    assert_eq!(retry.exterior_gas, reference_step.exterior_gas);
    assert_eq!(retry.work_pressure_pa, reference_step.work_pressure_pa);
}

#[test]
fn g1_gibbs_pr_membrane_impulse_uses_work_pressure_and_true_endpoint_readback() {
    let (patch, reference, fixed) = pyramid();
    let source = real_gas();
    let exterior = source.pressure_response().unwrap().pressure_pa;
    let model = GasMembrane::new(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        exterior,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut velocity = [[0.0; 3]; 5];
    velocity[0][2] = 0.3;
    let loads = [[0.0; 3]; 5];
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
    let settings = GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 330.0,
            ..controls().gas
        },
        ..controls()
    };
    let old = motion.state().clone();
    let dt = 2e-5;
    let frame = motion
        .step_discrete_gradient(dt, &loads, 20.0, None, settings, || false)
        .unwrap();
    let native = source
        .clone()
        .into_isochoric()
        .discrete_gradient_heat_volume(
            motion.accepted().gas.volume_m3(),
            20.0,
            settings.gas,
            || false,
        )
        .unwrap();
    assert_eq!(frame.gas, *native.step());
    assert_eq!(frame.work_pressure_pa, native.step().work_pressure_pa);
    assert_eq!(
        frame.gas_sampling_pressure_pa,
        native.thermodynamic_mean().pressure_pa()
    );
    assert_eq!(
        frame.gas_sampling_temperature_k,
        native.thermodynamic_mean().temperature_k()
    );
    assert!(frame.exterior_gas.is_none());
    assert!(frame.exterior_sampling_pressure_pa.is_none());
    assert!(frame.exterior_sampling_temperature_k.is_none());
    assert_eq!(frame.exterior_work_pressure_pa, exterior);
    assert!((frame.work_pressure_pa - frame.gas_sampling_pressure_pa).abs() > 1e-4);
    let mut midpoint = reference.clone();
    midpoint[0][2] = 0.5 * (old.q[0] + motion.state().q[0]);
    let internal = patch
        .evaluate(&midpoint, || false)
        .unwrap()
        .internal_gradient_n[0][2];
    let gradient = patch
        .pressure_segment_load_n(&reference, &motion.accepted().displacement_m, 1.0, || false)
        .unwrap()[0][2];
    near(
        patch.nodal_mass_kg()[0] * (motion.state().v[0] - old.v[0]) / dt + internal
            - (frame.work_pressure_pa - exterior) * gradient,
        0.0,
        5e-7,
        0.0,
    );
    let accepted = motion.accepted();
    let actual_pressure = accepted.gas.pressure_response().unwrap().pressure_pa;
    assert_eq!(accepted.pressure.pressure_pa, actual_pressure);
    let endpoint_internal = patch
        .evaluate(&accepted.displacement_m, || false)
        .unwrap()
        .internal_gradient_n[0][2];
    let endpoint_gradient = patch
        .pressure_load_n(&accepted.displacement_m, 1.0, || false)
        .unwrap()[0][2];
    near(
        accepted.residual_n[0][2],
        endpoint_internal - (actual_pressure - exterior) * endpoint_gradient,
        1e-10,
        1e-12,
    );
    assert!(frame.pressure_work_defect_j.abs() < 1e-8);
}
