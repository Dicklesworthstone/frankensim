//! Actual endpoint rate stress on the original layer and mechanical owners.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionError, ModalFrictionPort,
    nonlinear::{MappedModalPotentialBinding, ModalCoordinate},
};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment};

#[test]
fn total_viscous_stress_and_force_are_width_refinement_invariant() {
    let gate = CancelGate::new_clock_free();
    let (axial_eta, bending_eta) = (30_000.0, 40_000.0);
    let (axial, normal, axial_speed, normal_speed) = (0.0001, 0.002, 0.025, -0.1);
    let k = core::f64::consts::PI / 0.6;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let mean = 40.0 / (180.0 * YOUNG * area) + 0.25 * (k * normal).powi(2);
    let mean_rate = 0.5 * k * k * normal * normal_speed;
    let axial_stress = YOUNG * mean
        + axial_eta * mean_rate
        + (YOUNG * k * axial + axial_eta * k * axial_speed).abs();
    let bending_stress = RADIUS * k * k * (YOUNG * normal + bending_eta * normal_speed).abs();
    let mut previous = None;
    for cells in [1, 3, 9] {
        let source = ribbon_with_viscosities(cells, 100e6, axial_eta, bending_eta, &gate)
            .with_tensile_limit(LIMIT)
            .unwrap();
        let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
        let assessment = source
            .strain()
            .unwrap()
            .tensile_state_assessment(
                &[axial * scale, normal * scale],
                &[axial_speed * scale, normal_speed * scale],
            )
            .unwrap()
            .unwrap();
        near(
            assessment.governing_tensile_stress_bound_pa,
            axial_stress + bending_stress,
        );
        near(assessment.per_fiber_tension_bound_n, area * axial_stress);
        near(
            assessment.curvature_rate_bound_per_m_s,
            k * k * normal_speed.abs(),
        );
        if let Some(previous) = previous {
            near(assessment.governing_tensile_stress_bound_pa, previous);
        }
        previous = Some(assessment.governing_tensile_stress_bound_pa);
    }
}

#[test]
fn layered_rates_keep_original_viscosities_and_complete_geometric_harmonics() {
    let gate = CancelGate::new_clock_free();
    let radii = [32e-6, RADIUS];
    let young = [6e9, YOUNG];
    let eta_a = [1e3, 1e6];
    let eta_b = [8e5, 3e5];
    let limits = [200e6, 60e6];
    let section = ConcentricElasticSection::try_new(
        0.0,
        (0..2)
            .map(|i| ConcentricMaterialLayer {
                outer_radius_m: radii[i],
                material: fiber_material_with_limit(
                    1300.0,
                    young[i],
                    eta_a[i],
                    eta_b[i],
                    Some(limits[i]),
                ),
            })
            .collect(),
    )
    .unwrap();
    let preload_strain = 40.0 / (180.0 * section.axial_stiffness_n());
    let source = ribbon_with_section(3, section, 2, 2, &gate)
        .with_tensile_limit(LIMIT)
        .unwrap();
    let k = core::f64::consts::PI / 0.6;
    let u = [1e-4, -2e-5];
    let y = [0.001, 0.0004];
    let udot = [0.4, -0.2];
    let ydot = [-0.1, 0.05];
    let a = [k * y[0], 2.0 * k * y[1]];
    let b = [k * ydot[0], 2.0 * k * ydot[1]];
    let mean = preload_strain + 0.25 * (a[0] * a[0] + a[1] * a[1]);
    let mean_rate = 0.5 * (a[0] * b[0] + a[1] * b[1]);
    let strain = [
        k * u[0] + 0.5 * a[0] * a[1],
        2.0 * k * u[1] + 0.25 * a[0] * a[0],
    ];
    let rates = [
        k * udot[0] + 0.5 * (a[0] * b[1] + b[0] * a[1]),
        2.0 * k * udot[1] + 0.5 * a[0] * b[0],
    ];
    let stress = |i: usize| {
        young[i] * mean
            + eta_a[i] * mean_rate
            + (0..2)
                .map(|j| (young[i] * strain[j] + eta_a[i] * rates[j]).abs())
                .sum::<f64>()
            + radii[i]
                * (0..2)
                    .map(|j| {
                        ((j + 1) as f64 * k).powi(2) * (young[i] * y[j] + eta_b[i] * ydot[j]).abs()
                    })
                    .sum::<f64>()
    };
    assert!(stress(1) / limits[1] > stress(0) / limits[0]);
    let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
    let q = [u[0] * scale, u[1] * scale, y[0] * scale, y[1] * scale];
    let v = [
        udot[0] * scale,
        udot[1] * scale,
        ydot[0] * scale,
        ydot[1] * scale,
    ];
    let assessment = source
        .strain()
        .unwrap()
        .tensile_state_assessment(&q, &v)
        .unwrap()
        .unwrap();
    assert_eq!(assessment.governing_layer, 1);
    near(assessment.governing_tensile_stress_bound_pa, stress(1));
    near(assessment.utilization, stress(1) / limits[1]);
}

fn bind_state(
    source: &TensionedFiberRibbon,
    q: &[f64],
    velocity: &[f64],
    gate: &CancelGate,
) -> Result<FrictionModalSystem, ModalFrictionError> {
    let mut model = source.cell_model().clone();
    model
        .restore_states(
            &q.iter()
                .zip(velocity)
                .map(|(&q, &velocity)| ModalAcousticState {
                    displacement_m_sqrt_kg: q,
                    velocity_m_sqrt_kg_per_s: velocity,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let core = FrictionModalSystem::new(
        CoupledModalSystem::new(vec![model], vec![], options().network, gate).unwrap(),
        ModalFrictionPort {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 0.0],
            },
            right: None,
            law: law(),
        },
        options().friction,
        gate,
    )
    .unwrap();
    core.with_mapped_nonlinear_potentials(
        vec![MappedModalPotentialBinding {
            coordinates: (0..2)
                .map(|mode| ModalCoordinate { component: 0, mode })
                .collect(),
            replaces_linear: vec![0],
            potential: Box::new(source.strain().unwrap().clone()),
        }],
        strain_config(),
        gate,
    )
}

#[test]
fn signed_unloading_cancels_stress_and_cold_binding_uses_actual_velocity() {
    let gate = CancelGate::new_clock_free();
    let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
    let eta = 30_000.0;
    let source = ribbon_with_viscosities(3, preload + 100.0, eta, 0.0, &gate)
        .with_tensile_limit(LIMIT)
        .unwrap();
    let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
    let q = [0.0001 * scale, 0.0];
    let velocity = [-YOUNG / eta * q[0], 0.0];
    let strain = source.strain().unwrap();
    assert!(
        strain
            .tensile_assessment(&q)
            .unwrap()
            .unwrap()
            .exceeds_limit()
    );
    strain.validate_coordinates(&q).unwrap();
    strain.validate_state(&q, &velocity).unwrap();
    let assessment = strain
        .tensile_state_assessment(&q, &velocity)
        .unwrap()
        .unwrap();
    near(assessment.governing_tensile_stress_bound_pa, preload);
    bind_state(&source, &q, &velocity, &gate).unwrap();
    assert!(bind_state(&source, &q, &[-velocity[0], 0.0], &gate).is_err());
    assert!(strain.tensile_state_assessment(&q, &[0.0]).is_err());
    assert!(
        strain
            .tensile_state_assessment(&q, &[f64::NAN, 0.0])
            .is_err()
    );

    // Bending loss remains in the native linear owner, but its signed local
    // stress must still be combined before the fiber's radius/span bound.
    let bend = ribbon_with_viscosities(3, 100e6, 0.0, eta, &gate)
        .with_tensile_limit(LIMIT)
        .unwrap();
    let q = [0.0, 0.002 * scale];
    let velocity = [0.0, -YOUNG / eta * q[1]];
    assert!(!bend.strain().unwrap().has_constitutive_loss());
    let assessment = bend
        .strain()
        .unwrap()
        .tensile_state_assessment(&q, &velocity)
        .unwrap()
        .unwrap();
    let k = core::f64::consts::PI / 0.6;
    near(
        assessment.governing_tensile_stress_bound_pa,
        preload + 0.25 * YOUNG * (k * 0.002).powi(2),
    );
}

#[test]
fn static_reinitialization_cannot_bypass_viscous_material_strength() {
    let gate = CancelGate::new_clock_free();
    let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
    let source = ribbon_with_viscosities(3, preload + 100.0, 30_000.0, 0.0, &gate)
        .with_tensile_limit(LIMIT)
        .unwrap();
    let mut system = bind_state(&source, &[0.0; 2], &[0.0; 2], &gate)
        .unwrap()
        .with_additional_ports(vec![], joint_limits(), &gate)
        .unwrap();
    let before = system.components()[0].states().to_vec();
    let energy = system.total_energy_j().unwrap();
    let error = system
        .initialize_released_equilibrium(&[1.0, 0.0], &gate)
        .unwrap_err();
    assert!(format!("{error:?}").contains("tensile"), "{error:?}");
    assert_eq!(system.components()[0].states(), before.as_slice());
    assert_eq!(system.total_energy_j().unwrap(), energy);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.nonlinear_preload_diagnostics().is_none());
    system
        .initialize_released_equilibrium(&[0.0; 2], &gate)
        .unwrap();
}

#[test]
fn cold_restore_checks_viscous_and_elastic_strength_before_publication() {
    let gate = CancelGate::new_clock_free();
    let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
    for eta in [30_000.0, 0.0] {
        let source = ribbon_with_viscosities(3, preload + 100.0, eta, 0.0, &gate)
            .with_tensile_limit(LIMIT)
            .unwrap();
        let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
        let mut system = bind_state(&source, &[0.0; 2], &[0.0; 2], &gate)
            .unwrap()
            .with_additional_ports(vec![], joint_limits(), &gate)
            .unwrap();
        let before = system.components()[0].states().to_vec();
        let energy = system.total_energy_j().unwrap();
        let mut rejected = [ModalAcousticState::default(); 2];
        let mut accepted = rejected;
        if eta > 0.0 {
            // Zero displacement leaves elastic stress exactly at preload;
            // only the original axial rate stress crosses the 100 Pa margin.
            rejected[0].velocity_m_sqrt_kg_per_s = 0.01 * scale;
            accepted[0].velocity_m_sqrt_kg_per_s = 0.0001 * scale;
        } else {
            // A lossless guard does not request velocity validation, but its
            // coordinate/strength domain must still admit a cold restoration.
            rejected[0].displacement_m_sqrt_kg = 1e-6 * scale;
            accepted[0].displacement_m_sqrt_kg = 1e-9 * scale;
        }
        let error = system.restore_component_states(0, &rejected).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(system.components()[0].states(), before.as_slice());
        assert_eq!(system.total_energy_j().unwrap(), energy);
        assert_eq!(system.samples_rendered(), 0);
        assert!(system.last_frame().is_none());
        assert!(system.nonlinear_preload_diagnostics().is_none());

        system.restore_component_states(0, &accepted).unwrap();
        assert_eq!(system.components()[0].states(), accepted.as_slice());
        assert!(system.total_energy_j().unwrap() > energy);
        assert_eq!(system.samples_rendered(), 0);
        assert!(system.last_frame().is_none());
    }
}

#[test]
fn native_endpoint_rate_crossing_rolls_back_without_changing_admitted_mechanics() {
    with_cx(|cx, gate| {
        for (axial_eta, bending_eta, mode, amplitude) in
            [(30_000.0, 0.0, 0, 1.0), (0.0, 30_000.0, 1, 1000.0)]
        {
            let original = ribbon_with_viscosities(3, 100e6, axial_eta, bending_eta, gate);
            let guarded = original.clone().with_tensile_limit(LIMIT).unwrap();
            let mut control = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(guarded.clone(), 0.4, ribbon_limits(), gate)
                .unwrap();
            let mut unchecked = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(original, 0.4, ribbon_limits(), gate)
                .unwrap();
            let first = control.source().fiber_ribbon().unwrap().component_start();
            let start: usize = control.system().components()[..first]
                .iter()
                .map(|c| c.modes().len())
                .sum();
            let mut force = vec![0.0; control.system().mode_count()];
            force[start + mode] = amplitude;
            let drives = vec![drive(0.0, 0.0); 3];
            let frame = control.step(&force, &drives).unwrap().clone();
            assert_eq!(unchecked.step(&force, &drives).unwrap(), &frame);
            assert_eq!(states(&control), states(&unchecked));
            assert_eq!(
                control.system().total_energy_j().unwrap(),
                unchecked.system().total_energy_j().unwrap()
            );
            let endpoint = control.system().components()[first].states();
            let q: Vec<_> = endpoint.iter().map(|s| s.displacement_m_sqrt_kg).collect();
            let velocity: Vec<_> = endpoint
                .iter()
                .map(|s| s.velocity_m_sqrt_kg_per_s)
                .collect();
            let secant: Vec<_> = q
                .iter()
                .map(|q| q / control.system().sample_period_s())
                .collect();
            let strain = guarded.strain().unwrap();
            let actual = strain
                .tensile_state_assessment(&q, &velocity)
                .unwrap()
                .unwrap()
                .governing_tensile_stress_bound_pa;
            let wrong_rate = strain
                .tensile_state_assessment(&q, &secant)
                .unwrap()
                .unwrap()
                .governing_tensile_stress_bound_pa;
            let elastic = strain
                .tensile_assessment(&q)
                .unwrap()
                .unwrap()
                .governing_tensile_stress_bound_pa;
            assert!(actual > wrong_rate && wrong_rate > elastic);
            let source_limit = f64::midpoint(actual, wrong_rate);
            let source = ribbon_with_viscosities(3, source_limit, axial_eta, bending_eta, gate)
                .with_tensile_limit(LIMIT)
                .unwrap();
            let mut rejected = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(source, 0.4, ribbon_limits(), gate)
                .unwrap();
            let before = states(&rejected);
            let energy = rejected.system().total_energy_j().unwrap();
            let error = rejected.step(&force, &drives).unwrap_err();
            assert!(format!("{error:?}").contains("tensile"), "{error:?}");
            assert_eq!(states(&rejected), before);
            assert_eq!(rejected.system().total_energy_j().unwrap(), energy);
            assert_eq!(rejected.system().samples_rendered(), 0);
            assert!(rejected.system().last_frame().is_none());
            force.fill(0.0);
            rejected.step(&force, &drives).unwrap();
            assert_eq!(rejected.system().samples_rendered(), 1);
        }
    });
}
