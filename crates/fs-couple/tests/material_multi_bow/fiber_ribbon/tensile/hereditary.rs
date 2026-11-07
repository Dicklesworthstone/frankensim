//! Real axial/bending Maxwell memories participate in the original transaction.
use super::*;

fn hereditary(
    cells: usize,
    limit: f64,
    eta: f64,
    selected: bool,
    gate: &CancelGate,
) -> TensionedFiberRibbon {
    let source = ribbon_with_viscosities(cells, limit, eta, eta, gate);
    let spectrum = source
        .axial_material()
        .elastic()
        .clone()
        .with_prony_bending(&[LayerRelaxationProperties {
            equilibrium_modulus: "young_modulus".into(),
            branches: vec![("delta".into(), "tau".into())],
        }])
        .unwrap();
    let source = source
        .with_material_relaxation(spectrum, 128, gate)
        .unwrap();
    if selected {
        source.with_tensile_limit(LIMIT).unwrap()
    } else {
        source
    }
}

// Independent original-fiber reconstruction for this one-arm, one-harmonic
// manufactured fixture. Both memory slices come from the native accepted owner.
fn stresses(bank: &MultiBowStringBank, eta: f64) -> (f64, f64, f64) {
    let fiber = bank.source().fiber_ribbon().unwrap();
    let source = fiber.model();
    let state = bank.system().components()[fiber.component_start()].states();
    let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
    let (u, y) = (
        state[0].displacement_m_sqrt_kg / scale,
        state[1].displacement_m_sqrt_kg / scale,
    );
    let (udot, ydot) = (
        state[0].velocity_m_sqrt_kg_per_s / scale,
        state[1].velocity_m_sqrt_kg_per_s / scale,
    );
    let k = core::f64::consts::PI / 0.6;
    let mean = 0.25 * (k * y).powi(2);
    let harmonic = k * u;
    let base = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
    let mean_stress = base + YOUNG * mean + eta * 0.5 * k * k * y * ydot;
    let axial_stress = YOUNG * harmonic + eta * k * udot;
    let bend = YOUNG * y + eta * ydot;
    let channels = source.material_relaxation().unwrap().channels();
    assert_eq!(channels.len(), 1);
    assert_eq!(
        bank.system().material_memory().len(),
        source.resolution().cells
    );
    let histories = bank.system().nonlinear_material_memory();
    assert_eq!(histories.len(), source.resolution().cells);
    assert_eq!(histories[0].len(), 2);
    let relaxed_y = bank.system().material_memory()[0] / channels[0].coupling_per_s / scale;
    let complete = mean_stress
        + YOUNG * (mean - histories[0][0])
        + (axial_stress + YOUNG * (harmonic - histories[0][1])).abs()
        + RADIUS * k * k * (bend + YOUNG * (y - relaxed_y)).abs();
    let no_memory = mean_stress + axial_stress.abs() + RADIUS * k * k * bend.abs();
    let zero_memory = mean_stress
        + YOUNG * mean
        + (axial_stress + YOUNG * harmonic).abs()
        + RADIUS * k * k * (bend + YOUNG * y).abs();
    (complete, no_memory, zero_memory)
}

#[test]
fn hereditary_tensile_crossing_uses_both_actual_histories_and_rolls_back_every_owner() {
    with_cx(|cx, gate| {
        for (mode, amplitude) in [(0, 1.0), (1, 1000.0)] {
            let eta = 30_000.0;
            let mut control = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    hereditary(3, 100e6, eta, true, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap();
            let mut unchecked = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    hereditary(3, 100e6, eta, false, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
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
                control.system().material_memory(),
                unchecked.system().material_memory()
            );
            assert_eq!(
                control.system().nonlinear_material_memory(),
                unchecked.system().nonlinear_material_memory()
            );
            assert_eq!(
                control.system().total_energy_j().unwrap(),
                unchecked.system().total_energy_j().unwrap()
            );
            let (actual, missing, zero) = stresses(&control, eta);
            assert!(
                actual > missing && zero > actual,
                "mode {mode}: {missing:e} {actual:e} {zero:e}"
            );

            // An equilibrium/Kelvin-only gate would admit this same endpoint.
            let mut refused = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    hereditary(3, f64::midpoint(actual, missing), eta, true, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap();
            let before = states(&refused);
            let linear = refused.system().material_memory().to_vec();
            let nonlinear: Vec<_> = refused
                .system()
                .nonlinear_material_memory()
                .iter()
                .map(|z| z.to_vec())
                .collect();
            let energy = refused.system().total_energy_j().unwrap();
            let error = refused.step(&force, &drives).unwrap_err();
            assert!(format!("{error:?}").contains("tensile"), "{error:?}");
            assert_eq!(states(&refused), before);
            assert_eq!(refused.system().material_memory(), linear.as_slice());
            assert_eq!(
                refused.system().nonlinear_material_memory(),
                nonlinear.iter().map(Vec::as_slice).collect::<Vec<_>>()
            );
            assert_eq!(refused.system().total_energy_j().unwrap(), energy);
            assert_eq!(refused.system().samples_rendered(), 0);
            assert!(refused.system().last_frame().is_none());

            // A fabricated zero-memory endpoint would wrongly reject this
            // source. Actual staged relaxation must admit it unchanged.
            let mut admitted = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    hereditary(3, f64::midpoint(actual, zero), eta, true, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap();
            assert_eq!(admitted.step(&force, &drives).unwrap(), &frame);
            force.fill(0.0);
            refused.step(&force, &drives).unwrap();
            assert_eq!(refused.system().samples_rendered(), 1);
        }
    });
}

#[test]
fn hereditary_cold_restore_retains_memory_and_checks_axial_and_bending_stress() {
    with_cx(|cx, gate| {
        let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
        for mode in [0, 1] {
            let source = hereditary(3, preload + 100.0, 0.0, true, gate);
            let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
            let bank = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(source, 0.4, ribbon_limits(), gate)
                .unwrap();
            let component = bank.source().fiber_ribbon().unwrap().component_start();
            let (_, mut system) = bank.into_parts();
            let before = system.components()[component].states().to_vec();
            let energy = system.total_energy_j().unwrap();
            let mut states = [ModalAcousticState::default(); 2];
            let amplitude = if mode == 0 {
                100.0 / (1.5 * YOUNG * core::f64::consts::PI / 0.6)
            } else {
                2e-5
            };
            states[mode].displacement_m_sqrt_kg = amplitude * scale;
            let error = system
                .restore_component_states(component, &states)
                .unwrap_err();
            assert!(format!("{error:?}").contains("tensile"), "{error:?}");
            assert_eq!(system.components()[component].states(), before.as_slice());
            assert_eq!(system.total_energy_j().unwrap(), energy);
            assert!(system.material_memory().iter().all(|z| *z == 0.0));
            assert!(
                system
                    .nonlinear_material_memory()
                    .iter()
                    .all(|z| z.iter().all(|z| *z == 0.0))
            );
            assert_eq!(system.samples_rendered(), 0);
            states[mode].displacement_m_sqrt_kg *= 0.1;
            system.restore_component_states(component, &states).unwrap();
            assert_eq!(system.components()[component].states(), states.as_slice());
            assert!(system.material_memory().iter().all(|z| *z == 0.0));
            assert!(
                system
                    .nonlinear_material_memory()
                    .iter()
                    .all(|z| z.iter().all(|z| *z == 0.0))
            );
        }
    });
}

#[test]
fn hereditary_admission_combines_signed_elastic_viscous_and_memory_stress_before_bounding() {
    with_cx(|cx, gate| {
        let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
        let eta = 30_000.0;
        let source = hereditary(3, preload + 100.0, eta, true, gate);
        let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
        let bank = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(source, 0.4, ribbon_limits(), gate)
            .unwrap();
        let component = bank.source().fiber_ribbon().unwrap().component_start();
        let (_, mut system) = bank.into_parts();
        let mut state = [ModalAcousticState::default(); 2];
        state[0].displacement_m_sqrt_kg = 1e-5 * scale;
        // The original zero material memory supplies +DeltaE*epsilon.
        // The actual unloading rate cancels E*epsilon + DeltaE*epsilon.
        state[0].velocity_m_sqrt_kg_per_s = -2.0 * YOUNG / eta * state[0].displacement_m_sqrt_kg;
        system.restore_component_states(component, &state).unwrap();
        assert_eq!(system.components()[component].states(), state.as_slice());
        let energy = system.total_energy_j().unwrap();
        let mut incomplete = state;
        // Cancelling elastic stress alone leaves real tensile Maxwell stress.
        incomplete[0].velocity_m_sqrt_kg_per_s *= 0.5;
        let error = system
            .restore_component_states(component, &incomplete)
            .unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(system.components()[component].states(), state.as_slice());
        assert_eq!(system.total_energy_j().unwrap(), energy);
        assert!(system.material_memory().iter().all(|z| *z == 0.0));
        assert!(
            system
                .nonlinear_material_memory()
                .iter()
                .all(|z| z.iter().all(|z| *z == 0.0))
        );
    });
}

#[test]
fn hereditary_tensile_stress_is_invariant_under_physical_width_partition() {
    with_cx(|cx, gate| {
        let eta = 30_000.0;
        let mut previous = None;
        for cells in [1, 3, 9] {
            let mut bank = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    hereditary(cells, 100e6, eta, true, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap();
            let first = bank.source().fiber_ribbon().unwrap().component_start();
            let start: usize = bank.system().components()[..first]
                .iter()
                .map(|c| c.modes().len())
                .sum();
            let mut force = vec![0.0; bank.system().mode_count()];
            // Same physical load per original fiber. Each strip's modal mass
            // and force scale together, rather than holding strip force fixed.
            force[start] = 1.0 / (cells as f64).sqrt();
            force[start + 1] = 1000.0 / (cells as f64).sqrt();
            bank.step(&force, &vec![drive(0.0, 0.0); cells]).unwrap();
            let stress = stresses(&bank, eta).0;
            if let Some(previous) = previous {
                assert!((stress / previous - 1.0).abs() < 2e-12);
            }
            previous = Some(stress);
        }
    });
}

#[test]
fn hereditary_layers_retain_independent_moduli_viscosities_memories_and_limits() {
    with_cx(|cx, gate| {
        let radii = [32e-6_f64, RADIUS];
        let young = [6e9, YOUNG];
        let eta = [60_000.0, 300_000.0];
        let area = [
            core::f64::consts::PI * radii[0].powi(2),
            core::f64::consts::PI * (radii[1].powi(2) - radii[0].powi(2)),
        ];
        let preload_strain = 40.0 / (180.0 * (young[0] * area[0] + young[1] * area[1]));
        let u = 1e-5;
        let speed = -2.0 * young[0] / eta[0] * u;
        let k = core::f64::consts::PI / 0.6;
        let outer_increment = (2.0 * young[1] * k * u + eta[1] * k * speed).abs();
        assert!(outer_increment > 1e6);
        for admitted in [true, false] {
            let limits = [
                young[0] * preload_strain + 100.0,
                young[1] * preload_strain + if admitted { 1.1 } else { 0.9 } * outer_increment,
            ];
            let section = ConcentricElasticSection::try_new(
                0.0,
                (0..2)
                    .map(|layer| ConcentricMaterialLayer {
                        outer_radius_m: radii[layer],
                        material: fiber_material_with_limit(
                            1300.0,
                            young[layer],
                            eta[layer],
                            0.0,
                            Some(limits[layer]),
                        ),
                    })
                    .collect(),
            )
            .unwrap();
            let spectrum = section
                .clone()
                .with_prony_bending(&vec![
                    LayerRelaxationProperties {
                        equilibrium_modulus: "young_modulus".into(),
                        branches: vec![("delta".into(), "tau".into())],
                    };
                    2
                ])
                .unwrap();
            let source = ribbon_with_section(3, section, 1, 1, gate)
                .with_tensile_limit(LIMIT)
                .unwrap()
                .with_material_relaxation(spectrum, 128, gate)
                .unwrap();
            let scale = (0.5 * source.cell_line_density_kg_m() * 0.6).sqrt();
            let bank = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(source, 0.4, ribbon_limits(), gate)
                .unwrap();
            let component = bank.source().fiber_ribbon().unwrap().component_start();
            let (_, mut system) = bank.into_parts();
            let energy = system.total_energy_j().unwrap();
            let mut state = [ModalAcousticState::default(); 2];
            state[0].displacement_m_sqrt_kg = u * scale;
            state[0].velocity_m_sqrt_kg_per_s = speed * scale;
            let result = system.restore_component_states(component, &state);
            if admitted {
                result.unwrap();
                assert_eq!(system.components()[component].states(), state.as_slice());
            } else {
                let error = result.unwrap_err();
                assert!(format!("{error:?}").contains("tensile"), "{error:?}");
                assert_eq!(system.total_energy_j().unwrap(), energy);
                assert!(
                    system.components()[component]
                        .states()
                        .iter()
                        .all(|state| *state == ModalAcousticState::default())
                );
            }
        }
    });
}
