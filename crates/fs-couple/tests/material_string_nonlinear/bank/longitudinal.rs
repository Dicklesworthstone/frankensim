//! Dynamic axial modes join the actual source/bow/body/fluid owner.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;

fn selection(modes: &[usize]) -> LongitudinalStringResolution {
    LongitudinalStringResolution {
        axial_modes: modes.to_vec(),
        maximum_slope: 0.2,
        maximum_strain: 0.01,
    }
}

#[test]
fn longitudinal_extension_preserves_original_components_and_excites_real_axial_inertia() {
    with_cx(|cx, gate| {
        let original = physical_bank(cx, gate, 200e9, gas(), 2);
        let modes: Vec<_> = original
            .system()
            .components()
            .iter()
            .map(|c| c.modes().to_vec())
            .collect();
        let states: Vec<_> = original
            .system()
            .components()
            .iter()
            .map(|c| c.states().to_vec())
            .collect();
        let weights = original
            .system()
            .acceleration_pressure_weights()
            .unwrap()
            .to_vec();
        let old_count = original.system().mode_count();
        let mut bank = original
            .with_longitudinal_waves(vec![selection(&[1, 2]); 2], nonlinear_config(), gate)
            .unwrap();
        assert_eq!(
            bank.system().mode_count(),
            old_count + 4 * (1 + resolution().intervals)
        );
        assert_eq!(
            bank.source().longitudinal_waves().unwrap().strings().len(),
            2
        );
        assert!(bank.source().cylinder_gas_loading().is_some());
        assert!(bank.source().material_plate().is_some());
        assert_eq!(
            &bank.system().acceleration_pressure_weights().unwrap()[..old_count],
            weights.as_slice()
        );
        assert!(
            bank.system().acceleration_pressure_weights().unwrap()[old_count..]
                .iter()
                .all(|w| *w == 0.0)
        );
        for (index, (mode, state)) in modes.iter().zip(&states).enumerate() {
            assert_eq!(bank.system().components()[index].modes(), mode);
            assert_eq!(bank.system().components()[index].states(), state);
        }
        let axial = bank.source().longitudinal_waves().unwrap().strings()[0].axial_component();
        let fluid = bank.source().longitudinal_waves().unwrap().strings()[0]
            .axial_fluid_component()
            .expect("selected gas must load retained axial motion");
        let zero = vec![0.0; bank.system().mode_count()];
        let mut driven_axial = 0.0_f64;
        let mut axial_gas_energy = 0.0_f64;
        let mut sympathetic = 0.0_f64;
        for _ in 0..256 {
            let frame = *bank.step_under_gate(&zero, active(), gate).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            let report = bank.longitudinal_observation(0).unwrap();
            assert!(report.minimum_tension_n > 0.0 && report.maximum_abs_dynamic_strain <= 0.01);
            driven_axial = driven_axial.max(
                bank.system().components()[axial].states()[1]
                    .displacement_m_sqrt_kg
                    .abs(),
            );
            axial_gas_energy = axial_gas_energy.max(
                bank.system().components()[fluid]
                    .states()
                    .iter()
                    .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2))
                    .sum(),
            );
            sympathetic = sympathetic.max(
                bank.system().components()[1].states()[0]
                    .displacement_m_sqrt_kg
                    .abs(),
            );
        }
        assert!(
            driven_axial > 1e-18,
            "nonlinear transverse strain must drive the retained axial harmonic"
        );
        assert!(
            axial_gas_energy > 0.0,
            "bowed finite strain must also drive axial gas history"
        );
        assert!(
            sympathetic > 1e-15,
            "unplayed strings must stay coupled through the body"
        );
        let before = bank
            .system()
            .nonlinear_frame()
            .unwrap()
            .friction
            .network_energy_j;
        for _ in 0..64 {
            bank.step_under_gate(&zero, released(), gate).unwrap();
        }
        assert!(
            bank.system()
                .nonlinear_frame()
                .unwrap()
                .friction
                .network_energy_j
                <= before + 1e-12
        );
    });
}

#[test]
fn elastic_material_scales_axial_frequencies_while_gas_keeps_the_solid_axial_inertia() {
    with_cx(|cx, gate| {
        let build = |young, medium| {
            physical_bank(cx, gate, young, medium, 1)
                .with_longitudinal_waves(vec![selection(&[1, 2])], nonlinear_config(), gate)
                .unwrap()
        };
        let a = build(200e9, gas());
        let softer = build(50e9, gas());
        let denser_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap();
        let b = build(200e9, denser_gas);
        let law = |bank: &ReciprocalStringBank| {
            bank.source().longitudinal_waves().unwrap().strings()[0].clone()
        };
        let (la, ls, lb) = (law(&a), law(&softer), law(&b));
        let expected = core::f64::consts::PI / 0.5 * (200e9_f64 / 7800.0).sqrt();
        assert!((la.model().axial_omegas()[0] / expected - 1.0).abs() < 1e-14);
        assert!((ls.model().axial_omegas()[0] / la.model().axial_omegas()[0] - 0.5).abs() < 1e-14);
        assert_eq!(la.model().axial_omegas(), lb.model().axial_omegas());
        assert_eq!(
            la.solid_linear_density_kg_m(),
            lb.solid_linear_density_kg_m()
        );
        assert!(lb.transverse_linear_density_kg_m() > la.transverse_linear_density_kg_m());
        assert!(la.transverse_linear_density_kg_m() > la.solid_linear_density_kg_m());
    });
}

#[test]
fn axial_linear_limit_is_the_existing_exact_oscillator_without_double_stiffness() {
    with_cx(|cx, gate| {
        let mut bank = bank_without_causal_gas(cx, gate, 200e9, gas(), 1)
            .with_longitudinal_waves(vec![selection(&[1, 2])], nonlinear_config(), gate)
            .unwrap();
        let component = bank.source().longitudinal_waves().unwrap().strings()[0].axial_component();
        assert!(
            bank.source().longitudinal_waves().unwrap().strings()[0]
                .axial_fluid_component()
                .is_none()
        );
        let mut exact = bank.system().components()[component].clone();
        let offset: usize = bank.system().components()[..component]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let mut forces = vec![0.0; bank.system().mode_count()];
        for step in 0..96 {
            let load = if step == 0 { 1.0 } else { 0.0 };
            forces[offset] = load;
            exact.step(&[load, 0.0]).unwrap();
            let frame = *bank.step_under_gate(&forces, released(), gate).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            for (actual, expected) in bank.system().components()[component]
                .states()
                .iter()
                .zip(exact.states())
            {
                assert!(
                    (actual.displacement_m_sqrt_kg - expected.displacement_m_sqrt_kg).abs() < 1e-18
                );
                assert!(
                    (actual.velocity_m_sqrt_kg_per_s - expected.velocity_m_sqrt_kg_per_s).abs()
                        < 1e-14
                );
            }
            assert!(
                bank.system().components()[0]
                    .states()
                    .iter()
                    .all(|s| s.displacement_m_sqrt_kg == 0.0)
            );
        }
    });
}

#[test]
fn axial_bank_matches_independent_physical_force_and_momentum_coordinates() {
    with_cx(|cx, gate| {
        let mut bank = physical_bank(cx, gate, 200e9, gas(), 1)
            .with_longitudinal_waves(vec![selection(&[1])], nonlinear_config(), gate)
            .unwrap();
        let law = &bank.source().longitudinal_waves().unwrap().strings()[0];
        let component = law.axial_component();
        let fluid_component = law.axial_fluid_component().unwrap();
        let omega = law.model().axial_omegas()[0];
        // Integral sin²(k*x) dx = L/2. Compare the actual bank with a native
        // pHS in physical amplitude/momentum coordinates, with independently
        // derived effective length, dry mass, spring and work-conjugate force.
        let length = 0.5 * bank.cards()[0].length_m;
        let mass = length * law.solid_linear_density_kg_m();
        let root_mass = mass.sqrt();
        let flow = bank.source().cylinder_gas_loading().unwrap().strings()[0]
            .flow()
            .clone();
        let body = flow
            .axial_body(
                length,
                mass,
                mass * omega * omega,
                resolution().intervals + 2,
                gate,
            )
            .unwrap();
        assert_eq!(body.moving_mass_kg(), mass);
        let offset: usize = bank.system().components()[..component]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let mut forces = vec![0.0; bank.system().mode_count()];
        let mut state = vec![0.0; resolution().intervals + 2];
        for sample in 0..16 {
            let load = if sample == 0 { 1.0 } else { 0.0 };
            forces[offset] = load;
            bank.step_under_gate(&forces, released(), gate).unwrap();
            // The reference uses the distinct native pHS midpoint solver at
            // 128 substeps; its O(dt²) phase error bounds this comparison.
            for _ in 0..128 {
                state = fs_phs::step(
                    body.system(),
                    &state,
                    &[root_mass * load],
                    1.0 / (48_000.0 * 128.0),
                )
                .unwrap()
                .x;
            }
            let actual = bank.system().components()[component].states()[0];
            let q = root_mass * state[0];
            let v = state[1] / root_mass;
            let scale = q.hypot(v / omega);
            let error = (actual.displacement_m_sqrt_kg - q)
                .hypot((actual.velocity_m_sqrt_kg_per_s - v) / omega);
            assert!(
                error < 4e-5 * scale,
                "sample={sample}, physical axial error={error:e}/{scale:e}"
            );
            for (index, actual) in bank.system().components()[fluid_component]
                .states()
                .iter()
                .enumerate()
            {
                let expected = state[index + 2] / body.memory_masses_kg()[index].sqrt();
                assert!((actual.velocity_m_sqrt_kg_per_s - expected).abs() < 4e-5 * scale * omega);
            }
        }
    });
}

#[test]
fn axial_gas_viscosity_changes_retained_motion_and_closes_the_original_work_ledger() {
    with_cx(|cx, gate| {
        let mut specification = GasSpec::dry_air_ussa1976();
        specification.sutherland_beta *= 4.0;
        let viscous = GasState::try_new(&specification, 293.15, 101_325.0).unwrap();
        assert_eq!(viscous.density, gas().density);
        let build = |medium| {
            physical_bank(cx, gate, 200e9, medium, 1)
                .with_longitudinal_waves(vec![selection(&[1])], nonlinear_config(), gate)
                .unwrap()
        };
        let (mut a, mut b) = (build(gas()), build(viscous));
        let law = &a.source().longitudinal_waves().unwrap().strings()[0];
        let component = law.axial_component();
        let fluid_component = law.axial_fluid_component().unwrap();
        assert_eq!(
            law.model().axial_omegas(),
            b.source().longitudinal_waves().unwrap().strings()[0]
                .model()
                .axial_omegas()
        );
        let offset: usize = a.system().components()[..component]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let mut forces = vec![0.0; a.system().mode_count()];
        let (mut work, mut loss, mut memory_energy, mut difference) = (0.0, 0.0, 0.0_f64, 0.0_f64);
        for sample in 0..192 {
            forces[offset] = if sample == 0 { 1.0 } else { 0.0 };
            let fa = *a.step_under_gate(&forces, released(), gate).unwrap();
            let fb = *b.step_under_gate(&forces, released(), gate).unwrap();
            work += fa.external_work_j;
            loss += fa.network_dissipation_j;
            assert_eq!(fa.drive_work_j, 0.0);
            assert_eq!(fa.friction_dissipation_j, 0.0);
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
            memory_energy = memory_energy.max(
                a.system().components()[fluid_component]
                    .states()
                    .iter()
                    .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2))
                    .sum(),
            );
            difference = difference.max(
                (a.system().components()[component].states()[0].velocity_m_sqrt_kg_per_s
                    - b.system().components()[component].states()[0].velocity_m_sqrt_kg_per_s)
                    .abs(),
            );
        }
        assert!(
            memory_energy > 0.0,
            "axial gas must receive real kinetic energy"
        );
        assert!(
            difference > 1e-9,
            "source viscosity must affect the actual axial ringdown"
        );
        assert!(loss > 0.0);
        assert!((a.system().total_energy_j().unwrap() + loss - work).abs() < 2e-12);
    });
}

#[test]
fn complete_axial_memory_is_admitted_before_any_modes_can_be_silently_omitted() {
    with_cx(|cx, gate| {
        let admitted = physical_bank(cx, gate, 200e9, gas(), 2)
            .with_longitudinal_waves(vec![selection(&[1, 2, 3]); 2], nonlinear_config(), gate)
            .unwrap();
        assert_eq!(admitted.system().mode_count(), 62);
        // All four axial frequencies fit the clock, but the full gas history
        // would need 72 coordinates, exceeding the caller's original cap 64.
        assert!(
            physical_bank(cx, gate, 200e9, gas(), 2)
                .with_longitudinal_waves(
                    vec![selection(&[1, 2, 3, 4]); 2],
                    nonlinear_config(),
                    gate
                )
                .is_err()
        );
    });
}

#[test]
fn invalid_selection_clock_and_duplicate_strain_laws_refuse_without_changing_physics() {
    with_cx(|cx, gate| {
        for bad in [vec![], vec![0], vec![2, 1], vec![1, 1], vec![1025]] {
            assert!(selection(&bad).validate().is_err());
        }
        // Retained mode five exceeds this original 48 kHz clock. It must not be
        // silently deleted or retuned even though its material band is valid.
        assert!(
            physical_bank(cx, gate, 200e9, gas(), 1)
                .with_longitudinal_waves(vec![selection(&[5])], nonlinear_config(), gate)
                .is_err()
        );
        let bank = physical_bank(cx, gate, 200e9, gas(), 1)
            .with_geometric_stretching(0.2, nonlinear_config(), gate)
            .unwrap();
        assert!(
            bank.with_longitudinal_waves(vec![selection(&[2])], nonlinear_config(), gate)
                .is_err()
        );
        let bank = physical_bank(cx, gate, 200e9, gas(), 1)
            .with_longitudinal_waves(vec![selection(&[2])], nonlinear_config(), gate)
            .unwrap();
        assert!(
            bank.with_geometric_stretching(0.2, nonlinear_config(), gate)
                .is_err()
        );
    });
}
