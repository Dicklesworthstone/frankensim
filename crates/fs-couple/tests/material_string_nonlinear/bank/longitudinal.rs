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
        assert_eq!(bank.system().mode_count(), old_count + 4);
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
        let zero = vec![0.0; bank.system().mode_count()];
        let mut driven_axial = 0.0_f64;
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
fn elastic_material_scales_axial_frequencies_while_gas_changes_only_transverse_loading() {
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
        let mut bank = physical_bank(cx, gate, 200e9, gas(), 1)
            .with_longitudinal_waves(vec![selection(&[1, 2])], nonlinear_config(), gate)
            .unwrap();
        let component = bank.source().longitudinal_waves().unwrap().strings()[0].axial_component();
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
