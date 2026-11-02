//! The omitted string tail is not a thermal material or a heat source.
use super::*;

fn fiber_instrument(
    cx: &Cx<'_>,
    gate: &CancelGate,
    card: &MaterialCard,
    young: f64,
) -> MultiBowStringBank {
    base(cx, gate, 7800.0)
        .with_residual_flexibility(gate)
        .unwrap()
        .with_tensioned_fiber_ribbon(
            super::super::fiber_ribbon::fibers(3, 2, 1300.0, young, 50.0, 40.0, 48000, gate),
            0.4,
            ribbon_limits(),
            gate,
        )
        .unwrap()
        .with_bow_shear_layers(vec![(0, layer(card))], gate)
        .unwrap()
        .with_bow_residual_flexibility(gate)
        .unwrap()
}

#[test]
fn receiving_residual_fibers_and_finite_film_heat_change_the_same_pressure() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating_law_with_capacity(true, None, None, CP);
        let (other_card, other_caloric) = coating_law_with_capacity(true, None, None, 2.0 * CP);
        let mut bank = fiber_instrument(cx, gate, &card, 3e9);
        let mut stiffer = fiber_instrument(cx, gate, &card, 12e9);
        let mut larger_capacity = fiber_instrument(cx, gate, &other_card, 3e9);
        let original_receiver = bank.system().port(0).unwrap().right.clone();
        let first_fiber = bank.source().fiber_ribbon().unwrap().component_start();
        // The tail belongs to the string: changing fiber material cannot
        // secretly retune that tail or the original string/body spectrum.
        assert_eq!(
            bank.system().friction_residual_compliance(),
            stiffer.system().friction_residual_compliance()
        );
        for (a, b) in bank.system().components()[..first_fiber]
            .iter()
            .zip(&stiffer.system().components()[..first_fiber])
        {
            assert_eq!(a.modes(), b.modes());
        }
        let mut heat = thermal(caloric.curve(), false);
        let mut stiff_heat = thermal(caloric.curve(), false);
        let mut larger_heat = thermal(other_caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let larger_bound = binding(&other_card, other_caloric.curve());
        let forces = vec![0.0; bank.system().mode_count()];
        let (mut work, mut physical_heat, mut outside_loss, mut numerical) = (0.0, 0.0, 0.0, 0.0);
        let (mut fiber_pressure_change, mut caloric_pressure_change) = (0.0, 0.0);
        for sample in 0..128 {
            let drives = [
                drive(if sample < 64 { 0.001 } else { -0.001 }, 0.1),
                drive(0.0, 0.0),
                drive(0.0, 0.0),
            ];
            if sample == 8 {
                let before = states(&bank);
                let cells = heat.cells().clone();
                let history = bank.system().material_shear_state(0).unwrap().clone();
                let accepted = bank.system().last_frame().unwrap().clone();
                let mut limited = controls(1);
                limited.maximum_mechanical_trials = 1;
                assert!(
                    bank.step_with_layer_thermal_network(
                        &mut heat,
                        &bound,
                        limited,
                        &[0.0],
                        &forces,
                        &drives,
                        gate,
                        other
                    )
                    .is_err()
                );
                assert_eq!(states(&bank), before);
                assert_eq!(heat.cells(), &cells);
                assert_eq!(bank.system().material_shear_state(0), Some(&history));
                assert_eq!(bank.system().last_frame(), Some(&accepted));
            }
            let report = bank
                .step_with_layer_thermal_network(
                    &mut heat,
                    &bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let stiff_report = stiffer
                .step_with_layer_thermal_network(
                    &mut stiff_heat,
                    &bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let caloric_report = larger_capacity
                .step_with_layer_thermal_network(
                    &mut larger_heat,
                    &larger_bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let frame = &report.mechanics;
            work += frame.drive_work_j;
            physical_heat += report.deposited_heat_j;
            numerical += frame.friction_residual_loss_j;
            outside_loss += frame.network_dissipation_j + frame.friction_residual_loss_j;
            assert!(
                (report.deposited_heat_j - frame.ports[0].shear_layer.unwrap().dissipated_j).abs()
                    < 1e-15
            );
            assert!((total_heat(&heat) - physical_heat).abs() < 1e-11);
            assert!(
                (bank.system().total_energy_j().unwrap() + total_heat(&heat) + outside_loss - work)
                    .abs()
                    < 2e-10
            );
            assert!(report.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            assert_eq!(bank.system().port(0).unwrap().right, original_receiver);
            fiber_pressure_change +=
                (frame.observer_pressure_pa - stiff_report.mechanics.observer_pressure_pa).abs();
            caloric_pressure_change +=
                (frame.observer_pressure_pa - caloric_report.mechanics.observer_pressure_pa).abs();
        }
        assert!(
            bank.system().components()[first_fiber]
                .states()
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        );
        assert!(physical_heat > 0.0 && numerical > 0.0);
        assert!(fiber_pressure_change > 1e-15 && caloric_pressure_change > 1e-15);
        assert_eq!(
            bank.bow_shear_layer(0).unwrap().temperature_kelvin(),
            heat.temperature_k(FILM).unwrap()
        );
        assert_ne!(
            heat.temperature_k(FILM).unwrap(),
            larger_heat.temperature_k(FILM).unwrap()
        );
    });
}

fn instrument(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowStringBank {
    base(cx, gate, 7800.0)
        .with_residual_flexibility(gate)
        .unwrap()
        .with_additional_bows(vec![extra()], joint_limits(), gate)
        .unwrap()
        .with_bow_residual_flexibility(gate)
        .unwrap()
        .with_bow_shear_layers(vec![(0, layer(card))], gate)
        .unwrap()
}

#[test]
fn residual_film_heat_changes_pressure_without_heating_numerical_tail_loss() {
    with_cx(|cx, gate| {
        // Only finite caloric capacity differs. Both films use the same
        // manufactured temperature-dependent yield/viscosity, not a sound fit.
        let (card, caloric) = coating_law_with_capacity(true, None, None, CP);
        let (other_card, other_caloric) = coating_law_with_capacity(true, None, None, 2.0 * CP);
        let mut bank = instrument(cx, gate, &card);
        let mut changed = instrument(cx, gate, &other_card);
        let mut heat = thermal(caloric.curve(), false);
        let mut other_heat = thermal(other_caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let other_bound = binding(&other_card, other_caloric.curve());
        let forces = vec![0.0; bank.system().mode_count()];
        let (mut work, mut physical_heat, mut outside_loss, mut tail_loss, mut pressure_difference) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for sample in 0..256 {
            let drives = [
                drive(if sample < 128 { 0.001 } else { -0.001 }, 0.1),
                drive(0.0, 0.0),
            ];
            if sample == 8 {
                // Refuse after an actual mechanical/thermal candidate has been
                // evaluated. Neither accepted film nor omitted-force history
                // may be replaced by a failed Newton trial.
                let old_states = states(&bank);
                let old_frame = bank.system().last_frame().unwrap().clone();
                let old_history = bank.system().material_shear_state(0).unwrap().clone();
                let old_source = bank.bow_shear_layer(0).unwrap().clone();
                let old_energy = bank.system().total_energy_j().unwrap();
                let old_cells = heat.cells().clone();
                let old_time = heat.time_s();
                let mut limited = controls(1);
                limited.maximum_mechanical_trials = 1;
                assert!(
                    bank.step_with_layer_thermal_network(
                        &mut heat,
                        &bound,
                        limited,
                        &[0.0],
                        &forces,
                        &drives,
                        gate,
                        other
                    )
                    .is_err()
                );
                assert_eq!(states(&bank), old_states);
                assert_eq!(bank.system().last_frame(), Some(&old_frame));
                assert_eq!(bank.system().material_shear_state(0), Some(&old_history));
                assert_eq!(bank.bow_shear_layer(0), Some(&old_source));
                assert_eq!(bank.system().total_energy_j().unwrap(), old_energy);
                assert_eq!(heat.cells(), &old_cells);
                assert_eq!(heat.time_s(), old_time);
            }
            let report = bank
                .step_with_layer_thermal_network(
                    &mut heat,
                    &bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let compared = changed
                .step_with_layer_thermal_network(
                    &mut other_heat,
                    &other_bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let frame = &report.mechanics;
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            assert_eq!(report.routed_plastic_heat_j, frame.ports[0].dissipation_j);
            assert_eq!(report.routed_relaxation_heat_j, 0.0);
            assert!((report.deposited_heat_j - report.routed_plastic_heat_j).abs() < 1e-15);
            work += frame.drive_work_j;
            physical_heat += report.routed_plastic_heat_j;
            tail_loss += frame.friction_residual_loss_j;
            outside_loss += frame.network_dissipation_j + frame.friction_residual_loss_j;
            // The film's plastic heat changes its finite caloric state. Network
            // damping and the backward-Euler tail defect remain outside it.
            assert!((total_heat(&heat) - physical_heat).abs() < 1e-11);
            assert!(
                (bank.system().total_energy_j().unwrap() + total_heat(&heat) + outside_loss - work)
                    .abs()
                    < 2e-10
            );
            assert_eq!(
                bank.bow_shear_layer(0).unwrap().temperature_kelvin(),
                heat.temperature_k(FILM).unwrap()
            );
            pressure_difference +=
                (frame.observer_pressure_pa - compared.mechanics.observer_pressure_pa).abs();
        }
        assert!(physical_heat > 0.0 && tail_loss > 1e-14);
        assert!(heat.temperature_k(FILM).unwrap() > 300.0);
        assert!(heat.temperature_k(FILM).unwrap() != other_heat.temperature_k(FILM).unwrap());
        assert!(
            pressure_difference > 1e-15,
            "finite film heat must reach the original pressure solve"
        );
        assert_eq!(bank.system().samples_rendered(), 256);
    });
}
