//! G1/G3/G4: one simultaneous source transaction on the actual material owner.

use super::*;
use fs_lbm::material::free_surface::enthalpy::{FreeSurfaceMassInflow, FreeSurfaceMassOutflow};

fn streams(datum: f64) -> (FreeSurfaceMassInflow, FreeSurfaceMassOutflow) {
    (
        FreeSurfaceMassInflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: INTERFACE,
                mass_kg: 1e-7,
                specific_enthalpy_j_kg: datum + 250_000.0,
                heat_j: 0.0,
            },
            velocity_m_s: [-0.04, 0.05, -0.02],
        },
        FreeSurfaceMassOutflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: INTERFACE,
                mass_kg: -1e-7,
                specific_enthalpy_j_kg: datum + 150_000.0,
                heat_j: 0.0,
            },
            velocity_m_s: [0.04, -0.03, 0.01],
        },
    )
}

#[test]
fn g1_equal_throughput_changes_actual_momentum_and_enthalpy_without_changing_fill_or_clock() {
    for datum in [0.0, 1e9] {
        with_liquid(datum, |liquid, gate| {
            let (incoming, outgoing) = streams(datum);
            let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            let mass = before.density_kg_m3 * before.fill_fraction * 1e-9;
            let energy = liquid.total_enthalpy_j().unwrap();
            let initial_h = liquid.specific_enthalpy_at(INTERFACE).unwrap();
            let report = liquid
                .exchange_interface_throughflow(&[incoming], &[outgoing], 1e-8, controls(), gate)
                .unwrap();
            let after = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            assert_eq!(report.exchange.mass_change_kg, 0.0);
            assert_eq!(report.samples[0].incoming_mass_kg, 1e-7);
            assert_eq!(report.samples[0].outgoing_mass_kg, 1e-7);
            assert_eq!(before.fill_fraction, after.fill_fraction);
            assert_eq!(before.density_kg_m3, after.density_kg_m3);
            assert_eq!(
                before.pressure_perturbation_pa,
                after.pressure_perturbation_pa
            );
            let delta: [f64; 3] = core::array::from_fn(|k| {
                1e-7 * (incoming.velocity_m_s[k] - outgoing.velocity_m_s[k]) / mass
            });
            for k in 0..3 {
                assert!((after.velocity_m_s[k] - before.velocity_m_s[k] - delta[k]).abs() < 2e-14);
                assert!(
                    (report.exchange.carried_momentum_kg_m_s[k] - mass * delta[k]).abs() < 1e-21
                );
            }
            // At fixed mass u(s) is linear; independently integrate both
            // relative-speed quadratics rather than inferring one from balance.
            let loss = |velocity: [f64; 3]| {
                0.5e-7
                    * (0..3)
                        .map(|k| {
                            let relative = velocity[k] - before.velocity_m_s[k];
                            relative * relative - relative * delta[k] + delta[k] * delta[k] / 3.0
                        })
                        .sum::<f64>()
            };
            assert!((report.mixing_heat_j - loss(incoming.velocity_m_s)).abs() < 1e-21);
            assert!((report.mechanical_work_j - loss(outgoing.velocity_m_s)).abs() < 1e-21);
            assert!(report.mixing_heat_j > 0.0 && report.mechanical_work_j > 0.0);
            let carried_h = 1e-7
                * (incoming.exchange.specific_enthalpy_j_kg
                    - outgoing.exchange.specific_enthalpy_j_kg);
            near(report.exchange.carried_enthalpy_j, carried_h, 2e-13);
            near(
                liquid.specific_enthalpy_at(INTERFACE).unwrap(),
                initial_h + (carried_h + report.mixing_heat_j) / mass,
                2e-13,
            );
            assert!(report.absolute_total_energy_defect_j <= controls().energy_tolerance_j);
            assert!(
                report.exchange.absolute_momentum_defect_kg_m_s
                    <= controls().momentum_tolerance_kg_m_s
            );
            assert_eq!(report.exchange.physical_time_s, 0.0);
            assert_eq!(liquid.physical_time_s(), 0.0);
            let next = liquid.step(gate).unwrap();
            near(
                next.enthalpy_before_j,
                energy + report.exchange.enthalpy_change_j,
                2e-12,
            );
            assert_eq!(liquid.physical_time_s(), 0.001);
        });
    }
}

#[test]
fn g3_single_direction_preserves_complete_original_native_results() {
    for outgoing in [false, true] {
        let original = with_liquid(0.0, |liquid, gate| {
            let (input, output) = streams(0.0);
            let (exchange, carried, kinetic, heat, work, defect) = if outgoing {
                let report = liquid
                    .eject_interface_outflow(&[output], 1e-8, controls(), gate)
                    .unwrap();
                (
                    report.exchange,
                    report.carried_kinetic_energy_j,
                    report.kinetic_energy_change_j,
                    0.0,
                    report.mechanical_work_j,
                    report.absolute_total_energy_defect_j,
                )
            } else {
                let report = liquid
                    .mix_interface_inflow(&[input], controls(), gate)
                    .unwrap();
                (
                    report.exchange,
                    report.carried_kinetic_energy_j,
                    report.kinetic_energy_change_j,
                    report.mixing_heat_j,
                    0.0,
                    report.absolute_total_energy_defect_j,
                )
            };
            (
                exchange,
                carried,
                kinetic,
                heat,
                work,
                defect,
                liquid.flow().lattice().grid().state_bits(),
                liquid.specific_enthalpy_at(INTERFACE),
                liquid.material_at(INTERFACE).unwrap().clone(),
            )
        });
        with_liquid(0.0, |liquid, gate| {
            let (input, output) = streams(0.0);
            let inputs = if outgoing { vec![] } else { vec![input] };
            let outputs = if outgoing { vec![output] } else { vec![] };
            let report = liquid
                .exchange_interface_throughflow(&inputs, &outputs, 1e-8, controls(), gate)
                .unwrap();
            assert_eq!(report.exchange, original.0);
            assert_eq!(report.carried_kinetic_energy_j, original.1);
            assert_eq!(report.kinetic_energy_change_j, original.2);
            assert_eq!(report.mixing_heat_j, original.3);
            assert_eq!(report.mechanical_work_j, original.4);
            assert_eq!(report.absolute_total_energy_defect_j, original.5);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), original.6);
            assert_eq!(liquid.specific_enthalpy_at(INTERFACE), original.7);
            assert_eq!(liquid.material_at(INTERFACE).unwrap(), &original.8);
        });
    }
}

#[test]
fn g4_work_source_capacity_cancellation_and_dropped_proposal_preserve_all_accepted_state() {
    with_liquid(0.0, |liquid, gate| {
        let (input, output) = streams(0.0);
        let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let bits = liquid.flow().lattice().grid().state_bits();
        let energy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(INTERFACE).unwrap().clone();
        assert!(
            liquid
                .prepare_interface_throughflow(&[input], &[output], 0.0, controls(), gate)
                .is_err()
        );
        let mut bad = [output; 4];
        bad[0].exchange.mass_kg = -1e-6;
        bad[1].exchange.specific_enthalpy_j_kg = 600_001.0;
        bad[2].velocity_m_s = [1e6, 0.0, 0.0];
        bad[3].exchange.mass_kg = 1e-7;
        for invalid in bad {
            assert!(
                liquid
                    .prepare_interface_throughflow(&[input], &[invalid], 1.0, controls(), gate)
                    .is_err()
            );
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .prepare_interface_throughflow(&[input], &[output], 1e-8, controls(), &cancelled)
                .is_err()
        );
        let expected = {
            let proposal = liquid
                .prepare_interface_throughflow(&[input], &[output], 1e-8, controls(), gate)
                .unwrap();
            proposal.report().clone()
        };
        assert_eq!(
            liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap(),
            before
        );
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
        assert_eq!(liquid.material_at(INTERFACE).unwrap(), &material);
        assert_eq!(liquid.physical_time_s(), 0.0);
        let accepted = liquid
            .prepare_interface_throughflow(&[input], &[output], 1e-8, controls(), gate)
            .unwrap()
            .publish();
        assert_eq!(accepted, expected);
    });
}
