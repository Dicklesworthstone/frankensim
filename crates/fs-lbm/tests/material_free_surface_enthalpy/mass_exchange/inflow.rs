use super::*;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceMassInflow;

fn inflow(datum: f64) -> FreeSurfaceMassInflow {
    FreeSurfaceMassInflow {
        exchange: FreeSurfaceMassExchange {
            coordinates: INTERFACE,
            mass_kg: 1e-7,
            specific_enthalpy_j_kg: datum + 100_000.0,
            heat_j: 0.02,
        },
        velocity_m_s: [-0.04, 0.05, -0.02],
    }
}

#[test]
fn g1_well_mixed_inflow_fills_and_accelerates_actual_liquid_with_closed_total_energy() {
    with_liquid(0.0, |liquid, gate| {
        let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let input = inflow(0.0);
        let mass = before.density_kg_m3 * before.fill_fraction * 1e-9;
        let energy = liquid.total_enthalpy_j().unwrap();
        let report = liquid
            .mix_interface_inflow(&[input], controls(), gate)
            .unwrap();
        let after = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let dm = report.exchange.mass_change_kg;
        let next_mass = mass + dm;
        let expected_heat = 0.5 * mass * dm / next_mass
            * (0..3)
                .map(|d| (input.velocity_m_s[d] - before.velocity_m_s[d]).powi(2))
                .sum::<f64>();
        let old_kinetic = 0.5 * mass * before.velocity_m_s.iter().map(|v| v * v).sum::<f64>();
        let next_kinetic = 0.5 * next_mass * after.velocity_m_s.iter().map(|v| v * v).sum::<f64>();
        let carried_kinetic = 0.5 * dm * input.velocity_m_s.iter().map(|v| v * v).sum::<f64>();
        assert_eq!(before.density_kg_m3, after.density_kg_m3);
        assert_eq!(
            before.pressure_perturbation_pa,
            after.pressure_perturbation_pa
        );
        near(after.fill_fraction, 0.5, 2e-14);
        for d in 0..3 {
            near(
                next_mass * after.velocity_m_s[d],
                mass * before.velocity_m_s[d] + dm * input.velocity_m_s[d],
                2e-13,
            );
            near(
                report.exchange.carried_momentum_kg_m_s[d],
                dm * input.velocity_m_s[d],
                2e-14,
            );
        }
        near(report.mixing_heat_j, expected_heat, 2e-13);
        near(report.carried_kinetic_energy_j, carried_kinetic, 2e-13);
        near(
            report.kinetic_energy_change_j,
            next_kinetic - old_kinetic,
            2e-13,
        );
        near(
            liquid.specific_enthalpy_at(INTERFACE).unwrap(),
            220_000.0 + expected_heat / next_mass,
            2e-13,
        );
        assert!(
            (liquid.total_enthalpy_j().unwrap() - energy + next_kinetic
                - old_kinetic
                - dm * input.exchange.specific_enthalpy_j_kg
                - input.exchange.heat_j
                - carried_kinetic)
                .abs()
                < controls().energy_tolerance_j
        );
        assert_eq!(report.exchange.heat_j, input.exchange.heat_j);
        assert!(
            report.exchange.absolute_momentum_defect_kg_m_s < controls().momentum_tolerance_kg_m_s
        );
        assert!(report.absolute_total_energy_defect_j < controls().energy_tolerance_j);
        assert_eq!(liquid.physical_time_s(), 0.0);
        // The original native timeline consumes the resulting filled, moving,
        // heated state; inflow did not advance a separate clock.
        let advanced = liquid.step(gate).unwrap();
        near(
            advanced.enthalpy_before_j,
            energy + report.exchange.enthalpy_change_j,
            2e-12,
        );
        assert_eq!(liquid.physical_time_s(), 0.001);
    });
}

#[test]
fn g3_comoving_inflow_preserves_the_original_exchange_arithmetic() {
    let original = with_liquid(0.0, |liquid, gate| {
        let report = liquid
            .exchange_interface_mass(&[inflow(0.0).exchange], controls(), gate)
            .unwrap();
        (
            report,
            liquid.flow().lattice().grid().state_bits(),
            liquid.specific_enthalpy_at(INTERFACE),
        )
    });
    with_liquid(0.0, |liquid, gate| {
        let velocity = liquid
            .flow()
            .sample_at_si(3, 3, 3)
            .unwrap()
            .unwrap()
            .velocity_m_s;
        let report = liquid
            .mix_interface_inflow(
                &[FreeSurfaceMassInflow {
                    velocity_m_s: velocity,
                    ..inflow(0.0)
                }],
                controls(),
                gate,
            )
            .unwrap();
        assert_eq!(report.exchange, original.0);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original.1);
        assert_eq!(liquid.specific_enthalpy_at(INTERFACE), original.2);
        assert_eq!(report.mixing_heat_j, 0.0);
    });
}

#[test]
fn g4_inflow_refusals_do_not_publish_fill_momentum_heat_or_sources() {
    with_liquid(0.0, |liquid, gate| {
        let bits = liquid.flow().lattice().grid().state_bits();
        let sample = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let energy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(INTERFACE).unwrap().clone();
        let valid = inflow(0.0);
        let mut cases = vec![valid; 7];
        cases[0].exchange.mass_kg = -1e-8;
        cases[1].exchange.mass_kg = 1e-6;
        cases[2].exchange.mass_kg = 1e-30;
        cases[3].exchange.heat_j = 1.0;
        // Even a tiny inflow cannot bypass the original speed domain by
        // producing a slow mass-weighted endpoint.
        cases[4].exchange.mass_kg = 1e-12;
        cases[4].velocity_m_s = [1e6, 0.0, 0.0];
        // Both mixed endpoints would fit the caloric chart; the incoming
        // liquid itself must also have an admitted state on that chart.
        cases[5].exchange.specific_enthalpy_j_kg = -1.0;
        cases[6].exchange.specific_enthalpy_j_kg = 600_001.0;
        for input in cases {
            assert!(
                liquid
                    .mix_interface_inflow(&[input], controls(), gate)
                    .is_err()
            );
        }
        // Each source fits by itself; the shared cell cannot store their sum.
        let individually_storable = FreeSurfaceMassInflow {
            exchange: FreeSurfaceMassExchange {
                mass_kg: 4e-7,
                ..valid.exchange
            },
            ..valid
        };
        assert!(
            liquid
                .mix_interface_inflow(
                    &[individually_storable, individually_storable],
                    controls(),
                    gate,
                )
                .is_err()
        );
        // An average inside the chart cannot admit an invalid source stream.
        assert!(
            liquid
                .mix_interface_inflow(
                    &[
                        valid,
                        FreeSurfaceMassInflow {
                            exchange: FreeSurfaceMassExchange {
                                specific_enthalpy_j_kg: -1.0,
                                ..valid.exchange
                            },
                            ..valid
                        }
                    ],
                    controls(),
                    gate,
                )
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .mix_interface_inflow(&[valid], controls(), &cancelled)
                .is_err()
        );
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(
            liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap(),
            sample
        );
        assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
        assert_eq!(liquid.material_at(INTERFACE).unwrap(), &material);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g1_opposing_inflows_share_native_mass_projection_and_retain_all_mixing_heat() {
    with_liquid(0.0, |liquid, gate| {
        let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let mass = before.density_kg_m3 * before.fill_fraction * 1e-9;
        let energy_before = liquid.total_enthalpy_j().unwrap();
        let sources = [
            FreeSurfaceMassInflow {
                exchange: FreeSurfaceMassExchange {
                    mass_kg: 5e-8,
                    specific_enthalpy_j_kg: 100_000.0,
                    heat_j: 0.0,
                    ..inflow(0.0).exchange
                },
                velocity_m_s: [0.06, 0.0, 0.0],
            },
            FreeSurfaceMassInflow {
                exchange: FreeSurfaceMassExchange {
                    mass_kg: 5e-8,
                    specific_enthalpy_j_kg: 300_000.0,
                    heat_j: 0.0,
                    ..inflow(0.0).exchange
                },
                velocity_m_s: [-0.06, 0.0, 0.0],
            },
        ];
        let report = liquid
            .mix_interface_inflow(&sources, controls(), gate)
            .unwrap();
        let after = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let dm = report.exchange.mass_change_kg;
        let represented = dm / sources.iter().map(|s| s.exchange.mass_kg).sum::<f64>();
        let remaining = mass + dm;
        let incoming_kinetic: f64 = sources
            .iter()
            .map(|s| {
                0.5 * represented
                    * s.exchange.mass_kg
                    * s.velocity_m_s.iter().map(|u| u * u).sum::<f64>()
            })
            .sum();
        let incoming_enthalpy: f64 = sources
            .iter()
            .map(|s| represented * s.exchange.mass_kg * s.exchange.specific_enthalpy_j_kg)
            .sum();
        let old_kinetic = 0.5 * mass * before.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
        let new_kinetic = 0.5 * remaining * after.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
        for d in 0..3 {
            let carried: f64 = sources
                .iter()
                .map(|s| represented * s.exchange.mass_kg * s.velocity_m_s[d])
                .sum();
            near(
                remaining * after.velocity_m_s[d],
                mass * before.velocity_m_s[d] + carried,
                2e-13,
            );
            near(report.exchange.carried_momentum_kg_m_s[d], carried, 2e-14);
        }
        assert_eq!(report.exchange.samples.len(), 1);
        assert_eq!(report.exchange.heat_j, 0.0);
        assert_eq!(before.density_kg_m3, after.density_kg_m3);
        near(after.fill_fraction, 0.5, 2e-14);
        near(report.carried_kinetic_energy_j, incoming_kinetic, 1e-16);
        near(report.exchange.carried_enthalpy_j, incoming_enthalpy, 1e-15);
        near(
            report.mixing_heat_j,
            incoming_kinetic + old_kinetic - new_kinetic,
            1e-16,
        );
        assert!(report.mixing_heat_j >= 0.99 * incoming_kinetic);
        assert!(
            (liquid.total_enthalpy_j().unwrap() - energy_before + new_kinetic
                - old_kinetic
                - incoming_enthalpy
                - incoming_kinetic)
                .abs()
                < controls().energy_tolerance_j
        );
        assert!(report.absolute_total_energy_defect_j < controls().energy_tolerance_j);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
