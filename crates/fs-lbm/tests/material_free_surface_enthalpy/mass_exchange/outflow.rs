use super::*;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceMassOutflow;

fn outflow(datum: f64) -> FreeSurfaceMassOutflow {
    FreeSurfaceMassOutflow {
        exchange: FreeSurfaceMassExchange {
            coordinates: INTERFACE,
            mass_kg: -1e-7,
            specific_enthalpy_j_kg: datum + 200_000.0,
            heat_j: 0.0,
        },
        velocity_m_s: [-0.04, 0.05, -0.02],
    }
}

#[test]
fn g1_powered_outflow_drains_actual_liquid_and_recoils_without_inventing_heat() {
    for datum in [0.0, 1e9] {
        with_liquid(datum, |liquid, gate| {
            let input = outflow(datum);
            let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            let mass = before.density_kg_m3 * before.fill_fraction * 1e-9;
            let energy = liquid.total_enthalpy_j().unwrap();
            let h = liquid.specific_enthalpy_at(INTERFACE).unwrap();
            let report = liquid
                .eject_interface_outflow(&[input], 1e-8, controls(), gate)
                .unwrap();
            let after = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            let dm = report.exchange.mass_change_kg;
            let remaining = mass + dm;
            let expected_work = -0.5 * dm * mass / remaining
                * (0..3)
                    .map(|d| (input.velocity_m_s[d] - before.velocity_m_s[d]).powi(2))
                    .sum::<f64>();
            let kinetic = 0.5 * remaining * after.velocity_m_s.iter().map(|u| u * u).sum::<f64>()
                - 0.5 * mass * before.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
            let carried = 0.5 * dm * input.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
            assert_eq!(before.density_kg_m3, after.density_kg_m3);
            assert_eq!(
                before.pressure_perturbation_pa,
                after.pressure_perturbation_pa
            );
            near(after.fill_fraction, 0.3, 2e-14);
            for d in 0..3 {
                near(
                    remaining * after.velocity_m_s[d],
                    mass * before.velocity_m_s[d] + dm * input.velocity_m_s[d],
                    2e-13,
                );
                near(
                    report.exchange.carried_momentum_kg_m_s[d],
                    dm * input.velocity_m_s[d],
                    2e-14,
                );
            }
            near(report.mechanical_work_j, expected_work, 2e-13);
            near(report.kinetic_energy_change_j, kinetic, 2e-13);
            near(report.carried_kinetic_energy_j, carried, 2e-13);
            assert_eq!(liquid.specific_enthalpy_at(INTERFACE).unwrap(), h);
            assert_eq!(report.exchange.heat_j, 0.0);
            near(report.exchange.enthalpy_change_j, dm * h, 2e-13);
            assert!((kinetic - carried - report.mechanical_work_j).abs() < 1e-20);
            assert!(report.absolute_total_energy_defect_j < controls().energy_tolerance_j);
            assert!(
                report.exchange.absolute_momentum_defect_kg_m_s
                    < controls().momentum_tolerance_kg_m_s
            );
            assert_eq!(liquid.physical_time_s(), 0.0);
            // The resulting volume, velocity and original h advance through
            // the actual native thermal/flow owner on its retained timeline.
            let advanced = liquid.step(gate).unwrap();
            near(
                advanced.enthalpy_before_j,
                energy + report.exchange.enthalpy_change_j,
                2e-12,
            );
            assert_eq!(liquid.physical_time_s(), 0.001);
        });
    }
}

#[test]
fn g3_comoving_outflow_preserves_original_signed_exchange_exactly_without_work() {
    let original = with_liquid(0.0, |liquid, gate| {
        let report = liquid
            .exchange_interface_mass(&[outflow(0.0).exchange], controls(), gate)
            .unwrap();
        (
            report,
            liquid.flow().lattice().grid().state_bits(),
            liquid.specific_enthalpy_at(INTERFACE),
            liquid.material_at(INTERFACE).unwrap().clone(),
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
            .eject_interface_outflow(
                &[FreeSurfaceMassOutflow {
                    velocity_m_s: velocity,
                    ..outflow(0.0)
                }],
                0.0,
                controls(),
                gate,
            )
            .unwrap();
        assert_eq!(report.exchange, original.0);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original.1);
        assert_eq!(liquid.specific_enthalpy_at(INTERFACE), original.2);
        assert_eq!(liquid.material_at(INTERFACE).unwrap(), &original.3);
        assert_eq!(report.mechanical_work_j, 0.0);
    });
}

#[test]
fn g4_outflow_work_speed_capacity_and_source_refusals_preserve_all_accepted_state() {
    with_liquid(0.0, |liquid, gate| {
        let bits = liquid.flow().lattice().grid().state_bits();
        let sample = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let energy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(INTERFACE).unwrap().clone();
        let valid = outflow(0.0);
        for maximum_work in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                liquid
                    .eject_interface_outflow(&[valid], maximum_work, controls(), gate)
                    .is_err()
            );
        }
        let mut cases = vec![valid; 8];
        cases[0].exchange.mass_kg = 1e-7;
        cases[1].exchange.mass_kg = -1e-6;
        cases[2].exchange.mass_kg = -1e-30;
        cases[3].exchange.heat_j = 1.0;
        // An admitted outlet speed cannot hide an excessive remaining recoil
        // as the actual native liquid inventory approaches exhaustion.
        cases[4].exchange.mass_kg = -3.999e-7;
        // Tiny mass cannot bypass the outlet's own velocity domain.
        cases[5].exchange.mass_kg = -1e-12;
        cases[5].velocity_m_s = [1e6, 0.0, 0.0];
        cases[6].exchange.specific_enthalpy_j_kg = -1.0;
        cases[7].exchange.specific_enthalpy_j_kg = 600_001.0;
        for input in cases {
            assert!(
                liquid
                    .eject_interface_outflow(&[input], 1.0, controls(), gate)
                    .is_err()
            );
        }
        assert!(
            liquid
                .eject_interface_outflow(&[valid, valid], 1.0, controls(), gate)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .eject_interface_outflow(&[valid], 1.0, controls(), &cancelled)
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
