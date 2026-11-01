//! Independent physical modes of a coupled capacitor block at real SI scales.
//! A large stationary store must not hide either small current or small work.
use fs_phs::circuit::capacitance::{CapacitanceBlock, assemble_circuit_with_capacitance_blocks};
use fs_phs::circuit::{Branch, CircuitGraph};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn coupled_small_capacitance_modes_retain_distinct_decay_and_work_beside_large_bias() {
    // C=c[[2,-1],[-1,2]], with equal R on both physical terminals. The common
    // and differential voltage modes have time constants R*c and 3*R*c.
    // For dt=.2ms and R*c=1ms, the native midpoint map has exact modal factors
    // 9/11 and 29/31. These formulas do not query the constitutive implementation.
    let dt = 2e-4;
    let common = 2.0;
    let differential = 3.0;
    let initial_voltage = [5.0, -1.0];
    let common_after = common * 9.0 / 11.0;
    let differential_after = differential * 29.0 / 31.0;
    let expected_voltage = [
        common_after + differential_after,
        common_after - differential_after,
    ];
    for capacitance in [1e-9, 1e-12, 1e-15] {
        let resistance = 1e-3 / capacitance;
        let expected_charge = [
            capacitance * (2.0 * expected_voltage[0] - expected_voltage[1]),
            capacitance * (2.0 * expected_voltage[1] - expected_voltage[0]),
        ];
        let initial_charge = [11.0 * capacitance, -7.0 * capacitance];
        let expected_energy_before = 31.0 * capacitance;
        let expected_energy_after = capacitance
            * (common_after * common_after + 3.0 * differential_after * differential_after);
        let expected_heat: [f64; 2] = core::array::from_fn(|i| {
            let midpoint = (initial_voltage[i] + expected_voltage[i]) / 2.0;
            dt * midpoint * midpoint / resistance
        });
        near(
            expected_energy_after - expected_energy_before,
            -expected_heat.iter().sum::<f64>(),
            2e-14 * capacitance,
        );
        for bias in [0.0, 1e8] {
            let block = CapacitanceBlock::try_new(
                vec![0, 2],
                vec![
                    2.0 * capacitance,
                    -capacitance,
                    -capacitance,
                    2.0 * capacitance,
                ],
            )
            .unwrap();
            let graph = CircuitGraph {
                node_count: 4,
                branches: vec![
                    (
                        1,
                        0,
                        Branch::Capacitor {
                            farads: 2.0 * capacitance,
                        },
                    ),
                    (1, 0, Branch::Resistor { ohms: resistance }),
                    (
                        2,
                        0,
                        Branch::Capacitor {
                            farads: 2.0 * capacitance,
                        },
                    ),
                    (2, 0, Branch::Resistor { ohms: resistance }),
                    (3, 0, Branch::Capacitor { farads: 1.0 }),
                    (3, 0, Branch::VoltageSource { port: 0 }),
                ],
                transformers: vec![],
            };
            let dae = assemble_circuit_with_capacitance_blocks(&graph, &[block.clone()], &[], &[])
                .unwrap();
            let input = [bias];
            let mut initial = vec![0.0; dae.system.state_dim()];
            for (index, charge) in initial_charge.into_iter().chain([bias]).enumerate() {
                initial[dae.charge_index[index]] = charge;
            }
            let initial = dae.consistent_initial_state(&initial, &input).unwrap();
            let initial_effort = dae.capacitor_voltages(&initial).unwrap();
            for i in 0..2 {
                near(initial_effort[i], initial_voltage[i], 2e-12);
            }
            near(
                block.stored_energy_j(&initial_charge).unwrap(),
                expected_energy_before,
                2e-13 * capacitance,
            );
            let (record, native_defect) = dae.step_audited(&initial, &input, dt).unwrap();
            let actual_effort = dae.capacitor_voltages(&record.x).unwrap();
            let actual_charge: [f64; 2] = core::array::from_fn(|i| record.x[dae.charge_index[i]]);
            let mut actual_heat = 0.0;
            for i in 0..2 {
                near(actual_charge[i], expected_charge[i], 2e-10 * capacitance);
                near(actual_effort[i], expected_voltage[i], 2e-10);
                let midpoint = (initial_effort[i] + actual_effort[i]) / 2.0;
                // KCL checks each physical terminal, not just total energy.
                near(
                    actual_charge[i] - initial_charge[i],
                    -dt * midpoint / resistance,
                    2e-10 * capacitance,
                );
                let heat = dt * midpoint * midpoint / resistance;
                near(heat, expected_heat[i], 2e-10 * capacitance);
                actual_heat += heat;
            }
            let delta_energy = expected_energy_after - expected_energy_before;
            near(
                block
                    .energy_change_j(&initial_charge, &actual_charge)
                    .unwrap(),
                delta_energy,
                2e-10 * capacitance,
            );
            near(record.delta_h, delta_energy, 2e-10 * capacitance);
            near(record.dissipated, actual_heat, 2e-10 * capacitance);
            near(record.supplied, 0.0, 2e-10 * capacitance);
            near(record.balance_residual(), 0.0, 2e-10 * capacitance);
            assert!(
                native_defect <= 2e-10 * capacitance,
                "native defect={native_defect:e}"
            );
            assert!(record.delta_h < 0.0 && record.dissipated > 0.0);
            if bias != 0.0 {
                assert_eq!(
                    record.x[dae.charge_index[2]].to_bits(),
                    initial[dae.charge_index[2]].to_bits(),
                );
                // Neither endpoint sum can resolve this actual capacitor work.
                assert_eq!(
                    dae.system.hamiltonian(&initial),
                    dae.system.hamiltonian(&record.x),
                );
            }
        }
    }
}
