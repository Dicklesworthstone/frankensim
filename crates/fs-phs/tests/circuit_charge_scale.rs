//! Charge continuity must remain an equation when capacitor charges are pC.

use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};

#[test]
fn picofarad_rc_does_not_accept_a_frozen_charge_at_voltage_scale_tolerance() {
    for capacitance in [1.0_f64, 1e-9, 1e-12, 1e-15] {
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Resistor {
                        ohms: capacitance.recip(),
                    },
                ),
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: capacitance,
                    },
                ),
            ],
            transformers: vec![],
        };
        let circuit = assemble_circuit(&graph).unwrap();
        // This already-consistent physical state avoids an unrelated algebraic
        // initialization transient masking the missing charge convergence test.
        let before = vec![capacitance, 1.0, -capacitance];
        assert_eq!(
            circuit.consistent_initial_state(&before, &[]).unwrap(),
            before
        );
        let dt = 0.01;
        let (step, defect) = circuit.step_audited(&before, &[], dt).unwrap();
        let expected = capacitance * (1.0 - dt / 2.0) / (1.0 + dt / 2.0);
        assert!(
            (step.x[0] / expected - 1.0).abs() < 2e-12,
            "C={capacitance:e}: q1={:e}, expected={expected:e}",
            step.x[0]
        );
        assert!(
            defect / capacitance < 2e-13,
            "C={capacitance:e}, defect={defect:e}"
        );
    }
}

#[test]
fn a_farad_scale_branch_cannot_hide_a_picocoulomb_branch_in_the_same_dae() {
    let small = 1e-12_f64;
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: small }),
            (2, 0, Branch::Capacitor { farads: 1.0 }),
            (
                1,
                0,
                Branch::Resistor {
                    ohms: small.recip(),
                },
            ),
            (2, 0, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    let circuit = assemble_circuit(&graph).unwrap();
    let before = vec![small, 1.0, 1.0, 1.0, -small, -1.0];
    let (step, _) = circuit.step_audited(&before, &[], 0.01).unwrap();
    let expected = 0.995 / 1.005;
    assert!((step.x[0] / small - expected).abs() < 2e-12);
    assert!((step.x[1] - expected).abs() < 2e-12);
}

#[test]
fn a_femtoampere_source_is_not_discarded_at_a_volt_biased_capacitor_node() {
    let capacitance = 1e-12;
    let current = 1e-15;
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (
                1,
                0,
                Branch::Capacitor {
                    farads: capacitance,
                },
            ),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let circuit = assemble_circuit(&graph).unwrap();
    // Capacitor/source were at equilibrium before the source switched on.
    let before = vec![capacitance, 1.0, 0.0];
    let (step, defect) = circuit.step_audited(&before, &[current], 1.0).unwrap();
    assert!((step.x[0] / (capacitance + current) - 1.0).abs() < 2e-13);
    let work = current * (1.0 + 0.5 * current / capacitance);
    assert!((step.supplied / work - 1.0).abs() < 2e-13);
    assert!(defect < 1e-25);
}
