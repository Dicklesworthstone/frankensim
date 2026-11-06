//! Native Kirchhoff/charge integration with retained dielectric switching.
//! Coefficients are scaled analytical fixtures, not measured ferroelectric fits.
#![cfg(feature = "experimental-dielectric-hysteresis")]

use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitSolve, MaterialResistor,
    dielectric_hysteretic::{
        DielectricCircuitInitial, DielectricCircuitSolve, DielectricHystereticCircuit,
    },
    hysteretic::HystereticCircuitError,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor},
    field_storage::dielectric::hysteresis::{
        DielectricHysteresisBranch, ScalarDielectricHysteresis, UniformHystereticCapacitor,
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: resistance,
                dims: ELECTRICAL_RESISTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic dielectric-circuit resistor".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "analytical circuit fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}

fn linear(capacitance: f64, limit: f64) -> UniformHystereticCapacitor {
    UniformHystereticCapacitor::new(
        ScalarDielectricHysteresis::new(capacitance, limit, vec![]).unwrap(),
        1.0,
        1.0,
    )
    .unwrap()
}

fn switching(area: f64, gap: f64) -> UniformHystereticCapacitor {
    UniformHystereticCapacitor::new(
        ScalarDielectricHysteresis::new(
            0.01,
            30.0,
            vec![DielectricHysteresisBranch {
                susceptibility_f_m: 0.2,
                pinning_field_v_m: 2.0,
                saturation_polarization_c_m2: 1.0,
            }],
        )
        .unwrap(),
        area,
        gap,
    )
    .unwrap()
}

fn initial(
    capacitor: &UniformHystereticCapacitor,
    voltage: f64,
    fluxes: Vec<f64>,
    charges: Vec<f64>,
    inputs: Vec<f64>,
) -> DielectricCircuitInitial {
    DielectricCircuitInitial {
        dielectric: capacitor
            .candidate_voltage(&capacitor.initial_state(), voltage)
            .unwrap()
            .material
            .state,
        fluxes_wb: fluxes,
        remaining_charges_c: charges,
        external_inputs: inputs,
        time_s: 0.0,
    }
}

fn controls() -> DielectricCircuitSolve {
    DielectricCircuitSolve {
        maximum_circuit_trials: 256,
        maximum_newton_iterations: 32,
        difference_step_c: 1e-6,
        terminal_current_tolerance_a: 1e-8,
        energy_tolerance_j: 1e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            linear_restart: 2,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn rc_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Capacitor { farads: 0.2 }),
        ],
        transformers: vec![],
    }
}

#[test]
fn zero_branch_limit_matches_native_rc_and_rlc_with_other_charge_and_flux() {
    for extra_storage in [false, true] {
        let capacitor = linear(0.2, 100.0);
        let (graph, bound, resistor_branch, fluxes, remaining, inputs) = if extra_storage {
            // Original branch order interleaves the other C, external port gap,
            // sourced R and L before the capacitor being replaced.
            (
                CircuitGraph {
                    node_count: 4,
                    branches: vec![
                        (3, 0, Branch::Capacitor { farads: 0.3 }),
                        (1, 0, Branch::VoltageSource { port: 2 }),
                        (1, 2, Branch::Resistor { ohms: 999.0 }),
                        (2, 3, Branch::Inductor { henries: 0.5 }),
                        (2, 0, Branch::Capacitor { farads: 0.2 }),
                    ],
                    transformers: vec![],
                },
                4,
                2,
                vec![0.02],
                vec![0.015],
                vec![0.0, 0.0, 1.0],
            )
        } else {
            (rc_graph(), 2, 1, vec![], vec![], vec![1.0])
        };
        let sources = vec![resistor(resistor_branch, 0.8)];
        let start = initial(
            &capacitor,
            0.1,
            fluxes.clone(),
            remaining.clone(),
            inputs.clone(),
        );
        let mut native_charges = remaining;
        native_charges.push(0.02);
        let mut native = MaterialCircuit::new(
            graph.clone(),
            sources.clone(),
            fluxes,
            native_charges,
            &inputs,
        )
        .unwrap();
        let mut history =
            DielectricHystereticCircuit::new(graph, bound, capacitor, sources.clone(), start)
                .unwrap();
        assert_eq!(
            history.native_circuit().fluxes_wb().len(),
            usize::from(extra_storage)
        );
        assert_eq!(
            history.native_circuit().charges_c().len(),
            usize::from(extra_storage)
        );
        for tick in 1..=20 {
            let a = native
                .step_candidate(sources.clone(), &inputs, 0.005, controls().circuit, || {
                    false
                })
                .unwrap();
            let b = history
                .step_candidate(
                    f64::from(tick) * 0.005,
                    sources.clone(),
                    &inputs,
                    controls(),
                    || false,
                )
                .unwrap();
            near(
                b.circuit().charge_c(),
                *a.circuit().charges_c().last().unwrap(),
                1e-10,
            );
            for (&actual, &expected) in b
                .circuit()
                .native_circuit()
                .fluxes_wb()
                .iter()
                .zip(a.circuit().fluxes_wb())
            {
                near(actual, expected, 1e-10);
            }
            for (&actual, &expected) in b
                .circuit()
                .native_circuit()
                .charges_c()
                .iter()
                .zip(a.circuit().charges_c())
            {
                near(actual, expected, 1e-10);
            }
            near(b.step().stored_after_j, a.step().stored_after_j, 1e-10);
            near(b.step().supplied_work_j, a.step().native.supplied, 1e-10);
            near(b.step().heat_j, a.step().heat_j, 1e-10);
            near(b.step().dielectric.switching_heat_j, 0.0, 0.0);
            near(b.step().energy_residual_j, 0.0, 1e-9);
            assert_eq!(b.step().circuit.resistors[0].material, sources[0]);
            native = a.into_parts().0;
            history = b.into_parts().0;
        }
    }
}

#[test]
fn unrelated_storage_bias_does_not_erase_native_transfer_work() {
    let capacitor = linear(0.2, 100.0);
    // The isolated fourth node stores a large, unchanged energy at just 1 V.
    // A small real energy transfer into the other C/L must retain its own
    // native delta_h even though subtracting total endpoint energies loses it.
    let graph = CircuitGraph {
        node_count: 5,
        branches: vec![
            (3, 0, Branch::Capacitor { farads: 0.3 }),
            (4, 0, Branch::Capacitor { farads: 1e16 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 3, Branch::Inductor { henries: 0.5 }),
            (2, 0, Branch::Capacitor { farads: 0.2 }),
        ],
        transformers: vec![],
    };
    let sources = vec![resistor(3, 0.8)];
    let start = initial(&capacitor, 0.1, vec![0.02], vec![0.015, 1e16], vec![1.0]);
    let native = MaterialCircuit::new(
        graph.clone(),
        sources.clone(),
        vec![0.02],
        vec![0.015, 1e16, 0.02],
        &[1.0],
    )
    .unwrap();
    let circuit =
        DielectricHystereticCircuit::new(graph, 5, capacitor, sources.clone(), start).unwrap();
    let expected = native
        .step_candidate(sources.clone(), &[1.0], 0.005, controls().circuit, || false)
        .unwrap();
    let actual = circuit
        .step_candidate(0.005, sources, &[1.0], controls(), || false)
        .unwrap();
    let record = actual.step();
    assert_eq!(
        record.circuit.stored_before_j,
        record.circuit.stored_after_j
    );
    assert!(record.circuit.native.delta_h > 1e-6);
    near(record.stored_change_j, expected.step().native.delta_h, 1e-9);
    near(
        record.supplied_work_j,
        expected.step().native.supplied,
        1e-9,
    );
    near(record.heat_j, expected.step().heat_j, 1e-9);
    near(record.energy_residual_j, 0.0, 1e-9);
    near(
        actual.circuit().charge_c(),
        expected.circuit().charges_c()[2],
        1e-10,
    );
    assert_eq!(actual.circuit().native_circuit().charges_c()[1], 1e16);
}

#[test]
fn current_driven_major_loop_closes_exact_switching_heat_and_external_work() {
    let capacitor = switching(0.3, 0.02);
    let volume = capacitor.volume_m3();
    let start = initial(&capacitor, 0.2, vec![], vec![], vec![0.0, 0.0, 0.0]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 999.0 }),
            (1, 0, Branch::CurrentSource { port: 2 }),
        ],
        transformers: vec![],
    };
    let mut circuit = DielectricHystereticCircuit::new(graph, 0, capacitor, vec![], start).unwrap();
    let top = circuit.dielectric_state().clone();
    let initial_energy = circuit.dielectric_stored_energy_j();
    let (mut heat, mut supplied) = (0.0, 0.0);
    // q=±0.33 C at E=±10 V/m. Each reversal is an explicit segment;
    // exact work integration resolves its pinning and saturation crossings.
    for (target, current) in [(1.0, -0.66), (2.0, 0.66)] {
        let next = circuit
            .step_candidate(target, vec![], &[0.0, 0.0, current], controls(), || false)
            .unwrap();
        near(next.step().terminal_current_a, current, 1e-12);
        for probe in next.step().current_start_mid_end_a {
            near(probe, current, 1e-12);
        }
        heat += next.step().dielectric.switching_heat_j;
        supplied += next.step().supplied_work_j;
        near(
            next.step().heat_j,
            next.step().dielectric.switching_heat_j,
            1e-14,
        );
        near(next.step().circuit.heat_j, 0.0, 0.0);
        near(next.step().numerical_dissipation_j, 0.0, 0.0);
        near(next.step().energy_residual_j, 0.0, 1e-10);
        circuit = next.into_parts().0;
    }
    near(heat, 4.0 * 2.0 * 1.0 * volume, 1e-10);
    near(supplied, heat, 1e-10);
    near(circuit.charge_c(), 0.33, 1e-10);
    near(circuit.voltage_v(), 0.2, 1e-10);
    near(circuit.dielectric_stored_energy_j(), initial_energy, 1e-10);
    assert_eq!(
        circuit.dielectric_state().polarizations_c_m2(),
        top.polarizations_c_m2()
    );
}

#[test]
fn active_switching_rc_segment_is_second_order_against_the_independent_ode() {
    // At E=3 the branch is switching and remains below saturation at E=7.
    // A=d=R=1: (.01+.2) dV/dt = 4-V, V(0)=3.
    let expected = 4.0 - (-0.2_f64 / 0.21).exp();
    let mut previous_error = None;
    for count in [10, 20, 40] {
        let capacitor = switching(1.0, 1.0);
        let start = initial(&capacitor, 3.0, vec![], vec![], vec![4.0]);
        let sources = vec![resistor(1, 1.0)];
        let mut circuit =
            DielectricHystereticCircuit::new(rc_graph(), 2, capacitor, sources.clone(), start)
                .unwrap();
        let mut switching_heat = 0.0;
        for tick in 1..=count {
            let next = circuit
                .step_candidate(
                    0.2 * f64::from(tick) / f64::from(count),
                    sources.clone(),
                    &[4.0],
                    controls(),
                    || false,
                )
                .unwrap();
            switching_heat += next.step().dielectric.switching_heat_j;
            circuit = next.into_parts().0;
        }
        let error = (circuit.voltage_v() - expected).abs();
        if let Some(previous) = previous_error {
            let ratio: f64 = previous / error;
            assert!((3.8..4.2).contains(&ratio), "second-order ratio {ratio}");
        }
        near(
            switching_heat,
            2.0 * 0.2 * (circuit.voltage_v() - 3.0),
            2e-10,
        );
        previous_error = Some(error);
    }
    assert!(previous_error.unwrap() < 2e-5);
}

#[test]
fn cancellation_trial_budget_and_field_domain_refusals_leave_history_untouched() {
    let capacitor = linear(0.2, 100.0);
    let start = initial(&capacitor, 0.0, vec![], vec![], vec![1.0]);
    let sources = vec![resistor(1, 1.0)];
    let circuit =
        DielectricHystereticCircuit::new(rc_graph(), 2, capacitor, sources.clone(), start).unwrap();
    let saved = circuit.clone();
    assert!(matches!(
        circuit.step_candidate(0.1, sources.clone(), &[1.0], controls(), || true),
        Err(HystereticCircuitError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        circuit.step_candidate(0.1, sources.clone(), &[1.0], controls(), || {
            polls += 1;
            polls >= 6
        }),
        Err(HystereticCircuitError::Cancelled)
    ));
    assert!(polls >= 6);
    let mut bounded = controls();
    bounded.maximum_circuit_trials = 4;
    assert!(matches!(
        circuit.step_candidate(0.1, sources.clone(), &[1.0], bounded, || false),
        Err(HystereticCircuitError::TrialBudget)
    ));
    assert_eq!(circuit, saved);
    let a = circuit
        .step_candidate(0.1, sources.clone(), &[1.0], controls(), || false)
        .unwrap();
    let b = circuit
        .step_candidate(0.1, sources, &[1.0], controls(), || false)
        .unwrap();
    assert_eq!(
        a.circuit(),
        b.circuit(),
        "refused trials do not alter a retry"
    );

    let limited = linear(0.2, 2.0);
    let start = initial(&limited, 0.0, vec![], vec![], vec![0.0]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 0.2 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let limited = DielectricHystereticCircuit::new(graph, 0, limited, vec![], start).unwrap();
    let saved = limited.clone();
    assert!(matches!(
        limited.step_candidate(1.0, vec![], &[1.0], controls(), || false),
        Err(HystereticCircuitError::Material(_))
    ));
    assert_eq!(limited, saved);
}

#[test]
fn overlong_lc_interval_detects_current_reversal_before_publishing_history() {
    let capacitor = linear(1.0, 100.0);
    let start = initial(&capacitor, 0.0, vec![1.0], vec![], vec![]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    };
    let circuit = DielectricHystereticCircuit::new(graph, 1, capacitor, vec![], start).unwrap();
    let saved = circuit.clone();
    // Midpoint gives q=-.8 and phi=-.6 at dt=4: capacitor current starts
    // negative and ends positive, so this cannot be one monotone segment.
    match circuit.step_candidate(4.0, vec![], &[], controls(), || false) {
        Err(HystereticCircuitError::ElectricReversal {
            charge_change_c,
            current_start_mid_end_a,
        }) => {
            near(charge_change_c, -0.8, 1e-10);
            near(current_start_mid_end_a[0], -1.0, 1e-10);
            near(current_start_mid_end_a[2], 0.6, 1e-10);
        }
        Err(other) => panic!("expected a detected LC reversal, got {other}"),
        Ok(_) => panic!("an interval crossing the charge turn must refuse"),
    }
    assert_eq!(circuit, saved);
    let accepted = circuit
        .step_candidate(0.1, vec![], &[], controls(), || false)
        .unwrap();
    near(accepted.step().stored_after_j, 0.5, 1e-10);
    near(accepted.step().heat_j, 0.0, 0.0);
}
