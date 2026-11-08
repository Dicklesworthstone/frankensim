//! Original native networks closed on source-bound 3D harmonic windings.
use super::winding::{basis, end_model, fixture, impressed, winding_controls};
use super::*;
use fs_em::spatial::vector_potential::induction::harmonic::{
    circuit::{HarmonicCircuitControls, HarmonicSourceKind},
    winding::{SpatialHarmonicWindingError, SpatialHarmonicWindings},
};
use fs_phs::{
    circuit::{Branch, CircuitGraph},
    harmonic::HarmonicBudget,
};
use std::cell::Cell;

fn circuit_controls() -> HarmonicCircuitControls {
    HarmonicCircuitControls {
        native: HarmonicBudget {
            max_dimension: 128,
            max_dense_bytes: 4_000_000,
            max_dense_work: 10_000_000,
            relative_residual_tolerance: 1e-11,
        },
        max_interface_scalars: 100_000,
        reference_resistance_scale: 1.0,
        relative_terminal_residual_tolerance: 1e-11,
        voltage_tolerance_v: 3e-8,
        current_tolerance_a: 3e-8,
        power_tolerance_w: 8e-8,
    }
}

// Deliberately reverse physical winding order relative to original branches,
// and leave a zero source-control gap at index one.
fn voltage_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 0.7 }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 2 }),
            (3, 4, Branch::Resistor { ohms: 1.1 }),
            (4, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    }
}

fn power(voltage: C64, current: C64) -> C64 {
    (current.conj() * voltage).scale(0.5)
}

#[test]
fn loaded_3d_circuit_matches_independent_two_loop_equations_and_source_phase() {
    let f = fixture(2, false);
    let boundary = trace(&f.problem, [0.0; 3], false);
    // Three original static energy solves determine L by polarization. This
    // reference never reads the harmonic impedance or winding-basis L matrix.
    let energy = |current: &[C64]| {
        let source: Vec<_> = impressed(&f.profiles, current)
            .into_iter()
            .map(|j| j.map(|v| v.re))
            .collect();
        f.problem
            .solve(&source, &boundary, static_controls(), &mut || false)
            .unwrap()
            .stored_energy_j()
    };
    let e0 = energy(&[C64::ONE, C64::ZERO]);
    let e1 = energy(&[C64::ZERO, C64::ONE]);
    let es = energy(&[C64::ONE, C64::ONE]);
    // Authored diagonal ends plus 0.07 H * [1,-0.4][1,-0.4]^T.
    let l = [
        2.0 * e0 + 0.02 + 0.07,
        es - e0 - e1 - 0.028,
        es - e0 - e1 - 0.028,
        2.0 * e1 + 0.03 + 0.0112,
    ];
    assert!(l[1].abs() > 1e-5);
    let mut resistance = [0.0; 2];
    for source in &f.sources {
        let squared: f64 = f.profiles[source.winding][source.cell]
            .iter()
            .map(|v| v * v)
            .sum();
        resistance[source.winding] +=
            source.conductor.resistivity_ohm_m() * squared / (48.0 * source.fill_fraction);
    }
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &boundary,
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let graph = voltage_graph();
    let inputs = [C64::new(1.3, -0.2), C64::ZERO, C64::new(-0.6, 0.4)];
    for omega in [0.7, 9.0] {
        let frequency = prepared
            .at_frequency(omega, &f.velocity, &mut || false)
            .unwrap();
        let a = C64::new(resistance[0] + 1.1, -omega * l[0]);
        let b = C64::new(0.0, -omega * l[1]);
        let c = C64::new(0.0, -omega * l[2]);
        let d = C64::new(resistance[1] + 0.7, -omega * l[3]);
        let determinant = a * d - b * c;
        let expected = [
            (d * inputs[2] - b * inputs[0]) / determinant,
            (a * inputs[0] - c * inputs[2]) / determinant,
        ];
        let expected_supply = power(inputs[2], expected[0]) + power(inputs[0], expected[1]);
        let expected_resistor_heat =
            0.5 * (1.1 * expected[0].norm_sq() + 0.7 * expected[1].norm_sq());
        let expected_copper =
            0.5 * (resistance[0] * expected[0].norm_sq() + resistance[1] * expected[1].norm_sq());
        for scale in [0.25, 1.0, 4.0] {
            for phase in [C64::ONE, C64::new(0.6, 0.8)] {
                let drive = inputs.map(|v| v * phase);
                let response = frequency
                    .circuit_response(
                        &graph,
                        &[5, 2],
                        &[],
                        &drive,
                        HarmonicCircuitControls {
                            reference_resistance_scale: scale,
                            ..circuit_controls()
                        },
                        &mut || false,
                    )
                    .unwrap();
                for w in 0..2 {
                    complex_near(response.field.currents_a[w], expected[w] * phase, 3e-8);
                }
                complex_near(response.branches[5].current_a, expected[0] * phase, 3e-8);
                complex_near(response.branches[2].current_a, expected[1] * phase, 3e-8);
                complex_near(response.source_outputs[2], expected[0] * phase, 3e-8);
                complex_near(response.source_outputs[0], expected[1] * phase, 3e-8);
                assert_eq!(response.source_outputs[1], C64::ZERO);
                assert!(
                    response
                        .sources
                        .iter()
                        .all(|s| s.kind == HarmonicSourceKind::Voltage)
                );
                complex_near(response.supplied_power_va, expected_supply, 5e-8);
                near(response.mean_resistor_power_w, expected_resistor_heat, 5e-8);
                near(response.mean_field_power_w, expected_copper, 5e-8);
                near(response.mean_field_mechanical_power_w, 0.0, 1e-15);
                near(response.mean_network_stored_energy_j, 0.0, 1e-15);
                complex_near(response.power_residual_va, C64::ZERO, 8e-8);
                near(
                    response.supplied_power_va.im,
                    -2.0 * omega * response.field.mean_magnetic_energy_j,
                    8e-8,
                );
                assert!(
                    response.maximum_kirchhoff_current_defect_a
                        < circuit_controls().current_tolerance_a
                );
                assert!(
                    response.maximum_terminal_voltage_defect_v
                        < circuit_controls().voltage_tolerance_v
                );
            }
        }
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn moving_eddy_circuit_keeps_signed_mechanical_work_separate_from_physical_heat() {
    let f = fixture(3, true);
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let omega = 3.0;
    let frequency = prepared
        .at_frequency(omega, &f.velocity, &mut || false)
        .unwrap();
    let currents = [C64::new(0.7, -0.2), C64::new(-0.3, 0.4)];
    // Ideal current sources establish both currents independently of the field
    // impedance; their reported voltage and signed work remain unknown outputs.
    let graph = CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 2 }),
            (1, 2, Branch::Resistor { ohms: 0.7 }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::CurrentSource { port: 0 }),
            (3, 4, Branch::Resistor { ohms: 1.1 }),
            (4, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    };
    let direct = state
        .solve_harmonic(
            &impressed(&f.profiles, &currents),
            &f.velocity,
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let (mut heat, mut motion, mut reactive_motion) = (0.0, 0.0, 0.0);
    for record in direct.cells() {
        let tet = f.problem.mesh().tets[record.cell];
        let rho = f.eddy[record.cell].as_ref().unwrap().resistivity_ohm_m();
        let b = direct.flux_density_t()[record.cell];
        let br = b.map(|v| v.re);
        let bi = b.map(|v| v.im);
        for a in 0..4 {
            for b in 0..4 {
                // Independent exact tetrahedral P1 product quadrature.
                let mass = record.volume_m3 * if a == b { 0.1 } else { 0.05 };
                let ja = record.current_density_a_m2[a];
                let jb = record.current_density_a_m2[b];
                heat += 0.5
                    * mass
                    * rho
                    * (0..3)
                        .map(|k| ja[k].re * jb[k].re + ja[k].im * jb[k].im)
                        .sum::<f64>();
                let rr = cross(jb.map(|v| v.re), br);
                let ii = cross(jb.map(|v| v.im), bi);
                let ri = cross(jb.map(|v| v.re), bi);
                let ir = cross(jb.map(|v| v.im), br);
                for k in 0..3 {
                    motion += 0.5 * mass * f.velocity[tet[a] as usize][k] * (rr[k] + ii[k]);
                    reactive_motion +=
                        0.5 * mass * f.velocity[tet[a] as usize][k] * (ri[k] - ir[k]);
                }
            }
        }
    }
    assert!(heat > 1e-7 && motion.abs() > 1e-12);
    near(heat, direct.balance().mean_ohmic_power_w, 2e-8);
    near(motion, direct.balance().mean_mechanical_power_w, 2e-8);
    near(
        reactive_motion,
        direct.balance().reactive_mechanical_power_var,
        2e-8,
    );
    let copper: f64 = f
        .sources
        .iter()
        .map(|source| {
            let squared: f64 = f.profiles[source.winding][source.cell]
                .iter()
                .map(|v| v * v)
                .sum();
            0.5 * source.conductor.resistivity_ohm_m()
                * squared
                * currents[source.winding].norm_sq()
                / (162.0 * source.fill_fraction)
        })
        .sum();
    let resistor_heat = 0.5 * (0.7 * currents[0].norm_sq() + 1.1 * currents[1].norm_sq());
    let end_energy = 0.25
        * (0.02 * currents[0].norm_sq()
            + 0.03 * currents[1].norm_sq()
            + 0.07 * (currents[0] - currents[1].scale(0.4)).norm_sq());
    for phase in [C64::ONE, C64::new(-0.8, 0.6)] {
        let inputs = [currents[1] * phase, C64::ZERO, currents[0] * phase];
        let response = frequency
            .circuit_response(
                &graph,
                &[2, 5],
                &[],
                &inputs,
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        for w in 0..2 {
            complex_near(response.field.currents_a[w], currents[w] * phase, 2e-9);
        }
        for (actual, original) in response
            .field
            .field
            .flux_density_t()
            .iter()
            .flatten()
            .zip(direct.flux_density_t().iter().flatten())
        {
            complex_near(*actual, *original * phase, 2e-8);
        }
        assert!(
            response
                .sources
                .iter()
                .all(|s| s.kind == HarmonicSourceKind::Current)
        );
        near(response.mean_resistor_power_w, resistor_heat, 3e-8);
        near(response.mean_field_power_w, copper + heat, 5e-8);
        near(response.mean_field_mechanical_power_w, motion, 3e-8);
        near(
            response.supplied_power_va.re,
            resistor_heat + copper + heat + motion,
            8e-8,
        );
        near(
            response.supplied_power_va.im,
            reactive_motion - 2.0 * omega * (direct.balance().mean_magnetic_energy_j + end_energy),
            8e-8,
        );
        complex_near(response.power_residual_va, C64::ZERO, 8e-8);
        assert_eq!(response.source_outputs[1], C64::ZERO);
        assert!(
            response.maximum_kirchhoff_current_defect_a < circuit_controls().current_tolerance_a
        );
    }
    // This authored moving shield tests signed work, not a calibrated motor or
    // a claim that its TOTAL electrical supply is negative (generation).
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn circuit_mapping_work_budget_and_late_cancel_preserve_field_and_allow_identical_retry() {
    let f = fixture(2, false);
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let frequency = prepared
        .at_frequency(3.0, &f.velocity, &mut || false)
        .unwrap();
    let graph = voltage_graph();
    let inputs = [C64::new(1.3, -0.2), C64::ZERO, C64::new(-0.6, 0.4)];
    for mapping in [&[2usize, 2][..], &[2][..]] {
        assert!(
            frequency
                .circuit_response(
                    &graph,
                    mapping,
                    &[],
                    &inputs,
                    circuit_controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let mut limited = circuit_controls();
    limited.native.max_dense_work = 0;
    assert!(
        frequency
            .circuit_response(&graph, &[5, 2], &[], &inputs, limited, &mut || false)
            .is_err()
    );
    let mut wrong = inputs;
    wrong[1] = C64::ONE;
    assert!(
        frequency
            .circuit_response(
                &graph,
                &[5, 2],
                &[],
                &wrong,
                circuit_controls(),
                &mut || false
            )
            .is_err()
    );
    let calls = Cell::new(0);
    let expected = frequency
        .circuit_response(
            &graph,
            &[5, 2],
            &[],
            &inputs,
            circuit_controls(),
            &mut || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    for stop in [1, calls.get()] {
        let count = Cell::new(0);
        assert!(matches!(
            frequency.circuit_response(
                &graph,
                &[5, 2],
                &[],
                &inputs,
                circuit_controls(),
                &mut || {
                    count.set(count.get() + 1);
                    count.get() == stop
                }
            ),
            Err(SpatialHarmonicWindingError::Cancelled)
        ));
    }
    let replay = frequency
        .circuit_response(
            &graph,
            &[5, 2],
            &[],
            &inputs,
            circuit_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(replay.field.currents_a, expected.field.currents_a);
    assert_eq!(replay.node_voltages_v, expected.node_voltages_v);
    assert_eq!(
        replay.field.field.potential_wb(),
        expected.field.field.potential_wb()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}
