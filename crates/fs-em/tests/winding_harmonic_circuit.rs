//! Original source-driven networks coupled to actual spatial magnetic losses.
//! These authored fixtures test circuit/field identities, not measured devices.

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use fs_em::{
    planar::{
        eddy::{
            EddyConductor,
            harmonic::{
                HarmonicWindingControls, PlanarHarmonicWindings,
                circuit::{HarmonicCircuitControls, HarmonicCircuitError, HarmonicSourceKind},
            },
        },
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        WindingConfig,
        incremental::{IncrementalWindingConfig, PlanarIncrementalWindingBasis},
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_math::c64::C64;
use fs_phs::{
    circuit::{Branch, CircuitGraph, MutualInductance, TransformerLink},
    harmonic::{HarmonicBudget, HarmonicError},
};

fn basis(problem: &PlanarNonlinearMagnetostatics) -> PlanarIncrementalWindingBasis<'_> {
    PlanarIncrementalWindingBasis::new(
        problem,
        vec![
            vec![24.0, 0.0, 0.0, 0.0, -24.0, 0.0, 0.0, 0.0],
            vec![0.0, -12.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0],
        ],
        &common::boundary(),
        vec![2.0, 0.0],
        vec![0.02, 0.03],
        &[0.0; 8],
        IncrementalWindingConfig {
            winding: WindingConfig {
                max_windings: 4,
                max_basis_values: 10_000,
                pairing_tolerance_h: 1e-10,
                current_balance_tolerance_turns: 1e-10,
                inverse_residual_tolerance: 1e-10,
            },
            minimum_smooth_field_margin_t: 1e-8,
        },
        common::material::controls().nonlinear,
        &mut || false,
    )
    .unwrap()
}
fn prepared<'b, 'm>(
    basis: &'b PlanarIncrementalWindingBasis<'m>,
    source: &MaterialCard,
) -> PlanarHarmonicWindings<'b, 'm> {
    let conductor =
        common::material::conductor(source, 300.0, ElectricalPropertyBasis::Resistivity);
    let eddies = [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 9 },
            conductor: conductor.clone(),
        })
        .collect();
    let copper = common::conductors(source, 300.0)
        .into_iter()
        .filter(|source| [0, 1, 4, 5].contains(&source.triangle))
        .collect();
    PlanarHarmonicWindings::new(
        basis,
        eddies,
        &[0, 1, 4, 5],
        copper,
        HarmonicWindingControls {
            relative_tolerance: 1e-11,
            absolute_residual_a_m: 1e-12,
            restart: 8,
            max_cycles: 50,
            max_krylov_scalars: 10_000,
            max_field_scalars: 10_000,
            impedance_tolerance_ohm: 1e-8,
            power_tolerance_w: 1e-8,
            current_balance_tolerance_a: 1e-10,
            max_cell_skin_ratio: 1000.0,
        },
        &mut || false,
    )
    .unwrap()
}
fn controls() -> HarmonicCircuitControls {
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
        voltage_tolerance_v: 1e-8,
        current_tolerance_a: 1e-8,
        power_tolerance_w: 1e-8,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn complex_near(actual: C64, expected: C64) {
    near(actual.re, expected.re);
    near(actual.im, expected.im);
}
fn power(v: C64, i: C64) -> C64 {
    (i.conj() * v).scale(0.5)
}
fn loaded() -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: f64::NAN }),
            (2, 0, Branch::Resistor { ohms: 1.7 }),
            (3, 0, Branch::Inductor { henries: f64::NAN }),
            (0, 3, Branch::Resistor { ohms: 4.2 }),
        ],
        transformers: vec![],
    }
}
fn two_port(z: &[C64], drive: C64, first_load: C64, second_load: C64) -> [C64; 2] {
    let a = z[0] + first_load;
    let d = z[3] + second_load;
    let determinant = a * d - z[1] * z[2];
    [drive * d / determinant, -drive * z[2] / determinant]
}

#[test]
fn voltage_source_and_secondary_load_follow_full_lossy_mutual_impedance() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.7);
    let basis = basis(&problem);
    let model = prepared(&basis, &source);
    let drive = C64::new(2.0, -0.7);
    for omega in [0.7, 13.0, 60.0] {
        let frequency = model.at_frequency(omega, &mut || false).unwrap();
        let expected = two_port(
            frequency.impedance_ohm(),
            drive,
            C64::from_re(1.7),
            C64::from_re(4.2),
        );
        let mut previous = None;
        for scale in [0.25, 1.0, 4.0] {
            let mut limits = controls();
            limits.reference_resistance_scale = scale;
            let response = frequency
                .circuit_response(&loaded(), &[1, 3], &[], &[drive], limits, &mut || false)
                .unwrap();
            for (actual, want) in response.field.currents_a.iter().zip(expected) {
                complex_near(*actual, want);
            }
            complex_near(response.branches[1].current_a, expected[0]);
            complex_near(response.branches[3].current_a, expected[1]);
            complex_near(response.branches[4].current_a, expected[1]);
            complex_near(response.source_outputs[0], expected[0]);
            complex_near(response.supplied_power_va, power(drive, expected[0]));
            near(
                response.mean_resistor_power_w,
                0.5 * (1.7 * expected[0].norm_sq() + 4.2 * expected[1].norm_sq()),
            );
            near(
                response.supplied_power_va.re,
                response.mean_resistor_power_w + response.field.mean_dissipated_power_w(),
            );
            assert!(response.field.mean_eddy_power_w > 0.0);
            assert!(response.field.mean_winding_power_w > 0.0);
            assert!(response.supplied_power_va.im < 0.0);
            near(response.mean_network_stored_energy_j, 0.0);
            complex_near(response.power_residual_va, C64::ZERO);
            assert!(response.maximum_kirchhoff_current_defect_a < 1e-8);
            for cell in &response.field.eddy_cells {
                assert_eq!(
                    cell.source.conductor.material().card_identity(),
                    source.content_hash()
                );
            }
            if let Some(old) = previous {
                complex_near(response.field.currents_a[0], old);
            }
            previous = Some(response.field.currents_a[0]);
        }
    }
}

#[test]
fn ideal_current_and_voltage_clamps_preserve_original_source_power_conventions() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.7);
    let basis = basis(&problem);
    let model = prepared(&basis, &source);
    let frequency = model.at_frequency(17.0, &mut || false).unwrap();
    let current = C64::new(0.2, -0.1);
    let voltage = C64::new(-0.3, 0.8);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 0.0 }),
            (2, 0, Branch::Inductor { henries: 0.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
            (2, 0, Branch::VoltageSource { port: 1 }),
        ],
        transformers: vec![],
    };
    let z = frequency.impedance_ohm();
    let secondary = (voltage - z[2] * current) / z[3];
    let primary_voltage = z[0] * current + z[1] * secondary;
    let response = frequency
        .circuit_response(
            &graph,
            &[0, 1],
            &[],
            &[current, voltage],
            controls(),
            &mut || false,
        )
        .unwrap();
    complex_near(response.field.currents_a[0], current);
    complex_near(response.field.currents_a[1], secondary);
    complex_near(response.source_outputs[0], primary_voltage);
    complex_near(response.source_outputs[1], secondary);
    assert_eq!(response.sources[0].kind, HarmonicSourceKind::Current);
    assert_eq!(response.sources[1].kind, HarmonicSourceKind::Voltage);
    complex_near(
        response.sources[0].supplied_power_va,
        power(primary_voltage, current),
    );
    complex_near(
        response.sources[1].supplied_power_va,
        power(voltage, secondary),
    );
    complex_near(response.branches[2].current_a, -current);
    complex_near(response.branches[3].current_a, -secondary);
    complex_near(response.power_residual_va, C64::ZERO);
}

#[test]
fn shared_voltage_control_retains_individual_branch_currents_and_total_output() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.7);
    let basis = basis(&problem);
    let model = prepared(&basis, &source);
    let frequency = model.at_frequency(12.0, &mut || false).unwrap();
    let drive = C64::new(1.0, 0.4);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 0.0 }),
            (2, 0, Branch::Inductor { henries: 0.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::VoltageSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let z = frequency.impedance_ohm();
    let determinant = z[0] * z[3] - z[1] * z[2];
    let currents = [
        drive * (z[3] - z[1]) / determinant,
        drive * (z[0] - z[2]) / determinant,
    ];
    let response = frequency
        .circuit_response(&graph, &[0, 1], &[], &[drive], controls(), &mut || false)
        .unwrap();
    for (index, &current) in currents.iter().enumerate() {
        complex_near(response.sources[index].output, current);
        complex_near(
            response.sources[index].supplied_power_va,
            power(drive, current),
        );
        assert_eq!(response.sources[index].port, 0);
    }
    complex_near(response.source_outputs[0], currents[0] + currents[1]);
    complex_near(response.power_residual_va, C64::ZERO);
}

#[test]
fn transformer_and_mutually_coupled_external_rlc_survive_exact_series_resonance() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.7);
    let basis = basis(&problem);
    let model = prepared(&basis, &source);
    let omega = 2.0;
    let frequency = model.at_frequency(omega, &mut || false).unwrap();
    let drive = C64::new(3.0, -1.0);
    let ratio = 2.5;
    let total_l = 0.7;
    let graph = CircuitGraph {
        node_count: 8,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 3, Branch::Inductor { henries: 0.0 }),
            (3, 4, Branch::Resistor { ohms: 1.7 }),
            (4, 5, Branch::Inductor { henries: 0.25 }),
            (5, 6, Branch::Inductor { henries: 0.25 }),
            (
                6,
                0,
                Branch::Capacitor {
                    farads: 1.0 / (omega * omega * total_l),
                },
            ),
            (7, 0, Branch::Inductor { henries: 0.0 }),
            (0, 7, Branch::Resistor { ohms: 4.2 }),
        ],
        transformers: vec![TransformerLink {
            primary: (1, 0),
            secondary: (2, 0),
            ratio,
        }],
    };
    let expected = two_port(
        frequency.impedance_ohm(),
        drive.scale(1.0 / ratio),
        C64::from_re(1.7),
        C64::from_re(4.2),
    );
    let response = frequency
        .circuit_response(
            &graph,
            &[1, 6],
            &[MutualInductance {
                branches: (3, 4),
                henries: 0.1,
            }],
            &[drive],
            controls(),
            &mut || false,
        )
        .unwrap();
    complex_near(response.field.currents_a[0], expected[0]);
    complex_near(response.field.currents_a[1], expected[1]);
    for branch in [2, 3, 4, 5] {
        complex_near(response.branches[branch].current_a, expected[0]);
    }
    complex_near(response.node_voltages_v[2], drive.scale(1.0 / ratio));
    complex_near(
        response.transformers[0].primary_current_a,
        expected[0].scale(1.0 / ratio),
    );
    complex_near(response.transformers[0].secondary_current_a, -expected[0]);
    complex_near(response.transformers[0].power_residual_va, C64::ZERO);
    complex_near(response.source_outputs[0], expected[0].scale(1.0 / ratio));
    near(
        response.mean_network_stored_energy_j,
        0.5 * total_l * expected[0].norm_sq(),
    );
    near(
        response.supplied_power_va.im,
        response.field.reactive_power_var,
    );
    complex_near(response.power_residual_va, C64::ZERO);
}

#[test]
fn unsupported_bindings_budgets_and_cancelled_work_never_publish_a_response() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.7);
    let basis = basis(&problem);
    let model = prepared(&basis, &source);
    let frequency = model.at_frequency(1.0, &mut || false).unwrap();
    let graph = loaded();
    let original = frequency.impedance_ohm().to_vec();
    for mapping in [vec![], vec![1], vec![1, 1], vec![1, 2], vec![1, 99]] {
        assert!(
            frequency
                .circuit_response(&graph, &mapping, &[], &[C64::ONE], controls(), &mut || {
                    false
                })
                .is_err()
        );
    }
    assert!(
        frequency
            .circuit_response(
                &graph,
                &[1, 3],
                &[MutualInductance {
                    branches: (1, 3),
                    henries: 0.1
                }],
                &[C64::ONE],
                controls(),
                &mut || false
            )
            .is_err()
    );
    for inputs in [
        vec![],
        vec![C64::ONE, C64::ONE],
        vec![C64::new(f64::NAN, 0.0)],
    ] {
        assert!(
            frequency
                .circuit_response(&graph, &[1, 3], &[], &inputs, controls(), &mut || false)
                .is_err()
        );
    }
    for limit in [0, 1, 2] {
        let mut limits = controls();
        match limit {
            0 => limits.max_interface_scalars = 1,
            1 => limits.native.max_dense_bytes = 1,
            _ => limits.native.max_dense_work = 1,
        }
        assert!(
            frequency
                .circuit_response(&graph, &[1, 3], &[], &[C64::ONE], limits, &mut || false)
                .is_err()
        );
    }
    let mut mixed = graph.clone();
    mixed
        .branches
        .push((3, 0, Branch::CurrentSource { port: 0 }));
    assert!(
        frequency
            .circuit_response(&mixed, &[1, 3], &[], &[C64::ONE], controls(), &mut || false)
            .is_err()
    );
    for stop in [1, 10, 50] {
        let mut calls = 0;
        assert!(matches!(
            frequency.circuit_response(&graph, &[1, 3], &[], &[C64::ONE], controls(), &mut || {
                calls += 1;
                calls == stop
            }),
            Err(HarmonicCircuitError::Cancelled)
        ));
    }
    // A disconnected lossless mode at its true pole is not damped by the
    // reference termination on an entirely different magnetic terminal.
    let mut pole = graph.clone();
    pole.node_count += 1;
    pole.branches.extend([
        (4, 0, Branch::Inductor { henries: 1.0 }),
        (4, 0, Branch::Capacitor { farads: 1.0 }),
    ]);
    assert!(matches!(
        frequency.circuit_response(&pole, &[1, 3], &[], &[C64::ONE], controls(), &mut || false),
        Err(HarmonicCircuitError::Native(HarmonicError::Singular { .. }))
    ));
    assert_eq!(frequency.impedance_ohm(), original);
    assert!(
        frequency
            .circuit_response(&graph, &[1, 3], &[], &[C64::ONE], controls(), &mut || false)
            .is_ok()
    );
}
