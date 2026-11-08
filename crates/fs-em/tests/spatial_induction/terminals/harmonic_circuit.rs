//! Original AC networks drive physical electrodes and the resolved 3D current.
use super::*;
use crate::harmonic::{complex_near, harmonic_controls};
use fs_em::spatial::vector_potential::induction::{
    harmonic::{
        circuit::HarmonicCircuitControls,
        terminals::{SpatialHarmonicTerminalControls, SpatialHarmonicTerminalFrequency},
    },
    terminals::current::SpatialCurrentInduction,
};
use fs_math::c64::C64;
use fs_phs::{
    circuit::{Branch, CircuitGraph},
    harmonic::HarmonicBudget,
};

fn terminal_controls() -> SpatialHarmonicTerminalControls {
    SpatialHarmonicTerminalControls {
        field: harmonic_controls(),
        maximum_ports: 4,
        maximum_response_scalars: 100_000,
        current_probe_a: 0.25,
    }
}

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

fn current_source<'m>(
    problem: &'m SpatialVectorMagnetostatics,
    sources: Vec<Option<ResolvedConductorState>>,
    faces: &[Vec<usize>],
) -> SpatialCurrentInduction<'m> {
    let initial = problem
        .solve(
            &vec![[0.0; 3]; problem.mesh().tets.len()],
            &trace(problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    SpatialTerminalInduction::new(native, faces, preparation(), &mut || false)
        .unwrap()
        .into_current_driven(preparation(), &mut || false)
        .unwrap()
}

fn series_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 5,
        branches: vec![
            (4, 0, Branch::VoltageSource { port: 0 }),
            (4, 3, Branch::Resistor { ohms: 0.4 }),
            (3, 2, Branch::Inductor { henries: 0.03 }),
            (2, 1, Branch::Capacitor { farads: 0.7 }),
            (1, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    }
}

#[test]
fn series_rlc_matches_independent_diffusion_impedance_and_physical_orientation() {
    // Independently integrated six-tet cube: int W_x=1/6, int |W|^2=1/5,
    // K=2 and sigma=5. Eliminating electrode voltage under imposed current
    // gives D=K-i*w*sigma*(1/5-1/36), Z=1/sigma-i*w/(36*D).
    // The native field is inductive and dispersive, not its DC resistance.
    let (problem, sources) = cube(1, None);
    let graph = series_graph();
    let drive = C64::new(1.2, -0.7);
    let velocity = [[0.0; 3]; 8];
    for reverse in [false, true] {
        let mut faces = electrodes(&problem);
        if reverse {
            faces.reverse();
        }
        let source = current_source(&problem, sources.clone(), &faces);
        for omega in [0.3, 3.0, 40.0] {
            let frequency = SpatialHarmonicTerminalFrequency::new(
                &source,
                omega,
                &velocity,
                terminal_controls(),
                &mut || false,
            )
            .unwrap();
            assert_eq!(frequency.ports().len(), 1);
            assert_eq!(frequency.ports()[0].terminal, 1);
            assert_eq!(frequency.ports()[0].reference_terminal, 0);
            assert_eq!(frequency.unit_reports().len(), 1);
            let d = C64::new(2.0, -omega * 5.0 * (1.0 / 5.0 - 1.0 / 36.0));
            let z = C64::from_re(0.2) + C64::new(0.0, -omega / 36.0) / d;
            complex_near(frequency.impedance_ohm()[0], z, 3e-9);
            let expected = drive / (z + C64::new(0.4, 1.0 / (omega * 0.7) - omega * 0.03));
            let a = expected.scale(if reverse { 1.0 / 6.0 } else { -1.0 / 6.0 }) / d;
            for scale in [0.25, 4.0] {
                let response = frequency
                    .circuit_response(
                        &graph,
                        &[4],
                        &[],
                        &[drive],
                        HarmonicCircuitControls {
                            reference_resistance_scale: scale,
                            ..circuit_controls()
                        },
                        &mut || false,
                    )
                    .unwrap();
                complex_near(response.field.currents_a[0], expected, 3e-8);
                complex_near(response.field.terminal_voltages_v[0], z * expected, 3e-8);
                complex_near(response.branches[4].current_a, expected, 3e-8);
                complex_near(response.source_outputs[0], expected, 3e-8);
                complex_near(
                    response.field.field.potential_wb()[diagonal(&problem)],
                    a,
                    3e-8,
                );
                near(
                    response.mean_field_power_w,
                    0.5 * z.re * expected.norm_sq(),
                    3e-8,
                );
                near(
                    response.mean_resistor_power_w,
                    0.2 * expected.norm_sq(),
                    3e-8,
                );
                near(response.mean_field_mechanical_power_w, 0.0, 1e-15);
                near(
                    response.mean_network_stored_energy_j,
                    0.25 * expected.norm_sq() * (0.03 + 1.0 / (omega * omega * 0.7)),
                    3e-8,
                );
                complex_near(
                    response.supplied_power_va,
                    (drive * expected.conj()).scale(0.5),
                    3e-8,
                );
                complex_near(response.power_residual_va, C64::ZERO, 8e-8);
                assert!(response.maximum_kirchhoff_current_defect_a < 3e-8);
            }
        }
        assert!(source.potential_wb().iter().all(|&a| a == 0.0));
        assert_eq!(source.time_s(), 0.0);
    }
}

#[test]
fn moving_electrode_circuit_preserves_independently_integrated_signed_mechanical_work() {
    let (problem, sources) = cube(2, None);
    let source = current_source(&problem, sources, &electrodes(&problem));
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 0.4 }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    };
    let drive = C64::new(0.7, -0.2);
    let mut signed_work = Vec::new();
    for speed in [-0.03, 0.03] {
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| p.map(|x| speed * x))
            .collect();
        let frequency = SpatialHarmonicTerminalFrequency::new(
            &source,
            3.0,
            &velocity,
            terminal_controls(),
            &mut || false,
        )
        .unwrap();
        let response = frequency
            .circuit_response(&graph, &[2], &[], &[drive], circuit_controls(), &mut || {
                false
            })
            .unwrap();
        complex_near(response.field.currents_a[0], drive, 3e-8);
        let field = &response.field.field;
        let (mut heat, mut motion) = (0.0, 0.0);
        for cell in field.cells() {
            let tet = problem.mesh().tets[cell.cell];
            let b = field.flux_density_t()[cell.cell];
            let rho = source.conductors()[cell.cell]
                .as_ref()
                .unwrap()
                .resistivity_ohm_m();
            for a in 0..4 {
                for j in 0..4 {
                    // Exact tetrahedral P1 products, independently evaluated
                    // from the retained current, magnetic field and motion.
                    let weight = 0.5 * cell.volume_m3 * if a == j { 0.1 } else { 0.05 };
                    let ja = cell.current_density_a_m2[a];
                    let jj = cell.current_density_a_m2[j];
                    heat += weight * rho * (0..3).map(|k| (ja[k] * jj[k].conj()).re).sum::<f64>();
                    let rr = cross(jj.map(|v| v.re), b.map(|v| v.re));
                    let ii = cross(jj.map(|v| v.im), b.map(|v| v.im));
                    motion += weight
                        * (0..3)
                            .map(|k| velocity[tet[a] as usize][k] * (rr[k] + ii[k]))
                            .sum::<f64>();
                }
            }
        }
        assert!(heat > 1e-6 && motion.abs() > 1e-9);
        near(response.mean_field_power_w, heat, 3e-8);
        near(response.mean_field_mechanical_power_w, motion, 3e-8);
        near(response.mean_resistor_power_w, 0.2 * drive.norm_sq(), 3e-8);
        near(
            response.supplied_power_va.re,
            heat + motion + 0.2 * drive.norm_sq(),
            8e-8,
        );
        complex_near(response.power_residual_va, C64::ZERO, 8e-8);
        near(field.balance().mean_impressed_power_w, 0.0, 1e-15);
        signed_work.push(motion);
    }
    assert!(
        signed_work[0] * signed_work[1] < 0.0,
        "signed motion: {signed_work:?}"
    );
    assert_eq!(source.time_s(), 0.0);
}

#[test]
fn three_electrode_ports_share_one_physical_reference_and_match_voltage_control() {
    let (problem, sources) = cube(2, None);
    let planes = electrodes(&problem);
    let first = planes[1][0];
    let second = *planes[1]
        .iter()
        .find(|&&face| {
            problem.mesh().faces[face]
                .iter()
                .all(|v| !problem.mesh().faces[first].contains(v))
        })
        .unwrap();
    let source = current_source(
        &problem,
        sources,
        &[planes[0].clone(), vec![first], vec![second]],
    );
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let frequency = SpatialHarmonicTerminalFrequency::new(
        &source,
        3.0,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        frequency
            .ports()
            .iter()
            .map(|p| (p.terminal, p.reference_terminal))
            .collect::<Vec<_>>(),
        vec![(1, 0), (2, 0)]
    );
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::VoltageSource { port: 1 }),
            (1, 0, Branch::Inductor { henries: f64::NAN }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    };
    let drive = [C64::new(0.9, -0.4), C64::new(-0.2, 0.3)];
    let voltage_source = source
        .clone()
        .into_voltage_driven(preparation(), &mut || false)
        .unwrap();
    let voltage = voltage_source
        .solve_harmonic_voltage(
            &[C64::ZERO, drive[0], drive[1]],
            &vec![[C64::ZERO; 3]; problem.mesh().tets.len()],
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let response = frequency
        .circuit_response(
            &graph,
            &[2, 3],
            &[],
            &drive,
            circuit_controls(),
            &mut || false,
        )
        .unwrap();
    for port in 0..2 {
        complex_near(
            response.field.currents_a[port],
            voltage.terminals()[port + 1].incoming_current_a,
            3e-8,
        );
        complex_near(response.field.terminal_voltages_v[port], drive[port], 3e-8);
    }
    for (&actual, &expected) in response
        .field
        .field
        .potential_wb()
        .iter()
        .zip(voltage.potential_wb())
    {
        complex_near(actual, expected, 3e-8);
    }
    near(
        response.mean_field_power_w,
        voltage.balance().mean_ohmic_power_w,
        8e-8,
    );
    let mut inconsistent = graph.clone();
    inconsistent.node_count = 4;
    inconsistent.branches[3].1 = 3;
    inconsistent
        .branches
        .push((3, 0, Branch::Resistor { ohms: 1.0 }));
    assert!(
        frequency
            .circuit_response(
                &inconsistent,
                &[2, 3],
                &[],
                &drive,
                circuit_controls(),
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn terminal_preparation_and_cancelled_actual_drive_respect_budgets_and_keep_source() {
    let (problem, sources) = cube(1, None);
    let source = current_source(&problem, sources, &electrodes(&problem));
    let before = source.potential_wb().to_vec();
    let velocity = [[0.0; 3]; 8];
    for limited in [
        SpatialHarmonicTerminalControls {
            maximum_ports: 0,
            ..terminal_controls()
        },
        SpatialHarmonicTerminalControls {
            maximum_response_scalars: 0,
            ..terminal_controls()
        },
        SpatialHarmonicTerminalControls {
            current_probe_a: 0.0,
            ..terminal_controls()
        },
    ] {
        assert!(
            SpatialHarmonicTerminalFrequency::new(&source, 3.0, &velocity, limited, &mut || false)
                .is_err()
        );
    }
    assert!(
        SpatialHarmonicTerminalFrequency::new(
            &source,
            3.0,
            &velocity,
            terminal_controls(),
            &mut || true
        )
        .is_err()
    );
    let frequency = SpatialHarmonicTerminalFrequency::new(
        &source,
        3.0,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    let drive = [C64::new(0.8, -0.3)];
    let mut polls = 0;
    let accepted = frequency
        .current_response(&drive, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 2);
    let mut seen = 0;
    assert!(
        frequency
            .current_response(&drive, &mut || {
                seen += 1;
                seen >= polls / 2
            })
            .is_err()
    );
    let replay = frequency.current_response(&drive, &mut || false).unwrap();
    assert_eq!(replay.field.potential_wb(), accepted.field.potential_wb());
    assert_eq!(replay.terminal_voltages_v, accepted.terminal_voltages_v);
    assert!(frequency.current_response(&[], &mut || false).is_err());
    assert!(
        frequency
            .current_response(&[C64::new(f64::NAN, 0.0)], &mut || false)
            .is_err()
    );
    let limited = HarmonicCircuitControls {
        max_interface_scalars: 0,
        ..circuit_controls()
    };
    assert!(
        frequency
            .circuit_response(&series_graph(), &[4], &[], &drive, limited, &mut || false)
            .is_err()
    );
    assert!(
        frequency
            .circuit_response(
                &series_graph(),
                &[4],
                &[],
                &drive,
                circuit_controls(),
                &mut || true
            )
            .is_err()
    );
    assert_eq!(source.potential_wb(), before);
    assert_eq!(source.time_s(), 0.0);
}
