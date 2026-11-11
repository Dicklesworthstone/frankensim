//! G1/G3/G4: actual grouped-electrode solves, without terminal probe fields.
//! Manufactured layers test current transfer and local work, not measured data.
#[path = "terminals/motion.rs"]
mod motion;

use super::*;
use fs_em::conduction::{
    axisymmetric::AxisymmetricConduction,
    terminals::{ConductionTerminalControl, ConductionTerminalLimits},
};

pub(super) const TERMINALS: [&[usize]; 3] = [&[0, 1], &[2, 3, 4, 5], &[6, 7]];

pub(super) fn terminal_limits() -> ConductionTerminalLimits {
    ConductionTerminalLimits {
        maximum_terminals: 4,
        maximum_workspace_values: 4096,
    }
}

pub(super) fn disconnected_layers(metric: Metric2, layers: usize) -> TriComplex2 {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    for layer in 0..layers {
        let z = layer as f64;
        let first = vertices.len() as u32;
        vertices.extend([[0.0, z], [1.0, z], [1.0, z + 1.0], [0.0, z + 1.0]]);
        faces.extend([[first, first + 1, first + 2], [first, first + 2, first + 3]]);
    }
    TriComplex2::from_indexed_triangles("separate conducting layers", vertices, faces, metric)
        .unwrap()
}

fn conductors() -> Vec<ResolvedConductorState> {
    [1.0, 1.0, 1.0 / 3.0, 1.0 / 3.0]
        .map(|sigma| conductor(sigma, ElectricalPropertyBasis::Conductivity))
        .to_vec()
}

#[test]
fn actual_planar_terminal_constraints_recover_floating_transfer_and_source_work() {
    let area = 0.4;
    let original = PlanarConduction::new(
        disconnected_layers(Metric2::planar(area).unwrap(), 2),
        conductors(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let source_snapshot = original.conductors().to_vec();
    for sign in [-1.0, 1.0] {
        // Two unit lengths: rho=[1,3], |J|=2, V=[0,2,8], P/A=16.
        let currents = [-sign * 2.0 * area, 0.0, sign * 2.0 * area];
        let voltages = [0.0, sign * 2.0, sign * 8.0];
        let ends = [0.0, sign * 8.0];
        let reversed_currents = [sign * 2.0 * area, 0.0];
        for control in [
            ConductionTerminalControl::Currents(&currents),
            ConductionTerminalControl::Voltages(&voltages),
            ConductionTerminalControl::Mixed {
                current_terminals: &[2, 1],
                currents_a: &reversed_currents,
                voltages_v: &[1e250],
            },
            ConductionTerminalControl::Mixed {
                current_terminals: &[1],
                currents_a: &[0.0],
                voltages_v: &ends,
            },
        ] {
            let (field, terminals) = original
                .solve_terminals(
                    &TERMINALS,
                    control,
                    terminal_limits(),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for i in 0..3 {
                near(terminals.voltage_v[i], voltages[i]);
                near(terminals.current_a[i], currents[i]);
            }
            near(terminals.supplied_power_w, 16.0 * area);
            near(terminals.power_residual_w, 0.0);
            assert!(terminals.maximum_current_residual_a <= controls().current_tolerance_a);
            near(field.power_w(), 16.0 * area);
            for triangle in 0..4 {
                let e = if triangle < 2 { 2.0 } else { 6.0 };
                near(field.electric_field_v_m()[triangle][0], 0.0);
                near(field.electric_field_v_m()[triangle][1], -sign * e);
                near(field.current_density_a_m2()[triangle][1], -sign * 2.0);
                near(field.cell_power_w()[triangle], area * e);
            }
            let reaction = |vertices: &[usize]| {
                vertices
                    .iter()
                    .map(|&v| field.nodal_boundary_current_a()[v].unwrap())
                    .sum::<f64>()
            };
            near(reaction(&[2, 3]), sign * 2.0 * area);
            near(reaction(&[4, 5]), -sign * 2.0 * area);
            near(reaction(TERMINALS[1]), 0.0);
        }
        let boundary: Vec<_> = TERMINALS
            .iter()
            .zip(voltages)
            .flat_map(|(vertices, voltage)| vertices.iter().map(move |&v| (v, voltage)))
            .collect();
        let direct = original
            .solve(&boundary, controls(), &mut || false)
            .unwrap();
        let (grouped, _) = original
            .solve_terminals(
                &TERMINALS,
                ConductionTerminalControl::Voltages(&voltages),
                terminal_limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        for (&a, &b) in direct.cell_power_w().iter().zip(grouped.cell_power_w()) {
            near(a, b);
        }
    }
    let reordered = [TERMINALS[1], TERMINALS[2], TERMINALS[0]];
    let (_, report) = original
        .solve_terminals(
            &reordered,
            ConductionTerminalControl::Currents(&[0.0, 2.0 * area, -2.0 * area]),
            terminal_limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&actual, expected) in report.voltage_v.iter().zip([0.0, 6.0, -2.0]) {
        near(actual, expected);
    }
    near(report.supplied_power_w, 16.0 * area);
    assert_eq!(original.conductors(), source_snapshot);
}

#[test]
fn swept_terminal_constraints_keep_exact_disk_reactions_and_local_heat() {
    for angle in [0.7, 2.0, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let original = AxisymmetricConduction::new(
            disconnected_layers(Metric2::axisymmetric(angle).unwrap(), 2),
            conductors(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let (field, report) = original
            .solve_terminals(
                &TERMINALS,
                ConductionTerminalControl::Mixed {
                    current_terminals: &[2, 1],
                    currents_a: &[2.0 * area, 0.0],
                    voltages_v: &[1e250],
                },
                terminal_limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        near(report.voltage_v[1], 2.0);
        near(report.voltage_v[2], 8.0);
        near(report.current_a[0], -2.0 * area);
        near(report.current_a[1], 0.0);
        near(report.current_a[2], 2.0 * area);
        near(report.supplied_power_w, 16.0 * area);
        near(field.power_w(), 16.0 * area);
        // Axis and outer-rim reactions integrate N_i*r, not equal nodal shares.
        near(
            field.nodal_boundary_current_a()[0].unwrap(),
            -2.0 * area / 3.0,
        );
        near(
            field.nodal_boundary_current_a()[1].unwrap(),
            -4.0 * area / 3.0,
        );
        near(
            field.nodal_boundary_current_a()[6].unwrap(),
            4.0 * area / 3.0,
        );
        near(
            field.nodal_boundary_current_a()[7].unwrap(),
            2.0 * area / 3.0,
        );
        for triangle in 0..4 {
            let density = if triangle < 2 { 4.0 } else { 12.0 };
            near(
                field.cell_power_w()[triangle],
                density * original.mesh().face_measure(triangle).unwrap(),
            );
        }
    }
}

#[test]
fn direct_terminal_budget_trace_and_cancellation_refusals_preserve_the_owner() {
    let original = PlanarConduction::new(
        disconnected_layers(Metric2::planar(1.0).unwrap(), 2),
        conductors(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let currents = [-2.0, 0.0, 2.0];
    let mut calls = 0;
    let (baseline, _) = original
        .solve_terminals(
            &TERMINALS,
            ConductionTerminalControl::Currents(&currents),
            terminal_limits(),
            controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut attempt = 0;
    assert!(matches!(
        original.solve_terminals(
            &TERMINALS,
            ConductionTerminalControl::Currents(&currents),
            terminal_limits(),
            controls(),
            &mut || {
                attempt += 1;
                attempt >= calls - 1
            },
        ),
        Err(ConductionError::Cancelled)
    ));
    for control in [
        ConductionTerminalControl::Currents(&[-2.0, 0.0, 3.0]),
        ConductionTerminalControl::Voltages(&[0.0, 8.0]),
        ConductionTerminalControl::Mixed {
            current_terminals: &[1, 1],
            currents_a: &[0.0, 0.0],
            voltages_v: &[0.0],
        },
    ] {
        assert!(
            original
                .solve_terminals(
                    &TERMINALS,
                    control,
                    terminal_limits(),
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    for bound in [
        ConductionTerminalLimits {
            maximum_terminals: 2,
            ..terminal_limits()
        },
        ConductionTerminalLimits {
            maximum_workspace_values: 1,
            ..terminal_limits()
        },
    ] {
        assert!(
            original
                .solve_terminals(
                    &TERMINALS,
                    ConductionTerminalControl::Currents(&currents),
                    bound,
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let axis = AxisymmetricConduction::new(
        disconnected_layers(Metric2::axisymmetric(2.0).unwrap(), 1),
        conductors()[..2].to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        axis.solve_terminals(
            &[&[0, 3], &[1, 2]],
            ConductionTerminalControl::Voltages(&[0.0, 1.0]),
            terminal_limits(),
            controls(),
            &mut || false,
        )
        .is_err()
    );
    let (retry, _) = original
        .solve_terminals(
            &TERMINALS,
            ConductionTerminalControl::Currents(&currents),
            terminal_limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.voltage_v(), baseline.voltage_v());
    assert_eq!(retry.cell_power_w(), baseline.cell_power_w());
    assert_eq!(original.conductors(), conductors());
}
