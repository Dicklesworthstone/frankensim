//! G1/G3/G4: nonlinear terminal constraints at their actual physical amplitude.
//! Manufactured roots independently fix current, local heat and source work.
use super::super::terminals::{TERMINALS, disconnected_layers, terminal_limits};
use super::*;
use fs_em::conduction::{
    ConductionContact,
    axisymmetric::response::AxisymmetricConductionResponse,
    response::{
        ConductionResponseReport,
        terminals::{ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalReport},
    },
};

#[path = "../../common/electrical_interface_response.rs"]
mod interface_response;

fn laws(coordinate: ConductorMagnitude) -> Vec<ScalarConductorResponse> {
    [0, 0, 1, 1]
        .map(|i| {
            let (basis, slope) = if coordinate == ConductorMagnitude::ElectricField {
                (ElectricalPropertyBasis::Conductivity, [0.2, 0.4][i])
            } else {
                (ElectricalPropertyBasis::Resistivity, [0.1, 0.2][i])
            };
            response(basis, coordinate, [1.0, 3.0][i], slope, 20.0)
        })
        .to_vec()
}

fn fields(j: f64, coordinate: ConductorMagnitude) -> [f64; 2] {
    if coordinate == ConductorMagnitude::ElectricField {
        // J=aE+bE^2: use the stable positive root.
        [0, 1].map(|i| {
            let a: f64 = [1.0, 3.0][i];
            let b = [0.2, 0.4][i];
            2.0 * j / (a + (a * a + 4.0 * b * j).sqrt())
        })
    } else {
        [0, 1].map(|i| j * ([1.0, 3.0][i] + [0.1, 0.2][i] * j))
    }
}

fn drive<'a>(control: ConductionTerminalControl<'a>) -> ConductionTerminalDrive<'a> {
    ConductionTerminalDrive {
        terminals: &TERMINALS,
        control,
        limits: terminal_limits(),
    }
}

fn check(
    mesh: &TriComplex2,
    before: &[ResolvedConductorState],
    field: &ConductionResponseReport,
    terminal: &ConductionTerminalReport,
    coordinate: ConductorMagnitude,
    j: f64,
    sign: f64,
    area: f64,
) {
    let e = fields(j, coordinate);
    let voltage = e[0] + e[1];
    for (i, exact) in [0.0, e[0], voltage].into_iter().enumerate() {
        near(terminal.voltage_v[i], sign * exact);
    }
    for (i, exact) in [-j * area, 0.0, j * area].into_iter().enumerate() {
        near(terminal.current_a[i], sign * exact);
    }
    near(terminal.supplied_power_w, area * j * voltage);
    near(terminal.power_residual_w, 0.0);
    assert!(terminal.maximum_current_residual_a <= controls().current_tolerance_a);
    near(field.power_w, area * j * voltage);
    near(field.power_balance_residual_w, 0.0);
    for triangle in 0..4 {
        let expected_e = e[triangle / 2];
        near(field.electric_field_v_m[triangle][0], 0.0);
        near(field.electric_field_v_m[triangle][1], -sign * expected_e);
        near(field.current_density_a_m2[triangle][1], -sign * j);
        near(
            field.cell_power_w[triangle],
            j * expected_e * mesh.face_measure(triangle).unwrap(),
        );
        let sampled = &field.evaluated[triangle].conductor;
        near(
            sampled.material().query().axes()["operating"],
            if coordinate == ConductorMagnitude::ElectricField {
                expected_e
            } else {
                j
            },
        );
        assert_eq!(
            sampled.source_property().answer().receipt.selected,
            before[triangle].source_property().answer().receipt.selected,
        );
        assert!(field.evaluated[triangle].current_density_defect_a_m2[1].abs() < 1e-9);
    }
    let reaction = |vertices: &[usize]| {
        vertices
            .iter()
            .map(|&v| field.nodal_boundary_current_a[v].unwrap())
            .sum::<f64>()
    };
    near(reaction(&[2, 3]), sign * area * j);
    near(reaction(&[4, 5]), -sign * area * j);
    near(reaction(TERMINALS[1]), 0.0);
}

#[test]
fn nonlinear_planar_terminals_close_each_drive_and_continue_at_a_new_amplitude() {
    let coordinate = ConductorMagnitude::ElectricField;
    let area = 0.4;
    let model = PlanarConductionResponse::new(
        disconnected_layers(Metric2::planar(area).unwrap(), 2),
        laws(coordinate),
        limits(),
        &mut || false,
    )
    .unwrap();
    let before = model.conductors().to_vec();
    let e = fields(2.0, coordinate);
    for sign in [-1.0, 1.0] {
        let current = [-sign * 2.0 * area, 0.0, sign * 2.0 * area];
        let voltage = [0.0, sign * e[0], sign * (e[0] + e[1])];
        let endpoints = [0.0, voltage[2]];
        let reversed = [sign * 2.0 * area, 0.0];
        for control in [
            ConductionTerminalControl::Currents(&current),
            ConductionTerminalControl::Voltages(&voltage),
            ConductionTerminalControl::Mixed {
                current_terminals: &[1],
                currents_a: &[0.0],
                voltages_v: &endpoints,
            },
            ConductionTerminalControl::Mixed {
                current_terminals: &[2, 1],
                currents_a: &reversed,
                voltages_v: &[1e250],
            },
        ] {
            let (candidate, terminal) = model
                .solve_terminals_candidate(
                    drive(control),
                    limits(),
                    controls(),
                    policy(),
                    128,
                    &mut || false,
                )
                .unwrap();
            assert!(candidate.report().field_solves > 1);
            check(
                model.mesh(),
                &before,
                candidate.report(),
                &terminal,
                coordinate,
                2.0,
                sign,
                area,
            );
        }
    }
    let (first, _) = model
        .solve_terminals_candidate(
            drive(ConductionTerminalControl::Currents(&[
                -2.0 * area,
                0.0,
                2.0 * area,
            ])),
            limits(),
            controls(),
            policy(),
            128,
            &mut || false,
        )
        .unwrap();
    let (continued, _) = first.into_parts();
    let reordered = [TERMINALS[1], TERMINALS[2], TERMINALS[0]];
    let (next, terminal) = continued
        .solve_terminals_candidate(
            ConductionTerminalDrive {
                terminals: &reordered,
                control: ConductionTerminalControl::Currents(&[0.0, 3.0 * area, -3.0 * area]),
                limits: terminal_limits(),
            },
            limits(),
            controls(),
            policy(),
            128,
            &mut || false,
        )
        .unwrap();
    let e3 = fields(3.0, coordinate);
    for (i, expected) in [0.0, e3[1], -e3[0]].into_iter().enumerate() {
        near(terminal.voltage_v[i], expected);
    }
    near(terminal.supplied_power_w, area * 3.0 * (e3[0] + e3[1]));
    near(next.report().electric_field_v_m[0][1], -e3[0]);
    near(next.report().electric_field_v_m[2][1], -e3[1]);
    assert!((e3[0] + e3[1] - 1.5 * (e[0] + e[1])).abs() > 0.1);
    assert_eq!(model.conductors(), before);
}

#[test]
fn nonlinear_swept_terminals_keep_original_e_and_j_laws_and_sector_heat() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        for angle in [0.7, 2.0, core::f64::consts::TAU] {
            let area = angle / 2.0;
            let model = AxisymmetricConductionResponse::new(
                disconnected_layers(Metric2::axisymmetric(angle).unwrap(), 2),
                laws(coordinate),
                limits(),
                &mut || false,
            )
            .unwrap();
            let before = model.conductors().to_vec();
            let (candidate, terminal) = model
                .solve_terminals_candidate(
                    drive(ConductionTerminalControl::Mixed {
                        current_terminals: &[2, 1],
                        currents_a: &[2.0 * area, 0.0],
                        voltages_v: &[1e250],
                    }),
                    limits(),
                    controls(),
                    policy(),
                    128,
                    &mut || false,
                )
                .unwrap();
            check(
                model.mesh(),
                &before,
                candidate.report(),
                &terminal,
                coordinate,
                2.0,
                1.0,
                area,
            );
            assert!(candidate.report().field_solves > 1);
            assert_eq!(model.conductors(), before);
        }
    }
}

#[test]
fn actual_terminal_fields_stay_inside_source_domains_without_artificial_probe_voltages() {
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.1,
        0.002,
    );
    let model = PlanarConductionResponse::new(
        disconnected_layers(Metric2::planar(1.0).unwrap(), 2),
        vec![law; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    // No unit-voltage basis query is supported. The physical E and every
    // actual source iteration remain below the source's 0.002 V/m upper bound.
    let j: f64 = 0.001;
    let e = 2.0 * j / (1.0 + (1.0 + 0.4 * j).sqrt());
    let (candidate, terminal) = model
        .solve_terminals_candidate(
            drive(ConductionTerminalControl::Currents(&[-j, 0.0, j])),
            limits(),
            controls(),
            policy(),
            128,
            &mut || false,
        )
        .unwrap();
    near(terminal.voltage_v[1], e);
    near(terminal.voltage_v[2], 2.0 * e);
    near(terminal.supplied_power_w, 2.0 * j * e);
    for sampled in &candidate.report().evaluated {
        near(sampled.conductor.material().query().axes()["operating"], e);
    }
}

#[test]
fn nonlinear_terminal_field_budget_and_late_cancel_allow_exact_retry() {
    let model = PlanarConductionResponse::new(
        disconnected_layers(Metric2::planar(1.0).unwrap(), 2),
        laws(ConductorMagnitude::ElectricField),
        limits(),
        &mut || false,
    )
    .unwrap();
    let before = model.conductors().to_vec();
    let currents = [-2.0, 0.0, 2.0];
    let actual = drive(ConductionTerminalControl::Currents(&currents));
    let mut calls = 0;
    let (baseline, terminal) = model
        .solve_terminals_candidate(actual, limits(), controls(), policy(), 128, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let needed = baseline.report().field_solves;
    assert!(needed > 1);
    for budget in [0, 1, needed - 1] {
        assert!(matches!(
            model.solve_terminals_candidate(
                actual,
                limits(),
                controls(),
                policy(),
                budget,
                &mut || false,
            ),
            Err(ConductionResponseError::FieldBudget { .. })
        ));
    }
    let mut attempted = 0;
    assert!(matches!(
        model.solve_terminals_candidate(
            actual,
            limits(),
            controls(),
            policy(),
            needed,
            &mut || {
                attempted += 1;
                attempted >= calls - 1
            },
        ),
        Err(ConductionResponseError::Cancelled)
    ));
    let (retry, retry_terminal) = model
        .solve_terminals_candidate(actual, limits(), controls(), policy(), needed, &mut || {
            false
        })
        .unwrap();
    assert_eq!(retry.report().field_solves, needed);
    assert_eq!(retry.report().voltage_v, baseline.report().voltage_v);
    assert_eq!(retry.report().cell_power_w, baseline.report().cell_power_w);
    assert_eq!(retry_terminal.voltage_v, terminal.voltage_v);
    assert_eq!(
        retry_terminal.supplied_power_w.to_bits(),
        terminal.supplied_power_w.to_bits()
    );
    assert_eq!(model.conductors(), before);
}

#[test]
fn nonlinear_bulk_and_surface_laws_close_before_floating_terminal_heat_is_reported() {
    let mut bulk = laws(ConductorMagnitude::ElectricField);
    let third = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        2.0,
        0.1,
        20.0,
    );
    bulk.extend([third.clone(), third]);
    let law = interface_response::law(bulk[0].original(), bulk[2].original(), 0.15);
    let contact = ConductionContact {
        first_vertices: [2, 3],
        second_vertices: [5, 4],
        source: law.original().clone(),
    };
    let terminals: [&[usize]; 3] = [&[0, 1], &[6, 7, 8, 9], &[10, 11]];
    let j: f64 = 2.0;
    let e = fields(j, ConductorMagnitude::ElectricField);
    let e3 = 2.0 * j / (2.0 + (4.0 + 0.4 * j).sqrt());
    let rc = 0.5 + 0.15 * j;
    let middle = e[0] + e[1] + j * rc;
    let top = middle + e3;
    for angle in [0.7, 2.0] {
        let area = angle / 2.0;
        let model = AxisymmetricConductionResponse::new_with_contact_response(
            disconnected_layers(Metric2::axisymmetric(angle).unwrap(), 3),
            bulk.clone(),
            vec![contact.clone()],
            vec![law.clone()],
            interface_response::controls(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let currents = [-j * area, 0.0, j * area];
        let voltages = [0.0, middle, top];
        for control in [
            ConductionTerminalControl::Currents(&currents),
            ConductionTerminalControl::Voltages(&voltages),
        ] {
            let (next, terminal) = model
                .solve_terminals_candidate(
                    ConductionTerminalDrive {
                        terminals: &terminals,
                        control,
                        limits: terminal_limits(),
                    },
                    limits(),
                    controls(),
                    policy(),
                    128,
                    &mut || false,
                )
                .unwrap();
            let field = next.report();
            near(terminal.voltage_v[1], middle);
            near(terminal.voltage_v[2], top);
            near(terminal.current_a[1], 0.0);
            near(terminal.supplied_power_w, area * j * top);
            near(field.power_w, area * j * top);
            near(field.contact_response[0].current_a, -j * area);
            near(field.contact_response[0].power_w, area * j * j * rc);
            near(
                field.cell_power_w.iter().sum::<f64>() + field.contact_response[0].power_w,
                terminal.supplied_power_w,
            );
            near(next.problem().contacts()[0].source.resistance_ohm_m2(), rc);
            near(
                next.problem().contacts()[0]
                    .source
                    .resolved()
                    .state_coordinate("normal_current")
                    .unwrap(),
                j,
            );
            for (triangle, &power) in field.cell_power_w.iter().enumerate() {
                let expected_e = [e[0], e[1], e3][triangle / 2];
                near(
                    power,
                    area * j
                        * expected_e
                        * if triangle % 2 == 0 {
                            2.0 / 3.0
                        } else {
                            1.0 / 3.0
                        },
                );
            }
            let surface = &next.problem().contacts()[0].source;
            assert_eq!(
                surface.resolved().properties()[0].answer().receipt.selected,
                contact.source.resolved().properties()[0]
                    .answer()
                    .receipt
                    .selected,
            );
            let (continued, _) = next.into_parts();
            let (again, _) = continued
                .solve_terminals_candidate(
                    ConductionTerminalDrive {
                        terminals: &terminals,
                        control: ConductionTerminalControl::Currents(&currents),
                        limits: terminal_limits(),
                    },
                    limits(),
                    controls(),
                    policy(),
                    128,
                    &mut || false,
                )
                .unwrap();
            near(again.report().power_w, area * j * top);
        }
        near(model.contacts()[0].source.resistance_ohm_m2(), 0.5);
    }
}
