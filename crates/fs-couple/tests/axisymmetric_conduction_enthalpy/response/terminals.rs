//! G1/G3/G4: nonlinear ideal electrodes share the original two-body heat solve.
//! The source laws and independently eliminated balances are manufactured data.
use super::*;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyReport, CircuitEnthalpySolve,
    conduction::response::ConductionResponseThermalStep, material::MaterialNetworkSample,
};
use fs_em::conduction::response::{
    PlanarConductionResponse,
    terminals::{ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits},
};

const TERMINALS: [&[usize]; 3] = [&[0, 1], &[2, 3, 4, 5], &[6, 7]];

pub(super) fn terminal_drive<'a>(
    terminals: &'a [&'a [usize]],
    control: ConductionTerminalControl<'a>,
) -> ConductionTerminalDrive<'a> {
    ConductionTerminalDrive {
        terminals,
        control,
        limits: ConductionTerminalLimits {
            maximum_terminals: 4,
            maximum_workspace_values: 4096,
        },
    }
}

enum Problem {
    Planar(PlanarConductionResponse),
    Axisymmetric(AxisymmetricConductionResponse),
}
type Report = CircuitEnthalpyReport<MaterialNetworkSample, ConductionResponseThermalStep>;
type Endpoint<'t> = (DepositionSimulation<'t>, Problem, Report);

impl Problem {
    fn new(
        cards: &[MaterialCard; 2],
        coordinate: ConductorMagnitude,
        angle: f64,
        planar: bool,
    ) -> Self {
        let mesh = TriComplex2::from_indexed_triangles(
            "source-distinct layers joined only by a floating ideal electrode",
            vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
            vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
            if planar {
                Metric2::planar(angle / 2.0).unwrap()
            } else {
                Metric2::axisymmetric(angle).unwrap()
            },
        )
        .unwrap();
        let sources = bulk_sources(cards, coordinate);
        if planar {
            Self::Planar(
                PlanarConductionResponse::new(mesh, sources, limits(), &mut || false).unwrap(),
            )
        } else {
            Self::Axisymmetric(
                AxisymmetricConductionResponse::new(mesh, sources, limits(), &mut || false)
                    .unwrap(),
            )
        }
    }
    fn mesh(&self) -> &TriComplex2 {
        match self {
            Self::Planar(problem) => problem.mesh(),
            Self::Axisymmetric(problem) => problem.mesh(),
        }
    }
    fn conductors(&self) -> &[ResolvedConductorState] {
        match self {
            Self::Planar(problem) => problem.conductors(),
            Self::Axisymmetric(problem) => problem.conductors(),
        }
    }
    fn step<'t>(
        &self,
        layout: &ConductionThermalNetwork<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: ConductionTerminalDrive<'_>,
        controls: CircuitEnthalpySolve,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Endpoint<'t>, CircuitEnthalpyError<ConductionThermalError>> {
        match self {
            Self::Planar(problem) => layout
                .step_response_terminal_candidate(
                    thermal,
                    target_s,
                    problem,
                    drive,
                    &[],
                    limits(),
                    field_controls(),
                    response_controls(),
                    controls,
                    &[0.0; 2],
                    cancel,
                )
                .map(|candidate| {
                    let (thermal, problem, report) = candidate.into_parts();
                    (thermal, Self::Planar(problem), report)
                }),
            Self::Axisymmetric(problem) => layout
                .step_axisymmetric_response_terminal_candidate(
                    thermal,
                    target_s,
                    problem,
                    drive,
                    &[],
                    limits(),
                    field_controls(),
                    response_controls(),
                    controls,
                    &[0.0; 2],
                    cancel,
                )
                .map(|candidate| {
                    let (thermal, problem, report) = candidate.into_parts();
                    (thermal, Self::Axisymmetric(problem), report)
                }),
        }
    }
}

#[test]
fn nonlinear_floating_electrodes_close_local_heat_and_continue_in_both_geometries() {
    let (j, dt) = (2.0, 0.1);
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let field_law = coordinate == ConductorMagnitude::ElectricField;
        let (a, ca) = nonlinear_material(0, coordinate, None);
        let (b, cb) = nonlinear_material(1, coordinate, None);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = response_point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        // Independent BE balances retain unequal Cp, original G=A and rho(E/J,T).
        let (theta, e) = endpoint(j, dt, [0.0; 2], field_law);
        for (planar, angle) in [
            (false, 0.7),
            (false, 2.0),
            (false, core::f64::consts::TAU),
            (true, 2.0),
        ] {
            let area = angle / 2.0;
            let thermal = initial(&curves, angle);
            let original = Problem::new(&cards, coordinate, angle, planar);
            let currents = [-j * area, 0.0, j * area];
            let voltages = [0.0, e[0], e[0] + e[1]];
            let mixed_currents = [j * area, 0.0];
            let endpoints = [0.0, e[0] + e[1]];
            for control in [
                ConductionTerminalControl::Currents(&currents),
                ConductionTerminalControl::Mixed {
                    current_terminals: &[2, 1],
                    currents_a: &mixed_currents,
                    voltages_v: &[1e250],
                },
                ConductionTerminalControl::Mixed {
                    current_terminals: &[1],
                    currents_a: &[0.0],
                    voltages_v: &endpoints,
                },
                ConductionTerminalControl::Voltages(&voltages),
            ] {
                let (next, problem, report) = original
                    .step(
                        &layout,
                        &thermal,
                        dt,
                        terminal_drive(&TERMINALS, control),
                        config(),
                        &mut || false,
                    )
                    .unwrap();
                let step = &report.electrical;
                let terminals = step.terminals.as_ref().unwrap();
                for terminal in 0..3 {
                    near(terminals.current_a[terminal], currents[terminal]);
                    near(terminals.voltage_v[terminal], voltages[terminal]);
                }
                let face_current = |vertices: &[usize]| {
                    vertices
                        .iter()
                        .map(|&v| step.field.nodal_boundary_current_a[v].unwrap())
                        .sum::<f64>()
                };
                near(face_current(&[2, 3]), j * area);
                near(face_current(&[4, 5]), -j * area);
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(next.temperature_k(cell).unwrap(), 300.0 + theta[i]);
                    for triangle in 2 * i..2 * i + 2 {
                        near(step.field.current_density_a_m2[triangle][0], 0.0);
                        near(step.field.current_density_a_m2[triangle][1], -j);
                        near(step.field.electric_field_v_m[triangle][1], -e[i]);
                        near(
                            step.field.cell_power_w[triangle],
                            j * e[i] * original.mesh().face_measure(triangle).unwrap(),
                        );
                        let source = &problem.conductors()[triangle];
                        near(source.material().query().axes()["T"], 300.0 + theta[i]);
                        near(
                            source.material().query().axes()["operating"],
                            if field_law { e[i] } else { j },
                        );
                        assert_eq!(
                            source.source_property().answer().receipt.selected,
                            original.conductors()[triangle]
                                .source_property()
                                .answer()
                                .receipt
                                .selected
                        );
                    }
                }
                near(step.heat_j, dt * area * j * (e[0] + e[1]));
                near(step.supplied_work_j, dt * terminals.supplied_power_w);
                near(step.supplied_work_j, step.heat_j);
                near(next.energy_j() - thermal.energy_j(), step.heat_j);
                near(report.deposited_joule_j, step.heat_j);
                near(report.source.links[0].conductance_w_k, area);
                assert!(report.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
                assert!(step.field.field_solves > 1);
                assert!(report.circuit_trials > step.field.field_solves);
                // Continue from BOTH accepted owners, now with the middle floating
                // electrode as terminal zero and a different control ordering.
                let reordered = [TERMINALS[1], TERMINALS[2], TERMINALS[0]];
                let reordered_currents = [0.0, j * area, -j * area];
                let drive2 = terminal_drive(
                    &reordered,
                    ConductionTerminalControl::Currents(&reordered_currents),
                );
                let (continued, _, report2) = problem
                    .step(&layout, &next, 2.0 * dt, drive2, config(), &mut || false)
                    .unwrap();
                let (theta2, e2) = endpoint(j, dt, theta, field_law);
                near(
                    continued.temperature_k(CellId(30)).unwrap(),
                    300.0 + theta2[0],
                );
                near(
                    continued.temperature_k(CellId(10)).unwrap(),
                    300.0 + theta2[1],
                );
                let terminal2 = report2.electrical.terminals.as_ref().unwrap();
                near(terminal2.voltage_v[0], 0.0);
                near(terminal2.voltage_v[1], e2[1]);
                near(terminal2.voltage_v[2], -e2[0]);
                near(
                    report2.electrical.supplied_work_j,
                    dt * j * area * (e2[0] + e2[1]),
                );
            }
            assert_eq!(thermal.time_s(), 0.0);
            near(
                original.conductors()[0].material().query().axes()["operating"],
                0.0,
            );
        }
    }
}

#[test]
fn nonlinear_terminal_budget_and_late_cancel_publish_neither_thermal_nor_electrical_state() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = nonlinear_material(0, coordinate, None);
    let (b, cb) = nonlinear_material(1, coordinate, None);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = response_point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    for planar in [false, true] {
        let original = Problem::new(&cards, coordinate, 2.0, planar);
        let thermal = initial(&curves, 2.0);
        let drive = terminal_drive(
            &TERMINALS,
            ConductionTerminalControl::Currents(&[-2.0, 0.0, 2.0]),
        );
        let mut polls = 0;
        let (good, _, report) = original
            .step(&layout, &thermal, 0.1, drive, config(), &mut || {
                polls += 1;
                false
            })
            .unwrap();
        assert!(report.electrical.field.field_solves > 1);
        for cap in [1, report.circuit_trials - 1] {
            let mut controls = config();
            controls.maximum_circuit_trials = cap;
            assert!(matches!(
                original.step(&layout, &thermal, 0.1, drive, controls, &mut || false),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
        }
        let mut late = 0;
        assert!(matches!(
            original.step(&layout, &thermal, 0.1, drive, config(), &mut || {
                late += 1;
                late >= polls - 1
            }),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        let invalid = terminal_drive(
            &TERMINALS,
            ConductionTerminalControl::Currents(&[-2.0, 0.0, 3.0]),
        );
        assert!(
            original
                .step(&layout, &thermal, 0.1, invalid, config(), &mut || false)
                .is_err()
        );
        let mut exact = config();
        exact.maximum_circuit_trials = report.circuit_trials;
        let (retry, _, repeated) = original
            .step(&layout, &thermal, 0.1, drive, exact, &mut || false)
            .unwrap();
        assert_eq!(retry.energy_j().to_bits(), good.energy_j().to_bits());
        assert_eq!(
            repeated.electrical.field.voltage_v,
            report.electrical.field.voltage_v
        );
        assert_eq!(repeated.circuit_trials, report.circuit_trials);
        assert_eq!(thermal.time_s(), 0.0);
        near(
            original.conductors()[0].material().query().axes()["T"],
            300.0,
        );
    }
}
