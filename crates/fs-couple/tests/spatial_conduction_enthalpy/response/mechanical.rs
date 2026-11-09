//! Independent coupled current, native momentum and finite-heat endpoints.
use super::*;
use fs_couple::circuit_enthalpy::{
    conduction::spatial::{
        mechanical::{
            SpatialConductionKickControls, response::SpatialConductionResponseKickControls,
        },
        response::mechanical::step_spatial_conduction_response_mechanical_candidate,
    },
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_em::{
    conduction::terminals::{ConductionTerminalControl, ConductionTerminalLimits},
    spatial::conduction::{motion::SpatialLorentzLoadLimits, response::ConductionTerminalDrive},
};
use fs_exec::{CancelGate, Cx};
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::FreeFloatingBaseState;
use fs_mbd_material::{
    RegionBinding, ResolvedBody, VolumeMoments, deformation::TetrahedralMass, resolve_body,
};

fn mechanical_controls() -> SpatialConductionResponseKickControls {
    SpatialConductionResponseKickControls {
        physical: SpatialConductionKickControls {
            field: field_controls(),
            loads: SpatialLorentzLoadLimits {
                max_vertices: 100,
                max_cells: 100,
            },
            maximum_field_solves: 10_000,
            angular_velocity_tolerance_rad_s: 1e-8,
            linear_velocity_tolerance_m_s: 1e-8,
            energy_tolerance_j: 1e-7,
        },
        response: controls(),
        maximum_newton_iterations: 16,
        angular_velocity_scale_rad_s: 1.0,
        linear_velocity_scale_m_s: 1.0,
        difference_step: 1e-4,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
            linear_restart: 6,
            ..config().newton
        },
    }
}

fn translating(speed: f64) -> FreeFloatingBaseState {
    FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 0.), Vec3::new(0., speed, 0.)),
    )
}

fn centered_problem(cards: &[MaterialCard; 2]) -> SpatialConductionResponse {
    let original = problem(cards, ConductorMagnitude::ElectricField);
    let positions = original
        .positions_m()
        .iter()
        .map(|p| [p[0] - 1.0, p[1] - 0.5, p[2] - 0.5])
        .collect();
    SpatialConductionResponse::new(
        positions,
        original.mesh().tets.clone(),
        original.sources().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap()
}

fn body(problem: &SpatialConductionResponse, cards: &[MaterialCard; 2]) -> ResolvedBody {
    let context = query(300.0, ConductorMagnitude::ElectricField);
    let keys: Vec<_> = (0..problem.mesh().tets.len())
        .map(|i| format!("original-tet-{i}"))
        .collect();
    let regions: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .map(|(i, tet)| RegionBinding {
            key: &keys[i],
            material: &cards[i / 6],
            point: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::tetrahedron(tet.map(|v| {
                let p = problem.positions_m()[v as usize];
                Vec3::new(p[0], p[1], p[2])
            }))
            .unwrap(),
        })
        .collect();
    resolve_body(&regions, regions.len(), || false).unwrap()
}

fn bind_mass<'a>(
    problem: &SpatialConductionResponse,
    body: &'a ResolvedBody,
) -> TetrahedralMass<'a> {
    let keys: Vec<_> = (0..problem.mesh().tets.len())
        .map(|i| format!("original-tet-{i}"))
        .collect();
    let names: Vec<_> = keys.iter().map(String::as_str).collect();
    let tets: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .map(|tet| tet.map(|v| v as usize))
        .collect();
    TetrahedralMass::bind(body, problem.positions_m(), &tets, &names, 100, || false).unwrap()
}

fn faces(problem: &SpatialConductionResponse) -> [Vec<usize>; 2] {
    [-1.0, 1.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect()
    })
}

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 89,
                kernel_id: 53,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}

#[test]
fn current_momentum_and_each_original_finite_body_close_at_one_endpoint_and_continue() {
    let (a, ca) = nonlinear_material(1., CP[0], ALPHA[0], 0.05, ConductorMagnitude::ElectricField);
    let (b, cb) = nonlinear_material(2., CP[1], ALPHA[1], 0.08, ConductorMagnitude::ElectricField);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = query(300.0, ConductorMagnitude::ElectricField);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let original = centered_problem(&cards);
    let body = body(&original, &cards);
    let mass = bind_mass(&original, &body);
    let ends = faces(&original);
    let terminals = ends.each_ref().map(Vec::as_slice);
    let regions = regions();
    let mut solve = config();
    solve.maximum_circuit_trials = 100_000;
    let j = 2.0;
    let dt = 0.1;
    for control in [
        ConductionTerminalControl::Currents(&[-2.0, 2.0]),
        ConductionTerminalControl::Mixed {
            current_terminals: &[1],
            currents_a: &[2.0],
            voltages_v: &[13.0],
        },
    ] {
        let drive = ConductionTerminalDrive {
            terminals: &terminals,
            control,
            limits: ConductionTerminalLimits {
                maximum_terminals: 2,
                maximum_workspace_values: 20_000,
            },
        };
        for (speed, external_y) in [(0.5, 0.0), (-5.0, 0.6)] {
            with_cx(|cx, _| {
                let mut field = original.clone();
                let mut thermal = initial(&curves, 999.0, None);
                let mut mechanical = translating(speed);
                let mut old = [0.0; 2];
                for step in 1..=2 {
                    let before_heat = thermal.energy_j();
                    let before_speed = mechanical.twist_body.linear.y;
                    let (theta, ep) = endpoint(j, dt, old, [1., 2.], [0.05, 0.08], true);
                    // Uniform J and B give F_y=2*J in total volume two.
                    // Original sourced mass is two, so this is an independent
                    // native momentum and kinetic-work reference.
                    let after = before_speed + dt * (2.0 * j + external_y) / 2.0;
                    let midpoint = 0.5 * (before_speed + after);
                    let voltage = ep[0] + ep[1] + 2.0 * midpoint;
                    let next = step_spatial_conduction_response_mechanical_candidate(
                        &heat,
                        cx,
                        &thermal,
                        step as f64 * dt,
                        &field,
                        &mass,
                        drive,
                        &[[0., 0., 1.]; 12],
                        mechanical,
                        Wrench::new(Vec3::new(0., 0., 0.), Vec3::new(0., external_y, 0.)),
                        &regions,
                        1e-12,
                        limits(),
                        mechanical_controls(),
                        solve,
                        &[0.; 2],
                    )
                    .unwrap();
                    let report = next.report();
                    let native = &report.electrical.field;
                    let kick = &native.field.kick;
                    let electrodes = native.terminals.as_ref().unwrap();
                    near(electrodes.voltage_v[1], voltage);
                    near(electrodes.current_a[1], j);
                    near(next.mechanical().base_after.twist_body.linear.y, after);
                    near(kick.midpoint_twist_body.linear.y, midpoint);
                    let angular = kick.mechanical.base_after.twist_body.angular;
                    for component in [angular.x, angular.y, angular.z] {
                        near(component, 0.0);
                    }
                    near(kick.lorentz_wrench_body.force.y, 2.0 * j);
                    near(
                        kick.mechanical.kinetic_energy_change,
                        (after - before_speed) * (after + before_speed),
                    );
                    near(kick.external_work_j, dt * external_y * midpoint);
                    near(kick.mechanical_work_j, dt * 2.0 * j * midpoint);
                    near(report.electrical.supplied_work_j, dt * j * voltage);
                    near(kick.terminal_work_j, report.electrical.supplied_work_j);
                    near(report.deposited_joule_j, dt * j * (ep[0] + ep[1]));
                    near(
                        report.electrical.net_supplied_work_j,
                        report.deposited_joule_j,
                    );
                    near(
                        next.thermal().energy_j() - before_heat,
                        report.deposited_joule_j,
                    );
                    near(
                        next.thermal().energy_j() - before_heat
                            + kick.mechanical.kinetic_energy_change
                            - report.electrical.supplied_work_j
                            - kick.external_work_j
                            - report.non_electrical_heat_j,
                        report.total_energy_residual_j,
                    );
                    near(report.total_energy_residual_j, 0.0);
                    assert!(native.field_solves > native.iterations);
                    assert!(report.circuit_trials > native.field_solves);
                    assert_eq!(kick.mechanical.base_after.world_from_base, Se3::identity());
                    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                        near(
                            next.thermal().temperature_k(cell).unwrap(),
                            300.0 + theta[i],
                        );
                        let transfer = dt * (theta[1 - i] - theta[i]);
                        near(CP[i] * (theta[i] - old[i]), dt * j * ep[i] + transfer);
                        for tet in 6 * i..6 * i + 6 {
                            near(kick.field.current_density_a_m2()[tet][0], -j);
                            near(kick.field.conductor_field_v_m()[tet][0], -ep[i]);
                            near(kick.field.cell_joule_power_w()[tet], j * ep[i] / 6.0);
                            let source = &native.evaluated[tet].conductor;
                            near(source.material().query().axes()["T"], 300.0 + theta[i]);
                            near(source.material().query().axes()["operating"], ep[i]);
                            assert!(source.has_same_electrical_source(&field.conductors()[tet]));
                        }
                    }
                    if speed < 0.0 {
                        assert!(report.electrical.supplied_work_j < 0.0);
                        assert!(kick.mechanical.kinetic_energy_change < 0.0);
                    } else {
                        assert!(report.electrical.supplied_work_j > 0.0);
                        assert!(kick.mechanical.kinetic_energy_change > 0.0);
                    }
                    mechanical = next.mechanical().base_after;
                    (thermal, field, _) = next.into_parts();
                    old = theta;
                }
            });
        }
    }
    assert_eq!(context.axes()["operating"], 0.0);
    assert!(
        original
            .conductors()
            .iter()
            .all(|c| c.material().query().axes()["T"] == 300.0)
    );
}

#[test]
fn coupled_budget_and_cancellation_preserve_all_three_owners_for_retry() {
    let (a, ca) = nonlinear_material(1., CP[0], ALPHA[0], 0.05, ConductorMagnitude::ElectricField);
    let (b, cb) = nonlinear_material(2., CP[1], ALPHA[1], 0.08, ConductorMagnitude::ElectricField);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = query(300.0, ConductorMagnitude::ElectricField);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let field = centered_problem(&cards);
    let body = body(&field, &cards);
    let mass = bind_mass(&field, &body);
    let ends = faces(&field);
    let terminals = ends.each_ref().map(Vec::as_slice);
    let drive = ConductionTerminalDrive {
        terminals: &terminals,
        control: ConductionTerminalControl::Currents(&[-2., 2.]),
        limits: ConductionTerminalLimits {
            maximum_terminals: 2,
            maximum_workspace_values: 20_000,
        },
    };
    let thermal = initial(&curves, 999.0, None);
    let mechanical = translating(0.5);
    let regions = regions();
    let zero = Wrench::new(Vec3::new(0., 0., 0.), Vec3::new(0., 0., 0.));
    with_cx(|cx, gate| {
        let mut solve = config();
        solve.maximum_circuit_trials = 1;
        let attempt = |controls| {
            step_spatial_conduction_response_mechanical_candidate(
                &heat,
                cx,
                &thermal,
                0.1,
                &field,
                &mass,
                drive,
                &[[0., 0., 1.]; 12],
                mechanical,
                zero,
                &regions,
                1e-12,
                limits(),
                mechanical_controls(),
                controls,
                &[0.; 2],
            )
        };
        assert!(matches!(
            attempt(solve),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(mechanical.twist_body.linear.y, 0.5);
        assert_eq!(
            field.conductors()[0].material().query().axes()["operating"],
            0.0
        );
        solve.maximum_circuit_trials = 100_000;
        let next = attempt(solve).unwrap();
        let (theta, _) = endpoint(2.0, 0.1, [0.; 2], [1., 2.], [0.05, 0.08], true);
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta[0],
        );
        near(next.mechanical().base_after.twist_body.linear.y, 0.7);
        gate.request();
        assert!(matches!(
            attempt(solve),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(field.conductors()[0].material().query().axes()["T"], 300.0);
    });
}
