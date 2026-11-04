//! Native moving-current heat changes the original finite material temperatures.
use super::*;
use fs_em::spatial::conduction::motion::SpatialConductionMotion;
use std::cell::Cell;

fn model(
    cards: &[MaterialCard; 2],
    context: &QueryPoint,
    refinement: usize,
) -> (SpatialConduction, Vec<SpatialConductionRegion>) {
    let sources = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            context,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let xs = (0..=2 * refinement)
        .map(|i| i as f64 / refinement as f64)
        .collect::<Vec<_>>();
    let (positions, tets) = slab(&xs);
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..12 * refinement)
            .map(|i| sources[i / (6 * refinement)].clone())
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let regions = (0..12 * refinement)
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: [CellId(30), CellId(10)][tetrahedron / (6 * refinement)],
        })
        .collect();
    (problem, regions)
}

#[test]
fn prescribed_current_motion_keeps_local_heat_and_material_feedback_separate_from_work() {
    for (rho0, refinement, speed) in [
        (1.0, 1, 3.0),
        (1.0, 2, -10.0),
        (2.0, 1, -10.0),
        (2.0, 2, 0.0),
    ] {
        let (a, ca) = material(rho0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let mut thermal = initial(&curves, 999.0, None);
        let (mut problem, regions) = model(&cards, &context, refinement);
        let network = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let velocities = vec![[0.0, speed, 0.0]; regions.len()];
        let fields = vec![[0.0, 0.0, 1.0]; regions.len()];
        let motion = SpatialConductionMotion {
            velocity_m_s: &velocities,
            flux_density_t: &fields,
        };
        let (dt, current) = (0.5, 2.0);
        let mut theta = [0.0; 2];
        for step in 1..=2 {
            let before = thermal.energy_j();
            let drive = boundary(&problem, 2.0, current, None);
            let target = step as f64 * dt;
            let next = network
                .step_moving_candidate(
                    &thermal,
                    target,
                    &problem,
                    &drive,
                    motion,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            // Independent two-cell BE elimination, C=[1,2], G=1, I=2.
            // Prescribed current fixes local I²rho(T) heat, independent of u.
            let aa = 1.5 - 2.0 * rho0 * 0.002;
            let rhs = [theta[0] + 2.0 * rho0, 2.0 * theta[1] + 6.0];
            let det = 2.5 * aa - 0.25;
            theta = [
                (2.5 * rhs[0] + 0.5 * rhs[1]) / det,
                (0.5 * rhs[0] + aa * rhs[1]) / det,
            ];
            let rho = [rho0 * (1.0 + 0.002 * theta[0]), 3.0];
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + theta[0],
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + theta[1],
            );
            let report = next.report();
            let electrical = &report.electrical;
            let field = &electrical.field;
            for tet in 0..regions.len() {
                let resistivity = rho[tet / (6 * refinement)];
                let volume = 1.0 / (6 * refinement) as f64;
                near(field.current_density_a_m2()[tet][0], -current);
                near(field.conductor_field_v_m()[tet][0], -current * resistivity);
                near(
                    field.electric_field_v_m()[tet][0],
                    -current * resistivity - speed,
                );
                near(
                    field.cell_joule_power_w()[tet],
                    current * current * resistivity * volume,
                );
                near(
                    field.cell_mechanical_power_w()[tet],
                    current * speed * volume,
                );
                assert!(
                    problem.conductors()[tet]
                        .has_same_electrical_source(&next.problem().conductors()[tet])
                );
            }
            let voltage = current * (rho[0] + rho[1]) + 2.0 * speed;
            near(electrical.terminal_work_j, dt * current * voltage);
            near(electrical.mechanical_work_j, dt * 2.0 * current * speed);
            near(
                electrical.heat_j,
                dt * current * current * (rho[0] + rho[1]),
            );
            near(electrical.net_supplied_work_j, electrical.heat_j);
            near(next.thermal().energy_j() - before, electrical.heat_j);
            near(report.deposited_joule_j, electrical.heat_j);
            near(report.non_electrical_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            if speed < 0.0 {
                assert!(electrical.terminal_work_j < 0.0);
                assert!(electrical.mechanical_work_j < electrical.terminal_work_j);
            }
            if speed == 0.0 {
                let stationary = network
                    .step_candidate(
                        &thermal,
                        target,
                        &problem,
                        &drive,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                near(stationary.thermal().energy_j(), next.thermal().energy_j());
                near(
                    stationary.report().electrical.supplied_work_j,
                    electrical.terminal_work_j,
                );
                for (&a, &b) in stationary
                    .report()
                    .electrical
                    .field
                    .cell_power_w
                    .iter()
                    .zip(field.cell_joule_power_w())
                {
                    near(a, b);
                }
            }
            let (next_thermal, next_problem, _) = next.into_parts();
            thermal = next_thermal;
            problem = next_problem;
        }
    }
}

#[test]
fn short_circuit_braking_heats_exact_cells_while_open_circuit_screens_motion() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let (problem, regions) = model(&cards, &context, 1);
    let network = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let velocities = vec![[0.0, 2.0, 0.0]; regions.len()];
    let fields = vec![[0.0, 0.0, 1.0]; regions.len()];
    let motion = SpatialConductionMotion {
        velocity_m_s: &velocities,
        flux_density_t: &fields,
    };
    let short = boundary(&problem, 2.0, 0.0, Some(0.0));
    let next = network
        .step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            &short,
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // Motional EMF=4 V, series resistance=4 ohm, J=+1 A/m².
    // Unit-cell powers [1,3] W; native G=1 W/K; C=[1,2] J/K.
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + 4.0 / 7.0,
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + 5.0 / 7.0,
    );
    let report = next.report();
    near(report.electrical.terminal_work_j, 0.0);
    near(report.electrical.mechanical_work_j, -2.0);
    near(report.deposited_joule_j, 2.0);
    near(next.thermal().energy_j() - thermal.energy_j(), 2.0);
    near(report.total_energy_residual_j, 0.0);
    for tet in 0..regions.len() {
        near(report.electrical.field.current_density_a_m2()[tet][0], 1.0);
        near(
            report.electrical.field.cell_joule_power_w()[tet],
            [1.0, 3.0][tet / 6] / 6.0,
        );
    }
    let open = network
        .step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            &SpatialConductionBoundary::default(),
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(open.thermal().energy_j(), thermal.energy_j());
    near(open.report().deposited_joule_j, 0.0);
    near(open.report().electrical.mechanical_work_j, 0.0);
    for e in open.report().electrical.field.electric_field_v_m() {
        near(e[0], -2.0);
    }
}

#[test]
fn moving_heat_keeps_source_budget_and_late_cancellation_transactional() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let before = thermal.energy_j();
    let (problem, regions) = model(&cards, &context, 1);
    let network = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let velocities = vec![[0.0, 2.0, 0.0]; regions.len()];
    let fields = vec![[0.0, 0.0, 1.0]; regions.len()];
    let motion = SpatialConductionMotion {
        velocity_m_s: &velocities,
        flux_density_t: &fields,
    };
    let drive = boundary(&problem, 2.0, 2.0, None);
    let calls = Cell::new(0usize);
    let accepted = network
        .step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    assert!(accepted.thermal().energy_j() > before);
    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        network.step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                calls.get() >= last
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut limited = config();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            motion,
            limits(),
            field_controls(),
            limited,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    // An out-of-chart thermal proposal cannot extrapolate electrical sources.
    assert!(
        network
            .step_moving_candidate(
                &thermal,
                0.5,
                &problem,
                &drive,
                motion,
                limits(),
                field_controls(),
                config(),
                &[10000.0; 2],
                || false,
            )
            .is_err()
    );
    let (stale, _) = model(&cards, &point(400.0), 1);
    assert!(
        network
            .step_moving_candidate(
                &thermal,
                0.5,
                &stale,
                &drive,
                motion,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    near(thermal.time_s(), 0.0);
    near(thermal.energy_j(), before);
    near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
    let unchanged = problem
        .solve_with_motion(&drive, motion, field_controls(), &mut || false)
        .unwrap();
    near(unchanged.joule_power_w(), 16.0);
}
