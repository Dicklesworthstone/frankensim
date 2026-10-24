//! G1/G3/G4: retained spatial induction, endpoint material feedback and the
//! original enthalpy/conduction owner. Coefficients are manufactured test data.
#[path = "winding_enthalpy/common.rs"]
mod common;
#[path = "../../fs-em/tests/winding_ohmic/common.rs"]
mod field;
#[path = "eddy_enthalpy/thermal.rs"]
mod thermal;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    magnetic::{MagneticRegion, MagneticThermalNetwork},
    material::MaterialCell,
};
use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{EddyConductor, EddyControls, PlanarEddy},
};
use fs_matdb::{MaterialCard, PropertyKey};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    phase::ResolvedHeatCapacityCurve,
    state_point::THERMAL_CONDUCTIVITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_thermal::enthalpy::*;

fn field_controls() -> EddyControls {
    EddyControls {
        linear: field::linear(),
        energy_tolerance_j: 1e-7,
        current_balance_tolerance_a: 1e-8,
    }
}
fn sources(problem: &PlanarMagnetostatics) -> Vec<EddyConductor> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, material)| EddyConductor {
            triangle,
            region: 0,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                material.clone(),
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap(),
        })
        .collect()
}
fn square(card: &MaterialCard, depth: f64) -> PlanarMagnetostatics {
    let mesh = TriComplex2::from_indexed_triangles(
        "eddy-enthalpy-single-independent-mode",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    PlanarMagnetostatics::new(
        mesh,
        vec![field::bundle(card, 300.0, ElectricalPropertyBasis::Resistivity); 4],
        field::limits(),
        &mut || false,
    )
    .unwrap()
}
fn initial_field(problem: &PlanarMagnetostatics, mu: f64) -> PlanarEddy<'_> {
    let bc: Vec<_> = (0..4).map(|node| (node, 0.0)).collect();
    // The only free coordinate is A_center=3; K44=4*depth/mu.
    let initial = problem
        .solve(&[36.0 / mu; 4], &bc, field::linear(), &mut || false)
        .unwrap();
    PlanarEddy::new(
        &initial,
        sources(problem),
        &[],
        &bc,
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap()
}
fn single_thermal(curve: &ResolvedHeatCapacityCurve, volume: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(9),
                mass_kg: volume,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: common::boundary(0.0),
            }),
        }],
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn single_routes() -> Vec<MagneticRegion> {
    (0..4)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CellId(7),
        })
        .collect()
}

#[test]
fn same_card_temperature_changes_diffusion_and_deposits_its_spatial_heat_once() {
    let depth = 0.5;
    let mu = 0.1;
    let rho0 = 1.0 / 18.0;
    let dt = 0.01;
    let mut outcomes = Vec::new();
    for alpha in [0.0, 0.005] {
        let card = field::card(rho0, alpha, mu, ElectricalPropertyBasis::Resistivity);
        let curves = common::curves(&[card.clone(), card.clone()]);
        let thermal = single_thermal(&curves[0], depth);
        let problem = square(&card, depth);
        let state = initial_field(&problem, mu);
        let context = field::point(300.0);
        let cells = [MaterialCell {
            cell: CellId(7),
            material_id: MaterialId(9),
            curve: curves[0].curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }];
        let routes = single_routes();
        let layout = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &routes,
            volume_tolerance_m3: 1e-12,
        };
        let next = layout
            .step_eddy_candidate(
                &thermal,
                dt,
                &state,
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let temperature = next.thermal().temperature_k(CellId(7)).unwrap();
        // Independent scalar system: projected P1 mass is sigma*V/18,
        // the magnetic stiffness is 4*V/mu, and Cp*mass is 4*V.
        // Eliminate the field amplitude; solve its monotone thermal balance
        // by a separate bracketed scalar method, not the production Newton.
        let k = 4.0 * depth / mu;
        let amplitude = |t: f64| {
            let c = depth / (18.0 * rho0 * (1.0 + alpha * (t - 300.0)));
            3.0 * (c - dt * k / 2.0) / (c + dt * k / 2.0)
        };
        let mut lo = 300.0;
        let mut hi = 300.0 + 0.5 * k * 9.0 / (4.0 * depth);
        for _ in 0..80 {
            let mid = (lo + hi) / 2.0;
            let a = amplitude(mid);
            if 4.0 * depth * (mid - 300.0) > 0.5 * k * (9.0 - a * a) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        common::near(temperature, (lo + hi) / 2.0);
        common::near(next.field().potential_wb_m()[4], amplitude(temperature));
        let report = next.report();
        common::near(
            report.electrical.heat_j,
            4.0 * depth * (temperature - 300.0),
        );
        common::near(report.deposited_joule_j, report.electrical.heat_j);
        common::near(
            report.electrical.stored_before_j - report.electrical.stored_after_j,
            report.thermal.stored_energy_change_j,
        );
        common::near(report.electrical.supplied_work_j, 0.0);
        assert!(report.total_energy_residual_j.abs() <= common::config().total_energy_tolerance_j);
        for loss in &report.electrical.cells {
            common::near(loss.heat_j, report.electrical.heat_j / 4.0);
            let source = loss.source.conductor.material();
            assert_eq!(source.card_identity(), card.content_hash());
            common::near(
                source
                    .query_point()
                    .iter()
                    .find(|(n, _)| n == "T")
                    .unwrap()
                    .1,
                temperature,
            );
            common::near(
                loss.source.conductor.resistivity_ohm_m(),
                rho0 * (1.0 + alpha * (temperature - 300.0)),
            );
            // Spatial mean J is zero in every face of this mode. Heating
            // remains positive: averaging J before squaring would lose it.
            common::near(loss.current_density_a_m2.iter().sum::<f64>(), 0.0);
            assert!(loss.heat_j > 0.0);
        }
        let continued = layout
            .step_eddy_candidate(
                next.thermal(),
                2.0 * dt,
                next.field(),
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        common::near(
            continued.report().electrical.stored_before_j,
            next.field().energy_j(),
        );
        common::near(
            continued.field().energy_j()
                + 4.0 * depth * (continued.thermal().temperature_k(CellId(7)).unwrap() - 300.0),
            state.energy_j(),
        );
        assert_eq!(continued.field().time_s(), 2.0 * dt);
        assert_eq!(continued.thermal().time_s(), 2.0 * dt);
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        common::near(thermal.temperature_k(CellId(7)).unwrap(), 300.0);
        outcomes.push((temperature, next.field().potential_wb_m()[4]));
    }
    assert!(outcomes[1].0 > outcomes[0].0);
    assert!(outcomes[1].1 < outcomes[0].1);
}

#[test]
fn distributed_induction_conduction_and_authored_sources_share_original_cells_and_clock() {
    let cards = common::cards(1.0, 0.002);
    let curves = common::curves(&cards);
    let thermal = common::initial(&curves, 0.25, None);
    let (problem, profiles, boundary_nodes) =
        field::problem(&cards, ElectricalPropertyBasis::Resistivity, 0.5);
    let bc: Vec<_> = boundary_nodes.iter().map(|&node| (node, 0.0)).collect();
    let current: Vec<_> = profiles[0].iter().map(|j| 20.0 * j).collect();
    let initial = problem
        .solve(&current, &bc, field::linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        sources(&problem),
        &[],
        &bc,
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap();
    let routes: Vec<_> = problem
        .mesh()
        .faces()
        .iter()
        .enumerate()
        .map(|(triangle, face)| {
            let y: f64 = face
                .iter()
                .map(|&v| problem.mesh().vertices()[v as usize][1] / 3.0)
                .sum();
            MagneticRegion::Cell {
                triangle,
                cell: if y < 0.5 { CellId(30) } else { CellId(10) },
            }
        })
        .collect();
    let context = field::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = common::cells(&cards, &curves, &context, &key);
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &common::LINKS,
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let dt = 0.05;
    let zero = vec![0.0; current.len()];
    let next = layout
        .step_eddy_candidate(
            &thermal,
            dt,
            &state,
            &zero,
            common::config(),
            field_controls(),
            &[0.0; 3],
            || false,
        )
        .unwrap();
    let t0 = next.thermal().temperature_k(CellId(30)).unwrap();
    let t1 = next.thermal().temperature_k(CellId(10)).unwrap();
    let report = next.report();
    let mut heat = [0.0; 2];
    for loss in &report.electrical.cells {
        let cell = match routes[loss.source.triangle] {
            MagneticRegion::Cell { cell, .. } => cell,
            _ => unreachable!(),
        };
        let owner = usize::from(cell == CellId(10));
        heat[owner] += loss.heat_j;
        assert_eq!(
            loss.source.conductor.material().card_identity(),
            cards[owner].content_hash()
        );
        let source_t = loss
            .source
            .conductor
            .material()
            .query_point()
            .iter()
            .find(|(name, _)| name == "T")
            .unwrap()
            .1;
        common::near(source_t, if owner == 0 { t0 } else { t1 });
    }
    let g = report.source.network.links[0].conductance_w_k;
    common::near(g, 2.0);
    common::near(t0 - 300.0, heat[0] + dt * g * (t1 - t0));
    common::near(t1 - 300.0, heat[1] + dt * g * (t0 - t1));
    common::near(
        0.4 * (next.thermal().temperature_k(CellId(50)).unwrap() - 300.0),
        dt * 0.7,
    );
    common::near(report.deposited_joule_j, heat.iter().sum());
    common::near(report.non_electrical_heat_j, dt * 0.7);
    common::near(
        report.electrical.stored_after_j - report.electrical.stored_before_j
            + report.thermal.stored_energy_change_j,
        dt * 0.7,
    );
    assert!(report.maximum_cell_residual_j <= common::config().cell_energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= common::config().total_energy_tolerance_j);
    assert_eq!(next.thermal().time_s(), dt);
    assert_eq!(next.field().time_s(), dt);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn geometry_source_clock_budget_and_late_cancellation_preserve_both_histories() {
    let card = field::card(1.0 / 18.0, 0.005, 0.1, ElectricalPropertyBasis::Resistivity);
    let curves = common::curves(&[card.clone(), card.clone()]);
    let thermal = single_thermal(&curves[0], 0.5);
    let problem = square(&card, 0.5);
    let state = initial_field(&problem, 0.1);
    let context = field::point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curves[0].curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let routes = single_routes();
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let wrong_volume = single_thermal(&curves[0], 0.4);
    assert!(
        layout
            .step_eddy_candidate(
                &wrong_volume,
                0.01,
                &state,
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let mut missing = routes.clone();
    missing[0] = MagneticRegion::Fixed { triangle: 0 };
    let incomplete = MagneticThermalNetwork {
        magnetic_regions: &missing,
        ..layout
    };
    assert!(
        incomplete
            .step_eddy_candidate(
                &thermal,
                0.01,
                &state,
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let cloned = curves[0].clone();
    let foreign_cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: cloned.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let foreign = MagneticThermalNetwork {
        cells: &foreign_cells,
        ..layout
    };
    assert!(
        foreign
            .step_eddy_candidate(
                &thermal,
                0.01,
                &state,
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let mut budget = common::config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_eddy_candidate(
            &thermal,
            0.01,
            &state,
            &[0.0; 4],
            budget,
            field_controls(),
            &[0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut calls = 0;
    let expected = layout
        .step_eddy_candidate(
            &thermal,
            0.01,
            &state,
            &[0.0; 4],
            common::config(),
            field_controls(),
            &[0.0],
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        layout.step_eddy_candidate(
            &thermal,
            0.01,
            &state,
            &[0.0; 4],
            common::config(),
            field_controls(),
            &[0.0],
            || {
                seen += 1;
                seen == calls
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        layout
            .step_eddy_candidate(
                &thermal,
                0.02,
                expected.field(),
                &[0.0; 4],
                common::config(),
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let replay = layout
        .step_eddy_candidate(
            &thermal,
            0.01,
            &state,
            &[0.0; 4],
            common::config(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(
        expected.field().potential_wb_m(),
        replay.field().potential_wb_m()
    );
    assert_eq!(
        expected.report().electrical.cells,
        replay.report().electrical.cells
    );
    assert_eq!(
        expected.thermal().cells()[&CellId(7)].specific_h_j_kg,
        replay.thermal().cells()[&CellId(7)].specific_h_j_kg
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(state.time_s(), 0.0);
}

