//! G1/G3/G4: actual spatial current, local heat and endpoint material feedback.
//! Synthetic source laws verify the coupled model, not empirical metal data.
#[path = "circuit_enthalpy/common.rs"]
mod common;

use common::{config, initial, material, material_with_electrical_basis, near, point};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve,
    conduction::{ConductionDrive, ConductionRegion, ConductionThermalNetwork},
    material::{ConductionLink, MaterialCell},
};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::{
    conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis},
    phase::ResolvedHeatCapacityCurve,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use fs_thermal::enthalpy::{CellId, MaterialId};

const REGIONS: [ConductionRegion; 4] = [
    ConductionRegion {
        triangle: 3,
        cell: CellId(10),
    },
    ConductionRegion {
        triangle: 0,
        cell: CellId(30),
    },
    ConductionRegion {
        triangle: 2,
        cell: CellId(10),
    },
    ConductionRegion {
        triangle: 1,
        cell: CellId(30),
    },
];
const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 16,
        max_dofs: 32,
        max_triplets: 144,
    }
}
fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 100,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn mesh(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "electrothermal layered specimen",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, 1.0],
        ],
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
fn problem(
    cards: &[MaterialCard; 2],
    bases: [ElectricalPropertyBasis; 2],
    t: f64,
    depth: f64,
) -> PlanarConduction {
    let sources = [0, 0, 1, 1]
        .map(|i| {
            resolve_conductor_state_with_basis(
                &cards[i],
                &point(t),
                MaterialPropertySelection::SingleClaimOnly,
                bases[i],
            )
            .unwrap()
        })
        .to_vec();
    PlanarConduction::new(mesh(depth), sources, limits(), &mut || false).unwrap()
}
fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> Vec<MaterialCell<'a>> {
    [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect()
}
fn network<'a>(cells: &'a [MaterialCell<'a>]) -> ConductionThermalNetwork<'a> {
    ConductionThermalNetwork {
        cells,
        links: &LINKS,
        regions: &REGIONS,
        volume_tolerance_m3: 1e-12,
    }
}
fn ends(voltage: f64, offset: f64) -> [(usize, f64); 4] {
    [
        (0, offset),
        (3, offset),
        (2, voltage + offset),
        (5, voltage + offset),
    ]
}

#[test]
fn spatial_series_current_heats_its_actual_cells_and_conducts_between_them() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    // Same Ohmic law from the other independent source basis, not a fallback.
    let (b, cb) = material_with_electrical_basis(
        ElectricalPropertyBasis::Conductivity,
        1.0 / 3.0,
        2.0,
        0.0,
        1.0,
    );
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let original = problem(
        &cards,
        [
            ElectricalPropertyBasis::Resistivity,
            ElectricalPropertyBasis::Conductivity,
        ],
        300.0,
        1.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let next = network(&bindings)
        .step_candidate(
            &thermal,
            0.5,
            &original,
            &ends(8.0, 0.0),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // Independent exact layered-conduction solution: I=8/(1+3)=2 A.
    // Actual cell heat powers are I²R: 4 and 12 W. The 2-cell BE heat system is
    // [1.5,-.5;-.5,2.5] dT=[2,6], with sourced G=1 W/K (not the authored 999).
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + 8.0 / 3.5,
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + 10.0 / 3.5,
    );
    let r = next.report();
    let f = &r.electrical;
    near(r.source.links[0].conductance_w_k, 1.0);
    near(f.voltage_v[1], 2.0);
    near(f.voltage_v[4], 2.0);
    for (triangle, &power) in f.cell_power_w.iter().enumerate() {
        near(power, if triangle < 2 { 2.0 } else { 6.0 });
        near(f.current_density_a_m2[triangle][0], -2.0);
        near(f.current_density_a_m2[triangle][1], 0.0);
        near(
            f.electric_field_v_m[triangle][0],
            if triangle < 2 { -2.0 } else { -6.0 },
        );
    }
    near(
        f.nodal_boundary_current_a[2].unwrap() + f.nodal_boundary_current_a[5].unwrap(),
        2.0,
    );
    near(f.heat_j, 8.0);
    near(f.supplied_work_j, 8.0);
    near(r.deposited_joule_j, 8.0);
    near(r.non_electrical_heat_j, 0.0);
    near(next.thermal().energy_j() - thermal.energy_j(), 8.0);
    assert!(r.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
    assert!(f.linear.as_ref().unwrap().converged_euclidean());
    assert_eq!(
        next.problem().conductors()[2].source_basis(),
        ElectricalPropertyBasis::Conductivity
    );
    for (i, conductor) in next.problem().conductors().iter().enumerate() {
        assert_eq!(
            conductor.material().card_identity(),
            cards[i / 2].content_hash()
        );
        assert_eq!(
            conductor.source_property().answer().receipt.selected,
            original.conductors()[i]
                .source_property()
                .answer()
                .receipt
                .selected
        );
    }
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
}

#[test]
fn nonlinear_source_feedback_and_material_substitution_change_the_continuing_field() {
    let mut results = Vec::new();
    for rho in [1.0, 2.0] {
        let (a, ca) = material(rho, 2.0, 0.005);
        let (b, cb) = material(rho, 2.0, 0.005);
        let cards = [a, b];
        let curves = [ca, cb];
        let thermal = initial(&curves, 999.0, None);
        let original = problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        );
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        let first = layout
            .step_candidate(
                &thermal,
                0.2,
                &original,
                &ends(10.0, 0.0),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        // Identical cells have no intercell heat. C dT=dt V²/[4rho(1+alpha dT)].
        let theta = ((1.0_f64 + 0.005 * 0.2 * 100.0 / (rho * 2.0)).sqrt() - 1.0) / 0.01;
        let t = first.thermal().temperature_k(CellId(30)).unwrap();
        near(t, 300.0 + theta);
        near(first.thermal().temperature_k(CellId(10)).unwrap(), t);
        let final_rho = rho * (1.0 + 0.005 * theta);
        near(
            first.problem().conductors()[0].resistivity_ohm_m(),
            final_rho,
        );
        let current = -first.report().electrical.current_density_a_m2[0][0];
        near(current, 10.0 / (2.0 * final_rho));
        assert!(current < 10.0 / (2.0 * rho));
        let continued = layout
            .step_candidate(
                first.thermal(),
                0.4,
                first.problem(),
                &ends(10.0, 0.0),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let warmer = continued.thermal().temperature_k(CellId(30)).unwrap();
        assert!(warmer > t);
        assert!(-continued.report().electrical.current_density_a_m2[0][0] < current);
        near(
            continued.problem().conductors()[0].resistivity_ohm_m(),
            rho * (1.0 + 0.005 * (warmer - 300.0)),
        );
        near(
            continued.thermal().energy_j() - first.thermal().energy_j(),
            continued.report().electrical.supplied_work_j,
        );
        near(first.report().electrical.start_s, 0.0);
        near(continued.report().electrical.start_s, 0.2);
        assert_eq!(first.thermal().time_s(), 0.2);
        results.push((t, current));
    }
    assert!(results[1].0 < results[0].0 && results[1].1 < results[0].1);
}

#[test]
fn electrode_reversal_and_common_gauge_preserve_heat_and_source_pins() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 0.0, None);
    let original = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        1.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let run = |v, offset| {
        layout
            .step_candidate(
                &thermal,
                0.1,
                &original,
                &ends(v, offset),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap()
    };
    let forward = run(8.0, 0.0);
    let reverse = run(-8.0, 64.0);
    for id in [CellId(30), CellId(10)] {
        near(
            forward.thermal().temperature_k(id).unwrap(),
            reverse.thermal().temperature_k(id).unwrap(),
        );
    }
    for i in 0..4 {
        near(
            forward.report().electrical.cell_power_w[i],
            reverse.report().electrical.cell_power_w[i],
        );
        near(
            forward.report().electrical.current_density_a_m2[i][0],
            -reverse.report().electrical.current_density_a_m2[i][0],
        );
        assert_eq!(
            forward.problem().conductors()[i]
                .source_property()
                .answer()
                .receipt
                .selected,
            reverse.problem().conductors()[i]
                .source_property()
                .answer()
                .receipt
                .selected
        );
    }
}

#[test]
fn source_geometry_route_budget_and_late_cancellation_refuse_without_advancing() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 0.0, None);
    let original = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        1.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let saved = thermal.cells().clone();
    let sources = original.conductors().to_vec();
    let duplicate = [REGIONS[0]; 4];
    for regions in [&REGIONS[..3], &duplicate[..]] {
        let bad = ConductionThermalNetwork {
            regions,
            ..network(&bindings)
        };
        assert!(
            bad.step_candidate(
                &thermal,
                0.1,
                &original,
                &ends(8.0, 0.0),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
        );
    }
    for inconsistent in [
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            301.0,
            1.0,
        ),
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            2.0,
        ),
    ] {
        assert!(
            layout
                .step_candidate(
                    &thermal,
                    0.1,
                    &inconsistent,
                    &ends(8.0, 0.0),
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
    }
    // Unsupported temperature must not silently freeze the original conductor.
    assert!(
        layout
            .step_candidate(
                &thermal,
                10.0,
                &original,
                &ends(100.0, 0.0),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    let success = layout
        .step_candidate(
            &thermal,
            0.1,
            &original,
            &ends(8.0, 0.0),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let limited = CircuitEnthalpySolve {
        maximum_circuit_trials: success.report().circuit_trials - 1,
        ..config()
    };
    assert!(matches!(
        layout.step_candidate(
            &thermal,
            0.1,
            &original,
            &ends(8.0, 0.0),
            limits(),
            field_controls(),
            limited,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut actual = 0;
    assert!(matches!(
        layout.step_candidate(
            &thermal,
            0.1,
            &original,
            &ends(8.0, 0.0),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                actual += 1;
                actual == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.conductors(), &sources);
}

#[test]
fn current_drive_rebuilds_material_voltage_and_positive_thermal_feedback() {
    let mut results = Vec::new();
    for rho in [1.0, 2.0] {
        let (a, ca) = material(rho, 2.0, 0.005);
        let (b, cb) = material(rho, 2.0, 0.005);
        let cards = [a, b];
        let curves = [ca, cb];
        let mut thermal = initial(&curves, 999.0, None);
        let mut field = problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        );
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        let current = 3.0;
        let dt = 0.2;
        let mut theta = 0.0;
        for step in 1..=3 {
            let next = layout
                .step_driven_candidate(
                    &thermal,
                    step as f64 * dt,
                    &field,
                    ConductionDrive::Current {
                        positive: &[2, 5],
                        negative: &[0, 3],
                        current_a: current,
                    },
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            // Independent BE solution: C*(theta1-theta0) =
            // dt*I^2*rho0*(1+alpha*theta1), with C=2 J/K per slab.
            theta = (2.0 * theta + dt * current * current * rho)
                / (2.0 - dt * current * current * rho * 0.005);
            let t = next.thermal().temperature_k(CellId(30)).unwrap();
            near(t, 300.0 + theta);
            near(next.thermal().temperature_k(CellId(10)).unwrap(), t);
            let final_rho = rho * (1.0 + 0.005 * theta);
            let voltage = 2.0 * current * final_rho;
            let report = next.report();
            let electrical = &report.electrical;
            near(electrical.voltage_v[2] - electrical.voltage_v[0], voltage);
            near(
                electrical.nodal_boundary_current_a[2].unwrap()
                    + electrical.nodal_boundary_current_a[5].unwrap(),
                current,
            );
            near(electrical.supplied_work_j, dt * voltage * current);
            near(electrical.heat_j, dt * current * current * 2.0 * final_rho);
            near(report.deposited_joule_j, electrical.heat_j);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                electrical.supplied_work_j,
            );
            assert!(
                voltage > 2.0 * current * rho,
                "heating must change the required source voltage"
            );
            for i in 0..4 {
                near(electrical.current_density_a_m2[i][0], -current);
                near(
                    electrical.cell_power_w[i],
                    0.5 * current * current * final_rho,
                );
                near(
                    next.problem().conductors()[i].resistivity_ohm_m(),
                    final_rho,
                );
                assert_eq!(
                    next.problem().conductors()[i]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                    field.conductors()[i]
                        .source_property()
                        .answer()
                        .receipt
                        .selected
                );
            }
            near(report.total_energy_residual_j, 0.0);
            let (t, p, _) = next.into_parts();
            thermal = t;
            field = p;
        }
        results.push(theta);
    }
    assert!(
        results[1] > results[0],
        "fixed current must heat the higher-resistivity material more"
    );
}

#[test]
fn current_reversal_zero_and_matching_voltage_share_native_heat_equations() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.001);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let field = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        1.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let run = |current| {
        layout
            .step_driven_candidate(
                &thermal,
                0.2,
                &field,
                ConductionDrive::Current {
                    positive: &[2, 5],
                    negative: &[0, 3],
                    current_a: current,
                },
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap()
    };
    let forward = run(2.0);
    let reverse = run(-2.0);
    let zero = run(0.0);
    let voltage = forward.report().electrical.voltage_v[2];
    // Same BE endpoint equations when the voltage source is set to the
    // current-source solution's endpoint voltage; no continuous-time identity.
    let matched = layout
        .step_candidate(
            &thermal,
            0.2,
            &field,
            &ends(voltage, 0.0),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    for id in [CellId(30), CellId(10)] {
        near(
            forward.thermal().temperature_k(id).unwrap(),
            reverse.thermal().temperature_k(id).unwrap(),
        );
        near(
            forward.thermal().temperature_k(id).unwrap(),
            matched.thermal().temperature_k(id).unwrap(),
        );
        near(zero.thermal().temperature_k(id).unwrap(), 300.0);
    }
    near(
        forward.report().electrical.heat_j,
        reverse.report().electrical.heat_j,
    );
    near(
        forward.report().electrical.supplied_work_j,
        reverse.report().electrical.supplied_work_j,
    );
    near(zero.report().deposited_joule_j, 0.0);
    for triangle in 0..4 {
        near(
            forward.report().electrical.current_density_a_m2[triangle][0],
            -reverse.report().electrical.current_density_a_m2[triangle][0],
        );
        near(zero.report().electrical.cell_power_w[triangle], 0.0);
    }
}

#[test]
fn invalid_current_electrodes_and_late_cancellation_return_no_partial_state() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.001);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 1.0, None);
    let field = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        1.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    for current in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            layout
                .step_driven_candidate(
                    &thermal,
                    0.2,
                    &field,
                    ConductionDrive::Current {
                        positive: &[2, 5],
                        negative: &[0, 3],
                        current_a: current
                    },
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
    }
    // Even an unpowered electrode must retain physically admitted trace geometry.
    for (positive, negative) in [(&[2][..], &[0, 3][..]), (&[2, 5][..], &[2, 5][..])] {
        assert!(
            layout
                .step_driven_candidate(
                    &thermal,
                    0.2,
                    &field,
                    ConductionDrive::Current {
                        positive,
                        negative,
                        current_a: 0.0
                    },
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
    }
    let mut polls = 0;
    layout
        .step_driven_candidate(
            &thermal,
            0.2,
            &field,
            ConductionDrive::Current {
                positive: &[2, 5],
                negative: &[0, 3],
                current_a: 2.0,
            },
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        layout.step_driven_candidate(
            &thermal,
            0.2,
            &field,
            ConductionDrive::Current {
                positive: &[2, 5],
                negative: &[0, 3],
                current_a: 2.0
            },
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
    near(field.conductors()[0].resistivity_ohm_m(), 1.0);
}
