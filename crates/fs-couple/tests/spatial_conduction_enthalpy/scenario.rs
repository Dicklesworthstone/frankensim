//! G1/G3/G4: canonical electrical histories reach real current and finite heat.
use super::*;
use fs_couple::scenario_electric::{
    ElectricBoundary3, PositiveCurrent, ScenarioElectricError, conduction_boundary_at,
    solve_conduction3,
};
use fs_qty::QtyAny;
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, ScalarPayload, TableInterpolation,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Environment, FrameTree, Physics, Scenario, WORLD,
};

const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const CURRENT: Dims = Dims([-2, 0, 0, 0, 1, 0]);
fn row(kind: BcKind, values: [f64; 2]) -> BoundaryCondition {
    let dims = if kind == BcKind::ElectricPotential {
        VOLTAGE
    } else {
        CURRENT
    };
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(dims),
        PayloadId::new("basis/world-cartesian").unwrap(),
        WORLD,
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    let value = ScalarPayload::new(
        meta,
        SampleSource::table(
            [0.0, 1.0]
                .map(|v| QtyAny::new(v, Dims([0, 0, 1, 0, 0, 0])))
                .to_vec(),
            values.map(|v| QtyAny::new(v, dims)).to_vec(),
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .unwrap(),
    )
    .unwrap();
    BoundaryCondition {
        region: format!("electrode-{kind:?}"),
        physics: Physics::Electrics,
        kind,
        value: Some(BcValue::Typed(Payload::Scalar(value))),
        compatibility: None,
        frame: 0,
    }
}

#[test]
fn canonical_voltage_history_drives_heterogeneous_current_and_local_joule_heat() {
    let cards = [material(2.0, 1.0, 0.0).0, material(5.0, 2.0, 0.0).0];
    let sources = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..12).map(|i| sources[i / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let [right, left] = electrodes(&problem, 2.0);
    let mut scenario = Scenario::new("electrical-history", 19, Environment::earth_lab());
    scenario.base_bcs = vec![
        row(BcKind::ElectricPotential, [0.0; 2]),
        row(BcKind::ElectricPotential, [0.0, 7.0]),
    ];
    scenario.base_bcs[0].region = "ground".into();
    assert!(scenario.validate().is_empty());
    let decoded = parse_ir(&write_ir(&scenario)).unwrap();
    let restored = decoded.scenario();
    let rows = [
        ElectricBoundary3::Potential {
            faces: &left,
            row: &restored.base_bcs[0],
        },
        ElectricBoundary3::Potential {
            faces: &right,
            row: &restored.base_bcs[1],
        },
    ];
    for time in [0.0, 0.5, 1.0] {
        let result = solve_conduction3(
            &problem,
            &rows,
            &restored.frames,
            time,
            field_controls(),
            &mut || false,
        )
        .unwrap();
        // R=2+5 ohms at area 1 m². The authored voltage is 7t volts.
        near(result.time_s, time);
        near(result.field.power_w, 7.0 * time * time);
        near(result.field.supplied_power_w, result.field.power_w);
        for (i, j) in result.field.current_density_a_m2.iter().enumerate() {
            near(j[0], -time);
            near(j[1], 0.0);
            near(j[2], 0.0);
            near(
                result.field.cell_power_w[i],
                time * time * [2.0, 5.0][i / 6] / 6.0,
            );
        }
    }
    assert!(matches!(
        solve_conduction3(
            &problem,
            &rows,
            &restored.frames,
            2.0,
            field_controls(),
            &mut || false
        ),
        Err(ScenarioElectricError::Sample(_))
    ));
    assert!(matches!(
        solve_conduction3(
            &problem,
            &rows,
            &restored.frames,
            0.5,
            field_controls(),
            &mut || true
        ),
        Err(ScenarioElectricError::Cancelled)
    ));
    let overlap = [rows[0], rows[0]];
    assert!(
        conduction_boundary_at(&problem, &overlap, &restored.frames, 0.5, &mut || false).is_err()
    );
    let mut bad = restored.base_bcs[1].clone();
    bad.frame = 3;
    assert!(
        conduction_boundary_at(
            &problem,
            &[ElectricBoundary3::Potential {
                faces: &right,
                row: &bad
            }],
            &restored.frames,
            0.5,
            &mut || false
        )
        .is_err()
    );
    // Failure and retry do not alter original material queries or geometry.
    near(
        solve_conduction3(
            &problem,
            &rows,
            &restored.frames,
            0.5,
            field_controls(),
            &mut || false,
        )
        .unwrap()
        .field
        .power_w,
        1.75,
    );
}

#[test]
fn canonical_current_history_enters_original_temperature_feedback_once() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let sources = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            &context,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..12).map(|i| sources[i / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let regions = (0..12)
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: cells[tetrahedron / 6].cell,
        })
        .collect::<Vec<_>>();
    let network = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let [right, left] = electrodes(&problem, 2.0);
    let ground = row(BcKind::ElectricPotential, [0.0; 2]);
    let current = row(BcKind::NormalCurrentDensity, [0.0, 4.0]);
    for positive in [PositiveCurrent::IntoBody, PositiveCurrent::OutOfBody] {
        let rows = [
            ElectricBoundary3::Potential {
                faces: &left,
                row: &ground,
            },
            ElectricBoundary3::NormalCurrent {
                faces: &right,
                row: &current,
                positive,
            },
        ];
        // Explicit endpoint sampling holds I=2 A over this BE interval.
        let drive =
            conduction_boundary_at(&problem, &rows, &FrameTree::new(), 0.5, &mut || false).unwrap();
        let candidate = network
            .step_candidate(
                &thermal,
                0.5,
                &problem,
                &drive,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        // Independent BE elimination with the actual rho(T), C=[1,2], G=1.
        let aa = 1.5 - 0.004;
        let bb = 2.5;
        let determinant = aa * bb - 0.25;
        let theta = [
            (2.0 * bb + 3.0) / determinant,
            (6.0 * aa + 1.0) / determinant,
        ];
        near(
            candidate.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta[0],
        );
        near(
            candidate.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + theta[1],
        );
        let report = candidate.report();
        near(report.deposited_joule_j, report.electrical.heat_j);
        near(
            candidate.thermal().energy_j() - thermal.energy_j(),
            report.electrical.supplied_work_j,
        );
        near(report.total_energy_residual_j, 0.0);
        let expected = if positive == PositiveCurrent::IntoBody {
            -2.0
        } else {
            2.0
        };
        near(report.electrical.field.current_density_a_m2[0][0], expected);
    }
    assert_eq!(thermal.time_s(), 0.0);
}
