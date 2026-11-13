//! G1/G3/G4: canonical magnetic histories reach the original material fields.
#[path = "../../../fs-em/tests/common/mod.rs"]
mod common;
#[path = "scenario_magnetic/vector.rs"]
mod vector;

use common::{material, near};
use fs_couple::scenario_magnetic::{
    MagneticFluxBoundary3, PositiveMagneticFlux, ScenarioMagneticError, magnetic_flux_boundary_at,
    solve_magnetic_flux3,
};
use fs_em::spatial::{SpatialMagneticControls, SpatialMagneticMaterial, SpatialMagnetostatics};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, ScalarPayload, TableInterpolation,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Environment, FrameTree, Physics, Scenario, WORLD,
};

const TESLA: Dims = Dims([0, 1, -2, 0, -1, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

fn problem() -> SpatialMagnetostatics {
    SpatialMagnetostatics::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![[0, 1, 2, 3]],
        vec![SpatialMagneticMaterial::Linear(material(2.0))],
        common::limits(),
        &mut || false,
    )
    .unwrap()
}

fn controls() -> SpatialMagneticControls {
    SpatialMagneticControls {
        cg: common::config(),
        flux_tolerance_wb: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}

fn flux_row(value: f64, parity: OrientationParity) -> BoundaryCondition {
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(TESLA),
        PayloadId::new("basis/scalar").unwrap(),
        WORLD,
        parity,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    BoundaryCondition {
        region: "magnetic-exterior".into(),
        physics: Physics::Magnetics,
        kind: BcKind::NormalMagneticFluxDensity,
        value: Some(BcValue::Typed(Payload::Scalar(
            ScalarPayload::new(
                meta,
                SampleSource::table(
                    [0.0, 1.0].map(|x| QtyAny::new(x, TIME)).to_vec(),
                    [0.0, value].map(|x| QtyAny::new(x, TESLA)).to_vec(),
                    TableInterpolation::Linear,
                    OutsideDomainPolicy::Refuse,
                )
                .unwrap(),
            )
            .unwrap(),
        ))),
        compatibility: None,
        frame: 0,
    }
}

fn history(problem: &SpatialMagnetostatics) -> Scenario {
    let mut scenario = Scenario::new("normal-flux-history", 73, Environment::earth_lab());
    for face in problem.exterior_faces() {
        let n = problem.exterior_normal(face).unwrap();
        let value = n[0] * 0.3 - n[1] * 0.2 + n[2] * 0.4;
        let mut row = flux_row(value, OrientationParity::Odd);
        row.region = format!("original-face-{face}");
        scenario.base_bcs.push(row);
    }
    scenario
}

#[test]
fn canonical_normal_flux_history_recovers_analytical_field_and_energy() {
    let problem = problem();
    let scenario = history(&problem);
    assert!(scenario.validate().is_empty());
    let decoded = parse_ir(&write_ir(&scenario)).unwrap();
    let restored = decoded.scenario();
    let faces: Vec<_> = problem
        .exterior_faces()
        .into_iter()
        .map(|face| [face])
        .collect();
    let rows: Vec<_> = faces
        .iter()
        .zip(&restored.base_bcs)
        .map(|(faces, row)| MagneticFluxBoundary3 {
            faces,
            row,
            positive: PositiveMagneticFlux::OutOfBody,
        })
        .collect();
    for time in [0.0, 0.25, 0.5, 1.0] {
        let result = solve_magnetic_flux3(
            &problem,
            &rows,
            &restored.frames,
            time,
            controls(),
            &mut || false,
        )
        .unwrap();
        let b = [0.3 * time, -0.2 * time, 0.4 * time];
        assert_eq!(result.time_s, time);
        assert_eq!(result.boundary.flux.len(), 4);
        assert_eq!(result.field.gauge_vertices().len(), 1);
        for axis in 0..3 {
            near(result.field.flux_density_t()[0][axis], b[axis]);
            near(result.field.field_strength_a_m()[0][axis], b[axis] / 2.0);
        }
        near(
            result.field.stored_energy_j(),
            b.iter().map(|x| x * x).sum::<f64>() / 24.0,
        );
        assert!(result.field.maximum_weak_residual_wb() <= controls().flux_tolerance_wb);
    }
    // Only the authored normal convention changes sign; geometry already owns
    // the actual outward normal, independently of face winding or frame axes.
    let inward_rows: Vec<_> = rows
        .iter()
        .map(|row| MagneticFluxBoundary3 {
            positive: PositiveMagneticFlux::IntoBody,
            ..*row
        })
        .collect();
    let reversed = solve_magnetic_flux3(
        &problem,
        &inward_rows,
        &restored.frames,
        0.5,
        controls(),
        &mut || false,
    )
    .unwrap();
    for (actual, expected) in reversed.field.flux_density_t()[0]
        .iter()
        .zip([-0.15, 0.1, -0.2])
    {
        near(*actual, expected);
    }
}

#[test]
fn flux_compatibility_metadata_and_cancellation_keep_the_native_problem_reusable() {
    let problem = problem();
    let scenario = history(&problem);
    let faces: Vec<_> = problem
        .exterior_faces()
        .into_iter()
        .map(|face| [face])
        .collect();
    let rows: Vec<_> = faces
        .iter()
        .zip(&scenario.base_bcs)
        .map(|(faces, row)| MagneticFluxBoundary3 {
            faces,
            row,
            positive: PositiveMagneticFlux::OutOfBody,
        })
        .collect();
    // A gauge must not supply the omitted physical flux.
    assert!(matches!(
        solve_magnetic_flux3(
            &problem,
            &rows[1..],
            &scenario.frames,
            0.5,
            controls(),
            &mut || false
        ),
        Err(ScenarioMagneticError::Native(_))
    ));
    let duplicate = [rows[0], rows[0]];
    assert!(
        magnetic_flux_boundary_at(&problem, &duplicate, &scenario.frames, 0.5, &mut || false)
            .is_err()
    );
    let invalid = [MagneticFluxBoundary3 {
        faces: &[usize::MAX],
        ..rows[0]
    }];
    assert!(
        magnetic_flux_boundary_at(&problem, &invalid, &scenario.frames, 0.5, &mut || false)
            .is_err()
    );
    let even = flux_row(1.0, OrientationParity::Even);
    let invalid = [MagneticFluxBoundary3 {
        row: &even,
        ..rows[0]
    }];
    assert!(
        magnetic_flux_boundary_at(&problem, &invalid, &scenario.frames, 0.5, &mut || false)
            .is_err()
    );
    let mut mismatched_frame = scenario.base_bcs[0].clone();
    mismatched_frame.frame = 12;
    let invalid = [MagneticFluxBoundary3 {
        row: &mismatched_frame,
        ..rows[0]
    }];
    assert!(
        magnetic_flux_boundary_at(&problem, &invalid, &scenario.frames, 0.5, &mut || false)
            .is_err()
    );
    assert!(
        magnetic_flux_boundary_at(&problem, &rows, &scenario.frames, 2.0, &mut || false).is_err()
    );
    assert!(
        magnetic_flux_boundary_at(&problem, &[], &FrameTree::new(), f64::NAN, &mut || false)
            .is_err()
    );
    for stop in [1, 10, 25] {
        let mut polls = 0;
        assert!(matches!(
            solve_magnetic_flux3(
                &problem,
                &rows,
                &scenario.frames,
                0.5,
                controls(),
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(ScenarioMagneticError::Cancelled)
        ));
    }
    let result = solve_magnetic_flux3(
        &problem,
        &rows,
        &scenario.frames,
        0.5,
        controls(),
        &mut || false,
    )
    .unwrap();
    near(result.field.flux_density_t()[0][0], 0.15);
}
