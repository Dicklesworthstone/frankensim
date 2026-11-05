use super::*;
use crate::compressible::spatial::{State, StepControl};
use fs_atmo::mann::{MANN_TARGET_V1, MannParams};
use fs_atmo::mann_field::MannQuadrature;
use fs_qty::Time;

const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: super::super::BoundaryCondition::Transmissive,
    upper: super::super::BoundaryCondition::Transmissive,
};
const BOUNDARIES: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};
const R: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

fn field() -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        391,
        MannParams {
            length_m: 2.0,
            alpha_eps23: 0.01,
            ..MANN_TARGET_V1
        },
        MannQuadrature {
            n_modes: 12,
            min_kl: 0.1,
            max_kl: 8.0,
        },
        [0.3, -0.2, 0.1],
    )
    .unwrap();
    let mut field = HomogeneousMannField::at_anchor(model, 719, vec![1.0; 12]).unwrap();
    field.advance_to(0).unwrap();
    field
}

fn config() -> MannReservoirConfig3 {
    MannReservoirConfig3 {
        grid_origin_in_field_m: [1.0, 2.0, 3.0],
        field_to_grid: R,
        maximum_mode_samples: 10_000,
    }
}

fn grid() -> CartesianGrid {
    CartesianGrid::new(3, 2, 2, 0.1, 0.2, 0.3).unwrap()
}

fn frozen(source: &HomogeneousMannField) -> MannReservoirSource3<'_> {
    MannReservoirSource3::Frozen {
        source,
        reference_time_s: 0.0,
        convection_velocity_m_s: [0.5, 0.0, 0.0],
    }
}

fn flow(grid: CartesianGrid) -> Euler3d {
    let gas = crate::compressible::tests::gas(0.2);
    let state = State::from_primitive(&gas, 1.0, [0.0; 3], 1.0).unwrap();
    Euler3d::new(
        gas,
        grid,
        vec![state; grid.shape().iter().product()],
        BOUNDARIES,
        || false,
    )
    .unwrap()
}

#[test]
fn frozen_native_faces_keep_original_realization_frame_and_conservation() {
    let source = field();
    let grid = grid();
    let inlet = MannReservoirInflow3::new(
        source.model(),
        grid,
        &[BoundarySide::XMin, BoundarySide::YMax, BoundarySide::ZMin],
        config(),
        MannReservoirObservation3::Point,
        0,
        || false,
    )
    .unwrap();
    let mut flow = flow(grid);
    let sample_before = source.sample([0.7, 0.8, 0.9]).unwrap();
    let baseline = |_, _, _| {
        Ok(ReservoirFaceState {
            density: 1.0,
            velocity: [0.2, -0.1, 0.05],
            pressure_pa: 1.0,
        })
    };
    let report = inlet
        .apply_frozen(
            &mut flow,
            frozen(&source),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            baseline,
            || false,
        )
        .unwrap();
    assert_eq!(report.face_count, 16);
    assert_eq!(report.mode_samples, 16 * 12);
    assert_eq!(report.sampled_time_s, 0.0);
    assert_eq!(inlet.workspace_bytes(), 0);
    for &side in inlet.sides() {
        for face in 0..grid.face_count(side) {
            let p = grid.face_center(side, face).unwrap();
            // Independent quarter-turn map: field=(1+y, 2-x, 3+z).
            let u = source
                .sample([1.0 + p[1], 2.0 - p[0], 3.0 + p[2]])
                .unwrap()
                .u;
            let boundary = flow.reservoir_face(side, face).unwrap();
            assert_eq!(boundary.velocity, [0.2 - u[1], -0.1 + u[0], 0.05 + u[2]]);
            assert_eq!(boundary.density, 1.0);
            assert_eq!(boundary.pressure_pa, 1.0);
        }
    }
    let bytes = flow.reservoir_profile_workspace_bytes();
    let step = flow
        .step(Time::new(0.0001), StepControl::default(), || false)
        .unwrap();
    for residual in step.conservation_residual.array() {
        assert!(residual.abs() < 2e-12);
    }
    let report = inlet
        .apply_frozen(
            &mut flow,
            frozen(&source),
            0.001,
            MannReservoirAction3::Update,
            baseline,
            || false,
        )
        .unwrap();
    assert_eq!(report.sampled_time_s, flow.time().value());
    let p = grid.face_center(BoundarySide::XMin, 0).unwrap();
    let u = source
        .sample([
            1.0 + p[1],
            2.0 - p[0] + 0.5 * flow.time().value(),
            3.0 + p[2],
        ])
        .unwrap()
        .u;
    assert_eq!(
        flow.reservoir_face(BoundarySide::XMin, 0).unwrap().velocity,
        [0.2 - u[1], -0.1 + u[0], 0.05 + u[2]]
    );
    assert_eq!(bytes, flow.reservoir_profile_workspace_bytes());
    assert_eq!(source.tick(), 0);
    assert_eq!(source.sample([0.7, 0.8, 0.9]).unwrap(), sample_before);
}

#[test]
fn rectangular_face_average_matches_spatial_quadrature_without_rms_restoration() {
    let source = field();
    let grid = grid();
    let inlet = MannReservoirInflow3::new(
        source.model(),
        grid,
        &[BoundarySide::XMin, BoundarySide::XMax],
        config(),
        MannReservoirObservation3::FaceAverage,
        12 * 8,
        || false,
    )
    .unwrap();
    assert_eq!(inlet.workspace_bytes(), 12 * 8); // Opposite faces share gains.
    let trace = inlet
        .trace(&grid, BOUNDARIES, 0.0, 0.001, frozen(&source))
        .unwrap();
    let center = grid.face_center(BoundarySide::XMin, 1).unwrap();
    let mut independent = [0.0; 3];
    let n = 80;
    for a in 0..n {
        for b in 0..n {
            let y = center[1] + ((a as f64 + 0.5) / n as f64 - 0.5) * 0.2;
            let z = center[2] + ((b as f64 + 0.5) / n as f64 - 0.5) * 0.3;
            let u = source
                .sample([1.0 + y, 2.0 - center[0], 3.0 + z])
                .unwrap()
                .u;
            for (sum, value) in independent.iter_mut().zip([-u[1], u[0], u[2]]) {
                *sum += value / (n * n) as f64;
            }
        }
    }
    let actual = trace.fluctuation(BoundarySide::XMin, 1).unwrap();
    for i in 0..3 {
        assert!((actual[i] - independent[i]).abs() < 2e-7);
    }
    // Original full-tensor covariance remains queryable using the SAME filter.
    let filter = inlet.spatial_filter(BoundarySide::XMin).unwrap();
    assert_eq!(
        filter.half_edges_m(),
        [[0.0; 3], [0.1, 0.0, 0.0], [0.0, 0.0, 0.15]]
    );
}

#[test]
fn mixture_preserves_original_species_axis_and_late_refusal_is_atomic() {
    let source = field();
    let grid = grid();
    let gas = crate::compressible::mixture::tests::standard(0.2);
    let count: usize = grid.shape().iter().product();
    let mut states = vec![0.0; count * 7];
    for row in states.chunks_exact_mut(7) {
        super::super::mixture::state_into(&gas, &[0.6, 0.4, 0.0], [0.0; 3], 1.0, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, states, BOUNDARIES, bytes, || false).unwrap();
    let inlet = MannReservoirInflow3::new(
        source.model(),
        grid,
        &[BoundarySide::XMin],
        config(),
        MannReservoirObservation3::Point,
        0,
        || false,
    )
    .unwrap();
    let mean = |_, face: usize, _: [f64; 3], rho: &mut [f64]| {
        let first = 0.3 + 0.05 * face as f64;
        rho.copy_from_slice(&[first, 1.0 - first, 0.0]);
        Ok(MixtureReservoirFaceState {
            velocity: [0.2, 0.0, 0.0],
            pressure_pa: 1.0,
        })
    };
    inlet
        .apply_mixture(
            &mut flow,
            frozen(&source),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            mean,
            || false,
        )
        .unwrap();
    let saved: Vec<_> = (0..4)
        .map(|face| {
            let row = flow.reservoir_face(BoundarySide::XMin, face).unwrap();
            assert_eq!(row.partial_densities[2], 0.0);
            assert_eq!(row.partial_densities[0], 0.3 + 0.05 * face as f64);
            (
                row.partial_densities.to_vec(),
                row.velocity,
                row.pressure_pa,
            )
        })
        .collect();
    let bad = inlet.apply_mixture(
        &mut flow,
        frozen(&source),
        0.001,
        MannReservoirAction3::Update,
        |side, face, point, rho| {
            let mut value = mean(side, face, point, rho)?;
            value.pressure_pa = if face == 3 { -1.0 } else { 2.0 };
            Ok(value)
        },
        || false,
    );
    assert!(bad.is_err());
    for (face, expected) in saved.iter().enumerate() {
        let row = flow.reservoir_face(BoundarySide::XMin, face).unwrap();
        assert_eq!(
            (row.partial_densities, row.velocity, row.pressure_pa),
            (expected.0.as_slice(), expected.1, expected.2)
        );
    }
    let report = flow
        .step(Time::new(0.0001), StepControl::default(), || false)
        .unwrap();
    for residual in report.conservation_residual.species_masses {
        assert!(residual.abs() < 2e-12);
    }
    for residual in report.conservation_residual.momentum {
        assert!(residual.abs() < 2e-12);
    }
    assert!(report.conservation_residual.energy.abs() < 2e-12);
}

#[test]
fn source_time_band_capacity_and_mismatched_model_refuse_before_publication() {
    let source = field();
    let grid = grid();
    let sides = [BoundarySide::XMin];
    let mut small_work = config();
    small_work.maximum_mode_samples = 47;
    assert!(matches!(
        MannReservoirInflow3::new(
            source.model(),
            grid,
            &sides,
            small_work,
            MannReservoirObservation3::Point,
            0,
            || false
        ),
        Err(MannReservoirError3::Work {
            required: 48,
            maximum: 47
        })
    ));
    assert!(
        MannReservoirInflow3::new(
            source.model(),
            grid,
            &sides,
            config(),
            MannReservoirObservation3::FaceAverage,
            95,
            || false
        )
        .is_err()
    );
    let coarse = CartesianGrid::new(2, 2, 2, 1.0, 0.2, 0.3).unwrap();
    assert!(matches!(
        MannReservoirInflow3::new(
            source.model(),
            coarse,
            &sides,
            config(),
            MannReservoirObservation3::Point,
            0,
            || false
        ),
        Err(MannReservoirError3::SpatialNyquist { axis: 0, .. })
    ));
    let inlet = MannReservoirInflow3::new(
        source.model(),
        grid,
        &sides,
        config(),
        MannReservoirObservation3::Point,
        0,
        || false,
    )
    .unwrap();
    assert!(matches!(
        inlet.trace(&grid, BOUNDARIES, 0.0, 10.0, frozen(&source)),
        Err(MannReservoirError3::AdvectiveNyquist { .. })
    ));
    let periodic = Boundaries3d {
        y: AxisBoundaries::Periodic,
        ..BOUNDARIES
    };
    assert!(
        inlet
            .trace(&grid, periodic, 0.0, 0.001, frozen(&source))
            .is_err()
    );
    let evolving = FixedStepMannField::from_snapshot(source.clone(), 119, 0.001, 10_000).unwrap();
    assert!(
        inlet
            .trace(
                &grid,
                BOUNDARIES,
                0.001,
                0.001,
                MannReservoirSource3::Evolving(&evolving)
            )
            .is_err()
    );
    assert!(
        inlet
            .trace(
                &grid,
                BOUNDARIES,
                0.0,
                0.001,
                MannReservoirSource3::Evolving(&evolving)
            )
            .is_ok()
    );
    let mut flow = flow(grid);
    let baseline = |_, _, _| {
        Ok(ReservoirFaceState {
            density: 1.0,
            velocity: [0.0; 3],
            pressure_pa: 1.0,
        })
    };
    assert!(
        inlet
            .apply_frozen(
                &mut flow,
                frozen(&source),
                0.001,
                MannReservoirAction3::Prepare {
                    maximum_profile_bytes: 100_000
                },
                baseline,
                || true
            )
            .is_err()
    );
    assert_eq!(flow.prepared_reservoir_profile_sides(), [false; 6]);
    inlet
        .apply_frozen(
            &mut flow,
            frozen(&source),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            baseline,
            || false,
        )
        .unwrap();
    let before = flow.reservoir_face(BoundarySide::XMin, 0).unwrap();
    let model = MannSpectralModel::build(
        392,
        source.model().parameters(),
        source.model().quadrature(),
        [0.3, -0.2, 0.1],
    )
    .unwrap();
    let different = HomogeneousMannField::at_anchor(model, 719, vec![1.0; 12]).unwrap();
    assert!(
        inlet
            .apply_frozen(
                &mut flow,
                frozen(&different),
                0.001,
                MannReservoirAction3::Update,
                baseline,
                || false
            )
            .is_err()
    );
    assert_eq!(flow.reservoir_face(BoundarySide::XMin, 0).unwrap(), before);
}
