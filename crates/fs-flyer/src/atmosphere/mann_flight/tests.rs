use super::*;
use fs_atmo::mann::MANN_TARGET_V1;
use fs_atmo::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use std::cell::Cell;

const MODES: usize = 8;
const BUDGET: usize = 3 * 88 * MODES;
const WORK: usize = 8 * MODES;

fn mean() -> FlatSiteLogLaw {
    FlatSiteLogLaw {
        scenario_effective_z0_m: 0.001,
        displacement_height_m: 0.0,
        reference_height_m: 3.0,
        reference_speed_mps: 4.0,
    }
}

fn source() -> FixedStepMannField {
    let model = MannSpectralModel::build(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: MODES,
            min_kl: 0.1,
            max_kl: 12.0,
        },
        [4.0, -0.3, 0.0],
    )
    .unwrap();
    let snapshot = HomogeneousMannField::at_anchor(model, 29, vec![0.04; MODES]).unwrap();
    FixedStepMannField::from_snapshot(snapshot, 91, 0.0079, 88 * MODES).unwrap()
}

fn wind() -> MannFlightWind {
    MannFlightWind::prepare(
        source(),
        mean(),
        MannFlightBoundary::GroundPlane,
        2.5,
        BUDGET,
    )
    .unwrap()
}

// Direct source observation, with an independently written signed tensor map.
fn reference(field: &FixedStepMannField, point: [f64; 3]) -> FlowSample {
    let height = 2.5 - point[2];
    let flow = field
        .sample_wall_blocked(
            [-point[0], point[1], height],
            MannPlaneWall::new(0.0).unwrap(),
        )
        .unwrap();
    let g = flow.grad;
    FlowSample {
        u: [-(flow.u[0] + mean().speed(height)), flow.u[1], -flow.u[2]],
        grad: [
            [g[0][0], -g[0][1], g[0][2] + mean().dspeed_dh(height)],
            [-g[1][0], g[1][1], -g[1][2]],
            [g[2][0], -g[2][1], g[2][2]],
        ],
    }
}

#[test]
fn moving_material_points_share_the_actual_sequential_half_step_path() {
    let mut owner = wind();
    let mut oracle = source();
    let origin = owner.time_s();
    let capacity = owner.workspace_bytes();
    for _ in 0..5 {
        let initial = oracle.clone();
        oracle
            .advance_to(oracle.step() + 1, WORK / 2, || false)
            .unwrap();
        let midpoint = oracle.clone();
        oracle
            .advance_to(oracle.step() + 1, WORK / 2, || false)
            .unwrap();
        let fields = [&initial, &midpoint, &oracle];
        owner
            .with_stages(
                WORK,
                || false,
                |stages| {
                    // Reordered/repeated observations cannot draw fresh amplitudes.
                    for k in [2, 0, 1, 1, 0, 2] {
                        let elapsed = stages.times[k] - origin;
                        let point = [
                            0.7 + 8.0 * elapsed,
                            -0.4 + 0.3 * elapsed,
                            0.2 - 0.4 * elapsed,
                        ];
                        let observed = stages.sample(stages.times[k], point)?;
                        assert_eq!(observed, reference(fields[k], point));
                        assert_eq!(observed, stages.sample(stages.times[k], point)?);
                        assert!(
                            (observed.grad[0][0] + observed.grad[1][1] + observed.grad[2][2]).abs()
                                < 2e-12
                        );
                    }
                    assert_ne!(
                        stages.sample(stages.times[0], [0.7, -0.4, 0.2])?.u,
                        stages.sample(stages.times[1], [0.7, -0.4, 0.2])?.u
                    );
                    Ok(())
                },
            )
            .unwrap();
        assert_eq!(owner.checkpoint().field, oracle.checkpoint());
        assert_eq!(owner.workspace_bytes(), capacity);
    }
}

#[test]
fn world_gradient_and_ground_plane_have_the_same_physical_frame() {
    let owner = wind();
    let point = [0.7, -0.4, 0.2];
    let flow = owner.sample(point).unwrap();
    let epsilon = 1e-5;
    for j in 0..3 {
        let mut positive = point;
        let mut negative = point;
        positive[j] += epsilon;
        negative[j] -= epsilon;
        let up = owner.sample(positive).unwrap().u;
        let um = owner.sample(negative).unwrap().u;
        for i in 0..3 {
            assert!(((up[i] - um[i]) / (2.0 * epsilon) - flow.grad[i][j]).abs() < 2e-7);
        }
    }
    assert_eq!(owner.sample([0.7, -0.4, 2.5]).unwrap().u[2], 0.0);
    assert_eq!(
        owner.sample([0.7, -0.4, 2.5001]).unwrap_err().code,
        "mann-flight-applicability"
    );
}

#[test]
fn failed_stage_queries_work_and_late_cancellation_leave_the_future_path_unchanged() {
    let mut owner = wind();
    let before = owner.checkpoint();
    let called = Cell::new(false);
    assert_eq!(
        owner
            .with_stages(
                WORK - 1,
                || false,
                |_| {
                    called.set(true);
                    Ok(())
                }
            )
            .unwrap_err()
            .code,
        "mann-flight-work"
    );
    assert!(!called.get());
    assert_eq!(owner.checkpoint(), before);

    let mut polls = 0;
    assert!(
        owner
            .with_stages(
                WORK,
                || {
                    polls += 1;
                    polls == 8
                },
                |_| Ok(())
            )
            .is_err()
    );
    assert_eq!(owner.checkpoint(), before);
    for bad_time in [
        f64::NAN,
        owner.time_s() + 0.25 * owner.dt_s(),
        owner.time_s() - owner.dt_s(),
    ] {
        assert_eq!(
            owner
                .with_stages(WORK, || false, |stages| stages.sample(bad_time, [0.0; 3]))
                .unwrap_err()
                .code,
            "mann-flight-stage-time"
        );
        assert_eq!(owner.checkpoint(), before);
    }
    assert_eq!(
        owner
            .with_stages(
                WORK,
                || false,
                |stages| stages.sample(stages.times[1], [0.0, 0.0, 3.0])
            )
            .unwrap_err()
            .code,
        "mann-flight-applicability"
    );
    assert_eq!(owner.checkpoint(), before);

    let finished = Cell::new(false);
    assert_eq!(
        owner
            .with_stages(
                WORK,
                || finished.get(),
                |stages| {
                    stages.sample(stages.times[1], [0.0; 3])?;
                    finished.set(true);
                    Ok(())
                }
            )
            .unwrap_err()
            .code,
        "mann-flight-cancelled"
    );
    assert_eq!(owner.checkpoint(), before);

    let mut restored = MannFlightWind::resume(before, BUDGET).unwrap();
    for _ in 0..3 {
        let observe = |stages: &Stages<'_>| {
            stages
                .times
                .map(|t| stages.sample(t, [0.1, 0.2, -0.3]))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()
        };
        assert_eq!(
            owner.with_stages(WORK, || false, observe).unwrap(),
            restored.with_stages(WORK, || false, observe).unwrap()
        );
        assert_eq!(owner.checkpoint(), restored.checkpoint());
    }
}

#[test]
fn source_and_domain_admission_rejects_invalid_or_underfunded_models() {
    assert_eq!(
        MannFlightWind::prepare(
            source(),
            mean(),
            MannFlightBoundary::GroundPlane,
            0.0,
            BUDGET - 1
        )
        .unwrap_err()
        .code,
        "mann-flight-budget"
    );
    for minimum_height_m in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            MannFlightWind::prepare(
                source(),
                mean(),
                MannFlightBoundary::HomogeneousAbove { minimum_height_m },
                0.0,
                BUDGET
            )
            .unwrap_err()
            .code,
            "mann-flight-invalid"
        );
    }
    let homogeneous = MannFlightWind::prepare(
        source(),
        mean(),
        MannFlightBoundary::HomogeneousAbove {
            minimum_height_m: 1.0,
        },
        2.5,
        BUDGET,
    )
    .unwrap();
    assert!(homogeneous.sample([0.0, 0.0, 1.5]).is_ok());
    assert_eq!(
        homogeneous.sample([0.0, 0.0, 1.5001]).unwrap_err().code,
        "mann-flight-applicability"
    );
    let mut checkpoint = wind().checkpoint();
    checkpoint.version += 1;
    assert_eq!(
        MannFlightWind::resume(checkpoint, BUDGET).unwrap_err().code,
        "mann-flight-checkpoint"
    );
}
