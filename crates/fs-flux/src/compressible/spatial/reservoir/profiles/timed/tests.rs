use super::*;
use crate::compressible::tests::{close, gas};

const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const DUCT: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};

fn primitive(velocity: [f64; 3]) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0,
        velocity,
        pressure_pa: 1.0,
    }
}

fn field(
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    velocity: [f64; 3],
    sides: &[BoundarySide],
    sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
) -> Euler3d {
    let source = gas(0.4);
    let q = State::from_primitive(&source, 1.0, velocity, 1.0).unwrap();
    let mut flow = Euler3d::new(source, grid, vec![q; grid.len()], boundaries, || false).unwrap();
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(sides)
        .unwrap();
    flow.prepare_reservoir_profiles(sides, bytes, sample, || false)
        .unwrap();
    flow
}

fn balances(report: &StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 8e-12);
    }
}

#[test]
fn g0_constant_stage_source_preserves_native_cells_and_face_quadrature_bits() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let boundaries = Boundaries3d {
        x: OPEN,
        y: OPEN,
        z: OPEN,
    };
    let value = |side: BoundarySide, face: usize| ReservoirFaceState {
        density: 1.0 + 0.01 * face as f64,
        velocity: [0.2, -0.1, 0.03 * side.sign()],
        pressure_pa: 1.0 + 0.02 * face as f64,
    };
    let make = || {
        field(
            grid,
            boundaries,
            [0.2, -0.1, 0.0],
            &BoundarySide::ALL,
            |side, face, _| Ok(value(side, face)),
        )
    };
    let face_count: usize = BoundarySide::ALL
        .into_iter()
        .map(|side| grid.face_count(side))
        .sum();
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut held = make();
        let mut timed = make();
        let bytes = timed.reservoir_profile_workspace_bytes();
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        for _ in 0..2 {
            let start = timed.time().value();
            let mut visits = Vec::new();
            let actual = timed
                .step_with_reservoir_profiles(
                    Time::new(1e-4),
                    controls,
                    2 * face_count,
                    |time, side, face, point| {
                        assert_eq!(point, grid.face_center(side, face).unwrap());
                        visits.push((time, side, face));
                        Ok(value(side, face))
                    },
                    || false,
                )
                .unwrap();
            let expected = held.step(Time::new(1e-4), controls, || false).unwrap();
            assert_eq!(actual.elapsed, expected.elapsed);
            assert_eq!(actual.halvings, 0);
            assert_eq!(actual.hll_fallback_faces, expected.hll_fallback_faces);
            assert_eq!(actual.limited_cells, expected.limited_cells);
            assert_eq!(actual.before, expected.before);
            assert_eq!(actual.after, expected.after);
            assert_eq!(actual.boundary_change, expected.boundary_change);
            assert_eq!(actual.conservation_residual, expected.conservation_residual);
            assert_eq!(timed.states(), held.states());
            assert_eq!(timed.boundary_impulses, held.boundary_impulses);
            let expected_visits: Vec<_> = [start, start + actual.elapsed.value()]
                .into_iter()
                .flat_map(|time| {
                    BoundarySide::ALL.into_iter().flat_map(move |side| {
                        (0..grid.face_count(side)).map(move |face| (time, side, face))
                    })
                })
                .collect();
            assert_eq!(visits, expected_visits);
            assert_eq!(timed.reservoir_profile_workspace_bytes(), bytes);
        }
    }
}

fn varying(time: f64) -> ReservoirFaceState {
    primitive([5.0 + 2.0 * time, 0.3 + 0.4 * time, -0.2 + 0.1 * time])
}

// Analytic integrals of u, u*v and u*v² for affine prescribed velocities.
// The cubic energy flux uses the original gas's referenced internal energy.
fn linear_integral(a: f64, b: f64, time: f64) -> f64 {
    a * time + 0.5 * b * time.powi(2)
}

fn product_integral(a: f64, b: f64, c: f64, d: f64, time: f64) -> f64 {
    a * c * time + 0.5 * (a * d + b * c) * time.powi(2) + b * d * time.powi(3) / 3.0
}

fn cubic_integral(a: f64, b: f64, c: f64, d: f64, time: f64) -> f64 {
    a * c * c * time
        + 0.5 * (2.0 * a * c * d + b * c * c) * time.powi(2)
        + (a * d * d + 2.0 * b * c * d) * time.powi(3) / 3.0
        + 0.25 * b * d * d * time.powi(4)
}

#[test]
fn g1_unsteady_supersonic_inflow_has_second_order_physical_boundary_exchange() {
    let grid = CartesianGrid::new(2, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let end = 0.125;
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut errors = Vec::new();
        for steps in [4, 8, 16] {
            let mut flow = field(
                grid,
                DUCT,
                varying(0.0).velocity,
                &[BoundarySide::XMin],
                |_, _, _| Ok(varying(0.0)),
            );
            let internal = State::from_primitive(flow.gas(), 1.0, [0.0; 3], 1.0)
                .unwrap()
                .energy;
            let u_integral = linear_integral(5.0, 2.0, end);
            let exact = [
                u_integral,
                product_integral(5.0, 2.0, 5.0, 2.0, end) + end,
                product_integral(5.0, 2.0, 0.3, 0.4, end),
                product_integral(5.0, 2.0, -0.2, 0.1, end),
                (internal + 1.0) * u_integral
                    + 0.5
                        * (cubic_integral(5.0, 2.0, 5.0, 2.0, end)
                            + cubic_integral(5.0, 2.0, 0.3, 0.4, end)
                            + cubic_integral(5.0, 2.0, -0.2, 0.1, end)),
            ];
            let mut received = [0.0; 5];
            for _ in 0..steps {
                let report = flow
                    .step_with_reservoir_profiles(
                        Time::new(end / steps as f64),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        2,
                        |time, _, _, _| Ok(varying(time)),
                        || false,
                    )
                    .unwrap();
                assert_eq!(report.halvings, 0);
                balances(&report);
                for (total, outward) in received
                    .iter_mut()
                    .zip(flow.boundary_impulses(BoundarySide::XMin)[0].array())
                {
                    *total -= outward;
                }
                assert_eq!(
                    flow.reservoir_face(BoundarySide::XMin, 0),
                    Some(varying(flow.time().value()))
                );
            }
            assert_eq!(flow.time().value(), end);
            close(received[0], exact[0], 5e-13);
            assert!(flow.states().iter().any(|q| q.momentum[0] > 5.0));
            errors.push(std::array::from_fn::<_, 4, _>(|k| {
                (received[k + 1] - exact[k + 1]).abs()
            }));
        }
        for pair in errors.windows(2) {
            for component in 0..4 {
                let ratio = pair[0][component] / pair[1][component];
                assert!(
                    (3.98..4.02).contains(&ratio),
                    "component={component}, ratio={ratio}"
                );
            }
        }
    }
}

fn retry_field() -> Euler3d {
    let grid = CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
    field(grid, DUCT, [0.0; 3], &[BoundarySide::XMin], |_, _, _| {
        Ok(primitive([-0.2, 0.0, 0.0]))
    })
}

type Snapshot = (
    Vec<State>,
    Time,
    Option<StepReport>,
    Vec<BoundaryImpulse>,
    Option<ReservoirFaceState>,
);
fn snapshot(flow: &Euler3d) -> Snapshot {
    (
        flow.states.clone(),
        flow.time(),
        flow.last_step,
        flow.boundary_impulses.clone(),
        flow.reservoir_face(BoundarySide::XMin, 0),
    )
}

#[test]
fn g1_endpoint_cfl_retries_at_the_shortened_time_with_a_call_wide_budget() {
    let mut flow = retry_field();
    let mut visits = Vec::new();
    let report = flow
        .step_with_reservoir_profiles(
            Time::new(0.05),
            StepControl::default(),
            4,
            |time, _, _, _| {
                visits.push(time);
                Ok(primitive([200.0 * time, 0.0, 0.0]))
            },
            || false,
        )
        .unwrap();
    assert_eq!(visits, [0.0, 0.05, 0.0, 0.025]);
    assert_eq!(report.halvings, 1);
    assert!(!report.first_order_retry);
    assert_eq!(report.elapsed.value(), 0.025);
    assert_eq!(
        flow.reservoir_face(BoundarySide::XMin, 0),
        Some(primitive([5.0, 0.0, 0.0]))
    );
    balances(&report);

    let mut refused = retry_field();
    let before = snapshot(&refused);
    let mut visits = Vec::new();
    assert!(matches!(
        refused.step_with_reservoir_profiles(
            Time::new(0.05),
            StepControl::default(),
            3,
            |time, _, _, _| {
                visits.push(time);
                Ok(primitive([200.0 * time, 0.0, 0.0]))
            },
            || false,
        ),
        Err(Error::Invalid(
            "timed reservoir profile face evaluation budget exhausted"
        ))
    ));
    assert_eq!(visits, [0.0, 0.05, 0.0]);
    assert_eq!(snapshot(&refused), before);
}

fn accepted_field() -> Euler3d {
    let grid = CartesianGrid::new(2, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let mut flow = field(
        grid,
        DUCT,
        varying(0.0).velocity,
        &[BoundarySide::XMin],
        |_, _, _| Ok(varying(0.0)),
    );
    flow.step_with_reservoir_profiles(
        Time::new(0.01),
        StepControl::default(),
        2,
        |time, _, _, _| Ok(varying(time)),
        || false,
    )
    .unwrap();
    flow
}

#[test]
fn g4_endpoint_refusals_and_last_poll_cancellation_preserve_the_accepted_transaction() {
    let mut probe = accepted_field();
    let mut total_polls = 0;
    probe
        .step_with_reservoir_profiles(
            Time::new(0.01),
            StepControl::default(),
            2,
            |time, _, _, _| Ok(varying(time)),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = accepted_field();
    let before = snapshot(&flow);
    let mut polls = 0;
    assert_eq!(
        flow.step_with_reservoir_profiles(
            Time::new(0.01),
            StepControl::default(),
            2,
            |time, _, _, _| Ok(varying(time)),
            || {
                polls += 1;
                polls == total_polls
            }
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(polls, total_polls);
    assert_eq!(snapshot(&flow), before);

    for invalid_source in [false, true] {
        let mut calls = 0;
        let result = flow.step_with_reservoir_profiles(
            Time::new(0.01),
            StepControl::default(),
            20,
            |time, _, _, _| {
                calls += 1;
                if time > before.1.value() {
                    if invalid_source {
                        Ok(ReservoirFaceState {
                            pressure_pa: f64::NAN,
                            ..varying(time)
                        })
                    } else {
                        Err(Error::Invalid("source has no replayable endpoint"))
                    }
                } else {
                    Ok(varying(time))
                }
            },
            || false,
        );
        assert!(result.is_err());
        assert_eq!(calls, 2, "provider/source refusal must not be retried");
        assert_eq!(snapshot(&flow), before);
    }
    // Both source stages and the physical state can still advance after failed
    // attempts left partially overwritten, but strictly private, scratch.
    flow.step_with_reservoir_profiles(
        Time::new(0.01),
        StepControl::default(),
        2,
        |time, _, _, _| Ok(varying(time)),
        || false,
    )
    .unwrap();
    assert_eq!(snapshot(&flow), snapshot(&probe));
}
