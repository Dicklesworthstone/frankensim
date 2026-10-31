use super::*;
use crate::helmholtz::bc::{ForceDipole, Monopole};
use crate::helmholtz::{IterativeOptions, Medium};
use crate::near_field::Options;
use crate::observation_transfers::{
    IterativeBoundaryObservation, IterativeObservationOptions, LocalBoundary, LocalWall,
    PointSourceSpec, ReceiverSpec,
};
use crate::panel3d::SpherePanels;
use core::f64::consts::TAU;
use std::cell::Cell;

fn options() -> IterativeObservationOptions {
    IterativeObservationOptions {
        highest_omega: TAU * 400.0,
        minimum_clearance_m: 1e-3,
        receiver_quadrature: Options::default(),
        solver: IterativeOptions {
            fmm_order: 6,
            leaf_capacity: 16,
            tol: 1e-9,
            restart: 60,
            max_iters: 600,
            ..IterativeOptions::default()
        },
        lower_window: None,
    }
}

/// G2/G3: the outgoing total field agrees with a distant finite receiver,
/// and coherent source cancellation survives the complete angular projector.
#[test]
fn total_far_field_matches_receiver_limit_and_preserves_coherent_interference() {
    let surface = SpherePanels::icosphere(0.025, 2).unwrap();
    let medium = Medium {
        density: 1.2,
        sound_speed: 340.0,
    };
    let gate = CancelGate::new_clock_free();
    // Independent coordinates may drive the same physical location. They
    // must cancel as complex fields, even though each radiates positive power.
    let sources = [
        PointSourceSpec {
            id: "first".into(),
            position_m: [0.02, 0.01, 0.09],
        },
        PointSourceSpec {
            id: "second".into(),
            position_m: [0.02, 0.01, 0.09],
        },
    ];
    let receiver = [ReceiverSpec {
        id: "distant".into(),
        position_m: [0.0, 0.0, 100.0],
    }];
    let boundary = LocalBoundary::try_new_iterative(
        medium,
        vec![LocalWall::Rigid],
        vec![0; surface.areas().len()],
    )
    .unwrap();
    let observer = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium,
        &sources,
        &receiver,
        boundary,
        options(),
        &gate,
    )
    .unwrap();
    let omega = TAU * 300.0;
    let response = observer
        .solve_frequency(omega, 4 * surface.areas().len(), &gate)
        .unwrap();
    let directions = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let required = 2 * directions.len() * (16 * surface.areas().len() + 1);
    let far = response.far_field(&directions, required, &gate).unwrap();
    assert_eq!(far.kernel_evaluations, required);
    assert_eq!(far.amplitudes[..3], far.amplitudes[3..]);
    assert_eq!(far.sources.as_slice(), response.sources());
    let received = response.observe_receivers(&gate).unwrap();
    let phase = -omega / medium.sound_speed * 100.0;
    let from_finite =
        received.transfer[0].scale(100.0) * C64::new(det::cos(phase), det::sin(phase));
    let relative =
        det::sqrt((from_finite - far.amplitudes[0]).norm_sq() / far.amplitudes[0].norm_sq());
    assert!(relative < 2e-3, "finite/outgoing limit: {relative}");

    let projection = IterativeDirectivityOptions {
        l_max: 6,
        angular: DirectivityProjectionOptions {
            maximum_directions: 50_000,
            relative_tolerance: 1e-8,
        },
        maximum_kernel_evaluations: 100_000_000,
        maximum_combined_entries: 2 * surface.areas().len(),
    };
    let quadrature = response
        .coherent_directivity(&[C64::ONE, C64::new(0.0, 1.0)], projection, &gate)
        .unwrap();
    assert!(quadrature.table.captured_fraction > 1.0 - 1e-6);
    for (index, &direction) in directions.iter().enumerate() {
        let expected = far.amplitudes[index] * C64::new(1.0, 1.0);
        let relative = det::sqrt(
            (quadrature.table.evaluate(direction) - expected).norm_sq() / expected.norm_sq(),
        );
        assert!(relative < 2e-4, "coherent SH reconstruction: {relative}");
    }
    let cancelled_field = response
        .coherent_directivity(&[C64::ONE, C64::from_re(-1.0)], projection, &gate)
        .unwrap();
    assert!(
        cancelled_field
            .table
            .coefficients
            .iter()
            .all(|value| value.norm_sq() == 0.0)
    );
    assert_eq!(cancelled_field.table.captured_fraction, 0.0);

    assert!(matches!(
        response.far_field(&directions, required - 1, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Helmholtz(HelmholtzError::FarFieldBudget { .. })
        ))
    ));
    assert!(response.far_field(&[[0.0; 3]], required, &gate).is_err());
    assert!(
        response
            .coherent_directivity(&[C64::ONE], projection, &gate)
            .is_err()
    );
    assert!(
        response
            .coherent_directivity(
                &[C64::ONE; 2],
                IterativeDirectivityOptions {
                    maximum_combined_entries: 2 * surface.areas().len() - 1,
                    ..projection
                },
                &gate
            )
            .is_err()
    );
    assert!(matches!(
        response.coherent_directivity(
            &[C64::ONE; 2],
            IterativeDirectivityOptions {
                maximum_kernel_evaluations: 1,
                ..projection
            },
            &gate
        ),
        Err(IterativeObservationError::Transfer(
            TransferError::Helmholtz(HelmholtzError::FarFieldBudget { .. })
        ))
    ));

    // The checked seam must retain the legacy numerical owner exactly, while
    // stopping inside the first 320-panel direction when its gate is raised.
    let solution = &response.solutions()[0];
    let legacy = helmholtz::far_field_integrated(&surface, solution, medium, &directions);
    let checked = helmholtz::far_field_integrated_cancellable(
        &surface,
        solution,
        medium,
        &directions,
        directions.len() * 16 * surface.areas().len(),
        || false,
    )
    .unwrap();
    assert_eq!(legacy, checked);
    let polls = Cell::new(0);
    assert_eq!(
        helmholtz::far_field_integrated_cancellable(
            &surface,
            solution,
            medium,
            &directions,
            directions.len() * 16 * surface.areas().len(),
            || {
                polls.set(polls.get() + 1);
                polls.get() >= 4
            }
        ),
        Err(HelmholtzError::Cancelled)
    );
    assert_eq!(polls.get(), 4);
    gate.request();
    assert!(matches!(
        response.far_field(&directions, required, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
}

/// A tiny admitted direction-norm discrepancy can otherwise become a large
/// direct/scattered phase error for a source far from the coordinate origin.
#[test]
fn direct_incident_phase_uses_the_boundary_direction_normalization() {
    let gate = CancelGate::new_clock_free();
    let sources = [
        IncidentSource::Monopole(Monopole {
            position_m: [1e9, 0.0, 0.0],
            volume_velocity_m3_s: C64::new(0.3, -0.8),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [1e9, 0.0, 0.0],
            force_on_fluid_n: [C64::new(0.3, -0.8), C64::ONE, C64::ZERO],
        }),
    ];
    for source in sources {
        let mut exact = [C64::ZERO];
        let mut rounded = [C64::ZERO];
        add_compact_sources(
            &mut exact,
            &[[1.0, 0.0, 0.0]],
            &[source],
            6.0,
            2500.0,
            &gate,
        )
        .unwrap();
        add_compact_sources(
            &mut rounded,
            &[[1.0 + 5e-11, 0.0, 0.0]],
            &[source],
            6.0,
            2500.0,
            &gate,
        )
        .unwrap();
        assert_eq!(exact, rounded);
    }
}
