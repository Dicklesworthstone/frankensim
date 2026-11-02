use super::*;
use crate::modal_mouth::tests::{air, budget, fixture};
use fs_bem::near_field::FirstOrder;
use fs_duct::modal::{ModeCount, PreparedModalDuct};
use fs_duct::{Duct, LossModel, Segment};

fn native_budget() -> Budget {
    Budget {
        maximum_unknowns: 1000,
        maximum_cache_bytes: 64 << 20,
        maximum_work: 4_000_000_000,
        angular_refinement: 0,
    }
}
fn field_budget() -> ModalMouthObservationBudget {
    ModalMouthObservationBudget {
        maximum_receivers: 8,
        maximum_observation_bytes: 64 << 20,
        maximum_contraction_terms: 100_000_000,
    }
}
fn quadrature() -> Options {
    Options {
        relative_tolerance: 1e-7,
        maximum_depth: 12,
        maximum_kernel_evaluations: 2_000_000,
    }
}
fn close(a: C64, b: C64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-24),
        "{a:?} != {b:?}"
    );
}

#[test]
fn g1_coherent_modal_sound_and_reaction_match_one_original_mixed_source() {
    let (surface, patch) = fixture(12, 3);
    let gate = CancelGate::new();
    let mouth = BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate).unwrap();
    let prepared = mouth
        .prepare_galerkin_field(2400., native_budget(), &gate)
        .unwrap();
    let positions = [[0.08, 0.013, 0.009], [0.1, -0.025, 0.01]];
    let receiver = prepared
        .receivers(&positions, 0.005, field_budget(), &gate)
        .unwrap();
    let response = receiver.solve(2300., quadrature(), &gate).unwrap();
    let work = response.kernel_evaluations();
    assert!(work >= surface.areas().len() * positions.len() * 80);
    assert_eq!(response.source().ports(), 2);
    let flows = [C64::new(2e-6, -0.3e-6), C64::new(-1.1e-6, 0.9e-6)];
    let observed = response.observe(&flows, &gate).unwrap();

    // A separately prepared native boundary receives the ACTUAL combined
    // physical velocity, not columns from the new observation bank.
    let original = Prepared::new(&surface, 2400. / air().sound_speed, native_budget(), || {
        false
    })
    .unwrap();
    let velocity: Vec<_> = (0..surface.areas().len())
        .map(|panel| {
            flows[0].scale(mouth.velocity_trace(0).unwrap()[panel])
                + flows[1].scale(mouth.velocity_trace(1).unwrap()[panel])
        })
        .collect();
    let mixed = original
        .solve_batch(
            2300. / air().sound_speed,
            prepared.medium(),
            &[&velocity],
            || false,
        )
        .unwrap();
    let rows = original
        .receivers(&positions, 0.005)
        .unwrap()
        .prepare(
            2300. / air().sound_speed,
            prepared.medium(),
            quadrature(),
            true,
            || false,
        )
        .unwrap();
    let actual = rows.evaluate_velocity(&mixed[0]).unwrap();
    close(
        observed.complex_mouth_power_w,
        mixed[0].surface_power,
        2e-10,
    );
    for r in 0..positions.len() {
        close(
            observed.field.scalar.pressure[r],
            actual.scalar.pressure[r],
            2e-10,
        );
        for c in 0..3 {
            close(
                observed.field.particle_velocity_m_s[r][c],
                actual.particle_velocity_m_s[r][c],
                2e-10,
            );
            let intensity =
                0.5 * (actual.scalar.pressure[r] * actual.particle_velocity_m_s[r][c].conj()).re;
            assert!(
                (observed.mean_intensity_w_m2[r][c] - intensity).abs()
                    < 2e-10 * intensity.abs().max(1e-24)
            );
        }
        assert!(observed.field.scalar.quadrature_error_estimate_pa[r] >= 0.);
    }
    for mode in 0..2 {
        let direct = sum_complex(
            mixed[0]
                .pressure
                .iter()
                .zip(surface.areas())
                .zip(mouth.velocity_trace(mode).unwrap())
                .map(|((&p, &a), &b)| p.scale(a * b)),
        );
        close(observed.modal_pressure_pa[mode], direct, 2e-10);
    }

    // Physical zero-net-flow modes remain audible; phase changes move both
    // interference and source work, while all boundary/receiver work is reused.
    let higher = response.observe(&[C64::ZERO, flows[1]], &gate).unwrap();
    assert!(higher.field.scalar.pressure[0].abs() > 1e-5);
    let opposite = response.observe(&[flows[0], -flows[1]], &gate).unwrap();
    assert!(
        (opposite.field.scalar.pressure[0] - observed.field.scalar.pressure[0]).abs()
            > 0.1 * higher.field.scalar.pressure[0].abs()
    );
    assert_eq!(response.kernel_evaluations(), work);
    let repeat = response.observe(&flows, &gate).unwrap();
    assert_eq!(repeat.field.scalar.pressure, observed.field.scalar.pressure);
    assert_eq!(repeat.complex_mouth_power_w, observed.complex_mouth_power_w);
}

#[test]
fn g1_loaded_duct_flows_drive_matched_pressure_velocity_and_directional_receiver() {
    let (surface, patch) = fixture(12, 3);
    let gate = CancelGate::new();
    let gas = air();
    let mouth = BemModalMouth::new(&surface, patch, gas, ModeCount(2), budget(), &gate).unwrap();
    let prepared = mouth
        .prepare_galerkin_field(2300., native_budget(), &gate)
        .unwrap();
    let positions = [[0.08, 0.013, 0.009]];
    let receiver = prepared
        .receivers(&positions, 0.005, field_budget(), &gate)
        .unwrap();
    let response = receiver.solve(2300., quadrature(), &gate).unwrap();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            length: 0.06,
            radius: 0.02,
        }],
    };
    let bore = PreparedModalDuct::new(&duct, 2, 1).unwrap();
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.02_f64.powi(2));
    let transfer = bore
        .mouth_waves_matrix(&gas, 2300., LossModel::Lossless, response.load())
        .unwrap();
    let input = [C64::new(1.2, -0.3), C64::new(-0.4, 0.7)];
    let reflected = transfer.excite(&[reference; 2], &input).unwrap();
    let observed = response.observe(&reflected.mouth_flow, &gate).unwrap();
    let input_work: f64 = input
        .iter()
        .zip(&reflected.reflected_pressure)
        .map(|(a, b)| 0.5 * (a.norm_sq() - b.norm_sq()) / reference)
        .sum();
    assert!(
        (input_work - observed.complex_mouth_power_w.re).abs() < 2e-8 * input_work.abs().max(1e-24)
    );
    let p = observed.field.scalar.pressure[0];
    let v = observed.field.particle_velocity_m_s[0];
    let front = FirstOrder::new(0.5, [-1., 0., 0.])
        .unwrap()
        .observe(p, v, prepared.medium())
        .unwrap();
    let back = FirstOrder::new(0.5, [1., 0., 0.])
        .unwrap()
        .observe(p, v, prepared.medium())
        .unwrap();
    close(front + back, p, 2e-14);
    assert!((front - back).abs() > 0.01 * p.abs());
    assert!(response.unit_fields()[1].scalar.pressure[0].abs() > 1.0);
}

#[test]
fn g4_original_band_receiver_budget_and_cancel_refusals_preserve_reusable_source() {
    let (surface, patch) = fixture(12, 3);
    let gate = CancelGate::new();
    let mouth = BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate).unwrap();
    assert!(
        mouth
            .prepare_galerkin_field(f64::NAN, native_budget(), &gate)
            .is_err()
    );
    let prepared = mouth
        .prepare_galerkin_field(2300., native_budget(), &gate)
        .unwrap();
    let positions = [[0.08, 0.013, 0.009]];
    assert!(
        prepared
            .receivers(
                &positions,
                0.005,
                ModalMouthObservationBudget {
                    maximum_observation_bytes: 1,
                    ..field_budget()
                },
                &gate
            )
            .is_err()
    );
    assert!(
        prepared
            .receivers(
                &positions,
                0.005,
                ModalMouthObservationBudget {
                    maximum_contraction_terms: 1,
                    ..field_budget()
                },
                &gate
            )
            .is_err()
    );
    assert!(
        prepared
            .receivers(&[[0.02, 0., 0.]], 0.005, field_budget(), &gate)
            .is_err()
    );
    let receiver = prepared
        .receivers(&positions, 0.005, field_budget(), &gate)
        .unwrap();
    assert!(receiver.solve(2301., quadrature(), &gate).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        receiver.solve(2300., quadrature(), &cancelled),
        Err(ModalMouthObservationError::Cancelled)
    ));
    let response = receiver.solve(2300., quadrature(), &gate).unwrap();
    let flows = [C64::from_re(1e-6), C64::new(0., 0.2e-6)];
    let before = response.observe(&flows, &gate).unwrap();
    assert!(response.observe(&flows[..1], &gate).is_err());
    assert!(
        response
            .observe(&[C64::new(f64::NAN, 0.); 2], &gate)
            .is_err()
    );
    assert!(matches!(
        response.observe(&flows, &cancelled),
        Err(ModalMouthObservationError::Cancelled)
    ));
    let after = response.observe(&flows, &gate).unwrap();
    assert_eq!(before.field.scalar.pressure, after.field.scalar.pressure);
    assert_eq!(
        before.field.particle_velocity_m_s,
        after.field.particle_velocity_m_s
    );
    assert_eq!(before.complex_mouth_power_w, after.complex_mouth_power_w);
}
