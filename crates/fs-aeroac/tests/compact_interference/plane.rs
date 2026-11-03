//! G0/G1/G3 half-space boundary, path-interference and physical-source statistics.
use super::*;
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary, PreparedHalfSpaceObserver};
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;

fn wall(boundary: IdealPlaneBoundary) -> IdealAcousticPlane {
    IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 2.0], boundary).unwrap()
}

fn half_space(
    frequency: f64,
    receiver: [f64; 3],
    positions: &[[f64; 3]],
    boundary: IdealPlaneBoundary,
) -> PreparedHalfSpaceObserver {
    observe(1.25, 340.0, frequency, receiver, positions)
        .with_plane(wall(boundary), 2 * positions.len())
        .unwrap()
}

fn source() -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: C64::new(2e-4, -1e-4),
        force_on_fluid_n: [C64::new(0.2, -0.3), C64::new(-0.1, 0.4), C64::new(0.6, 0.2)],
    }
}

#[test]
fn ideal_boundaries_hold_for_every_volume_and_force_component() {
    let positions = [[0.3, -0.2, 0.7]];
    let receiver = [1.5, 0.4, 0.0];
    for frequency in [10.0, 1000.0, 12000.0] {
        for component in 0..4 {
            let mut input = CompactSourcePhasor::default();
            if component == 0 {
                input.volume_flow_m3_s = C64::new(1e-4, 2e-4);
            } else {
                input.force_on_fluid_n[component - 1] = C64::new(0.4, -0.7);
            }
            let direct = observe(1.25, 340.0, frequency, receiver, &positions)
                .evaluate(&[input])
                .unwrap();
            let rigid = half_space(frequency, receiver, &positions, IdealPlaneBoundary::Rigid)
                .evaluate(&[input])
                .unwrap();
            close_complex(rigid.pressure_pa, direct.pressure_pa.scale(2.0));
            for axis in 0..2 {
                close_complex(
                    rigid.particle_velocity_m_s[axis],
                    direct.particle_velocity_m_s[axis].scale(2.0),
                );
            }
            assert_eq!(rigid.particle_velocity_m_s[2], C64::ZERO);
            assert_eq!(rigid.active_intensity_w_m2[2], 0.0);
            let soft = half_space(
                frequency,
                receiver,
                &positions,
                IdealPlaneBoundary::PressureRelease,
            )
            .evaluate(&[input])
            .unwrap();
            assert_eq!(soft.pressure_pa, C64::ZERO);
            assert_eq!(soft.particle_velocity_m_s[..2], [C64::ZERO; 2]);
            close_complex(
                soft.particle_velocity_m_s[2],
                direct.particle_velocity_m_s[2].scale(2.0),
            );
            assert_eq!(soft.active_intensity_w_m2, [0.0; 3]);
        }
    }
}

#[test]
fn reflection_retains_path_delay_spreading_and_frequency_notches() {
    let flow = C64::new(2e-4, -3e-4);
    let input = [CompactSourcePhasor::monopole(flow)];
    for (frequency, ratio) in [(85.0, 0.25), (170.0, 2.25)] {
        let prepared = half_space(
            frequency,
            [0.0, 0.0, 3.0],
            &[[0.0, 0.0, 1.0]],
            IdealPlaneBoundary::Rigid,
        );
        let direct = prepared.direct().evaluate(&input).unwrap();
        let field = prepared.evaluate(&input).unwrap();
        let omega = TAU * frequency;
        let green = phase(omega * 2.0 / 340.0).scale(1.0 / (8.0 * PI))
            + phase(omega * 4.0 / 340.0).scale(1.0 / (16.0 * PI));
        close_complex(
            field.pressure_pa,
            flow * C64::new(0.0, -omega * 1.25) * green,
        );
        close(
            field.pressure_pa.norm_sq() / direct.pressure_pa.norm_sq(),
            ratio,
        );
        assert_eq!(prepared.image_position_m(0), Some([0.0, 0.0, -1.0]));
        assert_eq!(prepared.reflected_delay_s(0), Some(4.0 / 340.0));
        assert_eq!(prepared.direct().source_delay_s(0), Some(2.0 / 340.0));
        assert_eq!(
            prepared.reflected_phase_delay_s(0),
            prepared.reflected_delay_s(0)
        );
        assert_eq!(prepared.image_position_m(1), None);
        assert_eq!(prepared.scope(), field.scope);
    }
}

#[test]
fn full_near_field_obeys_momentum_and_rigid_coordinate_changes() {
    let frequency = 217.0;
    let receiver = [0.17, -0.08, 0.05];
    let positions = [[0.02, 0.03, 0.08]];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let baseline = half_space(frequency, receiver, &positions, boundary)
            .evaluate(&[source()])
            .unwrap();
        for axis in 0..3 {
            let h = 1e-6;
            let mut left = receiver;
            let mut right = receiver;
            left[axis] -= h;
            right[axis] += h;
            let pl = half_space(frequency, left, &positions, boundary)
                .evaluate(&[source()])
                .unwrap()
                .pressure_pa;
            let pr = half_space(frequency, right, &positions, boundary)
                .evaluate(&[source()])
                .unwrap()
                .pressure_pa;
            let gradient = (pr - pl).scale(0.5 / h);
            let momentum =
                baseline.particle_velocity_m_s[axis] * C64::new(0.0, TAU * frequency * 1.25);
            assert!((gradient - momentum).abs() < 2e-7 * momentum.abs().max(1e-8));
        }
        let rotate = |v: [f64; 3]| [0.8 * v[0] + 0.6 * v[2], v[1], -0.6 * v[0] + 0.8 * v[2]];
        let offset = [3.0, -2.0, 4.0];
        let move_point = |p| {
            let v = rotate(p);
            core::array::from_fn(|i| v[i] + offset[i])
        };
        let force = source().force_on_fluid_n;
        let input = CompactSourcePhasor {
            volume_flow_m3_s: source().volume_flow_m3_s,
            force_on_fluid_n: [
                force[0].scale(0.8) + force[2].scale(0.6),
                force[1],
                force[0].scale(-0.6) + force[2].scale(0.8),
            ],
        };
        let plane = IdealAcousticPlane::new(offset, [0.6, 0.0, 0.8], boundary).unwrap();
        let transformed = observe(
            1.25,
            340.0,
            frequency,
            move_point(receiver),
            &[move_point(positions[0])],
        )
        .with_plane(plane, 2)
        .unwrap()
        .evaluate(&[input])
        .unwrap();
        close_complex(transformed.pressure_pa, baseline.pressure_pa);
        let v = baseline.particle_velocity_m_s;
        let expected = [
            v[0].scale(0.8) + v[2].scale(0.6),
            v[1],
            v[0].scale(-0.6) + v[2].scale(0.8),
        ];
        for axis in 0..3 {
            close_complex(transformed.particle_velocity_m_s[axis], expected[axis]);
        }
    }
}

#[test]
fn source_independence_never_erases_its_own_reflection_interference() {
    let frequency = 85.0;
    let receiver = [0.0, 0.0, 3.0];
    let position = [0.0, 0.0, 1.0];
    let a = CompactSourcePhasor::monopole(C64::from_re(1e-4));
    let one = half_space(frequency, receiver, &[position], IdealPlaneBoundary::Rigid);
    let field = one.evaluate(&[a]).unwrap();
    let independent = one.evaluate_independent_sources(&[a]).unwrap();
    close(
        independent.mean_square_pressure_pa2,
        0.5 * field.pressure_pa.norm_sq(),
    );
    for axis in 0..3 {
        close(
            independent.active_intensity_w_m2[axis],
            field.active_intensity_w_m2[axis],
        );
    }

    let two = half_space(
        frequency,
        receiver,
        &[position; 2],
        IdealPlaneBoundary::Rigid,
    );
    let negative = CompactSourcePhasor::monopole(-a.volume_flow_m3_s);
    assert_eq!(two.evaluate(&[a, negative]).unwrap().pressure_pa, C64::ZERO);
    close(
        two.evaluate_independent_sources(&[a, negative])
            .unwrap()
            .mean_square_pressure_pa2,
        2.0 * independent.mean_square_pressure_pa2,
    );
    let gamma = C64::new(0.3, 0.4);
    let b = CompactSourcePhasor::monopole(a.volume_flow_m3_s * gamma);
    let c = CompactSourcePhasor::monopole(a.volume_flow_m3_s.scale((1.0 - gamma.norm_sq()).sqrt()));
    let modes = [&[a, b][..], &[CompactSourcePhasor::default(), c][..]];
    let statistics = two
        .evaluate_coherence_modes(
            &modes,
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 8,
            },
        )
        .unwrap();
    close(
        statistics.mean_square_pressure_pa2,
        (2.0 + 2.0 * gamma.re) * independent.mean_square_pressure_pa2,
    );
    assert!(
        two.evaluate_coherence_modes(
            &modes,
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 7
            }
        )
        .is_err()
    );
}

#[test]
fn molecular_reflection_preserves_original_absorption_and_dispersion() {
    let state = relaxing_state(true);
    let frequency = 900.0;
    let omega = TAU * frequency;
    let input = [CompactSourcePhasor::monopole(C64::from_re(1e-4))];
    let prepared = molecular_observer(&state, frequency, [0.0, 0.0, 3.0], &[[0.0, 0.0, 1.0]])
        .with_plane(wall(IdealPlaneBoundary::Rigid), 2)
        .unwrap();
    assert_eq!(prepared.direct().molecular_state(), Some(&state));
    let k = prepared.direct().wavenumber_per_m();
    let green = |r: f64| phase(k.re * r).scale((-k.im * r).exp() / (4.0 * PI * r));
    let expected = (green(2.0) + green(4.0)) * C64::new(0.0, -state.gas().density * omega * 1e-4);
    close_complex(prepared.evaluate(&input).unwrap().pressure_pa, expected);
    close(
        prepared.reflected_phase_delay_s(0).unwrap(),
        4.0 * k.re / omega,
    );
    assert_ne!(
        prepared.reflected_phase_delay_s(0),
        prepared.reflected_delay_s(0)
    );
    let boundary = molecular_observer(&state, frequency, [0.2, 0.3, 0.0], &[[0.0, 0.0, 1.0]])
        .with_plane(wall(IdealPlaneBoundary::Rigid), 2)
        .unwrap()
        .evaluate(&[source()])
        .unwrap();
    assert_eq!(boundary.particle_velocity_m_s[2], C64::ZERO);
}

#[test]
fn reflected_lines_feed_existing_finite_window_beats() {
    let frequencies = [170.0, 183.0];
    let input = [CompactSourcePhasor::monopole(C64::new(1e-4, -2e-4))];
    let lines = frequencies.map(|frequency| {
        half_space(
            frequency,
            [0.0, 0.0, 3.0],
            &[[0.0, 0.0, 1.0]],
            IdealPlaneBoundary::Rigid,
        )
        .evaluate_harmonic_line(&input)
        .unwrap()
    });
    let prepared = PreparedHarmonicField::new(
        &lines,
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap();
    let start = 0.013;
    let duration = 0.007;
    let window = HarmonicWindow::new(start, duration).unwrap();
    let observed = prepared.observe_window(window).unwrap();
    let samples = 20000;
    let mut mean_square = 0.0;
    for i in 0..samples {
        let time = start + duration * (f64::from(i) + 0.5) / f64::from(samples);
        let pressure = frequencies
            .iter()
            .map(|frequency| {
                let omega = TAU * frequency;
                let g = phase(omega * (2.0 / 340.0 - time)).scale(1.0 / (8.0 * PI))
                    + phase(omega * (4.0 / 340.0 - time)).scale(1.0 / (16.0 * PI));
                (g * input[0].volume_flow_m3_s * C64::new(0.0, -1.25 * omega)).re
            })
            .sum::<f64>();
        mean_square += pressure * pressure / f64::from(samples);
    }
    assert!((observed.mean_square_pressure_pa2 - mean_square).abs() < 1e-7 * mean_square);
}

#[test]
fn half_space_admission_cancellation_and_dynamic_range_are_explicit() {
    assert!(IdealAcousticPlane::new([0.0; 3], [0.0; 3], IdealPlaneBoundary::Rigid).is_err());
    assert!(
        IdealAcousticPlane::new(
            [f64::NAN, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            IdealPlaneBoundary::Rigid
        )
        .is_err()
    );
    assert_eq!(
        IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1e308], IdealPlaneBoundary::Rigid)
            .unwrap()
            .inward_normal(),
        [0.0, 0.0, 1.0]
    );
    for source_z in [-1.0, 0.0] {
        assert!(
            observe(1.25, 340.0, 100.0, [1.0, 0.0, 2.0], &[[0.0, 0.0, source_z]])
                .with_plane(wall(IdealPlaneBoundary::Rigid), 2)
                .is_err()
        );
    }
    assert!(
        observe(1.25, 340.0, 100.0, [1.0, 0.0, -2.0], &[[0.0, 0.0, 1.0]])
            .with_plane(wall(IdealPlaneBoundary::Rigid), 2)
            .is_err()
    );
    let direct = || observe(1.25, 340.0, 100.0, [1.0, 0.0, 2.0], &[[0.0, 0.0, 1.0]]);
    assert!(
        direct()
            .with_plane(wall(IdealPlaneBoundary::Rigid), 1)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        direct().with_plane_cancellable(wall(IdealPlaneBoundary::Rigid), 2, || {
            polls += 1;
            polls == 4
        }),
        Err(fs_aeroac::AeroacError::Cancelled { .. })
    ));
    let prepared = half_space(
        100.0,
        [1.0, 0.0, 2.0],
        &[[0.0, 0.0, 1.0]; 3],
        IdealPlaneBoundary::Rigid,
    );
    let inputs = [1e16, 1.0, -1e16].map(|q| CompactSourcePhasor::monopole(C64::from_re(q)));
    let expected = direct()
        .with_plane(wall(IdealPlaneBoundary::Rigid), 2)
        .unwrap()
        .evaluate(&[inputs[1]])
        .unwrap();
    let accepted_pressure = prepared.evaluate(&inputs).unwrap().pressure_pa;
    close_complex(accepted_pressure, expected.pressure_pa);
    assert!(matches!(
        prepared.evaluate_cancellable(&inputs, || true),
        Err(fs_aeroac::AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        prepared.evaluate(&inputs).unwrap().pressure_pa,
        accepted_pressure
    );
    assert!(prepared.evaluate(&inputs[..2]).is_err());
    let empty = half_space(100.0, [0.0; 3], &[], IdealPlaneBoundary::Rigid);
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
}
