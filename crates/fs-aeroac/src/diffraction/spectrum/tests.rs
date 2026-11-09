//! G1/G3 physical-frame, coherent-field and angular-quadrature regressions.
use super::*;

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:?}, expected {expected:?}, error {:e}, tolerance {tolerance:e}",
        (actual - expected).abs(),
    );
}

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}

fn frequency() -> f64 {
    7.0 * medium().sound_speed_m_s() / core::f64::consts::TAU
}

fn canonical() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}

fn prepare(
    frame: HalfPlaneFrame,
    reference: [f64; 3],
    boundary: IdealPlaneBoundary,
    waves: &[HalfPlaneIncidentWave],
) -> PreparedHalfPlaneSpectrum {
    PreparedHalfPlaneSpectrum::new(
        medium(),
        frequency(),
        frame,
        boundary,
        reference,
        waves,
        HalfPlaneSpectrumBudget {
            max_components: 4096,
        },
        || false,
    )
    .unwrap()
}

fn wave(direction: [f64; 3], amplitude: C64) -> HalfPlaneIncidentWave {
    HalfPlaneIncidentWave {
        propagation_direction: direction,
        pressure_pa: amplitude,
    }
}

fn world(frame: HalfPlaneFrame, local: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|j| {
        frame.origin_m[j] + (0..3).map(|i| local[i] * frame.axes[i][j]).sum::<f64>()
    })
}

#[test]
fn placed_wave_preserves_absolute_phase_and_rotates_actual_velocity() {
    let origin = [1.25, -2.0, 0.5];
    let reference = [-0.125, 0.5, 0.25];
    let frame = HalfPlaneFrame::new(origin, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let local_direction = [-0.48, -0.64, 0.6];
    let direction = [0.6, -0.48, -0.64];
    let amplitude = C64::new(0.4, -0.7);
    let local_point = [0.25, -0.375, 0.5];
    let point = world(frame, local_point);
    let phase: f64 = 7.0
        * (0..3)
            .map(|i| direction[i] * (origin[i] - reference[i]))
            .sum::<f64>();
    let at_edge = amplitude * C64::new(phase.cos(), phase.sin());
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let field = prepare(frame, reference, boundary, &[wave(direction, amplitude)])
            .evaluate(HalfPlaneReceiver::WorldPosition(point), || false)
            .unwrap();
        let expected = HalfPlaneDiffraction::new(medium(), frequency(), local_direction, boundary)
            .unwrap()
            .evaluate(local_point, at_edge, || false)
            .unwrap();
        close(field.pressure_pa, expected.pressure_pa, 8e-15);
        for (w, l) in [(0, 2), (1, 0), (2, 1)] {
            close(
                field.particle_velocity_m_s[w],
                expected.particle_velocity_m_s[l],
                3e-17,
            );
            assert!(
                (field.active_intensity_w_m2[w] - expected.active_intensity_w_m2[l]).abs() < 3e-17
            );
            assert!(
                (field.reactive_intensity_w_m2[w] - expected.reactive_intensity_w_m2[l]).abs()
                    < 3e-17
            );
        }
        // Translate the entire world, including the actual phase reference.
        // Binary coordinates make this a bitwise metrology check, not a test
        // confounded by rounding a tiny displacement off a huge coordinate.
        let shift = [8.0, 16.0, -32.0];
        let moved_origin = core::array::from_fn(|i| origin[i] + shift[i]);
        let moved_reference = core::array::from_fn(|i| reference[i] + shift[i]);
        let moved_point = core::array::from_fn(|i| point[i] + shift[i]);
        let moved_frame =
            HalfPlaneFrame::new(moved_origin, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
        let moved = prepare(
            moved_frame,
            moved_reference,
            boundary,
            &[wave(direction, amplitude)],
        )
        .evaluate(HalfPlaneReceiver::WorldPosition(moved_point), || false)
        .unwrap();
        assert_eq!(moved.pressure_pa, field.pressure_pa);
        assert_eq!(moved.particle_velocity_m_s, field.particle_velocity_m_s);
    }
}

#[test]
fn source_reference_changes_retain_the_same_physical_incident_field() {
    let direction = [-0.48, -0.64, 0.6];
    let amplitude = C64::new(0.2, 0.7);
    let reference = [0.25, -0.5, 0.125];
    let phase: f64 = 7.0
        * (direction[0] * reference[0] + direction[1] * reference[1] + direction[2] * reference[2]);
    let shifted_amplitude = amplitude * C64::new(phase.cos(), phase.sin());
    let first = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &[wave(direction, amplitude)],
    );
    let second = prepare(
        canonical(),
        reference,
        IdealPlaneBoundary::Rigid,
        &[wave(direction, shifted_amplitude)],
    );
    for point in [[0.25, 0.3, 0.1], [0.5, -0.25, -0.125], [-0.3, 0.1, 0.2]] {
        let a = first
            .evaluate(HalfPlaneReceiver::WorldPosition(point), || false)
            .unwrap();
        let b = second
            .evaluate(HalfPlaneReceiver::WorldPosition(point), || false)
            .unwrap();
        close(a.pressure_pa, b.pressure_pa, 3e-15);
        for axis in 0..3 {
            close(
                a.particle_velocity_m_s[axis],
                b.particle_velocity_m_s[axis],
                2e-17,
            );
        }
    }
}

#[test]
fn coherent_fields_are_summed_before_power_and_retain_quiet_components() {
    let direction = [-0.48, -0.64, 0.6];
    let point = HalfPlaneReceiver::WorldPosition([0.2, -0.3, 0.1]);
    let single = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &[wave(direction, C64::ONE)],
    )
    .evaluate(point, || false)
    .unwrap();
    let doubled = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &[wave(direction, C64::ONE); 2],
    )
    .evaluate(point, || false)
    .unwrap();
    close(doubled.pressure_pa, single.pressure_pa.scale(2.0), 0.0);
    for axis in 0..3 {
        close(
            doubled.particle_velocity_m_s[axis],
            single.particle_velocity_m_s[axis].scale(2.0),
            0.0,
        );
        assert_eq!(
            doubled.active_intensity_w_m2[axis],
            4.0 * single.active_intensity_w_m2[axis]
        );
        assert_eq!(
            doubled.reactive_intensity_w_m2[axis],
            4.0 * single.reactive_intensity_w_m2[axis]
        );
    }
    let cancelled = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &[wave(direction, C64::ONE), wave(direction, -C64::ONE)],
    )
    .evaluate(point, || false)
    .unwrap();
    assert_eq!(cancelled.pressure_pa, C64::ZERO);
    assert_eq!(cancelled.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(cancelled.active_intensity_w_m2, [0.0; 3]);

    // The weak middle component survives cancellation of two strong fields.
    // A plain f64 running sum would round it away before the last component.
    let weak = C64::from_re(1e-20);
    let coherent = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &[
            wave(direction, C64::ONE),
            wave(direction, weak),
            wave(direction, -C64::ONE),
        ],
    )
    .evaluate(point, || false)
    .unwrap();
    close(coherent.pressure_pa, single.pressure_pa * weak, 1e-35);
    for axis in 0..3 {
        close(
            coherent.particle_velocity_m_s[axis],
            single.particle_velocity_m_s[axis] * weak,
            1e-37,
        );
    }
}

#[test]
fn explicit_rotated_sheet_faces_preserve_the_boundary_conditions() {
    let frame = HalfPlaneFrame::new([0.3, -0.4, 0.2], [2.0, -1.0, 0.0], [1.0, 2.0, 2.0]).unwrap();
    let local_direction = [-0.48, -0.64, 0.6];
    let direction =
        core::array::from_fn(|i| (0..3).map(|j| frame.axes[j][i] * local_direction[j]).sum());
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let spectrum = prepare(
            frame,
            frame.origin_m,
            boundary,
            &[wave(direction, C64::ONE)],
        );
        let canonical =
            HalfPlaneDiffraction::new(medium(), frequency(), local_direction, boundary).unwrap();
        for (face, signed_zero) in [
            (HalfPlaneFace::PositiveNormal, 0.0),
            (HalfPlaneFace::NegativeNormal, -0.0),
        ] {
            let field = spectrum
                .evaluate(
                    HalfPlaneReceiver::SheetFace {
                        distance_from_edge_m: 0.25,
                        along_edge_m: 0.125,
                        face,
                    },
                    || false,
                )
                .unwrap();
            let expected = canonical
                .evaluate([0.25, signed_zero, 0.125], C64::ONE, || false)
                .unwrap();
            close(field.pressure_pa, expected.pressure_pa, 2e-15);
            let normal_velocity = (0..3).fold(C64::ZERO, |sum, i| {
                sum + field.particle_velocity_m_s[i].scale(frame.axes[1][i])
            });
            match boundary {
                IdealPlaneBoundary::Rigid => close(normal_velocity, C64::ZERO, 2e-18),
                IdealPlaneBoundary::PressureRelease => assert_eq!(field.pressure_pa, C64::ZERO),
            }
        }
    }
}

#[test]
fn world_velocity_matches_the_gradient_of_the_coherent_pressure() {
    let frame = HalfPlaneFrame::new([0.3, -0.4, 0.2], [2.0, -1.0, 0.0], [1.0, 2.0, 2.0]).unwrap();
    let spectrum = prepare(
        frame,
        [-0.1, 0.2, 0.3],
        IdealPlaneBoundary::Rigid,
        &[
            wave([0.6, -0.8, 0.0], C64::new(0.4, -0.2)),
            wave([-0.48, -0.64, 0.6], C64::new(-0.1, 0.7)),
        ],
    );
    let point = world(frame, [0.3, -0.2, 0.1]);
    let line = spectrum
        .evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(point), || false)
        .unwrap();
    let field = spectrum
        .evaluate(HalfPlaneReceiver::WorldPosition(point), || false)
        .unwrap();
    assert_eq!(line.frequency_hz, frequency());
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    let h = 1e-5;
    for axis in 0..3 {
        let pressure = |offset: f64| {
            let mut receiver = point;
            receiver[axis] += offset;
            spectrum
                .evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(receiver), || false)
                .unwrap()
                .pressure_pa
        };
        let gradient = (pressure(-2.0 * h) - pressure(-h).scale(8.0) + pressure(h).scale(8.0)
            - pressure(2.0 * h))
        .scale(1.0 / (12.0 * h));
        let expected = C64::new(gradient.im, -gradient.re)
            .scale(1.0 / (core::f64::consts::TAU * frequency() * medium().density_kg_m3()));
        close(field.particle_velocity_m_s[axis], expected, 4e-12);
    }
}

#[test]
fn smooth_angular_spectrum_has_a_measured_refinement_ladder() {
    let beam = |count: usize| {
        let step = 1.1 / count as f64;
        let waves: Vec<_> = (0..count)
            .map(|i| {
                let angle = -0.55 + (i as f64 + 0.5) * step;
                wave(
                    [angle.sin(), -angle.cos(), 0.0],
                    C64::from_re(step * (-0.5 * (angle / 0.2).powi(2)).exp()),
                )
            })
            .collect();
        prepare(
            canonical(),
            [-0.1, 0.3, 0.0],
            IdealPlaneBoundary::Rigid,
            &waves,
        )
    };
    let reference = beam(1024);
    let points = [[0.3, -0.2, 0.1], [-0.4, -0.3, 0.1], [0.1, 0.2, -0.1]];
    let mut previous = f64::INFINITY;
    for count in [16, 32, 64] {
        let spectrum = beam(count);
        let mut error = 0.0_f64;
        for point in points {
            let r = reference
                .evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(point), || false)
                .unwrap();
            let a = spectrum
                .evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(point), || false)
                .unwrap();
            error = error.max((a.pressure_pa - r.pressure_pa).abs());
            for axis in 0..3 {
                error = error.max(
                    (a.particle_velocity_m_s[axis] - r.particle_velocity_m_s[axis]).abs()
                        * medium().characteristic_impedance_pa_s_m(),
                );
            }
        }
        assert!(
            error < previous * 0.3,
            "count {count}, error {error:e}, previous {previous:e}"
        );
        previous = error;
    }
    assert!(
        previous < 3e-5,
        "finest angular-quadrature error {previous:e}"
    );
}

#[test]
fn geometry_spectrum_budget_faces_and_cancellation_refuse_explicitly() {
    assert!(HalfPlaneFrame::new([f64::NAN, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).is_err());
    assert!(HalfPlaneFrame::new([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]).is_err());
    assert!(HalfPlaneFrame::new([0.0; 3], [0.0, 0.1, 1.0], [0.0, 1.0, 0.0]).is_err());
    let rotated = HalfPlaneFrame::new([0.0; 3], [2.0, -1.0, 0.0], [1.0, 2.0, 2.0]).unwrap();
    for direction in [rotated.axes[0], rotated.axes[2]] {
        assert!(
            PreparedHalfPlaneSpectrum::new(
                medium(),
                frequency(),
                rotated,
                IdealPlaneBoundary::Rigid,
                [0.0; 3],
                &[wave(direction, C64::ONE)],
                HalfPlaneSpectrumBudget { max_components: 1 },
                || false,
            )
            .is_err()
        );
    }
    let input = wave([0.0, -1.0, 0.0], C64::ONE);
    let spectrum = prepare(canonical(), [0.0; 3], IdealPlaneBoundary::Rigid, &[input]);
    for point in [
        [0.25, 0.0, 0.0],
        [0.25, -0.0, 0.0],
        [0.0; 3],
        [f64::NAN, 0.1, 0.0],
    ] {
        assert!(
            spectrum
                .evaluate(HalfPlaneReceiver::WorldPosition(point), || false)
                .is_err()
        );
    }
    assert!(
        spectrum
            .evaluate(HalfPlaneReceiver::WorldPosition([-0.25, 0.0, 0.0]), || {
                false
            })
            .is_ok()
    );
    assert!(
        spectrum
            .evaluate(
                HalfPlaneReceiver::SheetFace {
                    distance_from_edge_m: 0.0,
                    along_edge_m: 0.0,
                    face: HalfPlaneFace::PositiveNormal
                },
                || false
            )
            .is_err()
    );
    for (input, budget) in [
        (vec![], 1),
        (vec![input; 2], 1),
        (vec![input; MAX_HALF_PLANE_COMPONENTS + 1], usize::MAX),
    ] {
        assert!(
            PreparedHalfPlaneSpectrum::new(
                medium(),
                frequency(),
                canonical(),
                IdealPlaneBoundary::Rigid,
                [0.0; 3],
                &input,
                HalfPlaneSpectrumBudget {
                    max_components: budget
                },
                || false
            )
            .is_err()
        );
    }
    let maximum = prepare(
        canonical(),
        [0.0; 3],
        IdealPlaneBoundary::Rigid,
        &vec![input; MAX_HALF_PLANE_COMPONENTS],
    );
    assert_eq!(maximum.component_count(), MAX_HALF_PLANE_COMPONENTS);
    for after in [0, 1, 3] {
        let mut polls = 0;
        let result = maximum.evaluate(HalfPlaneReceiver::WorldPosition([0.3, -0.2, 0.0]), || {
            polls += 1;
            polls > after
        });
        assert!(matches!(result, Err(AeroacError::Cancelled { .. })));
    }
    let mut polls = 0;
    let result = PreparedHalfPlaneSpectrum::new(
        medium(),
        frequency(),
        canonical(),
        IdealPlaneBoundary::Rigid,
        [0.0; 3],
        &[input; 4],
        HalfPlaneSpectrumBudget { max_components: 4 },
        || {
            polls += 1;
            polls > 2
        },
    );
    assert!(matches!(result, Err(AeroacError::Cancelled { .. })));
}
