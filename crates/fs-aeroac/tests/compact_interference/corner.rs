//! G1/G3 finite image groups, physical wall conditions and source coherence.
use super::*;
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};

fn walls(boundaries: &[IdealPlaneBoundary]) -> Vec<IdealAcousticPlane> {
    boundaries
        .iter()
        .enumerate()
        .map(|(axis, &boundary)| {
            let mut normal = [0.0; 3];
            normal[axis] = 1.0;
            IdealAcousticPlane::new([0.0; 3], normal, boundary).unwrap()
        })
        .collect()
}

fn drive() -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: C64::new(2e-4, -1e-4),
        force_on_fluid_n: [C64::new(0.2, -0.3), C64::new(-0.1, 0.4), C64::new(0.6, 0.2)],
    }
}

// Independent Cartesian image positions and physical source vectors. This
// constructs no transfer matrices and uses the scalar retarded-pressure oracle.
fn images(
    position: [f64; 3],
    source: CompactSourcePhasor,
    boundaries: &[IdealPlaneBoundary],
) -> (Vec<[f64; 3]>, Vec<CompactSourcePhasor>) {
    let mut positions = Vec::new();
    let mut sources = Vec::new();
    for mask in 0..(1 << boundaries.len()) {
        let mut p = position;
        let mut s = source;
        for (axis, boundary) in boundaries.iter().enumerate() {
            if mask & (1 << axis) != 0 {
                p[axis] = -p[axis];
                s.force_on_fluid_n[axis] = -s.force_on_fluid_n[axis];
                if *boundary == IdealPlaneBoundary::PressureRelease {
                    s.volume_flow_m3_s = -s.volume_flow_m3_s;
                    s.force_on_fluid_n = s.force_on_fluid_n.map(|v| -v);
                }
            }
        }
        positions.push(p);
        sources.push(s);
    }
    (positions, sources)
}

#[test]
fn complete_image_groups_match_independent_pressure_and_near_field_momentum() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let source = [0.02, 0.03, 0.08];
    let receiver = [0.17, 0.08, 0.05];
    let boundaries = [Rigid, PressureRelease, Rigid];
    for count in 1..=3 {
        for frequency in [17.0, 1300.0] {
            let prepared = observe(1.25, 340.0, frequency, receiver, &[source])
                .with_orthogonal_planes(&walls(&boundaries[..count]), 1 << count)
                .unwrap();
            let field = prepared.evaluate(&[drive()]).unwrap();
            let (positions, inputs) = images(source, drive(), &boundaries[..count]);
            let pressure = |x| reference_pressure(1.25, 340.0, frequency, x, &positions, &inputs);
            close_complex(field.pressure_pa, pressure(receiver));
            for axis in 0..3 {
                let h = 2e-6;
                let mut left = receiver;
                let mut right = receiver;
                left[axis] -= h;
                right[axis] += h;
                let expected = (pressure(right) - pressure(left)).scale(0.5 / h)
                    / C64::new(0.0, TAU * frequency * 1.25);
                let actual = field.particle_velocity_m_s[axis];
                assert!((actual - expected).abs() < 3e-7 * actual.abs().max(1e-8));
            }
        }
    }
}

#[test]
fn all_wall_conditions_hold_for_monopoles_dipoles_edges_and_corners() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let position = [0.3, 0.5, 0.7];
    for boundaries in [
        [Rigid; 3],
        [PressureRelease; 3],
        [Rigid, PressureRelease, Rigid],
    ] {
        for count in [2, 3] {
            let planes = walls(&boundaries[..count]);
            for axis in 0..count {
                let mut receiver = [0.8, 0.9, 1.1];
                receiver[axis] = 0.0;
                for component in 0..4 {
                    let mut input = CompactSourcePhasor::default();
                    if component == 0 {
                        input.volume_flow_m3_s = C64::new(1e-4, 2e-4);
                    } else {
                        input.force_on_fluid_n[component - 1] = C64::new(0.4, -0.7);
                    }
                    let field = observe(1.25, 340.0, 317.0, receiver, &[position])
                        .with_orthogonal_planes(&planes, 1 << count)
                        .unwrap()
                        .evaluate(&[input])
                        .unwrap();
                    match boundaries[axis] {
                        Rigid => assert_eq!(field.particle_velocity_m_s[axis], C64::ZERO),
                        PressureRelease => assert_eq!(field.pressure_pa, C64::ZERO),
                    }
                    assert_eq!(field.active_intensity_w_m2[axis], 0.0);
                }
            }
            let receiver = if count == 2 {
                [0.0, 0.0, 1.0]
            } else {
                [0.0; 3]
            };
            let direct = observe(1.25, 340.0, 317.0, receiver, &[position]);
            let original = direct.evaluate(&[drive()]).unwrap();
            let field = direct
                .with_orthogonal_planes(&planes, 1 << count)
                .unwrap()
                .evaluate(&[drive()])
                .unwrap();
            if boundaries[..count].iter().all(|b| *b == Rigid) {
                close_complex(
                    field.pressure_pa,
                    original.pressure_pa.scale((1 << count) as f64),
                );
            } else {
                assert_eq!(field.pressure_pa, C64::ZERO);
            }
        }
    }
}

#[test]
fn reflection_group_is_invariant_under_rigid_motion_and_plane_order() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let position = [0.3, 0.5, 0.7];
    let receiver = [0.8, 0.9, 1.1];
    let planes = walls(&[Rigid, PressureRelease, Rigid]);
    let field = observe(1.25, 340.0, 317.0, receiver, &[position])
        .with_orthogonal_planes(&planes, 8)
        .unwrap()
        .evaluate(&[drive()])
        .unwrap();
    let rotate = |v: [f64; 3]| [0.8 * v[0] + 0.6 * v[2], v[1], -0.6 * v[0] + 0.8 * v[2]];
    let offset = [3.0, -2.0, 4.0];
    let move_point = |p| {
        let v = rotate(p);
        core::array::from_fn(|i| v[i] + offset[i])
    };
    let f = drive().force_on_fluid_n;
    let input = CompactSourcePhasor {
        volume_flow_m3_s: drive().volume_flow_m3_s,
        force_on_fluid_n: [
            f[0].scale(0.8) + f[2].scale(0.6),
            f[1],
            f[0].scale(-0.6) + f[2].scale(0.8),
        ],
    };
    let transformed: Vec<_> = planes
        .iter()
        .rev()
        .map(|p| IdealAcousticPlane::new(offset, rotate(p.inward_normal()), p.boundary()).unwrap())
        .collect();
    let actual = observe(
        1.25,
        340.0,
        317.0,
        move_point(receiver),
        &[move_point(position)],
    )
    .with_orthogonal_planes(&transformed, 8)
    .unwrap()
    .evaluate(&[input])
    .unwrap();
    close_complex(actual.pressure_pa, field.pressure_pa);
    let v = field.particle_velocity_m_s;
    let expected = [
        v[0].scale(0.8) + v[2].scale(0.6),
        v[1],
        v[0].scale(-0.6) + v[2].scale(0.8),
    ];
    for axis in 0..3 {
        close_complex(actual.particle_velocity_m_s[axis], expected[axis]);
    }
}

#[test]
fn folded_images_keep_physical_source_coherence_and_the_original_clock() {
    let planes = walls(&[IdealPlaneBoundary::Rigid; 3]);
    let position = [0.3, 0.5, 0.7];
    let original = observe(1.25, 340.0, 317.0, [0.0; 3], &[position; 2]);
    let direct = original
        .evaluate(&[drive(), CompactSourcePhasor::default()])
        .unwrap();
    let prepared = original
        .clone()
        .with_orthogonal_planes(&planes, 16)
        .unwrap();
    assert_eq!(prepared.reflection_planes(), planes);
    assert_eq!(prepared.source_position_m(0), Some(position));
    assert_eq!(prepared.source_delay_s(0), original.source_delay_s(0));
    assert_eq!(prepared.acoustic(), original.acoustic());
    let negative = CompactSourcePhasor {
        volume_flow_m3_s: -drive().volume_flow_m3_s,
        force_on_fluid_n: drive().force_on_fluid_n.map(|v| -v),
    };
    let input = [drive(), negative];
    assert_eq!(prepared.evaluate(&input).unwrap().pressure_pa, C64::ZERO);
    let independent = prepared.evaluate_independent_sources(&input).unwrap();
    close(
        independent.mean_square_pressure_pa2,
        64.0 * direct.pressure_pa.norm_sq(),
    );
    let one = [drive(), CompactSourcePhasor::default()];
    let field = prepared.evaluate(&one).unwrap();
    let statistics = prepared
        .evaluate_coherence_modes(
            &[&one],
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: 2,
            },
        )
        .unwrap();
    close(
        statistics.mean_square_pressure_pa2,
        0.5 * field.pressure_pa.norm_sq(),
    );
    assert_eq!(statistics.scope, field.scope);
    let line = prepared.evaluate_harmonic_line(&one).unwrap();
    assert_eq!(line.frequency_hz, 317.0);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert!(
        prepared
            .evaluate_coherence_modes(
                &[&one],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1,
                }
            )
            .is_err()
    );
}

#[test]
fn molecular_corner_retains_each_paths_original_absorption_and_phase() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let state = relaxing_state(true);
    let boundaries = [Rigid, PressureRelease, Rigid];
    let receiver = [1.8, 0.9, 2.1];
    let position = [0.3, 0.5, 0.7];
    let frequency = 900.0;
    let input = CompactSourcePhasor::monopole(C64::new(1e-4, 2e-4));
    let prepared = molecular_observer(&state, frequency, receiver, &[position])
        .with_orthogonal_planes(&walls(&boundaries), 8)
        .unwrap();
    assert_eq!(prepared.molecular_state(), Some(&state));
    let k = prepared.wavenumber_per_m();
    let (positions, sources) = images(position, input, &boundaries);
    let mut expected = C64::ZERO;
    for (p, source) in positions.iter().zip(sources) {
        let r = (0..3)
            .map(|i| (receiver[i] - p[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        let green = phase(k.re * r).scale((-k.im * r).exp() / (4.0 * PI * r));
        expected = expected
            + green
                * source.volume_flow_m3_s
                * C64::new(0.0, -state.gas().density * TAU * frequency);
    }
    close_complex(prepared.evaluate(&[input]).unwrap().pressure_pa, expected);
}

#[test]
fn boundary_preparation_refuses_invalid_groups_budgets_and_cancellation() {
    use IdealPlaneBoundary::Rigid;
    let planes = walls(&[Rigid; 3]);
    let original = || observe(1.25, 340.0, 317.0, [0.8, 0.9, 1.1], &[[0.3, 0.5, 0.7]]);
    assert!(original().with_orthogonal_planes(&[], 8).is_err());
    assert!(
        original()
            .with_orthogonal_planes(&[planes[0]; 4], 16)
            .is_err()
    );
    assert!(
        original()
            .with_orthogonal_planes(&[planes[0]; 2], 4)
            .is_err()
    );
    assert!(original().with_orthogonal_planes(&planes, 7).is_err());
    let oblique = IdealAcousticPlane::new([0.0; 3], [1.0, 1.0, 0.0], Rigid).unwrap();
    assert!(
        original()
            .with_orthogonal_planes(&[planes[0], oblique], 4)
            .is_err()
    );
    assert!(
        observe(1.25, 340.0, 317.0, [-0.1, 0.9, 1.1], &[[0.3, 0.5, 0.7]])
            .with_orthogonal_planes(&planes, 8)
            .is_err()
    );
    assert!(
        observe(1.25, 340.0, 317.0, [0.8, 0.9, 1.1], &[[0.3, 0.0, 0.7]])
            .with_orthogonal_planes(&planes, 8)
            .is_err()
    );
    for stop in [1, 8, 25] {
        let mut polls = 0;
        assert!(matches!(
            original().with_orthogonal_planes_cancellable(&planes, 8, || {
                polls += 1;
                polls == stop
            }),
            Err(fs_aeroac::AeroacError::Cancelled { .. })
        ));
    }
    let prepared = original().with_orthogonal_planes(&planes, 8).unwrap();
    assert!(prepared.clone().with_orthogonal_planes(&planes, 8).is_err());
    assert!(prepared.with_plane(planes[0], 2).is_err());
    let single = original().with_orthogonal_planes(&planes[..1], 2).unwrap();
    let established = original().with_plane(planes[0], 2).unwrap();
    close_complex(
        single.evaluate(&[drive()]).unwrap().pressure_pa,
        established.evaluate(&[drive()]).unwrap().pressure_pa,
    );
}
