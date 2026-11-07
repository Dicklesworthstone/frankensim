//! Physical stress-cell coherence and bounded reflection at impedance walls.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    ImpedanceWallFieldKind, PlanarImpedanceWall, PreparedImpedanceWallQuadrupole,
    PreparedImpedanceWallQuadrupoleObserver as Observer,
};
use fs_aeroac::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const HZ: f64 = 700.0;
const Z0: f64 = 1.2 * 343.0;

fn limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 524_288,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 1e-8,
            absolute_velocity_tolerance_m_s_per_n_m: 1e-11,
            relative_tolerance: 1e-10,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        },
    }
}
fn wall(normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        normalized.map(|z| z.scale(Z0)),
    )
    .unwrap()
}
fn stress() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.0007, 0.0001),
            C64::new(-0.0003, 0.0002),
            C64::new(0.0009, -0.0004),
            C64::new(0.00023, -0.00015),
            C64::new(-0.00011, 0.0003),
            C64::new(0.0004, 0.0002),
        ],
    }
}
fn scaled(source: QuadrupolePhasor, scale: C64) -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: source.stress_moment_n_m.map(|z| z * scale),
    }
}
fn components(value: &CoherentCompactField) -> [C64; 4] {
    [
        value.pressure_pa,
        value.particle_velocity_m_s[0],
        value.particle_velocity_m_s[1],
        value.particle_velocity_m_s[2],
    ]
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 1e-20 + 3e-12 * expected.abs(),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn statistics(actual: &CompactFieldStatistics, independent_fields: &[[C64; 4]]) {
    let moment = |i: usize, j: usize| {
        independent_fields
            .iter()
            .fold(C64::ZERO, |sum, f| sum + f[i].scale(0.5) * f[j].conj())
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for i in 0..3 {
        for j in 0..3 {
            close(actual.velocity_covariance_m2_s2[i][j], moment(i + 1, j + 1));
        }
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            moment(0, i + 1),
        );
    }
}

#[test]
fn finite_wall_arrays_match_native_cells_and_preserve_robin_absorption() {
    let wall = wall(Some(C64::new(1.7, 0.4)));
    let receiver = [0.13, -0.05, 0.0];
    let positions = [[0.0, 0.0, -0.21], [-0.07, 0.03, -0.14]];
    let sources = [stress(), scaled(stress(), C64::new(-0.4, 0.7))];
    let native: Vec<_> = positions
        .iter()
        .map(|&position| {
            PreparedImpedanceWallQuadrupole::new(
                wall,
                HZ,
                position,
                receiver,
                limits().quadrupole,
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    for total in [false, true] {
        let prepared = if total {
            Observer::new_total(wall, HZ, receiver, &positions, limits(), &mut || false)
        } else {
            Observer::new_reflected(wall, HZ, receiver, &positions, limits(), &mut || false)
        }
        .unwrap();
        let fields: Vec<_> = native
            .iter()
            .zip(sources)
            .map(|(cell, source)| {
                components(
                    &if total {
                        cell.evaluate_total(source)
                    } else {
                        cell.evaluate_reflected(source)
                    }
                    .unwrap(),
                )
            })
            .collect();
        let expected: [C64; 4] = core::array::from_fn(|i| fields[0][i] + fields[1][i]);
        let value = prepared.evaluate(&sources).unwrap();
        for (actual, expected) in components(&value).into_iter().zip(expected) {
            close(actual, expected);
        }
        let line = prepared.evaluate_harmonic_line(&sources).unwrap();
        assert_eq!(line.frequency_hz, HZ);
        assert_eq!(line.pressure_pa, value.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, value.particle_velocity_m_s);
        assert_eq!(prepared.source_count(), 2);
        assert_eq!(prepared.observer_position_m(), receiver);
        assert_eq!(prepared.frequency_hz(), HZ);
        assert_eq!(prepared.wall().impedance_pa_s_m(), wall.impedance_pa_s_m());
        assert_eq!(
            prepared.field_kind(),
            if total {
                ImpedanceWallFieldKind::Total
            } else {
                ImpedanceWallFieldKind::Reflected
            }
        );
        assert_eq!(prepared.source_position_m(1), Some(positions[1]));
        assert!(prepared.source_position_m(2).is_none());
        assert!(prepared.source_diagnostics(2).is_none());
        assert_eq!(
            prepared.evaluations(),
            native
                .iter()
                .map(|v| v.diagnostics().evaluations)
                .sum::<usize>()
        );
        assert!(prepared.evaluations() > 0);
        statistics(
            &prepared.evaluate_independent(&sources, || false).unwrap(),
            &fields,
        );
        statistics(
            &prepared
                .evaluate_coherence_modes(
                    &[&sources],
                    CoherenceBudget {
                        max_modes: 1,
                        max_source_mode_pairs: 2,
                    },
                    || false,
                )
                .unwrap(),
            &[expected],
        );
        if total {
            let z = wall.impedance_pa_s_m().unwrap();
            let residual = value.pressure_pa - z * value.particle_velocity_m_s[2];
            assert!(residual.abs() < 1e-9 * value.pressure_pa.abs() + 1e-12);
            let absorbed = 0.5 * z.re * value.particle_velocity_m_s[2].abs().powi(2);
            assert!(absorbed > 0.0);
            assert!((value.active_intensity_w_m2[2] - absorbed).abs() < 1e-9 * absorbed + 1e-15);
        }
    }
}

#[test]
fn stress_cell_cancellation_and_partial_coherence_keep_original_phases() {
    let position = [0.0, 0.0, -0.2];
    let prepared = Observer::new_total(
        wall(None),
        HZ,
        [0.1, -0.05, -0.08],
        &[position; 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    let source = stress();
    let opposite = [source, scaled(source, -C64::ONE)];
    let zero = prepared.evaluate(&opposite).unwrap();
    for component in components(&zero) {
        close(component, C64::ZERO);
    }
    let one = components(
        &prepared
            .evaluate(&[source, QuadrupolePhasor::default()])
            .unwrap(),
    );
    assert!(one[0].abs() > 0.0);
    statistics(
        &prepared.evaluate_independent(&opposite, || false).unwrap(),
        &[one, one],
    );
    let budget = CoherenceBudget {
        max_modes: 2,
        max_source_mode_pairs: 4,
    };
    statistics(
        &prepared
            .evaluate_coherence_modes(&[&opposite], budget, || false)
            .unwrap(),
        &[[C64::ZERO; 4]],
    );
    let correlated = [source, scaled(source, C64::from_re(-0.5))];
    let residual = [
        QuadrupolePhasor::default(),
        scaled(source, C64::from_re(0.75_f64.sqrt())),
    ];
    statistics(
        &prepared
            .evaluate_coherence_modes(&[&correlated, &residual], budget, || false)
            .unwrap(),
        &[one],
    );
}

#[test]
fn independent_cells_retain_ideal_wall_direct_reflected_cross_terms() {
    let position = [0.0, 0.0, -0.2];
    let receiver = [0.15, 0.0, 0.0];
    let mut source = QuadrupolePhasor::default();
    source.stress_moment_n_m[2] = C64::new(0.002, 0.001);
    for normalized in [None, Some(C64::ZERO)] {
        let wall = wall(normalized);
        let mut controls = limits();
        controls.maximum_total_evaluations = 1;
        controls.quadrupole.maximum_evaluations = 1;
        let reflected =
            Observer::new_reflected(wall, HZ, receiver, &[position], controls, &mut || false)
                .unwrap();
        let total =
            Observer::new_total(wall, HZ, receiver, &[position], controls, &mut || false).unwrap();
        let returned = reflected.evaluate(&[source]).unwrap();
        let full = total.evaluate(&[source]).unwrap();
        assert_eq!(total.evaluations(), 0);
        assert!(total.source_diagnostics(0).unwrap().exact_ideal_boundary);
        let independent = total.evaluate_independent(&[source], || false).unwrap();
        statistics(&independent, &[components(&full)]);
        let path_power = 0.5 * returned.pressure_pa.abs().powi(2);
        assert!(path_power > 0.0);
        if normalized.is_none() {
            close(full.pressure_pa, returned.pressure_pa.scale(2.0));
            assert!(
                (independent.mean_square_pressure_pa2 - 4.0 * path_power).abs()
                    < 1e-12 * path_power
            );
            assert!(full.particle_velocity_m_s[2].abs() < 1e-15);
        } else {
            assert!(full.pressure_pa.abs() < 1e-13 * returned.pressure_pa.abs());
            assert!(independent.mean_square_pressure_pa2 < 1e-24 * path_power);
            assert!(full.particle_velocity_m_s[2].abs() > 0.0);
        }
    }
}

#[test]
fn aggregate_work_cancellation_self_fields_and_empty_arrays_are_bounded() {
    let position = [0.0, 0.0, -0.2];
    let receiver = [0.1, 0.02, -0.07];
    let finite = wall(Some(C64::new(1.4, 0.3)));
    let one = Observer::new_reflected(finite, HZ, receiver, &[position], limits(), &mut || false)
        .unwrap();
    let mut controls = limits();
    controls.maximum_total_evaluations = one.evaluations();
    assert!(controls.maximum_total_evaluations > 0);
    assert!(
        Observer::new_reflected(finite, HZ, receiver, &[position], controls, &mut || false).is_ok()
    );
    assert!(
        Observer::new_reflected(finite, HZ, receiver, &[position; 2], controls, &mut || {
            false
        })
        .is_err()
    );
    controls = limits();
    controls.max_sources = 1;
    assert!(
        Observer::new_total(
            wall(None),
            HZ,
            receiver,
            &[position; 2],
            controls,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        Observer::new_total(wall(None), HZ, receiver, &[position], limits(), &mut || {
            true
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let self_field =
        Observer::new_reflected(wall(None), HZ, position, &[position], limits(), &mut || {
            false
        })
        .unwrap();
    assert!(
        components(&self_field.evaluate(&[stress()]).unwrap())
            .iter()
            .all(|v| v.abs().is_finite())
    );
    assert!(
        Observer::new_total(wall(None), HZ, position, &[position], limits(), &mut || {
            false
        })
        .is_err()
    );
    assert!(matches!(
        one.evaluate_cancellable(&[stress()], || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        one.evaluate_independent(&[stress()], || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mode = [stress()];
    let budget = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 1,
    };
    assert!(matches!(
        one.evaluate_coherence_modes(&[&mode], budget, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        one.evaluate_coherence_modes(&[&mode, &mode], budget, || false)
            .is_err()
    );
    assert!(one.evaluate(&[]).is_err());
    let mut bad = stress();
    bad.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(one.evaluate(&[bad]).is_err());
    let empty =
        Observer::new_total(wall(None), HZ, receiver, &[], limits(), &mut || false).unwrap();
    assert_eq!(components(&empty.evaluate(&[]).unwrap()), [C64::ZERO; 4]);
    assert!(Observer::new_total(wall(None), 0.0, receiver, &[], limits(), &mut || false).is_err());
    assert!(
        Observer::new_total(wall(None), HZ, [0.0, 0.0, 0.1], &[], limits(), &mut || {
            false
        })
        .is_err()
    );
}
