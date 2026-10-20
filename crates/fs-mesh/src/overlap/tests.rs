use super::*;

const TET: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

fn with_cx<R>(cancel: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    if cancel {
        gate.request();
    }
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 28,
                kernel_id: 28,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-14 * (1.0 + b.abs()), "{a} != {b}");
}
fn compare(a: VolumeMoments, b: VolumeMoments) {
    close(a.volume, b.volume);
    for i in 0..3 {
        close(a.first[i], b.first[i]);
        for j in 0..3 {
            close(a.second[i][j], b.second[i][j]);
        }
    }
}
fn right_tet() -> VolumeMoments {
    // Simplex monomial integral: ∫ x^a y^b z^c = a!b!c!/(a+b+c+3)!.
    VolumeMoments {
        volume: 1.0 / 6.0,
        first: [1.0 / 24.0; 3],
        second: core::array::from_fn(|i| {
            core::array::from_fn(|j| if i == j { 1.0 / 60.0 } else { 1.0 / 120.0 })
        }),
    }
}

#[test]
fn full_simplex_has_analytic_mass_centroid_and_inertia_moments() {
    with_cx(false, |cx| {
        let actual = tetrahedron_box_moments(TET, [0.; 3], [1.; 3], [0.; 3], cx).unwrap();
        compare(actual, right_tet());
        let mut reversed = TET;
        reversed.swap(1, 2);
        compare(
            tetrahedron_box_moments(reversed, [0.; 3], [1.; 3], [0.; 3], cx).unwrap(),
            actual,
        );
        assert_eq!(
            actual,
            tetrahedron_box_moments(TET, [0.; 3], [1.; 3], [0.; 3], cx).unwrap()
        );
    });
}

#[test]
fn clipped_corner_box_has_independent_rectangular_integrals() {
    with_cx(false, |cx| {
        // This entire cube lies inside x+y+z <= 1. All six clipping planes
        // contribute; no tetrahedron vertex lies in the interior of the box.
        let actual = tetrahedron_box_moments(TET, [0.05; 3], [0.25; 3], [0.; 3], cx).unwrap();
        let volume = 0.2_f64.powi(3);
        let mean = 0.15;
        let diagonal = (0.25_f64.powi(3) - 0.05_f64.powi(3)) / (3.0 * 0.2);
        compare(
            actual,
            VolumeMoments {
                volume,
                first: [volume * mean; 3],
                second: core::array::from_fn(|i| {
                    core::array::from_fn(|j| volume * if i == j { diagonal } else { mean * mean })
                }),
            },
        );
    });
}

#[test]
fn halfspace_cut_preserves_the_translated_small_simplex() {
    with_cx(false, |cx| {
        let actual = tetrahedron_box_moments(TET, [0.5, 0., 0.], [1.; 3], [0.; 3], cx).unwrap();
        let original = right_tet();
        let shift = [0.5, 0., 0.];
        let volume = original.volume / 8.0;
        let first = core::array::from_fn(|i| original.first[i] / 16.0 + shift[i] * volume);
        let second = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                original.second[i][j] / 32.0
                    + (shift[i] * original.first[j] + shift[j] * original.first[i]) / 16.0
                    + shift[i] * shift[j] * volume
            })
        });
        compare(
            actual,
            VolumeMoments {
                volume,
                first,
                second,
            },
        );
    });
}

#[test]
fn voxel_partition_and_world_translation_preserve_all_moments() {
    with_cx(false, |cx| {
        let mut sum = VolumeMoments::default();
        for x in 0..4 {
            for y in 0..4 {
                for z in 0..4 {
                    let lo = [x as f64 / 4., y as f64 / 4., z as f64 / 4.];
                    let hi = lo.map(|v| v + 0.25);
                    let m = tetrahedron_box_moments(TET, lo, hi, [0.; 3], cx).unwrap();
                    sum.volume += m.volume;
                    for i in 0..3 {
                        sum.first[i] += m.first[i];
                        for j in 0..3 {
                            sum.second[i][j] += m.second[i][j];
                        }
                    }
                }
            }
        }
        compare(sum, right_tet());
        let shift = [1024., -32., 16.];
        let shifted = TET.map(|p| core::array::from_fn(|i| p[i] + shift[i]));
        compare(
            tetrahedron_box_moments(shifted, shift, shift.map(|x| x + 1.), shift, cx).unwrap(),
            right_tet(),
        );
        assert_eq!(
            tetrahedron_box_moments(TET, [1.; 3], [2.; 3], [0.; 3], cx).unwrap(),
            VolumeMoments::default()
        );
    });
}

#[test]
fn invalid_geometry_and_cancellation_do_not_return_partial_overlap() {
    with_cx(false, |cx| {
        assert_eq!(
            tetrahedron_box_moments([[0.; 3]; 4], [0.; 3], [1.; 3], [0.; 3], cx),
            Err(OverlapError::DegenerateTetrahedron)
        );
        assert_eq!(
            tetrahedron_box_moments(TET, [0.; 3], [0.; 3], [0.; 3], cx),
            Err(OverlapError::InvalidInput)
        );
        assert_eq!(
            tetrahedron_box_moments(TET, [f64::NAN; 3], [1.; 3], [0.; 3], cx),
            Err(OverlapError::InvalidInput)
        );
    });
    with_cx(true, |cx| {
        assert_eq!(
            tetrahedron_box_moments(TET, [0.; 3], [1.; 3], [0.; 3], cx),
            Err(OverlapError::Cancelled)
        );
    });
}
