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

fn fan(centre: Point, face: [Point; 3]) -> ([Point; 3], f64) {
    let edges = face.map(|p| sub(p, centre));
    (edges, determinant(edges[0], edges[1], edges[2]))
}

#[test]
fn the_two_fans_that_refused_the_cut_tetrahedron_fixtures_are_zero_volume_by_the_error_model() {
    // Printed from overlap.rs (the `det < 0.0` refusal of `integrate`) on ts2, 2026-11-11, while the
    // free-surface transfer fixture a_cut_tetrahedron_transfers_true_volume_and_affine_motion refused
    // with Unrepresentable: the clipped face has three vertices that agree to 1e-16. G0 (algebraic
    // law: the model bounds the evaluated determinant).
    let (edges, det) = fan(
        [0.3823529411764705, 0.1411764705882353, 0.1588235294117647],
        [
            [0.5, 0.23529411764705868, 0.26470588235294124],
            [0.5, 0.2352941176470588, 0.2647058823529411],
            [0.4999999999999999, 0.2352941176470588, 0.26470588235294124],
        ],
    );
    assert_eq!(det, -1.2245106889247525e-19);
    // Inside the bound by a factor of at least ten (68 measured).
    assert!(det < 0.0 && -det * 10.0 < fan_det_noise(edges));
    // The same, for the copy of the fixture with every coordinate a multiple of 2^-12.
    let (edges, det) = fan(
        [0.5999999999999999, 0.1285714285714286, 0.042857142857142844],
        [
            [0.7142857142857143, 0.21428571428571422, 0.07142857142857142],
            [0.7142857142857143, 0.21428571428571427, 0.07142857142857137],
            [0.7142857142857143, 0.21428571428571427, 0.07142857142857142],
        ],
    );
    assert_eq!(det, -1.2390881971261612e-20);
    // Inside the bound by a factor of at least ten (332 measured).
    assert!(det < 0.0 && -det * 10.0 < fan_det_noise(edges));
}

#[test]
fn a_negative_determinant_beyond_the_error_model_and_an_inside_out_polyhedron_are_still_refused() {
    // A planted fan e = (0.1 x, 0.1 y, -z): det = -0.01 z, bound about 96 eps x 0.01 = 2.1e-16.
    // z = 1e-14 gives det -1e-16 (inside, zero volume); z = 1e-13 gives -1e-15 (beyond, refused).
    for (z, refused) in [(1e-14, false), (1e-13, true)] {
        let edges = [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, -z]];
        let det = determinant(edges[0], edges[1], edges[2]);
        assert!(det < 0.0);
        assert_eq!(det < -fan_det_noise(edges), refused, "z = {z}");
    }
    // The faces of TET with both windings: the outward one integrates, the inside-out one is
    // refused by the fan.
    let faces: Vec<Face> = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
        .map(|face| face.map(|i| TET[i]).to_vec())
        .into();
    let inside_out: Vec<Face> = faces
        .iter()
        .map(|face| face.iter().rev().copied().collect())
        .collect();
    with_cx(false, |cx| {
        let a = integrate(&faces, [0.; 3], 1.0, [0.; 3], cx);
        let b = integrate(&inside_out, [0.; 3], 1.0, [0.; 3], cx);
        let (outward, refused) = if a.is_ok() { (a, b) } else { (b, a) };
        compare(outward.unwrap(), right_tet());
        assert_eq!(refused, Err(OverlapError::Unrepresentable));
    });
}
