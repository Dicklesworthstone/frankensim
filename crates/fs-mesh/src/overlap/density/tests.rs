use super::super::tetrahedron_box_moments;
use super::*;

const TET: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const RHO: [f64; 4] = [1.0, 2.0, 3.0, 4.0];

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
                seed: 43,
                kernel_id: 43,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-13 * a.abs().max(b.abs()).max(1e-300),
        "{a:e} != {b:e}"
    );
}

fn compare(a: DensityMoments, b: DensityMoments) {
    close(a.volume, b.volume);
    close(a.mass, b.mass);
    for j in 0..3 {
        close(a.first[j], b.first[j]);
        for k in 0..3 {
            close(a.second[j][k], b.second[j][k]);
        }
    }
}

fn whole_tet() -> DensityMoments {
    // rho = 1 + x + 2y + 3z on the unit simplex. Independent rational
    // integrals use integral(x^a y^b z^c) = a! b! c! / (a+b+c+3)!.
    DensityMoments {
        volume: 1.0 / 6.0,
        mass: 5.0 / 12.0,
        first: [1.0 / 10.0, 13.0 / 120.0, 7.0 / 60.0],
        second: [
            [7.0 / 180.0, 1.0 / 48.0, 1.0 / 45.0],
            [1.0 / 48.0, 2.0 / 45.0, 17.0 / 720.0],
            [1.0 / 45.0, 17.0 / 720.0, 1.0 / 20.0],
        ],
    }
}

fn rectangular_monomial(lower: [f64; 3], upper: [f64; 3], powers: [i32; 3]) -> f64 {
    (0..3)
        .map(|j| {
            // Factored difference of powers keeps the independent thin-slab oracle
            // from subtracting nearly equal endpoint powers.
            let n = powers[j];
            let sum: f64 = (0..=n)
                .map(|k| upper[j].powi(n - k) * lower[j].powi(k))
                .sum();
            (upper[j] - lower[j]) * sum / (n + 1) as f64
        })
        .product()
}

fn rectangle(lower: [f64; 3], upper: [f64; 3], rho: [f64; 4]) -> DensityMoments {
    let integrate = |powers: [i32; 3]| {
        let mut value = rho[0] * rectangular_monomial(lower, upper, powers);
        for j in 0..3 {
            let mut linear = powers;
            linear[j] += 1;
            value += (rho[j + 1] - rho[0]) * rectangular_monomial(lower, upper, linear);
        }
        value
    };
    DensityMoments {
        volume: rectangular_monomial(lower, upper, [0; 3]),
        mass: integrate([0; 3]),
        first: core::array::from_fn(|j| {
            let mut powers = [0; 3];
            powers[j] = 1;
            integrate(powers)
        }),
        second: core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                let mut powers = [0; 3];
                powers[j] += 1;
                powers[k] += 1;
                integrate(powers)
            })
        }),
    }
}

#[test]
fn g1_original_density_changes_centroid_inertia_and_affine_kinetic_energy() {
    with_cx(false, |cx| {
        let m =
            tetrahedron_box_density_moments(TET, RHO, [0.0; 3], [1.0; 3], [0.0; 3], cx).unwrap();
        compare(m, whole_tet());
        assert!((m.first[0] / m.mass - 0.25).abs() > 0.009);
        let mut reversed = TET;
        let mut reversed_rho = RHO;
        reversed.swap(1, 2);
        reversed_rho.swap(1, 2);
        compare(
            tetrahedron_box_density_moments(
                reversed,
                reversed_rho,
                [0.0; 3],
                [1.0; 3],
                [0.0; 3],
                cx,
            )
            .unwrap(),
            m,
        );
        // An affine velocity [0, -z, y] has K = (second_yy+second_zz)/2.
        // Equal total mass at uniform density gives a different rotational energy.
        close(0.5 * (m.second[1][1] + m.second[2][2]), 17.0 / 360.0);
        assert!((0.5 * (m.second[1][1] + m.second[2][2]) - m.mass / 10.0).abs() > 0.005);
    });
}

#[test]
fn g1_six_plane_cuts_retain_density_in_thin_and_ordinary_rectangles() {
    with_cx(false, |cx| {
        for (lower, upper, rho) in [
            ([0.05; 3], [0.25; 3], RHO),
            (
                [0.2, 0.1, 0.05],
                [0.2 + 2.0_f64.powi(-20), 0.2, 0.1],
                [1e-12, 1e12, 3.0, 4.0],
            ),
        ] {
            let actual =
                tetrahedron_box_density_moments(TET, rho, lower, upper, [0.0; 3], cx).unwrap();
            compare(actual, rectangle(lower, upper, rho));
        }
    });
}

#[test]
fn g3_partition_and_world_translation_preserve_weighted_storage() {
    with_cx(false, |cx| {
        let mut sum = DensityMoments::default();
        for x in 0..4 {
            for y in 0..4 {
                for z in 0..4 {
                    let lower = [x as f64 / 4.0, y as f64 / 4.0, z as f64 / 4.0];
                    let upper = lower.map(|x| x + 0.25);
                    let m = tetrahedron_box_density_moments(TET, RHO, lower, upper, [0.0; 3], cx)
                        .unwrap();
                    sum.volume += m.volume;
                    sum.mass += m.mass;
                    for j in 0..3 {
                        sum.first[j] += m.first[j];
                        for k in 0..3 {
                            sum.second[j][k] += m.second[j][k];
                        }
                    }
                }
            }
        }
        compare(sum, whole_tet());
        let shift = [1024.0, -32.0, 16.0];
        let translated = TET.map(|p| core::array::from_fn(|j| p[j] + shift[j]));
        compare(
            tetrahedron_box_density_moments(
                translated,
                RHO,
                shift,
                shift.map(|x| x + 1.0),
                shift,
                cx,
            )
            .unwrap(),
            whole_tet(),
        );
        assert_eq!(
            tetrahedron_box_density_moments(TET, RHO, [1.0; 3], [2.0; 3], [0.0; 3], cx).unwrap(),
            DensityMoments::default()
        );
    });
}

#[test]
fn g3_uniform_field_and_physical_unit_scaling_match_original_geometry() {
    with_cx(false, |cx| {
        let geometry = tetrahedron_box_moments(TET, [0.05; 3], [0.25; 3], [0.0; 3], cx).unwrap();
        let rho = 713.0;
        let m = tetrahedron_box_density_moments(TET, [rho; 4], [0.05; 3], [0.25; 3], [0.0; 3], cx)
            .unwrap();
        compare(
            m,
            DensityMoments {
                volume: geometry.volume,
                mass: rho * geometry.volume,
                first: geometry.first.map(|x| rho * x),
                second: geometry.second.map(|row| row.map(|x| rho * x)),
            },
        );
        let scaled = tetrahedron_box_density_moments(
            TET.map(|p| p.map(|x| 3.0 * x)),
            RHO.map(|r| 7.0 * r),
            [0.0; 3],
            [3.0; 3],
            [0.0; 3],
            cx,
        )
        .unwrap();
        let original = whole_tet();
        compare(
            scaled,
            DensityMoments {
                volume: 27.0 * original.volume,
                mass: 7.0 * 27.0 * original.mass,
                first: original.first.map(|x| 7.0 * 81.0 * x),
                second: original.second.map(|row| row.map(|x| 7.0 * 243.0 * x)),
            },
        );
    });
}

#[test]
fn g4_invalid_density_geometry_and_cancellation_return_no_partial_moments() {
    with_cx(false, |cx| {
        for rho in [[0.0; 4], [-1.0; 4], [f64::NAN; 4], [f64::INFINITY; 4]] {
            assert_eq!(
                tetrahedron_box_density_moments(TET, rho, [0.0; 3], [1.0; 3], [0.0; 3], cx),
                Err(OverlapError::InvalidInput)
            );
        }
        assert_eq!(
            tetrahedron_box_density_moments(
                TET,
                [f64::from_bits(1), f64::MAX, 1.0, 1.0],
                [0.0; 3],
                [1.0; 3],
                [0.0; 3],
                cx,
            ),
            Err(OverlapError::Unrepresentable)
        );
        assert_eq!(
            tetrahedron_box_density_moments([[0.0; 3]; 4], RHO, [0.0; 3], [1.0; 3], [0.0; 3], cx),
            Err(OverlapError::DegenerateTetrahedron)
        );
    });
    with_cx(true, |cx| {
        assert_eq!(
            tetrahedron_box_density_moments(TET, RHO, [0.0; 3], [1.0; 3], [0.0; 3], cx),
            Err(OverlapError::Cancelled)
        );
    });
}
