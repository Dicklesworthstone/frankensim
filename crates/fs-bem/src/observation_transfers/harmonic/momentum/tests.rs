use super::*;

const BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn observer(
    medium: Medium,
    omega: f64,
    positions: &[[f64; 3]],
    columns: &[Vec<[C64; 4]>],
) -> PreparedHarmonicObservation {
    let sources: Vec<_> = (0..columns.len())
        .map(|s| SourceId::StructuralMode(format!("source-{s}")))
        .collect();
    let receivers: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(i, &position_m)| ReceiverSpec {
            id: format!("receiver-{i}"),
            position_m,
        })
        .collect();
    let mut samples = Vec::new();
    for column in columns {
        assert_eq!(column.len(), positions.len());
        for field in column {
            samples.push(field[0]);
            samples.extend(
                field[1..]
                    .iter()
                    .map(|&v| C64::new(0.0, omega * medium.density) * v),
            );
        }
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &sources,
        &receivers,
        ContentHash([9; 32]),
        samples,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn traveling_and_standing_waves_use_full_stress_and_keep_orientation() {
    let medium = Medium {
        density: 2.0,
        sound_speed: 3.0,
    };
    let position = [0.0, 2.0, -1.0];
    let right = [C64::from_re(6.0), C64::ONE, C64::ZERO, C64::ZERO];
    let left = [C64::from_re(6.0), -C64::ONE, C64::ZERO, C64::ZERO];
    let owner = observer(medium, 7.0, &[position], &[vec![right], vec![left]]);
    let q = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [2.0, 0.0, 0.0],
    }];
    let surface = owner.momentum_surface(&q, [0.0; 3], || false).unwrap();
    let traveling = surface
        .evaluate_coherent(&[C64::ONE, C64::ZERO], BUDGET, || false)
        .unwrap();
    close(traveling.linear_momentum_rate_n[0], 2.0, 1e-14);
    let expected_torque = cross(position, [2.0, 0.0, 0.0]);
    for (actual, expected) in traveling
        .angular_momentum_rate_nm
        .into_iter()
        .zip(expected_torque)
    {
        close(actual, expected, 2e-14);
    }
    // At the pressure antinode particle velocity vanishes, but radiation
    // stress does not. A power/c or rho*u*u-only observer would return zero.
    let standing = surface
        .evaluate_coherent(&[C64::ONE; 2], BUDGET, || false)
        .unwrap();
    close(standing.linear_momentum_rate_n[0], 4.0, 2e-14);
    let transverse = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [0.0, 2.0, 0.0],
    }];
    let surface_y = owner
        .momentum_surface(&transverse, [0.0; 3], || false)
        .unwrap();
    let antinode = surface_y
        .evaluate_coherent(&[C64::ONE; 2], BUDGET, || false)
        .unwrap();
    let node = surface_y
        .evaluate_coherent(&[C64::ONE, -C64::ONE], BUDGET, || false)
        .unwrap();
    close(antinode.linear_momentum_rate_n[1], 4.0, 2e-14);
    close(node.linear_momentum_rate_n[1], -4.0, 2e-14);
    let reversed = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [-2.0, 0.0, 0.0],
    }];
    let reverse = owner
        .momentum_surface(&reversed, [0.0; 3], || false)
        .unwrap()
        .evaluate_coherent(&[C64::ONE; 2], BUDGET, || false)
        .unwrap();
    assert_eq!(
        reverse.linear_momentum_rate_n[0],
        -standing.linear_momentum_rate_n[0]
    );
}

#[test]
fn coherence_modes_match_independent_time_integrated_stress_and_torque() {
    let medium = Medium {
        density: 1.7,
        sound_speed: 4.2,
    };
    let positions = [[-0.3, 0.8, 1.1], [0.7, -0.2, 0.4]];
    let columns = [
        vec![
            [
                C64::new(2.0, -0.5),
                C64::new(0.3, 0.1),
                C64::new(-0.1, 0.4),
                C64::new(0.2, -0.3),
            ],
            [
                C64::new(-1.0, 2.0),
                C64::new(0.2, 0.4),
                C64::new(0.5, -0.2),
                C64::new(-0.1, -0.3),
            ],
        ],
        vec![
            [
                C64::new(0.5, 1.0),
                C64::new(-0.4, 0.2),
                C64::new(0.2, 0.3),
                C64::new(0.1, 0.5),
            ],
            [
                C64::new(2.0, 0.3),
                C64::new(0.1, -0.3),
                C64::new(-0.2, 0.1),
                C64::new(0.4, -0.1),
            ],
        ],
    ];
    let owner = observer(medium, 8.0, &positions, &columns);
    let q = [
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.4, -0.2, 0.7],
        },
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [-0.3, 0.8, 0.2],
        },
    ];
    let origin = [-0.1, 0.2, 0.5];
    let surface = owner.momentum_surface(&q, origin, || false).unwrap();
    let factors = [
        C64::new(0.3, 0.7),
        C64::new(-0.1, 0.6),
        C64::new(0.2, -0.5),
        C64::new(0.8, 0.1),
    ];
    let actual = surface
        .evaluate_coherence_modes(&factors, 2, BUDGET, || false)
        .unwrap();
    let mut force = [0.0; 3];
    let mut torque = [0.0; 3];
    let count = 2048;
    // An independent real-time oracle: evaluate the instantaneous quadratic
    // inviscid stress, then integrate a complete cycle for each independent bank.
    for factor in factors.chunks_exact(2) {
        for (i, sample) in q.iter().enumerate() {
            let field: [C64; 4] = core::array::from_fn(|c| {
                columns[0][i][c] * factor[0] + columns[1][i][c] * factor[1]
            });
            for t in 0..count {
                let (sin, cos) =
                    det::sin_cos(core::f64::consts::TAU * (t as f64 + 0.5) / count as f64);
                let f = field.map(|z| z.re * cos + z.im * sin);
                let v2 = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
                let isotropic = 0.5
                    * (f[0] * f[0] / (medium.density * medium.sound_speed * medium.sound_speed)
                        - medium.density * v2);
                let velocity_area: f64 =
                    (0..3).map(|j| f[j + 1] * sample.oriented_area_m2[j]).sum();
                let traction: [f64; 3] = core::array::from_fn(|j| {
                    (medium.density * f[j + 1] * velocity_area
                        + isotropic * sample.oriented_area_m2[j])
                        / count as f64
                });
                let arm = core::array::from_fn(|j| positions[i][j] - origin[j]);
                let moment = cross(arm, traction);
                for j in 0..3 {
                    force[j] += traction[j];
                    torque[j] += moment[j];
                }
            }
        }
    }
    for i in 0..3 {
        close(actual.linear_momentum_rate_n[i], force[i], 4e-14);
        close(actual.angular_momentum_rate_nm[i], torque[i], 4e-14);
    }
    let independent = surface
        .evaluate_independent_sources(&[0.25, 4.0], BUDGET, || false)
        .unwrap();
    let diagonal = surface
        .evaluate_coherence_modes(
            &[C64::from_re(0.5), C64::ZERO, C64::ZERO, C64::from_re(2.0)],
            2,
            BUDGET,
            || false,
        )
        .unwrap();
    for i in 0..3 {
        close(
            independent.linear_momentum_rate_n[i],
            diagonal.linear_momentum_rate_n[i],
            1e-14,
        );
        close(
            independent.angular_momentum_rate_nm[i],
            diagonal.angular_momentum_rate_nm[i],
            1e-14,
        );
    }
    let phase = C64::new(det::cos(0.83), det::sin(0.83));
    let shifted: Vec<_> = factors.iter().map(|&z| z * phase).collect();
    let shifted = surface
        .evaluate_coherence_modes(&shifted, 2, BUDGET, || false)
        .unwrap();
    for i in 0..3 {
        close(
            shifted.linear_momentum_rate_n[i],
            actual.linear_momentum_rate_n[i],
            2e-14,
        );
        close(
            shifted.angular_momentum_rate_nm[i],
            actual.angular_momentum_rate_nm[i],
            2e-14,
        );
    }
}

// Exact outgoing monopole + three Cartesian dipoles. p=G[a+(1+i/kr)b.n],
// G=exp(ikr)/r; differentiate analytically, independently of the stress owner.
fn multipole(
    medium: Medium,
    omega: f64,
    radius: f64,
    n: [f64; 3],
    a: C64,
    b: [C64; 3],
) -> [C64; 4] {
    let k = omega / medium.sound_speed;
    let (sin, cos) = det::sin_cos(k * radius);
    let g = C64::new(cos, sin).scale(1.0 / radius);
    let bn: C64 = (0..3).fold(C64::ZERO, |v, i| v + b[i].scale(n[i]));
    let radial = C64::new(1.0, 1.0 / (k * radius));
    let shape = a + radial * bn;
    let euler = C64::new(0.0, -1.0 / (omega * medium.density));
    let mut out = [C64::ZERO; 4];
    out[0] = g * shape;
    for i in 0..3 {
        let gradient = C64::new(-1.0 / radius, k) * shape.scale(n[i])
            + radial * (b[i] - bn.scale(n[i])).scale(1.0 / radius)
            + C64::new(0.0, -1.0 / (k * radius * radius)) * bn.scale(n[i]);
        out[i + 1] = euler * g * gradient;
    }
    out
}

fn sphere(
    radius: f64,
    center: [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<HarmonicSurfaceQuadrature>) {
    let nodes = [
        (-0.8611363115940526, 0.3478548451374538),
        (-0.3399810435848563, 0.6521451548625461),
        (0.3399810435848563, 0.6521451548625461),
        (0.8611363115940526, 0.3478548451374538),
    ];
    let mut positions = Vec::new();
    let mut directions = Vec::new();
    let mut quadrature = Vec::new();
    for (z, weight) in nodes {
        for longitude in 0..12 {
            let (sin, cos) = det::sin_cos(core::f64::consts::TAU * (longitude as f64 + 0.3) / 12.0);
            let radial = det::sqrt(1.0 - z * z);
            let n = [radial * cos, radial * sin, z];
            positions.push(core::array::from_fn(|i| center[i] + radius * n[i]));
            directions.push(n);
            quadrature.push(HarmonicSurfaceQuadrature {
                receiver: positions.len() - 1,
                oriented_area_m2: n
                    .map(|v| v * radius * radius * weight * core::f64::consts::TAU / 12.0),
            });
        }
    }
    (positions, directions, quadrature)
}

#[test]
fn near_and_far_spheres_retain_radiation_recoil_and_vortex_angular_momentum() {
    let medium = Medium {
        density: 1.3,
        sound_speed: 2.7,
    };
    let omega = 5.4;
    let a = C64::new(0.7, -0.2);
    let b = [C64::new(1.0, 0.0), C64::new(0.0, 1.0), C64::new(0.2, -0.3)];
    let force: [f64; 3] = b.map(|v| {
        (a.conj() * v).re * 4.0 * core::f64::consts::PI
            / (3.0 * medium.density * medium.sound_speed * medium.sound_speed)
    });
    let torque: [f64; 3] = core::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        (b[j].conj() * b[k] - b[k].conj() * b[j]).im * 2.0 * core::f64::consts::PI
            / (3.0 * medium.density * medium.sound_speed * omega)
    });
    let center = [0.3, -0.5, 0.8];
    let origin = [-0.2, 0.4, 0.1];
    let extrinsic = cross(core::array::from_fn(|i| center[i] - origin[i]), force);
    for radius in [0.17, 0.8, 12.0] {
        let (positions, directions, q) = sphere(radius, center);
        let fields = directions
            .iter()
            .map(|&n| multipole(medium, omega, radius, n, a, b))
            .collect();
        let owner = observer(medium, omega, &positions, &[fields]);
        let surface = owner.momentum_surface(&q, origin, || false).unwrap();
        let actual = surface
            .evaluate_coherent(&[C64::ONE], BUDGET, || false)
            .unwrap();
        for i in 0..3 {
            close(actual.linear_momentum_rate_n[i], force[i], 5e-12);
            close(
                actual.angular_momentum_rate_nm[i],
                torque[i] + extrinsic[i],
                5e-12,
            );
        }
        // A circular dipole has m=+1: angular momentum rate = radiated P/omega.
        let vortex = [C64::ONE, C64::new(0.0, 1.0), C64::ZERO];
        let fields = directions
            .iter()
            .map(|&n| multipole(medium, omega, radius, n, C64::ZERO, vortex))
            .collect();
        let owner = observer(medium, omega, &positions, &[fields]);
        let actual = owner
            .momentum_surface(&q, center, || false)
            .unwrap()
            .evaluate_coherent(&[C64::ONE], BUDGET, || false)
            .unwrap();
        let power = 4.0 * core::f64::consts::PI / (3.0 * medium.density * medium.sound_speed);
        close(actual.angular_momentum_rate_nm[2], power / omega, 5e-12);
        for f in actual.linear_momentum_rate_n {
            close(f, 0.0, 5e-12);
        }
    }
}

#[test]
fn quiet_source_and_patch_cancellation_retain_physical_residuals() {
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let wave = |amplitude| {
        [
            C64::from_re(amplitude),
            C64::from_re(amplitude),
            C64::ZERO,
            C64::ZERO,
        ]
    };
    let owner = observer(
        medium,
        1.0,
        &[[0.0; 3]],
        &[vec![wave(1e16)], vec![wave(1.0)], vec![wave(-1e16)]],
    );
    let q = [
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [2e16, 0.0, 0.0],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [2.0, 0.0, 0.0],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [-2e16, 0.0, 0.0],
        },
    ];
    let surface = owner.momentum_surface(&q, [0.0; 3], || false).unwrap();
    let result = surface
        .evaluate_coherent(&[C64::ONE; 3], BUDGET, || false)
        .unwrap();
    assert_eq!(result.linear_momentum_rate_n, [1.0, 0.0, 0.0]);
    let cancel = surface
        .evaluate_coherent(&[C64::ONE, C64::ZERO, C64::ONE], BUDGET, || false)
        .unwrap();
    assert_eq!(cancel, HarmonicMomentumFlux::default());
}

#[test]
fn admission_budgets_empty_surfaces_and_cancellation_are_atomic() {
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let fields = vec![[C64::ONE, C64::ONE, C64::ZERO, C64::ZERO]];
    let owner = observer(medium, 1.0, &[[0.0; 3]], &[fields]);
    let q = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    let surface = owner.momentum_surface(&q, [0.0; 3], || false).unwrap();
    assert_eq!(surface.sources(), owner.sources());
    assert_eq!(surface.physical_receivers(), owner.receivers());
    assert_eq!(surface.producer_hash(), owner.producer_hash());
    assert_eq!(surface.angular_frequency_rad_s(), 1.0);
    assert_eq!(surface.quadrature(), &q);
    assert_eq!(surface.moment_origin_m(), [0.0; 3]);
    let exact = HarmonicWorkBudget {
        maximum_products: 11,
    };
    assert!(
        surface
            .evaluate_coherent(&[C64::ONE], exact, || false)
            .is_ok()
    );
    assert!(
        surface
            .evaluate_coherent(
                &[C64::ONE],
                HarmonicWorkBudget {
                    maximum_products: 10
                },
                || false
            )
            .is_err()
    );
    assert!(
        surface
            .evaluate_independent_sources(
                &[1.0],
                HarmonicWorkBudget {
                    maximum_products: 7
                },
                || false
            )
            .is_ok()
    );
    assert!(
        surface
            .evaluate_independent_sources(&[-1.0], BUDGET, || false)
            .is_err()
    );
    assert!(surface.evaluate_coherent(&[], BUDGET, || false).is_err());
    assert!(
        surface
            .evaluate_coherent(&[C64::new(f64::NAN, 0.0)], BUDGET, || false)
            .is_err()
    );
    assert!(
        owner
            .momentum_surface(&q, [f64::NAN, 0.0, 0.0], || false)
            .is_err()
    );
    for bad in [
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [0.0; 3],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [f64::INFINITY, 0.0, 0.0],
        },
    ] {
        assert!(owner.momentum_surface(&[bad], [0.0; 3], || false).is_err());
    }
    for stop in 0..=3 {
        let mut calls = 0;
        let result = surface.evaluate_coherent(&[C64::ONE], BUDGET, || {
            calls += 1;
            calls > stop
        });
        assert!(matches!(result, Err(TransferError::Cancelled)));
    }
    let zero = owner
        .momentum_surface(&[], [0.0; 3], || false)
        .unwrap()
        .evaluate_coherent(&[C64::ONE], BUDGET, || false)
        .unwrap();
    assert_eq!(zero, HarmonicMomentumFlux::default());
    assert_eq!(
        surface
            .evaluate_coherence_modes(&[], 0, BUDGET, || false)
            .unwrap(),
        zero
    );
    assert_eq!(
        surface
            .evaluate_independent_sources(&[0.0], BUDGET, || false)
            .unwrap(),
        zero
    );
    assert_eq!(
        surface
            .evaluate_coherent(&[C64::ONE], BUDGET, || false)
            .unwrap()
            .linear_momentum_rate_n,
        [0.5, 0.0, 0.0]
    );
}

#[test]
fn relaxing_medium_readout_calibration_cannot_authorize_radiation_stress() {
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let owner = observer(
        medium,
        1.0,
        &[[0.0; 3]],
        &[vec![[C64::ONE, C64::ONE, C64::ZERO, C64::ZERO]]],
    )
    .without_radiation_momentum();
    // Admission is physical, even when the requested observation is empty.
    assert!(owner.momentum_surface(&[], [0.0; 3], || false).is_err());
    assert_eq!(owner.angular_frequency_rad_s(), 1.0);
    assert_eq!(owner.field(&[C64::ONE], 0).unwrap()[0], C64::ONE);
}
