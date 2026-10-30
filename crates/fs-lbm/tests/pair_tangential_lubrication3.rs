//! Source scalar coefficients, objective pair torques, and native moving fluid.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3, LubricatedPlanesExchange3,
    LubricationPlane3, MovingPartialSaturationBody3, PartialSaturationBody3,
    PartialSaturationError3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const MU: f64 = 0.11;
const CUT: f64 = 0.004;
const ROUGH: f64 = 1e-7;
const NORMAL_CUT: f64 = 0.04;
const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; error={:.3e}",
        a - b
    );
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn grid(walls: bool) -> BoundaryGrid3 {
    BoundaryGrid3::new(
        8,
        8,
        8,
        0.83,
        [0.0; 3],
        if walls {
            BoundarySpec3::new([FaceBoundary3::Wall { velocity: [0.0; 3] }; 6])
        } else {
            BoundarySpec3::periodic()
        },
    )
}
fn body(center: [f64; 3], mass: f64) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        inertia_world: [[0.6, 0.04, -0.02], [0.04, 0.8, 0.03], [-0.02, 0.03, 0.9]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    }
}
fn model(mu: f64, maximum_pairs: usize) -> SphereLubrication3 {
    SphereLubrication3::new(mu, NORMAL_CUT, ROUGH, maximum_pairs)
        .unwrap()
        .with_tangential_pairs(CUT, ROUGH)
        .unwrap()
        .with_resistance_mode_limit(32)
        .unwrap()
}
fn frozen<const N: usize>(
    mut bodies: [PartialSaturationBody3; N],
    radii: [f64; N],
    external: [[f64; 6]; N],
    model: SphereLubrication3,
    planes: &[LubricationPlane3],
) -> ([PartialSaturationBody3; N], LubricatedPlanesExchange3<N>) {
    // Zero occupied volume isolates the real resistance/native inertia operator.
    // The sampled moving-body test below exercises the coupled fluid mechanism.
    let fields = core::array::from_fn::<_, N, _>(|_| vec![0.0; 512]);
    let report = grid(false)
        .step_with_lubricated_partial_saturation_planes(
            fields.each_ref().map(Vec::as_slice),
            radii,
            &mut bodies,
            external,
            model,
            planes,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    (bodies, report)
}
fn speeds<const N: usize>(bodies: &[PartialSaturationBody3; N]) -> [[f64; 6]; N] {
    bodies.map(|b| {
        [
            b.velocity[0],
            b.velocity[1],
            b.velocity[2],
            b.angular_velocity[0],
            b.angular_velocity[1],
            b.angular_velocity[2],
        ]
    })
}
fn normal_drag(mu: f64, radius: f64, gap: f64) -> f64 {
    6.0 * core::f64::consts::PI
        * mu
        * radius
        * radius
        * (1.0 / gap.max(ROUGH) - 1.0 / NORMAL_CUT).max(0.0)
}

// Independent direct YA/YB/YC scalar normalization from Townsend (2023),
// https://arxiv.org/abs/1802.08226, and the author's G2ya/G2yb/G2yc/G4yc.
// No production modal vectors, eigenvalues or Schur response are used here.
// Ordering: translational drag, the two translation/rotation cross terms,
// rotational self coefficients, and reciprocal rotational cross coefficient.
fn scalar_coefficients(a: f64, b: f64, gap: f64, mu: f64) -> [f64; 6] {
    if gap.max(ROUGH) >= CUT {
        return [0.0; 6];
    }
    let l = b / a;
    let log = (CUT / gap.max(ROUGH)).ln();
    let p = core::f64::consts::PI * mu * log;
    let ya = (6.0 * p * a) * (4.0 * l * (2.0 + l + 2.0 * l * l) / (15.0 * (1.0 + l).powi(3)));
    let ba = -(4.0 * p * a * a) * (-l * (4.0 + l) / (5.0 * (1.0 + l).powi(2)));
    let inv = 1.0 / l;
    let bb = -(4.0 * p * b * b) * (-inv * (4.0 + inv) / (5.0 * (1.0 + inv).powi(2)));
    let ca = (8.0 * p * a.powi(3)) * (2.0 * l / (5.0 * (1.0 + l)));
    let cb = (8.0 * p * b.powi(3)) * (2.0 * inv / (5.0 * (1.0 + inv)));
    let cross = (8.0 * p * a.powi(3)) * (l * l / (10.0 * (1.0 + l)));
    [ya, ba, bb, ca, cb, cross]
}

fn physical_loads<const N: usize>(
    bodies: &[PartialSaturationBody3; N],
    radii: [f64; N],
    velocity: &[[f64; 6]; N],
    mu: f64,
    planes: &[LubricationPlane3],
) -> [[f64; 6]; N] {
    let mut result = [[0.0; 6]; N];
    for a in 0..N {
        for b in a + 1..N {
            let delta = sub(bodies[b].center, bodies[a].center);
            let distance = delta[0].hypot(delta[1]).hypot(delta[2]);
            let n = delta.map(|x| x / distance);
            let gap = (distance - radii[a]) - radii[b];
            let relative = core::array::from_fn(|i| velocity[a][i] - velocity[b][i]);
            let vn = dot(n, relative);
            let tangent = core::array::from_fn(|i| relative[i] - vn * n[i]);
            let wa = [velocity[a][3], velocity[a][4], velocity[a][5]];
            let wb = [velocity[b][3], velocity[b][4], velocity[b][5]];
            let wan = dot(wa, n);
            let wbn = dot(wb, n);
            let cat = cross(wa, n);
            let cbt = cross(wb, n);
            let nt = cross(n, tangent);
            let [aa, ba, bb, ca, cb, cc] = scalar_coefficients(radii[a], radii[b], gap, mu);
            let completion = (radii[a] + radii[b]) / distance;
            let normal = normal_drag(mu, radii[a] * radii[b] / (radii[a] + radii[b]), gap);
            for i in 0..3 {
                let f = -normal * vn * n[i]
                    - aa * completion * completion * tangent[i]
                    - ba * completion * cat[i]
                    - bb * completion * cbt[i];
                result[a][i] += f;
                result[b][i] -= f;
                result[a][i + 3] -=
                    ba * completion * nt[i] + ca * (wa[i] - wan * n[i]) + cc * (wb[i] - wbn * n[i]);
                result[b][i + 3] -=
                    bb * completion * nt[i] + cc * (wa[i] - wan * n[i]) + cb * (wb[i] - wbn * n[i]);
            }
        }
    }
    for plane in planes {
        let n = plane.normal();
        for a in 0..N {
            let gap = dot(n, sub(bodies[a].center, plane.point())) - radii[a];
            let relative = core::array::from_fn(|i| velocity[a][i] - plane.velocity()[i]);
            let vn = dot(relative, n);
            let omega = [velocity[a][3], velocity[a][4], velocity[a][5]];
            let wn = dot(omega, n);
            let spin = cross(omega, n);
            let translation = cross(n, relative);
            let zeta = normal_drag(mu, radii[a], gap);
            let k = if gap.max(ROUGH) < CUT {
                4.0 * core::f64::consts::PI * mu * radii[a] / 5.0 * (CUT / gap.max(ROUGH)).ln()
            } else {
                0.0
            };
            for i in 0..3 {
                result[a][i] += -zeta * vn * n[i] - 4.0 * k * (relative[i] - vn * n[i])
                    + k * radii[a] * spin[i];
                result[a][i + 3] += k * radii[a] * translation[i]
                    - 4.0 * k * radii[a] * radii[a] * (omega[i] - wn * n[i]);
            }
        }
    }
    result
}

fn direct_step<const N: usize>(
    bodies: &[PartialSaturationBody3; N],
    radii: [f64; N],
    external: [[f64; 6]; N],
    planes: &[LubricationPlane3],
) -> Vec<f64> {
    let dim = 6 * N;
    let before = speeds(bodies);
    let affine = physical_loads(bodies, radii, &[[0.0; 6]; N], MU, planes);
    let initial_load = physical_loads(bodies, radii, &before, MU, planes);
    let mut m = vec![vec![0.0; dim + 1]; dim];
    for a in 0..N {
        for i in 0..6 {
            let row = 6 * a + i;
            if i < 3 {
                m[row][row] = bodies[a].mass;
            } else {
                for j in 0..3 {
                    m[row][6 * a + 3 + j] = bodies[a].inertia_world[i - 3][j];
                }
            }
            m[row][dim] = (0..6)
                .map(|j| m[row][6 * a + j] * before[a][j])
                .sum::<f64>()
                + external[a][i]
                + 0.5 * (initial_load[a][i] + affine[a][i]);
        }
    }
    for col in 0..dim {
        let mut unit = [[0.0; 6]; N];
        unit[col / 6][col % 6] = 1.0;
        let load = physical_loads(bodies, radii, &unit, MU, planes);
        for a in 0..N {
            for i in 0..6 {
                m[6 * a + i][col] -= 0.5 * (load[a][i] - affine[a][i]);
            }
        }
    }
    // Independent physical velocity/spin Gaussian elimination with pivoting.
    for col in 0..dim {
        let pivot = (col..dim)
            .max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))
            .unwrap();
        m.swap(col, pivot);
        let diagonal = m[col][col];
        assert!(diagonal.abs() > 1e-12);
        for j in col..=dim {
            m[col][j] /= diagonal;
        }
        for row in 0..dim {
            if row != col {
                let coefficient = m[row][col];
                for j in col..=dim {
                    m[row][j] -= coefficient * m[col][j];
                }
            }
        }
    }
    m.into_iter().map(|row| row[dim]).collect()
}

#[test]
fn unequal_scalar_prefactors_and_cross_signs_match_and_approach_contact_limits() {
    for radii in [[0.5, 0.5], [0.2, 2.0], [2.0, 0.2], [0.4, 1.2], [1.2, 0.4]] {
        let mut convergence = Vec::new();
        for gap in [0.002, 0.0002, 0.00002] {
            let mut initial = [
                body([3.0, 3.0, 3.0], 1.3),
                body([3.0, 3.0, 3.0 + radii[0] + radii[1] + gap], 2.1),
            ];
            let actual_gap = ((initial[1].center[2] - initial[0].center[2]) - radii[0]) - radii[1];
            let coefficients = scalar_coefficients(radii[0], radii[1], actual_gap, MU);
            let mut worst = 0.0_f64;
            for column in 0..3 {
                for b in &mut initial {
                    b.velocity = [0.0; 3];
                    b.angular_velocity = [0.0; 3];
                }
                let speed = 0.001;
                match column {
                    0 => initial[0].velocity[0] = speed,
                    1 => initial[0].angular_velocity[1] = speed,
                    _ => initial[1].angular_velocity[1] = speed,
                }
                let v = speeds(&initial);
                let load = physical_loads(&initial, radii, &v, MU, &[]);
                // Exact physical external counter-load holds the independent
                // coordinate probe stationary despite the real inertia solve.
                let external = load.map(|w| w.map(|x| -x));
                let (after, result) = frozen(initial, radii, external, model(MU, 1), &[]);
                let film = result.coupling.lubrication;
                assert_eq!(film.active_pairs, 1);
                assert_eq!(film.active_tangential_pairs, 1);
                assert_eq!(film.active_resistance_modes, 5);
                let final_v = speeds(&after);
                for a in 0..2 {
                    for i in 0..6 {
                        close(final_v[a][i], v[a][i], 2e-16);
                        let observed = if i < 3 {
                            film.body_impulses[a][i]
                        } else {
                            film.body_spin_impulses[a][i - 3]
                        };
                        close(observed, load[a][i], 4e-17);
                    }
                }
                let observed = [
                    -film.body_impulses[0][0] / speed,
                    -film.body_spin_impulses[0][1] / speed,
                    -film.body_spin_impulses[1][1] / speed,
                ];
                let contact = match column {
                    0 => [coefficients[0], coefficients[1], coefficients[2]],
                    1 => [coefficients[1], coefficients[3], coefficients[5]],
                    _ => [coefficients[2], coefficients[5], coefficients[4]],
                };
                for i in 0..3 {
                    worst = worst.max((observed[i] / contact[i] - 1.0).abs());
                }
                let work = load
                    .into_iter()
                    .zip(v)
                    .map(|(f, v)| f.into_iter().zip(v).map(|(f, v)| f * v).sum::<f64>())
                    .sum::<f64>();
                close(film.dissipated_energy, -work, 3e-20);
                close(film.work_residual, 0.0, 3e-20);
                for i in 0..3 {
                    close(
                        film.body_impulses[0][i] + film.body_impulses[1][i],
                        0.0,
                        0.0,
                    );
                    close(
                        film.body_angular_impulses[0][i] + film.body_angular_impulses[1][i],
                        0.0,
                        2e-17,
                    );
                }
            }
            convergence.push(worst);
        }
        assert!(
            convergence.windows(2).all(|v| v[0] / v[1] > 9.8),
            "{radii:?}: {convergence:?}"
        );
        eprintln!("YA/YB/YC contact-limit relative errors {radii:?}: {convergence:?}");
    }
}

#[test]
fn counterrotation_with_zero_surface_slip_dissipates_and_refines_at_second_order() {
    let radius = 0.5;
    let mut initial = [body([3.0, 4.0, 4.0], 2.0), body([4.002, 4.0, 4.0], 2.0)];
    let gap = (initial[1].center[0] - initial[0].center[0]) - 1.0;
    let k = 4.0 * core::f64::consts::PI * MU * (radius / 2.0) / 5.0 * (CUT / gap).ln();
    for (a, b) in initial.iter_mut().enumerate() {
        b.inertia_world = [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]];
        b.angular_velocity = [0.006, 0.0, if a == 0 { 0.01 } else { -0.01 }];
    }
    close(
        radius * (initial[0].angular_velocity[2] + initial[1].angular_velocity[2]),
        0.0,
        0.0,
    );
    let rate = 3.0 * k * radius * radius / 0.2;
    let (after, result) = frozen(initial, [radius; 2], [[0.0; 6]; 2], model(MU, 1), &[]);
    let gain = (1.0 - 0.5 * rate) / (1.0 + 0.5 * rate);
    for a in 0..2 {
        assert_eq!(after[a].velocity, [0.0; 3]);
        close(
            after[a].angular_velocity[2],
            gain * initial[a].angular_velocity[2],
            2e-17,
        );
        close(
            after[a].angular_velocity[0],
            initial[a].angular_velocity[0],
            0.0,
        );
    }
    assert!(result.coupling.lubrication.dissipated_energy > 0.0);
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let dt = 1.0 / f64::from(count);
        let mut state = initial.map(|b| PartialSaturationBody3 {
            angular_velocity: b.angular_velocity.map(|x| x * dt),
            ..b
        });
        for _ in 0..count {
            state = frozen(state, [radius; 2], [[0.0; 6]; 2], model(MU * dt, 1), &[]).0;
        }
        errors.push((state[0].angular_velocity[2] / dt - 0.01 * (-rate).exp()).abs());
    }
    assert!(
        errors
            .windows(2)
            .all(|e| (3.9..4.1).contains(&(e[0] / e[1]))),
        "{errors:?}"
    );
    eprintln!("pair counter-rotation refinement errors: {errors:?}");
}

#[test]
fn rigid_motion_nullspace_scene_covariance_and_endpoint_exchange_are_physical() {
    let radii = [0.4, 0.7];
    let mut initial = [body([3.0, 3.0, 3.0], 1.3), body([3.0, 3.0, 4.102], 2.1)];
    initial[0].velocity = [0.009, -0.003, 0.001];
    initial[1].velocity = [-0.002, 0.004, -0.001];
    initial[0].angular_velocity = [0.005, 0.006, -0.002];
    initial[1].angular_velocity = [-0.004, 0.003, 0.001];
    let (after, report) = frozen(initial, radii, [[0.0; 6]; 2], model(MU, 1), &[]);
    let omega = [0.001, -0.002, 0.003];
    let boost = [0.004, -0.003, 0.002];
    let shifted = initial.map(|b| {
        let orbital = cross(omega, b.center);
        PartialSaturationBody3 {
            velocity: core::array::from_fn(|i| b.velocity[i] + boost[i] + orbital[i]),
            angular_velocity: core::array::from_fn(|i| b.angular_velocity[i] + omega[i]),
            ..b
        }
    });
    let (changed, changed_report) = frozen(shifted, radii, [[0.0; 6]; 2], model(MU, 1), &[]);
    for a in 0..2 {
        for i in 0..3 {
            close(
                changed[a].velocity[i] - shifted[a].velocity[i],
                after[a].velocity[i] - initial[a].velocity[i],
                3e-17,
            );
            close(
                changed[a].angular_velocity[i] - shifted[a].angular_velocity[i],
                after[a].angular_velocity[i] - initial[a].angular_velocity[i],
                3e-17,
            );
            close(
                changed_report.coupling.lubrication.body_spin_impulses[a][i],
                report.coupling.lubrication.body_spin_impulses[a][i],
                8e-17,
            );
        }
    }
    let common = initial.map(|b| PartialSaturationBody3 {
        velocity: core::array::from_fn(|i| boost[i] + cross(omega, b.center)[i]),
        angular_velocity: omega,
        ..b
    });
    let (same, null) = frozen(common, radii, [[0.0; 6]; 2], model(MU, 1), &[]);
    for a in 0..2 {
        for i in 0..3 {
            close(same[a].velocity[i], common[a].velocity[i], 3e-18);
            close(same[a].angular_velocity[i], omega[i], 3e-18);
            close(null.coupling.lubrication.body_impulses[a][i], 0.0, 8e-18);
            close(
                null.coupling.lubrication.body_spin_impulses[a][i],
                0.0,
                8e-18,
            );
        }
    }
    assert!(null.coupling.lubrication.dissipated_energy < 1e-33);
    let (swapped, reversed) = frozen(
        [initial[1], initial[0]],
        [radii[1], radii[0]],
        [[0.0; 6]; 2],
        model(MU, 1),
        &[],
    );
    for a in 0..2 {
        for i in 0..3 {
            close(swapped[1 - a].velocity[i], after[a].velocity[i], 4e-17);
            close(
                swapped[1 - a].angular_velocity[i],
                after[a].angular_velocity[i],
                4e-17,
            );
            close(
                reversed.coupling.lubrication.body_spin_impulses[1 - a][i],
                report.coupling.lubrication.body_spin_impulses[a][i],
                1e-16,
            );
        }
    }
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.31, -0.27, 0.19))).unwrap();
    let rotated = initial.map(|b| {
        let r = rotation.matrix();
        let inertia = r
            .compose(Mat3 {
                m: core::array::from_fn(|i| b.inertia_world[i / 3][i % 3]),
            })
            .compose(r.transpose());
        let center = array(rotation.rotate(vector(b.center)).unwrap());
        PartialSaturationBody3 {
            center: core::array::from_fn(|i| center[i] + [0.4, -0.3, 0.2][i]),
            velocity: array(rotation.rotate(vector(b.velocity)).unwrap()),
            angular_velocity: array(rotation.rotate(vector(b.angular_velocity)).unwrap()),
            inertia_world: core::array::from_fn(|i| {
                core::array::from_fn(|j| inertia.m[3 * i.max(j) + i.min(j)])
            }),
            ..b
        }
    });
    let (changed, rr) = frozen(rotated, radii, [[0.0; 6]; 2], model(MU, 1), &[]);
    for a in 0..2 {
        let v = array(rotation.rotate(vector(after[a].velocity)).unwrap());
        let w = array(rotation.rotate(vector(after[a].angular_velocity)).unwrap());
        let j = array(
            rotation
                .rotate(vector(report.coupling.lubrication.body_spin_impulses[a]))
                .unwrap(),
        );
        for i in 0..3 {
            close(changed[a].velocity[i], v[i], 8e-16);
            close(changed[a].angular_velocity[i], w[i], 8e-16);
            close(
                rr.coupling.lubrication.body_spin_impulses[a][i],
                j[i],
                8e-16,
            );
        }
    }
    close(
        rr.coupling.lubrication.dissipated_energy,
        report.coupling.lubrication.dissipated_energy,
        2e-17,
    );
}

#[test]
fn shared_unequal_pairs_and_driven_wall_match_full_physical_inertia_system() {
    let radii = [0.4, 0.55, 0.7];
    let x = 3.0 + (0.952_f64.powi(2) - 0.15_f64.powi(2)).sqrt();
    let y = 3.0 + (1.253_f64.powi(2) - 0.15_f64.powi(2)).sqrt();
    let mut initial = [
        body([3.0, 3.0, 2.402], 1.3),
        body([x, 3.0, 2.552], 2.1),
        body([x, y, 2.702], 1.7),
    ];
    for (a, b) in initial.iter_mut().enumerate() {
        b.velocity = [0.004 - 0.003 * a as f64, 0.001 * a as f64, -0.0004];
        b.angular_velocity = [-0.001, 0.003 - 0.002 * a as f64, 0.002];
    }
    let planes = [LubricationPlane3::new(
        [0.0, 0.0, 2.0],
        [0.0, 0.0, 1.0],
        [0.001, -0.002, 0.0],
        NORMAL_CUT,
        ROUGH,
    )
    .unwrap()
    .with_tangential(CUT, ROUGH)
    .unwrap()];
    let external = core::array::from_fn(|a| {
        [
            0.0002,
            -0.0001 * a as f64,
            0.0001,
            0.00003,
            -0.00002,
            0.00001,
        ]
    });
    let expected = direct_step(&initial, radii, external, &planes);
    let (after, result) = frozen(initial, radii, external, model(MU, 5), &planes);
    let film = result.coupling.lubrication;
    assert_eq!((film.active_pairs, film.active_tangential_pairs), (2, 2));
    assert_eq!(
        (film.active_plane_pairs, film.active_tangential_plane_pairs),
        (3, 3)
    );
    assert_eq!(film.active_resistance_modes, 25);
    let v = speeds(&after);
    let old = speeds(&initial);
    let midpoint = core::array::from_fn(|a| core::array::from_fn(|i| 0.5 * (v[a][i] + old[a][i])));
    let actual = physical_loads(&initial, radii, &midpoint, MU, &planes);
    for a in 0..3 {
        for i in 0..6 {
            close(v[a][i], expected[6 * a + i], 3e-15);
            close(
                if i < 3 {
                    film.body_impulses[a][i]
                } else {
                    film.body_spin_impulses[a][i - 3]
                },
                actual[a][i],
                3e-15,
            );
        }
    }
    for i in 0..3 {
        close(
            film.body_impulses.iter().map(|j| j[i]).sum::<f64>() + film.plane_impulse[i],
            0.0,
            3e-17,
        );
        close(
            film.body_angular_impulses.iter().map(|j| j[i]).sum::<f64>()
                + film.plane_angular_impulse[i],
            0.0,
            3e-17,
        );
    }
    assert!(film.dissipated_energy > 0.0);
    close(film.work_residual, 0.0, 1e-17);
    close(
        result.planes[0].work,
        dot(result.planes[0].impulse, planes[0].velocity()),
        0.0,
    );
}

fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                result.push(if g.is_solid(x, y, z) {
                    [0.0; Q3]
                } else {
                    g.populations(x, y, z)
                });
            }
        }
    }
    result
}
fn fluid_delta(before: &[[f64; Q3]], g: &BoundaryGrid3) -> ([f64; 3], [f64; 3]) {
    let (mut p, mut h) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if g.is_solid(x, y, z) {
                    continue;
                }
                let f = g.populations(x, y, z);
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    for (i, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[i] += f64::from(e) * (f[q] - before[(z * 8 + y) * 8 + x][q]);
                    }
                }
                let orbital = cross(
                    [
                        x as f64 + 0.5 - ORIGIN[0],
                        y as f64 + 0.5 - ORIGIN[1],
                        z as f64 + 0.5 - ORIGIN[2],
                    ],
                    local,
                );
                for i in 0..3 {
                    p[i] += local[i];
                    h[i] += orbital[i];
                }
            }
        }
    }
    (p, h)
}
fn moving(a: usize) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: if a == 0 {
            [2.502, 3.4, 4.0]
        } else {
            [
                2.602,
                3.4 + (1.102_f64.powi(2) - 0.1_f64.powi(2)).sqrt(),
                4.0,
            ]
        },
        orientation: So3::identity(),
        mass: 1.3 + 0.2 * a as f64,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
            ],
        })
        .unwrap(),
        velocity: [0.0, 0.000002 - 0.000004 * a as f64, -0.000001],
        angular_momentum_body: Vec3::new(0.000002, -0.000001, 0.000001),
    }
}
fn momentum(b: MovingPartialSaturationBody3) -> ([f64; 3], [f64; 3]) {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(sub(b.center, ORIGIN), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}

#[test]
fn sampled_unequal_spinning_spheres_close_global_momentum_and_actual_wrench_power() {
    let mut g = grid(true);
    g.voxelize_sdf(|p| p[0] - 2.0);
    let force = [-2e-9, 1e-9, 0.5e-9];
    g.set_force_field(|_| force);
    let (mut forcing, mut moment) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if !g.is_solid(x, y, z) {
                    let torque = cross(
                        [
                            x as f64 + 0.5 - ORIGIN[0],
                            y as f64 + 0.5 - ORIGIN[1],
                            z as f64 + 0.5 - ORIGIN[2],
                        ],
                        force,
                    );
                    for i in 0..3 {
                        forcing[i] += force[i];
                        moment[i] += torque[i];
                    }
                }
            }
        }
    }
    let planes = [LubricationPlane3::new(
        [2.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0; 3],
        NORMAL_CUT,
        ROUGH,
    )
    .unwrap()
    .with_tangential(CUT, ROUGH)
    .unwrap()];
    let radii = [0.5, 0.6];
    let mut bodies = [moving(0), moving(1)];
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 4).unwrap());
    let external = [[-1e-8, 4e-9, -2e-9, 1e-9, -2e-9, 1e-9]; 2];
    let mut maximum = 0.0_f64;
    let mut loss = 0.0;
    for _ in 0..24 {
        let before = snapshot(&g);
        let old = bodies;
        let result = g
            .step_with_lubricated_spheres_and_planes(
                &mut samplers,
                radii,
                &mut bodies,
                external,
                model(MU, 3),
                &planes,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        let report = result.motion;
        let film = report.coupling.lubrication;
        assert_eq!(film.active_tangential_pairs, 1);
        assert_eq!(film.active_resistance_modes, 15);
        let frozen = core::array::from_fn(|a| body(report.bodies[a].midpoint_center, old[a].mass));
        let midpoint = core::array::from_fn(|a| {
            let b = report.coupling.bodies[a];
            [
                b.velocity_midpoint[0],
                b.velocity_midpoint[1],
                b.velocity_midpoint[2],
                b.angular_velocity_midpoint[0],
                b.angular_velocity_midpoint[1],
                b.angular_velocity_midpoint[2],
            ]
        });
        let load = physical_loads(&frozen, radii, &midpoint, MU, &planes);
        let (mut dp, mut dh) = fluid_delta(&before, &g);
        let mut work = 0.0;
        for a in 0..2 {
            let p0 = momentum(old[a]);
            let p1 = momentum(bodies[a]);
            let center = core::array::from_fn(|i| old[a].center[i] + 0.5 * old[a].velocity[i]);
            let torque = cross(sub(center, ORIGIN), external[a][..3].try_into().unwrap());
            for i in 0..3 {
                dp[i] += p1.0[i] - p0.0[i] - external[a][i];
                dh[i] += p1.1[i] - p0.1[i] - external[a][i + 3] - torque[i];
                close(film.body_impulses[a][i], load[a][i], 5e-16);
                close(film.body_spin_impulses[a][i], load[a][i + 3], 5e-17);
                close(report.bodies[a].momentum_residual[i], 0.0, 2e-15);
                close(report.bodies[a].momentum_residual[i + 3], 0.0, 2e-15);
            }
            work += load[a]
                .into_iter()
                .zip(midpoint[a])
                .map(|(f, v)| f * v)
                .sum::<f64>();
            close(report.bodies[a].work_residual, 0.0, 3e-17);
        }
        close(work + film.dissipated_energy, 0.0, 2e-19);
        close(film.work_residual, 0.0, 2e-19);
        loss += film.dissipated_energy;
        for i in 0..3 {
            let p =
                dp[i] + report.coupling.walls.wall_impulse[i] + film.plane_impulse[i] - forcing[i];
            let h = dh[i]
                + report.coupling.walls.wall_angular_impulse[i]
                + film.plane_angular_impulse[i]
                - moment[i];
            close(p, 0.0, 3e-12);
            close(h, 0.0, 9e-12);
            maximum = maximum.max(p.abs()).max(h.abs());
        }
    }
    assert!(loss > 0.0);
    eprintln!(
        "24 unequal pair/wall spin ticks: maximum global momentum defect={maximum:e}; modeled loss={loss:e}"
    );
}

#[test]
fn separate_scales_budgets_and_late_two_spin_refusal_preserve_atomic_state() {
    let base = SphereLubrication3::new(MU, NORMAL_CUT, ROUGH, 1).unwrap();
    let mut initial = [body([3.0, 4.0, 4.0], 1.3), body([4.02, 4.0, 4.0], 1.5)];
    initial[0].velocity = [0.009, -0.003, 0.002];
    initial[1].angular_velocity = [0.001, -0.002, 0.003];
    let a = frozen(initial, [0.5; 2], [[0.0; 6]; 2], base, &[]);
    let b = frozen(
        initial,
        [0.5; 2],
        [[0.0; 6]; 2],
        base.with_tangential_pairs(CUT, ROUGH).unwrap(),
        &[],
    );
    assert_eq!(a, b);
    let near = [
        initial[0],
        PartialSaturationBody3 {
            center: [4.002, 4.0, 4.0],
            ..initial[1]
        },
    ];
    let mut g = grid(false);
    let before = snapshot(&g);
    let mut state = near;
    let fields = [vec![0.0; 512], vec![0.0; 512]];
    for m in [
        model(MU, 1).with_resistance_mode_limit(4).unwrap(),
        base.with_tangential_pairs(0.01001, ROUGH).unwrap(),
    ] {
        assert!(
            g.step_with_lubricated_partial_saturation_bodies(
                fields.each_ref().map(Vec::as_slice),
                [0.5; 2],
                &mut state,
                [[0.0; 6]; 2],
                m,
                ORIGIN,
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(state, near);
        assert_eq!(snapshot(&g), before);
    }
    let (_, admitted) = frozen(
        near,
        [0.5; 2],
        [[0.0; 6]; 2],
        model(MU, 1).with_resistance_mode_limit(5).unwrap(),
        &[],
    );
    assert_eq!(admitted.coupling.lubrication.active_resistance_modes, 5);
    let tangent_only = SphereLubrication3::new(MU, 0.0001, ROUGH, 1)
        .unwrap()
        .with_tangential_pairs(CUT, 0.003)
        .unwrap();
    let mut spin = near;
    for b in &mut spin {
        b.velocity = [0.0; 3];
    }
    spin[0].angular_velocity = [0.0, 0.0, 0.01];
    spin[1].angular_velocity = [0.0, 0.0, -0.01];
    let (_, rough) = frozen(spin, [0.5; 2], [[0.0; 6]; 2], tangent_only, &[]);
    assert_eq!(rough.coupling.lubrication.active_resistance_modes, 4);
    // Equal-radius pure counterrotation has no line-of-centers contribution,
    // so both actual gaps below declared roughness must give identical updates.
    let mut wider = spin;
    wider[1].center[0] = 4.0025;
    assert_eq!(
        frozen(spin, [0.5; 2], [[0.0; 6]; 2], tangent_only, &[])
            .0
            .map(|b| b.angular_velocity),
        frozen(wider, [0.5; 2], [[0.0; 6]; 2], tangent_only, &[])
            .0
            .map(|b| b.angular_velocity)
    );
    assert!(base.with_tangential_pairs(f64::NAN, ROUGH).is_err());
    assert!(base.with_tangential_pairs(CUT, CUT * 2.0).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        g.step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 2],
            &mut state,
            [[0.0; 6]; 2],
            model(MU, 1),
            ORIGIN,
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(state, near);
    assert_eq!(snapshot(&g), before);
    let mut moving = [moving(0), moving(1)];
    moving[0].center = [3.0, 4.0, 4.0];
    moving[1].center = [4.102, 4.0, 4.0];
    for b in &mut moving {
        b.velocity = [0.0; 3];
    }
    let old = moving;
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 3).unwrap());
    let weak = model(1e-10, 1);
    assert!(matches!(
        g.step_with_lubricated_spheres(
            &mut samplers,
            [0.5, 0.6],
            &mut moving,
            [
                [0.04, 0.0, 0.0, 0.0, 0.001, 0.0],
                [-0.04, 0.0, 0.0, 0.0, -0.001, 0.0]
            ],
            weak,
            ORIGIN,
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Invalid {
            reason: "moving assembly swept enclosing spheres meet"
        })
    ));
    assert_eq!(moving, old);
    assert_eq!(snapshot(&g), before);
    let actual = g
        .step_with_lubricated_spheres(
            &mut samplers,
            [0.5, 0.6],
            &mut moving,
            [[0.0; 6]; 2],
            weak,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual.coupling.lubrication.active_tangential_pairs, 1);
    let mut fresh = grid(false);
    let mut expected = old;
    let mut fresh_samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&fresh, 3).unwrap());
    let repeated = fresh
        .step_with_lubricated_spheres(
            &mut fresh_samplers,
            [0.5, 0.6],
            &mut expected,
            [[0.0; 6]; 2],
            weak,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual, repeated);
    assert_eq!(moving, expected);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
