//! Independent physical resistance blocks and native angular-momentum coupling.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3, LubricatedPlanesExchange3,
    LubricationPlane3, MovingPartialSaturationBody3, PartialSaturationBody3,
    PartialSaturationError3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
const MU: f64 = 0.11;
const NORMAL_CUTOFF: f64 = 0.04;
const TANGENT_CUTOFF: f64 = 0.01;
const ROUGHNESS: f64 = 0.0001;
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.17e} != {b:.17e}: {:.3e}", a - b);
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
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
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
        inertia_world: [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    }
}
fn wall(point: [f64; 3], normal: [f64; 3], velocity: [f64; 3]) -> LubricationPlane3 {
    LubricationPlane3::new(point, normal, velocity, NORMAL_CUTOFF, ROUGHNESS)
        .unwrap()
        .with_tangential(TANGENT_CUTOFF, ROUGHNESS)
        .unwrap()
}
fn normal_resistance(mu: f64, r: f64, gap: f64, cutoff: f64, roughness: f64) -> f64 {
    if gap.max(roughness) >= cutoff {
        0.0
    } else {
        6.0 * core::f64::consts::PI * mu * r * r * (1.0 / gap.max(roughness) - 1.0 / cutoff)
    }
}
fn logarithmic_k(mu: f64, r: f64, gap: f64, cutoff: f64, roughness: f64) -> f64 {
    if gap.max(roughness) >= cutoff {
        0.0
    } else {
        (4.0 / 5.0) * core::f64::consts::PI * mu * r * (cutoff / gap.max(roughness)).ln()
    }
}
fn frozen<const N: usize>(
    mut bodies: [PartialSaturationBody3; N],
    radii: [f64; N],
    external: [[f64; 6]; N],
    model: SphereLubrication3,
    planes: &[LubricationPlane3],
) -> ([PartialSaturationBody3; N], LubricatedPlanesExchange3<N>) {
    // This zero-occupancy operator limit isolates the actual hydrodynamic
    // resistance block from resolved fluid drag. Sampled trajectories follow.
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

#[test]
fn analytic_translation_rotation_block_retains_reciprocity_and_support_couple() {
    let radius = 0.5;
    let k = logarithmic_k(
        MU,
        radius,
        (2.505 - 2.0) - radius,
        TANGENT_CUTOFF,
        ROUGHNESS,
    );
    let mass = 2.0;
    let inertia = 0.2;
    let speed = 0.003;
    for (u, omega) in [(0.012, -0.009), (0.012, 0.0), (0.0, 0.009)] {
        let mut initial = body([4.0, 4.0, 2.505], mass);
        initial.velocity = [u, 0.0, 0.0];
        initial.angular_velocity = [0.0, omega, 0.006];
        let external = [0.001, 0.0, 0.0, 0.0, -0.00002, 0.0];
        let plane = wall([0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [speed, 0.0, 0.0]);
        // Solve physical velocity/spin coordinates directly, without production
        // mode factorization or the impulse-space Schur complement.
        let aa = mass + 2.0 * k;
        let dd = inertia + 2.0 * k * radius * radius;
        let ab = -0.5 * k * radius;
        let x = (mass - 2.0 * k) * u + 0.5 * k * radius * omega + external[0] + 4.0 * k * speed;
        let y = 0.5 * k * radius * u + (inertia - 2.0 * k * radius * radius) * omega + external[4]
            - k * radius * speed;
        let expected = [
            (dd * x - ab * y) / (aa * dd - ab * ab),
            (aa * y - ab * x) / (aa * dd - ab * ab),
        ];
        let (after, result) = frozen(
            [initial],
            [radius],
            [external],
            SphereLubrication3::new(MU, NORMAL_CUTOFF, ROUGHNESS, 1).unwrap(),
            &[plane],
        );
        close(after[0].velocity[0], expected[0], 7e-17);
        close(after[0].angular_velocity[1], expected[1], 7e-17);
        close(
            after[0].angular_velocity[2],
            initial.angular_velocity[2],
            0.0,
        );
        let film = result.coupling.lubrication;
        assert_eq!(film.active_plane_pairs, 1);
        assert_eq!(film.active_tangential_plane_pairs, 1);
        assert_eq!(film.active_resistance_modes, 5);
        assert_eq!(result.planes[0].tangential_spheres, 1);
        let vmid = 0.5 * (u + after[0].velocity[0]) - speed;
        let wmid = 0.5 * (omega + after[0].angular_velocity[1]);
        let force = -4.0 * k * vmid + k * radius * wmid;
        let torque = k * radius * vmid - 4.0 * k * radius * radius * wmid;
        close(film.body_impulses[0][0], force, 7e-17);
        close(film.body_spin_impulses[0][1], torque, 7e-17);
        let foot = 2.0;
        let arm = initial.center[2] - foot;
        close(
            result.planes[0].couple_impulse[1],
            -torque - arm * force,
            7e-17,
        );
        assert!(result.planes[0].couple_impulse[1].abs() > 1e-5);
        for i in 0..3 {
            close(
                film.body_impulses[0][i] + result.planes[0].impulse[i],
                0.0,
                3e-17,
            );
            close(
                film.body_angular_impulses[0][i] + result.planes[0].angular_impulse[i],
                0.0,
                1e-16,
            );
            close(result.coupling.bodies[0].momentum_residual[i], 0.0, 2e-17);
            close(
                result.coupling.bodies[0].momentum_residual[i + 3],
                0.0,
                2e-17,
            );
        }
        let loss = k
            * (4.0 * vmid * vmid - 2.0 * radius * vmid * wmid
                + 4.0 * radius * radius * wmid * wmid);
        close(film.dissipated_energy, loss, 2e-18);
        close(film.body_work[0] + result.planes[0].work + loss, 0.0, 2e-18);
        close(film.work_residual, 0.0, 2e-18);
    }
}

#[test]
fn rolling_and_opposed_modes_are_dissipative_and_converge_at_second_order() {
    let r = 0.5;
    let mass = 2.0;
    let k = logarithmic_k(MU, r, (2.505 - 2.0) - r, TANGENT_CUTOFF, ROUGHNESS);
    let plane = wall([0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [0.0; 3]);
    for sign in [-1.0, 1.0] {
        let mut initial = body([4.0, 4.0, 2.505], mass);
        // A thin circular ring about y has Iyy=mR², Ixx=Izz=mR²/2.
        initial.inertia_world = [[0.25, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.25]];
        initial.velocity = [0.01, 0.0, 0.0];
        initial.angular_velocity = [0.0, sign * 0.01 / r, 0.0];
        let eigenvalue = if sign > 0.0 { 3.0 * k } else { 5.0 * k };
        let gain = (mass - 0.5 * eigenvalue) / (mass + 0.5 * eigenvalue);
        let (after, result) = frozen(
            [initial],
            [r],
            [[0.0; 6]],
            SphereLubrication3::new(MU, NORMAL_CUTOFF, ROUGHNESS, 1).unwrap(),
            &[plane],
        );
        close(after[0].velocity[0], gain * initial.velocity[0], 3e-17);
        close(
            after[0].angular_velocity[1],
            gain * initial.angular_velocity[1],
            6e-17,
        );
        assert!(result.coupling.lubrication.dissipated_energy > 0.0);
        if sign > 0.0 {
            close(
                initial.velocity[0] - r * initial.angular_velocity[1],
                0.0,
                0.0,
            );
            assert!(
                after[0].velocity[0] < initial.velocity[0],
                "zero material-point slip still has pressure-driven rolling resistance"
            );
        }
    }
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = 1.0 / steps as f64;
        let mut initial = body([4.0, 4.0, 2.505], mass);
        initial.inertia_world = [[0.25, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.25]];
        // Lattice dt=1 represents physical dt here: velocities and dynamic
        // viscosity both scale with that time unit; length/mass units stay fixed.
        initial.velocity = [0.01 * dt, 0.0, 0.0];
        initial.angular_velocity = [0.0, 0.01 * dt / r, 0.0];
        let mut state = [initial];
        let model = SphereLubrication3::new(MU * dt, NORMAL_CUTOFF, ROUGHNESS, 1).unwrap();
        for _ in 0..steps {
            state = frozen(state, [r], [[0.0; 6]], model, &[plane]).0;
        }
        let expected = 0.01 * (-3.0 * k / mass).exp();
        errors.push((state[0].velocity[0] / dt - expected).abs());
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((3.9..4.1).contains(&ratio));
    }
    eprintln!("tangential rolling-mode refinement errors: {errors:?}");
}

// Independent direct physical force/torque law, with no modal factorization.
fn physical_loads<const N: usize>(
    bodies: &[PartialSaturationBody3; N],
    radii: [f64; N],
    planes: &[LubricationPlane3],
    velocity: &[[f64; 6]; N],
) -> [[f64; 6]; N] {
    let mut result = [[0.0; 6]; N];
    for a in 0..N {
        let v: [f64; 3] = velocity[a][..3].try_into().unwrap();
        let w: [f64; 3] = velocity[a][3..].try_into().unwrap();
        for plane in planes {
            let n = plane.normal();
            let gap = dot(n, sub(bodies[a].center, plane.point())) - radii[a];
            let u = sub(v, plane.velocity());
            let un = dot(u, n);
            let wn = dot(w, n);
            let ut = core::array::from_fn::<_, 3, _>(|i| u[i] - un * n[i]);
            let wt = core::array::from_fn::<_, 3, _>(|i| w[i] - wn * n[i]);
            let normal = normal_resistance(MU, radii[a], gap, NORMAL_CUTOFF, ROUGHNESS);
            let k = logarithmic_k(MU, radii[a], gap, TANGENT_CUTOFF, ROUGHNESS);
            let coupling_force = cross(w, n);
            let coupling_torque = cross(n, u);
            for i in 0..3 {
                result[a][i] +=
                    -normal * un * n[i] - 4.0 * k * ut[i] + k * radii[a] * coupling_force[i];
                result[a][i + 3] +=
                    k * radii[a] * coupling_torque[i] - 4.0 * k * radii[a] * radii[a] * wt[i];
            }
        }
    }
    for a in 0..N {
        for b in a + 1..N {
            let delta = sub(bodies[b].center, bodies[a].center);
            let distance = dot(delta, delta).sqrt();
            let n = delta.map(|x| x / distance);
            let effective = radii[a] * radii[b] / (radii[a] + radii[b]);
            let resistance = normal_resistance(
                MU,
                effective,
                distance - radii[a] - radii[b],
                NORMAL_CUTOFF,
                ROUGHNESS,
            );
            let relative = core::array::from_fn(|i| velocity[a][i] - velocity[b][i]);
            let force = -resistance * dot(n, relative);
            for i in 0..3 {
                result[a][i] += force * n[i];
                result[b][i] -= force * n[i];
            }
        }
    }
    result
}
fn solve_dense(mut a: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    // Test-only Gaussian elimination on the independently assembled physical
    // velocity equations, distinct from production's modal Cholesky system.
    let n = rhs.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap();
        a.swap(col, pivot);
        rhs.swap(col, pivot);
        assert!(a[col][col].abs() > 1e-12);
        for row in col + 1..n {
            let scale = a[row][col] / a[col][col];
            for j in col..n {
                a[row][j] -= scale * a[col][j];
            }
            rhs[row] -= scale * rhs[col];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] -= a[i][j] * rhs[j];
        }
        rhs[i] /= a[i][i];
    }
    rhs
}

#[test]
fn mixed_wall_and_sphere_modes_match_independent_full_tensor_system_in_rotated_frames() {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.31, -0.27, 0.19))).unwrap();
    for transformed in [false, true] {
        let rotate = |v| {
            if transformed {
                array(rotation.rotate(vector(v)).unwrap())
            } else {
                v
            }
        };
        let boost = if transformed {
            rotate([0.011, -0.007, 0.005])
        } else {
            [0.0; 3]
        };
        let initial = core::array::from_fn::<_, 2, _>(|a| {
            let mut b = body(
                rotate([2.505 + 1.01 * a as f64, 4.0, 2.505]),
                1.3 + 0.2 * a as f64,
            );
            let tensor = Mat3 {
                m: [
                    0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
                ],
            };
            let tensor = if transformed {
                rotation
                    .matrix()
                    .compose(tensor)
                    .compose(rotation.matrix().transpose())
            } else {
                tensor
            };
            b.inertia_world = core::array::from_fn(|i| {
                core::array::from_fn(|j| tensor.m[3 * i.max(j) + i.min(j)])
            });
            b.velocity = rotate([0.00002 - 0.00003 * a as f64, 0.00001, -0.00002]);
            for i in 0..3 {
                b.velocity[i] += boost[i];
            }
            b.angular_velocity = rotate([0.00006, -0.00004, 0.00002 + 0.00002 * a as f64]);
            b
        });
        let mut planes = [
            wall(rotate([0.0, 0.0, 2.0]), rotate([0.0, 0.0, 1.0]), boost),
            wall(rotate([2.0, 0.0, 0.0]), rotate([1.0, 0.0, 0.0]), boost),
        ];
        if transformed {
            planes.swap(0, 1);
        }
        let external = core::array::from_fn::<_, 2, _>(|a| {
            let force = rotate([0.000001, -0.000002, 0.000001 * (a + 1) as f64]);
            let torque = rotate([-0.0000001, 0.0000002, 0.0000001]);
            [
                force[0], force[1], force[2], torque[0], torque[1], torque[2],
            ]
        });
        let start = initial.map(|b| {
            [
                b.velocity[0],
                b.velocity[1],
                b.velocity[2],
                b.angular_velocity[0],
                b.angular_velocity[1],
                b.angular_velocity[2],
            ]
        });
        let zero = physical_loads(&initial, [0.5; 2], &planes, &[[0.0; 6]; 2]);
        let old_load = physical_loads(&initial, [0.5; 2], &planes, &start);
        let mut resistance = vec![vec![0.0; 12]; 12];
        for column in 0..12 {
            let mut unit = [[0.0; 6]; 2];
            unit[column / 6][column % 6] = 1.0;
            let load = physical_loads(&initial, [0.5; 2], &planes, &unit);
            for row in 0..12 {
                resistance[row][column] = zero[row / 6][row % 6] - load[row / 6][row % 6];
            }
        }
        let mut system = vec![vec![0.0; 12]; 12];
        let mut rhs = vec![0.0; 12];
        for a in 0..2 {
            for i in 0..6 {
                let row = 6 * a + i;
                let old_momentum = if i < 3 {
                    initial[a].mass * start[a][i]
                } else {
                    dot(
                        initial[a].inertia_world[i - 3],
                        start[a][3..].try_into().unwrap(),
                    )
                };
                rhs[row] = old_momentum + external[a][i] + 0.5 * (old_load[a][i] + zero[a][i]);
                for column in 0..12 {
                    system[row][column] = 0.5 * resistance[row][column];
                }
                if i < 3 {
                    system[row][row] += initial[a].mass;
                } else {
                    for j in 0..3 {
                        system[row][6 * a + 3 + j] += initial[a].inertia_world[i - 3][j];
                    }
                }
            }
        }
        let expected = solve_dense(system, rhs);
        let (after, result) = frozen(
            initial,
            [0.5; 2],
            external,
            SphereLubrication3::new(MU, NORMAL_CUTOFF, ROUGHNESS, 4).unwrap(),
            &planes,
        );
        assert_eq!(result.coupling.lubrication.active_resistance_modes, 16);
        assert_eq!(result.coupling.lubrication.active_pairs, 1);
        assert_eq!(result.coupling.lubrication.active_tangential_plane_pairs, 3);
        for a in 0..2 {
            for i in 0..3 {
                close(after[a].velocity[i], expected[6 * a + i], 3e-16);
                close(after[a].angular_velocity[i], expected[6 * a + i + 3], 5e-16);
            }
        }
        let film = result.coupling.lubrication;
        close(film.work_residual, 0.0, 2e-18);
        for i in 0..3 {
            close(
                film.body_impulses[0][i] + film.body_impulses[1][i] + film.plane_impulse[i],
                0.0,
                7e-17,
            );
            close(
                film.body_angular_impulses[0][i]
                    + film.body_angular_impulses[1][i]
                    + film.plane_angular_impulse[i],
                0.0,
                3e-16,
            );
        }
    }
}

fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut out = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                out.push(if g.is_solid(x, y, z) {
                    [0.0; Q3]
                } else {
                    g.populations(x, y, z)
                });
            }
        }
    }
    out
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
                    let d = f[q] - before[(z * 8 + y) * 8 + x][q];
                    for (i, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[i] += f64::from(e) * d;
                    }
                }
                let torque = cross(
                    [
                        x as f64 + 0.5 - ORIGIN[0],
                        y as f64 + 0.5 - ORIGIN[1],
                        z as f64 + 0.5 - ORIGIN[2],
                    ],
                    local,
                );
                for i in 0..3 {
                    p[i] += local[i];
                    h[i] += torque[i];
                }
            }
        }
    }
    (p, h)
}
fn moving(a: usize) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [2.505, 3.4 + 1.01 * a as f64, 4.0],
        orientation: So3::identity(),
        mass: 1.3 + 0.2 * a as f64,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
            ],
        })
        .unwrap(),
        velocity: [0.0, 0.00001 - 0.00002 * a as f64, -0.000005],
        angular_momentum_body: Vec3::new(0.000008, -0.000006, 0.000004),
    }
}
fn momentum(b: MovingPartialSaturationBody3) -> ([f64; 3], [f64; 3]) {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(sub(b.center, ORIGIN), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}

#[test]
fn sampled_rotating_spheres_near_actual_wall_close_physical_torque_and_power() {
    let mut g = grid(true);
    g.voxelize_sdf(|p| p[0] - 2.0);
    let force = [-2e-9, 1e-9, 0.5e-9];
    g.set_force_field(|_| force);
    let (mut forcing, mut moment) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if g.is_solid(x, y, z) {
                    continue;
                }
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
    let planes = [wall([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3])];
    let mut bodies = [moving(0), moving(1)];
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 4).unwrap());
    let model = SphereLubrication3::new(MU, NORMAL_CUTOFF, ROUGHNESS, 3).unwrap();
    let external = [[
        -0.00000001,
        0.000000004,
        -0.000000002,
        0.000000001,
        -0.000000002,
        0.000000001,
    ]; 2];
    let mut maximum = 0.0_f64;
    let mut loss = 0.0;
    for _ in 0..24 {
        let before = snapshot(&g);
        let old = bodies;
        let result = g
            .step_with_lubricated_spheres_and_planes(
                &mut samplers,
                [0.5; 2],
                &mut bodies,
                external,
                model,
                &planes,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        let report = result.motion;
        let film = report.coupling.lubrication;
        assert_eq!(film.active_resistance_modes, 11);
        assert_eq!(film.active_tangential_plane_pairs, 2);
        let (mut dp, mut dh) = fluid_delta(&before, &g);
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
        let loads = physical_loads(&frozen, [0.5; 2], &planes, &midpoint);
        let mut independent_work = 0.0;
        for a in 0..2 {
            let old_momentum = momentum(old[a]);
            let new_momentum = momentum(bodies[a]);
            let center = core::array::from_fn(|i| old[a].center[i] + 0.5 * old[a].velocity[i]);
            let torque = cross(sub(center, ORIGIN), external[a][..3].try_into().unwrap());
            for i in 0..3 {
                dp[i] += new_momentum.0[i] - old_momentum.0[i] - external[a][i];
                dh[i] += new_momentum.1[i] - old_momentum.1[i] - external[a][i + 3] - torque[i];
                close(film.body_impulses[a][i], loads[a][i], 5e-16);
                close(film.body_spin_impulses[a][i], loads[a][i + 3], 5e-17);
                close(report.bodies[a].momentum_residual[i], 0.0, 2e-15);
                close(report.bodies[a].momentum_residual[i + 3], 0.0, 2e-15);
            }
            independent_work += loads[a]
                .into_iter()
                .zip(midpoint[a])
                .map(|(f, v)| f * v)
                .sum::<f64>();
            close(report.bodies[a].work_residual, 0.0, 3e-17);
        }
        close(independent_work + film.dissipated_energy, 0.0, 3e-19);
        close(film.work_residual, 0.0, 3e-19);
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
        "24 full wall-resistance ticks: max world momentum defect={maximum:e}, modeled loss={loss:e}"
    );
}

#[test]
fn independent_tangent_cutoff_roughness_and_normal_only_fast_path_are_retained() {
    let model = SphereLubrication3::new(MU, 0.1, 0.0001, 1).unwrap();
    let mut initial = body([4.0, 4.0, 2.52], 2.0);
    initial.velocity = [0.01, 0.0, 0.0];
    initial.angular_velocity = [0.0, 0.02, 0.0];
    let normal =
        LubricationPlane3::new([0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [0.0; 3], 0.1, 0.0001).unwrap();
    let (a, ra) = frozen([initial], [0.5], [[0.0; 6]], model, &[normal]);
    let (b, rb) = frozen(
        [initial],
        [0.5],
        [[0.0; 6]],
        model,
        &[normal.with_tangential(0.01, 0.0001).unwrap()],
    );
    assert_eq!(a, b);
    assert_eq!(ra, rb);
    let mut initial = body([4.0, 4.0, 3.0], 2.0);
    initial.velocity = [0.01, 0.0, 0.0];
    initial.angular_velocity = [0.0, 0.02, 0.0];
    let mut previous = None;
    for gap in [1.0 / 512.0, 1.0 / 256.0] {
        let plane = LubricationPlane3::new(
            [0.0, 0.0, 2.0 - gap],
            [0.0, 0.0, 1.0],
            [0.0; 3],
            1.0 / 1024.0,
            1.0 / 2048.0,
        )
        .unwrap()
        .with_tangential(1.0 / 64.0, 1.0 / 128.0)
        .unwrap();
        let (after, report) = frozen([initial], [1.0], [[0.0; 6]], model, &[plane]);
        assert_eq!(report.coupling.lubrication.active_resistance_modes, 4);
        if let Some(old) = previous {
            assert_eq!(after, old);
        }
        previous = Some(after);
    }
}

#[test]
fn mode_budgets_asymptotic_admission_and_late_coupled_failure_are_atomic() {
    let mut g = grid(false);
    let before = snapshot(&g);
    let field = vec![0.0; 512];
    let initial = body([4.0, 4.0, 2.505], 2.0);
    let mut bodies = [initial];
    let base = SphereLubrication3::new(MU, NORMAL_CUTOFF, ROUGHNESS, 1).unwrap();
    let plane = wall([0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [0.0; 3]);
    assert!(
        g.step_with_lubricated_partial_saturation_planes(
            [&field],
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            base.with_resistance_mode_limit(4).unwrap(),
            &[plane],
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, [initial]);
    assert_eq!(snapshot(&g), before);
    let excessive = plane.with_tangential(0.011, ROUGHNESS).unwrap();
    assert!(
        g.step_with_lubricated_partial_saturation_planes(
            [&field],
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            base,
            &[excessive],
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, [initial]);
    assert_eq!(snapshot(&g), before);
    assert!(base.with_resistance_mode_limit(0).is_err());
    assert!(base.with_resistance_mode_limit(321).is_err());
    assert!(base.with_resistance_mode_limit(320).is_ok());
    assert!(plane.with_tangential(f64::NAN, ROUGHNESS).is_err());
    assert!(plane.with_tangential(0.01, 0.02).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        g.step_with_lubricated_partial_saturation_planes(
            [&field],
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            base,
            &[plane],
            ORIGIN,
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    let (_, accepted) = frozen(
        [initial],
        [0.5],
        [[0.0; 6]],
        base.with_resistance_mode_limit(5).unwrap(),
        &[plane],
    );
    assert_eq!(accepted.coupling.lubrication.active_resistance_modes, 5);
    let mut initial = moving(0);
    initial.velocity = [0.0; 3];
    initial.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    let mut bodies = [initial];
    let mut sampler = [BodyFrameSampler3::new(&g, 3).unwrap()];
    let weak = SphereLubrication3::new(1e-10, NORMAL_CUTOFF, ROUGHNESS, 1).unwrap();
    let driven = wall([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.007, 0.002, 0.0]);
    assert!(matches!(
        g.step_with_lubricated_spheres_and_planes(
            &mut sampler,
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            weak,
            &[driven],
            ORIGIN,
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Invalid {
            reason: "moving sphere meets prescribed lubrication plane"
        })
    ));
    assert_eq!(bodies, [initial]);
    assert_eq!(snapshot(&g), before);
    let driven = wall([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0001, 0.002, 0.0]);
    let accepted = g
        .step_with_lubricated_spheres_and_planes(
            &mut sampler,
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            weak,
            &[driven],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let mut fresh = grid(false);
    let mut expected = [initial];
    let mut fresh_sampler = [BodyFrameSampler3::new(&fresh, 3).unwrap()];
    let repeated = fresh
        .step_with_lubricated_spheres_and_planes(
            &mut fresh_sampler,
            [0.5],
            &mut expected,
            [[0.0; 6]],
            weak,
            &[driven],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(bodies, expected);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
