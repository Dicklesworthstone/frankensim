//! G0/G3: native Guo forcing, immersed recoil and full spatial balance.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3,
    MovingPartialSaturationBody3, PartialSaturationBody3, PartialSaturationError3, Q3, W3,
    equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const N: usize = 8;
const TAU: f64 = 0.83;
const ORIGIN: [f64; 3] = [-4.2, 7.1, 2.3];
const UNIFORM: [f64; 3] = [2e-5, -1.1e-5, 0.7e-5];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}: {}", a - b);
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn arr(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn index(c: [usize; 3]) -> usize {
    (c[2] * N + c[1]) * N + c[0]
}
fn field(c: [usize; 3]) -> [f64; 3] {
    [
        1.3e-6 * (c[1] as f64 - 3.5),
        -0.8e-6 * (c[2] as f64 - 3.5),
        1.7e-6 * (c[0] as f64 - 3.5),
    ]
}
fn grid(walls: bool) -> BoundaryGrid3 {
    let mut g = BoundaryGrid3::new(
        N,
        N,
        N,
        TAU,
        UNIFORM,
        if walls {
            BoundarySpec3::new([FaceBoundary3::Wall { velocity: [0.0; 3] }; 6])
        } else {
            BoundarySpec3::periodic()
        },
    );
    g.set_force_field(field);
    g
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut f = Vec::new();
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                f.push(g.populations(x, y, z));
            }
        }
    }
    f
}
fn fluid_delta(before: &[[f64; Q3]], g: &BoundaryGrid3) -> (f64, [f64; 3], [f64; 3]) {
    let (mut mass, mut p, mut h) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let mut local = [0.0; 3];
                let f = g.populations(x, y, z);
                for q in 0..Q3 {
                    let d = f[q] - before[index([x, y, z])][q];
                    mass += d;
                    for (i, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[i] += f64::from(e) * d;
                    }
                }
                let arm = [
                    x as f64 + 0.5 - ORIGIN[0],
                    y as f64 + 0.5 - ORIGIN[1],
                    z as f64 + 0.5 - ORIGIN[2],
                ];
                let angular = cross(arm, local);
                for i in 0..3 {
                    p[i] += local[i];
                    h[i] += angular[i];
                }
            }
        }
    }
    (mass, p, h)
}
fn imposed() -> ([f64; 3], [f64; 3]) {
    let (mut p, mut h) = ([0.0; 3], [0.0; 3]);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let f = field([x, y, z]);
                let f = core::array::from_fn(|i| f[i] + UNIFORM[i]);
                let arm = [
                    x as f64 + 0.5 - ORIGIN[0],
                    y as f64 + 0.5 - ORIGIN[1],
                    z as f64 + 0.5 - ORIGIN[2],
                ];
                let torque = cross(arm, f);
                for i in 0..3 {
                    p[i] += f[i];
                    h[i] += torque[i];
                }
            }
        }
    }
    (p, h)
}
fn body(a: usize) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [3.2 + 0.6 * a as f64, 3.7 - 0.3 * a as f64, 3.4],
        mass: 0.7 + 0.5 * a as f64,
        inertia_world: [[1.3, 0.12, -0.08], [0.12, 1.7, 0.06], [-0.08, 0.06, 2.0]],
        velocity: [-0.014 + 0.01 * a as f64, 0.008, 0.006],
        angular_velocity: [0.004, -0.003, 0.005 - 0.002 * a as f64],
    }
}
fn fields() -> [Vec<f64>; 2] {
    let mut f = core::array::from_fn(|_| vec![0.0; N * N * N]);
    for (c, a, b) in [
        ([4, 3, 3], 0.3, 0.6),
        ([3, 4, 3], 0.5, 0.2),
        ([2, 3, 4], 0.4, 0.0),
        ([4, 2, 2], 0.0, 0.7),
    ] {
        f[0][index(c)] = a;
        f[1][index(c)] = b;
    }
    f
}
// Independent Hermite source, with no call to native forcing/reframing helpers.
fn source(u: [f64; 3], force: [f64; 3]) -> [f64; Q3] {
    core::array::from_fn(|q| {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        W3[q] * (3.0 * dot(e, force) - 3.0 * dot(u, force) + 9.0 * dot(e, u) * dot(e, force))
    })
}

#[test]
fn forced_shared_cell_matches_independent_physical_population_formula() {
    let cell = [3, 3, 3];
    let force = [6e-4, -3e-4, 2e-4];
    let (rho, u) = (1.07, [0.024, -0.017, 0.011]);
    let mut g = BoundaryGrid3::new(N, N, N, TAU, [0.0; 3], BoundarySpec3::periodic());
    g.set_force_field(|c| if c == cell { force } else { [0.0; 3] });
    let eq = equilibrium3(rho, u);
    let s = source(u, force);
    // A real nonequilibrium stress with zero mass/momentum enters the oracle.
    let mut physical = eq;
    physical[1] += 2e-4;
    physical[2] += 2e-4;
    physical[0] -= 4e-4;
    let f = core::array::from_fn(|q| physical[q] - 0.5 * s[q]);
    g.set_populations(3, 3, 3, &f);
    for i in 0..3 {
        close(g.velocity(3, 3, 3)[i], u[i], 5e-17);
    }
    let mut fractions = [vec![0.0; N * N * N], vec![0.0; N * N * N]];
    fractions[0][index(cell)] = 0.3;
    fractions[1][index(cell)] = 0.6;
    let mut bodies = [body(0), body(1)];
    for b in &mut bodies {
        b.center = [3.5; 3];
        b.angular_velocity = [0.0; 3];
    }
    let before = bodies;
    let external = [
        [0.0003, -0.0002, 0.0001, 0.0, 0.0, 0.0],
        [-0.0001, 0.0004, 0.0, 0.0, 0.0, 0.0],
    ];
    let receipt = g
        .step_with_partial_saturation_bodies(
            fractions.each_ref().map(Vec::as_slice),
            &mut bodies,
            external,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let mut expected: [f64; Q3] =
        core::array::from_fn(|q| f[q] + (eq[q] - f[q]) / TAU + (1.0 - 0.5 / TAU) * s[q]);
    for a in 0..2 {
        let b = [0.3, 0.6][a] * (TAU - 0.5) / (1.0 - 0.9 + TAU - 0.5);
        let weight = b * rho;
        let after: [f64; 3] = core::array::from_fn(|i| {
            before[a].velocity[i]
                + (weight * (u[i] - before[a].velocity[i]) + external[a][i])
                    / (before[a].mass + 0.5 * weight)
        });
        let mid = core::array::from_fn(|i| 0.5 * (after[i] + before[a].velocity[i]));
        let solid = equilibrium3(rho, mid);
        for i in 0..3 {
            close(bodies[a].velocity[i], after[i], 2e-16);
        }
        for q in 0..Q3 {
            let opp = fs_lbm::d3q19::OPP3[q];
            let fluid = (eq[q] - physical[q]) / TAU;
            let solid = physical[opp] - physical[q] + solid[q] - eq[opp];
            expected[q] += b * (solid - fluid);
        }
        for d in receipt.bodies[a].momentum_residual {
            close(d, 0.0, 3e-16);
        }
    }
    let mut post_momentum = [0.0; 3];
    for q in 0..Q3 {
        let e = E3[q];
        let dest = [(3 + e.0) as usize, (3 + e.1) as usize, (3 + e.2) as usize];
        close(
            g.populations(dest[0], dest[1], dest[2])[q],
            expected[q],
            3e-16,
        );
        for (i, e) in [e.0, e.1, e.2].into_iter().enumerate() {
            post_momentum[i] += f64::from(e) * expected[q];
        }
    }
    let post_rho: f64 = expected.iter().sum();
    let mean =
        core::array::from_fn(|i| 0.5 * (u[i] + (post_momentum[i] + 0.5 * force[i]) / post_rho));
    close(receipt.forcing.midpoint_work, dot(force, mean), 2e-19);
    let torque = cross(core::array::from_fn(|i| 3.5 - ORIGIN[i]), force);
    for i in 0..3 {
        close(receipt.forcing.impulse[i], force[i], 0.0);
        close(receipt.forcing.angular_impulse[i], torque[i], 2e-18);
    }
}

#[test]
fn zero_occupancy_retains_native_guo_populations_bit_for_bit() {
    let (mut native, mut coupled) = (grid(false), grid(false));
    let empty = vec![0.0; N * N * N];
    let mut bodies = [body(0), body(1)];
    let initial = bodies;
    let expected = imposed();
    for _ in 0..12 {
        native.step();
        let report = coupled
            .step_with_partial_saturation_bodies(
                [&empty, &empty],
                &mut bodies,
                [[0.0; 6]; 2],
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(snapshot(&native), snapshot(&coupled));
        assert_eq!(bodies, initial);
        for i in 0..3 {
            close(report.forcing.impulse[i], expected.0[i], 4e-17);
            close(report.forcing.angular_impulse[i], expected.1[i], 3e-16);
        }
        assert_eq!(report.bodies[0].hydrodynamic.coupled_cells, 0);
    }
}

#[test]
fn forced_assembly_closes_independent_world_momentum_with_shared_cells_and_walls() {
    forced_assembly_balance(true);
}

#[test]
fn forced_periodic_assembly_closes_independent_world_angular_momentum() {
    forced_assembly_balance(false);
}

fn forced_assembly_balance(walls: bool) {
    let mut g = grid(walls);
    let f = fields();
    let mut bodies = [body(0), body(1)];
    let external = [
        [0.0005, -0.0003, 0.0002, -0.0004, 0.0006, 0.0003],
        [-0.0002, 0.0001, 0.0004, 0.0002, -0.0001, 0.0002],
    ];
    let imposed = imposed();
    let mut max_defect = 0.0_f64;
    let mut max_transport = 0.0_f64;
    for _ in 0..32 {
        let initial = snapshot(&g);
        let old = bodies;
        let report = g
            .step_with_partial_saturation_bodies(
                f.each_ref().map(Vec::as_slice),
                &mut bodies,
                external,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        let (mass, mut p, mut h) = fluid_delta(&initial, &g);
        close(mass, 0.0, 1e-12);
        for a in 0..2 {
            let dp = core::array::from_fn(|i| {
                old[a].mass * (bodies[a].velocity[i] - old[a].velocity[i])
            });
            let dw = core::array::from_fn(|i| {
                bodies[a].angular_velocity[i] - old[a].angular_velocity[i]
            });
            let dh = old[a].inertia_world.map(|row| dot(row, dw));
            let arm = core::array::from_fn(|i| old[a].center[i] - ORIGIN[i]);
            let orbital = cross(arm, dp);
            let ext_orbital = cross(arm, external[a][..3].try_into().unwrap());
            let mean_v =
                core::array::from_fn(|i| 0.5 * (old[a].velocity[i] + bodies[a].velocity[i]));
            let mean_w = core::array::from_fn(|i| {
                0.5 * (old[a].angular_velocity[i] + bodies[a].angular_velocity[i])
            });
            close(
                dot(dp, mean_v) + dot(dh, mean_w),
                report.bodies[a].external_work + report.bodies[a].hydrodynamic.body_work,
                4e-16,
            );
            for i in 0..3 {
                p[i] += dp[i] - external[a][i];
                h[i] += dh[i] + orbital[i] - external[a][i + 3] - ext_orbital[i];
            }
        }
        for i in 0..3 {
            let dp = p[i] + report.walls.wall_impulse[i] - imposed.0[i];
            let transport = report.periodic.angular_momentum_change[i];
            let dh = h[i] + report.walls.wall_angular_impulse[i] - imposed.1[i] - transport;
            max_transport = max_transport.max(transport.abs());
            max_defect = max_defect.max(dp.abs()).max(dh.abs());
            close(dp, 0.0, 2e-12);
            close(dh, 0.0, 9e-12);
            close(report.forcing.impulse[i], imposed.0[i], 4e-17);
            close(report.forcing.angular_impulse[i], imposed.1[i], 3e-16);
        }
        assert_eq!(report.walls.measured_links > 0, walls);
        assert_eq!(report.periodic.face_crossings > 0, !walls);
    }
    assert_eq!(max_transport > 1e-6, !walls);
    eprintln!("32 forced shared-cell ticks, maximum world balance defect {max_defect:e}");
}

fn moving(a: usize) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [3.8 + 0.72 * a as f64, 3.8, 4.2],
        orientation: So3::identity(),
        mass: 0.5 + 0.2 * a as f64,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.022, 0.003, -0.002, 0.003, 0.028, 0.001, -0.002, 0.001, 0.035,
            ],
        })
        .unwrap(),
        velocity: [0.001, 0.002, -0.001],
        angular_momentum_body: Vec3::new(0.00003, -0.00002, 0.00004),
    }
}
fn momentum(b: MovingPartialSaturationBody3) -> ([f64; 3], [f64; 3]) {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = arr(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| b.center[i] - ORIGIN[i]), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}

#[test]
fn moving_forced_full_tensor_bodies_keep_world_balance_and_single_force_ownership() {
    moving_forced_balance(true);
}

#[test]
fn moving_periodic_bodies_keep_world_angular_balance_and_single_transport_ownership() {
    moving_forced_balance(false);
}

fn moving_forced_balance(walls: bool) {
    let mut g = grid(walls);
    let shape = BodyFrameEllipsoid3::sphere(0.34).unwrap();
    let mut samplers = [
        BodyFrameSampler3::new(&g, 6).unwrap(),
        BodyFrameSampler3::new(&g, 6).unwrap(),
    ];
    let mut bodies = [moving(0), moving(1)];
    let external = [[0.00002, -0.00001, 0.00001, 0.000002, 0.000003, -0.000001]; 2];
    let imposed = imposed();
    let mut max_defect = 0.0_f64;
    let mut shared = false;
    for _ in 0..24 {
        let before = snapshot(&g);
        let old_bodies = bodies;
        let old = bodies.map(momentum);
        let receipt = g
            .step_with_moving_partial_saturation_bodies(
                &mut samplers,
                [&shape, &shape],
                &mut bodies,
                external,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        shared |= samplers[0]
            .fractions()
            .iter()
            .zip(samplers[1].fractions())
            .any(|(&a, &b)| a > 0.0 && b > 0.0);
        let (mass, mut p, mut h) = fluid_delta(&before, &g);
        close(mass, 0.0, 1e-12);
        for a in 0..2 {
            let after = momentum(bodies[a]);
            let midpoint_arm = core::array::from_fn(|i| {
                old_bodies[a].center[i] + 0.5 * old_bodies[a].velocity[i] - ORIGIN[i]
            });
            let external_orbital = cross(midpoint_arm, external[a][..3].try_into().unwrap());
            for i in 0..3 {
                p[i] += after.0[i] - old[a].0[i] - external[a][i];
                h[i] += after.1[i] - old[a].1[i] - external[a][i + 3] - external_orbital[i];
                close(
                    receipt.bodies[a].external_angular_impulse[i],
                    external[a][i + 3] + external_orbital[i],
                    2e-18,
                );
            }
            assert_eq!(receipt.bodies[a].coupling.forcing, Default::default());
            assert_eq!(receipt.bodies[a].coupling.periodic, Default::default());
            close(receipt.bodies[a].work_residual, 0.0, 2e-16);
        }
        for i in 0..3 {
            let dp = p[i] + receipt.coupling.walls.wall_impulse[i] - imposed.0[i];
            let dh = h[i] + receipt.coupling.walls.wall_angular_impulse[i]
                - imposed.1[i]
                - receipt.coupling.periodic.angular_momentum_change[i];
            max_defect = max_defect.max(dp.abs()).max(dh.abs());
            close(dp, 0.0, 2e-12);
            close(dh, 0.0, 9e-12);
            close(receipt.coupling.forcing.impulse[i], imposed.0[i], 4e-17);
            close(
                receipt.coupling.forcing.angular_impulse[i],
                imposed.1[i],
                3e-16,
            );
        }
    }
    assert!(shared);
    eprintln!("24 forced moving-body ticks, maximum world balance defect {max_defect:e}");
}

#[test]
fn forced_single_and_singleton_assembly_retain_identical_results() {
    let (mut single, mut assembly) = (grid(false), grid(false));
    let fractions = fields();
    let (mut body, mut bodies) = (body(0), [body(0)]);
    for _ in 0..12 {
        let external = [0.0003, -0.0002, 0.0001, -0.0001, 0.0002, 0.0];
        let a = single
            .step_with_partial_saturation_body(
                &fractions[0],
                &mut body,
                external,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        let b = assembly
            .step_with_partial_saturation_bodies(
                [fractions[0].as_slice()],
                &mut bodies,
                [external],
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(body, bodies[0]);
        assert_eq!(snapshot(&single), snapshot(&assembly));
        assert_eq!(a.hydrodynamic, b.bodies[0].hydrodynamic);
        assert_eq!(a.forcing, b.forcing);
        assert_eq!(a.walls, b.walls);
        assert_eq!(a.periodic, b.periodic);
        assert!(a.periodic.face_crossings > 0);
        assert_eq!(a.momentum_residual, b.bodies[0].momentum_residual);
    }
}

#[test]
fn late_force_overflow_and_cancellation_publish_no_fluid_or_body_state() {
    let (mut refused, mut clean) = (grid(true), grid(true));
    let f = fields();
    let mut bodies = [body(0), body(1)];
    let mut expected = bodies;
    let before = snapshot(&refused);
    refused.set_force_field(|c| {
        if c == [7, 7, 7] {
            [1e307, 0.0, 0.0]
        } else {
            field(c)
        }
    });
    assert!(matches!(
        refused.step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Collision { .. })
    ));
    assert_eq!(snapshot(&refused), before);
    assert_eq!(bodies, expected);
    refused.set_force_field(field);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        refused.step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            ORIGIN,
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    ));
    assert_eq!(snapshot(&refused), before);
    assert_eq!(bodies, expected);
    let a = refused
        .step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let b = clean
        .step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut expected,
            [[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(bodies, expected);
    assert_eq!(snapshot(&refused), snapshot(&clean));
}
