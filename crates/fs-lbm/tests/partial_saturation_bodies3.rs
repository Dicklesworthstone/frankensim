//! Native many-body impulse, work, shared-cell and trajectory regressions.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, BoundarySpec3, E3,
    FaceBoundary3, MovingPartialSaturationBody3, PartialSaturationBody3, PartialSaturationError3,
    Q3, equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.16e} != {b:.16e}, defect={:.4e}",
        a - b
    );
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
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn index(c: [usize; 3], n: usize) -> usize {
    (c[2] * n + c[1]) * n + c[0]
}
fn grid(n: usize, walls: bool) -> BoundaryGrid3 {
    BoundaryGrid3::new(
        n,
        n,
        n,
        0.83,
        [0.0; 3],
        if walls {
            BoundarySpec3::new([FaceBoundary3::Wall { velocity: [0.0; 3] }; 6])
        } else {
            BoundarySpec3::periodic()
        },
    )
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let d = g.dimensions();
    let mut f = Vec::new();
    for z in 0..d[2] {
        for y in 0..d[1] {
            for x in 0..d[0] {
                f.push(g.populations(x, y, z));
            }
        }
    }
    f
}
fn fluid_delta(
    before: &[[f64; Q3]],
    g: &BoundaryGrid3,
    origin: [f64; 3],
) -> (f64, [f64; 3], [f64; 3]) {
    let n = g.dimensions()[0];
    let (mut mass, mut p, mut h) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let f = g.populations(x, y, z);
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    let d = f[q] - before[index([x, y, z], n)][q];
                    mass += d;
                    for (a, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[a] += f64::from(e) * d;
                    }
                }
                let arm = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
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
fn body(a: usize) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [3.2 + 0.6 * a as f64, 3.7 - 0.3 * a as f64, 3.4],
        mass: 0.7 + 0.5 * a as f64,
        inertia_world: [[1.3, 0.12, -0.08], [0.12, 1.7, 0.06], [-0.08, 0.06, 2.0]],
        velocity: [-0.014 + 0.01 * a as f64, 0.008, 0.006],
        angular_velocity: [0.004, -0.003, 0.005 - 0.002 * a as f64],
    }
}
fn body_delta(
    before: PartialSaturationBody3,
    after: PartialSaturationBody3,
    origin: [f64; 3],
) -> ([f64; 3], [f64; 3], f64) {
    let p = core::array::from_fn(|i| before.mass * (after.velocity[i] - before.velocity[i]));
    let dw = core::array::from_fn(|i| after.angular_velocity[i] - before.angular_velocity[i]);
    let spin = before.inertia_world.map(|row| dot(row, dw));
    let orbital = cross(core::array::from_fn(|i| before.center[i] - origin[i]), p);
    let h = core::array::from_fn(|i| spin[i] + orbital[i]);
    let v = core::array::from_fn(|i| 0.5 * (after.velocity[i] + before.velocity[i]));
    let w =
        core::array::from_fn(|i| 0.5 * (after.angular_velocity[i] + before.angular_velocity[i]));
    (p, h, dot(p, v) + dot(spin, w))
}
fn fields() -> [Vec<f64>; 2] {
    let mut f = core::array::from_fn(|_| vec![0.0; 512]);
    for (c, a, b) in [
        ([4, 3, 3], 0.3, 0.6),
        ([3, 4, 3], 0.5, 0.2),
        ([2, 3, 4], 0.4, 0.0),
        ([4, 2, 2], 0.0, 0.7),
    ] {
        f[0][index(c, 8)] = a;
        f[1][index(c, 8)] = b;
    }
    f
}

#[test]
fn shared_cell_matches_independent_midpoint_blocks_and_population_formula() {
    let mut g = grid(8, false);
    let rho = 1.07;
    let u = [0.024, -0.017, 0.011];
    let initial = equilibrium3(rho, u);
    g.set_populations(3, 3, 3, &initial);
    let mut f = core::array::from_fn::<_, 2, _>(|_| vec![0.0; 512]);
    f[0][index([3, 3, 3], 8)] = 0.3;
    f[1][index([3, 3, 3], 8)] = 0.6;
    let mut bodies = [body(0), body(1)];
    for b in &mut bodies {
        b.center = [3.5; 3];
        b.angular_velocity = [0.0; 3];
    }
    let before = bodies;
    let impulses = [
        [0.0003, -0.0002, 0.0001, 0.0, 0.0, 0.0],
        [-0.0001, 0.0004, 0.0, 0.0, 0.0, 0.0],
    ];
    let report = g
        .step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut bodies,
            impulses,
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    let b = [0.3, 0.6].map(|e| e * (0.83 - 0.5) / (1.0 - 0.9 + 0.83 - 0.5));
    let mut expected = initial;
    for a in 0..2 {
        let weight = b[a] * rho;
        let after: [f64; 3] = core::array::from_fn(|i| {
            before[a].velocity[i]
                + (weight * (u[i] - before[a].velocity[i]) + impulses[a][i])
                    / (before[a].mass + 0.5 * weight)
        });
        let mid = core::array::from_fn(|i| 0.5 * (before[a].velocity[i] + after[i]));
        let solid = equilibrium3(rho, mid);
        for i in 0..3 {
            close(bodies[a].velocity[i], after[i], 2e-16);
        }
        close(
            report.bodies[a].hydrodynamic.coupling_weight_sum,
            b[a],
            2e-16,
        );
        for q in 0..Q3 {
            expected[q] += b[a] * (solid[q] - initial[q]);
        }
    }
    for q in 0..Q3 {
        let e = E3[q];
        let dest = [(3 + e.0) as usize, (3 + e.1) as usize, (3 + e.2) as usize];
        close(
            g.populations(dest[0], dest[1], dest[2])[q],
            expected[q],
            4e-16,
        );
    }
}

#[test]
fn simultaneous_full_tensor_bodies_close_world_momentum_and_work() {
    let origin = [-4.2, 7.1, 2.3];
    let mut g = grid(8, true);
    let f = fields();
    let mut bodies = [body(0), body(1)];
    let impulse = [
        [0.0005, -0.0003, 0.0002, -0.0004, 0.0006, 0.0003],
        [-0.0002, 0.0001, 0.0004, 0.0002, -0.0001, 0.0002],
    ];
    let mut max_defect = 0.0_f64;
    let mut wall_total = 0.0;
    for step in 0..48 {
        if step == 0 {
            g.set_populations(4, 3, 3, &equilibrium3(1.07, [0.024, -0.017, 0.011]));
        }
        let population = snapshot(&g);
        let before = bodies;
        let receipt = g
            .step_with_partial_saturation_bodies(
                f.each_ref().map(Vec::as_slice),
                &mut bodies,
                impulse,
                origin,
                &CancelGate::new(),
            )
            .unwrap();
        let (mass, mut p, mut h) = fluid_delta(&population, &g, origin);
        close(mass, 0.0, 8e-13);
        for a in 0..2 {
            let (dp, dh, energy) = body_delta(before[a], bodies[a], origin);
            let arm = core::array::from_fn(|i| before[a].center[i] - origin[i]);
            let orbital = cross(arm, impulse[a][..3].try_into().unwrap());
            for i in 0..3 {
                p[i] += dp[i] - impulse[a][i];
                h[i] += dh[i] - impulse[a][i + 3] - orbital[i];
            }
            let response = receipt.bodies[a];
            close(energy, response.kinetic_energy_change, 2e-17);
            close(
                energy,
                response.external_work + response.hydrodynamic.body_work,
                4e-16,
            );
            for d in response.momentum_residual {
                close(d, 0.0, 5e-14);
            }
        }
        for i in 0..3 {
            let dp = p[i] + receipt.walls.wall_impulse[i];
            let dh = h[i] + receipt.walls.wall_angular_impulse[i];
            max_defect = max_defect.max(dp.abs()).max(dh.abs());
            close(dp, 0.0, 1e-12);
            close(dh, 0.0, 8e-12);
            wall_total += receipt.walls.wall_impulse[i].abs();
        }
        close(receipt.walls.wall_work, 0.0, 0.0);
    }
    assert!(wall_total > 1e-4);
    eprintln!("48-tick multiple-body maximum world balance defect: {max_defect:e}");
}

#[test]
fn body_order_permutation_and_zero_fraction_members_preserve_one_tick() {
    let f = fields();
    let mut a = grid(8, false);
    let mut b = grid(8, false);
    let mut ab = [body(0), body(1)];
    let mut ba = [ab[1], ab[0]];
    let ra = a
        .step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut ab,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    let rb = b
        .step_with_partial_saturation_bodies(
            [&f[1], &f[0]],
            &mut ba,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(ab, [ba[1], ba[0]]);
    assert_eq!(ra.bodies, [rb.bodies[1], rb.bodies[0]]);
    for (fa, fb) in snapshot(&a).iter().zip(snapshot(&b)) {
        for q in 0..Q3 {
            close(fa[q], fb[q], 2e-16);
        }
    }
    let zeros = vec![0.0; 512];
    let mut ordinary = grid(8, false);
    let mut many = grid(8, false);
    ordinary.step();
    let mut body = [body(0), body(1)];
    let r = many
        .step_with_partial_saturation_bodies(
            [&zeros, &zeros],
            &mut body,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(snapshot(&ordinary), snapshot(&many));
    assert_eq!(body, [self::body(0), self::body(1)]);
    assert_eq!(r.population_impulse_residual, [0.0; 3]);
}

fn moving(center: [f64; 3], velocity: [f64; 3]) -> MovingPartialSaturationBody3 {
    let inertia = InertiaTensor3::new(Mat3 {
        m: [1.3, 0.12, -0.08, 0.12, 1.7, 0.06, -0.08, 0.06, 2.0],
    })
    .unwrap();
    MovingPartialSaturationBody3 {
        center,
        velocity,
        inertia,
        mass: 0.7,
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.15))).unwrap(),
        angular_momentum_body: inertia.momentum(Vec3::new(0.006, -0.004, 0.008)).unwrap(),
    }
}
fn moving_momenta(body: MovingPartialSaturationBody3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let p = body.velocity.map(|v| body.mass * v);
    let spin = array(body.orientation.rotate(body.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| body.center[i] - origin[i]), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}

#[test]
fn moving_assembly_preserves_actual_spatial_momentum_and_replays() {
    let mut a = grid(8, true);
    let mut b = grid(8, true);
    let mut aa = [
        moving([2.8, 3.7, 4.0], [0.013, -0.007, 0.006]),
        moving([4.4, 3.9, 4.0], [-0.008, 0.009, -0.003]),
    ];
    let mut bb = aa;
    let initial = aa;
    let shapes = [
        BodyFrameEllipsoid3::new(Vec3::new(0.68, 0.47, 0.51)).unwrap(),
        BodyFrameEllipsoid3::sphere(0.62).unwrap(),
    ];
    let mut sa = core::array::from_fn(|_| BodyFrameSampler3::new(&a, 4).unwrap());
    let mut sb = core::array::from_fn(|_| BodyFrameSampler3::new(&b, 4).unwrap());
    let impulse = [
        [0.00002, -0.00001, 0.00001, -0.00001, 0.00002, 0.00001],
        [-0.00001, 0.00002, 0.00001, 0.00001, -0.00002, 0.00001],
    ];
    let origin = [-3.1, 8.0, 1.7];
    let mut maximum = 0.0_f64;
    let mut shared = false;
    for _ in 0..24 {
        let population = snapshot(&a);
        let before = aa;
        let ra = a
            .step_with_moving_partial_saturation_bodies(
                &mut sa,
                [&shapes[0], &shapes[1]],
                &mut aa,
                impulse,
                origin,
                &CancelGate::new(),
            )
            .unwrap();
        let rb = b
            .step_with_moving_partial_saturation_bodies(
                &mut sb,
                [&shapes[0], &shapes[1]],
                &mut bb,
                impulse,
                origin,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(aa, bb);
        assert_eq!(ra, rb);
        assert_eq!(snapshot(&a), snapshot(&b));
        shared |= sa[0]
            .fractions()
            .iter()
            .zip(sa[1].fractions())
            .any(|(a, b)| *a > 0.0 && *b > 0.0);
        let (mass, mut p, mut h) = fluid_delta(&population, &a, origin);
        close(mass, 0.0, 8e-13);
        for i in 0..2 {
            let (p0, h0) = moving_momenta(before[i], origin);
            let (p1, h1) = moving_momenta(aa[i], origin);
            for k in 0..3 {
                p[k] += p1[k] - p0[k] - impulse[i][k];
                h[k] += h1[k] - h0[k] - ra.bodies[i].external_angular_impulse[k];
            }
            for r in ra.bodies[i].momentum_residual {
                close(r, 0.0, 8e-14);
            }
            close(ra.bodies[i].work_residual, 0.0, 5e-16);
            assert_eq!(ra.bodies[i].coupling.walls.measured_links, 0);
        }
        for i in 0..3 {
            maximum = maximum
                .max((p[i] + ra.coupling.walls.wall_impulse[i]).abs())
                .max((h[i] + ra.coupling.walls.wall_angular_impulse[i]).abs());
            close(p[i] + ra.coupling.walls.wall_impulse[i], 0.0, 1e-12);
            close(h[i] + ra.coupling.walls.wall_angular_impulse[i], 0.0, 8e-12);
        }
    }
    assert!(shared);
    assert_ne!(aa[0].center, initial[0].center);
    assert_ne!(aa[1].orientation, initial[1].orientation);
    eprintln!("24-tick moving assembly maximum world balance defect: {maximum:e}");
}

#[test]
fn single_member_preserves_original_body_and_populations_bit_for_bit() {
    let mut a = grid(8, false);
    let mut b = grid(8, false);
    let mut single = moving([3.7, 3.9, 4.1], [0.014, -0.008, 0.006]);
    let mut many = [single];
    let shape = BodyFrameEllipsoid3::new(Vec3::new(1.13, 0.72, 0.51)).unwrap();
    let mut sa = BodyFrameSampler3::new(&a, 3).unwrap();
    let mut sb = [BodyFrameSampler3::new(&b, 3).unwrap()];
    for _ in 0..16 {
        let ra = a
            .step_with_moving_partial_saturation_body(
                &mut sa,
                &shape,
                &mut single,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        let rb = b
            .step_with_moving_partial_saturation_bodies(
                &mut sb,
                [&shape],
                &mut many,
                [[0.0; 6]],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(single, many[0]);
        assert_eq!(snapshot(&a), snapshot(&b));
        assert_eq!(ra, rb.bodies[0]);
    }
}

#[test]
fn fraction_overfill_and_late_interbody_collision_roll_back_the_entire_assembly() {
    let mut g = grid(8, false);
    let old = snapshot(&g);
    let mut f = fields();
    f[1][index([4, 3, 3], 8)] = 0.8;
    let mut bodies = [body(0), body(1)];
    let before = bodies;
    assert!(
        g.step_with_partial_saturation_bodies(
            f.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, before);
    assert_eq!(snapshot(&g), old);
    let mut pair = [
        moving([2.9, 4.0, 4.0], [0.0; 3]),
        moving([4.2, 4.0, 4.0], [0.0; 3]),
    ];
    for b in &mut pair {
        b.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    }
    let initial = pair;
    let shape = BodyFrameEllipsoid3::sphere(0.6).unwrap();
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 4).unwrap());
    let force = [
        [0.11, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.11, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let error = g
        .step_with_moving_partial_saturation_bodies(
            &mut samplers,
            [&shape, &shape],
            &mut pair,
            force,
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap_err();
    assert_eq!(
        error,
        PartialSaturationError3::Invalid {
            reason: "moving assembly swept enclosing spheres meet"
        }
    );
    assert_eq!(pair, initial);
    assert_eq!(snapshot(&g), old);
    // Successful retry sees the original fluid, body and topology admission.
    g.step_with_moving_partial_saturation_bodies(
        &mut samplers,
        [&shape, &shape],
        &mut pair,
        [[0.0; 6]; 2],
        [0.0; 3],
        &CancelGate::new(),
    )
    .unwrap();
}

#[test]
fn cancellation_during_last_body_finalization_publishes_nothing() {
    use std::cell::Cell;
    struct Cancelling<'a> {
        calls: Cell<usize>,
        gate: &'a CancelGate,
        sphere: BodyFrameEllipsoid3,
    }
    impl BodyFrameShape3 for Cancelling<'_> {
        fn bounding_radius(&self) -> f64 {
            let n = self.calls.get() + 1;
            self.calls.set(n);
            if n == 4 {
                self.gate.request();
            }
            self.sphere.bounding_radius()
        }
        fn contains(&self, p: Vec3) -> bool {
            self.sphere.contains(p)
        }
    }
    let mut g = grid(8, false);
    let old = snapshot(&g);
    let gate = CancelGate::new();
    let sphere = BodyFrameEllipsoid3::sphere(0.6).unwrap();
    let last = Cancelling {
        calls: Cell::new(0),
        gate: &gate,
        sphere,
    };
    let mut pair = [
        moving([2.8, 4.0, 4.0], [0.01, 0.0, 0.0]),
        moving([4.4, 4.0, 4.0], [0.01, 0.0, 0.0]),
    ];
    let initial = pair;
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 3).unwrap());
    let error = g
        .step_with_moving_partial_saturation_bodies(
            &mut samplers,
            [&sphere, &last],
            &mut pair,
            [[0.0; 6]; 2],
            [0.0; 3],
            &gate,
        )
        .unwrap_err();
    assert_eq!(error, PartialSaturationError3::Cancelled);
    assert_eq!(last.calls.get(), 4);
    assert_eq!(pair, initial);
    assert_eq!(snapshot(&g), old);
}
