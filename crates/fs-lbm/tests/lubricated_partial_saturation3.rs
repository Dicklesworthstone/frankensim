//! Normal-film asymptotics and simultaneous native fluid/body momentum.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3,
    MovingPartialSaturationBody3, PartialSaturationBody3, PartialSaturationError3, Q3,
    SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}, delta={:.3e}",
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
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
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
fn body(x: f64, mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [x, 4.0, 4.0],
        mass,
        inertia_world: [[0.6, 0.04, -0.02], [0.04, 0.8, 0.03], [-0.02, 0.03, 0.9]],
        velocity,
        angular_velocity: [0.002, -0.001, 0.003],
    }
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                result.push(grid.populations(x, y, z));
            }
        }
    }
    result
}
fn resistance(mu: f64, ra: f64, rb: f64, gap: f64, cutoff: f64, roughness: f64) -> f64 {
    let r = ra * rb / (ra + rb);
    6.0 * core::f64::consts::PI * mu * r * r * (1.0 / gap.max(roughness) - 1.0 / cutoff).max(0.0)
}

#[test]
fn unequal_spheres_match_normal_stokes_impulse_and_physical_work() {
    // The zero-occupancy operator limit isolates the missing-film law from
    // resolved fluid drag. Actual sphere sampling is exercised below.
    let fields = [vec![0.0; 512], vec![0.0; 512]];
    let radii = [0.4, 0.7];
    let initial = [
        body(3.0, 2.0, [0.023, 0.004, -0.003]),
        body(4.3, 3.0, [-0.012, 0.004, -0.003]),
    ];
    let zeta = resistance(0.11, radii[0], radii[1], 0.2, 0.8, 0.01);
    let lambda = -zeta * (initial[0].velocity[0] - initial[1].velocity[0])
        / (1.0 + 0.5 * zeta * (1.0 / initial[0].mass + 1.0 / initial[1].mass));
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.31, -0.27, 0.19))).unwrap();
    for rotated in [false, true] {
        let mut before = initial;
        let n = if rotated {
            for b in &mut before {
                b.center = array(
                    rotation
                        .rotate(Vec3::new(b.center[0], b.center[1], b.center[2]))
                        .unwrap(),
                );
                b.velocity = array(
                    rotation
                        .rotate(Vec3::new(b.velocity[0], b.velocity[1], b.velocity[2]))
                        .unwrap(),
                );
                b.angular_velocity = array(
                    rotation
                        .rotate(Vec3::new(
                            b.angular_velocity[0],
                            b.angular_velocity[1],
                            b.angular_velocity[2],
                        ))
                        .unwrap(),
                );
                let tensor = Mat3 {
                    m: core::array::from_fn(|i| b.inertia_world[i / 3][i % 3]),
                };
                let r = rotation.matrix();
                let tensor = r.compose(tensor).compose(r.transpose());
                b.inertia_world = core::array::from_fn(|i| {
                    core::array::from_fn(|j| tensor.m[3 * i.max(j) + i.min(j)])
                });
            }
            array(rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap())
        } else {
            [1.0, 0.0, 0.0]
        };
        let mut state = before;
        let mut g = grid(false);
        let receipt = g
            .step_with_lubricated_partial_saturation_bodies(
                fields.each_ref().map(Vec::as_slice),
                radii,
                &mut state,
                [[0.0; 6]; 2],
                SphereLubrication3::new(0.11, 0.8, 0.01, 1).unwrap(),
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(receipt.lubrication.active_pairs, 1);
        let mut work = 0.0;
        for a in 0..2 {
            let sign = if a == 0 { 1.0 } else { -1.0 };
            for i in 0..3 {
                close(
                    state[a].velocity[i],
                    before[a].velocity[i] + sign * lambda * n[i] / before[a].mass,
                    3e-17,
                );
                close(
                    receipt.lubrication.body_impulses[a][i],
                    sign * lambda * n[i],
                    3e-17,
                );
                close(
                    state[a].angular_velocity[i],
                    before[a].angular_velocity[i],
                    3e-17,
                );
                close(receipt.bodies[a].momentum_residual[i], 0.0, 3e-17);
            }
            let dp = core::array::from_fn(|i| {
                before[a].mass * (state[a].velocity[i] - before[a].velocity[i])
            });
            let midpoint =
                core::array::from_fn(|i| 0.5 * (state[a].velocity[i] + before[a].velocity[i]));
            work += dot(dp, midpoint);
            close(receipt.lubrication.body_work[a], dot(dp, midpoint), 2e-18);
        }
        for i in 0..3 {
            close(
                receipt.lubrication.body_impulses[0][i] + receipt.lubrication.body_impulses[1][i],
                0.0,
                0.0,
            );
            close(
                receipt.lubrication.body_angular_impulses[0][i]
                    + receipt.lubrication.body_angular_impulses[1][i],
                0.0,
                4e-17,
            );
        }
        let slip = dot(
            n,
            core::array::from_fn(|i| {
                0.5 * ((before[0].velocity[i] + state[0].velocity[i])
                    - (before[1].velocity[i] + state[1].velocity[i]))
            }),
        );
        close(
            receipt.lubrication.dissipated_energy,
            zeta * slip * slip,
            3e-18,
        );
        close(work + receipt.lubrication.dissipated_energy, 0.0, 3e-18);
        close(receipt.lubrication.work_residual, 0.0, 3e-18);
    }
}

#[test]
fn shared_body_cluster_matches_both_analytic_modes_and_permutations() {
    let fields = [vec![0.0; 512], vec![0.0; 512], vec![0.0; 512]];
    let initial = [
        body(2.8, 1.7, [0.025, 0.003, 0.0]),
        body(4.0, 1.7, [-0.02, 0.003, 0.0]),
        body(5.2, 1.7, [-0.005, 0.003, 0.0]),
    ];
    // [0.015,0,-0.015] and [0.01,-0.02,0.01] are exact chain
    // Laplacian eigenmodes, with eigenvalues zeta and 3*zeta.
    let zeta = resistance(0.11, 0.5, 0.5, 0.2, 0.6, 0.01);
    let gain = |e: f64| (1.7 - 0.5 * e) / (1.7 + 0.5 * e);
    let expected = [
        0.015 * gain(zeta) + 0.01 * gain(3.0 * zeta),
        -0.02 * gain(3.0 * zeta),
        -0.015 * gain(zeta) + 0.01 * gain(3.0 * zeta),
    ];
    for order in [[0, 1, 2], [2, 0, 1]] {
        for boost in [0.0, 0.017] {
            let mut states = order.map(|i| {
                let mut b = initial[i];
                b.velocity[0] += boost;
                b
            });
            let before = states;
            let receipt = grid(false)
                .step_with_lubricated_partial_saturation_bodies(
                    fields.each_ref().map(Vec::as_slice),
                    [0.5; 3],
                    &mut states,
                    [[0.0; 6]; 3],
                    SphereLubrication3::new(0.11, 0.6, 0.01, 2).unwrap(),
                    ORIGIN,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(receipt.lubrication.active_pairs, 2);
            let mut actual_work = 0.0;
            for a in 0..3 {
                close(states[a].velocity[0], expected[order[a]] + boost, 4e-17);
                actual_work += 0.5
                    * before[a].mass
                    * (states[a].velocity[0] - before[a].velocity[0])
                    * (states[a].velocity[0] + before[a].velocity[0]);
            }
            close(
                actual_work + receipt.lubrication.dissipated_energy,
                0.0,
                4e-18,
            );
            close(receipt.lubrication.constitutive_residual, 0.0, 2e-17);
        }
    }
}

#[test]
fn roughness_saturation_cutoff_and_singleton_retain_declared_limits() {
    let fields = [vec![0.0; 512], vec![0.0; 512]];
    let mut outcomes = Vec::new();
    for gap in [0.08, 0.04] {
        let mut bodies = [
            body(3.0, 1.0, [0.01, 0.0, 0.0]),
            body(4.0 + gap, 1.0, [-0.01, 0.0, 0.0]),
        ];
        let receipt = grid(false)
            .step_with_lubricated_partial_saturation_bodies(
                fields.each_ref().map(Vec::as_slice),
                [0.5; 2],
                &mut bodies,
                [[0.0; 6]; 2],
                SphereLubrication3::new(0.11, 0.5, 0.1, 1).unwrap(),
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        outcomes.push((
            bodies.map(|b| b.velocity),
            receipt.lubrication.dissipated_energy,
        ));
    }
    assert_eq!(outcomes[0], outcomes[1]);
    let initial = [
        body(2.6, 1.0, [0.01, 0.0, 0.0]),
        body(5.0, 1.0, [-0.01, 0.0, 0.0]),
    ];
    let (mut a, mut b) = (initial, initial);
    let (mut ga, mut gb) = (grid(false), grid(false));
    let ra = ga
        .step_with_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            &mut a,
            [[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let rb = gb
        .step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 2],
            &mut b,
            [[0.0; 6]; 2],
            SphereLubrication3::new(0.11, 0.5, 0.1, 1).unwrap(),
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(ra, rb);
    assert_eq!(snapshot(&ga), snapshot(&gb));
    let (mut a, mut b) = ([initial[0]], [initial[0]]);
    let (mut ga, mut gb) = (grid(false), grid(false));
    let ra = ga
        .step_with_partial_saturation_bodies(
            [fields[0].as_slice()],
            &mut a,
            [[0.0; 6]],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let rb = gb
        .step_with_lubricated_partial_saturation_bodies(
            [fields[0].as_slice()],
            [0.5],
            &mut b,
            [[0.0; 6]],
            SphereLubrication3::new(0.11, 0.5, 0.1, 1).unwrap(),
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(ra, rb);
    assert_eq!(snapshot(&ga), snapshot(&gb));
    assert!(SphereLubrication3::new(0.11, 0.1, 0.2, 1).is_err());
    assert!(SphereLubrication3::new(0.0, 0.1, 0.01, 1).is_err());
    assert!(SphereLubrication3::new(0.11, 0.1, 0.01, 65).is_err());
}

fn moving(a: usize) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [2.8 + 1.2 * a as f64, 4.0, 4.0],
        orientation: So3::identity(),
        mass: 0.7 + 0.2 * a as f64,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
            ],
        })
        .unwrap(),
        velocity: [0.003 - 0.003 * a as f64, 0.001, -0.001],
        angular_momentum_body: Vec3::new(0.00008, -0.00006, 0.00004),
    }
}
fn momentum(b: MovingPartialSaturationBody3) -> ([f64; 3], [f64; 3]) {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| b.center[i] - ORIGIN[i]), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}
fn kinetic(b: MovingPartialSaturationBody3) -> f64 {
    let omega = b.inertia.velocity(b.angular_momentum_body).unwrap();
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5 * dot(array(omega), array(b.angular_momentum_body))
}
fn fluid_delta(before: &[[f64; Q3]], g: &BoundaryGrid3) -> ([f64; 3], [f64; 3]) {
    let (mut p, mut h) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
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

#[test]
fn moving_full_tensor_spheres_keep_forced_world_momentum_and_viscous_work() {
    let mut g = grid(true);
    let force = |c: [usize; 3]| [1e-6 * (1 + c[1]) as f64, -2e-6, 0.5e-6 * (1 + c[0]) as f64];
    g.set_force_field(force);
    let (mut imposed_p, mut imposed_h) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let f = force([x, y, z]);
                let h = cross(
                    [
                        x as f64 + 0.5 - ORIGIN[0],
                        y as f64 + 0.5 - ORIGIN[1],
                        z as f64 + 0.5 - ORIGIN[2],
                    ],
                    f,
                );
                for i in 0..3 {
                    imposed_p[i] += f[i];
                    imposed_h[i] += h[i];
                }
            }
        }
    }
    let mut bodies = core::array::from_fn::<_, 3, _>(moving);
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 4).unwrap());
    let external = [[0.00001, -0.00002, 0.00001, 0.000001, -0.000002, 0.000001]; 3];
    let model = SphereLubrication3::new(0.11, 0.6, 0.02, 3).unwrap();
    let mut max_defect = 0.0_f64;
    let mut loss = 0.0;
    for _ in 0..24 {
        let before = snapshot(&g);
        let original = bodies;
        let old = bodies.map(momentum);
        let report = g
            .step_with_lubricated_spheres(
                &mut samplers,
                [0.5; 3],
                &mut bodies,
                external,
                model,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(report.coupling.lubrication.active_pairs, 2);
        let (mut dp, mut dh) = fluid_delta(&before, &g);
        let mut body_work = 0.0;
        for a in 0..3 {
            let now = momentum(bodies[a]);
            let arm = core::array::from_fn(|i| {
                original[a].center[i] + 0.5 * original[a].velocity[i] - ORIGIN[i]
            });
            let moment = cross(arm, external[a][..3].try_into().unwrap());
            for i in 0..3 {
                dp[i] += now.0[i] - old[a].0[i] - external[a][i];
                dh[i] += now.1[i] - old[a].1[i] - external[a][i + 3] - moment[i];
                close(report.bodies[a].momentum_residual[i], 0.0, 2e-15);
                close(report.bodies[a].momentum_residual[i + 3], 0.0, 2e-15);
            }
            let midpoint =
                core::array::from_fn(|i| 0.5 * (original[a].velocity[i] + bodies[a].velocity[i]));
            let w = dot(report.coupling.lubrication.body_impulses[a], midpoint);
            close(report.coupling.lubrication.body_work[a], w, 2e-18);
            body_work += w;
            let body = &report.bodies[a];
            close(
                kinetic(bodies[a]) - kinetic(original[a]),
                body.coupling.external_work
                    + body.coupling.hydrodynamic.body_work
                    + w
                    + body.free_rotation_energy_change,
                3e-17,
            );
            close(body.work_residual, 0.0, 3e-17);
        }
        let lub = report.coupling.lubrication;
        let mut independent_impulses = [[0.0; 3]; 3];
        let mut independent_loss = 0.0;
        for a in 0..2 {
            let delta = core::array::from_fn(|i| {
                report.bodies[a + 1].midpoint_center[i] - report.bodies[a].midpoint_center[i]
            });
            let d = dot(delta, delta).sqrt();
            let n = delta.map(|x| x / d);
            let slip = dot(
                n,
                core::array::from_fn(|i| {
                    0.5 * ((original[a].velocity[i] + bodies[a].velocity[i])
                        - (original[a + 1].velocity[i] + bodies[a + 1].velocity[i]))
                }),
            );
            let zeta = resistance(0.11, 0.5, 0.5, d - 1.0, 0.6, 0.02);
            let lambda = -zeta * slip;
            independent_loss += zeta * slip * slip;
            for i in 0..3 {
                independent_impulses[a][i] += lambda * n[i];
                independent_impulses[a + 1][i] -= lambda * n[i];
            }
        }
        for a in 0..3 {
            for i in 0..3 {
                close(lub.body_impulses[a][i], independent_impulses[a][i], 2e-17);
            }
        }
        close(lub.dissipated_energy, independent_loss, 2e-18);
        close(body_work + lub.dissipated_energy, 0.0, 2e-18);
        loss += lub.dissipated_energy;
        for i in 0..3 {
            let p = dp[i] + report.coupling.walls.wall_impulse[i] - imposed_p[i];
            let h = dh[i] + report.coupling.walls.wall_angular_impulse[i] - imposed_h[i];
            close(p, 0.0, 3e-12);
            close(h, 0.0, 9e-12);
            max_defect = max_defect.max(p.abs()).max(h.abs());
        }
    }
    assert!(loss > 0.0);
    eprintln!(
        "24 forced lubricated moving-body ticks: max world momentum defect={max_defect:e}, viscous loss={loss:e}"
    );
}

#[test]
fn stiff_midpoint_relaxation_is_passive_with_explicit_unresolved_reversal() {
    let fields = [vec![0.0; 512], vec![0.0; 512]];
    let mut states = [
        body(3.0, 1.0, [0.01, 0.0, 0.0]),
        body(4.001, 1.0, [-0.01, 0.0, 0.0]),
    ];
    let before = states;
    let receipt = grid(false)
        .step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 2],
            &mut states,
            [[0.0; 6]; 2],
            SphereLubrication3::new(0.11, 0.5, 0.001, 1).unwrap(),
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    // Do not mistake A-stability/passivity for resolved near-contact motion.
    assert!(states[0].velocity[0] < 0.0);
    let before_energy = before
        .iter()
        .map(|b| 0.5 * b.mass * dot(b.velocity, b.velocity))
        .sum::<f64>();
    let after_energy = states
        .iter()
        .map(|b| 0.5 * b.mass * dot(b.velocity, b.velocity))
        .sum::<f64>();
    assert!(after_energy < before_energy);
    close(
        after_energy - before_energy + receipt.lubrication.dissipated_energy,
        0.0,
        3e-18,
    );
}

#[test]
fn pair_limits_overlap_late_motion_failure_and_cancellation_publish_nothing() {
    let fields = [vec![0.0; 512], vec![0.0; 512], vec![0.0; 512]];
    let initial = [
        body(2.8, 1.0, [0.01, 0.0, 0.0]),
        body(4.0, 1.0, [0.0; 3]),
        body(5.2, 1.0, [-0.01, 0.0, 0.0]),
    ];
    let mut g = grid(false);
    let before = snapshot(&g);
    let mut bodies = initial;
    assert!(
        g.step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 3],
            &mut bodies,
            [[0.0; 6]; 3],
            SphereLubrication3::new(0.11, 0.6, 0.02, 1).unwrap(),
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&g), before);
    assert!(
        g.step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.7; 3],
            &mut bodies,
            [[0.0; 6]; 3],
            SphereLubrication3::new(0.11, 0.6, 0.02, 3).unwrap(),
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&g), before);
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        g.step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 3],
            &mut bodies,
            [[0.0; 6]; 3],
            SphereLubrication3::new(0.11, 0.6, 0.02, 3).unwrap(),
            ORIGIN,
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&g), before);
    let mut moving = [moving(0), moving(1)];
    moving[0].center[0] = 3.0;
    moving[1].center[0] = 4.04;
    for b in &mut moving {
        b.velocity = [0.0; 3];
    }
    let old = moving;
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 3).unwrap());
    let model = SphereLubrication3::new(0.000001, 0.6, 0.01, 1).unwrap();
    let impulses = [
        [0.08, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.08, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    assert!(
        g.step_with_lubricated_spheres(
            &mut samplers,
            [0.5; 2],
            &mut moving,
            impulses,
            model,
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(moving, old);
    assert_eq!(snapshot(&g), before);
    let mut fresh = grid(false);
    let mut expected = old;
    let mut fresh_samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&fresh, 3).unwrap());
    let a = g
        .step_with_lubricated_spheres(
            &mut samplers,
            [0.5; 2],
            &mut moving,
            [[0.0; 6]; 2],
            model,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let b = fresh
        .step_with_lubricated_spheres(
            &mut fresh_samplers,
            [0.5; 2],
            &mut expected,
            [[0.0; 6]; 2],
            model,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(moving, expected);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
