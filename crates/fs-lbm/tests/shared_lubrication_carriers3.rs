//! One actual finite inertia across several original sphere/plane film ports.
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, CarrierLubricationFace3,
    E3, FiniteLubricationPlane3, LubricationPlane3, LubricationPlanePatch3, PartialSaturationBody3,
    PartialSaturationError3, Q3, SharedLubricationSurface3, SphereLubrication3,
};

const RADIUS: f64 = 0.7;
const MU: f64 = 0.11;
const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
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
fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(8, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        velocity,
        inertia_world: [
            [0.12, 0.009, -0.004],
            [0.009, 0.17, 0.006],
            [-0.004, 0.006, 0.21],
        ]
        .map(|r| r.map(|x| mass * x)),
        angular_velocity: [0.0; 3],
    }
}
fn patch(point: [f64; 3], normal: [f64; 3]) -> LubricationPlanePatch3 {
    LubricationPlanePatch3::new(point, normal, 0.04, 0.001).unwrap()
}
fn model(films: usize) -> SphereLubrication3 {
    SphereLubrication3::new(MU, 0.1, 0.001, films)
        .unwrap()
        .with_resistance_mode_limit(5 * films)
        .unwrap()
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..512)
        .map(|i| g.populations(i % 8, (i / 8) % 8, i / 64))
        .collect()
}
fn momentum(b: PartialSaturationBody3) -> [f64; 6] {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = b.inertia_world.map(|r| dot(r, b.angular_velocity));
    let orbital = cross(core::array::from_fn(|i| b.center[i] - ORIGIN[i]), p);
    core::array::from_fn(|i| {
        if i < 3 {
            p[i]
        } else {
            spin[i - 3] + orbital[i - 3]
        }
    })
}
fn external(b: PartialSaturationBody3, j: [f64; 6]) -> [f64; 6] {
    let orbital = cross(
        core::array::from_fn(|i| b.center[i] - ORIGIN[i]),
        core::array::from_fn(|i| j[i]),
    );
    core::array::from_fn(|i| if i < 3 { j[i] } else { j[i] + orbital[i - 3] })
}
fn energy(b: PartialSaturationBody3) -> f64 {
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                b.angular_velocity,
                b.inertia_world.map(|r| dot(r, b.angular_velocity)),
            )
}
fn opposing() -> (
    [PartialSaturationBody3; 2],
    PartialSaturationBody3,
    [LubricationPlanePatch3; 2],
) {
    (
        [
            body([3.21, 4.0, 4.0], 2.0, [0.012, 0.0, 0.0]),
            body([4.79, 4.0, 4.0], 3.0, [-0.006, 0.0, 0.0]),
        ],
        body([4.0, 4.0, 4.0], 5.0, [0.003, 0.0, 0.0]),
        [
            patch([2.5, 4.0, 4.0], [1.0, 0.0, 0.0]),
            patch([5.5, 4.0, 4.0], [-1.0, 0.0, 0.0]),
        ],
    )
}
fn finite(p: LubricationPlanePatch3, b: PartialSaturationBody3) -> FiniteLubricationPlane3 {
    let mut result =
        FiniteLubricationPlane3::new(p.point(), p.normal(), b, p.cutoff_gap(), p.roughness_gap())
            .unwrap();
    if let Some((cutoff, roughness)) = p.tangential_gaps() {
        result = result.with_tangential(cutoff, roughness).unwrap();
    }
    result
}

fn attached(carrier: usize, patch: LubricationPlanePatch3) -> SharedLubricationSurface3 {
    let mut face = CarrierLubricationFace3::new(
        carrier,
        patch.point(),
        patch.normal(),
        patch.cutoff_gap(),
        patch.roughness_gap(),
    )
    .unwrap();
    if let Some((cutoff, roughness)) = patch.tangential_gaps() {
        face = face.with_tangential(cutoff, roughness).unwrap();
    }
    SharedLubricationSurface3::Carried(face)
}

#[test]
fn two_faces_share_the_three_mass_midpoint_response_and_one_authored_load() {
    let (initial, carrier, patches) = opposing();
    let mut spheres = initial;
    let mut carriers = [carrier];
    let mut g = grid();
    let zero = vec![0.0; 512];
    let surfaces = patches.map(|patch| attached(0, patch));
    let loads = [
        [0.002, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let load = [0.0015, 0.0, 0.0, 0.0, 0.0, 0.0];
    let report = g
        .step_with_lubricated_partial_saturation_shared_carriers(
            [&zero, &zero],
            [RADIUS; 2],
            &mut spheres,
            loads,
            model(2),
            &surfaces,
            &mut carriers,
            &[load],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let k: [f64; 2] = core::array::from_fn(|i| {
        let gap = dot(
            patches[i].normal(),
            core::array::from_fn(|j| initial[i].center[j] - patches[i].point()[j]),
        ) - RADIUS;
        6.0 * core::f64::consts::PI
            * MU
            * RADIUS
            * RADIUS
            * (1.0 / gap - 1.0 / patches[i].cutoff_gap())
    });
    // Independent 2-by-2 elimination of the THREE physical x velocities.
    // These are Cartesian sphere impulses, so the shared-carrier off-diagonal
    // term is positive even though the geometric plane normals are opposite.
    let diagonal: [f64; 2] =
        core::array::from_fn(|i| 1.0 / k[i] + 0.5 * (1.0 / initial[i].mass + 1.0 / carrier.mass));
    let off = 0.5 / carrier.mass;
    let rhs: [f64; 2] = core::array::from_fn(|i| {
        -(initial[i].velocity[0] - carrier.velocity[0]
            + 0.5 * (loads[i][0] / initial[i].mass - load[0] / carrier.mass))
    });
    let determinant = diagonal[0] * diagonal[1] - off * off;
    let impulse = [
        (rhs[0] * diagonal[1] - off * rhs[1]) / determinant,
        (diagonal[0] * rhs[1] - off * rhs[0]) / determinant,
    ];
    for i in 0..2 {
        close(
            report.coupling.lubrication.body_impulses[i][0],
            impulse[i],
            3e-17,
        );
        close(
            spheres[i].velocity[0],
            initial[i].velocity[0] + (loads[i][0] + impulse[i]) / initial[i].mass,
            2e-17,
        );
        assert_eq!(report.surfaces[i].carrier, Some(0));
        close(report.surfaces[i].impulse[0], -impulse[i], 3e-17);
    }
    close(
        carriers[0].velocity[0],
        carrier.velocity[0] + (load[0] - impulse[0] - impulse[1]) / carrier.mass,
        2e-17,
    );
    assert_eq!(report.carriers.len(), 1);
    let shared = report.carriers[0];
    assert_eq!(shared.motion.body_after, carriers[0]);
    close(shared.impulse[0], -impulse[0] - impulse[1], 3e-17);
    let wall_mid = 0.5 * (carrier.velocity[0] + carriers[0].velocity[0]);
    close(shared.motion.external_work, load[0] * wall_mid, 3e-20);
    close(shared.work, shared.impulse[0] * wall_mid, 3e-20);
    close(
        report.surfaces.iter().map(|s| s.work).sum(),
        shared.work,
        3e-20,
    );
    let loss: f64 = (0..2)
        .map(|i| {
            k[i] * (0.5 * (initial[i].velocity[0] + spheres[i].velocity[0]) - wall_mid).powi(2)
        })
        .sum();
    close(report.coupling.lubrication.dissipated_energy, loss, 3e-18);
    let delta = energy(spheres[0]) + energy(spheres[1]) + energy(carriers[0])
        - energy(initial[0])
        - energy(initial[1])
        - energy(carrier);
    let authored = report
        .coupling
        .bodies
        .iter()
        .map(|b| b.external_work)
        .sum::<f64>()
        + shared.motion.external_work;
    close(delta - authored + loss, 0.0, 3e-18);
    for residual in shared.motion.momentum_residual {
        close(residual, 0.0, 3e-16);
    }
    close(shared.motion.work_residual, 0.0, 3e-18);
    // Two copied inertias do NOT describe the same mechanism.
    let mut independent = patches.map(|p| finite(p, carrier));
    let mut other_spheres = initial;
    let mut other = grid();
    other
        .step_with_lubricated_partial_saturation_finite_planes(
            [&zero, &zero],
            [RADIUS; 2],
            &mut other_spheres,
            loads,
            model(2),
            &mut independent,
            [load, [0.0; 6]],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert!((independent[0].body().velocity[0] - independent[1].body().velocity[0]).abs() > 1e-4);
    assert!((other_spheres[1].velocity[0] - spheres[1].velocity[0]).abs() > 1e-4);
}

#[test]
fn sampled_fluid_full_tensor_and_prescribed_work_close_with_two_shared_faces() {
    let mut g = grid();
    let force = [2e-7, -1e-7, 0.5e-7];
    g.set_force_field(|cell| {
        if cell.into_iter().all(|i| (2..6).contains(&i)) {
            force
        } else {
            [0.0; 3]
        }
    });
    let mut sphere = body([3.5, 4.5, 4.5], 3.0, [-0.008, 0.004, -0.002]);
    sphere.angular_velocity = [0.003, -0.004, 0.002];
    let mut carrier = body([2.6, 3.6, 4.8], 5.0, [0.002, -0.001, 0.003]);
    carrier.angular_velocity = [-0.002, 0.003, -0.001];
    let mut spheres = [sphere];
    let mut carriers = [carrier];
    let mut sampler = BodyFrameSampler3::new(&g, 4).unwrap();
    sampler
        .sample(
            sphere.center,
            So3::identity(),
            &BodyFrameEllipsoid3::sphere(RADIUS).unwrap(),
            &CancelGate::new(),
        )
        .unwrap();
    let fraction = sampler.fractions().to_vec();
    let tangent = |point, normal| patch(point, normal).with_tangential(0.013, 0.001).unwrap();
    let surfaces = [
        attached(0, tangent([2.79, 4.0, 4.0], [1.0, 0.0, 0.0])),
        attached(0, tangent([3.0, 3.79, 4.0], [0.0, 1.0, 0.0])),
        SharedLubricationSurface3::Prescribed(
            LubricationPlane3::new(
                [4.21, 4.3, 4.6],
                [-1.0, 0.0, 0.0],
                [0.001, -0.002, 0.001],
                0.04,
                0.001,
            )
            .unwrap()
            .with_tangential(0.013, 0.001)
            .unwrap()
            .with_axial_spin(0.003)
            .unwrap(),
        ),
    ];
    let load = [0.0002, -0.0001, 0.0003, 0.00001, -0.00002, 0.00001];
    let carrier_load = [-0.0003, 0.0002, -0.0001, -0.00002, 0.00001, 0.00003];
    for _ in 0..2 {
        let before = snapshot(&g);
        let old_sphere = spheres[0];
        let old_carrier = carriers[0];
        let report = g
            .step_with_lubricated_partial_saturation_shared_carriers(
                [&fraction],
                [RADIUS],
                &mut spheres,
                [load],
                model(3),
                &surfaces,
                &mut carriers,
                &[carrier_load],
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(report.coupling.lubrication.active_resistance_modes, 15);
        assert!(report.coupling.bodies[0].hydrodynamic.coupled_cells > 0);
        let after = snapshot(&g);
        let mut balance = [0.0; 6];
        for index in 0..512 {
            let cell = [index % 8, (index / 8) % 8, index / 64];
            let mut dp: [f64; 3] = core::array::from_fn(|axis| {
                (0..Q3)
                    .map(|q| {
                        let e = [E3[q].0, E3[q].1, E3[q].2];
                        f64::from(e[axis]) * (after[index][q] - before[index][q])
                    })
                    .sum()
            });
            if cell.into_iter().all(|i| (2..6).contains(&i)) {
                for i in 0..3 {
                    dp[i] -= force[i];
                }
            }
            let dl = cross(
                core::array::from_fn(|i| cell[i] as f64 + 0.5 - ORIGIN[i]),
                dp,
            );
            for i in 0..3 {
                balance[i] += dp[i];
                balance[i + 3] += dl[i];
            }
        }
        for (before, after, load) in [
            (old_sphere, spheres[0], load),
            (old_carrier, carriers[0], carrier_load),
        ] {
            let p0 = momentum(before);
            let p1 = momentum(after);
            let source = external(before, load);
            for i in 0..6 {
                balance[i] += p1[i] - p0[i] - source[i];
            }
        }
        for i in 0..3 {
            close(balance[i] + report.surfaces[2].impulse[i], 0.0, 3e-12);
            close(
                balance[i + 3] + report.surfaces[2].angular_impulse[i]
                    - report.coupling.periodic.angular_momentum_change[i],
                0.0,
                1e-11,
            );
        }
        let shared = report.carriers[0];
        let mid = shared.motion;
        for patch in &report.surfaces[..2] {
            let orbital = cross(
                core::array::from_fn(|i| old_carrier.center[i] - ORIGIN[i]),
                patch.impulse,
            );
            let spin = core::array::from_fn(|i| patch.angular_impulse[i] - orbital[i]);
            close(
                patch.work,
                dot(patch.impulse, mid.velocity_midpoint)
                    + dot(spin, mid.angular_velocity_midpoint),
                2e-18,
            );
        }
        assert_eq!(shared.motion.body_after, carriers[0]);
        close(
            report.surfaces[..2].iter().map(|p| p.work).sum(),
            shared.work,
            2e-18,
        );
        close(
            energy(carriers[0]) - energy(old_carrier),
            shared.motion.kinetic_energy_change,
            2e-18,
        );
        close(shared.motion.work_residual, 0.0, 2e-18);
        close(report.coupling.lubrication.work_residual, 0.0, 2e-18);
        assert!(report.coupling.lubrication.dissipated_energy >= 0.0);
    }
}

#[test]
fn one_face_per_registry_body_retains_legacy_bits_and_explicit_owner_order() {
    let (initial, left, patches) = opposing();
    let right = PartialSaturationBody3 {
        mass: 7.0,
        velocity: [-0.002, 0.0, 0.0],
        ..left
    };
    let mut carriers = [left, right];
    let mut legacy = [finite(patches[1], right), finite(patches[0], left)];
    let surfaces = [attached(1, patches[1]), attached(0, patches[0])];
    let loads = [
        [0.0002, 0.0001, -0.0001, 0.00001, -0.00002, 0.00001],
        [-0.0001, 0.0002, 0.0001, -0.00002, 0.00001, 0.00003],
    ];
    let mut a = grid();
    let mut b = grid();
    let mut spheres_a = initial;
    let mut spheres_b = initial;
    let zero = vec![0.0; 512];
    let old = a
        .step_with_lubricated_partial_saturation_finite_planes(
            [&zero, &zero],
            [RADIUS; 2],
            &mut spheres_a,
            [[0.0; 6]; 2],
            model(2),
            &mut legacy,
            [loads[1], loads[0]],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let new = b
        .step_with_lubricated_partial_saturation_shared_carriers(
            [&zero, &zero],
            [RADIUS; 2],
            &mut spheres_b,
            [[0.0; 6]; 2],
            model(2),
            &surfaces,
            &mut carriers,
            &loads,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(new.coupling, old.coupling);
    assert_eq!(spheres_a, spheres_b);
    assert_eq!(snapshot(&a), snapshot(&b));
    for owner in 0..2 {
        assert_eq!(carriers[owner], legacy[1 - owner].body());
        assert_eq!(
            new.carriers[owner].motion,
            old.planes[1 - owner].finite_body.unwrap()
        );
        assert_eq!(new.carriers[owner].impulse, old.planes[1 - owner].impulse);
        assert_eq!(new.carriers[owner].work, old.planes[1 - owner].work);
    }
    assert_eq!(new.surfaces[0].carrier, Some(1));
    assert_eq!(new.surfaces[1].carrier, Some(0));
}

#[test]
fn invalid_bindings_and_failed_shared_trials_publish_nothing_and_retry_fresh() {
    let (initial, carrier, patches) = opposing();
    let mut spheres = initial;
    let mut carriers = vec![carrier];
    let mut g = grid();
    let original = snapshot(&g);
    let zero = vec![0.0; 512];
    let gate = CancelGate::new();
    let good = patches.map(|patch| attached(0, patch));
    let run = |g: &mut BoundaryGrid3,
               spheres: &mut [PartialSaturationBody3; 2],
               carriers: &mut [PartialSaturationBody3],
               loads: &[[f64; 6]],
               surfaces: &[SharedLubricationSurface3],
               model,
               gate: &CancelGate| {
        g.step_with_lubricated_partial_saturation_shared_carriers(
            [&zero, &zero],
            [RADIUS; 2],
            spheres,
            [[0.0; 6]; 2],
            model,
            surfaces,
            carriers,
            loads,
            ORIGIN,
            gate,
        )
    };
    for (loads, surfaces, model) in [
        (vec![], good.to_vec(), model(2)),
        (vec![[0.0; 6]; 2], good.to_vec(), model(2)),
        (vec![[0.0; 6]], vec![attached(1, patches[0])], model(2)),
        (vec![[0.0; 6]], vec![good[0], good[0]], model(2)),
        (
            vec![[0.0; 6]],
            good.to_vec(),
            model(2).with_resistance_mode_limit(1).unwrap(),
        ),
        (
            vec![[1e4, 0.0, 0.0, 0.0, 0.0, 0.0]],
            good.to_vec(),
            model(2),
        ),
    ] {
        assert!(
            run(
                &mut g,
                &mut spheres,
                &mut carriers,
                &loads,
                &surfaces,
                model,
                &gate
            )
            .is_err()
        );
        assert_eq!(spheres, initial);
        assert_eq!(carriers, vec![carrier]);
        assert_eq!(snapshot(&g), original);
    }
    let mut unused = [carrier, carrier];
    assert!(
        run(
            &mut g,
            &mut spheres,
            &mut unused,
            &[[0.0; 6]; 2],
            &good,
            model(2),
            &gate
        )
        .is_err()
    );
    assert_eq!(unused, [carrier; 2]);
    assert_eq!(snapshot(&g), original);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        run(
            &mut g,
            &mut spheres,
            &mut carriers,
            &[[0.0; 6]],
            &good,
            model(2),
            &cancelled
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(spheres, initial);
    assert_eq!(carriers, vec![carrier]);
    assert_eq!(snapshot(&g), original);
    let accepted = run(
        &mut g,
        &mut spheres,
        &mut carriers,
        &[[0.0; 6]],
        &good,
        model(2),
        &gate,
    )
    .unwrap();
    let mut fresh = grid();
    let mut fresh_spheres = initial;
    let mut fresh_carriers = [carrier];
    let repeated = run(
        &mut fresh,
        &mut fresh_spheres,
        &mut fresh_carriers,
        &[[0.0; 6]],
        &good,
        model(2),
        &gate,
    )
    .unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(spheres, fresh_spheres);
    assert_eq!(carriers, Vec::from(fresh_carriers));
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
