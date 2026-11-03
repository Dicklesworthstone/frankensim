//! G0/G3/G4/G5: actual moving geometry, native free motion and atomic coupling.
use std::cell::Cell;

use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, BoundarySpec3, E3,
    FaceBoundary3, MovingPartialSaturationBody3, PartialSaturationError3, Q3, equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const N: usize = 8;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; difference={:.4e}",
        a - b
    );
}

fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn index(x: usize, y: usize, z: usize) -> usize {
    (z * N + y) * N + x
}
fn grid(boundaries: BoundarySpec3) -> BoundaryGrid3 {
    BoundaryGrid3::new(N, N, N, 0.83, [0.0; 3], boundaries)
}

fn body() -> MovingPartialSaturationBody3 {
    let inertia = InertiaTensor3::new(Mat3 {
        m: [1.3, 0.12, -0.08, 0.12, 1.7, 0.06, -0.08, 0.06, 2.0],
    })
    .unwrap();
    MovingPartialSaturationBody3 {
        center: [3.7, 3.9, 4.1],
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.15))).unwrap(),
        mass: 0.7,
        inertia,
        velocity: [0.014, -0.008, 0.006],
        angular_momentum_body: inertia.momentum(Vec3::new(0.008, -0.006, 0.01)).unwrap(),
    }
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut values = Vec::with_capacity(N * N * N);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                values.push(if grid.is_solid(x, y, z) {
                    [0.0; Q3]
                } else {
                    grid.populations(x, y, z)
                });
            }
        }
    }
    values
}

fn fluid_change(
    before: &[[f64; Q3]],
    after: &BoundaryGrid3,
    origin: [f64; 3],
) -> (f64, [f64; 3], [f64; 3]) {
    let (mut mass, mut linear, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let mut local = [0.0; 3];
                let after = after.populations(x, y, z);
                for q in 0..Q3 {
                    let change = after[q] - before[index(x, y, z)][q];
                    mass += change;
                    local[0] += f64::from(E3[q].0) * change;
                    local[1] += f64::from(E3[q].1) * change;
                    local[2] += f64::from(E3[q].2) * change;
                }
                let arm = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let local_angular = cross(arm, local);
                for i in 0..3 {
                    linear[i] += local[i];
                    angular[i] += local_angular[i];
                }
            }
        }
    }
    (mass, linear, angular)
}

fn body_momenta(body: MovingPartialSaturationBody3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let p = body.velocity.map(|v| body.mass * v);
    let spin = array(body.orientation.rotate(body.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| body.center[i] - origin[i]), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}

fn body_energy(body: MovingPartialSaturationBody3) -> f64 {
    0.5 * body.mass * dot(body.velocity, body.velocity)
        + 0.5
            * dot(
                array(body.angular_momentum_body),
                array(body.inertia.velocity(body.angular_momentum_body).unwrap()),
            )
}

#[test]
fn ellipsoid_sampling_rotates_in_body_frame_and_keeps_raw_counts() {
    let grid = grid(BoundarySpec3::periodic());
    let shape = BodyFrameEllipsoid3::new(Vec3::new(1.47, 0.83, 0.61)).unwrap();
    let gate = CancelGate::new();
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    let volume = sampler
        .sample([4.0; 3], So3::identity(), &shape, &gate)
        .unwrap();
    let initial = sampler.fractions().to_vec();
    let quarter = So3::exp(So3Tangent::new(Vec3::new(
        0.0,
        0.0,
        core::f64::consts::FRAC_PI_2,
    )))
    .unwrap();
    let rotated_volume = sampler.sample([4.0; 3], quarter, &shape, &gate).unwrap();
    assert_eq!(volume, rotated_volume);
    assert_ne!(initial, sampler.fractions());
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                assert_eq!(
                    initial[index(x, y, z)],
                    sampler.fractions()[index(N - 1 - y, x, z)]
                );
            }
        }
    }
    assert!(
        sampler
            .fractions()
            .iter()
            .all(|v| (64.0 * v).fract() == 0.0)
    );
    let measured: f64 = sampler.fractions().iter().sum();
    assert_eq!(measured, volume);
    let sphere = BodyFrameEllipsoid3::sphere(1.27).unwrap();
    let mut finer = BodyFrameSampler3::new(&grid, 8).unwrap();
    let estimate = finer
        .sample([4.07, 3.93, 4.11], quarter, &sphere, &gate)
        .unwrap();
    let exact = 4.0 * core::f64::consts::PI * 1.27_f64.powi(3) / 3.0;
    close(estimate, exact, 0.04 * exact);
    for n in [0, 9, usize::MAX] {
        assert!(BodyFrameSampler3::new(&grid, n).is_err());
    }
    assert!(BodyFrameSampler3::new(&grid, 1).is_ok());
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(BodyFrameEllipsoid3::sphere(radius).is_err());
    }
}

#[test]
fn translated_ellipsoid_refreshes_native_coupling_and_replays_exactly() {
    let mut first = grid(BoundarySpec3::periodic());
    let mut second = grid(BoundarySpec3::periodic());
    let velocity = [0.035, 0.0, 0.007];
    let equilibrium = equilibrium3(1.02, velocity);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                first.set_populations(x, y, z, &equilibrium);
                second.set_populations(x, y, z, &equilibrium);
            }
        }
    }
    let shape = BodyFrameEllipsoid3::new(Vec3::new(1.13, 0.72, 0.51)).unwrap();
    let mut a = body();
    a.center = [3.0, 3.7, 4.0];
    a.velocity = velocity;
    a.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    let mut b = a;
    let initial = a;
    let mut sa = BodyFrameSampler3::new(&first, 3).unwrap();
    let mut sb = BodyFrameSampler3::new(&second, 3).unwrap();
    let mut independent = BodyFrameSampler3::new(&first, 3).unwrap();
    sa.sample(a.center, a.orientation, &shape, &CancelGate::new())
        .unwrap();
    let old_fractions = sa.fractions().to_vec();
    let mut volume_variation = 0.0;
    for _ in 0..32 {
        let ra = first
            .step_with_moving_partial_saturation_body(
                &mut sa,
                &shape,
                &mut a,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        let rb = second
            .step_with_moving_partial_saturation_body(
                &mut sb,
                &shape,
                &mut b,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(ra, rb);
        assert_eq!(a, b);
        let volume = independent
            .sample(
                ra.midpoint_center,
                ra.midpoint_orientation,
                &shape,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(independent.fractions(), sa.fractions());
        assert_eq!(volume, ra.sampled_solid_volume[1]);
        volume_variation += (ra.sampled_solid_volume[2] - ra.sampled_solid_volume[0]).abs();
        for impulse in ra.coupling.hydrodynamic.body_impulse {
            close(impulse, 0.0, 2e-13);
        }
        for i in 0..3 {
            close(a.velocity[i], velocity[i], 2e-13);
        }
    }
    assert!(
        volume_variation > 0.0,
        "finite quadrature volume variation must be visible"
    );
    assert_ne!(old_fractions, sa.fractions());
    assert_eq!(snapshot(&first), snapshot(&second));
    for i in 0..3 {
        close(a.center[i], initial.center[i] + 32.0 * velocity[i], 3e-12);
    }
    for cell in snapshot(&first) {
        for q in 0..Q3 {
            close(cell[q], equilibrium[q], 3e-13);
        }
    }
}

#[test]
fn moving_rotating_tensor_body_closes_actual_fluid_wall_and_body_momenta() {
    let mut grid = grid(BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]));
    grid.perturb(127, 0.001);
    let mut body = body();
    let initial = body;
    let shape = BodyFrameEllipsoid3::new(Vec3::new(1.2, 0.8, 0.6)).unwrap();
    let mut sampler = BodyFrameSampler3::new(&grid, 3).unwrap();
    let origin = [-2.0, 1.3, 5.7];
    let external = [0.0003, 0.0002, -0.0001, 0.0002, -0.0001, 0.0003];
    let mut wall_activity = 0.0;
    for _ in 0..24 {
        let populations = snapshot(&grid);
        let before = body;
        let (p0, h0) = body_momenta(before, origin);
        let receipt = grid
            .step_with_moving_partial_saturation_body(
                &mut sampler,
                &shape,
                &mut body,
                external,
                origin,
                &CancelGate::new(),
            )
            .unwrap();
        let (mass, fluid_p, fluid_h) = fluid_change(&populations, &grid, origin);
        let (p1, h1) = body_momenta(body, origin);
        close(mass, 0.0, 5e-13);
        let center_arm = core::array::from_fn(|i| receipt.midpoint_center[i] - origin[i]);
        let orbital = cross(center_arm, external[..3].try_into().unwrap());
        for i in 0..3 {
            close(
                fluid_p[i] + p1[i] - p0[i] + receipt.coupling.walls.wall_impulse[i],
                external[i],
                2e-11,
            );
            close(
                fluid_h[i] + h1[i] - h0[i] + receipt.coupling.walls.wall_angular_impulse[i],
                external[i + 3] + orbital[i],
                3e-10,
            );
            close(
                receipt.external_angular_impulse[i],
                external[i + 3] + orbital[i],
                2e-18,
            );
            wall_activity += receipt.coupling.walls.wall_impulse[i].abs();
        }
        for residual in receipt.momentum_residual {
            close(residual, 0.0, 2e-14);
        }
        close(
            body_energy(body) - body_energy(before),
            receipt.kinetic_energy_change,
            2e-17,
        );
        close(receipt.work_residual, 0.0, 2e-14);
        close(
            receipt.kinetic_energy_change,
            receipt.coupling.hydrodynamic.body_work
                + receipt.coupling.external_work
                + receipt.free_rotation_energy_change,
            2e-14,
        );
    }
    assert_ne!(body.center, initial.center);
    assert_ne!(body.orientation, initial.orientation);
    assert_ne!(body.angular_momentum_body, initial.angular_momentum_body);
    assert!(wall_activity > 1e-4);
}

struct EmptyShape;
impl BodyFrameShape3 for EmptyShape {
    fn bounding_radius(&self) -> f64 {
        0.5
    }
    fn contains(&self, _: Vec3) -> bool {
        false
    }
}

#[test]
fn zero_occupancy_recovers_native_free_rotation_and_ordinary_fluid_bits() {
    let boundaries = BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]);
    let mut coupled = grid(boundaries);
    let mut reference = grid(boundaries);
    coupled.perturb(71, 0.001);
    reference.perturb(71, 0.001);
    let mut body = body();
    body.velocity = [0.0; 3];
    let initial_spin = body.orientation.rotate(body.angular_momentum_body).unwrap();
    let mut sampler = BodyFrameSampler3::new(&coupled, 1).unwrap();
    for _ in 0..64 {
        let mut expected = (body.orientation, body.angular_momentum_body);
        for _ in 0..2 {
            expected = body
                .inertia
                .momentum_step(expected.0, expected.1, Vec3::new(0.0, 0.0, 0.0), 0.5)
                .unwrap();
        }
        reference.step();
        let receipt = coupled
            .step_with_moving_partial_saturation_body(
                &mut sampler,
                &EmptyShape,
                &mut body,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!((body.orientation, body.angular_momentum_body), expected);
        assert_eq!(snapshot(&coupled), snapshot(&reference));
        assert_eq!(receipt.sampled_solid_volume, [0.0; 3]);
        close(
            receipt.kinetic_energy_change,
            receipt.free_rotation_energy_change,
            2e-18,
        );
    }
    let final_spin = body.orientation.rotate(body.angular_momentum_body).unwrap();
    for (a, b) in array(initial_spin).into_iter().zip(array(final_spin)) {
        close(a, b, 2e-14);
    }
}

#[test]
fn late_swept_support_refusal_keeps_fluid_body_and_retry_unchanged() {
    let mut grid = grid(BoundarySpec3::periodic());
    let shape = BodyFrameEllipsoid3::sphere(1.0).unwrap();
    let mut body = body();
    body.center = [2.001, 4.0, 4.0];
    body.velocity = [0.0; 3];
    body.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    let initial = body;
    let before = snapshot(&grid);
    let mut sampler = BodyFrameSampler3::new(&grid, 2).unwrap();
    let result = grid.step_with_moving_partial_saturation_body(
        &mut sampler,
        &shape,
        &mut body,
        [-0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0; 3],
        &CancelGate::new(),
    );
    assert!(matches!(
        result,
        Err(PartialSaturationError3::Invalid {
            reason: "moving body swept support reaches the domain halo"
        })
    ));
    assert_eq!(body, initial);
    assert_eq!(snapshot(&grid), before);
    // Topology admission must still be open after the candidate was streamed.
    grid.voxelize_sdf(|_| 1.0);
    let mut reference = crate::grid(BoundarySpec3::periodic());
    let mut reference_body = initial;
    let mut reference_sampler = BodyFrameSampler3::new(&reference, 2).unwrap();
    let receipt = grid
        .step_with_moving_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    let expected = reference
        .step_with_moving_partial_saturation_body(
            &mut reference_sampler,
            &shape,
            &mut reference_body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(receipt, expected);
    assert_eq!(body, reference_body);
    assert_eq!(snapshot(&grid), snapshot(&reference));
}

struct CancelOnEndpoint<'a> {
    gate: &'a CancelGate,
    calls: Cell<usize>,
}
impl BodyFrameShape3 for CancelOnEndpoint<'_> {
    fn bounding_radius(&self) -> f64 {
        1.0
    }
    fn contains(&self, point: Vec3) -> bool {
        let count = self.calls.get() + 1;
        self.calls.set(count);
        // At [4,4,4], each pose visits 8 cells * 8 subcell samples.
        if count == 129 {
            self.gate.request();
        }
        dot(array(point), array(point)) <= 1.0
    }
}

#[test]
fn endpoint_cancellation_and_fixed_obstacle_refuse_atomically() {
    let mut grid = grid(BoundarySpec3::periodic());
    let mut body = body();
    body.center = [4.0; 3];
    body.orientation = So3::identity();
    body.velocity = [0.0; 3];
    body.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    let initial = body;
    let before = snapshot(&grid);
    let mut sampler = BodyFrameSampler3::new(&grid, 2).unwrap();
    let gate = CancelGate::new();
    let shape = CancelOnEndpoint {
        gate: &gate,
        calls: Cell::new(0),
    };
    let result = grid.step_with_moving_partial_saturation_body(
        &mut sampler,
        &shape,
        &mut body,
        [0.0; 6],
        [0.0; 3],
        &gate,
    );
    assert_eq!(result, Err(PartialSaturationError3::Cancelled));
    assert!(
        shape.calls.get() >= 129,
        "must cancel only after native collision and stream preparation"
    );
    assert_eq!(body, initial);
    assert_eq!(snapshot(&grid), before);
    grid.voxelize_sdf(|p| if p == [3.5; 3] { -1.0 } else { 1.0 });
    let blocked = snapshot(&grid);
    let result = grid.step_with_moving_partial_saturation_body(
        &mut sampler,
        &BodyFrameEllipsoid3::sphere(1.0).unwrap(),
        &mut body,
        [0.0; 6],
        [0.0; 3],
        &CancelGate::new(),
    );
    assert!(matches!(
        result,
        Err(PartialSaturationError3::Cell {
            reason: "moving body swept support meets a fixed solid cell",
            ..
        })
    ));
    assert_eq!(body, initial);
    assert_eq!(snapshot(&grid), blocked);
}

fn sampled_first_moment(sampler: &BodyFrameSampler3, index: usize) -> [f64; 3] {
    sampler
        .cut_cell_first_moments()
        .iter()
        .find(|(cell, _)| *cell == index)
        .map_or([0.0; 3], |(_, moment)| *moment)
}

#[test]
fn cut_cell_first_moments_match_independent_sample_points_and_clear_old_support() {
    let grid = grid(BoundarySpec3::periodic());
    let shape = BodyFrameEllipsoid3::new(Vec3::new(1.47, 0.83, 0.61)).unwrap();
    let orientation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    let gate = CancelGate::new();
    for samples in [1, 3, 8] {
        let mut sampler = BodyFrameSampler3::new(&grid, samples).unwrap();
        for center in [[3.37, 4.11, 3.83], [5.1, 3.0, 4.5]] {
            let volume = sampler.sample(center, orientation, &shape, &gate).unwrap();
            let mut expected_volume = vec![0.0; N * N * N];
            let mut expected_moment = vec![[0.0; 3]; N * N * N];
            let weight = 1.0 / (samples * samples * samples) as f64;
            // Independently integrate the WORLD sample points, without using
            // the sampler's fractions, sparse cache or integer moment sums.
            for gz in 0..N * samples {
                for gy in 0..N * samples {
                    for gx in 0..N * samples {
                        let fine = [gx, gy, gz];
                        let point = fine.map(|v| (v as f64 + 0.5) / samples as f64);
                        let local = orientation
                            .inverse()
                            .rotate(Vec3::new(
                                point[0] - center[0],
                                point[1] - center[1],
                                point[2] - center[2],
                            ))
                            .unwrap();
                        if shape.contains(local) {
                            let cell = fine.map(|v| v / samples);
                            let index = index(cell[0], cell[1], cell[2]);
                            expected_volume[index] += weight;
                            for i in 0..3 {
                                expected_moment[index][i] +=
                                    weight * (point[i] - cell[i] as f64 - 0.5);
                            }
                        }
                    }
                }
            }
            close(volume, expected_volume.iter().sum(), 2e-13);
            let mut largest_offset = 0.0_f64;
            for index in 0..N * N * N {
                close(sampler.fractions()[index], expected_volume[index], 2e-14);
                let measured = sampled_first_moment(&sampler, index);
                for i in 0..3 {
                    close(measured[i], expected_moment[index][i], 2e-14);
                    largest_offset = largest_offset.max(measured[i].abs());
                }
                if sampler.fractions()[index] == 0.0 || sampler.fractions()[index] == 1.0 {
                    assert_eq!(measured, [0.0; 3]);
                }
            }
            if samples > 1 {
                assert!(
                    largest_offset > 0.02,
                    "off-center cuts must exercise the moment"
                );
            }
        }
    }

    // Two disjoint half-cell geometries retain opposite local moments even
    // though their total occupancy fills this one cell exactly.
    struct HalfCell;
    impl BodyFrameShape3 for HalfCell {
        fn bounding_radius(&self) -> f64 {
            0.75
        }
        fn contains(&self, p: Vec3) -> bool {
            p.x.abs() < 0.25 && p.y.abs() < 0.5 && p.z.abs() < 0.5
        }
    }
    let mut halves = [
        BodyFrameSampler3::new(&grid, 8).unwrap(),
        BodyFrameSampler3::new(&grid, 8).unwrap(),
    ];
    for (side, sampler) in halves.iter_mut().enumerate() {
        sampler
            .sample(
                [3.25 + 0.5 * side as f64, 3.5, 3.5],
                So3::identity(),
                &HalfCell,
                &gate,
            )
            .unwrap();
        let sign = if side == 0 { -1.0 } else { 1.0 };
        assert_eq!(sampler.fractions()[index(3, 3, 3)], 0.5);
        assert_eq!(
            sampled_first_moment(sampler, index(3, 3, 3)),
            [sign * 0.125, 0.0, 0.0]
        );
        let torque = cross(
            sampled_first_moment(sampler, index(3, 3, 3)),
            [0.0, 1.0, 0.0],
        );
        assert_eq!(torque, [0.0, 0.0, sign * 0.125]);
    }
}
