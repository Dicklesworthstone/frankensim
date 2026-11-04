//! Original rheology drives one native moving body and its resolved fluid.
use std::cell::Cell;

use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, BoundarySpec3,
    CollisionModel3, E3, MovingPartialSaturationBody3, PartialSaturationError3, Q3, W3,
    equilibrium3,
};
use fs_lbm::kinematics::KinematicsWorkspace3;
use fs_lbm::material::immersed::{MaterialImmersedError3, MaterialImmersedPoint3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};
use fs_time::lie::gyrostat::InertiaTensor3;

#[path = "material_immersed_rheology/multiple.rs"]
mod multiple;

const N: usize = 8;
const CELLS: usize = N * N * N;
const DX: f64 = 0.01;
const DT: f64 = 0.01;

fn transport(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1000.0, mu).unwrap()
}

fn model(tau: f64, family: usize) -> CollisionModel3 {
    match family {
        0 => CollisionModel3::Bgk { tau },
        1 => CollisionModel3::Trt { tau, magic: 0.137 },
        2 => CollisionModel3::CentralMoment {
            second_order_rate: 1.0 / tau,
            higher_order_rate: 0.73,
        },
        _ => CollisionModel3::ReducedCumulant {
            second_order_rate: 1.0 / tau,
            third_order_rate: 0.87,
            fourth_order_rate: 1.63,
        },
    }
}

fn grid(model: CollisionModel3, force: [f64; 3]) -> BoundaryGrid3 {
    BoundaryGrid3::with_collision_model(N, N, N, model, force, BoundarySpec3::periodic())
}

fn body() -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [4.0; 3],
        orientation: So3::identity(),
        mass: 0.7,
        inertia: InertiaTensor3::new(Mat3 {
            m: [1.3, 0.0, 0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 2.0],
        })
        .unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
    }
}

fn workspace() -> KinematicsWorkspace3 {
    KinematicsWorkspace3::new(
        [N; 3],
        KinematicsWorkspace3::required_bytes([N; 3]).unwrap(),
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..CELLS)
        .map(|i| grid.populations(i % N, i / N % N, i / (N * N)))
        .collect()
}

// Independently prescribe the physical Guo-frame velocity. The non-equilibrium
// perturbation adds stress/higher moments but exactly zero mass and momentum.
fn seed_velocity(grid: &mut BoundaryGrid3, velocity: impl Fn([usize; 3]) -> [f64; 3]) {
    for i in 0..CELLS {
        let cell = [i % N, i / N % N, i / (N * N)];
        let u = velocity(cell);
        let force = grid.force_density(cell[0], cell[1], cell[2]);
        let uf: f64 = (0..3).map(|a| u[a] * force[a]).sum();
        let mut f = equilibrium3(1.0, u);
        f[1] += 0.00012;
        f[2] += 0.00012;
        f[0] -= 0.00024;
        for q in 0..Q3 {
            let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
            let eu: f64 = (0..3).map(|a| e[a] * u[a]).sum();
            let ef: f64 = (0..3).map(|a| e[a] * force[a]).sum();
            f[q] -= 0.5 * W3[q] * (3.0 * ef - 3.0 * uf + 9.0 * eu * ef);
        }
        grid.set_populations(cell[0], cell[1], cell[2], &f);
    }
}

#[derive(Debug)]
struct SingleRead {
    point: MaterialImmersedPoint3,
    value: FluidTransportState,
    reads: Cell<usize>,
}
impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        assert_eq!(self.reads.replace(self.reads.get() + 1), 0);
        Ok(self.value)
    }
}

#[test]
fn constant_source_matches_original_moving_solver_for_every_collision_family() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let tau = scaling.flow().tau;
    let shape = BodyFrameEllipsoid3::sphere(1.13).unwrap();
    let gate = CancelGate::new_clock_free();
    for family in 0..4 {
        for force in [[0.0; 3], [2e-6, -1e-6, 5e-7]] {
            let mut actual = grid(model(1.37, family), force);
            let mut expected = grid(model(tau, family), force);
            seed_velocity(&mut actual, |_| [0.024, -0.017, 0.011]);
            seed_velocity(&mut expected, |_| [0.024, -0.017, 0.011]);
            let mut a = body();
            a.velocity = [-0.011, 0.007, -0.004];
            a.angular_momentum_body = Vec3::new(0.005, -0.003, 0.004);
            let mut b = a;
            let mut sa = BodyFrameSampler3::new(&actual, 3).unwrap();
            let mut sb = BodyFrameSampler3::new(&expected, 3).unwrap();
            let impulse = [1e-4, -2e-4, 3e-4, 2e-5, -1e-5, 3e-5];
            let origin = [0.2, -0.7, 0.4];
            let mut queries = 0;
            let report = scaling
                .step_moving_body_with_material_rheology(
                    &mut actual,
                    &mut sa,
                    &shape,
                    &mut a,
                    impulse,
                    origin,
                    &mut workspace(),
                    CELLS,
                    None,
                    &gate,
                    |point| {
                        queries += 1;
                        Ok(SingleRead {
                            point,
                            // Momentum-only sampling must not impose thermal tau.
                            value: FluidTransportState::try_new(1000.0, 1000.0, 1e-250, 1.0)?,
                            reads: Cell::new(0),
                        })
                    },
                )
                .unwrap();
            let expected_report = expected
                .step_with_moving_partial_saturation_body(
                    &mut sb, &shape, &mut b, impulse, origin, &gate,
                )
                .unwrap();
            assert_eq!(queries, CELLS);
            assert_eq!(report.momentum, expected_report);
            assert_eq!(a, b);
            assert_eq!(snapshot(&actual), snapshot(&expected));
            assert_eq!(actual.collision_model(), model(1.37, family));
            for sample in report.applied {
                assert_eq!(sample.point, sample.material.point);
                assert_eq!(sample.material.reads.get(), 1);
                assert_eq!(sample.shear_relaxation_time, tau);
            }
        }
    }
}

fn carreau(rate: f64, time: f64) -> f64 {
    0.2 + 1.4 * (1.0 + (time * rate).powi(2)).powf(-0.25)
}

#[test]
fn variable_carreau_viscosity_closes_independent_body_response_and_streamed_momentum() {
    let gate = CancelGate::new_clock_free();
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let shape = BodyFrameEllipsoid3::sphere(1.3).unwrap();
    let force = [2e-6, 0.0, 0.0];
    let speed = |cell: [usize; 3]| 0.001 * (cell[1] as f64 + 0.5 - 4.0).powi(2);
    let mut accepted_speeds = Vec::new();
    for time in [2.0, 12.0] {
        let mut fluid = grid(model(1.37, 2), force);
        seed_velocity(&mut fluid, |cell| [speed(cell), 0.0, 0.0]);
        let before = snapshot(&fluid);
        let mut carrier = body();
        let initial = carrier;
        let mut sampler = BodyFrameSampler3::new(&fluid, 4).unwrap();
        sampler
            .sample(carrier.center, carrier.orientation, &shape, &gate)
            .unwrap();
        let mut resistance = 0.0;
        let mut drive = 0.0;
        // The symmetric sphere and even-in-y quadratic velocity have zero
        // torque. Thus the independent spatial response reduces to one scalar.
        for (i, &epsilon) in sampler.fractions().iter().enumerate() {
            if epsilon == 0.0 {
                continue;
            }
            let cell = [i % N, i / N % N, i / (N * N)];
            let rate = (0.002 * (cell[1] as f64 + 0.5 - 4.0)).abs() / DT;
            let tau = 0.5 + 3.0 * (carreau(rate, time) / 1000.0) * DT / (DX * DX);
            let weight = epsilon * (tau - 0.5) / (1.0 - epsilon + tau - 0.5);
            resistance += weight;
            drive += weight * speed(cell);
        }
        let expected_speed = drive / (carrier.mass + 0.5 * resistance);
        let report = scaling
            .step_moving_body_with_material_rheology(
                &mut fluid,
                &mut sampler,
                &shape,
                &mut carrier,
                [0.0; 6],
                [0.0; 3],
                &mut workspace(),
                CELLS,
                None,
                &gate,
                |point| Ok(transport(carreau(point.shear_rate_per_s, time))),
            )
            .unwrap();
        close(carrier.velocity[0], expected_speed, 3e-14);
        close(
            carrier.center[0],
            initial.center[0] + 0.5 * expected_speed,
            2e-14,
        );
        for value in carrier.velocity[1..]
            .iter()
            .copied()
            .chain(carrier.angular_velocity_world().unwrap())
        {
            close(value, 0.0, 3e-14);
        }
        close(
            report.momentum.coupling.hydrodynamic.coupling_weight_sum,
            resistance,
            2e-12,
        );
        close(report.momentum.work_residual, 0.0, 2e-14);
        for sample in &report.applied {
            let y = sample.point.cell[1];
            if (1..N - 1).contains(&y) {
                close(
                    sample.point.shear_rate_per_s,
                    (0.002 * (y as f64 + 0.5 - 4.0)).abs() / DT,
                    4e-13,
                );
            }
        }
        let after = snapshot(&fluid);
        let mut mass_change = 0.0;
        let mut momentum_change = [0.0; 3];
        for (old, new) in before.iter().zip(&after) {
            for q in 0..Q3 {
                let delta = new[q] - old[q];
                mass_change += delta;
                let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
                for a in 0..3 {
                    momentum_change[a] += e[a] * delta;
                }
            }
        }
        close(mass_change, 0.0, 3e-12);
        for a in 0..3 {
            close(
                momentum_change[a] + initial.mass * carrier.velocity[a],
                CELLS as f64 * force[a],
                4e-12,
            );
            close(report.momentum.momentum_residual[a], 0.0, 3e-14);
        }
        accepted_speeds.push(carrier.velocity[0]);
    }
    assert!((accepted_speeds[0] - accepted_speeds[1]).abs() > 1e-6);
}

#[test]
fn source_rate_includes_every_symmetric_component_and_removes_rigid_rotation() {
    let gate = CancelGate::new_clock_free();
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let shape = BodyFrameEllipsoid3::sphere(1.0).unwrap();
    let symmetric = [
        [0.0007, 0.0015, -0.0009],
        [0.0015, -0.0004, 0.0008],
        [-0.0009, 0.0008, -0.0003],
    ];
    for strain in [false, true] {
        for spin in [0.0, 0.003] {
            let mut g = if strain { symmetric } else { [[0.0; 3]; 3] };
            g[0][1] -= spin;
            g[1][0] += spin;
            let expected_rate = if strain {
                (2.0_f64 * symmetric.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt() / DT
            } else {
                0.0
            };
            let mut fluid = grid(model(1.37, 1), [2e-6, -1e-6, 3e-6]);
            seed_velocity(&mut fluid, |cell| {
                g.map(|row| (0..3).map(|a| row[a] * (cell[a] as f64 + 0.5 - 4.0)).sum())
            });
            let mut carrier = body();
            let mut sampler = BodyFrameSampler3::new(&fluid, 2).unwrap();
            let report = scaling
                .step_moving_body_with_material_rheology(
                    &mut fluid,
                    &mut sampler,
                    &shape,
                    &mut carrier,
                    [0.0; 6],
                    [0.0; 3],
                    &mut workspace(),
                    CELLS,
                    None,
                    &gate,
                    |point| Ok(transport(0.4 + 0.2 * point.shear_rate_per_s)),
                )
                .unwrap();
            for sample in report.applied {
                if sample
                    .point
                    .cell
                    .into_iter()
                    .all(|c| (1..N - 1).contains(&c))
                {
                    close(sample.point.shear_rate_per_s, expected_rate, 6e-13);
                    close(
                        sample.transport.dynamic_viscosity_pa_s(),
                        0.4 + 0.2 * expected_rate,
                        2e-13,
                    );
                }
            }
        }
    }
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
        // Each of the first two stationary poses visits 8 cells * 8 samples.
        if count == 129 {
            self.gate.request();
        }
        point.x * point.x + point.y * point.y + point.z * point.z <= 1.0
    }
}

#[test]
fn source_budget_and_late_geometry_refusals_leave_one_accepted_state() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let mut fluid = grid(model(1.37, 3), [0.0; 3]);
    let mut carrier = body();
    let initial = carrier;
    let saved = snapshot(&fluid);
    let shape = BodyFrameEllipsoid3::sphere(1.0).unwrap();
    let mut sampler = BodyFrameSampler3::new(&fluid, 2).unwrap();
    let mut cache = workspace();
    for mode in 0..4 {
        let gate = CancelGate::new_clock_free();
        let mut queries = 0;
        let error = scaling
            .step_moving_body_with_material_rheology(
                &mut fluid,
                &mut sampler,
                &shape,
                &mut carrier,
                [0.0; 6],
                [0.0; 3],
                &mut cache,
                if mode == 0 { CELLS - 1 } else { CELLS },
                None,
                &gate,
                |point| {
                    queries += 1;
                    if point.cell == [N - 1; 3] {
                        if mode == 1 {
                            return Err(MaterialError::State {
                                what: "original source domain exhausted".into(),
                            });
                        }
                        if mode == 2 {
                            return FluidTransportState::try_new(999.0, 1000.0, 1000.0, 1.0);
                        }
                        if mode == 3 {
                            gate.request();
                        }
                    }
                    Ok(reference)
                },
            )
            .unwrap_err();
        if mode == 0 {
            assert_eq!(queries, 0);
            assert!(matches!(error, MaterialImmersedError3::Invalid(_)));
        } else if mode == 3 {
            assert!(matches!(
                error,
                MaterialImmersedError3::Flow(PartialSaturationError3::Cancelled)
            ));
        } else {
            assert!(
                matches!(error, MaterialImmersedError3::State { cell, .. } if cell == [N - 1; 3])
            );
        }
        assert_eq!(snapshot(&fluid), saved);
        assert_eq!(carrier, initial);
    }
    let gate = CancelGate::new_clock_free();
    let cancelling = CancelOnEndpoint {
        gate: &gate,
        calls: Cell::new(0),
    };
    let error = scaling
        .step_moving_body_with_material_rheology(
            &mut fluid,
            &mut sampler,
            &cancelling,
            &mut carrier,
            [0.0; 6],
            [0.0; 3],
            &mut cache,
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        MaterialImmersedError3::Flow(PartialSaturationError3::Cancelled)
    ));
    assert!(cancelling.calls.get() >= 129);
    assert_eq!(snapshot(&fluid), saved);
    assert_eq!(carrier, initial);
    // Reusing both native and material scratch after a late refusal is exact.
    let mut fresh = grid(model(1.37, 3), [0.0; 3]);
    let mut fresh_body = initial;
    let mut fresh_sampler = BodyFrameSampler3::new(&fresh, 2).unwrap();
    let gate = CancelGate::new_clock_free();
    let actual = scaling
        .step_moving_body_with_material_rheology(
            &mut fluid,
            &mut sampler,
            &shape,
            &mut carrier,
            [0.0; 6],
            [0.0; 3],
            &mut cache,
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap();
    let expected = scaling
        .step_moving_body_with_material_rheology(
            &mut fresh,
            &mut fresh_sampler,
            &shape,
            &mut fresh_body,
            [0.0; 6],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(carrier, fresh_body);
    assert_eq!(snapshot(&fluid), snapshot(&fresh));
}
