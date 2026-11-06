//! The real native flow and trajectory owners retain smaller physical changes.

use super::*;
use fs_lbm::d3q19::MovingPartialSaturationExchange3;

fn advance(
    grid: &mut BoundaryGrid3,
    sampler: &mut BodyFrameSampler3,
    shape: &impl BodyFrameShape3,
    body: &mut MovingPartialSaturationBody3,
    impulse: [f64; 6],
    assembly: bool,
) -> MovingPartialSaturationExchange3 {
    if assembly {
        grid.step_with_moving_partial_saturation_bodies(
            core::array::from_mut(sampler),
            [shape as &dyn BodyFrameShape3],
            core::array::from_mut(body),
            [impulse],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap()
        .bodies[0]
    } else {
        grid.step_with_moving_partial_saturation_body(
            sampler,
            shape,
            body,
            impulse,
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap()
    }
}

#[test]
fn original_couple_changes_native_spin_when_world_velocity_difference_loses_it() {
    for assembly in [false, true] {
        let mut grid = grid(BoundarySpec3::periodic());
        let mut body = body();
        body.center = [4.0; 3];
        body.velocity = [0.0; 3];
        body.inertia =
            InertiaTensor3::from_principal(Vec3::new(1.0, 1.0, 1.0), So3::identity()).unwrap();
        body.orientation = So3::exp(So3Tangent::new(Vec3::new(
            0.0,
            0.0,
            core::f64::consts::FRAC_PI_4,
        )))
        .unwrap();
        body.angular_momentum_body = Vec3::new(0.03125, 0.0, 0.0);
        let (midpoint, _) = body
            .inertia
            .momentum_step(
                body.orientation,
                body.angular_momentum_body,
                Vec3::new(0.0, 0.0, 0.0),
                0.5,
            )
            .unwrap();
        let impulse = 2.0_f64.powi(-62);
        let world = midpoint.rotate(Vec3::new(0.0, impulse, 0.0)).unwrap();
        let mut sampler = BodyFrameSampler3::new(&grid, 1).unwrap();
        let receipt = advance(
            &mut grid,
            &mut sampler,
            &EmptyShape,
            &mut body,
            [0.0, 0.0, 0.0, world.x, world.y, world.z],
            assembly,
        );
        close(body.angular_momentum_body.y / impulse, 1.0, 2e-15);
        assert_eq!(body.angular_momentum_body.x, 0.03125);
        close(
            receipt.kinetic_energy_change / (0.5 * impulse * impulse),
            1.0,
            4e-15,
        );
        assert_eq!(receipt.free_rotation_energy_change, 0.0);
        assert_eq!(receipt.coupling.hydrodynamic.body_work, 0.0);
    }
}

#[test]
fn moving_body_and_fluid_keep_transverse_energy_beside_large_common_translation() {
    for assembly in [false, true] {
        let mut grid = grid(BoundarySpec3::periodic());
        let mass = 2.0_f64.powi(80);
        let velocity = [2.0_f64.powi(-5), 0.0, 0.0];
        let f = equilibrium3(1.0, velocity);
        for z in 0..N {
            for y in 0..N {
                for x in 0..N {
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
        let mut body = MovingPartialSaturationBody3 {
            center: [4.0; 3],
            orientation: So3::identity(),
            mass,
            inertia: InertiaTensor3::from_principal(
                Vec3::new(0.4 * mass, 0.4 * mass, 0.4 * mass),
                So3::identity(),
            )
            .unwrap(),
            velocity,
            angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
        };
        let before = body;
        let populations = snapshot(&grid);
        let shape = BodyFrameEllipsoid3::sphere(1.0).unwrap();
        let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
        let transverse = 2.0_f64.powi(-32);
        let receipt = advance(
            &mut grid,
            &mut sampler,
            &shape,
            &mut body,
            [0.0, mass * transverse, 0.0, 0.0, 0.0, 0.0],
            assembly,
        );
        assert_eq!(body.velocity[0], velocity[0]);
        assert_eq!(body.velocity[1], transverse);
        close(body.velocity[2], 0.0, 1e-70);
        assert_eq!(body_energy(before), body_energy(body));
        close(
            receipt.kinetic_energy_change,
            0.5 * mass * transverse * transverse,
            0.0,
        );
        assert!(receipt.sampled_solid_volume[1] > 0.0);
        assert!(receipt.coupling.hydrodynamic.body_impulse[1] < 0.0);
        let (_, fluid, _) = fluid_change(&populations, &grid, [0.0; 3]);
        assert!(
            fluid[1] > 0.0,
            "the transverse kick must drive the actual fluid"
        );
        close(
            fluid[1] + receipt.coupling.hydrodynamic.body_impulse[1],
            0.0,
            5e-13,
        );
        close(receipt.work_residual, 0.0, 1e-12);
    }
}
