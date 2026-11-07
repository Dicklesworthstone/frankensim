use super::*;
use crate::compressible::spatial::molecular::sgs::{
    ChannelDynamicScalarConfig, ChannelDynamicSgsConfig, ChannelDynamicThermalSubgridConfig,
    EquilibriumWallConfig,
};
use fs_lbm::conjugate::turbulence::law_of_the_wall_ratios;

fn wall_config() -> EquilibriumWallConfig {
    EquilibriumWallConfig {
        maximum_wall_distance_plus: 300.0,
        maximum_mach_number: 0.3,
    }
}

fn wall(velocity: [f64; 3]) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity,
        thermal: WallThermalCondition::Adiabatic,
    }
}

fn channel(axis: usize, velocity: [f64; 3]) -> Euler3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[axis] = WALL;
    field(
        gas(0.0),
        CartesianGrid::new(4, 4, 4, 0.2, 0.3, 0.4).unwrap(),
        Boundaries3d {
            x: axes[0],
            y: axes[1],
            z: axes[2],
        },
        |_| [1.0, velocity[0], velocity[1], velocity[2], 2.0],
    )
}

fn operator(
    flow: &Euler3d,
    axis: usize,
    velocity: [f64; 3],
    dynamic: bool,
) -> MolecularTransport3d<ConstantTransport> {
    let source = law(flow.gas(), 1e-4, 0.0, 0.0);
    if dynamic {
        let cfg = ChannelDynamicThermalSubgridConfig {
            momentum: ChannelDynamicSgsConfig {
                maximum_coefficient_squared: 100.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            heat: ChannelDynamicScalarConfig {
                maximum_coefficient: 100.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        };
        let bytes = MolecularTransport3d::<ConstantTransport>::required_channel_dynamic_thermal_sgs_workspace_bytes(flow).unwrap();
        MolecularTransport3d::with_channel_dynamic_thermal_sgs(
            flow,
            source,
            cfg,
            [wall(velocity); 2],
            bytes,
            || false,
        )
        .unwrap()
    } else {
        MolecularTransport3d::with_sgs(
            flow,
            source,
            config(SgsModel::Smagorinsky { cs: 0.0 }),
            || false,
        )
        .unwrap()
        .with_wall(BoundarySide::ALL[2 * axis], wall(velocity))
        .unwrap()
        .with_wall(BoundarySide::ALL[2 * axis + 1], wall(velocity))
        .unwrap()
    }
}

fn prepare(op: &mut MolecularTransport3d<ConstantTransport>, flow: &Euler3d) -> f64 {
    spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        None,
        &flow.states,
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.viscosity,
        &mut op.work.flux,
        &mut op.work.binary,
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn equilibrium_sgs_wall_native_faces_use_secant_stress_wall_work_and_stability_on_every_axis() {
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[a] = 0.2;
        velocity[b] = -0.1;
        let wall_velocity = velocity.map(|u| 0.25 * u);
        let flow = channel(axis, velocity);
        let distance = 0.5 * flow.grid.spacing[axis];
        let speed = (0.15_f64.powi(2) + 0.075_f64.powi(2)).sqrt();
        let ratio = law_of_the_wall_ratios(speed * distance / 1e-4).1;
        assert!(ratio > 1.5);
        let mut base = operator(&flow, axis, wall_velocity, false);
        let base_rate = prepare(&mut base, &flow);
        for dynamic in [false, true] {
            let mut op = operator(&flow, axis, wall_velocity, dynamic);
            for upper in [false, true] {
                let side = BoundarySide::ALL[2 * axis + usize::from(upper)];
                op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
                assert_eq!(op.equilibrium_sgs_wall(side), Some(wall_config()));
            }
            assert!(prepare(&mut op, &flow) > base_rate);
            for upper in [false, true] {
                let sign = if upper { 1.0 } else { -1.0 };
                let mut face = [1; 3];
                face[axis] = if upper { 4 } else { 0 };
                let flux = op.work.flux[axis][op.grid.cell_face(axis, face)];
                close(flux.normal_momentum, 0.0, 1e-14);
                for (i, t) in [a, b].into_iter().enumerate() {
                    let momentum =
                        1e-4 * ratio * sign * (velocity[t] - wall_velocity[t]) / distance;
                    close(flux.transverse_momentum[i], momentum, 2e-14);
                }
                close(
                    flux.energy,
                    flux.transverse_momentum[0] * wall_velocity[a]
                        + flux.transverse_momentum[1] * wall_velocity[b],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn equilibrium_sgs_wall_braking_keeps_native_total_energy_and_moving_wall_work() {
    for wall_speed in [0.0, 0.1] {
        let mut base = channel(1, [0.4, 0.0, 0.0]);
        let mut active = channel(1, [0.4, 0.0, 0.0]);
        let mut base_op = operator(&base, 1, [wall_speed, 0.0, 0.0], false);
        let mut op = operator(&active, 1, [wall_speed, 0.0, 0.0], false);
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
        }
        let before = active.totals();
        let kinetic_before = partition(&active).0;
        base_op
            .advance(
                &mut base,
                Time::new(0.02),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        let report = op
            .advance(
                &mut active,
                Time::new(0.02),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        assert!(active.totals().momentum[0] < base.totals().momentum[0]);
        assert!(partition(&active).0 < kinetic_before);
        close(
            report.boundary_change.energy,
            wall_speed * report.boundary_change.momentum[0],
            1e-13,
        );
        close(
            active.totals().energy - before.energy,
            report.boundary_change.energy,
            2e-13,
        );
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 2e-13);
        }
        assert!(report.maximum_parabolic_number <= 0.40000000000001);
    }
}

#[test]
fn equilibrium_sgs_wall_late_domain_refusal_preserves_native_field_and_accepted_loads() {
    let mut flow = channel(1, [0.4, 0.0, 0.0]);
    let mut op = operator(&flow, 1, [0.0; 3], false);
    for side in [BoundarySide::YMin, BoundarySide::YMax] {
        op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
    }
    op.advance(
        &mut flow,
        Time::new(0.01),
        MolecularControl::default(),
        || false,
    )
    .unwrap();
    let accepted = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    // Lower wall faces prepare first; the upper wall refuses after those
    // transient fluxes have already been calculated.
    op = op
        .with_equilibrium_sgs_wall(
            BoundarySide::YMax,
            EquilibriumWallConfig {
                maximum_mach_number: 1e-6,
                ..wall_config()
            },
        )
        .unwrap();
    assert!(
        op.advance(
            &mut flow,
            Time::new(0.01),
            MolecularControl::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), accepted);
    assert_eq!(op.boundary_impulses, loads);
    assert!(
        op.with_wall(BoundarySide::YMin, MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
}
