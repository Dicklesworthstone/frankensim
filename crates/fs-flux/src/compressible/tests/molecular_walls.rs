//! Source-owned planar no-slip wall momentum, heat and prescribed mechanical work.
use super::molecular::{bank, law};
use super::*;
use crate::compressible::molecular::{
    MolecularControl, MolecularError, MolecularEuler2d, MolecularStepControl, MolecularTransport2d,
    MolecularWall2d, WallThermalCondition,
};
use crate::compressible::spatial::{
    Boundaries3d, BoundarySide as Side3d, CartesianGrid as Grid3d, Euler3d, State as State3d,
    molecular::{MolecularTransport3d, MolecularWall3d},
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn channel_boundaries(axis: usize) -> Boundaries2d {
    let mut boundaries = PERIODIC;
    if axis == 0 {
        boundaries.x = WALLS;
    } else {
        boundaries.y = WALLS;
    }
    boundaries
}

fn no_slip(axis: usize, speed: f64, thermal: WallThermalCondition) -> MolecularWall2d {
    let mut velocity = [0.0; 2];
    velocity[1 - axis] = speed;
    MolecularWall2d::NoSlip { velocity, thermal }
}

#[test]
fn all_four_sides_resolve_half_cell_shear_fourier_heat_and_prescribed_wall_power() {
    let grid = CartesianGrid::new(3, 4, 0.2, 0.3, 0.7).unwrap();
    let (mu, conductivity, duration) = (0.02, 0.05, 1e-8);
    for side in SIDES {
        let axis = side.axis();
        let tangent = 1 - axis;
        let mut velocity = [0.0; 2];
        velocity[tangent] = 0.15;
        let mut flow = field(gas(0.0), grid, channel_boundaries(axis), |_| {
            [1.0, velocity[0], velocity[1], 2.2]
        });
        let initial = flow.totals();
        let mut op =
            MolecularTransport2d::new(&flow, law(&flow, mu, 0.003, conductivity), || false)
                .unwrap()
                .with_wall(
                    side,
                    no_slip(
                        axis,
                        0.4,
                        WallThermalCondition::Isothermal(Temperature::new(3.0)),
                    ),
                )
                .unwrap();
        let report = op
            .advance(
                &mut flow,
                Time::new(duration),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        let distance = 0.5 * grid.spacing()[axis];
        // Outward wall exchange = -mu*(u_wall-u_cell)/distance on either end.
        // Wall energy is its exact tangential work plus outward Fourier heat.
        let expected_force = -mu * (0.4 - 0.15) / distance;
        let expected_heat = -conductivity * (3.0 - 2.2) / distance;
        let expected_power = 0.4 * expected_force + expected_heat;
        let area_time = grid.face_area(side) * duration;
        for impulse in op.boundary_impulses(side) {
            assert_eq!(impulse.mass, 0.0);
            close(impulse.momentum[axis] / area_time, 0.0, 1e-12);
            close(impulse.momentum[tangent] / area_time, expected_force, 5e-8);
            close(impulse.energy / area_time, expected_power, 5e-8);
        }
        assert!(flow.totals().momentum[tangent] > initial.momentum[tangent]);
        assert!(flow.totals().energy > initial.energy);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-13);
        }
        assert_eq!(flow.time().value(), 0.0);
    }
}

#[test]
fn complete_coupled_co_moving_walls_preserve_uniform_translation_on_both_axes() {
    for axis in 0..2 {
        let mut velocity = [0.0; 2];
        velocity[1 - axis] = 0.25;
        let flow = field(
            gas(0.0),
            CartesianGrid::new(3, 4, 0.2, 0.3, 0.7).unwrap(),
            channel_boundaries(axis),
            |_| [1.0, velocity[0], velocity[1], 2.2],
        );
        let original = flow.states().to_vec();
        let transport = law(&flow, 0.02, 0.003, 0.05);
        let wall = no_slip(axis, 0.25, WallThermalCondition::Adiabatic);
        let mut coupled = MolecularEuler2d::new(flow, transport, || false).unwrap();
        for side in &SIDES[2 * axis..2 * axis + 2] {
            coupled = coupled.with_wall(*side, wall).unwrap();
            assert_eq!(coupled.wall(*side), wall);
        }
        for _ in 0..3 {
            let report = coupled
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            for (current, initial) in coupled.flow().states().iter().zip(&original) {
                assert_eq!(
                    current.array().map(f64::to_bits),
                    initial.array().map(f64::to_bits)
                );
            }
            for residual in report.conservation_residual.array() {
                close(residual, 0.0, 3e-13);
            }
            assert_eq!(report.boundary_change.energy, 0.0);
            assert_eq!(report.boundary_change.momentum[1 - axis], 0.0);
        }
    }
}

#[test]
fn native_couette_flow_preserves_momentum_and_turns_wall_power_into_original_gas_heat() {
    let (mu, difference, height, area, duration) = (0.02, 0.3, 1.0, 0.7, 0.001);
    for axis in 0..2 {
        for boost in [0.0, 3.0] {
            let tangent = 1 - axis;
            let flow = field(
                gas(0.0),
                CartesianGrid::new(4, 4, 0.25, 0.25, 0.7).unwrap(),
                channel_boundaries(axis),
                |x| {
                    let mut velocity = [0.0; 2];
                    velocity[tangent] = boost + difference * x[axis] / height;
                    [1.0, velocity[0], velocity[1], 2.3]
                },
            );
            let initial = flow.totals();
            let velocities: Vec<_> = flow.states().iter().map(|q| q.momentum).collect();
            let transport = law(&flow, mu, 0.003, 0.0);
            let mut coupled = MolecularEuler2d::new(flow, transport, || false).unwrap();
            for (side, speed) in [
                (SIDES[2 * axis], boost),
                (SIDES[2 * axis + 1], boost + difference),
            ] {
                coupled = coupled
                    .with_wall(side, no_slip(axis, speed, WallThermalCondition::Adiabatic))
                    .unwrap();
            }
            let expected_power = mu * difference * difference / height * area;
            for _ in 0..8 {
                let report = coupled
                    .step(Time::new(duration), MolecularStepControl::default(), || {
                        false
                    })
                    .unwrap();
                close(report.boundary_change.momentum[tangent], 0.0, 2e-13);
                close(
                    report.boundary_change.energy,
                    expected_power * report.elapsed.value(),
                    2e-13,
                );
                for (side, speed) in [
                    (SIDES[2 * axis], boost),
                    (SIDES[2 * axis + 1], boost + difference),
                ] {
                    for impulse in coupled.boundary_impulses(side) {
                        close(impulse.energy, speed * impulse.momentum[tangent], 2e-13);
                    }
                }
            }
            let elapsed = coupled.flow().time().value();
            close(coupled.flow().totals().mass, initial.mass, 2e-13);
            close(
                coupled.flow().totals().momentum[tangent],
                initial.momentum[tangent],
                2e-12,
            );
            close(
                coupled.flow().totals().energy - initial.energy,
                expected_power * elapsed,
                2e-11,
            );
            // Viscous momentum divergence vanishes for this linear profile;
            // unchanged kinetic energy makes the admitted wall work internal heat.
            let expected_temperature = 2.3 + mu * (difference / height).powi(2) * elapsed / 2.5;
            for (state, original) in coupled.flow().states().iter().zip(&velocities) {
                for component in 0..2 {
                    close(state.momentum[component], original[component], 2e-12);
                }
                close(
                    state
                        .primitive(coupled.flow().gas(), || false)
                        .unwrap()
                        .temperature,
                    expected_temperature,
                    2e-11,
                );
            }
        }
    }
}

/// The compatible wall traction keeps the adjacent cell's tangential divergence in its normal
/// component: `stress.rs::no_slip_traction` states that both cross terms belong to the discrete
/// adjoint, even for constant viscosity, and that replacing them by the derivatives of a constant
/// boundary trace breaks the energy balance. The wall's own tangential derivative is zero, so the
/// term is the compatible operator's limitation, kept visible here (it approaches the physical
/// trace under refinement). The expected value is the formula, `-sign * (bulk - 2/3 mu) * d_t u_t`
/// of the adjacent cell, from the fixture's own numbers.
#[test]
fn constant_no_slip_wall_normal_traction_carries_the_adjacent_cells_tangential_divergence() {
    for side in SIDES {
        let axis = side.axis();
        let tangent = 1 - axis;
        let grid = CartesianGrid::new(3, 3, 0.25, 0.25, 0.7).unwrap();
        let mut flow = field(gas(0.0), grid, channel_boundaries(axis), |x| {
            let distance = if side.sign() > 0.0 {
                0.75 - x[axis]
            } else {
                x[axis]
            };
            let mut velocity = [0.0; 2];
            velocity[tangent] = 0.05 * x[tangent] * distance;
            [1.0, velocity[0], velocity[1], 2.0]
        });
        let (mu, bulk) = (0.02, 0.003);
        let mut op = MolecularTransport2d::new(&flow, law(&flow, mu, bulk, 0.0), || false)
            .unwrap()
            .with_wall(side, no_slip(axis, 0.0, WallThermalCondition::Adiabatic))
            .unwrap();
        let duration = 1e-9;
        op.advance(
            &mut flow,
            Time::new(duration),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
        let middle = op.boundary_impulses(side)[1];
        // The adjacent cell's d_t u_t = 0.05 * (its distance from the wall) is nonzero while the
        // constant wall's is zero; the compatible traction carries (bulk - 2/3 mu) times it. The
        // 1e-10 allowance is the O(dt) normal velocity the wall traction induces during the step
        // (about 1e-13 here), far below the term itself (6.5e-5).
        let adjacent_divergence = 0.05 * (0.5 * grid.spacing[axis]);
        let expected = -side.sign() * (bulk - (2.0 / 3.0) * mu) * adjacent_divergence;
        let traction = middle.momentum[axis] / (duration * grid.face_area(side));
        assert!(
            (traction - expected).abs() < 1e-10,
            "{traction:e} against the compatible formula {expected:e}"
        );
        assert!(middle.momentum[tangent] > 0.0);
        assert_eq!(middle.energy, 0.0);
    }
}

#[test]
fn planar_wall_evolution_matches_the_same_native_three_dimensional_extrusion() {
    wall_extrusion(None);
}

#[test]
fn navier_slip_and_interface_heat_match_the_native_three_dimensional_extrusion() {
    wall_extrusion(Some((0.17, 0.35)));
}

fn wall_extrusion(slip: Option<(f64, f64)>) {
    let grid = CartesianGrid::new(4, 3, 0.25, 0.5, 0.7).unwrap();
    let mut planar = field(gas(0.07), grid, channel_boundaries(1), |x| {
        [1.0, 0.1 + 0.05 * x[1], 0.0, 2.2 + 0.1 * x[0]]
    });
    let mut spatial = Euler3d::new(
        planar.gas().clone(),
        Grid3d::new(4, 3, 1, 0.25, 0.5, 0.7).unwrap(),
        planar
            .states()
            .iter()
            .map(|q| State3d {
                density: q.density,
                momentum: [q.momentum[0], q.momentum[1], 0.0],
                energy: q.energy,
            })
            .collect(),
        Boundaries3d {
            x: AxisBoundaries::Periodic,
            y: WALLS,
            z: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    let mut planar_op =
        MolecularTransport2d::new(&planar, law(&planar, 0.02, 0.003, 0.05), || false).unwrap();
    let mut spatial_op =
        MolecularTransport3d::new(&spatial, law(&planar, 0.02, 0.003, 0.05), || false).unwrap();
    for (side2, side3, speed, thermal) in [
        (
            BoundarySide::YMin,
            Side3d::YMin,
            -0.1,
            WallThermalCondition::Adiabatic,
        ),
        (
            BoundarySide::YMax,
            Side3d::YMax,
            0.3,
            WallThermalCondition::Isothermal(Temperature::new(3.2)),
        ),
    ] {
        planar_op = planar_op
            .with_wall(
                side2,
                match slip {
                    Some((slip_length_m, fluid_heat_fraction)) => MolecularWall2d::NavierSlip {
                        velocity: [speed, 0.0],
                        thermal,
                        slip_length_m,
                        fluid_heat_fraction,
                    },
                    None => no_slip(1, speed, thermal),
                },
            )
            .unwrap();
        spatial_op = spatial_op
            .with_wall(
                side3,
                match slip {
                    Some((slip_length_m, fluid_heat_fraction)) => MolecularWall3d::NavierSlip {
                        velocity: [speed, 0.0, 0.0],
                        thermal,
                        slip_length_m,
                        fluid_heat_fraction,
                    },
                    None => MolecularWall3d::NoSlip {
                        velocity: [speed, 0.0, 0.0],
                        thermal,
                    },
                },
            )
            .unwrap();
    }
    let planar_report = planar_op
        .advance(
            &mut planar,
            Time::new(0.001),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let spatial_report = spatial_op
        .advance(
            &mut spatial,
            Time::new(0.001),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(planar_report.accepted_steps, spatial_report.accepted_steps);
    for (a, b) in planar.states().iter().zip(spatial.states()) {
        close(a.density, b.density, 2e-13);
        close(a.momentum[0], b.momentum[0], 2e-13);
        close(a.momentum[1], b.momentum[1], 2e-13);
        assert_eq!(b.momentum[2], 0.0);
        close(a.energy, b.energy, 2e-13);
    }
    for (side2, side3) in [
        (BoundarySide::YMin, Side3d::YMin),
        (BoundarySide::YMax, Side3d::YMax),
    ] {
        for (a, b) in planar_op
            .boundary_impulses(side2)
            .iter()
            .zip(spatial_op.boundary_impulses(side3))
        {
            close(a.momentum[0], b.momentum[0], 2e-13);
            close(a.momentum[1], b.momentum[1], 2e-13);
            close(a.energy, b.energy, 2e-13);
        }
    }
}

fn coupled_channel() -> MolecularEuler2d<ConstantTransport> {
    let flow = field(
        gas(0.0),
        CartesianGrid::new(3, 4, 0.2, 0.3, 0.7).unwrap(),
        channel_boundaries(1),
        |_| [1.0, 0.1, 0.0, 2.2],
    );
    let transport = law(&flow, 0.02, 0.003, 0.05);
    MolecularEuler2d::new(flow, transport, || false)
        .unwrap()
        .with_wall(
            BoundarySide::YMin,
            no_slip(1, 0.0, WallThermalCondition::Adiabatic),
        )
        .unwrap()
        .with_wall(
            BoundarySide::YMax,
            no_slip(
                1,
                0.3,
                WallThermalCondition::Isothermal(Temperature::new(3.0)),
            ),
        )
        .unwrap()
}

#[test]
fn wall_history_and_late_cancellation_preserve_the_original_complete_interval() {
    let mut reference = coupled_channel();
    let mut polls = 0;
    let accepted = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut replay = coupled_channel();
    let original = replay.flow().states().to_vec();
    let mut calls = 0;
    assert!(matches!(
        replay.step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            calls == polls - 1
        }),
        Err(MolecularError::Cancelled)
    ));
    assert!(calls > polls / 2);
    assert_eq!(replay.flow().states(), original);
    assert_eq!(replay.flow().time().value(), 0.0);
    assert!(replay.last_step().is_none());
    let result = replay
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(result, accepted);
    assert_eq!(replay.flow().states(), reference.flow().states());
    for side in SIDES {
        assert_eq!(
            replay.boundary_impulses(side),
            reference.boundary_impulses(side)
        );
    }
    let side = BoundarySide::YMax;
    let same = replay.wall(side);
    let replay = replay.with_wall(side, same).unwrap();
    assert_eq!(replay.last_step(), Some(&accepted));
    let previous = replay.flow().states().to_vec();
    let changed = replay
        .with_wall(side, no_slip(1, 0.3, WallThermalCondition::Adiabatic))
        .unwrap();
    assert!(changed.last_step().is_none());
    assert_eq!(changed.flow().states(), previous);
    assert_eq!(changed.flow().time(), reference.flow().time());
}

#[test]
fn wall_admission_and_source_face_temperature_refuse_without_partial_flow() {
    for wall in [
        MolecularWall2d::NoSlip {
            velocity: [0.0, 0.1],
            thermal: WallThermalCondition::Adiabatic,
        },
        no_slip(1, f64::NAN, WallThermalCondition::Adiabatic),
        no_slip(
            1,
            0.0,
            WallThermalCondition::Isothermal(Temperature::new(60.0)),
        ),
    ] {
        assert!(
            coupled_channel()
                .with_wall(BoundarySide::YMax, wall)
                .is_err()
        );
    }
    assert!(
        coupled_channel()
            .with_wall(BoundarySide::XMin, MolecularWall2d::SlipAdiabatic)
            .is_err()
    );
    let mut flow = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 0.25, 0.25, 0.7).unwrap(),
        Boundaries2d { x: WALLS, y: OPEN },
        |_| [1.0, 0.0, 0.0, 2.0],
    );
    assert!(
        MolecularTransport2d::new(&flow, law(&flow, 0.02, 0.003, 0.05), || false)
            .unwrap()
            .with_wall(
                BoundarySide::YMin,
                no_slip(1, 0.0, WallThermalCondition::Adiabatic)
            )
            .is_err()
    );
    let transport = ConstantTransport::new(
        bank(flow.gas()),
        "wall-face-temperature-admission",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 2.5,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.02),
            bulk_viscosity: DynViscosity::new(0.003),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.05),
        },
        vec![0.0],
    )
    .unwrap();
    let mut op = MolecularTransport2d::new(&flow, transport, || false)
        .unwrap()
        .with_wall(
            BoundarySide::XMax,
            no_slip(
                0,
                0.0,
                WallThermalCondition::Isothermal(Temperature::new(4.0)),
            ),
        )
        .unwrap();
    let original = flow.states().to_vec();
    let boundary = op.boundary_impulses(BoundarySide::XMax).to_vec();
    // Cell T=2 is admitted, wall T=4 is inside the gas domain, but the actual
    // face mean T=3 exceeds this source law's declared range of validity.
    assert!(matches!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            MolecularControl::default(),
            || false
        ),
        Err(MolecularError::Law(_))
    ));
    assert_eq!(flow.states(), original);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(op.boundary_impulses(BoundarySide::XMax), boundary);
}
