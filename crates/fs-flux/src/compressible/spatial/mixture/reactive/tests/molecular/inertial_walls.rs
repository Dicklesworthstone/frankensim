use super::*;
mod body_forcing;
use crate::compressible::spatial::mixture::molecular::{MolecularWall3d, WallThermalCondition};
use crate::compressible::spatial::mixture::reactive::{
    InertialWallLoad, InertialWallsControl, Piston, ReactiveInertialWallsError,
    ReactiveInertialWallsStep,
};

fn controls() -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-2.0, 2.0]; 2],
        momentum_tolerance_n_s: 3e-10,
        energy_tolerance_j: 3e-9,
        maximum_iterations: 32,
    }
}

fn walls(velocity: [f64; 2]) -> [Piston; 2] {
    [
        Piston::new(0.3, velocity[0]).unwrap(),
        Piston::new(0.6, velocity[1]).unwrap(),
    ]
}

fn channel(axis: usize) -> Boundaries3d {
    let mut result = PERIODIC;
    match axis {
        0 => result.x = WALLS,
        1 => result.y = WALLS,
        2 => result.z = WALLS,
        _ => unreachable!(),
    }
    result
}

fn case(axis: usize, heat: f64, rate: f64) -> ReactingEuler3d {
    let (gas, mechanism) = mechanism(heat, rate, 20.0);
    let grid = CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap();
    coupled(
        field(gas, grid, channel(axis), |_| ([1.0, 0.0], [0.0; 3], 1.0)),
        mechanism,
    )
}

fn check(report: &ReactiveInertialWallsStep, gas: &MixtureGas, controls: InertialWallsControl) {
    let mechanics = &report.mechanics;
    let reactive = &mechanics.gas.reactive;
    balance(reactive);
    assert_eq!(reactive.elapsed, reactive.transport.elapsed);
    for residual in mechanics
        .momentum_residual_n_s
        .into_iter()
        .chain([mechanics.combined_momentum_residual_n_s])
    {
        assert!(residual.abs() <= controls.momentum_tolerance_n_s);
    }
    for residual in mechanics
        .mechanical_energy_residual_j
        .into_iter()
        .chain([mechanics.total_energy_residual_j])
    {
        assert!(residual.abs() <= controls.energy_tolerance_j);
    }
    let heat: f64 = mechanics.gas.wall_heat_j.iter().sum();
    let body_work: f64 = (0..2)
        .map(|end| {
            mechanics.kinetic_energy_change_j[end] + mechanics.exterior_pressure_work_j[end]
                - mechanics.external_force_work_j[end]
        })
        .sum();
    let delta = reactive.after.energy - reactive.before.energy;
    let chemical: f64 = (0..gas.species_count())
        .map(|k| {
            (reactive.after.species_masses[k] - reactive.before.species_masses[k])
                * gas.energy_references()[k]
        })
        .sum();
    assert!((delta + chemical + body_work - heat).abs() <= controls.energy_tolerance_j);
    assert!(
        (delta - reactive.chemistry.referenced_energy_change_j + body_work - heat).abs()
            <= controls.energy_tolerance_j
    );
    for axis in 0..3 {
        if axis != report.axis {
            close(
                reactive.after.momentum[axis] - reactive.before.momentum[axis],
                report.side_boundary_impulse_n_s[axis]
                    + report.guide_impulse_n_s[0][axis]
                    + report.guide_impulse_n_s[1][axis],
                3e-11,
            );
        }
    }
}

#[test]
fn g1_reaction_pressure_accelerates_both_finite_masses_beyond_thermoneutral_gas() {
    for axis in 0..3 {
        let mut reacting = case(axis, 24.0, 30.0);
        let mut neutral = case(axis, 0.0, 30.0);
        let mut hot_walls = walls([0.0; 2]);
        let mut cold_walls = hot_walls;
        for _ in 0..3 {
            let hot = reacting
                .step_with_inertial_walls(
                    axis,
                    &mut hot_walls,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    [InertialWallLoad::default(); 2],
                    controls(),
                    || false,
                )
                .unwrap();
            let cold = neutral
                .step_with_inertial_walls(
                    axis,
                    &mut cold_walls,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    [InertialWallLoad::default(); 2],
                    controls(),
                    || false,
                )
                .unwrap();
            check(&hot, reacting.flow().gas(), controls());
            check(&cold, neutral.flow().gas(), controls());
            assert_eq!(
                hot.mechanics.gas.reactive.elapsed,
                cold.mechanics.gas.reactive.elapsed
            );
            assert!(hot.mechanics.gas.reactive.chemistry.species_change_kg[0] < 0.0);
            assert!(hot.mechanics.gas.reactive.chemistry.species_change_kg[1] > 0.0);
            assert!(
                hot.mechanics
                    .gas
                    .reactive
                    .chemistry
                    .referenced_energy_change_j
                    > 0.0
            );
            assert_eq!(
                hot.mechanics.gas.reactive.chemistry.cell_solves,
                2 * reacting.flow().cell_count()
            );
        }
        assert!(hot_walls[0].velocity_m_s() < 1.05 * cold_walls[0].velocity_m_s());
        assert!(hot_walls[1].velocity_m_s() > 1.05 * cold_walls[1].velocity_m_s());
        assert!(reacting.flow().totals().species_masses[1] > 0.01);
    }
}

#[test]
fn g2_reactive_compression_and_expansion_preserve_species_and_reciprocal_momentum() {
    for axis in 0..3 {
        for velocity in [[-0.12, 0.08], [0.12, -0.08]] {
            let mut flow = case(axis, 16.0, 12.0);
            let original_grid = *flow.flow().grid();
            let mut pistons = walls(velocity);
            let original = pistons;
            let report = flow
                .step_with_inertial_walls(
                    axis,
                    &mut pistons,
                    Time::new(0.001),
                    ReactiveStepControl::default(),
                    [InertialWallLoad::default(); 2],
                    controls(),
                    || false,
                )
                .unwrap();
            check(&report, flow.flow().gas(), controls());
            let r = &report.mechanics.gas.reactive;
            assert!(r.chemistry.species_change_kg[1] > 0.0);
            close(r.after.mass(), r.before.mass(), 3e-12);
            for k in 0..2 {
                close(
                    r.after.species_masses[k] - r.before.species_masses[k],
                    r.chemistry.species_change_kg[k],
                    3e-12,
                );
            }
            let body_impulse: f64 = (0..2)
                .map(|end| {
                    pistons[end].mass_kg()
                        * (pistons[end].velocity_m_s() - original[end].velocity_m_s())
                })
                .sum();
            close(
                r.after.momentum[axis] - r.before.momentum[axis] + body_impulse,
                0.0,
                4e-10,
            );
            let spacing_change = flow.flow().grid().spacing()[axis] - original_grid.spacing()[axis];
            assert_eq!(
                spacing_change.is_sign_positive(),
                (velocity[1] - velocity[0]).is_sign_positive()
            );
            assert_eq!(flow.last_step(), Some(r));
        }
    }
}

fn charted(minimum: f64) -> ReactingEuler3d {
    // Unequal Cv makes changing Tmin a DIFFERENT reference shift per species,
    // while the exact underlying NASA bank and physical reaction remain fixed.
    let source_gas = make_gas(vec![
        source("A", R, 3.5, 0.0, 1.0e6),
        source("B", R, 4.5, 0.0, 1.0e6 - 24.0),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("reactive-piston-chart").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(12.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        source_gas.source().clone(),
        "same-physical-piston-mechanism-two-charts",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let gas = MixtureGas::new(
        source_gas.source().clone(),
        fs_gas::dynamics::CaloricDomain {
            minimum_temperature: minimum,
            ..source_gas.domain()
        },
        || false,
    )
    .unwrap();
    let flow = field(
        gas,
        CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap(),
        channel(0),
        |_| ([0.9, 0.1], [0.03, 0.0, 0.0], 2.0),
    );
    let chemistry = IsochoricKineticsSolve {
        minimum_temperature_k: 0.5,
        energy_tolerance_j: 1e-9,
        ..chemistry()
    };
    let reactor = CellChemistry::new(flow.gas().clone(), mechanism, chemistry, || false).unwrap();
    let bytes = ReactingEuler3d::required_workspace_bytes(&flow, &reactor).unwrap();
    ReactingEuler3d::with_chemistry(flow, reactor, bytes, || false).unwrap()
}

#[test]
fn g3_large_unequal_species_reference_charts_keep_actual_reaction_and_piston_dynamics() {
    let mut a = charted(0.001);
    let mut b = charted(0.5);
    assert_eq!(a.mechanism(), b.mechanism());
    let shifts: Vec<_> = a
        .flow()
        .gas()
        .energy_references()
        .iter()
        .zip(b.flow().gas().energy_references())
        .map(|(a, b)| b - a)
        .collect();
    assert!((shifts[0] - shifts[1]).abs() > 0.1);
    let mut wa = walls([-0.03, 0.05]);
    let mut wb = wa;
    let c = InertialWallsControl {
        energy_tolerance_j: 3e-8,
        ..controls()
    };
    let ra = a
        .step_with_inertial_walls(
            0,
            &mut wa,
            Time::new(0.001),
            ReactiveStepControl::default(),
            [InertialWallLoad::default(); 2],
            c,
            || false,
        )
        .unwrap();
    let rb = b
        .step_with_inertial_walls(
            0,
            &mut wb,
            Time::new(0.001),
            ReactiveStepControl::default(),
            [InertialWallLoad::default(); 2],
            c,
            || false,
        )
        .unwrap();
    for report in [&ra, &rb] {
        assert!(report.mechanics.gas.reactive.chemistry.species_change_kg[1] > 1e-4);
        assert!(
            report
                .mechanics
                .gas
                .reactive
                .chemistry
                .referenced_energy_change_j
                .abs()
                > 1e-3
        );
        assert!(report.mechanics.total_energy_residual_j.abs() <= c.energy_tolerance_j);
    }
    assert_eq!(
        ra.mechanics.gas.reactive.elapsed,
        rb.mechanics.gas.reactive.elapsed
    );
    for end in 0..2 {
        close(wa[end].velocity_m_s(), wb[end].velocity_m_s(), 3e-10);
    }
    for cell in 0..a.flow().cell_count() {
        let qa = a.flow().cell(cell).unwrap();
        let qb = b.flow().cell(cell).unwrap();
        let pa = qa.primitive(a.flow().gas(), || false).unwrap();
        let pb = qb.primitive(b.flow().gas(), || false).unwrap();
        for k in 0..2 {
            close(qa.partial_densities[k], qb.partial_densities[k], 3e-10);
        }
        close(pa.temperature, pb.temperature, 2e-8);
        close(pa.pressure, pb.pressure, 2e-8);
    }
}

fn viscous_case() -> (ReactingEuler3d, MolecularTransport3d<ConstantTransport>) {
    let (gas, mechanism) = mechanism(16.0, 10.0, 100.0);
    let flow = coupled(
        field(
            gas,
            CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap(),
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let first = 0.75 + 0.03 * x[0];
                (
                    [first, 1.0 - first],
                    [0.18 + 0.03 * x[1], 0.1 + 0.02 * x[2], -0.06],
                    1.0 + 0.15 * x[1] + 0.2 * x[2],
                )
            },
        ),
        mechanism,
    );
    let mut op = operator(flow.flow(), 0.015);
    for side in BoundarySide::ALL {
        let model = if side == BoundarySide::XMax {
            MolecularWall3d::NavierSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
                slip_length_m: 0.08,
                fluid_heat_fraction: 0.4,
            }
        } else {
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
            }
        };
        op = op.with_wall(side, model).unwrap();
    }
    (flow, op)
}

#[test]
fn g2_reacting_viscous_pistons_separate_wall_heat_and_retain_side_and_guide_loads() {
    let (mut flow, mut op) = viscous_case();
    let mut pistons = walls([-0.02, 0.03]);
    let bytes = budget(&flow, &op);
    let report = flow
        .step_with_dissipation_and_inertial_walls(
            &mut op,
            0,
            &mut pistons,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            [InertialWallLoad::default(); 2],
            controls(),
            || false,
        )
        .unwrap();
    check(&report, flow.flow().gas(), controls());
    let gas = &report.mechanics.gas;
    assert!(gas.reactive.chemistry.species_change_kg[1] > 0.0);
    assert!(gas.reactive.chemistry.referenced_energy_change_j > 0.0);
    assert!(gas.wall_heat_j.iter().sum::<f64>() > 1e-5);
    assert!(report.side_boundary_impulse_n_s[0] < -1e-6);
    assert!(
        report
            .guide_impulse_n_s
            .iter()
            .flatten()
            .any(|q| q.abs() > 1e-6)
    );
    assert!(
        report
            .gas_angular_impulse_n_m_s
            .iter()
            .flatten()
            .any(|q| q.abs() > 1e-6)
    );
    let split = gas.reactive.dissipation.unwrap();
    assert_eq!(op.last_step(), Some(&split.second));
    close(
        gas.wall_heat_j.iter().sum::<f64>() + gas.molecular_wall_work_j.iter().sum::<f64>(),
        split.first.boundary_change.energy + split.second.boundary_change.energy,
        3e-12,
    );
    let mut side = [0.0; 3];
    for wall in BoundarySide::ALL
        .into_iter()
        .filter(|side| side.axis() != 0)
    {
        for q in flow.boundary_impulses(wall) {
            for a in 0..3 {
                side[a] -= q.momentum[a];
            }
        }
    }
    for a in 0..3 {
        close(side[a], report.side_boundary_impulse_n_s[a], 3e-13);
    }
    // A subsequent stationary molecular/chemical split verifies that both
    // accepted owners now bind the actual moved geometry.
    flow.step_with_dissipation(
        &mut op,
        Time::new(0.0002),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        bytes,
        || false,
    )
    .unwrap();
}

#[test]
fn g4_reacting_piston_cancel_and_root_exhaustion_preserve_both_accepted_owners() {
    let (mut reference, mut reference_op) = viscous_case();
    let mut reference_walls = walls([-0.02, 0.03]);
    let bytes = budget(&reference, &reference_op);
    reference
        .step_with_dissipation_and_inertial_walls(
            &mut reference_op,
            0,
            &mut reference_walls,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            [InertialWallLoad::default(); 2],
            controls(),
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = reference
        .step_with_dissipation_and_inertial_walls(
            &mut reference_op,
            0,
            &mut reference_walls,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            [InertialWallLoad::default(); 2],
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [polls / 2, polls] {
        let (mut flow, mut op) = viscous_case();
        let mut pistons = walls([-0.02, 0.03]);
        flow.step_with_dissipation_and_inertial_walls(
            &mut op,
            0,
            &mut pistons,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            [InertialWallLoad::default(); 2],
            controls(),
            || false,
        )
        .unwrap();
        let state = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let previous = flow.last_step().copied();
        let operator_report = op.last_step().copied();
        let previous_walls = pistons;
        let wall_data = BoundarySide::ALL.map(|side| op.wall(side));
        let original_loads = loads(&flow);
        let original_molecular = molecular_loads(&op);
        let mut calls = 0;
        assert_eq!(
            flow.step_with_dissipation_and_inertial_walls(
                &mut op,
                0,
                &mut pistons,
                Time::new(0.001),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                [InertialWallLoad::default(); 2],
                controls(),
                || {
                    calls += 1;
                    calls == cutoff
                }
            ),
            Err(ReactiveInertialWallsError::Cancelled)
        );
        assert_eq!(flow.flow().states(), state);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), previous);
        assert_eq!(op.last_step().copied(), operator_report);
        assert_eq!(pistons, previous_walls);
        assert_eq!(loads(&flow), original_loads);
        assert_eq!(molecular_loads(&op), original_molecular);
        assert_eq!(BoundarySide::ALL.map(|side| op.wall(side)), wall_data);
        let recovered = flow
            .step_with_dissipation_and_inertial_walls(
                &mut op,
                0,
                &mut pistons,
                Time::new(0.001),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                [InertialWallLoad::default(); 2],
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(recovered, expected);
        assert_eq!(pistons, reference_walls);
        assert_eq!(flow.flow().states(), reference.flow().states());
    }
    let (mut flow, mut op) = viscous_case();
    let state = flow.flow().states().to_vec();
    let grid = *flow.flow().grid();
    let original_loads = loads(&flow);
    let original_molecular = molecular_loads(&op);
    let mut pistons = walls([0.0; 2]);
    let original = pistons;
    let refused = flow.step_with_dissipation_and_inertial_walls(
        &mut op,
        0,
        &mut pistons,
        Time::new(0.001),
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        DissipationControl::default(),
        bytes,
        [InertialWallLoad::default(); 2],
        InertialWallsControl {
            midpoint_velocity_bracket_m_s: [[-1e-9, 1e-9]; 2],
            maximum_iterations: 1,
            ..controls()
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(ReactiveInertialWallsError::StepRetries { attempts: 1, .. })
    ));
    assert_eq!(flow.flow().states(), state);
    assert_eq!(*flow.flow().grid(), grid);
    assert_eq!(pistons, original);
    assert_eq!(flow.last_step(), None);
    assert_eq!(op.last_step(), None);
    assert_eq!(loads(&flow), original_loads);
    assert_eq!(molecular_loads(&op), original_molecular);
}

#[test]
fn g2_fast_reaction_wave_shortening_retries_whole_finite_piston_interval() {
    let mut flow = case(0, 100.0, 500.0);
    let mut pistons = [
        Piston::new(5.0, 0.0).unwrap(),
        Piston::new(7.0, 0.0).unwrap(),
    ];
    let initial = pistons;
    let c = InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-0.5, 0.5]; 2],
        ..controls()
    };
    let report = flow
        .step_with_inertial_walls(
            0,
            &mut pistons,
            Time::new(0.02),
            ReactiveStepControl::default(),
            [InertialWallLoad::default(); 2],
            c,
            || false,
        )
        .unwrap();
    check(&report, flow.flow().gas(), c);
    assert!(report.mechanics.gas.reactive.coupling_retries > 0);
    assert!(report.mechanics.gas.reactive.elapsed.value() < 0.02);
    let mut direct = case(0, 100.0, 500.0);
    let mut direct_walls = initial;
    let recovered = direct
        .step_with_inertial_walls(
            0,
            &mut direct_walls,
            report.mechanics.gas.reactive.elapsed,
            ReactiveStepControl::default(),
            [InertialWallLoad::default(); 2],
            c,
            || false,
        )
        .unwrap();
    assert_eq!(flow.flow().states(), direct.flow().states());
    assert_eq!(flow.flow().grid(), direct.flow().grid());
    assert_eq!(pistons, direct_walls);
    assert_eq!(
        report.mechanics.gas.reactive.chemistry,
        recovered.mechanics.gas.reactive.chemistry
    );
}
