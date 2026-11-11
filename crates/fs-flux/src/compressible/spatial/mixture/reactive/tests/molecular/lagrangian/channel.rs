//! The original chemical and physical-wall owners share one pathline transaction.

use super::*;
use crate::compressible::spatial::mixture::molecular::{MolecularWall3d, WallThermalCondition};
use fs_lbm::sgs3::dynamic::channel::ScalarWall;

#[test]
fn g4_reacting_channel_history_retains_no_slip_heat_chemistry_and_atomic_replay() {
    fn noise(cell: usize, salt: u64) -> f64 {
        let mut h = ((cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
        h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        h ^= h >> 31;
        (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
    }
    let shape = [5, 4, 6];
    let spacing = [0.7, 1.1, 0.9];
    let walls = [[0.0, 0.0, 0.01], [0.0, 0.0, -0.02]];
    let temperatures = [2.2, 0.9];
    let make = || {
        let (gas, mechanism) = mechanism(5.0, 8.0, 0.0);
        let grid = CartesianGrid::new(
            shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
        )
        .unwrap();
        let flow = coupled(
            field(
                gas,
                grid,
                Boundaries3d {
                    y: WALLS,
                    ..PERIODIC
                },
                |point| {
                    let cell = grid.index(std::array::from_fn(|axis| {
                        (point[axis] / spacing[axis]).floor() as usize
                    }));
                    let eta = point[1] / (shape[1] as f64 * spacing[1]);
                    let envelope = 4.0 * eta * (1.0 - eta);
                    let density = 1.0 + 0.15 * noise(cell, 911);
                    let heat = 1.0 - eta + 0.05 * envelope * noise(cell, 101);
                    let first = 0.4 + 0.1 * heat + 0.02 * noise(cell, 291);
                    (
                        [density * first, density * (1.0 - first)],
                        std::array::from_fn(|axis| {
                            0.02 * envelope * noise(cell, 31 * axis as u64 + 7)
                        }),
                        density * (1.0 + heat),
                    )
                },
            ),
            mechanism,
        );
        let mut molecular = operator(flow.flow(), 0.02);
        let bytes = flow
            .required_lagrangian_dynamic_species_sgs_workspace_bytes(&molecular)
            .unwrap();
        let before = published(&flow, &molecular);
        assert_eq!(
            flow.attach_lagrangian_dynamic_channel_species_sgs(
                &mut molecular,
                config(),
                bytes,
                || false,
            ),
            Err(ReactiveError::Dissipation(DissipationError::Invalid(
                "Lagrangian reacting channel requires two configured no-slip molecular walls",
            )))
        );
        assert_eq!(published(&flow, &molecular), before);
        assert_eq!(molecular.lagrangian_channel_normal_axis(), None);
        for (side, velocity, temperature) in [
            (BoundarySide::YMin, walls[0], temperatures[0]),
            (BoundarySide::YMax, walls[1], temperatures[1]),
        ] {
            molecular = molecular
                .with_wall(
                    side,
                    MolecularWall3d::NoSlip {
                        velocity,
                        thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
                    },
                )
                .unwrap();
        }
        flow.attach_lagrangian_dynamic_channel_species_sgs(&mut molecular, config(), bytes, || {
            false
        })
        .unwrap();
        assert_eq!(molecular.lagrangian_channel_normal_axis(), Some(1));
        assert_eq!(budget(&flow, &molecular), bytes);
        (flow, molecular, bytes)
    };
    let (mut flow, mut molecular, bytes) = make();
    let (mut replay, mut replay_op, replay_bytes) = make();
    let mut history =
        LagrangianDynamicComponents::new_channel(shape, spacing, 1, 2, usize::MAX).unwrap();
    for duration in [0.001, 0.0005] {
        let before = published(&replay, &replay_op);
        let primitive: Vec<_> = (0..flow.flow().cell_count())
            .map(|i| {
                flow.flow()
                    .cell(i)
                    .unwrap()
                    .primitive(flow.flow().gas(), || false)
                    .unwrap()
            })
            .collect();
        let density: Vec<_> = primitive.iter().map(|p| p.density).collect();
        let velocity: Vec<_> = primitive.iter().map(|p| p.velocity).collect();
        let temperature: Vec<_> = primitive.iter().map(|p| p.temperature).collect();
        let trial = history
            .prepare_with_walls(
                config().momentum,
                config().heat,
                config().species,
                before.time.value(),
                &density,
                &velocity,
                walls,
                &temperature,
                temperatures.map(ScalarWall::Value),
                |component, cell| flow.flow().states()[6 * cell + component] / density[cell],
                || false,
            )
            .unwrap();
        let expected = (
            Some(trial.momentum_report()),
            Some(trial.scalar_report()),
            Some(trial.components_report()),
        );
        assert!(expected.2.unwrap().maximum_coefficient > 0.0);
        let mut polls = 0;
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(duration),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        balance(&report);
        assert_eq!(fits(&molecular), expected);
        assert_eq!(
            molecular.lagrangian_subgrid_sample_time(),
            Some(before.time)
        );
        assert_eq!(report.chemistry.cell_solves, 2 * flow.flow().cell_count());
        assert!(report.chemistry.species_change_kg[0] < 0.0);
        assert!(report.chemistry.referenced_energy_change_j > 0.0);
        assert!(report.boundary_change.energy.abs() > 1e-10);
        // The y walls are impermeable: their boundary rows are exactly 0.0. The periodic x and z
        // sides enter the outward sum as +f and -f for each face (bitwise equal pairs), which
        // cancel in real arithmetic but leave rounding in a floating-point sum of all the rows:
        // at most rows*eps*sum|f|. A face moves at most one cell's species mass per step (CFL < 1),
        // so sum|f| <= rows*m_cell with m_cell <= rho_max*V = 1.15*V (density is 1 + 0.15*noise).
        // A real wall flux would be of the face-flux size (about 1e-5), many orders above this.
        let rows = 2 * (shape[1] * shape[2] + shape[0] * shape[2] + shape[0] * shape[1]);
        let rounding = (rows * rows) as f64
            * f64::EPSILON
            * 1.15
            * (spacing[0] * spacing[1] * spacing[2]);
        assert!(
            report
                .boundary_change
                .species_masses
                .iter()
                .all(|mass| mass.abs() <= rounding),
            "species masses through the boundary {:?} exceed the rounding bound {rounding:e}",
            report.boundary_change.species_masses
        );
        trial.commit();

        // Cancel after the final reacting/molecular balance, immediately before
        // physical publication; all three old fits and wall loads survive.
        let mut calls = 0;
        assert_eq!(
            replay.step_with_dissipation(
                &mut replay_op,
                Time::new(duration),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                replay_bytes,
                || {
                    calls += 1;
                    calls == polls
                },
            ),
            Err(ReactiveError::Cancelled)
        );
        assert_eq!(published(&replay, &replay_op), before);
        assert_eq!(
            replay
                .step_with_dissipation(
                    &mut replay_op,
                    Time::new(duration),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    replay_bytes,
                    || false,
                )
                .unwrap(),
            report
        );
        assert_eq!(published(&replay, &replay_op), published(&flow, &molecular));
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            assert_eq!(replay_op.wall(side), molecular.wall(side));
        }
    }
}
