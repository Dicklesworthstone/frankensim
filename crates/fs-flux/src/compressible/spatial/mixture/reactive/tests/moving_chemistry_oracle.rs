//! Independent concentration-dependent reaction oracle on changing volumes.

use super::*;
use crate::compressible::spatial::MeshMotion;

#[test]
fn moving_bimolecular_source_uses_both_actual_volumes_on_every_axis() {
    // A + B -> 2 B has identical species calorics and molecular masses.
    // There is no chemical heat release, so temperature changes caused by
    // pressure work cannot alter this temperature-independent reaction rate.
    // Its *integrated* B mass obeys logistic growth at each fixed volume.
    // One impermeable cell retains both species masses through the ALE phase;
    // the two exact chemical half intervals therefore have a closed solution
    // involving BOTH physical volumes. A frozen volume gives a different
    // answer, even though total mass and chemical energy would still balance.
    let duration = 0.002;
    let mass_rate = 100.0;
    for axis in 0..3 {
        let mut converted = [0.0; 2];
        for (case, lower_speed) in [5.0, -5.0].into_iter().enumerate() {
            let gas = make_gas(vec![
                source("A", R, 3.5, 0.0, 0.0),
                source("B", R, 3.5, 0.0, 0.0),
            ]);
            let reaction = ElementaryReactionV1::new(
                ReactionId::new("moving-bimolecular-autocatalysis").unwrap(),
                vec![
                    (SpeciesId::new("A").unwrap(), 1),
                    (SpeciesId::new("B").unwrap(), 1),
                ],
                vec![(SpeciesId::new("B").unwrap(), 2)],
                ArrheniusRateV1::new(
                    Temperature::new(1.0),
                    ReactionProgressRateQuantityV1::new(mass_rate * R),
                    0.0,
                    MolarEnergyQuantityV1::new(0.0),
                )
                .unwrap(),
            )
            .unwrap()
            .irreversible();
            let kinetics = IdealGasKineticsV1::new(
                gas.source().clone(),
                "moving-bimolecular-volume-oracle",
                1,
                vec![reaction],
                1e-12,
            )
            .unwrap();
            let grid = CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
            let flow = field(
                gas,
                grid,
                Boundaries3d {
                    x: WALLS,
                    y: WALLS,
                    z: WALLS,
                },
                |_| ([0.8, 0.2], [0.0; 3], 1.0),
            );
            let source = CellChemistry::new(
                flow.gas().clone(),
                kinetics,
                IsochoricKineticsSolve {
                    relative_tolerance: 1e-9,
                    absolute_amount_tolerance_mol: 1e-15,
                    initial_step_s: 1e-5,
                    ..chemistry()
                },
                || false,
            )
            .unwrap();
            let bytes = ReactingEuler3d::required_workspace_bytes(&flow, &source).unwrap();
            let mut flow = ReactingEuler3d::with_chemistry(flow, source, bytes, || false).unwrap();
            let mut motion = MeshMotion::default();
            motion.lower_velocity[axis] = lower_speed;
            motion.upper_velocity[axis] = -lower_speed;
            let report = flow
                .step_moving(
                    Time::new(duration),
                    ReactiveStepControl {
                        maximum_coupling_retries: 0,
                        ..ReactiveStepControl::default()
                    },
                    motion,
                    || false,
                )
                .unwrap();
            assert_eq!(report.reactive.elapsed.value(), duration);
            let old_volume = grid.cell_volume();
            let new_volume = flow.flow().grid().cell_volume();
            close(new_volume, 1.0 - 2.0 * lower_speed * duration, 1e-14);
            let exposure = mass_rate * duration * 0.5 * (1.0 / old_volume + 1.0 / new_volume);
            let expected_b_mass = 1.0 / (1.0 + 4.0 * (-exposure).exp());
            let observed_b_mass = flow.flow().cell(0).unwrap().partial_densities[1] * new_volume;
            close(observed_b_mass, expected_b_mass, 2e-8);
            close(report.reactive.after.mass(), 1.0, 2e-13);
            close(
                report.reactive.chemistry.species_change_kg[1],
                expected_b_mass - 0.2,
                2e-8,
            );
            assert_eq!(report.reactive.chemistry.cell_solves, 2);
            assert!(report.reactive.chemistry.entropy_change_j_k > 0.0);
            close(
                report.reactive.chemistry.referenced_energy_change_j,
                0.0,
                1e-13,
            );
            balance(&report.reactive);
            converted[case] = observed_b_mass;
        }
        // Compression raises concentrations and conversion; expansion lowers
        // both. The gap is well above the declared native chemistry error.
        assert!(converted[0] - converted[1] > 5e-4);
    }
}
