//! G1/G3/G4 two finite baths and the original native diffusion field.
//! Source coefficients below are authored numerical fixtures, not calibration.

use super::*;

#[path = "paired/adaptive.rs"]
mod adaptive;
use fs_blake3::ContentHash;
use fs_couple::moisture::{
    cartesian::{MoistureCartesian, MoistureCartesianGeometry},
    radial::{MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry},
};
use fs_material::osmosis::solubility::kinetics::KineticSolubilitySolution;
use fs_qty::{Time, VolumetricFlowRate};
use fs_thermochem::solution::solubility::kinetics::{
    PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1,
};

fn pair_controls(volumes: [f64; 2]) -> PairedSolutionMoistureSolve {
    PairedSolutionMoistureSolve {
        solutions: volumes.map(controls),
        maximum_evaluations: 20_000,
        mass_tolerance_kg: 1e-12,
    }
}

fn pair_faces(slope: f64) -> [Option<SolutionMoistureFace>; 2] {
    [MoistureSolutionPort::First, MoistureSolutionPort::Second].map(|solution| {
        Some(SolutionMoistureFace {
            solution,
            sorption: atmosphere::sorption(slope),
        })
    })
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

#[test]
fn g1_distinct_finite_baths_retain_through_flow_with_zero_net_solid_uptake() {
    with_cx(|cx| {
        let activity = [0.997_f64, 0.992];
        let final_volume = [5e-5, 8e-5];
        let slope = 0.2;
        let diffusion = 0.01;
        let duration = 2.0;
        let m0 = slope * (activity[0] + activity[1]) * 0.5;
        let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![m0], cx);
        // Two actual half-cell resistances. The original stationary cell can
        // still pass water from one finite solution to the other.
        let transfer =
            activity.map(|a| initial.dry_mass_kg() * 2.0 * diffusion * duration * (slope * a - m0));
        let volumes: [f64; 2] =
            std::array::from_fn(|i| final_volume[i] + transfer[i] / WATER_DENSITY);
        let first_total = -activity[0].ln() * final_volume[0] / WATER_VOLUME;
        // The second bath has its own nonideal free-energy model. Invert its
        // independent solvent law -ln(a_w)/v = c + B*c²/2 at the target volume.
        let interaction = 0.001;
        let osmotic_concentration = -activity[1].ln() / WATER_VOLUME;
        let concentration = 2.0 * osmotic_concentration
            / (1.0 + (1.0 + 2.0 * interaction * osmotic_concentration).sqrt());
        let model = Arc::new(
            VirialSolutionModel::new(293.15, WATER_VOLUME, 1.0, vec![3000.0], vec![interaction])
                .unwrap(),
        );
        let baths = [
            reservoir(volumes[0], vec![0.3 * first_total, 0.7 * first_total]),
            FiniteMoistureSolution::new(
                AqueousSolution::nonideal(
                    model.clone(),
                    volumes[1],
                    vec![concentration * final_volume[1]],
                )
                .unwrap(),
                WATER_MOLAR_MASS,
            )
            .unwrap(),
        ];
        let faces = pair_faces(slope);
        let result = baths[0]
            .advance_paired(
                &baths[1],
                &initial,
                duration,
                &faces,
                config(),
                pair_controls(volumes),
                cx,
            )
            .unwrap();
        near(result.moisture().moisture()[0], m0, 2e-13);
        near(result.step().water_change_kg, 0.0, 2e-13);
        assert!(transfer[0] > 0.0 && transfer[1] < 0.0);
        assert!(result.step().field_evaluations > 1);
        for i in 0..2 {
            let bath = &result.solutions()[i];
            near(bath.solution().solvent_volume_m3(), final_volume[i], 2e-16);
            near(
                bath.solution().water_activity().unwrap(),
                activity[i],
                2e-12,
            );
            near(
                result.step().solutions[i].water_to_solid_kg,
                transfer[i],
                2e-13,
            );
            near(
                result.step().transport.boundary_water_kg[i],
                transfer[i],
                2e-13,
            );
            near(
                result.step().solutions[i].water_balance_residual_kg,
                0.0,
                1e-13,
            );
            assert_eq!(
                bath.solution().solute_moles(),
                baths[i].solution().solute_moles()
            );
            let source = result.step().transport.surface_sorption[i]
                .as_ref()
                .unwrap();
            near(
                source.query().axes()["relative_humidity"],
                activity[i],
                2e-12,
            );
            assert_eq!(baths[i].solution().solvent_volume_m3(), volumes[i]);
        }
        near(result.step().field_water_balance_residual_kg, 0.0, 1e-12);
        near(result.step().water_balance_residual_kg, 0.0, 1e-12);
        near(result.step().direct_water_balance_residual_kg, 0.0, 1e-12);
        assert!(std::ptr::eq(
            result.solutions()[1].solution().nonideal_model().unwrap(),
            model.as_ref()
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[m0]);

        // Rename the two ports without changing which physical bath touches
        // each face. Root nesting may differ; the physical endpoint must agree.
        let swapped_faces = [faces[1].clone(), faces[0].clone()];
        let swapped = baths[1]
            .advance_paired(
                &baths[0],
                &initial,
                duration,
                &swapped_faces,
                config(),
                pair_controls([volumes[1], volumes[0]]),
                cx,
            )
            .unwrap();
        near(
            swapped.moisture().moisture()[0],
            result.moisture().moisture()[0],
            3e-13,
        );
        for i in 0..2 {
            near(
                swapped.solutions()[1 - i].solution().solvent_volume_m3(),
                result.solutions()[i].solution().solvent_volume_m3(),
                3e-16,
            );
        }
    });
}

#[test]
fn g3_identical_baths_match_the_original_single_combined_inventory() {
    with_cx(|cx| {
        let volume = 1e-4;
        let each = reservoir(volume, vec![0.02]);
        let combined = reservoir(2.0 * volume, vec![0.04]);
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.17, 0.16, 0.17], cx);
        let law = atmosphere::sorption(0.2);
        let original = combined
            .advance(
                &initial,
                0.02,
                &[Some(law.clone()), Some(law)],
                config(),
                controls(2.0 * volume),
                cx,
            )
            .unwrap();
        let paired = each
            .advance_paired(
                &each,
                &initial,
                0.02,
                &pair_faces(0.2),
                config(),
                pair_controls([volume; 2]),
                cx,
            )
            .unwrap();
        for (&a, &b) in paired
            .moisture()
            .moisture()
            .iter()
            .zip(original.moisture().moisture())
        {
            near(a, b, 2e-12);
        }
        for bath in paired.solutions() {
            near(
                2.0 * bath.solution().solvent_volume_m3(),
                original.reservoir().solution().solvent_volume_m3(),
                3e-16,
            );
            assert_eq!(
                bath.solution().solute_moles(),
                each.solution().solute_moles()
            );
        }
        near(
            paired.step().solutions[0].water_to_solid_kg,
            paired.step().solutions[1].water_to_solid_kg,
            3e-13,
        );
    });
}

#[test]
fn g1_each_saturated_bath_keeps_its_own_finite_crystals_during_through_flow() {
    with_cx(|cx| {
        let saturation = [500.0, 1500.0];
        let activity = saturation.map(|c| (-WATER_VOLUME * c).exp());
        let final_volume: [f64; 2] = [6e-5, 8e-5];
        let slope = 0.2;
        let m0 = slope * (activity[0] + activity[1]) * 0.5;
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![m0], cx);
        let transfer = activity.map(|a| initial.dry_mass_kg() * 0.02 * (slope * a - m0));
        let volumes: [f64; 2] =
            std::array::from_fn(|i| final_volume[i] + transfer[i] / WATER_DENSITY);
        let totals: [f64; 2] =
            std::array::from_fn(|i| 2.0 * saturation[i] * volumes[i].max(final_volume[i]));
        let baths: [FiniteMoistureSolution; 2] =
            std::array::from_fn(|i| soluble_reservoir(volumes[i], totals[i], saturation[i]));
        let result = baths[0]
            .advance_paired(
                &baths[1],
                &initial,
                1.0,
                &pair_faces(slope),
                config(),
                pair_controls(volumes),
                cx,
            )
            .unwrap();
        near(result.moisture().moisture()[0], m0, 3e-13);
        for i in 0..2 {
            let before = baths[i].solubility().unwrap();
            let after = result.solutions()[i].solubility().unwrap();
            let phase = after.phase_partitions()[0];
            assert_eq!(before.components(), after.components());
            near(
                phase.dissolved.value(),
                saturation[i] * final_volume[i],
                3e-13,
            );
            near(
                phase.solid.value(),
                totals[i] - saturation[i] * final_volume[i],
                3e-13,
            );
            near(
                phase.solid.value() + phase.dissolved.value(),
                totals[i],
                1e-15,
            );
            near(
                result.step().solutions[i].water_to_solid_kg,
                transfer[i],
                3e-13,
            );
            near(
                result.step().solutions[i].water_activity[1],
                activity[i],
                2e-12,
            );
        }
        assert!(result.step().solutions[0].water_to_solid_kg > 0.0);
        assert!(result.step().solutions[1].water_to_solid_kg < 0.0);
    });
}

#[test]
fn g1_sealed_baths_keep_solvent_while_native_diffusion_and_finite_kinetics_continue() {
    with_cx(|cx| {
        let volume = 5e-5;
        let law = IdealSolubilityV1::new(
            SpeciesId::new("water").unwrap(),
            SpeciesId::new("authored-neutral-solute").unwrap(),
            Temperature::new(293.15),
            AmountConcentration::new(1000.0),
            0.0,
            [Temperature::new(293.15), Temperature::new(293.15)],
        )
        .unwrap();
        let kinetics = KineticSolubilitySolution::new(
            Temperature::new(293.15),
            WATER_VOLUME,
            Volume::new(volume),
            vec![(
                SolubilityKineticsV1::new(
                    law,
                    VolumetricFlowRate::new(2.5e-5),
                    PrecipitationPolicyV1::ExistingSolidOnly,
                )
                .unwrap(),
                SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.2)).unwrap(),
            )],
        )
        .unwrap();
        let mut expected = kinetics.clone();
        expected.advance(Time::new(2.0), || false).unwrap();
        let baths = [
            reservoir(volume, vec![0.02]),
            FiniteMoistureSolution::with_kinetic_solubility(kinetics, WATER_MOLAR_MASS).unwrap(),
        ];
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.25, 0.1, 0.05], cx);
        let result = baths[0]
            .advance_paired(
                &baths[1],
                &initial,
                2.0,
                &[None, None],
                config(),
                pair_controls([volume; 2]),
                cx,
            )
            .unwrap();
        assert_eq!(result.step().field_evaluations, 1);
        assert_eq!(result.step().transport.boundary_water_kg, [0.0; 2]);
        assert!(result.moisture().moisture()[0] < initial.moisture()[0]);
        assert!(result.moisture().moisture()[2] > initial.moisture()[2]);
        near(
            result.moisture().water_mass_kg(),
            initial.water_mass_kg(),
            1e-13,
        );
        for i in 0..2 {
            assert_eq!(result.solutions()[i].solution().solvent_volume_m3(), volume);
            assert_eq!(result.step().solutions[i].water_to_solid_kg, 0.0);
        }
        let actual = result.solutions()[1].kinetic_solubility().unwrap();
        assert_eq!(actual.last_steps(), expected.last_steps());
        assert!(actual.last_steps()[0].state().dissolved().value() > 0.0);
        assert_eq!(actual.last_steps()[0].state().total_solute().value(), 0.2);
        assert_eq!(
            baths[1].kinetic_solubility().unwrap().last_steps()[0]
                .state()
                .dissolved()
                .value(),
            0.0
        );
    });
}

#[test]
fn g1_cartesian_and_radial_pairs_use_native_areas_storage_and_resistances() {
    with_cx(|cx| {
        let activity = [0.995_f64, 0.98];
        let final_volume = [1e-4; 2];
        let slope = 0.2;
        let diffusion = 1e-5;
        let duration = 0.2;
        let length = 0.04_f64;
        let m0 = 0.1;
        let grid = MoistureCartesian::new(
            card([diffusion; 2], 0.3),
            point(0.1),
            ["moisture_diffusivity"; 3],
            MaterialPropertySelection::SingleClaimOnly,
            MoistureCartesianGeometry {
                lengths_m: [length; 3],
                cells: [1; 3],
                frame: ContentHash([4; 32]),
                maximum_moisture: 0.3,
                maximum_cells: 1,
            },
            vec![m0],
            cx,
        )
        .unwrap();
        let alpha = 6.0 * duration * diffusion / length.powi(2);
        let target = (m0 + alpha * slope * (activity[0] + activity[1])) / (1.0 + 2.0 * alpha);
        let transfers = activity.map(|a| grid.dry_mass_kg() * alpha * (slope * a - target));
        let volumes: [f64; 2] =
            std::array::from_fn(|i| final_volume[i] + transfers[i] / WATER_DENSITY);
        let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
            reservoir(
                volumes[i],
                vec![-activity[i].ln() * final_volume[i] / WATER_VOLUME],
            )
        });
        let base_faces = pair_faces(slope);
        let faces = std::array::from_fn(|i| base_faces[i % 2].clone());
        let result = baths[0]
            .advance_paired_cartesian(
                &baths[1],
                &grid,
                duration,
                &faces,
                config(),
                nonlinear(),
                pair_controls(volumes),
                cx,
            )
            .unwrap();
        near(result.moisture().moisture()[0], target, 3e-12);
        for i in 0..2 {
            near(
                result.step().solutions[i].water_to_solid_kg,
                transfers[i],
                3e-13,
            );
            near(
                result.solutions()[i].solution().solvent_volume_m3(),
                final_volume[i],
                3e-16,
            );
        }
        for side in 0..6 {
            near(
                result.step().transport.boundary_water_kg[side],
                transfers[side % 2] / 3.0,
                3e-13,
            );
        }

        let a = 0.01_f64;
        let b = 0.02_f64;
        let center = 0.5 * (a + b);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.08 },
        ] {
            let field = MoistureRadial::new(
                card([diffusion; 2], 0.3),
                point(0.1),
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                MoistureRadialGeometry {
                    symmetry,
                    faces_m: vec![a, b],
                    maximum_moisture: 0.3,
                    maximum_cells: 1,
                },
                vec![m0],
                cx,
            )
            .unwrap();
            let conductance = match symmetry {
                RadialMoistureSymmetry::Sphere => {
                    let factor = 500.0 * diffusion * 4.0 * core::f64::consts::PI;
                    [
                        factor / (1.0 / a - 1.0 / center),
                        factor / (1.0 / center - 1.0 / b),
                    ]
                }
                RadialMoistureSymmetry::Cylinder { length_m } => {
                    let factor = 500.0 * diffusion * 2.0 * core::f64::consts::PI * length_m;
                    [factor / (center / a).ln(), factor / (b / center).ln()]
                }
            };
            let dry = field.dry_mass_kg();
            let target = (dry * m0
                + duration * slope * (conductance[0] * activity[0] + conductance[1] * activity[1]))
                / (dry + duration * (conductance[0] + conductance[1]));
            let transfer: [f64; 2] =
                std::array::from_fn(|i| duration * conductance[i] * (slope * activity[i] - target));
            let volumes: [f64; 2] =
                std::array::from_fn(|i| final_volume[i] + transfer[i] / WATER_DENSITY);
            let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
                reservoir(
                    volumes[i],
                    vec![-activity[i].ln() * final_volume[i] / WATER_VOLUME],
                )
            });
            let result = baths[0]
                .advance_paired_radial(
                    &baths[1],
                    &field,
                    duration,
                    &base_faces,
                    config(),
                    nonlinear(),
                    pair_controls(volumes),
                    cx,
                )
                .unwrap();
            near(result.moisture().moisture()[0], target, 3e-12);
            for i in 0..2 {
                near(
                    result.step().transport.boundary_water_kg[i],
                    transfer[i],
                    3e-13,
                );
                near(
                    result.solutions()[i].solution().solvent_volume_m3(),
                    final_volume[i],
                    3e-16,
                );
            }
        }
    });
}

#[test]
fn g4_depletion_work_and_identity_refusals_preserve_both_baths_and_field() {
    with_cx(|cx| {
        let volume = 1e-4;
        let first = reservoir(volume, vec![0.02]);
        let second = reservoir(2.0 * volume, vec![0.05]);
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let faces = pair_faces(0.2);
        let complete = first
            .advance_paired(
                &second,
                &original,
                0.25,
                &faces,
                config(),
                pair_controls([volume, 2.0 * volume]),
                cx,
            )
            .unwrap();
        let mut late = pair_controls([volume, 2.0 * volume]);
        late.maximum_evaluations = complete.step().field_evaluations - 1;
        assert!(
            matches!(first.advance_paired(&second, &original, 0.25, &faces,
            config(), late, cx), Err(SolutionMoistureError::EvaluationBudget { evaluations, maximum })
                if evaluations == late.maximum_evaluations && maximum == evaluations)
        );
        let mut bounded = pair_controls([volume, 2.0 * volume]);
        bounded.maximum_evaluations = 1;
        assert!(matches!(
            first.advance_paired(&second, &original, 0.25, &faces, config(), bounded, cx),
            Err(SolutionMoistureError::EvaluationBudget { evaluations: 1, .. })
        ));
        let mut narrow = pair_controls([volume, 2.0 * volume]);
        narrow.solutions[0].solvent_volume_bracket_m3 = [0.999999 * volume, 1.000001 * volume];
        assert!(matches!(
            first.advance_paired(&second, &original, 0.25, &faces, config(), narrow, cx),
            Err(SolutionMoistureError::Unbracketed { .. })
        ));
        let tiny = reservoir(1e-7, vec![]);
        let pure = reservoir(volume, vec![]);
        assert!(matches!(
            tiny.advance_paired(
                &pure,
                &original,
                0.25,
                &faces,
                config(),
                pair_controls([1e-7, volume]),
                cx
            ),
            Err(SolutionMoistureError::Unbracketed { .. })
        ));
        let wrong = FiniteMoistureSolution::new(second.solution().clone(), 0.02).unwrap();
        assert!(matches!(
            first.advance_paired(
                &wrong,
                &original,
                0.25,
                &faces,
                config(),
                pair_controls([volume, 2.0 * volume]),
                cx
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert_eq!(first.solution().solvent_volume_m3(), volume);
        assert_eq!(second.solution().solvent_volume_m3(), 2.0 * volume);
        assert_eq!(first.solution().solute_moles(), &[0.02]);
        assert_eq!(second.solution().solute_moles(), &[0.05]);
        assert_eq!(tiny.solution().solvent_volume_m3(), 1e-7);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.moisture(), &[0.1]);
    });
}

#[test]
fn g4_cancelled_pair_publishes_no_partial_bath_or_specimen() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let volume = 1e-4;
        let baths = [reservoir(volume, vec![0.02]), reservoir(volume, vec![0.05])];
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            baths[0].advance_paired(
                &baths[1],
                &original,
                0.25,
                &pair_faces(0.2),
                config(),
                pair_controls([volume; 2]),
                &cx
            ),
            Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(baths[0].solution().solvent_volume_m3(), volume);
        assert_eq!(baths[1].solution().solvent_volume_m3(), volume);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
