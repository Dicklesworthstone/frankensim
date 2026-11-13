//! Curved-shell diffusion and the actual finite aqueous inventory.
//! Authored analytical fixtures; no measured material-validation claim.
use super::*;
use crate::{atmosphere as air, solution as bath};
use fs_couple::moisture::solution::{
    FiniteMoistureSolution, RadialMoistureParticipant, SolutionMoistureError,
};
use fs_material::osmosis::solubility::FiniteSolubilitySolution;
use fs_qty::{Amount, AmountConcentration, Temperature, Volume};
use fs_thermochem::{SpeciesId, solution::solubility::IdealSolubilityV1};

const DENSITY: f64 = 0.018 / bath::WATER_VOLUME;
const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

fn participant<'a>(
    id: u64,
    moisture: &'a MoistureRadial,
    faces: &'a [Option<fs_couple::moisture::atmosphere::DryBasisSorption>; 2],
) -> RadialMoistureParticipant<'a> {
    RadialMoistureParticipant {
        id,
        moisture,
        faces,
        transport: config(),
        nonlinear: nonlinear(),
    }
}

fn local_water(old: &MoistureRadial, new: &MoistureRadial) -> f64 {
    old.moisture()
        .iter()
        .zip(new.moisture())
        .zip(old.dry_masses_kg())
        .map(|((&a, &b), &dry)| dry * (b - a))
        .sum()
}

#[test]
fn radial_solution_be_matches_curved_resistance_and_finite_solvent_balance() {
    with_cx(|cx| {
        let radius: f64 = 0.01;
        let diffusion = 1e-5;
        let dt = 0.2;
        let slope = 0.25;
        let final_volume: f64 = 1e-4;
        let total_solute = 0.2;
        let activity = (-bath::WATER_VOLUME * total_solute / final_volume).exp();
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.08 },
        ] {
            let rate = match symmetry {
                RadialMoistureSymmetry::Sphere => 3.0 * diffusion / radius.powi(2),
                RadialMoistureSymmetry::Cylinder { .. } => {
                    2.0 * diffusion / (radius.powi(2) * 2.0_f64.ln())
                }
            };
            for m0 in [0.1, 0.28] {
                let initial = specimen(symmetry, vec![0.0, radius], vec![m0], [diffusion; 2], cx);
                let expected = (m0 + dt * rate * slope * activity) / (1.0 + dt * rate);
                let transfer = initial.dry_mass_kg() * (expected - m0);
                let initial_volume = final_volume + transfer / DENSITY;
                let reservoir =
                    bath::reservoir(initial_volume, vec![0.4 * total_solute, 0.6 * total_solute]);
                let original = reservoir.solution().clone();
                let faces = [None, Some(air::sorption(slope))];
                let result = reservoir
                    .advance_radial(
                        &initial,
                        dt,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(initial_volume),
                        cx,
                    )
                    .unwrap();
                near(result.moisture().moisture()[0], expected, 2e-12);
                near(
                    result.reservoir().solution().solvent_volume_m3(),
                    final_volume,
                    2e-15,
                );
                near(result.step().water_activity[1], activity, 2e-12);
                near(result.step().water_to_solid_kg, transfer, 2e-13);
                near(
                    result.step().transport.boundary_water_kg[1],
                    transfer,
                    2e-13,
                );
                assert_eq!(result.step().transport.boundary_water_kg[0], 0.0);
                near(result.step().water_balance_residual_kg, 0.0, 1e-13);
                near(
                    result.reservoir().solvent_mass_kg() + result.moisture().water_mass_kg(),
                    reservoir.solvent_mass_kg() + initial.water_mass_kg(),
                    2e-13,
                );
                assert_eq!(
                    result.reservoir().solution().solute_moles(),
                    original.solute_moles()
                );
                let boundary = result.step().transport.surface_sorption[1]
                    .as_ref()
                    .unwrap();
                near(
                    boundary.query().axes()["relative_humidity"],
                    activity,
                    2e-12,
                );
                assert_eq!(initial.moisture(), &[m0]);
                assert_eq!(reservoir.solution(), &original);
            }
        }
    });
}

#[test]
fn radial_solution_shared_exchange_is_canonical_and_uses_unequal_shell_masses() {
    with_cx(|cx| {
        let volume = 2e-4;
        let reservoir = bath::reservoir(volume, vec![0.8]);
        let a = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.003, 0.012],
            vec![0.03, 0.07],
            [2e-5; 2],
            cx,
        );
        let b = specimen(
            RadialMoistureSymmetry::Cylinder { length_m: 0.05 },
            vec![0.008, 0.01, 0.02],
            vec![0.27, 0.28],
            [5e-6; 2],
            cx,
        );
        let faces_a = [None, Some(air::sorption(0.25))];
        let faces_b = [Some(air::sorption(0.25)), Some(air::sorption(0.25))];
        for second_order in [false, true] {
            // Compare physical endpoints, not report type, across the two methods.
            let run = |reverse: bool| {
                let ordered = if reverse {
                    [participant(29, &b, &faces_b), participant(7, &a, &faces_a)]
                } else {
                    [participant(7, &a, &faces_a), participant(29, &b, &faces_b)]
                };
                if second_order {
                    let (r, solids, step) = reservoir
                        .advance_shared_radial_sdirk2(&ordered, 0.02, 4, bath::controls(volume), cx)
                        .unwrap()
                        .into_parts();
                    (
                        r,
                        solids
                            .into_iter()
                            .map(|s| (s.id, s.moisture))
                            .collect::<Vec<_>>(),
                        step,
                    )
                } else {
                    let (r, solids, step) = reservoir
                        .advance_shared_radial(&ordered, 0.02, 4, bath::controls(volume), cx)
                        .unwrap()
                        .into_parts();
                    (
                        r,
                        solids
                            .into_iter()
                            .map(|s| (s.id, s.moisture))
                            .collect::<Vec<_>>(),
                        step,
                    )
                }
            };
            let (r, solids, step) = run(false);
            let (reversed, reversed_solids, reversed_step) = run(true);
            assert_eq!(solids.iter().map(|s| s.0).collect::<Vec<_>>(), vec![7, 29]);
            assert_eq!(r.solution(), reversed.solution());
            for (x, y) in solids.iter().zip(&reversed_solids) {
                assert_eq!(x.1.moisture(), y.1.moisture());
            }
            assert_eq!(
                step.water_to_solids_kg.to_bits(),
                reversed_step.water_to_solids_kg.to_bits()
            );
            let uptake = local_water(&a, &solids[0].1);
            let release = local_water(&b, &solids[1].1);
            assert!(uptake > 0.0 && release < 0.0);
            near(step.water_to_solids_kg, uptake + release, 2e-14);
            near(
                r.solvent_mass_kg() - reservoir.solvent_mass_kg() + uptake + release,
                0.0,
                1e-13,
            );
            near(step.direct_water_balance_residual_kg, 0.0, 1e-13);
            assert_eq!(
                r.solution().solute_moles(),
                reservoir.solution().solute_moles()
            );
        }
    });
}

fn rk4(initial: f64, duration: f64, steps: usize, rate: impl Fn(f64) -> f64) -> f64 {
    let h = duration / steps as f64;
    let mut value = initial;
    for _ in 0..steps {
        let a = rate(value);
        let b = rate(value + 0.5 * h * a);
        let c = rate(value + 0.5 * h * b);
        let d = rate(value + h * c);
        value += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    value
}

#[test]
fn radial_solution_sdirk_refines_against_independent_finite_bath_ode() {
    with_cx(|cx| {
        let volume = 1e-4;
        let solute = 0.2;
        let slope = 0.25;
        let radius: f64 = 0.03;
        let diffusion = 0.001;
        let rate = 3.0 * diffusion / radius.powi(2);
        let duration = 0.06;
        let faces = [None, Some(air::sorption(slope))];
        for m0 in [0.1, 0.28] {
            let initial = specimen(
                RadialMoistureSymmetry::Sphere,
                vec![0.0, radius],
                vec![m0],
                [diffusion; 2],
                cx,
            );
            let original = bath::reservoir(volume, vec![solute]);
            let rhs = |m: f64| {
                let v = volume - initial.dry_mass_kg() * (m - m0) / DENSITY;
                rate * (slope * (-bath::WATER_VOLUME * solute / v).exp() - m)
            };
            let reference = rk4(m0, duration, 2048, rhs);
            near(reference, rk4(m0, duration, 4096, rhs), 3e-13);
            let errors = [4usize, 8, 16].map(|steps| {
                let mut field = initial.clone();
                let mut reservoir = original.clone();
                let dt = duration / steps as f64;
                for _ in 0..steps {
                    let (next_reservoir, next, step) = reservoir
                        .advance_radial_sdirk2(
                            &field,
                            dt,
                            &faces,
                            config(),
                            nonlinear(),
                            bath::controls(reservoir.solution().solvent_volume_m3()),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    assert_eq!(
                        step.transport.endpoint.interval_s,
                        [field.time_s(), field.time_s() + dt]
                    );
                    near(
                        step.transport.first_stage.interval_s[1],
                        field.time_s() + GAMMA * dt,
                        1e-15,
                    );
                    near(
                        step.transport.endpoint.boundary_water_kg[1],
                        step.water_to_solid_kg,
                        1e-12,
                    );
                    near(step.water_balance_residual_kg, 0.0, 1e-13);
                    field = next;
                    reservoir = next_reservoir;
                }
                near(
                    reservoir.solvent_mass_kg() + field.water_mass_kg(),
                    original.solvent_mass_kg() + initial.water_mass_kg(),
                    3e-12,
                );
                assert_eq!(
                    reservoir.solution().solute_moles(),
                    original.solution().solute_moles()
                );
                (field.moisture()[0] - reference).abs()
            });
            for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
                assert!(ratio > 3.8 && ratio < 4.5, "{errors:?}");
            }
        }
    });
}

#[test]
fn radial_solution_sdirk_retains_finite_crystal_buffer_and_full_water_account() {
    with_cx(|cx| {
        let volume = 1e-4;
        let saturation = 500.0;
        let total = 2.0 * saturation * volume;
        let law = IdealSolubilityV1::new(
            SpeciesId::new("water").unwrap(),
            SpeciesId::new("authored-solute").unwrap(),
            Temperature::new(293.15),
            AmountConcentration::new(saturation),
            0.0,
            [Temperature::new(293.15); 2],
        )
        .unwrap();
        let reservoir = FiniteMoistureSolution::with_solubility(
            FiniteSolubilitySolution::new(
                Temperature::new(293.15),
                bath::WATER_VOLUME,
                Volume::new(volume),
                vec![(law, Amount::new(total))],
            )
            .unwrap(),
            0.018,
        )
        .unwrap();
        let radius: f64 = 0.03;
        let diffusion = 1e-5;
        let dt = 0.2;
        let m0 = 0.1;
        let slope = 0.25;
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![m0],
            [diffusion; 2],
            cx,
        );
        let faces = [None, Some(air::sorption(slope))];
        let candidate = reservoir
            .advance_radial_sdirk2(
                &initial,
                dt,
                &faces,
                config(),
                nonlinear(),
                bath::controls(volume),
                cx,
            )
            .unwrap();
        let activity = (-bath::WATER_VOLUME * saturation).exp();
        let z = dt * 3.0 * diffusion / radius.powi(2);
        let decay = (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2);
        let expected = slope * activity + (m0 - slope * activity) * decay;
        near(candidate.moisture().moisture()[0], expected, 3e-12);
        let final_volume = candidate.reservoir().solution().solvent_volume_m3();
        let phase = candidate
            .reservoir()
            .solubility()
            .unwrap()
            .phase_partitions()[0];
        near(phase.dissolved.value(), saturation * final_volume, 2e-13);
        near(
            phase.solid.value(),
            total - saturation * final_volume,
            2e-13,
        );
        near(candidate.step().water_activity[1], activity, 2e-13);
        near(
            candidate.step().water_to_solid_kg,
            initial.dry_mass_kg() * (expected - m0),
            2e-13,
        );
        near(candidate.step().water_balance_residual_kg, 0.0, 1e-13);
    });
}

#[test]
fn radial_solution_refuses_invalid_topology_and_shared_admission_without_mutation() {
    with_cx(|cx| {
        let volume = 1e-4;
        let reservoir = bath::reservoir(volume, vec![0.2]);
        let before = reservoir.solution().clone();
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.01],
            vec![0.1],
            [1e-5; 2],
            cx,
        );
        let exposed_center = [Some(air::sorption(0.25)), None];
        let faces = [None, Some(air::sorption(0.25))];
        assert!(
            reservoir
                .advance_radial(
                    &initial,
                    0.1,
                    &exposed_center,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    cx
                )
                .is_err()
        );
        for second_order in [false, true] {
            let duplicate = [
                participant(1, &initial, &faces),
                participant(1, &initial, &faces),
            ];
            let invalid = if second_order {
                reservoir
                    .advance_shared_radial_sdirk2(&duplicate, 0.1, 2, bath::controls(volume), cx)
                    .err()
                    .unwrap()
            } else {
                reservoir
                    .advance_shared_radial(&duplicate, 0.1, 2, bath::controls(volume), cx)
                    .err()
                    .unwrap()
            };
            assert!(matches!(invalid, SolutionMoistureError::Invalid(_)));
        }
        assert!(
            reservoir
                .advance_shared_radial(&[], 0.1, 1, bath::controls(volume), cx)
                .is_err()
        );
        let distinct = [
            participant(1, &initial, &faces),
            participant(2, &initial, &faces),
        ];
        assert!(
            reservoir
                .advance_shared_radial_sdirk2(&distinct, 0.1, 1, bath::controls(volume), cx)
                .is_err()
        );
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(reservoir.solution(), &before);
    });
}
