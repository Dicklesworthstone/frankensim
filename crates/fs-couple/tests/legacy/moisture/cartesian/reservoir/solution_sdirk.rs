//! Finite-solution activity, phase buffering and RK water conservation.
//! All coefficients are authored reference data, not measured validation.
use super::*;
use fs_couple::moisture::{
    adaptive::{MoistureAdvanceError, MoistureTimeControls},
    solution::{FiniteMoistureSolution, SolutionMoistureError},
};
use fs_material::osmosis::solubility::{
    FiniteSolubilitySolution, kinetics::KineticSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Volume, VolumetricFlowRate};
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
const WATER_MOLAR_MASS: f64 = 0.018;
const WATER_DENSITY: f64 = WATER_MOLAR_MASS / bath::WATER_VOLUME;

fn decay(z: f64) -> f64 {
    (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2)
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

fn solubility_law(saturation: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-neutral-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap()
}

#[test]
fn g1_kirchhoff_finite_solution_retains_the_first_stage_source_envelope() {
    with_cx(|cx| {
        let volume = 1e-4;
        let reservoir = bath::reservoir(volume, vec![1.2, 0.8]);
        let before = reservoir.solution().clone();
        let initial = grid(
            card([0.01; 2], 0.3),
            grid_geometry([0.1; 3], [1; 3]),
            vec![0.04],
            cx,
        );
        let faces = std::array::from_fn(|_| Some(air::sorption(0.25)));
        let dt = 0.15;
        // With constant D the potential-integrated and ordinary face laws
        // coincide. The lower-volume root probe has a lower activity than the
        // first stage, but must retain the water delivered by that stage.
        let ordinary = reservoir
            .advance_cartesian_sdirk2(
                &initial,
                dt,
                &faces,
                config(),
                nonlinear(),
                bath::controls(volume),
                cx,
            )
            .unwrap();
        let integrated = initial.clone().with_kirchhoff_diffusivity();
        let candidate = reservoir
            .advance_cartesian_sdirk2(
                &integrated,
                dt,
                &faces,
                config(),
                nonlinear(),
                bath::controls(volume),
                cx,
            )
            .unwrap();
        near(
            candidate.moisture().moisture()[0],
            ordinary.moisture().moisture()[0],
            3e-12,
        );
        near(
            candidate.reservoir().solution().solvent_volume_m3(),
            ordinary.reservoir().solution().solvent_volume_m3(),
            3e-15,
        );
        near(
            candidate.reservoir().solvent_mass_kg() + candidate.moisture().water_mass_kg(),
            reservoir.solvent_mass_kg() + initial.water_mass_kg(),
            2e-12,
        );
        assert!(
            candidate
                .step()
                .transport
                .endpoint
                .kirchhoff_sources
                .is_some()
        );
        assert_eq!(
            candidate.reservoir().solution().solute_moles(),
            before.solute_moles()
        );
        assert_eq!(reservoir.solution(), &before);
        assert_eq!(integrated.moisture(), &[0.04]);
        assert_eq!(integrated.time_s(), 0.0);
    });
}

#[test]
fn g1_solution_sdirk_tracks_finite_activity_feedback_with_second_order_accuracy() {
    with_cx(|cx| {
        let volume = 1e-4;
        let solute = 0.02;
        let original = bath::reservoir(volume, vec![0.6 * solute, 0.4 * solute]);
        let before = original.solution().clone();
        let geometry = grid_geometry([0.1; 3], [1; 3]);
        let diffusion = 0.01;
        let slope = 0.25;
        let faces = std::array::from_fn(|_| Some(air::sorption(slope)));
        let rate = 12.0 * diffusion / geometry.lengths_m[0].powi(2);
        let duration = 0.08;
        for m0 in [0.1, 0.28] {
            let initial = grid(card([diffusion; 2], 0.3), geometry, vec![m0], cx);
            let dry = initial.dry_mass_kg();
            // Eliminate the finite solvent inventory. This independent scalar
            // ODE uses the ideal dilute activity, not the production integrator.
            let rhs = |m: f64| {
                let v = volume - dry * (m - m0) / WATER_DENSITY;
                let activity = (-bath::WATER_VOLUME * solute / v).exp();
                rate * (slope * activity - m)
            };
            let reference = rk4(m0, duration, 2048, rhs);
            near(reference, rk4(m0, duration, 4096, rhs), 5e-13);
            let errors = [4usize, 8, 16].map(|steps| {
                let dt = duration / steps as f64;
                let mut reservoir = original.clone();
                let mut field = initial.clone();
                for _ in 0..steps {
                    let (next_reservoir, next_field, step) = reservoir
                        .advance_cartesian_sdirk2(
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
                    near(step.water_balance_residual_kg, 0.0, 1e-13);
                    near(
                        step.transport.endpoint.boundary_water_kg.iter().sum(),
                        step.water_to_solid_kg,
                        1e-11,
                    );
                    near(
                        step.solvent_mass_kg[1] - step.solvent_mass_kg[0] + step.water_to_solid_kg,
                        0.0,
                        1e-13,
                    );
                    for source in step.transport.endpoint.surface_sorption.iter().flatten() {
                        near(
                            source.query().axes()["relative_humidity"],
                            step.water_activity[1],
                            2e-13,
                        );
                    }
                    reservoir = next_reservoir;
                    field = next_field;
                }
                let final_volume = volume - dry * (field.moisture()[0] - m0) / WATER_DENSITY;
                near(
                    reservoir.solution().solvent_volume_m3(),
                    final_volume,
                    3e-15,
                );
                near(
                    reservoir.solution().water_activity().unwrap(),
                    (-bath::WATER_VOLUME * solute / final_volume).exp(),
                    2e-11,
                );
                assert_eq!(
                    reservoir.solution().solute_moles(),
                    original.solution().solute_moles()
                );
                near(
                    reservoir.solvent_mass_kg() + field.water_mass_kg(),
                    original.solvent_mass_kg() + initial.water_mass_kg(),
                    3e-12,
                );
                near(field.time_s(), duration, 1e-15);
                assert_eq!(field.moisture()[0] > m0, m0 < slope);
                assert_eq!(
                    reservoir.solution().water_activity().unwrap()
                        < original.solution().water_activity().unwrap(),
                    m0 < slope,
                );
                (field.moisture()[0] - reference).abs()
            });
            for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
                assert!(ratio > 3.8 && ratio < 4.5, "{errors:?}");
            }
            assert_eq!(initial.moisture(), &[m0]);
            assert_eq!(initial.time_s(), 0.0);
        }
        assert_eq!(original.solution(), &before);
    });
}

#[test]
fn g1_solution_sdirk_matches_pure_solvent_and_saturated_phase_buffers() {
    with_cx(|cx| {
        let volume = 5e-5;
        let slope = 0.2;
        let diffusion = 1e-5;
        let geometry = grid_geometry([0.04; 3], [1; 3]);
        let rate = 12.0 * diffusion / geometry.lengths_m[0].powi(2);
        let dt = 2.0;
        let faces = std::array::from_fn(|_| Some(air::sorption(slope)));
        for saturation in [0.0, 1000.0] {
            let total = 2.0 * saturation * volume;
            let reservoir = if saturation == 0.0 {
                bath::reservoir(volume, Vec::new())
            } else {
                FiniteMoistureSolution::with_solubility(
                    FiniteSolubilitySolution::new(
                        Temperature::new(293.15),
                        bath::WATER_VOLUME,
                        Volume::new(volume),
                        vec![(solubility_law(saturation), Amount::new(total))],
                    )
                    .unwrap(),
                    WATER_MOLAR_MASS,
                )
                .unwrap()
            };
            let before = reservoir.solubility().cloned();
            let activity = (-bath::WATER_VOLUME * saturation).exp();
            let equilibrium = slope * activity;
            for m0 in [0.1, 0.27] {
                let initial = grid(card([diffusion; 2], 0.3), geometry, vec![m0], cx);
                let dry = initial.dry_mass_kg();
                let first = equilibrium + (m0 - equilibrium) / (1.0 + GAMMA * dt * rate);
                let expected = equilibrium + (m0 - equilibrium) * decay(rate * dt);
                let water = dry * (expected - m0);
                let expected_volume = volume - water / WATER_DENSITY;
                let candidate = reservoir
                    .advance_cartesian_sdirk2(
                        &initial,
                        dt,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(volume),
                        cx,
                    )
                    .unwrap();
                let report = candidate.step();
                near(candidate.moisture().moisture()[0], expected, 3e-13);
                near(
                    candidate.reservoir().solution().solvent_volume_m3(),
                    expected_volume,
                    2e-16,
                );
                near(report.water_to_solid_kg, water, 2e-13);
                near(report.water_balance_residual_kg, 0.0, 1e-13);
                near(
                    report.transport.endpoint.boundary_water_kg.iter().sum(),
                    water,
                    2e-13,
                );
                near(
                    report.transport.first_stage.transport_end[0].query().axes()["moisture_content"],
                    first,
                    3e-13,
                );
                near(
                    report.transport.first_stage.boundary_water_kg.iter().sum(),
                    dry * GAMMA * dt * rate * (equilibrium - first),
                    2e-13,
                );
                assert_eq!(report.transport.endpoint.interval_s, [0.0, dt]);
                for actual in report.water_activity {
                    near(actual, activity, 2e-14);
                }
                if let Some(phases) = candidate.reservoir().solubility() {
                    let phase = phases.phase_partitions()[0];
                    near(phase.dissolved.value(), saturation * expected_volume, 2e-13);
                    near(
                        phase.solid.value(),
                        total - saturation * expected_volume,
                        2e-13,
                    );
                    near(phase.dissolved.value() + phase.solid.value(), total, 2e-16);
                    let initial_solid =
                        before.as_ref().unwrap().phase_partitions()[0].solid.value();
                    near(
                        phase.solid.value() - initial_solid,
                        saturation * water / WATER_DENSITY,
                        2e-13,
                    );
                    assert!(phase.solid.value() > 0.0);
                    assert_eq!(phases.solution().solute_moles(), &[phase.dissolved.value()]);
                } else {
                    assert_eq!(report.water_activity, [1.0; 2]);
                    assert_eq!(report.iterations, 0);
                }
                assert_eq!(initial.moisture(), &[m0]);
            }
            assert_eq!(reservoir.solubility(), before.as_ref());
            assert_eq!(reservoir.solution().solvent_volume_m3(), volume);
        }
    });
}

#[test]
fn g1_shared_solution_sdirk_retains_opposing_transfers_when_net_water_is_zero() {
    with_cx(|cx| {
        let volume = 5e-5;
        let reservoir = bath::reservoir(volume, vec![0.1]);
        let before = reservoir.solution().clone();
        let activity = reservoir.solution().water_activity().unwrap();
        let equilibrium = 0.2 * activity;
        let ga = grid_geometry([0.04; 3], [1; 3]);
        let gb = grid_geometry([0.02; 3], [1; 3]);
        // Five active faces on a and six on b have the same decay rate. Their
        // dry-weighted excursions cancel at every time, not just the endpoint.
        let a = grid(card([1.2e-5; 2], 0.3), ga, vec![equilibrium - 0.004], cx);
        let volume_ratio =
            ga.lengths_m.iter().product::<f64>() / gb.lengths_m.iter().product::<f64>();
        let b0 = equilibrium + volume_ratio * 0.004;
        let b = grid(card([2.5e-6; 2], 0.3), gb, vec![b0], cx);
        let dry = [a.dry_mass_kg(), b.dry_mass_kg()];
        let old = [a.moisture()[0], b.moisture()[0]];
        let rate = 10.0 * 1.2e-5 / ga.lengths_m[0].powi(2);
        near(rate, 12.0 * 2.5e-6 / gb.lengths_m[0].powi(2), 1e-16);
        let dt = 2.0;
        let expected = old.map(|m| equilibrium + (m - equilibrium) * decay(rate * dt));
        let mut fa = std::array::from_fn(|_| Some(air::sorption(0.2)));
        fa[0] = Some(emission::law([0.0; 2]));
        let fb = std::array::from_fn(|_| Some(air::sorption(0.2)));
        let run = |reverse| {
            let entries = if reverse {
                [participant(3, &b, &fb), participant(7, &a, &fa)]
            } else {
                [participant(7, &a, &fa), participant(3, &b, &fb)]
            };
            reservoir
                .advance_shared_cartesian_sdirk2(&entries, dt, 2, bath::controls(volume), cx)
                .unwrap()
        };
        let accepted = run(false);
        let replay = run(true);
        assert_eq!(
            accepted
                .solids()
                .iter()
                .map(|body| body.id)
                .collect::<Vec<_>>(),
            vec![3, 7]
        );
        assert_eq!(
            accepted.reservoir().solution(),
            replay.reservoir().solution()
        );
        assert_eq!(accepted.reservoir().solution().solvent_volume_m3(), volume);
        near(accepted.step().water_to_solids_kg, 0.0, 1e-13);
        near(accepted.step().water_balance_residual_kg, 0.0, 1e-13);
        near(accepted.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(accepted.step().water_activity, [activity; 2]);
        for (body, repeated) in accepted.solids().iter().zip(replay.solids()) {
            let i = if body.id == 7 { 0 } else { 1 };
            let water = dry[i] * (expected[i] - old[i]);
            near(body.moisture.moisture()[0], expected[i], 3e-13);
            near(
                body.transport.endpoint.boundary_water_kg.iter().sum(),
                water,
                2e-13,
            );
            assert_eq!(water.is_sign_positive(), i == 0);
            assert!(water.abs() > 1e-6);
            assert_eq!(body.moisture.moisture(), repeated.moisture.moisture());
            assert_eq!(
                body.transport.endpoint.boundary_water_kg,
                repeated.transport.endpoint.boundary_water_kg
            );
            assert_eq!(body.transport.endpoint.interval_s, [0.0, dt]);
        }
        assert_eq!(
            accepted.solids()[1].transport.endpoint.boundary_water_kg[0],
            0.0
        );
        assert_eq!(
            accepted.solids()[1].transport.first_stage.boundary_water_kg[0],
            0.0
        );
        assert_eq!(reservoir.solution(), &before);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g4_solution_sdirk_refuses_unstaged_kinetics_and_cancellation_transactionally() {
    with_cx(|cx| {
        let volume = 5e-5;
        let kinetic = FiniteMoistureSolution::with_kinetic_solubility(
            KineticSolubilitySolution::new(
                Temperature::new(293.15),
                bath::WATER_VOLUME,
                Volume::new(volume),
                vec![(
                    SolubilityKineticsV1::new(
                        solubility_law(1000.0),
                        VolumetricFlowRate::new(2.5e-5),
                        PrecipitationPolicyV1::ExistingSolidOnly,
                    )
                    .unwrap(),
                    SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.2))
                        .unwrap(),
                )],
            )
            .unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        let before = kinetic.kinetic_solubility().unwrap().clone();
        let initial = grid(
            card([1e-5; 2], 0.3),
            grid_geometry([0.04; 3], [1; 3]),
            vec![0.1],
            cx,
        );
        let faces = std::array::from_fn(|_| Some(air::sorption(0.2)));
        for result in [
            kinetic
                .advance_cartesian_sdirk2(
                    &initial,
                    2.0,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    cx,
                )
                .map(|_| ()),
            kinetic
                .advance_shared_cartesian_sdirk2(
                    &[participant(7, &initial, &faces)],
                    2.0,
                    1,
                    bath::controls(volume),
                    cx,
                )
                .map(|_| ()),
        ] {
            let error = result.unwrap_err();
            assert!(
                matches!(
                    &error,
                    SolutionMoistureError::Invalid(message)
                        if *message == "second-order solution coupling requires algebraic phase partitioning"
                ),
                "{error:?}"
            );
        }
        assert_eq!(kinetic.kinetic_solubility(), Some(&before));
        let dissolved = bath::reservoir(volume, vec![0.02]);
        let before_dissolved = dissolved.solution().clone();
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
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
            assert!(matches!(
                dissolved.advance_cartesian_sdirk2(
                    &initial,
                    1.0,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    &cancelled,
                ),
                Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
            ));
        });
        assert_eq!(dissolved.solution(), &before_dissolved);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g1_solution_sdirk_adaptive_controls_error_and_only_debits_accepted_transfers() {
    with_cx(|cx| {
        let volume = 1e-4;
        let solute = 0.02;
        let reservoir = bath::reservoir(volume, vec![solute]);
        let before = reservoir.solution().clone();
        let geometry = grid_geometry([0.1; 3], [1; 3]);
        let initial = grid(card([0.01; 2], 0.3), geometry, vec![0.1], cx);
        let dry = initial.dry_mass_kg();
        let slope = 0.25;
        let rate = 12.0 * 0.01 / geometry.lengths_m[0].powi(2);
        let target = 0.08;
        let faces = std::array::from_fn(|_| Some(air::sorption(slope)));
        let rhs = |m: f64| {
            let v = volume - dry * (m - 0.1) / WATER_DENSITY;
            rate * (slope * (-bath::WATER_VOLUME * solute / v).exp() - m)
        };
        let reference = rk4(0.1, target, 4096, rhs);
        let time = |absolute| MoistureTimeControls {
            initial_step_s: target,
            minimum_step_s: target / 4096.0,
            maximum_step_s: target,
            maximum_evaluations: 20_000,
            moisture_absolute_tolerance: absolute,
            activity_absolute_tolerance: absolute,
            reservoir_water_absolute_tolerance_kg: dry * absolute,
            relative_tolerance: 0.0,
            water_balance_tolerance_kg: 1e-11,
        };
        let run = |absolute| {
            reservoir
                .advance_cartesian_sdirk2_adaptive_to(
                    &initial,
                    target,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    time(absolute),
                    cx,
                )
                .unwrap()
        };
        let (loose_bath, loose, loose_report) = run(1e-5);
        let (tight_bath, tight, tight_report) = run(1e-7);
        let loose_error = (loose.moisture()[0] - reference).abs();
        let tight_error = (tight.moisture()[0] - reference).abs();
        assert!(
            tight_error < loose_error / 4.0,
            "{loose_error} -> {tight_error}"
        );
        assert!(tight_report.accepted_substeps > loose_report.accepted_substeps);
        for (bath, field, report) in [
            (&loose_bath, &loose, &loose_report),
            (&tight_bath, &tight, &tight_report),
        ] {
            assert_eq!(field.time_s(), target);
            assert_eq!(report.interval_s, [0.0, target]);
            assert!(report.rejected_intervals > 0);
            assert!(report.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(report.evaluations % 2, 0);
            near(
                report.water_change_kg,
                field.water_mass_kg() - initial.water_mass_kg(),
                1e-12,
            );
            near(
                report.boundary_water_kg.iter().sum(),
                report.water_change_kg,
                1e-11,
            );
            near(
                reservoir.solvent_mass_kg() - bath.solvent_mass_kg(),
                report.water_change_kg,
                1e-11,
            );
            near(report.water_balance_residual_kg, 0.0, 1e-11);
            near(report.transport_balance_residual_kg, 0.0, 1e-11);
            assert_eq!(
                bath.solution().solute_moles(),
                reservoir.solution().solute_moles()
            );
            assert_eq!(
                report.reservoir_water_kg,
                Some([reservoir.solvent_mass_kg(), bath.solvent_mass_kg()])
            );
            assert_eq!(report.last_substep.transport.endpoint.interval_s[1], target);
        }
        let exhausted = MoistureTimeControls {
            maximum_evaluations: 5,
            ..time(1e-5)
        };
        assert!(matches!(
            reservoir.advance_cartesian_sdirk2_adaptive_to(
                &initial,
                target,
                &faces,
                config(),
                nonlinear(),
                bath::controls(volume),
                exhausted,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                accepted_intervals: 0,
                evaluations: 4,
                ..
            })
        ));
        assert_eq!(reservoir.solution(), &before);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
    });
}
