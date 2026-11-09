//! G0/G1/G3/G4: adaptive open molecular sources and finite material heat storage.
use super::*;
use crate::rocket_blowdown::effusion::{
    exchange::{EffusionAmbientPort, EffusionReservoirPort},
    thermal::trajectory::exchange::EffusionWallExchangeTrajectoryCandidate,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::effusion::drift::EffusiveDriftInterface;
use fs_thermochem::{
    Composition, CompositionBasis, Nasa9RegionV1, Nasa9StandardStateModelV1, SpeciesId,
    mixture::FrozenIdealGasMixtureModelV1,
};

fn supply<'a>(
    gas: &'a ClosedGasParcel,
    law: &'a EffusiveMomentumInterface,
) -> EffusionReservoirPort<'a> {
    EffusionReservoirPort {
        id: 9,
        source: gas,
        interface: law,
        area_m2: AREA,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn bath(model: &ThermochemicalGasModel, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model.clone(), source::point(t, p).unwrap(), VOLUME).unwrap()
}
fn ambient<'a>(
    gas: &'a ClosedGasParcel,
    law: &'a EffusiveDriftInterface,
    velocity: Vec3,
) -> EffusionAmbientPort<'a> {
    EffusionAmbientPort {
        id: 3,
        source: gas,
        interface: law,
        area_m2: AREA,
        velocity_world_m_s: velocity,
        inertial_frame: FRAME,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn accounts(
    before: &GasRocketState,
    input: UniformEnthalpyStepInput<'_>,
    result: &EffusionWallExchangeTrajectoryCandidate,
) {
    let (after, report) = (result.endpoint(), result.step());
    assert_eq!(report.duration_s, input.duration_s);
    assert_eq!(report.accepted_steps.len(), 2 * report.accepted_intervals);
    assert_eq!(
        report.attempted_intervals,
        report.accepted_intervals + report.rejected_intervals
    );
    assert!(report.maximum_accepted_error <= 1.0);
    near(
        after.reservoir().mass_kg() - before.reservoir().mass_kg(),
        report.net_mass_into_gas_kg,
        1e-12,
    );
    vector_near(
        after
            .vehicle()
            .velocity_m_s
            .scale(mass(after))
            .sub(before.vehicle().velocity_m_s.scale(mass(before)))
            .sub(report.net_momentum_into_vessel_world_kg_m_s)
            .sub(report.external_impulse_n_s),
        Vec3::ZERO,
        2e-9,
    );
    near(
        internal(after) - internal(before) + kinetic(after) - kinetic(before)
            + report.wall_storage_change_j
            - report.net_energy_into_vessel_world_j
            - report.external_work_j
            - report.external_thermal_heat_j,
        0.0,
        1e-4,
    );
    assert!(report.maximum_species_mass_residual_kg <= adaptive(1.0).window_mass_tolerance_kg);
    assert!(report.debit.energy_j <= adaptive(1.0).window_energy_tolerance_j);
}

#[test]
fn entering_species_controls_refinement_and_all_accepted_wall_source_steps_replay() {
    scope(|cx, _| {
        let catalog = source::load(MIXTURE, "effusive recoil analytic fixture").unwrap();
        let pure = |fractions| {
            catalog
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
                )
                .unwrap()
        };
        let before = GasRocketState::new(
            bath(&pure(vec![1.0, 0.0]), T0, P0),
            0.01,
            Vec3::ZERO,
            Vec3::ZERO,
            FRAME,
        )
        .unwrap();
        let source = bath(&pure(vec![0.0, 1.0]), 700.0, 1.3 * P0);
        let saved_source = source.clone();
        let b = SpeciesId::new("B").unwrap();
        let base = interface(&catalog);
        let law = EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                catalog.clone(),
                catalog
                    .thermodynamics()
                    .components()
                    .iter()
                    .map(|species| {
                        EffusiveSpeciesTransmission::try_new(
                            species.species().clone(),
                            if species.species() == &b { 1.0 } else { 0.0 },
                            "synthetic B-only aperture",
                        )
                        .unwrap()
                    })
                    .collect(),
                base.interface().support(),
                base.interface().regime().clone(),
            )
            .unwrap(),
            "synthetic uncollimated B-only transmission",
        )
        .unwrap();
        let drift = EffusiveDriftInterface::try_new(
            law.clone(),
            1000.0,
            "synthetic Maxwellian drift support",
        )
        .unwrap();
        let ports = [supply(&source, &law)];
        let ambients = [ambient(&source, &drift, Vec3::new(-30.0, 40.0, 0.0))];
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.3,
            internal_heat_j: 1.3,
            ..wall(&chart, 600.0)
        };
        let accuracy = |mass| {
            let mut a = adaptive(1e8);
            a.motion.relative_tolerance = 0.0;
            a.motion.species_mass_tolerance_kg = mass;
            a
        };
        let run = |mass| {
            before
                .advance_effusion_exchange_with_wall_adaptive(
                    cx,
                    free(),
                    &[],
                    &ports,
                    &ambients,
                    input,
                    &body(100.0),
                    0.4,
                    solve(),
                    accuracy(mass),
                )
                .unwrap()
        };
        let loose = run(2e-8);
        let tight = run(2e-10);
        assert_eq!(before.reservoir().species_mass_kg(&b).unwrap(), 0.0);
        assert!(tight.endpoint().reservoir().species_mass_kg(&b).unwrap() > 0.0);
        assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
        assert!(tight.step().rejected_intervals > 0);
        let (mut current, mut wall_state) = (before.clone(), input.initial);
        let (mut elapsed, mut supplied, mut discharged) = (0.0, 0.0, 0.0);
        for accepted in &tight.step().accepted_steps {
            let dt = accepted.duration_s;
            let next = current
                .advance_effusion_exchange_with_wall(
                    cx,
                    free(),
                    &[],
                    &ports,
                    &ambients,
                    UniformEnthalpyStepInput {
                        initial: wall_state,
                        duration_s: dt,
                        internal_heat_j: input.internal_heat_j * (dt / input.duration_s),
                        ..input
                    },
                    &body(100.0),
                    0.4 * (dt / input.duration_s),
                    solve(),
                )
                .unwrap();
            assert_eq!(next.step(), accepted);
            supplied += accepted.exchange.incoming.mass_kg;
            discharged += accepted.exchange.outgoing.mass_kg;
            (current, wall_state, _) = next.into_parts();
            elapsed += dt;
        }
        assert_eq!(&current, tight.endpoint());
        assert_eq!(&wall_state, tight.wall());
        near(elapsed, input.duration_s, 2e-15);
        near(supplied, tight.step().incoming.mass_kg, 2e-17);
        near(discharged, tight.step().outgoing.mass_kg, 2e-17);
        for residual in tight.step().prescribed_heat_residual_j {
            near(residual, 0.0, 2e-14);
        }
        accounts(&before, input, &tight);
        assert_eq!(source, saved_source);
    });
}

// Independent continuous monatomic gas/wall/translation balances with one
// comoving source and one vacuum aperture. No native sampler supplies this ODE.
fn continuous_source(duration: f64) -> [f64; 7] {
    let derivative = |y: [f64; 7]| {
        let [n, t, w, _, v, _, _] = y;
        let out = AREA * n / VOLUME * (R * t / (2.0 * core::f64::consts::PI * 0.028)).sqrt();
        let inc = AREA * 0.9 * P0 / (2.0 * core::f64::consts::PI * 0.028 * R * 300.0).sqrt();
        let dn = inc - 2.0 * out;
        let heat = 100.0 * WALL_AREA * (w - t);
        let molecular = 2.0 * R * (300.0 * inc - 2.0 * t * out);
        let impulse = -0.5 * AREA * (0.9 * P0 + 2.0 * n * R * t / VOLUME);
        [
            dn,
            (heat + molecular - 1.5 * R * t * dn) / (1.5 * R * n),
            -heat / (WALL_MASS * WALL_CP),
            v,
            impulse / (0.01 + 0.028 * n),
            0.028 * dn * v + impulse,
            molecular + v * impulse + 0.5 * 0.028 * dn * v * v,
        ]
    };
    let mut y = [P0 * VOLUME / (R * T0), T0, 600.0, 0.0, 0.0, 0.0, 0.0];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..7 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn tighter_wall_supply_discharge_histories_converge_to_the_independent_continuous_balance() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let source = bath(before.reservoir().model(), 300.0, 0.9 * P0);
        let law = interface(before.reservoir().model());
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.3,
            ..wall(&chart, 600.0)
        };
        let run = |scale| {
            before
                .advance_effusion_exchange_with_wall_adaptive(
                    cx,
                    free(),
                    &[outlet(&law, Vec3::new(1.0, 0.0, 0.0))],
                    &[supply(&source, &law)],
                    &[],
                    input,
                    &body(100.0),
                    0.0,
                    solve(),
                    adaptive(scale),
                )
                .unwrap()
        };
        let (loose, tight) = (run(1.0), run(0.04));
        let reference = continuous_source(input.duration_s);
        let discrepancy = |result: &EffusionWallExchangeTrajectoryCandidate| {
            let (end, step) = (result.endpoint(), result.step());
            [
                end.reservoir().mass_kg() / 0.028,
                end.reservoir().state().state().temperature,
                result.wall().temperature_k(),
                end.vehicle().position_m.x,
                end.vehicle().velocity_m_s.x,
                step.net_momentum_into_vessel_world_kg_m_s.x,
                step.net_energy_into_vessel_world_j,
            ]
            .into_iter()
            .zip(reference)
            .zip([0.03, 400.0, 600.0, 0.01, 0.1, 0.002, 2.0])
            .map(|((a, b), scale)| (a - b).abs() / scale)
            .sum::<f64>()
        };
        assert!(discrepancy(&tight) < 0.4 * discrepancy(&loose));
        assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
        accounts(&before, input, &tight);
    });
}

fn shifted_model(datum: f64) -> ThermochemicalGasModel {
    let model = source::load(MIXTURE, "effusive recoil analytic fixture").unwrap();
    let species = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            Nasa9StandardStateModelV1::new(
                species.species().clone(),
                species.molar_mass(),
                species.convention().clone(),
                species
                    .regions()
                    .iter()
                    .map(|region| {
                        let mut card = region.card().clone();
                        card.parameters.get_mut("a7").unwrap().value += datum;
                        Nasa9RegionV1::from_card(card).unwrap()
                    })
                    .collect(),
            )
            .unwrap()
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            model.thermodynamics().mole_composition().clone(),
        )
        .unwrap(),
        model.transport_materials().to_vec(),
        None,
    )
    .unwrap()
}

#[test]
fn both_gross_energy_estimates_ignore_formation_datums_and_obey_common_inertial_boosts() {
    scope(|cx, _| {
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..wall(&chart, 600.0)
        };
        let run = |datum, boost: Vec3, accuracy| {
            let model = shifted_model(datum);
            let before =
                GasRocketState::new(bath(&model, T0, P0), 0.01, Vec3::ZERO, boost, FRAME).unwrap();
            let source = bath(&model, 700.0, 1.3 * P0);
            let law = EffusiveDriftInterface::try_new(
                interface(&model),
                1000.0,
                "synthetic Maxwellian drift support",
            )
            .unwrap();
            before
                .advance_effusion_exchange_with_wall_adaptive(
                    cx,
                    free(),
                    &[],
                    &[],
                    &[ambient(
                        &source,
                        &law,
                        Vec3::new(-40.0, 70.0, 0.0).add(boost),
                    )],
                    input,
                    &body(100.0),
                    0.0,
                    solve(),
                    accuracy,
                )
                .unwrap()
        };
        let mut accuracy = adaptive(1e8);
        accuracy.motion.relative_tolerance = 0.0;
        accuracy.emitted_energy_tolerance_j = 1e-4;
        let a = run(0.0, Vec3::ZERO, accuracy);
        let shifted = run(1e12, Vec3::ZERO, accuracy);
        assert!(
            a.step().rejected_intervals > 0,
            "molecular energy alone must refine"
        );
        assert_eq!(
            a.step().accepted_intervals,
            shifted.step().accepted_intervals
        );
        near(
            a.endpoint().reservoir().state().state().temperature,
            shifted.endpoint().reservoir().state().state().temperature,
            1e-6,
        );
        near(
            a.step().incoming_reference_energy_world_j,
            shifted.step().incoming_reference_energy_world_j,
            2e-7,
        );
        near(
            a.step().outgoing_reference_energy_world_j,
            shifted.step().outgoing_reference_energy_world_j,
            2e-7,
        );
        near(
            a.step().net_reference_energy_into_vessel_world_j,
            shifted.step().net_reference_energy_into_vessel_world_j,
            2e-7,
        );

        // Identical fixed schedules isolate observer transformation from the
        // deliberately frame-specific local export tolerance.
        accuracy.emitted_energy_tolerance_j = 1e6;
        accuracy.motion.initial_step_s = 0.025;
        accuracy.motion.maximum_step_s = 0.025;
        let boost = Vec3::new(31.0, -7.0, 13.0);
        let (a, b) = (run(0.0, Vec3::ZERO, accuracy), run(0.0, boost, accuracy));
        assert_eq!(a.step().accepted_intervals, b.step().accepted_intervals);
        vector_near(
            b.endpoint().vehicle().velocity_m_s,
            a.endpoint().vehicle().velocity_m_s.add(boost),
            2e-10,
        );
        for (first, second, reference, boosted_reference) in [
            (
                a.step().incoming,
                b.step().incoming,
                a.step().incoming_reference_energy_world_j,
                b.step().incoming_reference_energy_world_j,
            ),
            (
                a.step().outgoing,
                b.step().outgoing,
                a.step().outgoing_reference_energy_world_j,
                b.step().outgoing_reference_energy_world_j,
            ),
        ] {
            near(first.mass_kg, second.mass_kg, 2e-16);
            let correction =
                boost.dot(first.momentum_world_kg_m_s) + 0.5 * first.mass_kg * boost.norm_squared();
            near(
                second.energy_world_j,
                first.energy_world_j + correction,
                2e-7,
            );
            near(boosted_reference, reference + correction, 2e-7);
        }
        near(
            b.step().net_reference_energy_into_vessel_world_j,
            a.step().net_reference_energy_into_vessel_world_j
                + boost.dot(a.step().net_momentum_into_vessel_world_kg_m_s)
                + 0.5 * a.step().net_mass_into_gas_kg * boost.norm_squared(),
            2e-7,
        );
    });
}

#[test]
fn fixed_window_budgets_and_original_source_work_cancellation_refusals_are_atomic() {
    scope(|cx, gate| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let source = bath(before.reservoir().model(), 700.0, P0);
        let saved_source = source.clone();
        let law = interface(before.reservoir().model());
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let ports = [supply(&source, &law)];
        let run = |accuracy| {
            before.advance_effusion_exchange_with_wall_adaptive(
                cx,
                free(),
                &[],
                &ports,
                &[],
                input,
                &body(100.0),
                0.0,
                solve(),
                accuracy,
            )
        };
        let accepted = run(adaptive(1.0)).unwrap();
        assert!(accepted.step().debit.energy_j > 0.0);
        let mut bounded = adaptive(1.0);
        bounded.window_energy_tolerance_j = 0.5 * accepted.step().debit.energy_j;
        assert!(matches!(
            run(bounded),
            Err(EffusionWallTrajectoryError::WindowBalance { .. })
        ));
        let mut exhausted = adaptive(1e-7);
        exhausted.motion.maximum_attempts = 1;
        assert!(matches!(
            run(exhausted),
            Err(EffusionWallTrajectoryError::Control(
                EffusionTrajectoryError::WorkBudget { .. }
            ))
        ));
        exhausted.motion.maximum_attempts = 4096;
        exhausted.motion.minimum_step_s = exhausted.motion.initial_step_s;
        assert!(matches!(
            run(exhausted),
            Err(EffusionWallTrajectoryError::Control(
                EffusionTrajectoryError::MinimumStep { .. }
            ))
        ));
        let narrow = EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                law.interface().gas_catalog().clone(),
                law.interface().transmissions().to_vec(),
                EffusiveSupport {
                    temperature_k: (410.0, 1800.0),
                    ..law.interface().support()
                },
                law.interface().regime().clone(),
            )
            .unwrap(),
            law.provenance(),
        )
        .unwrap();
        let mut calls = 0;
        assert!(matches!(
            before.advance_effusion_exchange_with_wall_convection_adaptive(
                cx,
                free(),
                &[],
                &[supply(&source, &narrow)],
                &[],
                input,
                &body(100.0),
                0.0,
                solve(),
                adaptive(1.0),
                |_| {
                    calls += 1;
                    Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                        coefficient_w_per_m2_k: 100.0,
                        source: (),
                    })
                },
            ),
            Err(EffusionWallTrajectoryError::Step {
                elapsed_s: 0.0,
                error: EffusionWallError::Effusion(EffusionRocketError::Source { .. }),
                ..
            })
        ));
        assert_eq!(calls, 0);
        let closed = before
            .advance_effusion_exchange_with_wall_adaptive(
                cx,
                free(),
                &[],
                &[EffusionReservoirPort {
                    area_m2: 0.0,
                    ..supply(&source, &narrow)
                }],
                &[],
                input,
                &body(100.0),
                0.0,
                solve(),
                adaptive(1.0),
            )
            .unwrap();
        assert_eq!(closed.step().incoming.mass_kg, 0.0);
        assert_eq!(closed.step().outgoing.mass_kg, 0.0);
        assert!(closed.wall().temperature_k() < input.initial.temperature_k());
        let cancelled = before.advance_effusion_exchange_with_wall_convection_adaptive(
            cx,
            free(),
            &[],
            &ports,
            &[],
            input,
            &body(100.0),
            0.0,
            solve(),
            adaptive(1.0),
            |_| {
                gate.request();
                Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
        );
        assert!(matches!(
            cancelled,
            Err(EffusionWallTrajectoryError::Cancelled)
        ));
        assert_eq!(before, saved);
        assert_eq!(source, saved_source);
        assert_eq!(input.initial.temperature_k(), 600.0);
    });
}
