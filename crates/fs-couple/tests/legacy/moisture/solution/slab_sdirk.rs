//! G1/G3/G4 native slab stages with original finite-solution thermodynamics.
//! Sources and coefficients are authored verification fixtures, not calibration.
use super::*;
use fs_couple::moisture::adaptive::{MoistureAdvanceError, MoistureTimeControls};
use fs_couple::moisture::atmosphere::DryBasisSorption;
use fs_material::osmosis::solubility::kinetics::KineticSolubilitySolution;
use fs_qty::VolumetricFlowRate;
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::kinetics::{
        PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1,
    },
};

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn participant<'a>(
    id: u64,
    slab: &'a MoistureSlab,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> MoistureParticipant<'a> {
    MoistureParticipant {
        id,
        slab,
        faces,
        transport: config(),
    }
}

fn decay(z: f64) -> f64 {
    (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2)
}

// Independent continuous scalar ODE oracle; no production endpoint is reused.
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

fn nonideal(
    volume: f64,
    amount: f64,
    maximum: f64,
) -> (FiniteMoistureSolution, Arc<VirialSolutionModel>) {
    let model = Arc::new(
        VirialSolutionModel::new(293.15, WATER_VOLUME, 1.0, vec![maximum], vec![0.002]).unwrap(),
    );
    let bath = FiniteMoistureSolution::new(
        AqueousSolution::nonideal(model.clone(), volume, vec![amount]).unwrap(),
        WATER_MOLAR_MASS,
    )
    .unwrap();
    (bath, model)
}

fn activity(volume: f64, amount: f64) -> f64 {
    let c = amount / volume;
    (-WATER_VOLUME * (c + 0.001 * c * c)).exp()
}

#[test]
fn g1_solution_slab_sdirk_matches_pure_solvent_and_algebraic_phase_buffers() {
    with_cx(|cx| {
        let volume = 1e-4;
        let dt = 2.0;
        let slope = 0.2;
        let rate = 0.04;
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        for saturation in [0.0, 1000.0] {
            let total = 2.0 * saturation * volume;
            let bath = if saturation == 0.0 {
                reservoir(volume, Vec::new())
            } else {
                soluble_reservoir(volume, total, saturation)
            };
            let original_solution = bath.solution().clone();
            let original_phases = bath.solubility().cloned();
            let activity = (-WATER_VOLUME * saturation).exp();
            let equilibrium = slope * activity;
            for (m0, kirchhoff) in [(0.1, false), (0.27, true)] {
                let mut original = slab(card([0.01; 2], 0.3), 1.0, vec![m0], cx);
                if kirchhoff {
                    original = original.with_kirchhoff_diffusivity(nonlinear()).unwrap();
                }
                let dry = original.dry_mass_kg();
                let first = equilibrium + (m0 - equilibrium) / (1.0 + GAMMA * dt * rate);
                let expected = equilibrium + (m0 - equilibrium) * decay(dt * rate);
                let water = dry * (expected - m0);
                let final_volume = volume - water / WATER_DENSITY;
                let accepted = bath
                    .advance_sdirk2(
                        &original,
                        dt,
                        &faces,
                        config(),
                        nonlinear(),
                        controls(volume),
                        cx,
                    )
                    .unwrap();
                let report = accepted.step();
                near(accepted.moisture().moisture()[0], expected, 3e-13);
                near(
                    accepted.reservoir().solution().solvent_volume_m3(),
                    final_volume,
                    2e-16,
                );
                near(
                    report.transport.first_stage.stored_water_kg[1] / dry,
                    first,
                    3e-13,
                );
                assert_eq!(report.transport.first_stage.interval_s, [0.0, GAMMA * dt]);
                assert_eq!(report.transport.endpoint.interval_s, [0.0, dt]);
                near(report.water_to_solid_kg, water, 2e-13);
                near(
                    report.solvent_mass_kg[0] - report.solvent_mass_kg[1],
                    water,
                    2e-13,
                );
                for face in report.transport.endpoint.boundary_water_kg {
                    near(face, water / 2.0, 2e-13);
                }
                for actual in report.water_activity {
                    near(actual, activity, 2e-14);
                }
                near(
                    accepted.reservoir().solvent_mass_kg() + accepted.moisture().water_mass_kg(),
                    bath.solvent_mass_kg() + original.water_mass_kg(),
                    2e-13,
                );
                assert_eq!(
                    report.transport.endpoint.kirchhoff_source.is_some(),
                    kirchhoff
                );
                let ordinary = accepted
                    .reservoir()
                    .advance(
                        accepted.moisture(),
                        dt,
                        &faces,
                        config(),
                        controls(volume),
                        cx,
                    )
                    .unwrap();
                assert_eq!(ordinary.step().transport.nonlinear.is_some(), kirchhoff);
                if let Some(phases) = accepted.reservoir().solubility() {
                    let phase = phases.phase_partitions()[0];
                    near(phase.dissolved.value(), saturation * final_volume, 2e-13);
                    near(
                        phase.solid.value(),
                        total - saturation * final_volume,
                        2e-13,
                    );
                    near(phase.dissolved.value() + phase.solid.value(), total, 2e-16);
                    assert_eq!(
                        phases.components(),
                        original_phases.as_ref().unwrap().components()
                    );
                    assert!(phase.solid.value() > 0.0);
                } else {
                    assert_eq!(report.water_activity, [1.0; 2]);
                    assert_eq!(report.iterations, 0);
                }
                assert_eq!(original.moisture(), &[m0]);
                assert_eq!(original.time_s(), 0.0);
            }
            assert_eq!(bath.solution(), &original_solution);
            assert_eq!(bath.solubility(), original_phases.as_ref());
        }
    });
}

#[test]
fn g1_solution_slab_sdirk_resolves_nonideal_feedback_with_second_order_accuracy() {
    with_cx(|cx| {
        let volume = 1e-4;
        let amount = 0.1;
        let duration = 8.0;
        let slope = 0.25;
        let (bath, model) = nonideal(volume, amount, 3000.0);
        let before = bath.solution().clone();
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        for m0 in [0.1, 0.28] {
            let initial = slab(card([0.01; 2], 0.3), 1.0, vec![m0], cx)
                .with_kirchhoff_diffusivity(nonlinear())
                .unwrap();
            let dry = initial.dry_mass_kg();
            // Eliminate solvent storage: V(M)=V0-m_dry*(M-M0)/rho_water.
            // The original virial chart gives -ln(a_w)/V_water=c+B*c²/2.
            let rhs = |m: f64| {
                let volume = volume - dry * (m - m0) / WATER_DENSITY;
                0.04 * (slope * activity(volume, amount) - m)
            };
            let exact = rk4(m0, duration, 2048, rhs);
            near(exact, rk4(m0, duration, 4096, rhs), 5e-13);
            let errors = [4usize, 8, 16].map(|steps| {
                let dt = duration / steps as f64;
                let mut reservoir = bath.clone();
                let mut moisture = initial.clone();
                for _ in 0..steps {
                    let (next_bath, next, report) = reservoir
                        .advance_sdirk2(
                            &moisture,
                            dt,
                            &faces,
                            config(),
                            nonlinear(),
                            controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    near(
                        report.transport.endpoint.boundary_water_kg.iter().sum(),
                        report.water_to_solid_kg,
                        1e-11,
                    );
                    near(
                        report.solvent_mass_kg[0] - report.solvent_mass_kg[1],
                        report.water_to_solid_kg,
                        2e-13,
                    );
                    let final_volume = next_bath.solution().solvent_volume_m3();
                    let aw = activity(final_volume, amount);
                    near(report.water_activity[1], aw, 2e-12);
                    let c = amount / final_volume;
                    near(
                        report.osmotic_pressure_pa[1],
                        R * 293.15 * (c + 0.001 * c * c),
                        2e-6,
                    );
                    for source in report.transport.endpoint.surface_sorption.iter().flatten() {
                        near(source.query().axes()["relative_humidity"], aw, 2e-12);
                    }
                    assert!(std::ptr::eq(
                        next_bath.solution().nonideal_model().unwrap(),
                        model.as_ref()
                    ));
                    assert_eq!(next_bath.solution().solute_moles(), &[amount]);
                    assert!(report.transport.endpoint.kirchhoff_source.is_some());
                    reservoir = next_bath;
                    moisture = next;
                }
                let expected_volume = volume - dry * (moisture.moisture()[0] - m0) / WATER_DENSITY;
                near(
                    reservoir.solution().solvent_volume_m3(),
                    expected_volume,
                    3e-15,
                );
                near(
                    reservoir.solvent_mass_kg() + moisture.water_mass_kg(),
                    bath.solvent_mass_kg() + initial.water_mass_kg(),
                    3e-12,
                );
                assert_eq!(moisture.moisture()[0] > m0, m0 < 0.25);
                assert_eq!(
                    reservoir.solution().water_activity().unwrap()
                        < before.water_activity().unwrap(),
                    m0 < 0.25
                );
                (moisture.moisture()[0] - exact).abs()
            });
            for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
                assert!(ratio > 3.8 && ratio < 4.5, "{errors:?}");
            }
            assert_eq!(initial.moisture(), &[m0]);
            assert_eq!(initial.time_s(), 0.0);
        }
        assert_eq!(bath.solution(), &before);
    });
}

#[test]
fn g1_shared_solution_slab_sdirk_keeps_original_nonideal_stores_at_both_stages() {
    with_cx(|cx| {
        let volume = 1e-4;
        let amount = 0.1;
        let (bath, model) = nonideal(volume, amount, 3000.0);
        let before_bath = bath.solution().clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let b = slab(card([0.03; 2], 0.3), 0.7, vec![0.28], cx)
            .with_kirchhoff_diffusivity(nonlinear())
            .unwrap();
        let faces_a = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let faces_b = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let fields = [participant(7, &a, &faces_a), participant(3, &b, &faces_b)];
        let before = [0.05, 0.28];
        let dry = [a.dry_mass_kg(), b.dry_mass_kg()];
        let slope = [0.2, 0.25];
        let rate = [0.04, 4.0 * 0.03 / 0.7_f64.powi(2)];
        let dt = 1.0;
        let alpha = rate.map(|rate| GAMMA * dt * rate);
        // Independent scalar elimination of two distinct one-cell slab laws.
        // Bisection is used only for this analytical stage equation.
        let stage = |known: [f64; 2]| {
            let cells = |v: f64| -> [f64; 2] {
                std::array::from_fn(|i| {
                    (before[i] + known[i] + alpha[i] * slope[i] * activity(v, amount))
                        / (1.0 + alpha[i])
                })
            };
            let (mut low, mut high) = (0.5 * volume, 1.5 * volume);
            for _ in 0..100 {
                let middle = low + 0.5 * (high - low);
                let residual = WATER_DENSITY * (middle - volume)
                    + (0..2)
                        .map(|i| dry[i] * (cells(middle)[i] - before[i]))
                        .sum::<f64>();
                if residual > 0.0 {
                    high = middle;
                } else {
                    low = middle;
                }
            }
            let v = low + 0.5 * (high - low);
            (v, cells(v))
        };
        let (first_volume, first) = stage([0.0; 2]);
        let first_activity = activity(first_volume, amount);
        let known = std::array::from_fn(|i| {
            (1.0 - GAMMA) * dt * rate[i] * (slope[i] * first_activity - first[i])
        });
        let (expected_volume, expected) = stage(known);
        let accepted = bath
            .advance_shared_sdirk2(&fields, dt, 2, nonlinear(), controls(volume), cx)
            .unwrap();
        let reversed = bath
            .advance_shared_sdirk2(
                &[participant(3, &b, &faces_b), participant(7, &a, &faces_a)],
                dt,
                2,
                nonlinear(),
                controls(volume),
                cx,
            )
            .unwrap();
        near(
            accepted.reservoir().solution().solvent_volume_m3(),
            expected_volume,
            3e-16,
        );
        near(
            accepted.step().water_activity[1],
            activity(expected_volume, amount),
            2e-12,
        );
        assert_eq!(
            accepted.reservoir().solution(),
            reversed.reservoir().solution()
        );
        assert_eq!(
            accepted.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![3, 7]
        );
        let mut water = 0.0;
        for (field, replay) in accepted.solids().iter().zip(reversed.solids()) {
            let i = if field.id == 7 { 0 } else { 1 };
            assert_eq!(field.moisture.moisture(), replay.moisture.moisture());
            near(field.moisture.moisture()[0], expected[i], 3e-12);
            near(
                field.transport.first_stage.stored_water_kg[1] / dry[i],
                first[i],
                3e-12,
            );
            near(
                field.transport.first_stage.reservoir_moisture[0].unwrap(),
                slope[i] * first_activity,
                2e-12,
            );
            near(
                field.transport.endpoint.reservoir_moisture[0].unwrap(),
                slope[i] * activity(expected_volume, amount),
                2e-12,
            );
            assert_eq!(field.transport.endpoint.interval_s, [0.0, dt]);
            let change = dry[i] * (expected[i] - before[i]);
            assert_eq!(change.is_sign_positive(), i == 0);
            for face in field.transport.endpoint.boundary_water_kg {
                near(face, change / 2.0, 2e-13);
            }
            water += change;
            assert_eq!(field.transport.endpoint.kirchhoff_source.is_some(), i == 1);
        }
        near(accepted.step().water_to_solids_kg, water, 2e-13);
        near(
            accepted.step().solvent_mass_kg[0] - accepted.step().solvent_mass_kg[1],
            water,
            2e-13,
        );
        near(accepted.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(accepted.reservoir().solution().solute_moles(), &[amount]);
        assert!(std::ptr::eq(
            accepted.reservoir().solution().nonideal_model().unwrap(),
            model.as_ref()
        ));
        assert_eq!(bath.solution(), &before_bath);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.28]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g4_solution_slab_sdirk_refuses_kinetics_exhaustion_sources_and_cancellation_atomically() {
    with_cx(|cx| {
        let volume = 1e-4;
        let law = IdealSolubilityV1::new(
            SpeciesId::new("water").unwrap(),
            SpeciesId::new("authored-neutral-solute").unwrap(),
            Temperature::new(293.15),
            AmountConcentration::new(1000.0),
            0.0,
            [Temperature::new(293.15), Temperature::new(293.15)],
        )
        .unwrap();
        let kinetic = FiniteMoistureSolution::with_kinetic_solubility(
            KineticSolubilitySolution::new(
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
                    SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.2))
                        .unwrap(),
                )],
            )
            .unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        let before_phase = kinetic.kinetic_solubility().unwrap().clone();
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        for refusal in [
            kinetic
                .advance_sdirk2(
                    &original,
                    1.0,
                    &faces,
                    config(),
                    nonlinear(),
                    controls(volume),
                    cx,
                )
                .map(|_| ()),
            kinetic
                .advance_shared_sdirk2(
                    &[participant(7, &original, &faces)],
                    1.0,
                    1,
                    nonlinear(),
                    controls(volume),
                    cx,
                )
                .map(|_| ()),
        ] {
            assert!(
                matches!(refusal, Err(SolutionMoistureError::Invalid(message))
                if message == "second-order solution coupling requires algebraic phase partitioning")
            );
        }
        let time = MoistureTimeControls {
            initial_step_s: 1.0,
            minimum_step_s: 1e-6,
            maximum_step_s: 1.0,
            maximum_evaluations: 100,
            moisture_absolute_tolerance: 1e-5,
            activity_absolute_tolerance: 1e-5,
            reservoir_water_absolute_tolerance_kg: 1e-7,
            relative_tolerance: 0.0,
            water_balance_tolerance_kg: 1e-10,
        };
        assert!(matches!(
            kinetic.advance_sdirk2_adaptive_to(
                &original,
                1.0,
                &faces,
                config(),
                nonlinear(),
                controls(volume),
                time,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                evaluations: 2,
                error: SolutionMoistureError::Invalid(_),
                ..
            })
        ));
        assert_eq!(kinetic.kinetic_solubility(), Some(&before_phase));

        let small = reservoir(1e-6, Vec::new());
        let before_small = small.solution().clone();
        assert!(matches!(
            small.advance_sdirk2(
                &original,
                1.0,
                &faces,
                config(),
                nonlinear(),
                controls(1e-6),
                cx,
            ),
            Err(SolutionMoistureError::Unbracketed { .. })
        ));
        assert!(
            small
                .advance_sdirk2(
                    &original,
                    0.01,
                    &faces,
                    config(),
                    nonlinear(),
                    controls(1e-6),
                    cx,
                )
                .is_ok()
        );
        assert_eq!(small.solution(), &before_small);

        let (bounded, _) = nonideal(volume, 0.1, 1100.0);
        let before_bounded = bounded.solution().clone();
        assert!(matches!(
            bounded.advance_sdirk2(
                &original,
                0.1,
                &faces,
                config(),
                nonlinear(),
                controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Solution(OsmosisError::Nonideal(
                NonidealSolutionError::OutsideDomain { .. }
            )))
        ));
        assert_eq!(bounded.solution(), &before_bounded);
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
                bounded.advance_shared_sdirk2(
                    &[participant(7, &original, &faces)],
                    0.1,
                    1,
                    nonlinear(),
                    controls(volume),
                    &cancelled,
                ),
                Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
            ));
        });
        assert_eq!(bounded.solution(), &before_bounded);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
