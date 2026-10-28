//! Direct coupled-chemistry behavior on the existing manufactured source set.
use super::*;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use fs_couple::cavity_phs::piston::reactive::subcycling::ReactionSubcycling;
use fs_solver::NewtonKrylovConfig;

fn policy(law: &GasReaction) -> ImplicitMechanismSolve {
    ImplicitMechanismSolve {
        extent_brackets_mol: law
            .channels()
            .map(|(k, _)| {
                (
                    k.reverse_rate_closure().equilibrium().reaction().clone(),
                    [-1.0, 1.0],
                )
            })
            .collect(),
        extent_tolerance_mol: 1e-11,
        interior_fraction: 0.05,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 1024,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_maximum: 1e-6,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn permissive() -> GasReactionSolve {
    GasReactionSolve {
        maximum_consumed_fraction: 0.99999,
        ..controls()
    }
}
fn implicit(law: GasReaction) -> GasReaction {
    let solve = policy(&law);
    law.with_implicit_mechanism(solve).unwrap()
}

#[test]
fn stiff_chain_solves_both_endpoint_rates_and_creates_then_consumes_an_intermediate() {
    let k = 10000.0;
    let dt = 1.0;
    let channels = laws(&[("ab", 0, 1, k), ("bc", 1, 2, k)], [0.0; 3]);
    let mechanism = implicit(GasReaction::combine(channels.clone()).unwrap());
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let saved = old.clone();
    // The original explicit path correctly refuses this unresolved interval.
    assert!(
        GasReaction::combine(channels.clone())
            .unwrap()
            .parcel_candidate(&old, dt, permissive(), || false)
            .is_err()
    );
    let result = mechanism
        .parcel_candidate(&old, dt, permissive(), || false)
        .unwrap();
    let total = old.mass_kg() / 0.028;
    let h = dt * k * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0 / 100000.0;
    assert!(h > 300.0);
    // Independent inverse of I - h Q for A <-> B <-> C, where all
    // standard states coincide and Q = [[-1,1,0],[1,-2,1],[0,1,-1]].
    let d = 1.0 + 4.0 * h + 3.0 * h * h;
    let expected = [
        total * (1.0 + 3.0 * h + h * h) / d,
        total * h * (1.0 + h) / d,
        total * h * h / d,
    ];
    for (sample, target) in result.step().species.iter().zip(expected) {
        near(sample.final_mol, target, 2e-11);
        assert!(sample.final_mol > 0.0);
    }
    near(result.step().extent_mol, total - expected[0], 2e-11);
    near(
        result.step().additional_reactions[0].extent_mol,
        expected[2],
        2e-11,
    );
    let solve = result.step().implicit_mechanism.as_ref().unwrap();
    assert!(solve.nonlinear.as_ref().unwrap().converged);
    assert_eq!(solve.equations.len(), 2);
    assert!(solve.maximum_residual_mol <= policy(&mechanism).extent_tolerance_mol);
    assert!(solve.evaluations <= policy(&mechanism).maximum_evaluations);
    assert!(result.step().implicit.is_none());
    assert!(result.step().entropy_change_j_k > 0.0);
    assert_eq!(result.gas().mass_kg(), old.mass_kg());
    assert_eq!(result.gas().volume_m3(), old.volume_m3());
    near(result.gas().state().state().temperature, 400.0, 1e-7);
    assert!(
        result.step().stored_energy_change_j.abs()
            <= old.mass_kg() * permissive().gas.energy_tolerance_j_kg
    );
    assert_eq!(old, saved);
    // Canonical source/channel ordering makes the complete solve replayable.
    let reversed = implicit(GasReaction::combine(channels.into_iter().rev().collect()).unwrap());
    assert_eq!(mechanism, reversed);
    assert_eq!(
        result,
        reversed
            .parcel_candidate(&old, dt, permissive(), || false)
            .unwrap()
    );
}

#[test]
fn every_material_rate_reads_one_actual_endpoint_and_formation_energy_drives_pressure_and_sound() {
    let channels = laws(
        &[("ab", 0, 1, 1.0), ("bc", 1, 2, 1.0)],
        [0.0, -3000.0, -6000.0],
    );
    let material = GasReaction::combine(
        channels
            .into_iter()
            .zip([1000.0, 2000.0])
            .map(|(c, k)| c.with_material_rate(rate_source(k)).unwrap())
            .collect(),
    )
    .unwrap();
    let mechanism = implicit(material);
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let dt = 0.01;
    let result = mechanism
        .parcel_candidate(&old, dt, permissive(), || false)
        .unwrap();
    let t = result.gas().state().state().temperature;
    assert!(t > 400.0 && t < 800.0);
    let n = old.mass_kg() / 0.028;
    let first = result.step().extent_mol;
    let second = result.step().additional_reactions[0].extent_mol;
    assert!(second > 0.0);
    near(
        t,
        400.0 + 3000.0 * (first + second) / (n * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        1e-7,
    );
    near(
        result.gas().state().state().pressure,
        100000.0 * t / 400.0,
        1e-5,
    );
    assert!(result.gas().state().state().sound_speed > old.state().state().sound_speed);
    for (sample, k) in [
        (result.step().rate_material.as_ref().unwrap(), 1000.0),
        (
            result.step().additional_reactions[0]
                .rate_material
                .as_ref()
                .unwrap(),
            2000.0,
        ),
    ] {
        near(
            sample
                .property(ACTIVITY_RATE_SCALE_PROPERTY)
                .unwrap()
                .value_si(),
            k * t / 400.0,
            1e-8,
        );
        assert!(sample.query_point().contains(&("catalyst-lot".into(), 7.0)));
        assert!(sample.query_point().contains(&("temperature".into(), t)));
    }
    for ((_, extent), rate) in result.step().reaction_extents().zip(
        std::iter::once(&result.step().rate)
            .chain(result.step().additional_reactions.iter().map(|r| &r.rate)),
    ) {
        near(
            extent,
            dt * old.volume_m3() * rate.net_progress_rate().unwrap().value(),
            1e-11,
        );
    }
    // Changing rate data, without a physics-code change, alters the coupled root.
    let constant = implicit(
        GasReaction::combine(laws(
            &[("ab", 0, 1, 1000.0), ("bc", 1, 2, 2000.0)],
            [0.0, -3000.0, -6000.0],
        ))
        .unwrap(),
    );
    let frozen = constant
        .parcel_candidate(&old, dt, permissive(), || false)
        .unwrap();
    assert!((result.step().extent_mol - frozen.step().extent_mol).abs() > 1e-7);
    assert_eq!(result.gas().point().axes().get("catalyst-lot"), Some(&7.0));
}

#[test]
fn an_unreachable_absent_submechanism_is_not_seeded_into_existence() {
    // Two distinct reaction columns, both confined to absent A and B; C is inert.
    let mechanism = implicit(
        GasReaction::combine(laws(
            &[("ab1", 0, 1, 10000.0), ("ab2", 0, 1, 5000.0)],
            [0.0; 3],
        ))
        .unwrap(),
    );
    let old = gas(&mechanism, [0.0, 0.0, 1.0], 400.0);
    let result = mechanism
        .parcel_candidate(&old, 1.0, permissive(), || false)
        .unwrap();
    assert_eq!(result.gas(), &old);
    assert!(
        result
            .step()
            .reaction_extents()
            .all(|(_, extent)| extent == 0.0)
    );
    let solve = result.step().implicit_mechanism.as_ref().unwrap();
    assert!(solve.nonlinear.is_none());
    assert_eq!(solve.evaluations, 1);
    assert_eq!(solve.maximum_residual_mol, 0.0);
}

#[test]
fn work_bounds_source_domains_and_cancellation_preserve_the_original_inventory() {
    let raw = GasReaction::combine(laws(
        &[("ab", 0, 1, 1000.0), ("bc", 1, 2, 1500.0)],
        [0.0; 3],
    ))
    .unwrap();
    let mechanism = implicit(raw.clone());
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let saved = old.clone();
    let mut polls = 0;
    mechanism
        .parcel_candidate(&old, 0.1, permissive(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            mechanism.parcel_candidate(&old, 0.1, permissive(), || {
                count += 1;
                count == stop
            }),
            Err(GasReactionError::Cancelled)
        );
        assert_eq!(old, saved);
    }
    let mut limited = policy(&raw);
    limited.maximum_evaluations = 1;
    assert!(matches!(
        raw.clone()
            .with_implicit_mechanism(limited)
            .unwrap()
            .parcel_candidate(&old, 0.1, permissive(), || false),
        Err(GasReactionError::ImplicitEvaluationLimit {
            evaluations: 1,
            maximum: 1
        })
    ));
    assert!(matches!(
        mechanism.parcel_candidate(&old, 0.1, controls(), || false),
        Err(GasReactionError::ConsumptionLimit { .. })
    ));
    let sourced = implicit(
        GasReaction::combine(
            laws(&[("ab", 0, 1, 1.0), ("bc", 1, 2, 1.0)], [0.0; 3])
                .into_iter()
                .map(|c| c.with_material_rate(rate_source(1000.0)).unwrap())
                .collect(),
        )
        .unwrap(),
    );
    let hot = gas(&sourced, [1.0, 0.0, 0.0], 900.0);
    assert!(matches!(
        sourced.parcel_candidate(&hot, 0.1, permissive(), || false),
        Err(GasReactionError::RateMaterial(_))
    ));
    assert_eq!(old, saved);
}

#[test]
fn mechanism_bounds_are_identity_bound_and_methods_cannot_be_mixed() {
    use fs_couple::cavity_phs::piston::reactive::implicit::ImplicitReactionSolve;
    let raw = GasReaction::combine(laws(
        &[("ab", 0, 1, 1000.0), ("bc", 1, 2, 1000.0)],
        [0.0; 3],
    ))
    .unwrap();
    let mut missing = policy(&raw);
    let key = ReactionId::new("ab").unwrap();
    missing.extent_brackets_mol.remove(&key);
    assert!(raw.clone().with_implicit_mechanism(missing).is_err());
    let mut foreign = policy(&raw);
    foreign.extent_brackets_mol.remove(&key);
    foreign
        .extent_brackets_mol
        .insert(ReactionId::new("foreign").unwrap(), [-1.0, 1.0]);
    assert!(raw.clone().with_implicit_mechanism(foreign).is_err());
    let mut bad = policy(&raw);
    bad.jacobian_relative_step = 0.0;
    assert!(raw.clone().with_implicit_mechanism(bad).is_err());
    let bound = implicit(raw.clone());
    assert!(bound.clone().with_implicit_mechanism(policy(&raw)).is_err());
    assert!(
        bound
            .clone()
            .with_implicit(ImplicitReactionSolve {
                extent_bracket_mol: [-1.0, 1.0],
                extent_tolerance_mol: 1e-11,
                maximum_iterations: 64
            })
            .is_err()
    );
    assert!(GasReaction::combine(vec![bound, raw]).is_err());
}

#[test]
fn original_subcycling_keeps_actual_coupled_equations_and_matches_manual_half_steps() {
    let mechanism = implicit(
        GasReaction::combine(laws(&[("ab", 0, 1, 100.0), ("bc", 1, 2, 100.0)], [0.0; 3])).unwrap(),
    );
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let dt = 0.1;
    let policy = ReactionSubcycling {
        absolute_amount_mol: 1.0,
        absolute_temperature_k: 1.0,
        absolute_pressure_pa: 100.0,
        relative_tolerance: 0.0,
        maximum_substeps: 2,
        maximum_evaluations: 3,
    };
    let refined = mechanism
        .clone()
        .with_subcycling(policy)
        .unwrap()
        .parcel_candidate(&old, dt, permissive(), || false)
        .unwrap();
    let first = mechanism
        .parcel_candidate(&old, 0.5 * dt, permissive(), || false)
        .unwrap();
    let second = mechanism
        .parcel_candidate(first.gas(), 0.5 * dt, permissive(), || false)
        .unwrap();
    assert_eq!(refined.gas(), second.gas());
    assert!(refined.step().implicit_mechanism.is_none());
    let substeps = &refined.step().subcycling.as_ref().unwrap().steps;
    assert_eq!(substeps.len(), 2);
    assert_eq!(&substeps[0].reaction, first.step());
    assert_eq!(&substeps[1].reaction, second.step());
    assert!(
        substeps
            .iter()
            .all(|s| s.reaction.implicit_mechanism.is_some())
    );
    // Independent exact relaxation of the same three-isomer linear system:
    // backward Euler refines toward, but does not equal, its exponential flow.
    let h = 100.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0 / 100000.0 * dt;
    let n = old.mass_kg() / 0.028;
    let exact_c = n * (1.0 / 3.0 - 0.5 * (-h).exp() + (-3.0 * h).exp() / 6.0);
    let coarse = mechanism
        .parcel_candidate(&old, dt, permissive(), || false)
        .unwrap();
    assert!(
        (refined.step().species[2].final_mol - exact_c).abs()
            < (coarse.step().species[2].final_mol - exact_c).abs()
    );
}

#[test]
fn original_piston_uses_the_solved_secondary_reaction_without_an_extra_heat_credit() {
    let mechanism = implicit(
        GasReaction::combine(laws(
            &[("ab", 0, 1, 1000.0), ("bc", 1, 2, 2000.0)],
            [0.0, -3000.0, -6000.0],
        ))
        .unwrap(),
    );
    // Resolve the native piston/compressibility clock over the fixture bounds.
    let dt = 0.0001;
    let piston = GasPiston::new(gas(&mechanism, [1.0, 0.0, 0.0], 400.0), 1.0, 0.01, 0.0).unwrap();
    let ch = source::chamber(10, &piston, dt);
    let reference = mechanism
        .parcel_candidate(piston.gas(), dt, permissive(), || false)
        .unwrap();
    let result = piston
        .reactive_candidate(&mechanism, ch.input, permissive(), ch.piston_solve, || {
            false
        })
        .unwrap();
    assert_eq!(&result.step().reaction, reference.step());
    assert!(result.step().reaction.additional_reactions[0].extent_mol > 0.0);
    assert!(result.step().reaction.implicit_mechanism.is_some());
    assert!(result.step().total_energy_residual_j.abs() <= ch.piston_solve.energy_tolerance_j);
    assert_eq!(piston.gas().state().state().temperature, 400.0);
}
