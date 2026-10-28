//! Constant-pressure chemistry through the original material and Newton owners.
//! Synthetic sources test equations, not empirical combustion accuracy.
#[path = "material_exchange/gas.rs"]
#[allow(dead_code)]
mod source;

use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionError, GasReactionSolve};
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::{GasEnergySolve, parcel::ClosedGasParcel}};
use fs_matdb::QueryPoint;
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::*;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn channels(columns: &[(&str, [i128; 3], f64)], atoms: [i128; 3], formation: [f64; 3]) -> Vec<GasReaction> {
    let components = ["A", "B", "C"].into_iter().enumerate()
        .map(|(i, name)| source::component(name, 0.028 * atoms[i] as f64, formation[i]))
        .collect::<Vec<_>>();
    let ids = components.iter().map(|(s, _)| s.species().clone()).collect::<Vec<_>>();
    let matrix = StoichiometricMatrix::new(ids.clone(),
        columns.iter().map(|(id, _, _)| ReactionId::new(*id).unwrap()).collect(),
        (0..3).map(|i| columns.iter().map(|(_, nu, _)| nu[i]).collect()).collect()).unwrap();
    let elements = ElementalMatrix::new(vec![ElementId::new("X").unwrap()], ids.clone(),
        vec![atoms.to_vec()]).unwrap();
    let certificate = verify_conservation(&elements, &matrix,
        &ChargeVector::new(ids, vec![0; 3]).unwrap()).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        components.iter().map(|(s, _)| s.clone()).collect(),
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5]).unwrap()).unwrap(),
        components.iter().map(|(_, t)| t.clone()).collect(), None).unwrap();
    columns.iter().map(|(name, nu, scale)| {
        let equilibrium = IdealGasReactionEquilibriumV1::new(matrix.clone(), certificate,
            ReactionId::new(*name).unwrap(), components.iter().enumerate()
                .filter(|(i, _)| nu[*i] != 0).map(|(_, (s, _))| s.clone()).collect()).unwrap();
        GasReaction::new(catalog.clone(), DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(equilibrium,
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(*scale)).unwrap()),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap()
    }).collect()
}
fn policy(law: &GasReaction) -> ImplicitMechanismSolve {
    ImplicitMechanismSolve {
        extent_brackets_mol: law.channels().map(|(k, _)|
            (k.reverse_rate_closure().equilibrium().reaction().clone(), [-1.0, 1.0])).collect(),
        extent_tolerance_mol: 1e-11, interior_fraction: 0.05, jacobian_relative_step: 1e-5,
        maximum_iterations: 32, maximum_evaluations: 1024,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-12, relative_tolerance: 1e-12,
            linear_restart: 8, max_linear_cycles: 8, forcing_maximum: 1e-6,
            ..NewtonKrylovConfig::default() },
    }
}
fn implicit(law: GasReaction) -> GasReaction {
    let policy = policy(&law);
    law.with_implicit_mechanism(policy).unwrap()
}
fn controls() -> GasReactionSolve {
    GasReactionSolve { gas: GasEnergySolve { minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0, energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 },
        maximum_consumed_fraction: 0.99999, species_tolerance_mol: 1e-11,
        entropy_tolerance_j_k: 1e-9 }
}
fn gas(law: &GasReaction, fractions: [f64; 3], t: f64) -> ClosedGasParcel {
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0).unwrap()
        .with("specimen-lot", 42.0).unwrap();
    ClosedGasParcel::new(law.catalog().with_composition(
        Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap()).unwrap(),
        point, 0.001).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a - b).abs() <= tolerance, "{a} != {b}"); }

#[test]
fn stiff_chain_matches_independent_backward_euler_and_canonical_order() {
    let k = 10000.0;
    let dt = 1.0;
    let laws = channels(&[("ab", [-1, 1, 0], k), ("bc", [0, -1, 1], k)], [1; 3], [0.0; 3]);
    let mechanism = implicit(GasReaction::combine(laws.clone()).unwrap());
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let result = mechanism.isobaric_implicit_candidate(&old, dt, 0.0, controls(), || false).unwrap();
    let reaction = &result.step().boundary.reaction;
    let total = old.mass_kg() / 0.028;
    let h = dt * k * R * 400.0 / 100000.0;
    let denominator = 1.0 + 4.0 * h + 3.0 * h * h;
    let expected = [total * (1.0 + 3.0 * h + h * h) / denominator,
        total * h * (1.0 + h) / denominator, total * h * h / denominator];
    for (row, expected) in reaction.species.iter().zip(expected) {
        near(row.final_mol, expected, 2e-11);
        assert!(row.final_mol > 0.0);
    }
    near(result.gas().state().state().temperature, 400.0, 1e-7);
    assert_eq!(result.gas().state().state().pressure, 100000.0);
    near(result.gas().volume_m3(), old.volume_m3(), 1e-12);
    let solve = reaction.implicit_mechanism.as_ref().unwrap();
    assert_eq!(solve.equations.len(), 2);
    assert!(solve.maximum_residual_mol <= policy(&mechanism).extent_tolerance_mol);
    assert!(solve.evaluations <= policy(&mechanism).maximum_evaluations);
    assert!(reaction.entropy_change_j_k > 0.0);
    let reversed = implicit(GasReaction::combine(laws.into_iter().rev().collect()).unwrap());
    assert_eq!(result, reversed.isobaric_implicit_candidate(&old, dt, 0.0, controls(), || false).unwrap());
    let rigid = mechanism.parcel_candidate(&old, dt, controls(), || false).unwrap();
    for (a, b) in reaction.species.iter().zip(&rigid.step().species) {
        near(a.final_mol, b.final_mol, 2e-11);
    }
}

#[test]
fn exothermic_isomer_rates_use_final_volume_and_only_external_heat() {
    let release = 10000.0;
    let law = implicit(channels(&[("ab", [-1, 1, 0], 100.0)], [1; 3],
        [0.0, -release, 0.0]).remove(0));
    let old = gas(&law, [0.75, 0.0, 0.25], 400.0);
    let saved = old.clone();
    let total = old.mass_kg() / 0.028;
    for heat in [-1.0, 0.0, 1.0] {
        let dt = 0.1;
        let result = law.isobaric_implicit_candidate(&old, dt, heat, controls(), || false).unwrap();
        let step = result.step();
        let reaction = &step.boundary.reaction;
        let extent = reaction.extent_mol;
        let expected_t = 400.0 + (release * extent + heat) / (total * 3.5 * R);
        let end = result.gas();
        near(end.state().state().temperature, expected_t, 1e-7);
        assert_eq!(end.mass_kg(), old.mass_kg());
        assert_eq!(end.state().state().pressure, 100000.0);
        near(end.volume_m3(), old.volume_m3() * expected_t / 400.0, 1e-12);
        let rate = reaction.rate.net_progress_rate().unwrap().value();
        near(extent, dt * end.volume_m3() * rate, policy(&law).extent_tolerance_mol);
        assert!((extent - dt * old.volume_m3() * rate).abs() > 1e-7,
            "fixed initial volume must not satisfy an expanding reactor's endpoint rate");
        let work = 100000.0 * (end.volume_m3() - old.volume_m3());
        assert_eq!(step.boundary.boundary.work_by_gas_j, work);
        near(reaction.stored_energy_change_j + work, heat, old.mass_kg() * controls().gas.energy_tolerance_j_kg);
        near(step.enthalpy_change_j, heat, old.mass_kg() * controls().gas.energy_tolerance_j_kg);
        assert!(step.boundary.entropy_residual_j_k >= -controls().entropy_tolerance_j_k);
        assert_eq!(end.point().axes().get("specimen-lot"), Some(&42.0));
        end.state().state().validate_coherence().unwrap();
    }
    assert_eq!(old, saved);
}

#[test]
fn two_to_one_reaction_retains_molecular_amount_and_pressure_work() {
    let release = 70000.0;
    let law = implicit(channels(&[("dimer", [-2, 1, 0], 10.0)], [1, 2, 1],
        [0.0, -release, 0.0]).remove(0));
    let old = gas(&law, [0.75, 0.0, 0.25], 400.0);
    let n0 = old.mass_kg() / 0.028;
    let result = law.isobaric_implicit_candidate(&old, 0.1, 0.0, controls(), || false).unwrap();
    let step = &result.step().boundary.reaction;
    let xi = step.extent_mol;
    assert!(xi > 0.0);
    for (row, expected) in step.species.iter().zip([0.75 * n0 - 2.0 * xi, xi, 0.25 * n0]) {
        near(row.final_mol, expected, controls().species_tolerance_mol);
    }
    let t = (n0 * 3.5 * R * 400.0 + release * xi) / ((n0 - xi) * 3.5 * R);
    near(result.gas().state().state().temperature, t, 1e-7);
    near(result.gas().volume_m3(), (n0 - xi) * R * t / 100000.0, 1e-12);
    near(step.species[0].final_mol + 2.0 * step.species[1].final_mol + step.species[2].final_mol,
        n0, 3.0 * controls().species_tolerance_mol);
    assert_eq!(result.gas().state().state().pressure, 100000.0);
    assert!(result.step().enthalpy_residual_j.abs() <= old.mass_kg() * controls().gas.energy_tolerance_j_kg);
}

#[test]
fn constant_pressure_chemical_trajectory_is_first_order() {
    // Equal isomers give the independent exact decay of n_A-n_B. This tests
    // the shared implicit trajectory, not an inferred unconditional theorem.
    let k = 20.0;
    let law = implicit(channels(&[("ab", [-1, 1, 0], k)], [1; 3], [0.0; 3]).remove(0));
    let initial = gas(&law, [1.0, 0.0, 0.0], 400.0);
    let n0 = initial.mass_kg() / 0.028;
    let exact_a = 0.5 * n0 * (1.0 + (-2.0 * k * R * 400.0 / 100000.0).exp());
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let mut state = initial.clone();
        let mut amount_a = n0;
        for _ in 0..count {
            let step = law.isobaric_implicit_candidate(&state, 1.0 / f64::from(count),
                0.0, controls(), || false).unwrap();
            amount_a = step.step().boundary.reaction.species[0].final_mol;
            state = step.into_parts().0;
        }
        errors.push((amount_a - exact_a).abs());
        assert_eq!(state.mass_kg(), initial.mass_kg());
        assert_eq!(state.state().state().pressure, 100000.0);
    }
    for pair in errors.windows(2) { assert!((1.8..2.2).contains(&(pair[0] / pair[1]))); }
}

#[test]
fn missing_policies_work_limits_and_cancellation_leave_inputs_unchanged() {
    let explicit = channels(&[("ab", [-1, 1, 0], 100.0)], [1; 3], [0.0, -10000.0, 0.0]).remove(0);
    let old = gas(&explicit, [1.0, 0.0, 0.0], 400.0);
    let saved = old.clone();
    assert!(explicit.isobaric_implicit_candidate(&old, 0.1, 0.0, controls(), || false).is_err());
    let mut limited = policy(&explicit);
    limited.maximum_evaluations = 1;
    let limited = explicit.clone().with_implicit_mechanism(limited).unwrap();
    assert!(matches!(limited.isobaric_implicit_candidate(&old, 0.1, 0.0, controls(), || false),
        Err(GasReactionError::ImplicitEvaluationLimit { .. })));
    let law = implicit(explicit);
    for (dt, heat) in [(0.0, 0.0), (0.1, f64::NAN), (0.1, f64::INFINITY), (0.1, 1e9)] {
        assert!(law.isobaric_implicit_candidate(&old, dt, heat, controls(), || false).is_err());
        assert_eq!(old, saved);
    }
    assert!(matches!(law.isobaric_implicit_candidate(&old, 0.1, 0.0,
        GasReactionSolve { maximum_consumed_fraction: 0.01, ..controls() }, || false),
        Err(GasReactionError::ConsumptionLimit { .. })));
    let mut count = 0;
    law.isobaric_implicit_candidate(&old, 0.1, 0.0, controls(), || { count += 1; false }).unwrap();
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(matches!(law.isobaric_implicit_candidate(&old, 0.1, 0.0, controls(), || {
            polls += 1; polls == stop
        }), Err(GasReactionError::Cancelled)));
        assert_eq!(old, saved);
    }
}
