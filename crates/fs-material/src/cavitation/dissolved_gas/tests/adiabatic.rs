//! Synthetic caloric/solubility fixtures; no empirical gas/solvent accuracy claim.
use super::*;
use crate::cavitation::dissolved_gas::adiabatic::{
    AdiabaticDissolutionControls, DissolutionSolventCalorics,
};

fn solvent(catalog: &ThermochemicalGasModel) -> DissolutionSolventCalorics {
    // Deliberately small authored thermal storage exposes the coupling. Neither
    // its mass nor its calorics is inferred from the liquid inventory's volume.
    DissolutionSolventCalorics {
        mass_kg: 0.02,
        heat_capacity: calorics(catalog).partial_enthalpy,
    }
}

fn thermal_controls() -> AdiabaticDissolutionControls {
    AdiabaticDissolutionControls {
        dissolution: controls(),
        temperature_tolerance_k: 1e-5,
        dissolved_moles_tolerance: 1e-9,
        minimum_step_s: 1e-5,
        maximum_trials: 4096,
    }
}

fn varying_solubility(catalog: &ThermochemicalGasModel) -> HenryDissolution {
    let mut h = claim(
        HENRY_SOLUBILITY_PROPERTY,
        HENRY_SOLUBILITY_DIMS,
        H,
        [1000.0, 1e7],
    );
    h.value = PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![
            (200.0, 2.0 * H),
            (300.0, H),
            (350.0, 0.5 * H),
            (800.0, 0.5 * H),
        ],
        dims: HENRY_SOLUBILITY_DIMS,
    };
    h.interpolation = InterpolationPolicy::LinearInside;
    let source = HenrySource {
        material: card(
            "synthetic temperature-dependent solubility",
            vec![
                h,
                claim(
                    DISSOLVED_DIFFUSIVITY_PROPERTY,
                    DIFFUSIVITY_DIMS,
                    D,
                    [1000.0, 1e7],
                ),
            ],
        ),
        point: context(T),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    HenryDissolution::new(catalog, id("A"), source, calorics(catalog)).unwrap()
}

#[test]
fn absorption_heats_exsolution_cools_and_all_original_stores_balance() {
    // Nonconstant NASA-9 Cp(T) prevents a constant-gamma energy surrogate.
    let catalog = model(vec![0.6, 0.4], 0.001);
    let solvent = solvent(&catalog);
    for n0 in [0.0, 0.006] {
        let mut law = varying_solubility(&catalog);
        let mut gas = gas(&catalog);
        let mut liquid = DissolvedGasInventory::new(VL, n0, 0.0).unwrap();
        let initial_gas = gas.clone();
        let initial_law = law.clone();
        let initial_u = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
            + n0 * MW * law.dissolved_enthalpy_j_kg;
        let step = law
            .advance_adiabatic(
                &mut gas,
                &mut liquid,
                &solvent,
                film(),
                200.0,
                thermal_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(liquid.time_s(), 200.0);
        assert_eq!(law.temperature_k(), gas.state().state().temperature);
        assert!((step.dissolved_moles_change > 0.0) == (n0 == 0.0));
        assert!((step.temperature_k > T) == (n0 == 0.0));
        assert!((step.temperature_k - T).abs() > 0.05);
        let solvent_change = solvent.mass_kg
            * solvent
                .heat_capacity
                .specific_enthalpy_change(T, step.temperature_k)
                .unwrap();
        let final_u = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
            + liquid.dissolved_moles() * MW * law.dissolved_enthalpy_j_kg;
        near(final_u - initial_u + solvent_change, 0.0, 3e-8);
        near(
            step.gas_internal_energy_change_j
                + step.dissolved_internal_energy_change_j
                + step.solvent_internal_energy_change_j,
            0.0,
            3e-8,
        );
        near(
            step.energy_residual_j,
            0.0,
            controls().gas.energy_tolerance_j,
        );
        near(
            gas.species_mass_kg(&id("A")).unwrap() + liquid.dissolved_moles() * MW,
            initial_gas.species_mass_kg(&id("A")).unwrap() + n0 * MW,
            1e-14,
        );
        near(
            gas.species_mass_kg(&id("B")).unwrap(),
            initial_gas.species_mass_kg(&id("B")).unwrap(),
            1e-15,
        );
        assert_eq!(law.source().material, initial_law.source().material);
        assert_eq!(law.calorics(), initial_law.calorics());
        assert_eq!(law.catalog(), initial_law.catalog());
        near(
            law.henry_mol_m3_pa(),
            H * (1.0 - 0.01 * (step.temperature_k - T)),
            1e-18,
        );
        assert!(step.estimated_temperature_error_k <= thermal_controls().temperature_tolerance_k);
        assert!(
            step.estimated_dissolved_moles_error <= thermal_controls().dissolved_moles_tolerance
        );
    }
}

#[test]
fn adaptive_thermal_feedback_converges_to_an_independent_closed_energy_ode() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let solvent = solvent(&catalog);
    let original_gas = gas(&catalog);
    let na = original_gas.species_mass_kg(&id("A")).unwrap() / MW;
    let nb = original_gas.species_mass_kg(&id("B")).unwrap() / 0.004;
    // Independent algebraic energy elimination for these constant-Cp fixtures,
    // then RK4 on dn/dt. This does not use the implementation's split update.
    let temperature = |n: f64| {
        let capacity = (na - n) * 2.5 * R + nb * 1.5 * R + n * MW * 500.0 + solvent.mass_kg * 500.0;
        T + n * (SOLUTION_HEAT - R * T) / capacity
    };
    let rate = |n: f64| {
        let t = temperature(n);
        let h = H * (1.0 - 0.01 * (t - T));
        film().area_m2 * D / film().thickness_m * (h * (na - n) * R * t / VG - n / VL)
    };
    let mut reference = 0.0;
    let dt = 0.1;
    for _ in 0..2000 {
        let k1 = rate(reference);
        let k2 = rate(reference + 0.5 * dt * k1);
        let k3 = rate(reference + 0.5 * dt * k2);
        let k4 = rate(reference + dt * k3);
        reference += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let run = |controls| {
        let mut law = varying_solubility(&catalog);
        let mut gas = original_gas.clone();
        let mut liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
        law.advance_adiabatic(
            &mut gas,
            &mut liquid,
            &solvent,
            film(),
            200.0,
            controls,
            || false,
        )
        .unwrap()
    };
    let mut loose = thermal_controls();
    loose.temperature_tolerance_k = 1e-2;
    loose.dissolved_moles_tolerance = 1e-6;
    let coarse = run(loose);
    let fine = run(thermal_controls());
    assert!(fine.accepted_steps > coarse.accepted_steps);
    assert!(
        (fine.state.dissolved_moles() - reference).abs()
            < 0.15 * (coarse.state.dissolved_moles() - reference).abs()
    );
    near(fine.state.dissolved_moles(), reference, 2e-9);
    near(
        fine.temperature_k,
        temperature(fine.state.dissolved_moles()),
        1e-8,
    );
    let mut isothermal_gas = original_gas.clone();
    let mut isothermal_liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
    varying_solubility(&catalog)
        .advance(
            &mut isothermal_gas,
            &mut isothermal_liquid,
            film(),
            200.0,
            controls(),
            || false,
        )
        .unwrap();
    assert!(fine.state.dissolved_moles() < isothermal_liquid.dissolved_moles() - 1e-8);
}

#[test]
fn thermal_refusal_work_exhaustion_and_final_cancellation_are_atomic() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let solvent = solvent(&catalog);
    let original = (
        varying_solubility(&catalog),
        gas(&catalog),
        DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap(),
    );
    let mut bracket = thermal_controls();
    bracket.dissolution.gas.gas.maximum_temperature_k = T + 1e-7;
    let mut exhausted = thermal_controls();
    exhausted.maximum_trials = 1;
    exhausted.temperature_tolerance_k = 1e-12;
    exhausted.dissolved_moles_tolerance = 1e-16;
    for rejected in [bracket, exhausted] {
        let (mut law, mut gas, mut liquid) = original.clone();
        assert!(
            law.advance_adiabatic(
                &mut gas,
                &mut liquid,
                &solvent,
                film(),
                200.0,
                rejected,
                || false,
            )
            .is_err()
        );
        assert_eq!((law, gas, liquid), original);
    }
    let (mut law, mut gas, mut liquid) = original.clone();
    let mut polls = 0;
    law.advance_adiabatic(
        &mut gas,
        &mut liquid,
        &solvent,
        film(),
        20.0,
        thermal_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let (mut law, mut gas, mut liquid) = original.clone();
    let mut calls = 0;
    let error = law
        .advance_adiabatic(
            &mut gas,
            &mut liquid,
            &solvent,
            film(),
            20.0,
            thermal_controls(),
            || {
                calls += 1;
                calls == polls
            },
        )
        .unwrap_err();
    assert_eq!(error, DissolvedGasError::Cancelled);
    assert_eq!((law, gas, liquid), original);
}
