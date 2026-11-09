use super::*;
use crate::chemistry::{IsochoricKineticsMethod, IsochoricKineticsSolve};
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::{CaloricDomain, MixtureGas};
use fs_qty::Temperature;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn reactor(heat: f64, reference: f64, method: IsochoricKineticsMethod) -> CellChemistry {
    let original = bank(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", R, 3.5, 0.0, reference - heat),
        source("absent", R, 4.0, 0.1, reference),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("tank-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(4.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        original.clone(),
        "finite-tank-isomers",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let gas = MixtureGas::new(
        original,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    CellChemistry::new(
        gas,
        mechanism,
        IsochoricKineticsSolve {
            method,
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            energy_tolerance_j: 1e-10,
            relative_tolerance: 2e-7,
            absolute_amount_tolerance_mol: 1e-14,
            entropy_tolerance_j_k: 1e-11,
            initial_step_s: 1e-4,
            maximum_step_s: 1.0,
            maximum_steps: 30000,
            ..IsochoricKineticsSolve::default()
        },
        || false,
    )
    .unwrap()
}

fn tank(chemistry: &CellChemistry) -> Plenum {
    Plenum::new(
        0.2,
        chemistry
            .gas()
            .conserved(&[1.0, 0.0, 0.0], 0.0, 1.0)
            .unwrap(),
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_compact_reactions_match_independent_reversible_decay_without_pressure_heating() {
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        for reference in [0.0, -1000.0, 1000.0] {
            let mut chemistry = reactor(0.0, reference, method);
            let original = tank(&chemistry);
            let mut candidate = original.clone();
            let bytes = chemistry.workspace_bytes();
            let buffer = candidate.state().partial_densities.as_ptr();
            let receipt = original
                .react_candidate_into(&mut chemistry, Time::new(0.03), 17, &mut candidate, || {
                    false
                })
                .unwrap();
            let exact_b = 0.5 * (1.0 - (-8.0_f64 * 0.03).exp());
            close(candidate.state().partial_densities[1], exact_b, 6e-5);
            close(candidate.state().partial_densities[0], 1.0 - exact_b, 6e-5);
            assert_eq!(candidate.state().partial_densities[2], 0.0);
            let p = chemistry
                .gas()
                .primitive(candidate.state(), || false)
                .unwrap();
            close(p.temperature, 1.0, 2e-10);
            close(p.pressure, 1.0, 2e-10);
            assert_eq!(candidate.state().momentum, 0.0);
            assert_eq!(candidate.volume_m3(), original.volume_m3());
            close(
                receipt.species_change_kg[1],
                candidate.state().partial_densities[1] * 0.2,
                1e-14,
            );
            close(receipt.referenced_energy_change_j, 0.0, 2e-12);
            assert_eq!(receipt.cell_solves, 1);
            assert!(receipt.entropy_change_j_k > 0.0);
            assert_eq!(candidate.state().partial_densities.as_ptr(), buffer);
            assert_eq!(chemistry.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g0_compact_chemical_energy_drives_pressure_once_under_formation_gauge_shifts() {
    let mut reference_state: Option<Vec<f64>> = None;
    for reference in [0.0, -1000.0, 1000.0] {
        let mut chemistry = reactor(8.0, reference, IsochoricKineticsMethod::Sdirk2);
        let original = tank(&chemistry);
        let mut candidate = original.clone();
        let receipt = original
            .react_candidate_into(&mut chemistry, Time::new(0.03), 2, &mut candidate, || false)
            .unwrap();
        let q = candidate.state();
        let converted = q.partial_densities[1];
        let p = chemistry.gas().primitive(q, || false).unwrap();
        assert!(converted > 0.0);
        close(p.temperature, 1.0 + 8.0 * converted / 2.5, 2e-10);
        close(p.pressure, p.temperature, 2e-11);
        let gauge_change: f64 = q
            .partial_densities
            .iter()
            .zip(original.state().partial_densities)
            .zip(chemistry.gas().energy_references())
            .map(|((after, before), energy)| (after - before) * energy)
            .sum();
        close(
            q.energy - original.state().energy + gauge_change,
            0.0,
            2e-10,
        );
        close(
            receipt.referenced_energy_change_j,
            8.0 * converted * 0.2,
            2e-10,
        );
        if let Some(previous) = &reference_state {
            for (a, b) in q.partial_densities.iter().zip(previous) {
                close(*a, *b, 3e-9);
            }
        } else {
            reference_state = Some(q.partial_densities.to_vec());
        }
    }
}

#[test]
fn g4_compact_chemical_candidates_keep_output_on_final_poll_shape_and_precision_refusal() {
    let make = || reactor(8.0, 0.0, IsochoricKineticsMethod::Sdirk2);
    let mut probe = make();
    let original = tank(&probe);
    let mut expected = original.clone();
    let mut polls = 0;
    original
        .react_candidate_into(&mut probe, Time::new(0.03), 9, &mut expected, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut chemistry = make();
    let mut output = original.clone();
    let buffer = output.state().partial_densities.as_ptr();
    let mut seen = 0;
    assert_eq!(
        original.react_candidate_into(&mut chemistry, Time::new(0.03), 9, &mut output, || {
            seen += 1;
            seen == polls
        },),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(seen, polls);
    assert_eq!(output, original);
    assert_eq!(output.state().partial_densities.as_ptr(), buffer);
    original
        .react_candidate_into(&mut chemistry, Time::new(0.03), 9, &mut output, || false)
        .unwrap();
    assert_eq!(output, expected);
    let mut wrong = Plenum::new(
        1.0,
        crate::multicomponent::Conserved {
            partial_densities: vec![1.0, 0.0],
            momentum: 0.0,
            energy: 1.0,
        },
    )
    .unwrap();
    let saved = wrong.clone();
    assert!(
        original
            .react_candidate_into(&mut chemistry, Time::new(0.03), 9, &mut wrong, || false)
            .is_err()
    );
    assert_eq!(wrong, saved);
    let mut shifted = reactor(0.0, 1e18, IsochoricKineticsMethod::Sdirk2);
    let lost_precision = tank(&shifted);
    let mut target = lost_precision.clone();
    assert!(matches!(
        lost_precision
            .react_candidate_into(&mut shifted, Time::new(0.03), 9, &mut target, || false,),
        Err(ReactiveError::EnergyConsistency { .. })
    ));
    assert_eq!(target, lost_precision);
}
