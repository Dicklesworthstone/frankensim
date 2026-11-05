//! G1/G3: evolving original pressure rates, conservative storage and time error.
use super::super::super::tests::{close, controls as flow_controls, energy, nozzle, parcel};
use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::energy::tests::bound_with_reference;
use crate::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::Pressure;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

fn controls() -> RealGasNozzleNetworkTimeSolve {
    RealGasNozzleNetworkTimeSolve {
        initial_dt_s: 0.1,
        minimum_dt_s: 1e-7,
        maximum_dt_s: 0.1,
        maximum_trials: 4096,
        maximum_steps: 8192,
        energy_absolute_tolerance_j: 1e-3,
        mass_absolute_tolerance_kg: 1e-8,
        temperature_absolute_tolerance_k: 1e-5,
        relative_tolerance: 0.0,
        cumulative_energy_tolerance_j: 1e-4,
        cumulative_mass_tolerance_kg: 1e-11,
    }
}
fn network() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.25,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 1e-7,
    }
}
fn nodes(reference: f64) -> BTreeMap<u64, IsochoricRealGasParcel> {
    let model = bound_with_reference(reference);
    BTreeMap::from([
        (10, parcel(model.clone(), 5e6, 0.01)),
        (20, parcel(model, 4.8e6, 0.01)),
    ])
}
fn endpoints(nodes: &BTreeMap<u64, IsochoricRealGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    nodes.keys().map(|&id| (id, energy())).collect()
}
fn links() -> Vec<RealGasNozzleNetworkLink> {
    vec![RealGasNozzleNetworkLink {
        id: 9,
        parcels: [20, 10],
        nozzle: nozzle(),
        solve: flow_controls(),
    }]
}
fn evolve(
    nodes: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasNozzleNetworkLink],
    duration: f64,
    controls: RealGasNozzleNetworkTimeSolve,
) -> RealGasNozzleNetworkEvolutionCandidate {
    advance_adaptive(
        nodes,
        links,
        duration,
        &endpoints(nodes),
        network(),
        controls,
        || false,
    )
    .unwrap()
}

#[test]
fn accepts_exact_original_half_trajectory_and_accounts_no_comparison_streams() {
    let original = nodes(0.0);
    let paths = links();
    let loose = RealGasNozzleNetworkTimeSolve {
        energy_absolute_tolerance_j: 100.0,
        mass_absolute_tolerance_kg: 0.1,
        temperature_absolute_tolerance_k: 10.0,
        ..controls()
    };
    let first = advance(
        &original,
        &paths,
        0.05,
        &endpoints(&original),
        network(),
        || false,
    )
    .unwrap();
    let second = advance(
        first.parcels(),
        &paths,
        0.05,
        &endpoints(&original),
        network(),
        || false,
    )
    .unwrap();
    let actual = evolve(&original, &paths, 0.1, loose);
    assert_eq!(actual.parcels(), second.parcels());
    assert_eq!(actual.report().comparisons, 1);
    assert_eq!(actual.report().evaluated_steps, 3);
    assert_eq!(actual.report().accepted_steps, 2);
    let path = &actual.report().paths[0];
    assert_eq!(path.parcels, [10, 20]);
    assert_eq!(path.mass_right_to_left_kg, 0.0);
    let mass =
        first.step().links[0].transferred_mass_kg + second.step().links[0].transferred_mass_kg;
    close(path.mass_left_to_right_kg, mass, 1e-18);
    let h = first.step().links[0].transferred_mass_kg
        * original[&10].state().specific_enthalpy_j_kg()
        + second.step().links[0].transferred_mass_kg
            * first.parcels()[&10].state().specific_enthalpy_j_kg();
    close(path.enthalpy_left_to_right_j, h, 1e-12);
    let mut total = 0.0;
    for (&id, before) in &original {
        let after = &actual.parcels()[&id];
        let du = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        close(du, if id == 10 { -h } else { h }, 1e-7);
        total += du;
    }
    close(total, actual.report().energy_residual_j, 1e-7);
}

#[test]
fn tighter_local_error_resolves_real_gas_blowdown_entropy_and_refined_trajectory() {
    let original = nodes(0.0);
    let paths = links();
    let loose = evolve(
        &original,
        &paths,
        0.1,
        RealGasNozzleNetworkTimeSolve {
            energy_absolute_tolerance_j: 0.2,
            mass_absolute_tolerance_kg: 1e-4,
            temperature_absolute_tolerance_k: 1e-2,
            ..controls()
        },
    );
    let fine = evolve(&original, &paths, 0.1, controls());
    assert!(fine.report().accuracy_rejections > 0);
    assert!(fine.report().accepted_steps > loose.report().accepted_steps);
    assert!(fine.report().maximum_accepted_normalized_error <= 1.0);
    let entropy_defect = |end: &IsochoricRealGasParcel| {
        original[&10]
            .model()
            .specific_entropy_change_between_densities(
                original[&10].density_kg_m3(),
                300.0,
                end.density_kg_m3(),
                end.state().temperature_k(),
            )
            .unwrap()
            .abs()
    };
    // A well-mixed adiabatic donor obeys ds=0. This independent constitutive
    // identity checks the evolving enthalpy path, not just controller counters.
    assert!(entropy_defect(&fine.parcels()[&10]) < 0.5 * entropy_defect(&loose.parcels()[&10]));
    let mut reference = original.clone();
    for _ in 0..256 {
        reference = advance(
            &reference,
            &paths,
            0.1 / 256.0,
            &endpoints(&original),
            network(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    }
    let deviation = |candidate: &RealGasNozzleNetworkEvolutionCandidate| {
        candidate
            .parcels()
            .iter()
            .map(|(&id, value)| {
                (value.state().temperature_k() - reference[&id].state().temperature_k()).abs()
            })
            .fold(0.0, f64::max)
    };
    assert!(deviation(&fine) < 0.5 * deviation(&loose));
    assert!(
        fine.report().accumulated_absolute_energy_residual_j
            <= controls().cumulative_energy_tolerance_j
    );
}

fn mixture(x_a: f64, references: [f64; 2]) -> ThermochemicalPengRobinsonModel {
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, references[0]),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, references[1]),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x_a, 1.0 - x_a]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap();
    ThermochemicalPengRobinsonModel::new(
        source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap()
}
fn mixed_nodes(references: [f64; 2]) -> BTreeMap<u64, IsochoricRealGasParcel> {
    BTreeMap::from([
        (10, parcel(mixture(0.75, references), 5e6, 0.01)),
        (20, parcel(mixture(0.35, references), 4.8e6, 0.01)),
    ])
}

#[test]
fn changing_mixtures_keep_time_estimates_under_independent_huge_formation_shifts() {
    let base = mixed_nodes([0.0; 2]);
    let shifted = mixed_nodes([1e24, -2e24]);
    let paths = links();
    let a = evolve(&base, &paths, 0.05, controls());
    let b = evolve(&shifted, &paths, 0.05, controls());
    assert_eq!(a.report().comparisons, b.report().comparisons);
    assert_eq!(a.report().accepted_steps, b.report().accepted_steps);
    assert_eq!(a.report().last_error, b.report().last_error);
    assert_eq!(a.report().suggested_dt_s, b.report().suggested_dt_s);
    for id in base.keys() {
        let x = &a.parcels()[id];
        let y = &b.parcels()[id];
        assert_eq!(x.mass_kg(), y.mass_kg());
        assert_eq!(x.state().temperature_k(), y.state().temperature_k());
        assert_eq!(x.state().pressure_pa(), y.state().pressure_pa());
        assert_eq!(
            x.model().mechanics().mole_fractions(),
            y.model().mechanics().mole_fractions()
        );
    }
    assert_ne!(
        a.parcels()[&20].model().mechanics().mole_fractions(),
        base[&20].model().mechanics().mole_fractions()
    );
    let coarse = advance(&base, &paths, 0.05, &endpoints(&base), network(), || false).unwrap();
    let first = advance(&base, &paths, 0.025, &endpoints(&base), network(), || false).unwrap();
    let second = advance(
        first.parcels(),
        &paths,
        0.025,
        &endpoints(&base),
        network(),
        || false,
    )
    .unwrap();
    let estimate = compare(
        &base,
        coarse.parcels(),
        second.parcels(),
        0.025,
        0.025,
        controls(),
        &mut || false,
    )
    .unwrap();
    // Independently subtract the explicit species gauge from absolute U for
    // the ordinary-reference case. The production estimator never does this.
    let gauged = |gas: &IsochoricRealGasParcel| {
        let inv = inventory(gas.model().source(), gas.mass_kg()).unwrap();
        let mut value = gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
        for species in gas.model().source().thermodynamics().components() {
            let u = species
                .evaluate_properties(Temperature::new(300.0))
                .unwrap()
                .u()
                .value()
                / species.molar_mass().value();
            value -= inv[species.species()] * u;
        }
        value
    };
    let expected = base
        .keys()
        .map(|id| (gauged(&coarse.parcels()[id]) - gauged(&second.parcels()[id])).abs())
        .fold(0.0, f64::max);
    close(estimate.maximum_energy_difference_j, expected, 1e-9);
}

#[test]
fn gross_limits_subdivide_and_path_order_orientation_preserve_evolution() {
    let original = nodes(0.0);
    let mut paths = links();
    paths[0].solve.maximum_donor_fraction = 5e-5;
    let mut second = paths[0];
    second.id = 2;
    second.nozzle.area_m2 *= 0.5;
    paths.push(second);
    let a = evolve(&original, &paths, 0.1, controls());
    assert!(a.report().feasibility_rejections > 0);
    paths.reverse();
    for path in &mut paths {
        path.parcels.swap(0, 1);
    }
    let b = evolve(&original, &paths, 0.1, controls());
    assert_eq!(a, b);
    assert_eq!(
        a.report().paths.iter().map(|p| p.id).collect::<Vec<_>>(),
        vec![2, 9]
    );
    close(
        a.report().paths[1].mass_left_to_right_kg,
        2.0 * a.report().paths[0].mass_left_to_right_kg,
        1e-15,
    );
}

#[test]
fn work_accuracy_source_and_late_cancellation_refuse_the_whole_window() {
    let original = nodes(0.0);
    let saved = original.clone();
    let paths = links();
    let input = endpoints(&original);
    let trial = |settings| {
        advance_adaptive(&original, &paths, 0.1, &input, network(), settings, || {
            false
        })
    };
    assert!(matches!(
        trial(RealGasNozzleNetworkTimeSolve {
            maximum_trials: 1,
            ..controls()
        }),
        Err(Error::Budget("adaptive nozzle comparisons"))
    ));
    assert!(matches!(
        trial(RealGasNozzleNetworkTimeSolve {
            maximum_steps: 2,
            ..controls()
        }),
        Err(Error::Budget("accepted nozzle half steps"))
    ));
    assert!(matches!(
        trial(RealGasNozzleNetworkTimeSolve {
            minimum_dt_s: 0.1,
            ..controls()
        }),
        Err(Error::Accuracy { .. })
    ));
    assert!(matches!(
        trial(RealGasNozzleNetworkTimeSolve {
            energy_absolute_tolerance_j: f64::NAN,
            ..controls()
        }),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        trial(RealGasNozzleNetworkTimeSolve {
            maximum_trials: 65_537,
            ..controls()
        }),
        Err(Error::Invalid(_))
    ));
    let mut wrong = input.clone();
    wrong.get_mut(&20).unwrap().maximum_temperature_k = 290.0;
    assert!(matches!(
        advance_adaptive(
            &original,
            &paths,
            0.1,
            &wrong,
            network(),
            controls(),
            || false
        ),
        Err(Error::Network(_))
    ));
    let mut polls = 0;
    advance_adaptive(
        &original,
        &paths,
        0.1,
        &input,
        network(),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let stop = polls - 1;
    let mut cancelled_polls = 0;
    assert!(matches!(
        advance_adaptive(
            &original,
            &paths,
            0.1,
            &input,
            network(),
            controls(),
            || {
                cancelled_polls += 1;
                cancelled_polls >= stop
            }
        ),
        Err(Error::Network(RealGasNozzleNetworkError::Cancelled))
    ));
    assert_eq!(original, saved);
}

#[test]
fn closed_paths_and_small_final_windows_keep_exact_states_without_rounding_time_away() {
    let original = nodes(0.0);
    let mut paths = links();
    paths[0].nozzle.area_m2 = 0.0;
    let result = evolve(
        &original,
        &paths,
        0.103,
        RealGasNozzleNetworkTimeSolve {
            minimum_dt_s: 0.01,
            ..controls()
        },
    );
    assert_eq!(result.parcels(), &original);
    assert_eq!(result.report().duration_s, 0.103);
    assert_eq!(result.report().accepted_steps, 4);
    assert_eq!(result.report().maximum_accepted_normalized_error, 0.0);
    assert_eq!(result.report().paths[0].mass_left_to_right_kg, 0.0);
    assert_eq!(result.report().paths[0].enthalpy_left_to_right_j, 0.0);
    let tiny = f64::from_bits(1);
    assert!(matches!(
        advance_adaptive(
            &original,
            &paths,
            tiny,
            &endpoints(&original),
            network(),
            controls(),
            || false
        ),
        Err(Error::Invalid(_))
    ));
}
