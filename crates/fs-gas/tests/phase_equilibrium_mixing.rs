//! Original enthalpy, species and entropy in common-pressure phase mixing.
#![cfg(feature = "material")]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{METHANE, MIXTURE_MASS, model, source};
use fs_gas::material::phase_mixing::{
    MixtureMixingControl as Control, MixtureMixingError as E, MixtureMixingInlet as Inlet,
    MixtureMixingStep as Step, mix_adiabatic_mixture_streams,
};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel as Model,
    flash::{
        ThermochemicalFlashError,
        equilibrium::{PhaseEnthalpyControls, mixing::PhaseMixingError as NativeError},
    },
};
use fs_qty::{MassFlowRate, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

const CH4_MASS: f64 = 0.016_042_46;
const C2H6_MASS: f64 = 0.030_069_04;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected
    );
}
fn inlet(model: &Model, t: f64, flow: f64, p: f64) -> Inlet<'_> {
    Inlet {
        model,
        temperature: Temperature::new(t),
        pressure: Pressure::new(p),
        mass_flow: MassFlowRate::new(flow),
    }
}
fn controls(reference: f64, bracket: [f64; 2]) -> Control {
    Control::new(
        Temperature::new(reference),
        PhaseEnthalpyControls::new(bracket.map(Temperature::new)),
    )
}
fn composition(base: &Model, methane: f64) -> Model {
    let joined = base
        .source()
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![1.0 - methane, methane]).unwrap(),
        )
        .unwrap();
    let old = base.source().thermodynamics().components();
    let indices: Vec<_> = joined
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            old.binary_search_by(|q| q.species().cmp(s.species()))
                .unwrap()
        })
        .collect();
    let bindings = joined
        .thermodynamics()
        .components()
        .iter()
        .zip(&indices)
        .map(|(s, &index)| {
            let critical = &base.mechanics().components()[index];
            PengRobinsonSpeciesBinding {
                species: s.species().clone(),
                critical_temperature: critical.critical_temperature,
                critical_pressure: critical.critical_pressure,
                acentric_factor: critical.acentric_factor,
            }
        })
        .collect();
    let mut interactions = Vec::new();
    for &i in &indices {
        for &j in &indices {
            interactions.push(base.mechanics().binary_interactions()[i * old.len() + j]);
        }
    }
    Model::new(joined, bindings, interactions).unwrap()
}
fn changed_pair(base: &Model, pair: f64) -> Model {
    let bindings = base
        .source()
        .thermodynamics()
        .components()
        .iter()
        .zip(base.mechanics().components())
        .map(|(s, c)| PengRobinsonSpeciesBinding {
            species: s.species().clone(),
            critical_temperature: c.critical_temperature,
            critical_pressure: c.critical_pressure,
            acentric_factor: c.acentric_factor,
        })
        .collect();
    Model::new(base.source().clone(), bindings, vec![0.0, pair, pair, 0.0]).unwrap()
}
fn cryogenic(
    left: &Model,
    right: &Model,
    control: Control,
    cancel: impl FnMut() -> bool,
) -> Result<Step, E> {
    mix_adiabatic_mixture_streams(
        [inlet(left, 205.0, 0.5, 3e6), inlet(right, 235.0, 1.0, 3e6)],
        None,
        control,
        cancel,
    )
}

#[test]
fn g1_distinct_liquid_and_vapor_feeds_form_the_independent_mixed_split_state() {
    let base = source();
    let left = composition(&base, 0.2);
    let right = composition(&base, 0.8);
    let result = cryogenic(&left, &right, controls(220.0, [195.0, 255.0]), || false).unwrap();
    let outlet = result.mixture.outlet.state();
    // Independent direct PR fugacity/Gibbs selection with conserved species
    // molar flows and analytic h/s departures, not the production PH owner.
    close(outlet.temperature_k(), 224.839_063_213_158_85, 3e-7);
    close(outlet.density_kg_m3(), 60.541_089_102_697_214, 3e-6);
    close(
        outlet.vapor_mass_fraction().unwrap(),
        0.650_950_996_948_399_6,
        3e-8,
    );
    close(
        outlet.phase().two_phase().unwrap().vapor_fraction_mol(),
        0.712_067_894_876_149_6,
        3e-8,
    );
    close(result.entropy_production_w_k, 13.569_438_151_808_754, 3e-6);
    close(
        result.mixture.model.mechanics().mole_fractions()[METHANE],
        0.645_878_963_176_998_2,
        2e-15,
    );
    close(
        result.mixture.model.mechanics().molar_mass().value(),
        0.021_009_567_052_680_78,
        1e-16,
    );
    assert!(
        result
            .mixture
            .inlets
            .iter()
            .all(|s| s.phase().two_phase().is_none())
    );
    assert_eq!(result.mass_flow.value(), 1.5);
    assert_eq!(result.species.len(), 2);
    assert!(result.mass_residual_kg_s.abs() < 1e-14);
    assert!(result.species.iter().all(|s| s.residual_kg_s.abs() < 1e-14));
    assert!(result.energy_residual_w.abs() < 3e-5);
    let methane = result
        .species
        .iter()
        .find(|s| s.species == SpeciesId::new("methane").unwrap())
        .unwrap();
    let m_left = 0.2 * CH4_MASS + 0.8 * C2H6_MASS;
    let m_right = 0.8 * CH4_MASS + 0.2 * C2H6_MASS;
    close(
        methane.inlet_mass_flow_kg_s[0],
        0.5 * 0.2 * CH4_MASS / m_left,
        1e-15,
    );
    close(
        methane.inlet_mass_flow_kg_s[1],
        0.8 * CH4_MASS / m_right,
        1e-15,
    );
    close(
        methane.outlet_mass_flow_kg_s,
        methane.inlet_mass_flow_kg_s.iter().sum(),
        1e-15,
    );
    // The result owns the mixed law and can continue through the original PH API.
    let next = result
        .mixture
        .model
        .phase_equilibrium_pressure_enthalpy_increment(
            Temperature::new(outlet.temperature_k()),
            Pressure::new(3e6),
            Pressure::new(3e6),
            1000.0,
            controls(220.0, [195.0, 255.0]).phase.equilibrium,
            || false,
        )
        .unwrap();
    close(next.accepted_increment_j_kg().unwrap(), 1000.0, 1e-5);
    assert!(next.state().temperature_k() > outlet.temperature_k());
}

#[test]
fn g1_real_heat_of_mixing_changes_temperature_even_when_both_feeds_have_equal_temperature() {
    let base = source();
    let left = composition(&base, 0.2);
    let right = composition(&base, 0.8);
    let result = mix_adiabatic_mixture_streams(
        [
            inlet(&left, 350.0, 0.5, 3e6),
            inlet(&right, 350.0, 1.0, 3e6),
        ],
        None,
        controls(350.0, [330.0, 370.0]),
        || false,
    )
    .unwrap();
    close(
        result.mixture.outlet.state().temperature_k(),
        349.174_247_846_999_44,
        3e-7,
    );
    close(
        result.mixture.outlet.state().density_kg_m3(),
        23.322_307_930_725_27,
        1e-6,
    );
    close(result.entropy_production_w_k, 87.008_588_948_907_35, 3e-6);
    assert!(result.mixture.outlet.state().phase().two_phase().is_none());
    assert!(result.mixture.outlet.state().temperature_k() < 349.2);
    assert!(result.energy_residual_w.abs() < 3e-5);
}

#[test]
fn g1_dilute_same_composition_recovers_exact_thermal_mixing_and_entropy() {
    let base = source();
    let result = mix_adiabatic_mixture_streams(
        [inlet(&base, 300.0, 0.5, 1.0), inlet(&base, 400.0, 1.0, 1.0)],
        None,
        controls(350.0, [295.0, 405.0]),
        || false,
    )
    .unwrap();
    let cp = 4.0 * R / MIXTURE_MASS;
    let ideal_t: f64 = (0.5 * 300.0 + 400.0) / 1.5;
    let ideal_entropy = cp * (0.5 * (ideal_t / 300.0).ln() + (ideal_t / 400.0).ln());
    close(
        result.mixture.outlet.state().temperature_k(),
        366.666_666_411_371_47,
        3e-7,
    );
    close(result.mixture.outlet.state().temperature_k(), ideal_t, 1e-6);
    close(result.entropy_production_w_k, 19.219_788_764_808_982, 1e-6);
    close(result.entropy_production_w_k, ideal_entropy, 1e-5);
    assert_eq!(result.ideal_mixing_entropy_flow_w_k, 0.0);
}

#[test]
fn g1_pure_feeds_mix_with_explicit_original_cross_pairs_and_keep_ideal_composition_entropy() {
    let base = source();
    let methane = composition(&base, 1.0);
    let ethane = composition(&base, 0.0);
    let feeds = [
        inlet(&methane, 350.0, 0.5, 1.0),
        inlet(&ethane, 350.0, 1.0, 1.0),
    ];
    let control = controls(350.0, [330.0, 370.0]);
    assert!(matches!(
        mix_adiabatic_mixture_streams(feeds, None, control, || false),
        Err(E::Material(NativeError::Catalog(_)))
    ));
    let result = mix_adiabatic_mixture_streams(feeds, Some(&base), control, || false).unwrap();
    let n = [0.5 / CH4_MASS, 1.0 / C2H6_MASS];
    let total = n[0] + n[1];
    let z = [n[0] / total, n[1] / total];
    let entropy = -R * (n[0] * z[0].ln() + n[1] * z[1].ln());
    close(
        result.mixture.model.mechanics().mole_fractions()[METHANE],
        z[0],
        1e-15,
    );
    close(result.ideal_mixing_entropy_flow_w_k, entropy, 1e-12);
    close(result.entropy_production_w_k, entropy, 1e-4);
    close(result.mixture.outlet.state().temperature_k(), 350.0, 1e-5);
    assert_eq!(
        result
            .mixture
            .model
            .source()
            .thermodynamics()
            .components()
            .len(),
        2
    );
    assert_eq!(result.mixture.model.source().transport_materials().len(), 2);
}

#[test]
fn g3_species_reference_shifts_reference_temperature_and_feed_order_preserve_physics() {
    let base = source();
    let left = composition(&base, 0.2);
    let right = composition(&base, 0.8);
    let control = controls(220.0, [195.0, 255.0]);
    let plain = cryogenic(&left, &right, control, || false).unwrap();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let shifted_left = composition(&shifted, 0.2);
    let shifted_right = composition(&shifted, 0.8);
    let gauged = cryogenic(&shifted_left, &shifted_right, control, || false).unwrap();
    assert_eq!(
        plain.mixture.outlet.state().temperature_k(),
        gauged.mixture.outlet.state().temperature_k()
    );
    assert_eq!(
        plain.mixture.outlet.state().phase(),
        gauged.mixture.outlet.state().phase()
    );
    assert_eq!(
        plain.mixture.target_enthalpy_increment_j_kg,
        gauged.mixture.target_enthalpy_increment_j_kg
    );
    assert_eq!(plain.energy_residual_w, gauged.energy_residual_w);
    assert_eq!(plain.entropy_production_w_k, gauged.entropy_production_w_k);
    assert_eq!(plain.species, gauged.species);
    for reference in [225.0, 230.0] {
        let moved =
            cryogenic(&left, &right, controls(reference, [195.0, 255.0]), || false).unwrap();
        close(
            moved.mixture.outlet.state().temperature_k(),
            plain.mixture.outlet.state().temperature_k(),
            3e-7,
        );
        close(
            moved.entropy_production_w_k,
            plain.entropy_production_w_k,
            3e-6,
        );
        assert_eq!(moved.species, plain.species);
    }
    let reverse = mix_adiabatic_mixture_streams(
        [
            inlet(&right, 235.0, 1.0, 3e6),
            inlet(&left, 205.0, 0.5, 3e6),
        ],
        None,
        control,
        || false,
    )
    .unwrap();
    close(
        reverse.mixture.outlet.state().temperature_k(),
        plain.mixture.outlet.state().temperature_k(),
        1e-10,
    );
    close(
        reverse.entropy_production_w_k,
        plain.entropy_production_w_k,
        1e-10,
    );
    assert_eq!(reverse.mixture.pt_evaluations, plain.mixture.pt_evaluations);
}

#[test]
fn g0_identical_feeds_reuse_the_mixed_reference_with_no_caloric_iteration() {
    let base = source();
    let mut control = controls(350.0, [330.0, 370.0]);
    control.phase.maximum_pt_evaluations = 3;
    let result =
        mix_adiabatic_mixture_streams([inlet(&base, 350.0, 1.0, 3e6); 2], None, control, || false)
            .unwrap();
    assert_eq!(result.mixture.outlet.state().temperature_k(), 350.0);
    assert_eq!(result.mixture.target_enthalpy_increment_j_kg, 0.0);
    assert_eq!(result.energy_residual_w, 0.0);
    assert_eq!(result.entropy_production_w_k, 0.0);
    assert_eq!(result.ideal_mixing_entropy_flow_w_k, 0.0);
    assert_eq!(result.mixture.pt_evaluations, 3);
    assert_eq!(result.mixture.outlet.pt_evaluations(), 0);
}

#[test]
fn g1_quadratic_mixing_entropy_survives_nearly_identical_compositions() {
    let base = source();
    let delta = 1e-8;
    let left = composition(&base, 0.5 - delta);
    let right = composition(&base, 0.5 + delta);
    // Equal molar streams. Their leading ideal mixing entropy is 4 R delta².
    // Separate temperatures and a reference away from their mean keep the PH
    // target resolved while this composition contribution remains quadratic.
    let result = mix_adiabatic_mixture_streams(
        [
            inlet(&left, 300.0, left.mechanics().molar_mass().value(), 1.0),
            inlet(&right, 400.0, right.mechanics().molar_mass().value(), 1.0),
        ],
        None,
        controls(340.0, [295.0, 405.0]),
        || false,
    )
    .unwrap();
    close(
        result.ideal_mixing_entropy_flow_w_k,
        4.0 * R * delta * delta,
        3e-22,
    );
    assert!(result.ideal_mixing_entropy_flow_w_k > 0.0);
}

#[test]
fn g0_incompatible_sources_pairs_and_pressures_are_native_refusals() {
    let base = source();
    let left = composition(&base, 0.2);
    let right = composition(&base, 0.8);
    let conflict = changed_pair(&right, 0.1);
    let control = controls(220.0, [195.0, 255.0]);
    assert!(matches!(
        cryogenic(&left, &conflict, control, || false),
        Err(E::Material(NativeError::Catalog(_)))
    ));
    let wrong_datum = composition(&model([1e6, 0.0], [0.0; 2], [100.0, 450.0]), 0.8);
    assert!(matches!(
        cryogenic(&left, &wrong_datum, control, || false),
        Err(E::Material(NativeError::Catalog(_)))
    ));
    let wrong_entropy = composition(&model([0.0; 2], [0.0, 100.0], [100.0, 450.0]), 0.8);
    assert!(matches!(
        cryogenic(&left, &wrong_entropy, control, || false),
        Err(E::Material(NativeError::Catalog(_)))
    ));
    assert!(matches!(
        mix_adiabatic_mixture_streams(
            [
                inlet(&left, 205.0, 0.5, 3e6),
                inlet(&right, 235.0, 1.0, 2e6)
            ],
            None,
            control,
            || false
        ),
        Err(E::Invalid(_))
    ));
}

#[test]
fn g0_original_domain_whole_call_work_and_final_cancellation_bound_every_attempt() {
    let base = source();
    let left = composition(&base, 0.2);
    let right = composition(&base, 0.8);
    let mut control = controls(220.0, [195.0, 255.0]);
    let mut polls = 0;
    let accepted = cryogenic(&left, &right, control, || {
        polls += 1;
        false
    })
    .unwrap();
    control.phase.maximum_pt_evaluations = accepted.mixture.pt_evaluations;
    let exact = cryogenic(&left, &right, control, || false).unwrap();
    assert_eq!(
        exact.mixture.outlet.state().temperature_k(),
        accepted.mixture.outlet.state().temperature_k()
    );
    assert_eq!(
        exact.mixture.pt_evaluations,
        accepted.mixture.pt_evaluations
    );
    control.phase.maximum_pt_evaluations -= 1;
    match cryogenic(&left, &right, control, || false).unwrap_err() {
        E::Material(NativeError::Phase {
            completed_pt_evaluations,
            source: ThermochemicalFlashError::WorkLimit { flash_evaluations },
            ..
        }) => assert_eq!(
            completed_pt_evaluations + flash_evaluations,
            control.phase.maximum_pt_evaluations
        ),
        E::Material(NativeError::WorkLimit { pt_evaluations }) => {
            assert_eq!(pt_evaluations, control.phase.maximum_pt_evaluations)
        }
        other => panic!("unexpected source work refusal: {other:?}"),
    }
    assert!(matches!(
        cryogenic(&left, &right, controls(220.0, [219.0, 221.0]), || false),
        Err(E::Material(NativeError::Phase {
            source: ThermochemicalFlashError::OutsideEnthalpyBracket,
            ..
        }))
    ));
    let mut replay = 0;
    assert_eq!(
        cryogenic(&left, &right, controls(220.0, [195.0, 255.0]), || {
            replay += 1;
            replay == polls
        })
        .unwrap_err(),
        E::Cancelled
    );
    assert_eq!(
        cryogenic(&left, &right, controls(220.0, [195.0, 255.0]), || false).unwrap(),
        accepted
    );
}
