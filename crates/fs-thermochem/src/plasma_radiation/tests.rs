//! G0/G2/G3: source normalization, original populations and applicability.

use super::*;
use crate::fusion::{NuclearSpeciesV1, plasma::DtFusionPlasmaV1};
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, InitialStatePolicy,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure};
use std::collections::BTreeMap;

const NAMES: [&str; 4] = ["ion_d", "ion_t", "ash_he", "electron"];
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn bank(electron_mass: f64, label: &str) -> IdealGasEquilibriumV1 {
    let masses = [
        NuclearSpeciesV1::Deuterium.mass_kg(),
        NuclearSpeciesV1::Tritium.mass_kg(),
        NuclearSpeciesV1::Helium4.mass_kg(),
        electron_mass,
    ];
    let models = NAMES
        .into_iter()
        .zip(masses)
        .map(|(name, mass)| {
            let mut parameters = BTreeMap::new();
            for ((name, dims), value) in required_nasa9_parameter_dims()
                .into_iter()
                .take(9)
                .zip([0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            {
                parameters.insert(name.into(), LawParameter { value, dims });
            }
            parameters.insert(
                "reference_pressure".into(),
                LawParameter {
                    value: 1e5,
                    dims: Pressure::DIMS,
                },
            );
            let card = ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 1e3, 2e9),
                sources: vec![],
                provenance: Provenance {
                    source: label.into(),
                    license: "analytic test fixture".into(),
                    artifact: None,
                },
            };
            Nasa9StandardStateModelV1::new(
                id(name),
                MolarMass::new(mass * AVOGADRO_PER_MOL),
                StandardStateConventionV1::new(
                    StandardStatePhaseV1::Gas,
                    ReferenceEquationOfStateV1::IdealGas,
                    Pressure::new(1e5),
                    ElementalReferenceIdV1::new("radiation-sensible-test").unwrap(),
                )
                .unwrap(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .collect();
    let names: Vec<_> = NAMES.map(id).into();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            names.clone(),
            vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(names, vec![1, 1, 2, -1]).unwrap(),
        models,
    )
    .unwrap()
}

fn controls() -> BremsstrahlungControlsV1 {
    BremsstrahlungControlsV1 {
        frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
        minimum_electron_temperature_k: 1e4,
        maximum_electron_temperature_k: 2e9,
        maximum_electron_thermal_rest_energy_ratio: 0.1,
        maximum_ion_to_electron_thermal_variance_ratio: 0.05,
        maximum_electron_thermal_phase_space_density: 0.01,
    }
}
fn fixture() -> ElectronIonBremsstrahlungV1 {
    ElectronIonBremsstrahlungV1::new(
        bank(ELECTRON_MASS_KG, "manufactured monatomic source"),
        controls(),
    )
    .unwrap()
}
fn amounts(law: &ElectronIonBremsstrahlungV1, densities: [f64; 3], volume: f64) -> [Amount; 4] {
    let mut result = [Amount::new(0.0); 4];
    for (name, density) in NAMES.into_iter().take(3).zip(densities) {
        let index = law
            .thermodynamics()
            .elemental_matrix()
            .species()
            .binary_search(&id(name))
            .unwrap();
        result[index] = Amount::new(density * volume / AVOGADRO_PER_MOL);
    }
    result[law.electron_index()] =
        Amount::new((densities[0] + densities[1] + 2.0 * densities[2]) * volume / AVOGADRO_PER_MOL);
    result
}
fn evaluate(
    law: &ElectronIonBremsstrahlungV1,
    densities: [f64; 3],
    volume: f64,
    te: f64,
) -> BremsstrahlungEvaluationV1 {
    law.evaluate_point(
        &amounts(law, densities, volume),
        Volume::new(volume),
        Temperature::new(1e8),
        Temperature::new(te),
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs(),
        "actual={a:.17e}, expected={b:.17e}"
    );
}

#[test]
fn g2_published_cgs_normalization_converts_to_si_once() {
    let law = fixture();
    let result = evaluate(&law, [5e19, 5e19, 0.0], 2.0, 1e8);
    // Itoh et al., equation (11): 1.7112e5 erg/(s cm³) at ne=ni=1e14 cm⁻³,
    // Te=1e8 K and explicitly declared g=1.2. Multiply by0.1 for W/m³.
    close(result.power_density_w_m3(), 17_112.0, 8.0 * f64::EPSILON);
    close(result.total_power_w(), 34_224.0, 8.0 * f64::EPSILON);
    close(
        result.total_power_over_sqrt_temperature_w_k_sqrt(),
        3.4224,
        8.0 * f64::EPSILON,
    );
    close(result.electron_density_per_m3(), 1e20, 4.0 * f64::EPSILON);
    close(
        result.charge_squared_ion_density_per_m3(),
        1e20,
        4.0 * f64::EPSILON,
    );
    close(
        result.derivative_w_per_k(),
        34_224.0 / (2.0e8),
        8.0 * f64::EPSILON,
    );
}

#[test]
fn g3_charge_square_weight_includes_newly_produced_helium_ash() {
    let law = fixture();
    let fuel = evaluate(&law, [5e19, 5e19, 0.0], 1.0, 1e8);
    let ash = evaluate(&law, [4e19, 4e19, 1e19], 1.0, 1e8);
    // D+ + T+ -> He++ preserves electron amount but increases sum(ni*Zi²).
    close(
        ash.electron_density_per_m3(),
        fuel.electron_density_per_m3(),
        4.0 * f64::EPSILON,
    );
    close(
        ash.charge_squared_ion_density_per_m3(),
        1.2e20,
        4.0 * f64::EPSILON,
    );
    close(
        ash.total_power_w() / fuel.total_power_w(),
        1.2,
        8.0 * f64::EPSILON,
    );
}

#[test]
fn g3_volume_density_temperature_and_gaunt_scaling() {
    let law = fixture();
    let p = evaluate(&law, [5e19, 5e19, 0.0], 1.0, 1e8);
    let extensive = evaluate(&law, [5e19, 5e19, 0.0], 4.0, 1e8);
    close(
        extensive.total_power_w(),
        4.0 * p.total_power_w(),
        4.0 * f64::EPSILON,
    );
    let diluted = law
        .evaluate_point(
            &amounts(&law, [5e19, 5e19, 0.0], 1.0),
            Volume::new(2.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || false,
        )
        .unwrap();
    close(
        diluted.total_power_w(),
        0.5 * p.total_power_w(),
        4.0 * f64::EPSILON,
    );
    close(
        diluted.power_density_w_m3(),
        0.25 * p.power_density_w_m3(),
        4.0 * f64::EPSILON,
    );
    let hotter = evaluate(&law, [5e19, 5e19, 0.0], 1.0, 4e8);
    close(
        hotter.total_power_w(),
        2.0 * p.total_power_w(),
        4.0 * f64::EPSILON,
    );
    assert_eq!(
        hotter.total_power_over_sqrt_temperature_w_k_sqrt(),
        p.total_power_over_sqrt_temperature_w_k_sqrt()
    );
    let mut c = controls();
    c.frequency_integrated_gaunt_factor = BremsstrahlungGauntFactorV1::Constant(2.4);
    let alternate = ElectronIonBremsstrahlungV1::new(law.thermodynamics().clone(), c).unwrap();
    close(
        evaluate(&alternate, [5e19, 5e19, 0.0], 1.0, 1e8).total_power_w(),
        2.0 * p.total_power_w(),
        4.0 * f64::EPSILON,
    );
    let h = 1e3;
    let plus = evaluate(&law, [5e19, 5e19, 0.0], 1.0, 1e8 + h).total_power_w();
    let minus = evaluate(&law, [5e19, 5e19, 0.0], 1.0, 1e8 - h).total_power_w();
    close((plus - minus) / (2.0 * h), p.derivative_w_per_k(), 1e-9);
}

#[test]
fn g0_original_bank_binding_and_original_source_refusals() {
    let original = bank(ELECTRON_MASS_KG, "manufactured monatomic source");
    let law = ElectronIonBremsstrahlungV1::new(original.clone(), controls()).unwrap();
    let fusion =
        DtFusionPlasmaV1::new(original.clone(), id(NAMES[0]), id(NAMES[1]), id(NAMES[2])).unwrap();
    assert_eq!(fusion.electron_index(), law.electron_index());
    assert!(core::ptr::eq(
        original.models().as_ptr(),
        law.thermodynamics().models().as_ptr()
    ));
    let ns = amounts(&law, [5e19, 5e19, 0.0], 1.0);
    let state = FrozenIonizedGasStateV1::new(
        &original,
        &ns,
        Volume::new(1.0),
        Temperature::new(1e8),
        Temperature::new(1e8),
        1e-12,
        || false,
    )
    .unwrap();
    let before = state.clone();
    law.evaluate(&state, || false).unwrap();
    assert_eq!(state, before);
    let other = ElectronIonBremsstrahlungV1::new(
        bank(ELECTRON_MASS_KG, "different original provenance"),
        controls(),
    )
    .unwrap();
    assert!(matches!(
        other.evaluate(&state, || false),
        Err(BremsstrahlungErrorV1::Invalid(_))
    ));
    assert!(matches!(
        law.evaluate_point(
            &ns,
            Volume::new(1.0),
            Temperature::new(3e9),
            Temperature::new(1e8),
            || false
        ),
        Err(BremsstrahlungErrorV1::Species { .. })
    ));
    assert!(
        ElectronIonBremsstrahlungV1::new(
            bank(2.0 * ELECTRON_MASS_KG, "wrong electron mass"),
            controls()
        )
        .is_err()
    );
}

#[test]
fn g0_nonrelativistic_nondegenerate_and_stationary_ion_domains_are_explicit() {
    let law = fixture();
    let ns = amounts(&law, [5e19, 5e19, 0.0], 1.0);
    for (th, te) in [(1e8, 1e9), (1e9, 1e5), (1e8, 1e3)] {
        assert!(matches!(
            law.evaluate_point(
                &ns,
                Volume::new(1.0),
                Temperature::new(th),
                Temperature::new(te),
                || false
            ),
            Err(BremsstrahlungErrorV1::Domain(_))
        ));
    }
    assert!(matches!(
        law.evaluate_point(
            &amounts(&law, [5e39, 5e39, 0.0], 1.0),
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || false
        ),
        Err(BremsstrahlungErrorV1::Domain(_))
    ));
    let mut c = controls();
    c.maximum_ion_to_electron_thermal_variance_ratio = 2e-4;
    let strict = ElectronIonBremsstrahlungV1::new(law.thermodynamics().clone(), c).unwrap();
    // Absent lighter fuel does not impose its ion-motion limit on helium alone.
    let helium = evaluate(&strict, [0.0, 0.0, 5e19], 1.0, 1e8);
    assert!(helium.maximum_ion_to_electron_thermal_variance_ratio() < 2e-4);
    assert!(
        strict
            .evaluate_point(
                &ns,
                Volume::new(1.0),
                Temperature::new(1e8),
                Temperature::new(1e8),
                || false
            )
            .is_err()
    );
}

#[test]
fn g0_absent_carriers_are_exact_zero_and_invalid_inputs_are_refused() {
    let law = fixture();
    let mut ns = amounts(&law, [5e19, 5e19, 0.0], 1.0);
    ns[law.electron_index()] = Amount::new(0.0);
    let absent = law
        .evaluate_point(
            &ns,
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || false,
        )
        .unwrap();
    assert_eq!(absent.total_power_w(), 0.0);
    assert_eq!(absent.derivative_w_per_k(), 0.0);
    assert_eq!(evaluate(&law, [0.0; 3], 1.0, 1e8).total_power_w(), 0.0);
    for volume in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            law.evaluate_point(
                &ns,
                Volume::new(volume),
                Temperature::new(1e8),
                Temperature::new(1e8),
                || false
            )
            .is_err()
        );
    }
    assert_eq!(
        law.evaluate_point(
            &ns,
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || true
        ),
        Err(BremsstrahlungErrorV1::Cancelled)
    );
    assert!(
        law.evaluate_point(
            &ns[..3],
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || false
        )
        .is_err()
    );
    ns[0] = Amount::new(-1.0);
    assert!(
        law.evaluate_point(
            &ns,
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(1e8),
            || false
        )
        .is_err()
    );
    for g in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut c = controls();
        c.frequency_integrated_gaunt_factor = BremsstrahlungGauntFactorV1::Constant(g);
        assert!(c.validate().is_err());
    }
    let mut c = controls();
    c.maximum_electron_thermal_phase_space_density = 1.0;
    assert!(c.validate().is_err());
}

#[test]
fn g2_itoh_table2_anchors_and_charge_temperature_scaling() {
    // Independently evaluated with 70-digit Decimal arithmetic, direct powers
    // rather than the implementation's simultaneous Horner operation tree.
    let anchors = [
        (1.0, 5e7, 1.1664502635097814, -6.216533811375891e-10),
        (2.0, 5e7, 1.2249777654893104, -1.0801979788864634e-9),
        (1.0, 1e8, 1.1487442167500536, -2.0746969464285134e-10),
        (2.0, 1e8, 1.1918104241032491, -4.20284717595417e-10),
        (1.0, 1579.0, 1.160573, 2.6939900060587048e-5),
        (1.0, 157900000.0, 1.139223, -1.5166916524559948e-10),
    ];
    for (z, temperature, expected, expected_derivative) in anchors {
        let (g, derivative_log_temperature) = itoh_nonrelativistic_gaunt(z, temperature).unwrap();
        close(g, expected, 2e-14);
        close(
            derivative_log_temperature / temperature,
            expected_derivative,
            2e-12,
        );
    }
    // At x=0 the value and slope select the two leading table coefficients.
    let (g, derivative_log_temperature) =
        itoh_nonrelativistic_gaunt(1.0, ITOH_RYDBERG_TEMPERATURE_K * sqrt(10.0)).unwrap();
    close(g, 1.42752, 4.0 * f64::EPSILON);
    close(derivative_log_temperature, -0.029649631715120525, 2e-14);
    let hydrogen = itoh_nonrelativistic_gaunt(1.0, 3e7).unwrap();
    let helium = itoh_nonrelativistic_gaunt(2.0, 1.2e8).unwrap();
    assert_eq!(hydrogen, helium); // Same Z²/Te, including the logarithmic slope.
}

#[test]
fn g3_itoh_mixed_ion_power_and_full_temperature_derivative() {
    let constant = fixture();
    let mut c = controls();
    c.frequency_integrated_gaunt_factor = BremsstrahlungGauntFactorV1::ItohNonrelativistic;
    let law = ElectronIonBremsstrahlungV1::new(constant.thermodynamics().clone(), c).unwrap();
    let densities = [4e19, 3e19, 1.5e19];
    let point = evaluate(&law, densities, 1.0, 1e8);
    // Independent 70-digit value: C*ne*sqrt(Te)*(7e19*g1 + 6e19*g2).
    close(point.power_density_w_m3(), 21663.894760226435, 2e-14);
    close(point.derivative_w_per_k(), 1.0265255526546084e-4, 2e-14);
    let h = 1e3;
    let plus = evaluate(&law, densities, 1.0, 1e8 + h).total_power_w();
    let minus = evaluate(&law, densities, 1.0, 1e8 - h).total_power_w();
    close((plus - minus) / (2.0 * h), point.derivative_w_per_k(), 1e-8);
    // A fitted coefficient varies with temperature; omitting dg/dTe would
    // yield P/(2Te), which disagrees with the physical source fit's slope.
    assert!((point.derivative_w_per_k() / (point.total_power_w() / (2.0e8)) - 1.0).abs() > 0.04);
    assert!(
        point.total_power_over_sqrt_temperature_w_k_sqrt()
            < evaluate(&law, densities, 1.0, 5e7).total_power_over_sqrt_temperature_w_k_sqrt()
    );
}

#[test]
fn g0_itoh_bounds_are_inclusive_and_apply_to_every_populated_ion() {
    for z in [1.0, 2.0] {
        let lower = z * z * ITOH_RYDBERG_TEMPERATURE_K / 100.0;
        let upper = z * z * ITOH_RYDBERG_TEMPERATURE_K / 1e-3;
        assert!(itoh_nonrelativistic_gaunt(z, lower).is_ok());
        assert!(itoh_nonrelativistic_gaunt(z, upper).is_ok());
        assert!(itoh_nonrelativistic_gaunt(z, lower.next_down()).is_err());
        assert!(itoh_nonrelativistic_gaunt(z, upper.next_up()).is_err());
    }
    let constant = fixture();
    let mut c = controls();
    c.frequency_integrated_gaunt_factor = BremsstrahlungGauntFactorV1::ItohNonrelativistic;
    c.minimum_electron_temperature_k = 1e3;
    let fit = ElectronIonBremsstrahlungV1::new(constant.thermodynamics().clone(), c).unwrap();
    let fuel = amounts(&fit, [5e15, 5e15, 0.0], 1.0);
    let with_ash = amounts(&fit, [4.9e15, 4.9e15, 1e14], 1.0);
    // At4000K hydrogen's gamma² lies in the fit, helium's does not. An absent
    // ash row does not constrain the law until it actually becomes populated.
    fit.evaluate_point(
        &fuel,
        Volume::new(1.0),
        Temperature::new(4000.0),
        Temperature::new(4000.0),
        || false,
    )
    .unwrap();
    assert!(matches!(
        fit.evaluate_point(
            &with_ash,
            Volume::new(1.0),
            Temperature::new(4000.0),
            Temperature::new(4000.0),
            || false
        ),
        Err(BremsstrahlungErrorV1::Domain(_))
    ));
    // Conversely, helium supports200MK but populated hydrogen does not.
    evaluate(&fit, [0.0, 0.0, 5e19], 1.0, 2e8);
    assert!(matches!(
        fit.evaluate_point(
            &amounts(&fit, [5e19, 5e19, 0.0], 1.0),
            Volume::new(1.0),
            Temperature::new(1e8),
            Temperature::new(2e8),
            || false
        ),
        Err(BremsstrahlungErrorV1::Domain(_))
    ));
    // The explicitly selected constant alternative remains available; a fit
    // refusal is never silently converted into that different source model.
    evaluate(&constant, [5e19, 5e19, 0.0], 1.0, 2e8);
}
