use super::*;
use crate::eos::peng_robinson::{
    PengRobinsonComponentV1, saturation::PengRobinsonSaturationOptionsV1,
};
use crate::{
    ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, SpeciesId, StandardStateConventionV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use std::collections::BTreeMap;

fn component() -> PengRobinsonComponentV1 {
    PengRobinsonComponentV1 {
        critical_temperature: Temperature::new(190.564),
        critical_pressure: Pressure::new(4.5992e6),
        acentric_factor: 0.01142,
        molar_mass: MolarMass::new(0.016_042_46),
    }
}
fn model() -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(vec![component()], vec![1.0], vec![]).unwrap()
}
fn source(formation: f64, entropy: f64) -> Nasa9StandardStateModelV1 {
    let mut a = [0.0; 9];
    a[2] = 4.0;
    a[7] = formation / R;
    a[8] = entropy / R;
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
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
        validity: ValidityDomain::unconstrained().with("T", 100.0, 500.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "constant-Cp pure-fluid caloric oracle".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new("CH4").unwrap(),
        component().molar_mass,
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1e5),
            ElementalReferenceIdV1::new("pure-fluid-test").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn saturation(model: &PengRobinsonMixtureV1) -> PengRobinsonSaturationV1 {
    model
        .pure_saturation_at_temperature(
            Temperature::new(150.0),
            [Pressure::new(1e5), Pressure::new(4.4e6)],
            PengRobinsonSaturationOptionsV1 {
                log_fugacity_tolerance: 1e-12,
                relative_pressure_tolerance: 1e-12,
                ..Default::default()
            },
        )
        .unwrap()
}
fn path<'a>(
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
) -> PengRobinsonIsobaricCaloricV1<'a> {
    PengRobinsonIsobaricCaloricV1::new(
        model,
        source,
        saturation(model),
        [Temperature::new(120.0), Temperature::new(180.0)],
        PureFluidIsobaricSolveV1::default(),
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute,
        "{a:.17e} != {b:.17e} (allowance {absolute})"
    );
}
fn absolute_h(model: &PengRobinsonMixtureV1, p: Pressure, t: f64, vapor: bool) -> f64 {
    let roots = model.evaluate_pt(Temperature::new(t), p).unwrap();
    (4.0 * R * t
        + roots
            .select(if vapor {
                PengRobinsonPhaseV1::Vapor
            } else {
                PengRobinsonPhaseV1::Liquid
            })
            .residual_enthalpy_j_mol())
        / component().molar_mass.value()
}

#[test]
fn g1_original_latent_heat_selects_quality_and_expansion_at_constant_temperature() {
    let model = model();
    let source = source(0.0, 0.0);
    let path = path(&model, &source);
    let mass = component().molar_mass.value();
    let latent = 6621.875_622_767_658 / mass; // Independent Maxwell-area coexistence oracle.
    close(path.latent_enthalpy_j_kg(), latent, 3e-6);
    let vl = path.saturation().liquid().molar_volume_m3_per_mol();
    let vv = path.saturation().vapor().molar_volume_m3_per_mol();
    for quality in [0.0, 0.1, 0.37, 0.8, 1.0] {
        let h = quality * path.latent_enthalpy_j_kg();
        let result = path.sample_specific_enthalpy(h, || false).unwrap();
        assert_eq!(result.evaluations, 0);
        assert_eq!(result.state.temperature, Temperature::new(150.0));
        assert_eq!(result.state.temperature_slope_kg_k_j, 0.0);
        close(result.state.vapor_fraction(), quality, 2e-16);
        close(
            mass / result.state.density_kg_m3,
            (1.0 - quality) * vl + quality * vv,
            4e-19,
        );
        close(
            result.state.specific_volume_slope_m3_j,
            (vv - vl) / path.saturation().latent_enthalpy_j_mol(),
            1e-22,
        );
        assert!(result.residual_j_kg.abs() <= path.controls().enthalpy_tolerance_j_kg);
    }
    assert!(matches!(
        path.sample_temperature(Temperature::new(150.0), || false),
        Err(E::CoexistenceAmountRequired)
    ));
    let change = path
        .change_between_specific_enthalpies(
            0.2 * path.latent_enthalpy_j_kg(),
            0.8 * path.latent_enthalpy_j_kg(),
            || false,
        )
        .unwrap();
    close(
        change.caloric_entropy_change_j_kg_k,
        0.6 * path.latent_enthalpy_j_kg() / 150.0,
        2e-12,
    );
    close(
        change.source_enthalpy_change_j_kg - 150.0 * change.source_entropy_change_j_kg_k,
        0.6 * path.latent_closure_residual_j_kg(),
        2e-10,
    );
}

#[test]
fn g1_sensible_inverse_and_tangents_use_the_original_eos_and_nasa_sources() {
    let model = model();
    let source = source(0.0, 0.0);
    let path = path(&model, &source);
    let h0 = absolute_h(&model, path.pressure(), 150.0, false);
    for t in [120.0, 133.7, 145.0, 155.0, 168.3, 180.0] {
        let vapor = t > 150.0;
        let point = path
            .sample_temperature(Temperature::new(t), || false)
            .unwrap();
        close(
            point.specific_enthalpy_j_kg,
            absolute_h(&model, path.pressure(), t, vapor) - h0,
            2e-7,
        );
        let inverse = path
            .sample_specific_enthalpy(point.specific_enthalpy_j_kg, || false)
            .unwrap();
        close(inverse.state.temperature.value(), t, 2e-9);
        let step = 0.001;
        let cp = (absolute_h(&model, path.pressure(), t + step, vapor)
            - absolute_h(&model, path.pressure(), t - step, vapor))
            / (2.0 * step);
        close(point.temperature_slope_kg_k_j, 1.0 / cp, 2e-12);
        let volume = |t| {
            let r = model
                .evaluate_pt(Temperature::new(t), path.pressure())
                .unwrap();
            r.select(if vapor {
                PengRobinsonPhaseV1::Vapor
            } else {
                PengRobinsonPhaseV1::Liquid
            })
            .molar_volume_m3_per_mol()
                / component().molar_mass.value()
        };
        close(
            point.specific_volume_slope_m3_j,
            (volume(t + step) - volume(t - step)) / (2.0 * step * cp),
            5e-15,
        );
    }
}

#[test]
fn g3_source_datums_cancel_before_sensible_and_latent_storage() {
    let model = model();
    let ordinary = source(0.0, 0.0);
    let shifted = source(1e14, -1e12);
    let base = path(&model, &ordinary);
    let offset = path(&model, &shifted);
    assert_eq!(
        base.specific_enthalpy_domain_j_kg(),
        offset.specific_enthalpy_domain_j_kg()
    );
    for t in [120.0, 145.0, 155.0, 180.0] {
        assert_eq!(
            base.sample_temperature(Temperature::new(t), || false)
                .unwrap(),
            offset
                .sample_temperature(Temperature::new(t), || false)
                .unwrap()
        );
    }
    let h = base
        .sample_temperature(Temperature::new(169.0), || false)
        .unwrap()
        .specific_enthalpy_j_kg;
    let tiny = base
        .change_between_specific_enthalpies(h, h + 0.001, || false)
        .unwrap();
    let other = offset
        .change_between_specific_enthalpies(h, h + 0.001, || false)
        .unwrap();
    assert_eq!(tiny, other);
    close(tiny.source_enthalpy_change_j_kg, 0.001, 2e-6);
    close(tiny.source_entropy_change_j_kg_k, 0.001 / 169.0, 2e-8);
}

#[test]
fn g1_cross_phase_energy_and_entropy_keep_latent_work_once_and_reverse() {
    let model = model();
    let source = source(0.0, 0.0);
    let path = path(&model, &source);
    let a = path
        .sample_temperature(Temperature::new(130.0), || false)
        .unwrap()
        .specific_enthalpy_j_kg;
    let b = path
        .sample_temperature(Temperature::new(170.0), || false)
        .unwrap()
        .specific_enthalpy_j_kg;
    let change = path
        .change_between_specific_enthalpies(a, b, || false)
        .unwrap();
    let reverse = path
        .change_between_specific_enthalpies(b, a, || false)
        .unwrap();
    close(change.source_enthalpy_change_j_kg, b - a, 2e-6);
    let residual_s = |t: f64, vapor| {
        let roots = model
            .evaluate_pt(Temperature::new(t), path.pressure())
            .unwrap();
        roots
            .select(if vapor {
                PengRobinsonPhaseV1::Vapor
            } else {
                PengRobinsonPhaseV1::Liquid
            })
            .residual_entropy_j_mol_k()
    };
    let expected_s = (4.0 * R * (170.0_f64 / 130.0).ln() + residual_s(170.0, true)
        - residual_s(130.0, false))
        / component().molar_mass.value();
    close(change.source_entropy_change_j_kg_k, expected_s, 2e-8);
    close(
        change.source_enthalpy_change_j_kg,
        -reverse.source_enthalpy_change_j_kg,
        0.0,
    );
    close(
        change.source_entropy_change_j_kg_k,
        -reverse.source_entropy_change_j_kg_k,
        0.0,
    );
    assert!(change.caloric_entropy_change_j_kg_k >= (b - a) / 170.0);
    assert!(change.caloric_entropy_change_j_kg_k <= (b - a) / 130.0);
}

#[test]
fn g4_domains_foreign_sources_and_bounded_work_do_not_change_the_path() {
    let model = model();
    let source = source(0.0, 0.0);
    let path = path(&model, &source);
    let [lo, hi] = path.specific_enthalpy_domain_j_kg();
    assert!(
        path.sample_specific_enthalpy(f64::from_bits(1), || false)
            .is_err()
    );
    for h in [lo.next_down(), hi.next_up(), f64::NAN, f64::INFINITY] {
        assert!(matches!(
            path.sample_specific_enthalpy(h, || false),
            Err(E::OutsideDomain)
        ));
    }
    assert!(
        PengRobinsonIsobaricCaloricV1::new(
            &model,
            &source,
            saturation(&model),
            [Temperature::new(120.0), Temperature::new(191.0)],
            Default::default(),
            || false
        )
        .is_err()
    );
    let mut foreign = component();
    foreign.molar_mass = MolarMass::new(0.02);
    let foreign = PengRobinsonMixtureV1::new(vec![foreign], vec![1.0], vec![]).unwrap();
    assert!(
        PengRobinsonIsobaricCaloricV1::new(
            &foreign,
            &source,
            saturation(&model),
            path.temperature_domain(),
            Default::default(),
            || false
        )
        .is_err()
    );
    let short = PengRobinsonIsobaricCaloricV1::new(
        &model,
        &source,
        saturation(&model),
        path.temperature_domain(),
        PureFluidIsobaricSolveV1 {
            maximum_evaluations: 1,
            ..Default::default()
        },
        || false,
    )
    .unwrap();
    let h = path
        .sample_temperature(Temperature::new(133.7), || false)
        .unwrap()
        .specific_enthalpy_j_kg;
    assert!(matches!(
        short.sample_specific_enthalpy(h, || false),
        Err(E::NotConverged { evaluations: 1, .. })
    ));
    let expected = path.sample_specific_enthalpy(h, || false).unwrap();
    let mut polls = 0;
    path.sample_specific_enthalpy(h, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut attempted = 0;
    assert!(matches!(
        path.sample_specific_enthalpy(h, || {
            attempted += 1;
            attempted == polls
        }),
        Err(E::Cancelled)
    ));
    assert_eq!(attempted, polls);
    assert_eq!(
        path.sample_specific_enthalpy(h, || false).unwrap(),
        expected
    );
}
