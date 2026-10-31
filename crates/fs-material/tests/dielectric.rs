//! G0/G1 passive Debye laws from synthetic source data, not measured validation.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{ELECTRIC_PERMITTIVITY_DIMS, dielectric::*},
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    Dims, QuantitySpec, Time,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

fn card(values: &[(&str, Dims, f64)], validity: ValidityDomain) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Debye law fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn properties() -> Vec<DebyePoleProperties> {
    vec![
        DebyePoleProperties {
            permittivity_increment: "delta1".into(),
            relaxation_time: "tau1".into(),
        },
        DebyePoleProperties {
            permittivity_increment: "delta2".into(),
            relaxation_time: "tau2".into(),
        },
    ]
}
fn resolve(source: &MaterialCard, poles: &[DebyePoleProperties]) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::resolve(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        poles,
        0.02,
        0.005,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn source_geometry_poles_and_energy_keep_their_physical_roles() {
    let source = card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
            ("tau1", Time::DIMS, 0.01),
            ("delta2", ELECTRIC_PERMITTIVITY_DIMS, 5e-6),
            ("tau2", Time::DIMS, 0.2),
        ],
        ValidityDomain::unconstrained(),
    );
    let names = properties();
    let law = resolve(&source, &names);
    assert_eq!(law.high_frequency_property(), "epsilon_inf");
    assert_eq!(law.pole_properties(), names);
    assert_eq!(law.material().properties().len(), 5);
    near(law.instantaneous_capacitance_f(), 8e-6);
    near(law.equilibrium_capacitance_f(), 40e-6);
    near(law.poles()[0].capacitance_f, 12e-6);
    near(law.poles()[1].capacitance_f, 20e-6);
    for pole in law.poles() {
        near(
            pole.resistance_ohm * pole.capacitance_f,
            pole.relaxation_time_s,
        );
    }
    let same = ResolvedDebyeCapacitor::try_from_material(
        law.material().clone(),
        "epsilon_inf",
        &names,
        0.02,
        0.005,
    )
    .unwrap();
    assert_eq!(same, law);
    let volts = [3.0, -2.0, 5.0];
    let charges = [12e-6 * volts[1], 20e-6 * volts[2]];
    near(
        law.stored_energy_j(8e-6 * volts[0], &charges).unwrap(),
        0.5 * (8e-6 * 9.0 + 12e-6 * 4.0 + 20e-6 * 25.0),
    );
    assert!(law.stored_energy_j(0.0, &[]).is_err());
}

#[test]
fn dc_half_pole_and_dispersive_limits_match_analytic_admittance() {
    let source = card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
            ("tau1", Time::DIMS, 0.01),
        ],
        ValidityDomain::unconstrained(),
    );
    let law = resolve(&source, &properties()[..1]);
    assert_eq!(
        law.admittance(0.0).unwrap(),
        DebyeAdmittance {
            conductance_s: 0.0,
            susceptance_s: 0.0
        }
    );
    let at_pole = law.admittance(100.0).unwrap();
    near(at_pole.conductance_s, 0.0006);
    near(at_pole.susceptance_s, 0.0014);
    let low = law.admittance(1e-6).unwrap();
    near(low.susceptance_s / 1e-6, law.equilibrium_capacitance_f());
    let high = law.admittance(1e10).unwrap();
    near(high.susceptance_s / 1e10, law.instantaneous_capacitance_f());
    near(high.conductance_s, 1.0 / law.poles()[0].resistance_ohm);
    for exponent in -9..=9 {
        let response = law.admittance(10.0_f64.powi(exponent)).unwrap();
        assert!(response.conductance_s >= 0.0 && response.susceptance_s > 0.0);
    }
    for omega in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(law.admittance(omega).is_err());
    }
}

#[test]
fn multiple_poles_add_and_zero_poles_preserve_the_lossless_limit() {
    let source = card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
            ("tau1", Time::DIMS, 0.01),
            ("delta2", ELECTRIC_PERMITTIVITY_DIMS, 5e-6),
            ("tau2", Time::DIMS, 0.2),
        ],
        ValidityDomain::unconstrained(),
    );
    let names = properties();
    let all = resolve(&source, &names);
    let first = resolve(&source, &names[..1]);
    let second = resolve(&source, &names[1..]);
    let empty = resolve(&source, &[]);
    near(
        empty.equilibrium_capacitance_f(),
        empty.instantaneous_capacitance_f(),
    );
    for omega in [0.0, 0.1, 10.0, 1000.0] {
        let total = all.admittance(omega).unwrap();
        let a = first.admittance(omega).unwrap();
        let b = second.admittance(omega).unwrap();
        let baseline = empty.admittance(omega).unwrap();
        near(total.conductance_s, a.conductance_s + b.conductance_s);
        near(
            total.susceptance_s,
            a.susceptance_s + b.susceptance_s - baseline.susceptance_s,
        );
        assert_eq!(baseline.conductance_s, 0.0);
    }
    near(
        empty.stored_energy_j(1e-5, &[]).unwrap(),
        0.5 * 1e-10 / 8e-6,
    );
}

#[test]
fn high_frequency_product_overflow_does_not_erase_physical_loss() {
    let source = card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 1e-200),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 1.0),
            ("tau1", Time::DIMS, 1e200),
        ],
        ValidityDomain::unconstrained(),
    );
    let law = resolve(&source, &properties()[..1]);
    let response = law.admittance(1e200).unwrap();
    near(response.conductance_s, 4e-200);
    near(response.susceptance_s, 4.0);
}

#[test]
fn invalid_units_passivity_duplicate_sources_and_frozen_frequency_refuse() {
    let names = properties();
    let attempt = |dims, high, delta, time| {
        let source = card(
            &[
                ("epsilon_inf", dims, high),
                ("delta1", ELECTRIC_PERMITTIVITY_DIMS, delta),
                ("tau1", Time::DIMS, time),
            ],
            ValidityDomain::unconstrained(),
        );
        ResolvedDebyeCapacitor::resolve(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &names[..1],
            1.0,
            1.0,
        )
    };
    assert!(attempt(Dims([0; 6]), 2.0, 1.0, 1.0).is_err());
    for (high, delta, time) in [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 0.0),
        (1.0, 1.0, -1.0),
    ] {
        assert!(attempt(ELECTRIC_PERMITTIVITY_DIMS, high, delta, time).is_err());
    }
    let law = attempt(ELECTRIC_PERMITTIVITY_DIMS, 2.0, 1.0, 1.0).unwrap();
    let duplicate = vec![names[0].clone(), names[0].clone()];
    assert!(
        ResolvedDebyeCapacitor::try_from_material(
            law.material().clone(),
            "epsilon_inf",
            &duplicate,
            1.0,
            1.0
        )
        .is_err()
    );
    assert!(
        ResolvedDebyeCapacitor::try_from_material(
            law.material().clone(),
            "epsilon_inf",
            &names[..1],
            1.0,
            0.0
        )
        .is_err()
    );
    let frequency = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Cyclic),
        ValueForm::Static,
    ));
    let source = card(
        &[("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2.0)],
        ValidityDomain::unconstrained().with_quantity("frequency", frequency, 1.0, 100.0),
    );
    let query = QueryPoint::new()
        .with_quantity("frequency", frequency, 10.0)
        .unwrap();
    assert!(
        ResolvedDebyeCapacitor::resolve(
            &source,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &[],
            1.0,
            1.0
        )
        .is_err()
    );
}
