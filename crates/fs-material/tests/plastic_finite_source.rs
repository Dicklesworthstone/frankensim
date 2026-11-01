//! Synthetic source-card transport and finite-strain mechanics, not material calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::source::{FINITE_J2_LAW_ID, ResolvedFiniteJ2};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Material, FiniteJ2Parameters};
use fs_qty::semantic::{
    FrequencyConvention, QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm,
};
use fs_qty::{Dims, QuantitySpec};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MU: f64 = 80.0e6;
const KAPPA: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const H: f64 = 8.0e6;

fn fixture() -> ConstitutiveModelCard {
    let source = fs_blake3::hash_bytes(
        b"synthetic finite J2 Hencky/Kirchhoff software fixture; no experimental calibration",
    );
    ConstitutiveModelCard {
        law: LawId(FINITE_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", MU, PA),
            ("bulk_modulus_pa", KAPPA, PA),
            ("yield_stress_pa", YIELD, PA),
            ("hardening_modulus_pa", H, PA),
            ("minimum_jacobian", 0.25, Dims::NONE),
            ("maximum_elastic_log_strain", 0.5, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained(),
        sources: vec![source],
        provenance: Provenance {
            source: "authored synthetic multiplicative plasticity fixture".into(),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    }
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn diagonal_deformation(volume_log: f64, axial_deviator: f64) -> [f64; 9] {
    [
        (volume_log + axial_deviator).exp(),
        0.0,
        0.0,
        0.0,
        (volume_log - 0.5 * axial_deviator).exp(),
        0.0,
        0.0,
        0.0,
        (volume_log - 0.5 * axial_deviator).exp(),
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-10 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn exact_source_coefficients_drive_kirchhoff_return_at_nonunit_jacobian() {
    let card = fixture();
    let source = ResolvedFiniteJ2::try_from_model_card(&card).unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.model_card().content_hash(), card.content_hash());
    assert_eq!(source.context(), &QueryPoint::new());
    assert_eq!(
        source.law().parameters(),
        &FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: KAPPA,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: H,
        }
    );
    assert_eq!(
        source.law().limits(),
        &FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        }
    );
    assert!(std::ptr::eq(source.finite_j2(), source.law()));

    // For commuting axisymmetric stretch, the exact log-strain radial return
    // reduces independently to this one-dimensional scalar equation.
    let (v, a) = (0.04, 0.20);
    let f = diagonal_deformation(v, a);
    let old = source.law().initial_state();
    assert_eq!(old.alpha, 0.0);
    assert_eq!(
        old.inverse_plastic_cauchy_green,
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
    let response = source.finite_j2().response(&f, &old).unwrap();
    let increment = (3.0 * MU * a - YIELD) / (3.0 * MU + H);
    let returned = a - increment;
    let pressure = 3.0 * KAPPA * v;
    close(response.jacobian, (3.0 * v).exp());
    close(response.plastic_increment, increment);
    close(response.state.alpha, increment);
    close(response.equivalent_kirchhoff_pa, YIELD + H * increment);
    assert!(response.jacobian > 1.1);
    // Cauchy yield is lower by J; it must not be silently identified with the
    // explicit Kirchhoff source value at this finite volume change.
    assert!(response.equivalent_kirchhoff_pa / response.jacobian < YIELD + H * increment);
    for (index, tau) in [
        (0, pressure + 2.0 * MU * returned),
        (4, pressure - MU * returned),
        (8, pressure - MU * returned),
    ] {
        close(response.kirchhoff_pa[index], tau);
        close(response.first_piola_pa[index], tau / f[index]);
    }
    close(
        response.state.inverse_plastic_cauchy_green[0],
        (-2.0 * increment).exp(),
    );
    close(
        response.state.inverse_plastic_cauchy_green[1],
        increment.exp(),
    );
    close(
        response.state.inverse_plastic_cauchy_green[2],
        increment.exp(),
    );
    close(response.plastic_dissipation_j_m3, YIELD * increment);
    let stored = 1.5 * MU * returned * returned
        + 0.5 * KAPPA * (3.0 * v).powi(2)
        + 0.5 * H * increment * increment;
    close(response.stored_energy_j_m3, stored);
    close(
        response.incremental_potential_j_m3,
        stored + YIELD * increment,
    );
    assert_eq!(old, source.law().initial_state());
}

#[test]
fn every_native_domain_limit_reaches_the_actual_response() {
    let plastic_f = diagonal_deformation(0.0, 0.20);
    for (parameter, limit, f) in [
        ("minimum_jacobian", 0.9, diagonal_deformation(-0.1, 0.0)),
        (
            "maximum_elastic_log_strain",
            0.01,
            diagonal_deformation(0.1, 0.0),
        ),
        ("maximum_plastic_log_strain", 0.01, plastic_f),
        ("maximum_accumulated_plastic_strain", 0.01, plastic_f),
    ] {
        let mut card = fixture();
        card.parameters.get_mut(parameter).unwrap().value = limit;
        let source = ResolvedFiniteJ2::try_from_model_card(&card).unwrap();
        let old = source.law().initial_state();
        assert!(source.law().response(&f, &old).is_err(), "{parameter}");
        assert_eq!(old, source.law().initial_state());
        assert_eq!(source.model_card(), &card);
    }
}

#[test]
fn all_parameters_are_required_dimensioned_and_finite() {
    let original = fixture();
    for name in original.parameters.keys() {
        let mut card = original.clone();
        card.parameters.remove(name);
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "missing {name}"
        );
        let mut card = original.clone();
        let parameter = card.parameters.get_mut(name).unwrap();
        parameter.dims = if parameter.dims == PA { Dims::NONE } else { PA };
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "units {name}"
        );
        let mut card = original.clone();
        card.parameters.get_mut(name).unwrap().value = f64::NAN;
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "finite {name}"
        );
    }
    for (name, invalid) in [
        ("shear_modulus_pa", 0.0),
        ("bulk_modulus_pa", 0.0),
        ("yield_stress_pa", -1.0),
        ("hardening_modulus_pa", -1.0),
        ("minimum_jacobian", 0.0),
        ("minimum_jacobian", 1.0),
        ("maximum_elastic_log_strain", 0.0),
        ("maximum_plastic_log_strain", 0.0),
        ("maximum_accumulated_plastic_strain", 0.0),
    ] {
        let mut card = original.clone();
        card.parameters.get_mut(name).unwrap().value = invalid;
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "range {name}"
        );
    }
}

#[test]
fn law_identity_history_and_extra_physics_cannot_be_reinterpreted() {
    for id in [
        "j2-plasticity",
        "finite-j2-hencky-cauchy",
        "finite-j2-hencky-kirchhoff-voce",
        "finite-j2-hencky-kirchhoff-thermoplastic",
    ] {
        let mut card = fixture();
        card.law = LawId(id.into());
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "{id}"
        );
    }
    let mut version = fixture();
    version.law_version = 2;
    let mut schema = fixture();
    schema.state_schema_version = 2;
    let mut declared = fixture();
    declared.initial_state = InitialStatePolicy::RequiresDeclaredState;
    for card in [version, schema, declared] {
        assert!(ResolvedFiniteJ2::try_from_model_card(&card).is_err());
    }
    for name in [
        "mu",
        "young_modulus_pa",
        "kinematic_hardening_modulus_pa",
        "voce_saturation_pa",
        "thermal_expansion_per_k",
        "temperature_evolution",
    ] {
        let mut card = fixture();
        card.parameters.insert(
            name.into(),
            LawParameter {
                value: 1.0,
                dims: PA,
            },
        );
        assert!(
            ResolvedFiniteJ2::try_from_model_card(&card).is_err(),
            "{name}"
        );
    }
}

#[test]
fn exact_typed_context_and_source_identity_are_retained() {
    let absolute = semantic(QuantityKind::AbsoluteTemperature);
    let pressure = semantic(QuantityKind::Pressure);
    let mut card = fixture();
    card.validity = card
        .validity
        .with_quantity("temperature", absolute, 290.0, 310.0)
        .with_quantity("environment_pressure", pressure, 1.0e5, 2.0e5);
    let context = QueryPoint::new()
        .with_quantity("temperature", absolute, 300.0)
        .unwrap()
        .with_quantity("environment_pressure", pressure, 1.5e5)
        .unwrap()
        .with("unconstrained_note_coordinate", 17.0)
        .unwrap();
    let source = ResolvedFiniteJ2::try_from_model_card_at(&card, &context).unwrap();
    assert_eq!(source.context(), &context);
    assert_eq!(source.model_card(), &card);
    assert!(ResolvedFiniteJ2::try_from_model_card(&card).is_err());
    for bad_context in [
        context.clone().with("temperature", 311.0).unwrap(),
        context.clone().with("environment_pressure", 0.9e5).unwrap(),
        context
            .clone()
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                300.0,
            )
            .unwrap(),
        context
            .clone()
            .with_quantity(
                "environment_pressure",
                semantic(QuantityKind::Stress),
                1.5e5,
            )
            .unwrap(),
    ] {
        assert!(ResolvedFiniteJ2::try_from_model_card_at(&card, &bad_context).is_err());
    }
    // Equal native coefficients do not collapse the complete source identity.
    let mut other = card.clone();
    other.provenance.source = "independent synthetic fixture provenance".into();
    let other = ResolvedFiniteJ2::try_from_model_card_at(&other, &context).unwrap();
    assert_eq!(source.law(), other.law());
    assert_ne!(
        source.model_card().content_hash(),
        other.model_card().content_hash()
    );
}

#[test]
fn mechanically_evolving_and_temperature_rate_axes_are_explicitly_refused() {
    for (axis, quantity, value) in [
        ("stress", semantic(QuantityKind::Stress), 1.0),
        (
            "strain",
            semantic(QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Normal,
            }),
            0.1,
        ),
        ("stretch", QuantitySpec::dimensional(Dims::NONE), 1.0),
        (
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            1.0,
        ),
        (
            "temperature_rate",
            QuantitySpec::dimensional(Dims([0, 0, -1, 1, 0, 0])),
            1.0,
        ),
        (
            "pressure_sample",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Instantaneous,
            )),
            1.0,
        ),
    ] {
        let mut card = fixture();
        card.validity = card.validity.with_quantity(axis, quantity, 0.0, 2.0);
        let context = QueryPoint::new()
            .with_quantity(axis, quantity, value)
            .unwrap();
        assert!(
            card.validity
                .contains_typed(context.axes(), context.axis_quantities())
        );
        let error = ResolvedFiniteJ2::try_from_model_card_at(&card, &context).unwrap_err();
        assert!(
            error.to_string().contains("fixed external condition"),
            "{axis}: {error}"
        );
    }
    let mut legacy = fixture();
    legacy.validity = legacy.validity.with("T", 290.0, 310.0);
    assert!(
        ResolvedFiniteJ2::try_from_model_card_at(
            &legacy,
            &QueryPoint::new().with("T", 300.0).unwrap(),
        )
        .is_ok()
    );
}
