//! G0/G1 source-driven hydrostatic compression and irreversible volume history.
//! Coefficients are authored numerical fixtures, not measured foam calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::porous::FinitePorousMaterial;
use fs_material::plastic::finite::porous::source::{FINITE_POROUS_LAW_ID, ResolvedFinitePorous};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const A: f64 = 1.2;

fn fixture() -> ConstitutiveModelCard {
    let source = fs_blake3::hash_bytes(b"authored finite porous compression software fixture");
    ConstitutiveModelCard {
        law: LawId(FINITE_POROUS_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", 80e6, PA),
            ("bulk_modulus_pa", K, PA),
            ("yield_stress_pa", Y, PA),
            ("hardening_modulus_pa", H, PA),
            ("pressure_sensitivity", A, Dims::NONE),
            ("minimum_jacobian", 0.2, Dims::NONE),
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
            source: "authored finite porous compression fixture".into(),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    }
}

fn dilation(stretch: f64) -> [f64; 9] {
    [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch]
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

#[test]
fn source_law_retains_hydrostatic_crush_storage_heat_and_stress_convention() {
    let card = fixture();
    let source = ResolvedFinitePorous::try_from_model_card(&card).unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.context(), &QueryPoint::new());
    assert!(std::ptr::eq(source.finite_porous(), source.law()));
    let law = source.finite_porous();
    let virgin = law.initial_state();
    let stretch = 0.8_f64;
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    // Independent hydrostatic return: beta*|tau_m|=Y+H*alpha,
    // ln(Jp)=-beta*alpha and tau_m=K*(ln(J)-ln(Jp)).
    let alpha = (-3.0 * K * beta * stretch.ln() - Y) / (H + K * beta * beta);
    let tau = -(Y + H * alpha) / beta;
    let response = law.response(&dilation(stretch), &virgin).unwrap();
    near(response.state.alpha, alpha, 1e-12);
    near(response.plastic_jacobian, (-beta * alpha).exp(), 1e-12);
    near(response.jacobian, stretch.powi(3), 1e-15);
    near(
        response.stored_energy_j_m3,
        0.5 * tau * tau / K + 0.5 * H * alpha * alpha,
        1e-6,
    );
    near(response.plastic_dissipation_j_m3, Y * alpha, 1e-6);
    for i in [0, 4, 8] {
        near(response.kirchhoff_pa[i], tau, 1e-6);
        near(response.first_piola_pa[i], tau / stretch, 1e-6);
    }
    // Cauchy stress differs substantially at this finite compression.
    assert!((tau / response.jacobian).abs() > 1.9 * tau.abs());
    assert_eq!(virgin, law.initial_state());

    let release = (-beta * alpha / 3.0).exp();
    let unloaded = law.response(&dilation(release), &response.state).unwrap();
    assert_eq!(unloaded.state, response.state);
    assert_eq!(unloaded.plastic_dissipation_j_m3, 0.0);
    for stress in unloaded.kirchhoff_pa {
        near(stress, 0.0, 1e-6);
    }
    assert!(
        unloaded.jacobian < 0.55,
        "permanent volume loss must survive unloading"
    );
    let reload = law
        .response(&dilation((stretch * release).sqrt()), &unloaded.state)
        .unwrap();
    assert_eq!(reload.state, unloaded.state);
    near(reload.kirchhoff_pa[0], 0.5 * tau, 1e-6);
}

#[test]
fn zero_sensitivity_and_zero_hardening_are_explicit_physical_choices() {
    let mut card = fixture();
    card.parameters
        .get_mut("pressure_sensitivity")
        .unwrap()
        .value = 0.0;
    let source = ResolvedFinitePorous::try_from_model_card(&card).unwrap();
    let law = source.law();
    let elastic = law.response(&dilation(0.8), &law.initial_state()).unwrap();
    assert_eq!(elastic.state, law.initial_state());
    near(elastic.kirchhoff_pa[0], 3.0 * K * 0.8_f64.ln(), 1e-6);
    assert_eq!(elastic.plastic_jacobian, 1.0);

    card.parameters
        .get_mut("pressure_sensitivity")
        .unwrap()
        .value = A;
    card.parameters
        .get_mut("hardening_modulus_pa")
        .unwrap()
        .value = 0.0;
    let source = ResolvedFinitePorous::try_from_model_card(&card).unwrap();
    let law = source.law();
    let plastic = law.response(&dilation(0.8), &law.initial_state()).unwrap();
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    near(plastic.kirchhoff_pa[0], -Y / beta, 1e-6);
    assert!(plastic.plastic_jacobian < 0.55);
}

#[test]
fn source_parameters_cannot_infer_units_law_history_or_missing_physics() {
    let original = fixture();
    for name in original.parameters.keys() {
        let mut missing = original.clone();
        missing.parameters.remove(name);
        assert!(
            ResolvedFinitePorous::try_from_model_card(&missing).is_err(),
            "missing {name}"
        );
        let mut wrong_units = original.clone();
        let value = wrong_units.parameters.get_mut(name).unwrap();
        value.dims = if value.dims == PA { Dims::NONE } else { PA };
        assert!(
            ResolvedFinitePorous::try_from_model_card(&wrong_units).is_err(),
            "units {name}"
        );
    }
    for id in [
        "porous-plasticity",
        "finite-j2-hencky-kirchhoff",
        "finite-porous-hencky-cauchy",
    ] {
        let mut card = original.clone();
        card.law = LawId(id.into());
        assert!(ResolvedFinitePorous::try_from_model_card(&card).is_err());
    }
    let mut version = original.clone();
    version.state_schema_version = 2;
    let mut history = original.clone();
    history.initial_state = InitialStatePolicy::RequiresDeclaredState;
    let mut extra = original.clone();
    extra.parameters.insert(
        "porosity".into(),
        LawParameter {
            value: 0.8,
            dims: Dims::NONE,
        },
    );
    let mut negative = original;
    negative
        .parameters
        .get_mut("pressure_sensitivity")
        .unwrap()
        .value = -1.0;
    for card in [version, history, extra, negative] {
        assert!(ResolvedFinitePorous::try_from_model_card(&card).is_err());
    }
}

#[test]
fn fixed_source_context_is_typed_and_evolving_stress_cannot_be_frozen() {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let stress = QuantitySpec::semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static));
    let mut card = fixture();
    card.validity = card
        .validity
        .with_quantity("temperature", temperature, 290.0, 310.0);
    let context = QueryPoint::new()
        .with_quantity("temperature", temperature, 300.0)
        .unwrap();
    let source = ResolvedFinitePorous::try_from_model_card_at(&card, &context).unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.context(), &context);
    assert!(ResolvedFinitePorous::try_from_model_card(&card).is_err());
    let outside = context.clone().with("temperature", 311.0).unwrap();
    assert!(ResolvedFinitePorous::try_from_model_card_at(&card, &outside).is_err());
    let wrong = context
        .clone()
        .with_quantity("temperature", QuantitySpec::dimensional(PA), 300.0)
        .unwrap();
    assert!(ResolvedFinitePorous::try_from_model_card_at(&card, &wrong).is_err());
    card.validity = card
        .validity
        .with_quantity("compressive_stress", stress, -1e9, 0.0);
    let dynamic = context
        .with_quantity("compressive_stress", stress, -1e6)
        .unwrap();
    assert!(ResolvedFinitePorous::try_from_model_card_at(&card, &dynamic).is_err());
}

#[test]
fn source_limits_reach_the_native_compression_and_history_owner() {
    for (name, value) in [
        ("minimum_jacobian", 0.6),
        ("maximum_accumulated_plastic_strain", 0.1),
    ] {
        let mut card = fixture();
        card.parameters.get_mut(name).unwrap().value = value;
        let source = ResolvedFinitePorous::try_from_model_card(&card).unwrap();
        let old = source.law().initial_state();
        assert!(
            source.law().response(&dilation(0.8), &old).is_err(),
            "{name}"
        );
        assert_eq!(old, source.law().initial_state());
    }
}

#[test]
fn explicit_voce_source_preserves_native_saturation_and_refuses_implicit_selection() {
    use fs_material::plastic::finite::porous::source::FINITE_POROUS_VOCE_LAW_ID;
    let mut card = fixture();
    card.parameters.insert(
        "voce_saturation_pa".into(),
        LawParameter {
            value: 12e6,
            dims: PA,
        },
    );
    card.parameters.insert(
        "voce_rate_per_strain".into(),
        LawParameter {
            value: 3.0,
            dims: Dims::NONE,
        },
    );
    assert!(ResolvedFinitePorous::try_from_model_card(&card).is_err());
    card.law = LawId(FINITE_POROUS_VOCE_LAW_ID.into());
    let resolved = ResolvedFinitePorous::try_from_model_card(&card).unwrap();
    let native = ResolvedFinitePorous::try_from_model_card(&fixture())
        .unwrap()
        .law()
        .clone()
        .with_voce_hardening(12e6, 3.0)
        .unwrap();
    assert_eq!(resolved.model_card(), &card);
    assert_eq!(resolved.law(), &native);
    assert_eq!(
        resolved
            .law()
            .response(&dilation(0.8), &native.initial_state())
            .unwrap(),
        native
            .response(&dilation(0.8), &native.initial_state())
            .unwrap()
    );
    for name in ["voce_saturation_pa", "voce_rate_per_strain"] {
        let mut missing = card.clone();
        missing.parameters.remove(name);
        assert!(ResolvedFinitePorous::try_from_model_card(&missing).is_err());
        let mut mistyped = card.clone();
        let parameter = mistyped.parameters.get_mut(name).unwrap();
        parameter.dims = if parameter.dims == PA { Dims::NONE } else { PA };
        assert!(ResolvedFinitePorous::try_from_model_card(&mistyped).is_err());
    }
}
