//! Original finite model, reference-area work and initiation source transport.
//! The fixture is authored synthetic data, not experimental material calibration.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::ductile::evolution::DuctileJ2Conditions;
use fs_material::ductile::evolution::source::{
    DuctileJ2SourceSettings, PLASTIC_WORK_PER_AREA_DIMS,
    POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
};
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
    ResolvedDuctileLocus, STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::source::{
    FINITE_POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY as WORK, ResolvedFiniteDuctileJ2,
    finite_ductile_requirements, resolve_finite_ductile_j2,
};
use fs_material::plastic::finite::source::{FINITE_J2_VOCE_LAW_ID, ResolvedFiniteJ2};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature(), t)
        .unwrap()
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("temperature", temperature(), 290.0, 310.0)
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic finite ductile verification".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn model() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(FINITE_J2_VOCE_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", 500.0, PA),
            ("bulk_modulus_pa", 1000.0, PA),
            ("yield_stress_pa", 3.0, PA),
            ("hardening_modulus_pa", 40.0, PA),
            ("minimum_jacobian", 0.1, Dims::NONE),
            ("maximum_elastic_log_strain", 0.5, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
            ("voce_saturation_pa", 15.0, PA),
            ("voce_rate_per_strain", 12.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: validity(),
        sources: vec![fs_blake3::hash_bytes(b"authored finite ductile test model")],
        provenance: provenance(),
    }
}
fn card(work_key: &str, work: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(work_key, PLASTIC_WORK_PER_AREA_DIMS),
            value: PropertyValue::Scalar {
                value: work,
                dims: PLASTIC_WORK_PER_AREA_DIMS,
            },
            validity: validity(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: STRESS_TRIAXIALITY_AXIS.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(-50.0, 0.002), (50.0, 0.002)],
                dims: Dims::NONE,
            },
            validity: validity()
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -50.0,
                    50.0,
                )
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(RATE),
                    1e-8,
                    100.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic finite ductile".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![model()],
    )
    .unwrap()
}
fn initiation(card: &MaterialCard) -> ResolvedDuctileLocus {
    resolve_ductile_locus(
        card,
        &point(300.0)
            .with_quantity(
                STRESS_TRIAXIALITY_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -50.0,
            )
            .unwrap()
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(RATE),
                0.01,
            )
            .unwrap(),
        &DuctileSourceRequest {
            temperature_axis: "temperature".into(),
            maximum_points: 2,
        },
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()
}
fn settings() -> DuctileJ2SourceSettings {
    DuctileJ2SourceSettings {
        characteristic_length_m: 0.02,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-10,
            yield_pa: 1e-8,
        },
    }
}
fn resolve(card: &MaterialCard) -> ResolvedFiniteDuctileJ2 {
    resolve_finite_ductile_j2(
        card,
        &card.models()[0],
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        initiation(card),
        settings(),
    )
    .unwrap()
}

#[test]
fn original_finite_voce_model_and_reference_work_drive_damage_at_nonunit_volume() {
    let card = card(WORK, 0.004);
    let source = resolve(&card);
    assert_eq!(source.material().card_identity(), card.content_hash());
    assert_eq!(source.plastic().model_card(), &card.models()[0]);
    assert_eq!(source.law().native(), source.plastic().law());
    assert_eq!(
        source.law().initiation_law().source().observations().len(),
        2
    );
    let law = source.law();
    let old = law.initial_state();
    let mut integrities = Vec::new();
    for dilation in [0.002_f64, 0.02] {
        let f = [
            (dilation + 0.04).exp(),
            0.0,
            0.0,
            0.0,
            (dilation - 0.02).exp(),
            0.0,
            0.0,
            0.0,
            (dilation - 0.02).exp(),
        ];
        let r = law
            .trial(
                &f,
                &old,
                DuctileJ2Conditions {
                    duration_s: 0.2,
                    temperature_k: 300.0,
                },
            )
            .unwrap();
        let original = source
            .plastic()
            .law()
            .response(&f, old.native_state())
            .unwrap();
        let alpha = original.state.alpha;
        let delta = alpha - 0.002;
        let work = delta * (3.0 + 40.0 * (0.002 + delta * 0.5))
            + 15.0 * (delta + (-0.024_f64).exp() * (-12.0 * delta).exp_m1() / 12.0);
        let g = (-0.02 * work / 0.004).exp();
        assert!((r.state.integrity() - g).abs() < 2e-12);
        assert!(r.state.damage() > 0.0);
        for i in 0..9 {
            assert!((r.first_piola_pa[i] - g * original.first_piola_pa[i]).abs() < 2e-10);
        }
        integrities.push(r.state.integrity());
    }
    // A changed current volume/area does not silently rescale reference Gp.
    assert!((integrities[0] - integrities[1]).abs() < 2e-12);
}

#[test]
fn small_strain_total_fracture_and_zero_work_cannot_alias_finite_reference_work() {
    for key in [
        POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
        "fracture_energy",
        WORK,
    ] {
        let source = card(key, if key == WORK { 0.0 } else { 0.004 });
        assert!(
            resolve_finite_ductile_j2(
                &source,
                &source.models()[0],
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
                initiation(&source),
                settings()
            )
            .is_err()
        );
    }
}

#[test]
fn numerically_identical_foreign_model_or_other_card_initiation_are_refused() {
    let source = card(WORK, 0.004);
    let mut foreign = source.models()[0].clone();
    foreign.provenance.source = "different original model, same coefficients".into();
    assert!(ResolvedFiniteJ2::try_from_model_card_at(&foreign, &point(300.0)).is_ok());
    assert!(
        resolve_finite_ductile_j2(
            &source,
            &foreign,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            initiation(&source),
            settings()
        )
        .is_err()
    );
    let other = card(WORK, 0.008);
    assert!(
        resolve_finite_ductile_j2(
            &source,
            &source.models()[0],
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            initiation(&other),
            settings()
        )
        .is_err()
    );
}

#[test]
fn complete_fixed_context_and_runtime_temperature_remain_original() {
    let source = card(WORK, 0.004);
    let material = resolve_material_state_point(
        &source,
        &point(300.0),
        &finite_ductile_requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    for context in [
        point(301.0),
        point(300.0)
            .with_quantity("pressure", QuantitySpec::dimensional(PA), 1e5)
            .unwrap(),
    ] {
        let plastic =
            ResolvedFiniteJ2::try_from_model_card_at(&source.models()[0], &context).unwrap();
        assert!(
            ResolvedFiniteDuctileJ2::try_from_finite_j2(
                &source,
                material.clone(),
                plastic,
                initiation(&source),
                settings()
            )
            .is_err()
        );
    }
    let resolved = resolve(&source);
    let old = resolved.law().initial_state();
    assert!(
        resolved
            .law()
            .trial(
                old.deformation(),
                &old,
                DuctileJ2Conditions {
                    duration_s: 0.2,
                    temperature_k: 301.0
                }
            )
            .is_err()
    );
    assert_eq!(old, resolved.law().initial_state());
}
