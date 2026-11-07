//! Same-source finite coefficient curvature reaches actual capacity and waves.

#[path = "finite_thermoelastic_coefficients/heating.rs"]
mod heating;

use fs_couple::finite_solid::FiniteSolidMaterial;
use fs_couple::finite_solid::thermal::thermodynamics::{
    CaloricReference, FiniteThermoelasticSourceError,
};
use fs_couple::finite_solid::thermal::{
    ExpansionThermalFrame, ExpansionThermalReference, ThermalReference,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::thermoelastic::ThermalReferenceJet;
use fs_material::hyper::thermoelastic::source::ResolvedNeoHookeanThermalCoefficients;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const F: [f64; 9] = [1.22, 0.16, 0.04, -0.05, 0.91, 0.12, 0.03, -0.06, 1.13];
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn card(slopes: [f64; 2], cp_scale: f64) -> MaterialCard {
    card_with_expansion(slopes, cp_scale, [-0.01, 0.0, 0.01])
}
fn card_with_expansion(slopes: [f64; 2], cp_scale: f64, expansion: [f64; 3]) -> MaterialCard {
    let artifact = fs_blake3::hash_bytes(b"authored original finite thermal coefficient consumer");
    let provenance = Provenance {
        source: "authored E/nu, Cp and engineering expansion; not measurements".into(),
        license: "test-only".into(),
        artifact: Some(artifact),
    };
    let validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        *point(300.0).axis_quantities().get("T").unwrap(),
        280.0,
        320.0,
    );
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            [3000.0 - 10.0 * slopes[0], 3000.0, 3000.0 + 10.0 * slopes[0]],
        ),
        (
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            [0.25 - 10.0 * slopes[1], 0.25, 0.25 + 10.0 * slopes[1]],
        ),
        ("density", Density::DIMS, [6.0; 3]),
        (
            "specific_heat_capacity",
            CP_DIMS,
            [10.0 * cp_scale, 12.0 * cp_scale, 16.0 * cp_scale],
        ),
        ("engineering_length_change", Dims::NONE, expansion),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [290.0, 300.0, 310.0].into_iter().zip(values).collect(),
                    dims,
                },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let model = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", 1200.0, Pressure::DIMS),
            ("lambda", 1200.0, Pressure::DIMS),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity,
        sources: vec![artifact],
        provenance,
    };
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored finite thermal solid".into(),
            phase: "solid".into(),
            process: "source-qualified analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap()
}
fn material(card: &MaterialCard) -> FiniteSolidMaterial {
    let expansion = PreparedFreeExpansion::try_new(
        card,
        FreeExpansionGeometryQuery {
            reference: &point(300.0),
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &PropertyKey::new("engineering_length_change", Dims::NONE),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(310.0),
        },
    )
    .unwrap();
    FiniteSolidMaterial::resolve_thermal(
        card,
        0,
        ThermalReference::Expansion(
            ExpansionThermalReference::try_new(
                card,
                expansion,
                300.0,
                ExpansionThermalFrame::Isotropic,
            )
            .unwrap(),
        ),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn coefficients(
    card: &MaterialCard,
    material: &FiniteSolidMaterial,
) -> ResolvedNeoHookeanThermalCoefficients {
    let source = resolve_material_state_point(
        card,
        material.law().context(),
        &[
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: -1.0,
                    upper: 0.5,
                },
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedNeoHookeanThermalCoefficients::try_new(card, material.law(), &source, [290.0, 310.0])
        .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(310.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        0.01,
    )
    .unwrap()
}
fn h_at_300(source: &ResolvedHeatCapacityCurve) -> f64 {
    source
        .curve()
        .knots()
        .iter()
        .find(|knot| knot.temperature_k == 300.0)
        .unwrap()
        .specific_enthalpy_j_kg
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn actual_finite_solid_uses_original_variable_cp_and_full_coefficient_response() {
    let card = card([-12.0, 0.004], 1.0);
    let material = material(&card);
    let coefficients = coefficients(&card, &material);
    let caloric = caloric(&card);
    let h = h_at_300(&caloric);
    let held = material
        .adiabatic_acoustic_state(
            &F,
            ReferenceDeformation::identity(),
            &caloric,
            h,
            CaloricReference::StressFreeZeroPressure,
        )
        .unwrap();
    let state = material
        .adiabatic_acoustic_state_with_coefficients(
            &F,
            ReferenceDeformation::identity(),
            &coefficients,
            &caloric,
            h,
            CaloricReference::StressFreeZeroPressure,
        )
        .unwrap();
    assert!(core::ptr::eq(state.material(), &material));
    assert!(core::ptr::eq(state.caloric_source(), &caloric));
    assert!(core::ptr::eq(
        state.thermal_coefficients().unwrap(),
        &coefficients
    ));
    assert!(held.thermal_coefficients().is_none());
    assert_eq!(
        state.thermodynamics().first_piola_pa(),
        held.thermodynamics().first_piola_pa()
    );
    assert_eq!(
        state.thermodynamics().isothermal_tangent_pa(),
        held.thermodynamics().isothermal_tangent_pa()
    );
    near(state.native().caloric_heat_capacity_j_m3_k(), 6.0 * 12.0);
    near(
        state.native().fixed_deformation_heat_capacity_j_m3_k(),
        6.0 * 12.0 + state.thermodynamics().mechanical_heat_capacity_j_m3_k(),
    );
    assert!(
        (state.native().fixed_deformation_heat_capacity_j_m3_k()
            - held.native().fixed_deformation_heat_capacity_j_m3_k())
        .abs()
            > 0.1
    );
    let jet = ThermalReferenceJet::try_new(
        ReferenceDeformation::identity(),
        I.map(|v| 0.001 * v),
        [0.0; 9],
    )
    .unwrap();
    let direct = coefficients
        .thermoelastic_state(&F, jet)
        .unwrap()
        .adiabatic_acoustics(6.0, 72.0)
        .unwrap();
    near(
        state.native().fixed_deformation_heat_capacity_j_m3_k(),
        direct.fixed_deformation_heat_capacity_j_m3_k(),
    );
    let waves = state.directional_waves([0.3, -0.4, 0.5]).unwrap();
    let expected = direct
        .acoustic()
        .directional_waves(waves.direction)
        .unwrap();
    let baseline = held.directional_waves(waves.direction).unwrap();
    for (actual, expected) in waves.modes.iter().zip(expected.modes) {
        near(
            actual.squared_phase_speed_m2_s2,
            expected.squared_phase_speed_m2_s2,
        );
    }
    assert!(
        (waves.modes[2].squared_phase_speed_m2_s2 - baseline.modes[2].squared_phase_speed_m2_s2)
            .abs()
            > 1e-3
    );
}

#[test]
fn actual_finite_solid_refuses_foreign_coefficient_receipts_and_nonpositive_capacity() {
    let card = card([-12.0, 0.004], 1.0);
    let material = material(&card);
    let coefficients = coefficients(&card, &material);
    let foreign_card = self::card([0.0; 2], 1.0);
    let foreign = self::material(&foreign_card);
    let caloric = caloric(&foreign_card);
    assert!(matches!(
        foreign.adiabatic_acoustic_state_with_coefficients(
            &F,
            ReferenceDeformation::identity(),
            &coefficients,
            &caloric,
            h_at_300(&caloric),
            CaloricReference::StressFreeZeroPressure
        ),
        Err(FiniteThermoelasticSourceError::Invalid(_))
    ));
    let low_card = self::card([-12.0, 0.004], 1e-9);
    let low = self::material(&low_card);
    let low_coefficients = self::coefficients(&low_card, &low);
    let low_caloric = self::caloric(&low_card);
    assert!(
        low.adiabatic_acoustic_state_with_coefficients(
            &F,
            ReferenceDeformation::identity(),
            &low_coefficients,
            &low_caloric,
            h_at_300(&low_caloric),
            CaloricReference::StressFreeZeroPressure
        )
        .is_err()
    );
    assert_eq!(material.law(), coefficients.model());
}
