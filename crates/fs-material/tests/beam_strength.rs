//! Beam-face strength from original thermal and material snapshots.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::strength::tsai_wu::{
    PreparedLaminateTsaiWu, TsaiWuInteractionSelection,
};
use fs_material::laminate::strength::{
    LaminaFace, LaminaFailureMode, LaminaStrengthProperties, LaminaStrengthSelection,
    LaminateStrengthError, PreparedLaminateStrength,
};
use fs_material::laminate::thermal::beam::PreparedLaminateBeam;
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const NAMES: [&str; 5] = ["Xt", "Xc", "Yt", "Yc", "S12"];
const LIMITS: [f64; 5] = [80e6, 160e6, 50e6, 100e6, 20e6];

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn claim(name: &str, dims: Dims, low: f64, high: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (400.0, high)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 400.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured beam strength fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(alpha: f64, alternate_at_hot: bool) -> Arc<MaterialCard> {
    let mut set = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 5000.0),
        ("young_modulus_1", Pressure::DIMS, 100e9),
        ("young_modulus_2", Pressure::DIMS, 100e9),
        ("poisson_ratio_12", Dims::NONE, 0.25),
        ("shear_modulus_12", Pressure::DIMS, 40e9),
        ("shear_modulus_31", Pressure::DIMS, 40e9),
        ("shear_modulus_23", Pressure::DIMS, 40e9),
        (CTE, INVERSE_TEMPERATURE_DIMS, alpha),
        ("interaction", Dims::NONE, 0.2),
    ] {
        set.insert_claim(claim(name, dims, value, value)).unwrap();
    }
    for (name, value) in NAMES.into_iter().zip(LIMITS) {
        set.insert_claim(claim(name, Pressure::DIMS, value, 0.5 * value))
            .unwrap();
    }
    if alternate_at_hot {
        let mut other = claim("Xt", Pressure::DIMS, 1e9, 1e9);
        other.validity =
            ValidityDomain::unconstrained().with_quantity("T", temperature(), 400.0, 400.0);
        other.provenance.source = "different strength experiment, unavailable at reference".into();
        set.insert_claim(other).unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic beam".into(),
                phase: "solid".into(),
                process: "numerical oracle".into(),
                revision: 0,
            },
            set,
            vec![],
        )
        .unwrap(),
    )
}
fn selection() -> LaminaStrengthSelection {
    LaminaStrengthSelection {
        properties: LaminaStrengthProperties {
            tensile: [NAMES[0].into(), NAMES[2].into()],
            compressive: [NAMES[1].into(), NAMES[3].into()],
            shear: NAMES[4].into(),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn prepare(
    cards: &[Arc<MaterialCard>],
    angles: &[f64],
    temperatures: &[f64],
) -> PreparedLaminateStrength {
    let inputs = cards
        .iter()
        .zip(angles)
        .map(|(card, angle)| LaminaSourceInput {
            card: Arc::clone(card),
            reference: point(300.0),
            thickness_m: 0.005,
            orientation_rad: *angle,
            elasticity: LaminaElasticSource::OrthotropicPlane12,
            expansion: LaminaExpansionSource::Isotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    let source =
        LaminateThermalSource::try_new(inputs, ReferenceSurface::Midsurface, || false).unwrap();
    PreparedLaminateStrength::resolve(
        source.at_temperatures(temperatures, || false).unwrap(),
        &vec![selection(); cards.len()],
        || false,
    )
    .unwrap()
}
fn beam(strength: &PreparedLaminateStrength, width: f64, offset: f64) -> PreparedLaminateBeam {
    PreparedLaminateBeam::try_new(
        Arc::clone(strength.thermal().prepared()),
        width,
        offset,
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_zero_resultant_bimetal_still_exceeds_hot_layer_strength() {
    let strength = prepare(
        &[card(1e-5, false), card(3e-5, false)],
        &[0.0; 2],
        &[400.0; 2],
    );
    let strip = beam(&strength, 0.02, 0.0);
    let [eps, kappa] = strip.free_generalized_strain();
    let response = strip.evaluate(eps, kappa, || false).unwrap();
    near(response.axial_force_n, 0.0, 100.0);
    near(response.bending_moment_n_m, 0.0, 1.0);
    near(response.strain_energy_j_m, 2.5, 2.5);
    let assessed = strength.assess_beam(&strip, eps, kappa, || false).unwrap();
    assert_eq!((assessed.layer, assessed.face), (0, LaminaFace::Top));
    assert_eq!(
        assessed.assessment.critical_mode(),
        Some(LaminaFailureMode::LongitudinalTension)
    );
    // Equal E and thickness give +/-E*(alpha1-alpha0)*deltaT/2 at the interface.
    // The actual current tensile strength has decreased from 80 to 40 MPa.
    near(assessed.assessment.material_stress_pa()[0], 100e6, 100e6);
    near(assessed.assessment.utilization(), 2.5, 2.5);
    assert!(assessed.assessment.reaches_limit());
    let upper = strength.assess_beam_layer(&strip, 1, eps, kappa).unwrap();
    near(upper.bottom.material_stress_pa()[0], -100e6, 100e6);
    assert_eq!(
        upper.bottom.critical_mode(),
        Some(LaminaFailureMode::LongitudinalCompression)
    );
    let shifted = beam(&strength, 0.09, 0.002);
    let invariant = strength
        .assess_beam(&shifted, eps - 0.002 * kappa, kappa, || false)
        .unwrap();
    assert_eq!(
        (invariant.layer, invariant.face),
        (assessed.layer, assessed.face)
    );
    near(
        invariant.assessment.utilization(),
        assessed.assessment.utilization(),
        2.5,
    );
}

#[test]
fn g3_rotated_material_axes_keep_signed_strength_modes_and_face_ties() {
    let strength = prepare(
        &[card(0.0, false)],
        &[core::f64::consts::FRAC_PI_2],
        &[300.0],
    );
    let strip = beam(&strength, 0.02, 0.0);
    for (strain, mode, expected) in [
        (0.0005, LaminaFailureMode::TransverseTension, 1.0),
        (-0.0005, LaminaFailureMode::TransverseCompression, 0.5),
    ] {
        let result = strength.assess_beam(&strip, strain, 0.0, || false).unwrap();
        assert_eq!((result.layer, result.face), (0, LaminaFace::Bottom));
        assert_eq!(result.assessment.critical_mode(), Some(mode));
        near(
            result.assessment.material_stress_pa()[1],
            100e9 * strain,
            50e6,
        );
        near(result.assessment.material_stress_pa()[0], 0.0, 1.0);
        near(result.assessment.utilization(), expected, 1.0);
    }
}

#[test]
fn g4_equal_mechanics_cannot_pair_cold_strength_with_hot_or_detached_snapshots() {
    let cold = prepare(&[card(0.0, true)], &[0.0], &[300.0]);
    let old_claim = cold
        .layer_strength(0)
        .unwrap()
        .resolved()
        .property("Xt")
        .unwrap()
        .answer()
        .receipt
        .selected;
    let hot = cold.at_temperatures(&[400.0], || false).unwrap();
    let cold_beam = beam(&cold, 0.02, 0.0);
    let hot_beam = beam(&hot, 0.02, 0.0);
    assert_eq!(cold_beam.tangent(), hot_beam.tangent());
    assert_eq!(
        cold_beam.thermal_resultants(),
        hot_beam.thermal_resultants()
    );
    assert!(matches!(
        cold.assess_beam(&hot_beam, 0.0005, 0.0, || false),
        Err(LaminateStrengthError::ThermalSnapshotMismatch)
    ));
    assert!(matches!(
        hot.assess_beam_layer(&cold_beam, 0, 0.0005, 0.0),
        Err(LaminateStrengthError::ThermalSnapshotMismatch)
    ));
    let detached = PreparedLaminateBeam::try_new(
        Arc::new(cold.thermal().prepared().as_ref().clone()),
        0.02,
        0.0,
        || false,
    )
    .unwrap();
    assert_eq!(detached.source(), cold_beam.source());
    assert!(matches!(
        cold.assess_beam(&detached, 0.0005, 0.0, || false),
        Err(LaminateStrengthError::ThermalSnapshotMismatch)
    ));
    let original = cold
        .clone()
        .assess_beam(&cold_beam, 0.0005, 0.0, || false)
        .unwrap();
    let current = hot.assess_beam(&hot_beam, 0.0005, 0.0, || false).unwrap();
    near(
        current.assessment.utilization(),
        2.0 * original.assessment.utilization(),
        2.0,
    );
    let selected = hot
        .layer_strength(0)
        .unwrap()
        .resolved()
        .property("Xt")
        .unwrap();
    assert_eq!(selected.answer().receipt.selected, old_claim);
    assert_eq!(selected.value_si(), 40e6);
    assert_eq!(
        hot.layer_strength(0).unwrap().resolved().query(),
        point(400.0)
    );
    assert!(
        PreparedLaminateStrength::resolve(hot.thermal().clone(), &[selection()], || false).is_err()
    );
}

#[test]
fn g1_tsai_wu_uses_the_same_rotated_beam_faces_and_original_interaction() {
    let maximum = prepare(
        &[card(0.0, false)],
        &[core::f64::consts::FRAC_PI_4],
        &[300.0],
    );
    let strip = beam(&maximum, 0.02, 0.0);
    let coupled = PreparedLaminateTsaiWu::from_maximum_stress(
        maximum,
        &[TsaiWuInteractionSelection {
            property: "interaction".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        || false,
    )
    .unwrap();
    let result = coupled.assess_beam(&strip, 0.0003, 0.0, || false).unwrap();
    assert_eq!((result.layer, result.face), (0, LaminaFace::Bottom));
    let stress: [f64; 3] = [15e6, 15e6, -15e6];
    for (actual, expected) in result
        .assessment
        .material_stress_pa()
        .into_iter()
        .zip(stress)
    {
        near(actual, expected, 15e6);
    }
    let [s1, s2, shear] = stress;
    let a = s1 * s1 / (80e6 * 160e6)
        + s2 * s2 / (50e6 * 100e6)
        + 0.4 * s1 * s2 / ((80e6_f64 * 160e6).sqrt() * (50e6_f64 * 100e6).sqrt())
        + (shear / 20e6).powi(2);
    let b = s1 * (1.0 / 80e6 - 1.0 / 160e6) + s2 * (1.0 / 50e6 - 1.0 / 100e6);
    near(
        result.assessment.utilization(),
        0.5 * (b + (b * b + 4.0 * a).sqrt()),
        1.0,
    );
    let faces = coupled.assess_beam_layer(&strip, 0, 0.0003, 0.0).unwrap();
    assert_eq!(faces.bottom, faces.top);
    let hot = coupled.at_temperatures(&[400.0], || false).unwrap();
    assert!(matches!(
        hot.assess_beam(&strip, 0.0003, 0.0, || false),
        Err(LaminateStrengthError::ThermalSnapshotMismatch)
    ));
    let hot_beam =
        PreparedLaminateBeam::try_new(Arc::clone(hot.thermal().prepared()), 0.02, 0.0, || false)
            .unwrap();
    near(
        hot.assess_beam(&hot_beam, 0.0003, 0.0, || false)
            .unwrap()
            .assessment
            .utilization(),
        2.0 * result.assessment.utilization(),
        2.0,
    );
}

#[test]
fn g4_beam_strength_cancellation_and_invalid_recovery_publish_no_partial_result() {
    let strength = prepare(
        &[card(1e-5, false), card(3e-5, false)],
        &[0.0; 2],
        &[400.0; 2],
    );
    let strip = beam(&strength, 0.02, 0.0);
    let mut polls = 0;
    assert!(matches!(
        strength.assess_beam(&strip, 0.0, 0.0, || {
            polls += 1;
            polls == 3
        }),
        Err(LaminateStrengthError::Cancelled)
    ));
    assert!(strength.assess_beam_layer(&strip, 2, 0.0, 0.0).is_err());
    assert!(
        strength
            .assess_beam(&strip, f64::NAN, 0.0, || false)
            .is_err()
    );
    let free = strip.free_generalized_strain();
    assert!(
        strength
            .assess_beam(&strip, free[0], free[1], || false)
            .unwrap()
            .assessment
            .reaches_limit()
    );
}
