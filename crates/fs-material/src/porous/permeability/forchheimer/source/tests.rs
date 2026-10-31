//! Synthetic observation-backed cards exercise the original material owner.
//! These fixtures establish implementation behavior, not empirical validation.
use super::*;
use crate::gas::{GasSpec, GasState};
use crate::porous::permeability::{
    DarcyGasBarrier, INTRINSIC_PERMEABILITY_DIMS, IntrinsicPermeability,
};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};

const BETA: &str = "calibrated_inertial_beta";
const K: &str = "intrinsic_permeability";

fn support() -> PorousGasSupport {
    PorousGasSupport {
        temperature_k: (250.0, 600.0),
        pressure_pa: (10000.0, 1e6),
        maximum_superficial_speed_m_s: 20.0,
    }
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            150000.0,
        )
        .unwrap()
        .with_quantity(
            SOURCE_SPEED_AXIS,
            QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
            0.1,
        )
        .unwrap()
}
fn card(beta_value: PropertyValue, beta_dims: Dims, speed_bounds: (f64, f64)) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let provenance = Provenance {
        source: "synthetic observation-shaped Darcy–Forchheimer fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic rigid plug, normal x, no empirical claim".into(),
            method: "synthetic fit of dp/L = mu*q/k + rho*beta*abs(q)*q".into(),
            artifact: fs_blake3::hash_domain("fixture", b"synthetic k beta joint fit"),
            caveats:
                "synthetic source-retention test; uncertainty is not propagated as independent data"
                    .into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let domain = |temperature: (f64, f64), pressure: (f64, f64)| {
        ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
                temperature.0,
                temperature.1,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                pressure.0,
                pressure.1,
            )
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(K, INTRINSIC_PERMEABILITY_DIMS),
            value: PropertyValue::Scalar {
                value: 1e-10,
                dims: INTRINSIC_PERMEABILITY_DIMS,
            },
            validity: domain((280.0, 350.0), (100000.0, 250000.0)),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.03,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance.clone(),
        })
        .unwrap();
    let interpolation = if matches!(&beta_value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BETA, beta_dims),
            value: beta_value,
            validity: domain((290.0, 320.0), (120000.0, 200000.0)).with_quantity(
                SOURCE_SPEED_AXIS,
                QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
                speed_bounds.0,
                speed_bounds.1,
            ),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.04,
                confidence: 0.95,
            },
            interpolation,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic pore skeleton".into(),
            phase: "solid".into(),
            process: "synthetic pressure-loss characterization".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn scalar_card(beta: f64) -> MaterialCard {
    card(
        PropertyValue::Scalar {
            value: beta,
            dims: ForchheimerCalibration::BETA_DIMS,
        },
        ForchheimerCalibration::BETA_DIMS,
        (0.0, 0.25),
    )
}
fn resolve(card: &MaterialCard, beta_dims: Dims, temperature: f64) -> ResolvedMaterialStatePoint {
    let requirements = [(K, INTRINSIC_PERMEABILITY_DIMS), (BETA, beta_dims)].map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative).unwrap()
    });
    resolve_material_state_point(
        card,
        &point(temperature),
        &requirements,
        MaterialPropertySelection::PreferObservationBacked,
    )
    .unwrap()
}
fn layer(source: ResolvedMaterialStatePoint) -> PorousGasLayer {
    PorousGasLayer::from_material_isotropic(source, K, 0.02, None).unwrap()
}
fn bind(layer: &PorousGasLayer) -> Result<PorousGasLayer, MaterialError> {
    layer.with_material_forchheimer(
        BETA,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
    )
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 4e-12 * a.abs().max(b.abs()).max(1e-300),
        "{a} != {b}"
    );
}

#[test]
fn actual_source_beta_changes_gas_and_continuum_resistance_without_losing_claims() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150000.0).unwrap();
    let mut previous = f64::INFINITY;
    for beta in [0.0, 2e6, 4e6] {
        let card = scalar_card(beta);
        let source = resolve(&card, ForchheimerCalibration::BETA_DIMS, 300.0);
        let original = layer(source.clone());
        let calibrated = bind(&original).unwrap();
        let calibration = calibrated.forchheimer_calibration().unwrap();
        assert_eq!(calibration.beta_per_m(), beta);
        assert_eq!(calibration.source(), Some(&source));
        assert_eq!(calibration.source_property(), Some(BETA));
        assert_eq!(calibrated.source(), Some(&source));
        assert_eq!(calibrated.permeability(), original.permeability());
        let barrier =
            DarcyGasBarrier::try_new(vec![calibrated], [1.0, 0.0, 0.0], support()).unwrap();
        assert_eq!(barrier.support().temperature_k, (290.0, 320.0));
        assert_eq!(barrier.support().pressure_pa, (120000.0, 200000.0));
        assert_eq!(barrier.support().maximum_superficial_speed_m_s, 0.25);
        let flow = barrier.mass_flow(&gas, 150100.0, 149900.0, 0.01).unwrap();
        let flux = flow.mass_flux_left_to_right_kg_m2_s;
        assert!(flux < previous);
        previous = flux;
        let rt = gas.specific_gas_constant * gas.temperature;
        near(
            200.0 * 150000.0 / rt,
            (gas.dynamic_viscosity * 0.02 / 1e-10) * flux + (beta * 0.02) * flux * flux,
        );
        let (a, b) = barrier
            .volume_resistance_coefficients(&gas, 0.01)
            .unwrap()
            .unwrap();
        near(a, gas.dynamic_viscosity * 0.02 / 1e-10 / 0.01);
        near(b, gas.density * beta * 0.02 / 0.01_f64.powi(2));
        let continuum = barrier.continuum_geometry().unwrap();
        near(continuum.inertial_geometric_resistance(), beta * 0.02);
        assert_eq!(
            continuum.barrier().layers()[0]
                .forchheimer_calibration()
                .unwrap()
                .source(),
            Some(&source)
        );
        assert!(
            continuum
                .admit_flow(300.0, 150100.0, 149900.0, 0.3)
                .is_err()
        );
        assert!(barrier.mass_flow(&gas, 201000.0, 150000.0, 0.01).is_err());
        let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 325.0, 150000.0).unwrap();
        assert!(barrier.mass_flow(&hot, 150100.0, 149900.0, 0.01).is_err());
    }
}

#[test]
fn source_beta_rejects_other_quantities_curves_and_unsupported_speed_domains() {
    for dims in [
        Dims::NONE,
        INTRINSIC_PERMEABILITY_DIMS,
        fs_qty::Pressure::DIMS,
    ] {
        let card = card(
            PropertyValue::Scalar { value: 2e6, dims },
            dims,
            (0.0, 0.25),
        );
        let source = resolve(&card, dims, 300.0);
        assert!(bind(&layer(source)).is_err());
    }
    let dims = ForchheimerCalibration::BETA_DIMS;
    let changing = card(
        PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(290.0, 1e6), (320.0, 2e6)],
            dims,
        },
        dims,
        (0.0, 0.25),
    );
    assert!(bind(&layer(resolve(&changing, dims, 300.0))).is_err());
    let limited = card(
        PropertyValue::Scalar { value: 2e6, dims },
        dims,
        (0.05, 0.25),
    );
    assert!(bind(&layer(resolve(&limited, dims, 300.0))).is_err());
    let card = scalar_card(2e6);
    let only_k = resolve_material_state_point(
        &card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            K,
            INTRINSIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(bind(&layer(only_k)).is_err());
}

#[test]
fn source_identity_context_and_direction_cannot_be_silently_substituted() {
    let card = scalar_card(2e6);
    let source = resolve(&card, ForchheimerCalibration::BETA_DIMS, 300.0);
    let original = layer(source.clone());
    let calibration = ForchheimerCalibration::from_material(
        source.clone(),
        BETA,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
    )
    .unwrap();
    let other = scalar_card(3e6);
    assert!(
        layer(resolve(&other, ForchheimerCalibration::BETA_DIMS, 300.0))
            .with_forchheimer_calibration(calibration.clone())
            .is_err()
    );
    assert!(
        layer(resolve(&card, ForchheimerCalibration::BETA_DIMS, 310.0))
            .with_forchheimer_calibration(calibration.clone())
            .is_err()
    );
    let authored = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(1e-10).unwrap(),
        0.02,
        0.0,
        "explicit authored skeleton",
    )
    .unwrap();
    assert!(authored.with_forchheimer_calibration(calibration).is_err());
    let sideways = ForchheimerCalibration::from_material(
        source.clone(),
        BETA,
        [0.0, 1.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
    )
    .unwrap();
    let sideways = original.with_forchheimer_calibration(sideways).unwrap();
    assert!(DarcyGasBarrier::try_new(vec![sideways], [1.0, 0.0, 0.0], support()).is_err());
    assert_eq!(original.source(), Some(&source));
    assert!(original.forchheimer_calibration().is_none());
    let authored_beta = ForchheimerCalibration::try_new(
        2e6,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
        "explicit model input",
    )
    .unwrap();
    assert!(authored_beta.source().is_none());
    assert!(authored.with_forchheimer_calibration(authored_beta).is_ok());
}
