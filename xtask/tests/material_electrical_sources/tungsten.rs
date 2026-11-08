//! Original NBS tungsten data on native electrical, caloric and surface owners.
//! This checks source transport and numerical response, not specimen validation.

use super::{close, compile};
use fs_conduction::SurfaceEmissivity;
use fs_matdb::{ClaimSelection, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, ElectricalPropertyBasis, OhmicDrive,
    resolve_conductor_state_with_property,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const SLUG: &str = "tungsten-polished-nbs-1971";
const TEMPERATURE: &str = "T_IPTS68";
// Nominal 1e-4 torr converted to Pa; these digits are not pressure accuracy.
const PRESSURE: f64 = 0.013332236842105263;
const MOLAR_MASS: f64 = 0.18385; // The paper's original calorimetry basis [kg/mol].
const POLICY: SelectionPolicy = SelectionPolicy::SingleClaimOnly;
const MOLAR_CP: [f64; 17] = [
    31.65, 32.49, 33.29, 34.08, 34.89, 35.72, 36.61, 37.57, 38.63, 39.81, 41.14, 42.62, 44.29,
    46.17, 48.27, 50.63, 53.25,
];

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), PRESSURE)
        .unwrap()
}

fn cp_key() -> PropertyKey {
    PropertyKey::new("specific-heat-capacity-cp", Dims([2, 0, -2, -1, 0, 0]))
}

fn reference_electrical_key() -> PropertyKey {
    PropertyKey::new(
        "electrical-resistivity-reference-geometry-298k",
        ELECTRICAL_RESISTIVITY_DIMS,
    )
}

fn reference_emittance_key() -> PropertyKey {
    PropertyKey::new(
        "hemispherical-total-emissivity-reference-area-298k",
        Dims::NONE,
    )
}

#[test]
fn source_fit_drives_native_reference_geometry_resistance_and_signed_joule_heat() {
    let pack = compile(SLUG, 5);
    let card = pack.card();
    let original = resolve_conductor_state_with_property(
        card,
        &point(2000.0),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
        &reference_electrical_key(),
    )
    .unwrap();
    // These are declared 298 K reference dimensions, as required by Eq.2.
    let node = MaterialOhmicNode::try_new(card, original, TEMPERATURE, [2000.0, 3600.0], 0.1, 1e-8)
        .unwrap();
    for t in [2000.0, 2300.0, 2500.5, 2856.0, 3000.0, 3600.0] {
        let conductor = node.resolve_at_temperature_k(t).unwrap();
        let resistance = (-14.08 + 0.03515 * t) * 1e-8 * 0.1 / 1e-8;
        close(conductor.resistance_ohm(), resistance);
        assert_eq!(conductor.material().query(), point(t));
        let receipt = &conductor.source_property().answer().receipt;
        assert!(receipt.observation_backed && !receipt.source_hashes.is_empty());
        card.claims().verify_receipt(receipt).unwrap();
        for drive in [
            OhmicDrive::Current(2.0),
            OhmicDrive::Current(-2.0),
            OhmicDrive::Voltage(12.0),
            OhmicDrive::Voltage(-12.0),
        ] {
            let current = conductor.current_a_for_drive(drive).unwrap();
            let expected_current = match drive {
                OhmicDrive::Current(i) => i,
                OhmicDrive::Voltage(v) => v / resistance,
            };
            close(current, expected_current);
            close(
                conductor.joule_power_w(current).unwrap(),
                expected_current * expected_current * resistance,
            );
        }
    }
    close(
        node.maximum_positive_joule_temperature_slope_w_k(OhmicDrive::Current(2.0))
            .unwrap(),
        4.0 * 0.03515 * 1e-8 * 0.1 / 1e-8,
    );
    for t in [293.0, 1999.0, 3601.0] {
        assert!(node.resolve_at_temperature_k(t).is_err());
    }
    assert!(
        resolve_conductor_state_with_property(
            card,
            &point(2800.0).with("p", 101325.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
            &reference_electrical_key(),
        )
        .is_err()
    );
    // Reference geometry must be an explicit source choice. Generic intrinsic
    // resistivity lookup cannot silently attach changing hot dimensions.
    assert!(
        card.claims()
            .query_typed(&super::electrical_key(), &point(2800.0), POLICY)
            .is_err()
    );
}

#[test]
fn source_cp_integrates_and_inverts_without_inventing_hot_density_or_conductivity() {
    let pack = compile(SLUG, 5);
    let card = pack.card();
    let key = cp_key();
    let lower = point(2000.0);
    let upper = point(3600.0);
    let query = IsobaricHeatCapacityQuery {
        specific_heat_cp: &key,
        temperature_axis: TEMPERATURE,
        pressure_axis: "p",
        lower: &lower,
        upper: &upper,
        selection: ClaimSelection::Policy(POLICY),
    };
    let caloric = IsobaricHeatCapacityLaw::try_new(card, query, 2000.0).unwrap();
    assert_eq!(caloric.source_line().knots.len(), 17);
    assert_eq!(caloric.material_card_identity(), card.content_hash());
    assert_eq!(caloric.pressure_pa(), PRESSURE);
    for (i, molar_cp) in MOLAR_CP.iter().enumerate() {
        let temperature = 2000.0 + 100.0 * i as f64;
        close(
            caloric.specific_heat_cp_j_kg_k(temperature).unwrap(),
            molar_cp / MOLAR_MASS,
        );
        let receipt = &caloric.source_line().knots[i].1.receipt;
        card.claims().verify_receipt(receipt).unwrap();
        assert!(receipt.observation_backed);
        assert_eq!(
            caloric.source_line().knots[i].1.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    // Independent sums of trapezoids of the printed molar rows, then /0.18385.
    // The native owner integrates the admitted affine-Cp law, not the source cubic.
    close(
        caloric
            .specific_enthalpy_change_j_kg(2000.0, 3600.0)
            .unwrap(),
        347381.01713353277,
    );
    let heat = caloric
        .specific_enthalpy_change_j_kg(2300.0, 3000.0)
        .unwrap();
    close(heat, 141876.52977971172);
    close(
        caloric
            .temperature_after_specific_heat(2300.0, heat)
            .unwrap(),
        3000.0,
    );
    close(
        caloric
            .temperature_after_specific_heat(3000.0, -heat)
            .unwrap(),
        2300.0,
    );
    assert!(
        caloric
            .temperature_after_specific_heat(3600.0, 1.0)
            .is_err()
    );
    assert!(caloric.specific_heat_cp_j_kg_k(1999.0).is_err());
    let wrong_pressure = lower.clone().with("p", 101325.0).unwrap();
    assert!(
        IsobaricHeatCapacityLaw::try_new(
            card,
            IsobaricHeatCapacityQuery {
                lower: &wrong_pressure,
                ..query
            },
            2000.0
        )
        .is_err()
    );
    let density = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
    let room = QueryPoint::new()
        .with_quantity(
            TEMPERATURE,
            quantity(QuantityKind::AbsoluteTemperature),
            293.0,
        )
        .unwrap();
    close(
        card.claims()
            .query_typed(&density, &room, POLICY)
            .unwrap()
            .evidence
            .value
            .value,
        19230.0,
    );
    assert!(
        card.claims()
            .query_typed(&density, &point(2800.0), POLICY)
            .is_err()
    );
    assert!(
        card.claims()
            .query_typed(&super::conductivity_key(), &point(2800.0), POLICY)
            .is_err()
    );
}

#[test]
fn total_and_normal_spectral_emittance_keep_their_distinct_domains_and_roles() {
    let pack = compile(SLUG, 5);
    let card = pack.card();
    for t in [2300.0, 2500.0, 2856.0, 3400.0] {
        let context = point(t);
        let answer = card
            .claims()
            .query_typed(&reference_emittance_key(), &context, POLICY)
            .unwrap();
        close(answer.evidence.value.value, 0.2627 + 2.770e-5 * t);
        card.claims().verify_receipt(&answer.receipt).unwrap();
        // Native generic surface ingress must not reinterpret this reference-area
        // coefficient as an intrinsic emissivity multiplying an expanded hot area.
        assert!(
            SurfaceEmissivity::from_condition(
                "NBS polished tube",
                ThermalMaterialPoint {
                    material: card,
                    context: &context,
                    temperature_axis: TEMPERATURE,
                    liquid_fraction_axis: None,
                },
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    // The paper explicitly labels the lower-temperature total-emittance entries
    // extrapolated. Neither they nor a normal sample extend the thermal source.
    for t in [2000.0, 2200.0, 3401.0] {
        let context = point(t);
        assert!(
            card.claims()
                .query_typed(&reference_emittance_key(), &context, POLICY)
                .is_err()
        );
    }
    let key = PropertyKey::new("normal-spectral-emissivity", Dims::NONE);
    let normal = point(2500.0)
        .with_quantity(
            "effective-wavelength",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            650e-9,
        )
        .unwrap()
        .with_quantity(
            "direction-cosine",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    let answer = card.claims().query_typed(&key, &normal, POLICY).unwrap();
    close(answer.evidence.value.value, 0.3804 - 5.060e-7 * 2500.0);
    assert!((answer.evidence.value.value - (0.2627 + 2.770e-5 * 2500.0)).abs() > 0.04);
    card.claims().verify_receipt(&answer.receipt).unwrap();
    for unsupported in [
        point(2500.0),
        normal.clone().with("effective-wavelength", 645e-9).unwrap(),
        normal.clone().with("direction-cosine", 0.5).unwrap(),
        normal.clone().with(TEMPERATURE, 3001.0).unwrap(),
    ] {
        assert!(
            card.claims()
                .query_typed(&key, &unsupported, POLICY)
                .is_err()
        );
    }
}
