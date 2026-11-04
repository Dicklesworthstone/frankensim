//! Analytic physics fixtures, not measured material qualifications.
use super::*;
use crate::thermal::transport::{
    SPECIFIC_HEAT,
    expansion::{FreeExpansionSourceQuery, FreeExpansionStrainSources},
    isobaric::IsobaricHeatCapacityQuery,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const PRESSURE: f64 = 100_000.0;
const CP: [(f64, f64); 3] = [(290.0, 2.0), (300.0, 2.0), (310.0, 4.0)];
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), pressure)
        .unwrap()
}
fn keys() -> [PropertyKey; 5] {
    [
        PropertyKey::new("cp", SPECIFIC_HEAT),
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("free-length-x", Dims::NONE),
        PropertyKey::new("free-length-y", Dims::NONE),
        PropertyKey::new("free-length-z", Dims::NONE),
    ]
}
fn claim(key: &PropertyKey, knots: &[(f64, f64)], pressure_qualified: bool) -> PropertyClaim {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        absolute(),
        knots[0].0,
        knots[knots.len() - 1].0,
    );
    if pressure_qualified {
        validity = validity.with_quantity(
            "p",
            QuantitySpec::dimensional(Pressure::DIMS),
            PRESSURE,
            2.0 * PRESSURE,
        );
    }
    PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: knots.to_vec(),
            dims: key.quantity().dims(),
        },
        validity,
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored Cp and free-expansion physics fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn fixture(
    capacity: &[(f64, f64)],
    strains: &[&[(f64, f64)]],
    pressure_qualified: bool,
    datum: f64,
) -> (
    MaterialCard,
    IsobaricHeatCapacityLaw,
    PreparedFreeExpansionDensity,
) {
    let keys = keys();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(&keys[0], capacity, true))
        .unwrap();
    let mut density = claim(
        &keys[1],
        &[(300.0, 1000.0), (300.0, 1000.0)],
        pressure_qualified,
    );
    density.value = PropertyValue::Scalar {
        value: 1000.0,
        dims: keys[1].quantity().dims(),
    };
    density.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claims.insert_claim(density).unwrap();
    for (i, knots) in strains.iter().enumerate() {
        claims
            .insert_claim(claim(&keys[2 + i], knots, pressure_qualified))
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "pure caloric specimen".into(),
            phase: "solid".into(),
            process: "maintained-pressure free expansion".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let (caloric, density) = prepare(
        &card,
        strains.len(),
        capacity[0].0,
        capacity[capacity.len() - 1].0,
        datum,
        PRESSURE,
    );
    (card, caloric, density)
}
fn prepare(
    card: &MaterialCard,
    directions: usize,
    low: f64,
    high: f64,
    datum: f64,
    pressure: f64,
) -> (IsobaricHeatCapacityLaw, PreparedFreeExpansionDensity) {
    let keys = keys();
    let (lower, upper, reference) = (
        point(low, pressure),
        point(high, pressure),
        point(300.0, pressure),
    );
    let caloric = IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &keys[0],
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &lower,
            upper: &upper,
            selection: POLICY,
        },
        datum,
    )
    .unwrap();
    let strains = if directions == 1 {
        FreeExpansionStrainSources::Isotropic {
            strain: &keys[2],
            selection: POLICY,
        }
    } else {
        FreeExpansionStrainSources::FixedPrincipal {
            frame: ContentHash([7; 32]),
            strains: [&keys[2], &keys[3], &keys[4]],
            selections: [POLICY; 3],
        }
    };
    let density = PreparedFreeExpansionDensity::try_new(
        card,
        FreeExpansionSourceQuery {
            density: &keys[1],
            density_selection: POLICY,
            reference: &reference,
            strains,
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
        },
    )
    .unwrap();
    (caloric, density)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn affine_cp_heating_and_cooling_close_enthalpy_work_internal_energy_and_entropy() {
    let (card, caloric, density) = fixture(&CP, &[&[(290.0, 0.0), (310.0, 0.02)]], true, 300.0);
    let law = IsobaricFreeExpansion::try_new(&card, &caloric, &density).unwrap();
    let step = law.advance_heat(295.0, 2.0, 45.0).unwrap();
    close(step.state.temperature_k(), 305.0, 1e-12);
    // Cp=2 below 300 K, and Cp=.2*T-58 above it; original reference strain is .01.
    let [s0, s1] = [1.005_f64 / 1.01, 1.015_f64 / 1.01];
    let volume_change = 0.002 * (s1.powi(3) - s0.powi(3));
    close(step.volume_change_m3, volume_change, 2e-18);
    close(step.enthalpy_change_j, 45.0, 2e-12);
    close(step.boundary_work_j, PRESSURE * volume_change, 2e-12);
    close(
        step.internal_energy_change_j + step.boundary_work_j,
        45.0,
        2e-12,
    );
    let entropy =
        2.0 * (2.0 * (300.0_f64 / 295.0).ln() + 0.2 * 5.0 - 58.0 * (305.0_f64 / 300.0).ln());
    close(step.entropy_change_j_k, entropy, 2e-14);
    assert!(step.heat_residual_j.abs() < 2e-12);
    close(
        step.state.density_kg_m3() * step.final_volume_m3,
        2.0,
        1e-15,
    );
    let cooled = law
        .advance_heat(step.state.temperature_k(), 2.0, -45.0)
        .unwrap();
    close(cooled.state.temperature_k(), 295.0, 1e-12);
    close(cooled.boundary_work_j, -step.boundary_work_j, 2e-12);
    close(cooled.entropy_change_j_k, -step.entropy_change_j_k, 2e-14);
    let unchanged = law.advance_heat(300.0, 2.0, 0.0).unwrap();
    assert_eq!(unchanged.initial, unchanged.state);
    assert_eq!(unchanged.internal_energy_change_j, 0.0);
    assert_eq!(unchanged.boundary_work_j, 0.0);
    assert_eq!(unchanged.entropy_change_j_k, 0.0);
    assert_eq!(unchanged.path_identity, law.identity());

    let (card, caloric, density) = fixture(&CP, &[&[(290.0, 0.01), (310.0, 0.01)]], true, 300.0);
    let rigid = IsobaricFreeExpansion::try_new(&card, &caloric, &density)
        .unwrap()
        .advance_heat(295.0, 2.0, 45.0)
        .unwrap();
    assert_eq!(rigid.boundary_work_j, 0.0);
    assert_eq!(rigid.internal_energy_change_j, rigid.enthalpy_change_j);
}

#[test]
fn anisotropic_original_lengths_and_source_breaks_determine_volume_work() {
    let (card, caloric, density) = fixture(
        &CP,
        &[
            &[(290.0, -0.01), (310.0, 0.03)],
            &[(290.0, 0.03), (310.0, -0.01)],
            &[(290.0, 0.0), (302.0, 0.024), (310.0, 0.016)],
        ],
        true,
        290.0,
    );
    let law = IsobaricFreeExpansion::try_new(&card, &caloric, &density).unwrap();
    let step = law.advance_heat(295.0, 2.0, 64.8).unwrap();
    close(step.state.temperature_k(), 308.0, 1e-12);
    let before = (1.0 / 1.01) * (1.02 / 1.01) * (1.01 / 1.02);
    let after = (1.026 / 1.01) * (0.994 / 1.01) * (1.018 / 1.02);
    close(step.volume_change_m3, 0.002 * (after - before), 2e-18);
    close(
        step.internal_energy_change_j + PRESSURE * step.volume_change_m3,
        64.8,
        2e-12,
    );
    assert_eq!(
        step.state,
        density.at_temperature(step.state.temperature_k()).unwrap()
    );
    assert_ne!(
        step.state.principal_stretches()[0],
        step.state.principal_stretches()[1]
    );
}

#[test]
fn tiny_heat_and_expansion_survive_a_large_enthalpy_datum_and_rounded_density() {
    let (card, caloric, density) = fixture(
        &[(1.0, 1e20), (290.0, 2.0), (300.0, 2.0), (310.0, 4.0)],
        &[&[(1.0, -0.289), (310.0, 0.02)]],
        true,
        1.0,
    );
    let law = IsobaricFreeExpansion::try_new(&card, &caloric, &density).unwrap();
    let final_k = 300.0_f64.next_up();
    let dt = final_k - 300.0;
    let heat = 4.0 * dt;
    assert_eq!(
        caloric.specific_enthalpy_j_kg(300.0).unwrap(),
        caloric.specific_enthalpy_j_kg(final_k).unwrap()
    );
    let step = law.advance_heat(300.0, 2.0, heat).unwrap();
    assert_eq!(step.state.temperature_k(), final_k);
    assert_eq!(step.initial.density_kg_m3(), step.state.density_kg_m3());
    let stretch_change = (0.001 / 1.01) * dt;
    let expected_volume =
        0.002 * (3.0 * stretch_change + 3.0 * stretch_change.powi(2) + stretch_change.powi(3));
    close(
        step.volume_change_m3,
        expected_volume,
        expected_volume.abs() * 2e-12,
    );
    assert!(step.boundary_work_j > 0.0);
    close(step.enthalpy_change_j, heat, heat.abs() * 2e-14);
    close(step.entropy_change_j_k, heat / 300.0, heat.abs() * 1e-15);
    close(
        step.internal_energy_change_j + step.boundary_work_j,
        heat,
        heat.abs() * 2e-14,
    );
    // Opposite principal expansions cancel at first order. The actual finite
    // volume still contracts quadratically, although both endpoint densities round equal.
    let (card, caloric, density) = fixture(
        &CP,
        &[
            &[(290.0, -0.01), (310.0, 0.01)],
            &[(290.0, 0.01), (310.0, -0.01)],
            &[(290.0, 0.0), (310.0, 0.0)],
        ],
        true,
        300.0,
    );
    let step = IsobaricFreeExpansion::try_new(&card, &caloric, &density)
        .unwrap()
        .advance_heat(300.0, 2.0, heat)
        .unwrap();
    let expected = -0.002 * (0.001 * dt).powi(2);
    close(step.volume_change_m3, expected, expected.abs() * 2e-12);
    assert!(step.boundary_work_j < 0.0);
}

#[test]
fn incompatible_or_pressure_unqualified_sources_cannot_supply_expansion_work() {
    let strains = &[&[(290.0, 0.0), (310.0, 0.02)][..]][..];
    let (card, caloric, density) = fixture(&CP, strains, true, 300.0);
    let (_, narrowed) = prepare(&card, 1, 295.0, 310.0, 300.0, PRESSURE);
    assert!(IsobaricFreeExpansion::try_new(&card, &caloric, &narrowed).is_err());
    let (_, other_pressure) = prepare(&card, 1, 290.0, 310.0, 300.0, 2.0 * PRESSURE);
    assert!(IsobaricFreeExpansion::try_new(&card, &caloric, &other_pressure).is_err());
    let (unqualified, unqualified_cp, unqualified_density) = fixture(&CP, strains, false, 300.0);
    assert!(
        IsobaricFreeExpansion::try_new(&unqualified, &unqualified_cp, &unqualified_density)
            .is_err()
    );
    assert!(IsobaricFreeExpansion::try_new(&unqualified, &caloric, &density).is_err());
    let law = IsobaricFreeExpansion::try_new(&card, &caloric, &density).unwrap();
    assert!(law.advance_heat(300.0, 0.0, 1.0).is_err());
    assert!(law.advance_heat(300.0, 2.0, 1e100).is_err());
    assert_eq!(
        law.advance_heat(300.0, 2.0, 0.0).unwrap().state,
        density.at_temperature(300.0).unwrap()
    );
}
