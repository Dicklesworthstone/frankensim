//! Synthetic source/native-history checks; no measured material qualification.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::{
    ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA, J2HardnessError,
    ResolvedJ2Hardness,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY, ResolvedJ2Plasticity,
    resolve_j2_plasticity,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec, Temperature};

const GATES: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-9,
    yield_pa: 1e-9,
};
const TEMPERATURE: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

fn fixture(
    value: PropertyValue,
    policy: InterpolationPolicy,
    hardening: f64,
) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: "synthetic hardness/plastic-history fixture, no measured data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let environment = ValidityDomain::unconstrained()
        .with_quantity("T", TEMPERATURE, 280.0, 380.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, hardening),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: environment.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic J2 specimen at declared accumulated plastic strain".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"hardness on J2 history"),
            caveats: "numerical regression only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic spherical indenter with radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic projected-area indentation protocol revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value,
            validity: environment.with_quantity(
                ALPHA,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.04,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: policy,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "synthetic plastic-history specimen".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (card, key)
}

fn curve(scale: f64) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: ALPHA.into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(0.0, scale * 100.0), (0.01, scale * 300.0), (0.04, scale * 150.0)],
        dims: Pressure::DIMS,
    }
}

fn plastic(card: &MaterialCard, temperature: f64) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &QueryPoint::new()
            .with_quantity("T", TEMPERATURE, temperature)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.05,
    )
    .unwrap()
}

fn prepare(
    card: &MaterialCard,
    key: &PropertyKey,
    temperature: f64,
    maximum_alpha: f64,
) -> ResolvedJ2Hardness {
    ResolvedJ2Hardness::new(
        plastic(card, temperature),
        card,
        key.clone(),
        MaterialPropertySelection::SingleClaimOnly,
        maximum_alpha,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_native_plastic_work_and_independent_hardness_share_one_accepted_history() {
    // Perfectly plastic mechanical response can still have a separately
    // measured hardness increase. No sigma_y-to-hardness factor is imposed.
    let (card, key) = fixture(curve(1.0), InterpolationPolicy::LinearInside, 0.0);
    let model = prepare(&card, &key, 300.0, 0.04);
    let mut state = model.initial_state();
    for axial in [0.006, 0.0058, -0.006, 0.006] {
        let strain = [axial, 0.0, 0.0, 0.0, 0.0, 0.0];
        let original = model
            .plastic()
            .law()
            .advance_with_energy(state.strain(), state.plastic(), &strain, GATES)
            .unwrap();
        let step = model.advance_candidate(&state, &strain, GATES).unwrap();
        assert_eq!(step.stress_pa, original.stress_pa);
        assert_eq!(step.state.plastic(), &original.state);
        assert_eq!(step.energy, original.energy);
        let exported = model.resolve_hardness_state(&card, &step.state).unwrap();
        near(step.state.hardness_pa(), exported.hardness_pa());
        assert_eq!(exported.key(), &key);
        assert_eq!(
            exported.property().answer().receipt.selected,
            model.hardness_source().property().answer().receipt.selected
        );
        assert_eq!(exported.resolved().state_coordinate("T"), Some(300.0));
        assert_eq!(
            exported.resolved().state_coordinate(ALPHA),
            Some(step.state.plastic().alpha)
        );
        state = step.state;
    }
    assert!(state.plastic().alpha > 0.01);
    assert_eq!(model.plastic().hardening_modulus_pa(), 0.0);
    assert_eq!(model.plastic().law().yield_stress, 1.0);
    assert!(state.hardness_pa() > 100.0);
    near(model.hardness_pa_at(0.025).unwrap(), 225.0);
    assert!(model.hardness_pa_at(0.04).unwrap() < model.hardness_pa_at(0.01).unwrap());
}

#[test]
fn g3_elastic_unloading_keeps_hardness_and_reversal_uses_accumulated_strain() {
    let (card, key) = fixture(curve(1.0), InterpolationPolicy::LinearInside, 100.0);
    let model = prepare(&card, &key, 300.0, 0.04);
    let initial = model.initial_state();
    let elastic = model
        .advance_candidate(&initial, &[1e-4, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
        .unwrap();
    assert_eq!(elastic.state.plastic().alpha, 0.0);
    assert_eq!(elastic.state.hardness_pa(), 100.0);
    let loaded = model
        .advance_candidate(&elastic.state, &[0.006, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
        .unwrap();
    assert!(loaded.state.plastic().alpha > 0.0);
    near(
        loaded.state.hardness_pa(),
        100.0 + 20000.0 * loaded.state.plastic().alpha,
    );
    let unloaded = model
        .advance_candidate(&loaded.state, &[0.0058, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
        .unwrap();
    assert_eq!(unloaded.state.plastic(), loaded.state.plastic());
    assert_eq!(unloaded.hardness_change_pa, 0.0);
    let reversed = model
        .advance_candidate(&unloaded.state, &[-0.006, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
        .unwrap();
    assert!(reversed.state.plastic().alpha > unloaded.state.plastic().alpha);
    assert!(reversed.state.hardness_pa() > unloaded.state.hardness_pa());
    // The accumulated history does not cancel when the signed loading reverses.
    assert!(reversed.state.plastic().plastic_strain[0] < 0.0);
}

fn surface_with_gap(gap: bool) -> PropertyValue {
    let spans: &[(f64, f64)] = if gap {
        &[(0.0, 0.01), (0.03, 0.04)]
    } else {
        &[(0.0, 0.04)]
    };
    let mut samples = Vec::new();
    let mut triangles = Vec::new();
    for &(lo, hi) in spans {
        let base = samples.len() as u32;
        for (alpha, temperature) in [(lo, 280.0), (hi, 280.0), (hi, 380.0), (lo, 380.0)] {
            samples.push(([alpha, temperature], 400.0 + 10000.0 * alpha - temperature));
        }
        triangles.extend([[base, base + 1, base + 2], [base, base + 2, base + 3]]);
    }
    PropertyValue::Surface2D {
        axes: [ALPHA.into(), "T".into()],
        axis_dims: [Dims::NONE, Temperature::DIMS],
        samples,
        triangles,
        dims: Pressure::DIMS,
    }
}

#[test]
fn g3_temperature_cuts_keep_original_cells_and_refuse_supported_endpoints_across_a_hole() {
    let (card, key) = fixture(
        surface_with_gap(false),
        InterpolationPolicy::LinearInside,
        100.0,
    );
    let cold = prepare(&card, &key, 300.0, 0.04);
    let warm = prepare(&card, &key, 350.0, 0.04);
    for alpha in [0.0, 0.004, 0.02, 0.03, 0.04] {
        near(cold.hardness_pa_at(alpha).unwrap(), 100.0 + 10000.0 * alpha);
        near(warm.hardness_pa_at(alpha).unwrap(), 50.0 + 10000.0 * alpha);
    }
    assert!(cold.source_line().knots.len() >= 3);
    assert!(
        warm.advance_candidate(
            &cold.initial_state(),
            &[0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
            GATES,
        )
        .is_err()
    );
    let (gapped, key) = fixture(
        surface_with_gap(true),
        InterpolationPolicy::LinearInside,
        100.0,
    );
    let result = ResolvedJ2Hardness::new(
        plastic(&gapped, 330.0),
        &gapped,
        key,
        MaterialPropertySelection::SingleClaimOnly,
        0.04,
    );
    assert!(matches!(result, Err(J2HardnessError::Support(_))));
}

#[test]
fn g0_source_range_and_native_refusals_leave_the_entire_history_unchanged() {
    let (card, key) = fixture(curve(1.0), InterpolationPolicy::LinearInside, 100.0);
    let model = prepare(&card, &key, 300.0, 0.001);
    let state = model.initial_state();
    let before = state.clone();
    for axial in [0.01, 0.1, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .advance_candidate(&state, &[axial, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
                .is_err()
        );
        assert_eq!(state, before);
    }
    for alpha in [-1.0, 0.0011, f64::NAN, f64::INFINITY] {
        assert!(model.hardness_pa_at(alpha).is_err());
    }
    let another = prepare(&card, &key, 300.0, 0.002);
    assert!(another.advance_candidate(&state, &[0.0; 6], GATES).is_err());
    let replacement = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    assert!(model.resolve_hardness_state(&replacement, &state).is_err());
    assert!(
        ResolvedJ2Hardness::new(
            plastic(&card, 300.0),
            &replacement,
            key,
            MaterialPropertySelection::SingleClaimOnly,
            0.01,
        )
        .is_err()
    );
    assert_eq!(state, before);
}

#[test]
fn g0_discrete_only_negative_interior_and_generic_strain_sources_cannot_evolve() {
    for (value, policy) in [
        (curve(1.0), InterpolationPolicy::TabulatedOnly),
        (
            PropertyValue::Curve {
                abscissa: ALPHA.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 100.0), (0.02, -1.0), (0.04, 100.0)],
                dims: Pressure::DIMS,
            },
            InterpolationPolicy::LinearInside,
        ),
        (
            PropertyValue::Curve {
                abscissa: "total-strain".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 100.0), (0.04, 200.0)],
                dims: Pressure::DIMS,
            },
            InterpolationPolicy::LinearInside,
        ),
    ] {
        let (card, key) = fixture(value, policy, 100.0);
        assert!(
            ResolvedJ2Hardness::new(
                plastic(&card, 300.0),
                &card,
                key,
                MaterialPropertySelection::SingleClaimOnly,
                0.04,
            )
            .is_err()
        );
    }
}
