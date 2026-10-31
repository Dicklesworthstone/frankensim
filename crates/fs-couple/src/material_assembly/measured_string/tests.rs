use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};

fn selectors() -> MeasuredStringProperties {
    MeasuredStringProperties {
        length: "L".into(),
        diameter: "d".into(),
        tension: "Tension".into(),
        linear_density: "mu".into(),
        bending_stiffness: "EI".into(),
    }
}

fn source(
    values: [f64; 5],
    band: Option<(f64, f64)>,
    wrong_dimensions: bool,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for ((name, mut dims), value) in selectors().roles().into_iter().zip(values) {
        if wrong_dimensions && name == "mu" {
            dims = fs_qty::Density::DIMS;
        }
        let mut validity = ValidityDomain::unconstrained().with("temperature_known", 0.0, 0.0);
        if let Some((low, high)) = band {
            validity = validity.with("omega", low, high);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.68,
                },
                provenance: Provenance {
                    source: "measured-section admission test".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
        // Deliberately permissive upstream admission tests the adapter itself.
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test whole-string specimen".into(),
            phase: "unresolved composite".into(),
            process: "frozen measured setup".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut query = QueryPoint::new().with("temperature_known", 0.0).unwrap();
    if let Some((low, high)) = band {
        query = query.with("omega", 0.5 * (low + high)).unwrap();
    }
    resolve_material_state_point(
        &card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn measured_cello_sections_keep_their_distinct_dispersion_without_inventing_bulk_moduli() {
    // These published mu/EI values are inferred effective coefficients. This
    // test preserves their stated continuum dispersion without adding fluid
    // inertia or claiming an independently weighed dry structural mass.
    // Lampis et al. JASA Express Letters 4:113201 (2024), Table 1;
    // values also reside in data/matdb/seed-v1/cello-g-strings-whole-lampis-2024.
    // This checks the transverse constitutive calculation, not cello fidelity.
    let specimens = [
        [0.70, 0.001092, 149.6, 0.00794, 0.000685],
        [0.70, 0.000877, 153.2, 0.00814, 0.000432],
        [0.70, 0.001150, 116.3, 0.00620, 0.000303],
        [0.70, 0.000895, 153.5, 0.00820, 0.000310],
    ];
    let mut ratios = Vec::new();
    for values in specimens {
        let state = source(values, None, false);
        let measured = ResolvedMeasuredStringSpecimen::try_new(30, &state, selectors()).unwrap();
        let string = measured.string();
        assert_eq!(string.length_m, values[0]);
        assert_eq!(2.0 * measured.radius_m(), values[1]);
        assert_eq!(string.tension_n, values[2]);
        assert_eq!(string.lin_density_kg_m, values[3]);
        assert_eq!(string.bending_stiffness_n_m2, values[4]);
        assert_eq!(string.axial_stiffness_n, 0.0);
        assert_eq!(string.damping_ratio, 0.0);
        assert!(string.kelvin_voigt_bending.is_none() && string.relaxation_bending.is_none());
        assert!(string.tension_loss_factor.is_none());
        assert!(measured.measurements().property("young_modulus").is_none());
        assert!(measured.measurements().property("density").is_none());
        assert!(measured.omega_band_rad_s().is_none());
        for mode in [1, 7, 30] {
            let k = (mode as f64) * core::f64::consts::PI / values[0];
            let expected = ((values[2] * k * k + values[4] * k * k * k * k) / values[3]).sqrt();
            let actual = fs_nlmodal::prestressed_beam_omega(
                values[0], values[2], values[3], values[4], mode,
            );
            assert!((actual / expected - 1.0).abs() < 2e-15);
        }
        let frequency = |mode| {
            fs_nlmodal::prestressed_beam_omega(
                string.length_m,
                string.tension_n,
                string.lin_density_kg_m,
                string.bending_stiffness_n_m2,
                mode,
            )
        };
        ratios.push(frequency(30) / (30.0 * frequency(1)));
    }
    assert!(ratios[0] > ratios[1] && ratios[1] > ratios[2] && ratios[2] > ratios[3]);
    assert!(
        ratios.iter().all(|ratio| *ratio > 1.0),
        "measured EI must round the stiff-string dispersion"
    );
}

#[test]
fn original_source_state_receipts_and_uncertainty_survive_exactly() {
    let state = source([0.7, 0.00115, 116.3, 0.0062, 0.000303], None, false);
    let measured = ResolvedMeasuredStringSpecimen::try_new(3, &state, selectors()).unwrap();
    assert_eq!(measured.measurements().identity(), state.identity());
    assert_eq!(
        measured.measurements().card_identity(),
        state.card_identity()
    );
    assert_eq!(
        measured.measurements().query_point(),
        &[("temperature_known".into(), 0.0)]
    );
    assert_eq!(measured.properties(), &selectors());
    for property in state.properties() {
        let preserved = measured
            .measurements()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            format!("{:?}", preserved.answer()),
            format!("{:?}", property.answer())
        );
    }
    // Material substitution changes the actual source values, with fixed T/L.
    let mut altered = [0.7, 0.00115, 116.3, 0.0062, 0.000303];
    altered[3] *= 4.0;
    let other =
        ResolvedMeasuredStringSpecimen::try_new(3, &source(altered, None, false), selectors())
            .unwrap();
    assert_ne!(
        other.measurements().identity(),
        measured.measurements().identity()
    );
    let a = measured.string();
    let b = other.string();
    let omega = |s: &PrestressedString| {
        fs_nlmodal::prestressed_beam_omega(
            s.length_m,
            s.tension_n,
            s.lin_density_kg_m,
            s.bending_stiffness_n_m2,
            3,
        )
    };
    assert!((omega(&b) / omega(&a) - 0.5).abs() < 2e-15);
}

#[test]
fn actual_source_frequency_domain_and_numerical_limits_are_enforced() {
    let values = [0.7, 0.00115, 116.3, 0.0062, 0.000303];
    let f1 = fs_nlmodal::prestressed_beam_omega(values[0], values[2], values[3], values[4], 1);
    let state = source(values, Some((0.0, 1.1 * f1)), false);
    let measured = ResolvedMeasuredStringSpecimen::try_new(1, &state, selectors()).unwrap();
    assert_eq!(measured.omega_band_rad_s(), Some((0.0, 1.1 * f1)));
    assert!(ResolvedMeasuredStringSpecimen::try_new(2, &state, selectors()).is_err());
    assert!(measured.admit_angular_frequency(1.2 * f1).is_err());
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(measured.admit_angular_frequency(omega).is_err());
    }
    for count in [0, fs_modal::time::MAX_TIME_DOMAIN_ACOUSTIC_MODES + 1] {
        assert!(ResolvedMeasuredStringSpecimen::try_new(count, &state, selectors()).is_err());
    }
}

#[test]
fn a_speaking_length_binds_nut_to_support_and_no_other_span_or_station() {
    // A guitar-like source measured over its 0.65 m scale (nut to saddle), on a 0.68 m pinned span.
    let values = [0.65, 0.0011, 71.6, 0.00624, 57e-6];
    let state = source(values, None, false);
    let measured = ResolvedMeasuredStringSpecimen::try_new(24, &state, selectors()).unwrap();
    assert_eq!(measured.speaking_length_m(), None);
    let station = 0.65 / 0.68;
    let bound = measured
        .clone()
        .on_span_from_speaking_length(0.68, station)
        .unwrap();
    let (original, string) = (measured.string(), bound.string());
    assert_eq!(string.length_m, 0.68);
    assert_eq!(bound.speaking_length_m(), Some(0.65));
    assert_eq!(string.tension_n, original.tension_n);
    assert_eq!(string.lin_density_kg_m, original.lin_density_kg_m);
    assert_eq!(string.bending_stiffness_n_m2, original.bending_stiffness_n_m2);
    assert_eq!(bound.radius_m(), measured.radius_m());
    assert_eq!(bound.measurements().identity(), measured.measurements().identity());
    // The measurement is never rescaled: a span or station that misses it by more than rounding refuses,
    // including the six-digit station 0.955882 (0.68 x 0.955882 = 0.64999976 m).
    for (span, fraction) in [
        (0.68, 0.955_882),
        (0.70, station),
        (0.65, station),
        (0.68, 1.0),
        (0.68, 0.0),
        (0.68, f64::NAN),
        (-0.68, station),
        (f64::INFINITY, station),
    ] {
        assert!(
            measured
                .clone()
                .on_span_from_speaking_length(span, fraction)
                .is_err(),
            "accepted span {span} station {fraction}"
        );
    }
    assert!(bound.on_span_from_speaking_length(0.68, station).is_err());
    // The longer span's modes are re-admitted against the source band: a band that holds the 0.65 m
    // fundamental but not the 4.4 percent lower 0.68 m one refuses the binding.
    let f1 = fs_nlmodal::prestressed_beam_omega(values[0], values[2], values[3], values[4], 1);
    let banded = source(values, Some((0.97 * f1, 1.1 * f1)), false);
    let narrow = ResolvedMeasuredStringSpecimen::try_new(1, &banded, selectors()).unwrap();
    assert!(narrow.on_span_from_speaking_length(0.68, station).is_err());
}

#[test]
fn missing_duplicate_wrong_quantity_and_nonphysical_measurements_refuse() {
    let values = [0.7, 0.00115, 116.3, 0.0062, 0.000303];
    let state = source(values, None, false);
    let mut missing = selectors();
    missing.bending_stiffness = "mean_Q".into();
    assert!(ResolvedMeasuredStringSpecimen::try_new(1, &state, missing).is_err());
    let mut duplicate = selectors();
    duplicate.length = duplicate.diameter.clone();
    assert!(duplicate.requirements().is_err());
    assert!(
        ResolvedMeasuredStringSpecimen::try_new(1, &source(values, None, true), selectors())
            .is_err()
    );
    for index in 0..values.len() {
        let mut zero = values;
        zero[index] = 0.0;
        assert!(
            ResolvedMeasuredStringSpecimen::try_new(1, &source(zero, None, false), selectors())
                .is_err()
        );
        let mut negative = values;
        negative[index] = -values[index];
        assert!(
            ResolvedMeasuredStringSpecimen::try_new(1, &source(negative, None, false), selectors())
                .is_err()
        );
    }
    let mut tiny = values;
    tiny[1] = f64::from_bits(1);
    assert!(
        ResolvedMeasuredStringSpecimen::try_new(1, &source(tiny, None, false), selectors())
            .is_err()
    );
}
