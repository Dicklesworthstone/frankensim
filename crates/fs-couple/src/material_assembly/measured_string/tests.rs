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

/// A row beside the five roles: name, value, dimensions and stated band.
type ExtraRow<'a> = (&'a str, f64, Dims, Option<(f64, f64)>);

fn source(
    values: [f64; 5],
    band: Option<(f64, f64)>,
    wrong_dimensions: bool,
) -> ResolvedMaterialStatePoint {
    source_with(values, band, wrong_dimensions, &[])
}

fn source_with(
    values: [f64; 5],
    band: Option<(f64, f64)>,
    wrong_dimensions: bool,
    extra: &[ExtraRow<'_>],
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let roles = selectors();
    let mut rows: Vec<ExtraRow<'_>> = roles
        .roles()
        .into_iter()
        .zip(values)
        .map(|((name, dims), value)| (name, dims, value))
        .map(|(name, dims, value)| (name, value, dims, band))
        .collect();
    rows.extend_from_slice(extra);
    // The query sits inside every stated band, so the narrowest one sets it.
    let query_omega = rows
        .iter()
        .filter_map(|row| row.3)
        .map(|(low, high)| 0.5 * (low + high))
        .reduce(f64::min);
    for (name, value, mut dims, band) in rows {
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
    if let Some(omega) = query_omega {
        query = query.with("omega", omega).unwrap();
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

/// Woodhouse 2004 Table I, E2, with the EJ46 diameter stand-in: L, d, T, mu, EI.
const E2: [f64; 5] = [0.65, 0.001118, 71.6, 0.00624, 57e-6];
const INVERSE_SECOND: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// (mode, omega, zeta) of Woodhouse's eq. (8), p. 956, friction and bending terms
/// only, as zeta = eta / 2, on the 0.68 m span. Worked from the printed form in
/// 50-digit decimal by a separate script, not by this crate.
const E2_LAW: [(usize, f64, f64); 5] = [
    (1, 4.948_902_138_948_980_65e2, 1.016_974_690_405_744_40e-5),
    (12, 5.945_893_126_065_888_00e3, 3.438_430_494_314_561_72e-5),
    (24, 1.193_524_631_775_046_69e4, 1.068_281_843_449_075_90e-4),
    (48, 2.421_505_335_623_435_26e4, 3.863_688_595_087_020_64e-4),
    (96, 5.109_373_231_935_070_41e4, 1.362_599_646_361_274_71e-3),
];

fn e2_source() -> ResolvedMaterialStatePoint {
    source_with(
        E2,
        None,
        false,
        &[
            ("eta_F_e2", 2e-5, Dims::NONE, None),
            ("eta_B_e2", 0.02, Dims::NONE, None),
            ("eta_F_d3", 5e-5, Dims::NONE, None),
            ("eta_A_e2", 1.2, INVERSE_SECOND, None),
        ],
    )
}

fn loss(friction: &str, bending: &str) -> MeasuredStringLossProperties {
    MeasuredStringLossProperties {
        friction: friction.into(),
        bending: bending.into(),
    }
}

fn e2_on_span() -> ResolvedMeasuredStringSpecimen {
    ResolvedMeasuredStringSpecimen::try_new(96, &e2_source(), selectors())
        .unwrap()
        .on_span_from_speaking_length(0.68, 0.65 / 0.68)
        .unwrap()
}

fn zeta(string: &PrestressedString, mode: usize) -> f64 {
    let omega = fs_nlmodal::prestressed_beam_omega(
        string.length_m,
        string.tension_n,
        string.lin_density_kg_m,
        string.bending_stiffness_n_m2,
        mode,
    );
    let wave_number = mode as f64 * core::f64::consts::PI / string.length_m;
    crate::acoustic_realize::intrinsic_string_mode_zeta(string, omega, wave_number).unwrap()
}

#[test]
fn woodhouse_e2_rows_give_his_printed_friction_and_bending_terms() {
    let plain = e2_on_span();
    let bound = plain.clone().with_loss_factors(loss("eta_F_e2", "eta_B_e2")).unwrap();
    let string = bound.string();
    assert_eq!(string.tension_loss_factor, Some(2e-5));
    assert_eq!(string.bending_loss_factor, Some(0.02));
    assert_eq!(bound.loss_properties(), Some(&loss("eta_F_e2", "eta_B_e2")));
    assert_eq!(bound.measurements().identity(), plain.measurements().identity());
    for (mode, omega, expected) in E2_LAW {
        let actual = fs_nlmodal::prestressed_beam_omega(
            string.length_m,
            string.tension_n,
            string.lin_density_kg_m,
            string.bending_stiffness_n_m2,
            mode,
        );
        assert!((actual / omega - 1.0).abs() < 1e-14, "mode {mode}: omega {actual:e}");
        let rate = zeta(&string, mode);
        assert!((rate / expected - 1.0).abs() < 1e-13, "mode {mode}: zeta {rate:e} vs {expected:e}");
    }
    // Only the loss moved: no other law arrived and no coefficient was rescaled.
    let same = PrestressedString {
        tension_loss_factor: None,
        bending_loss_factor: None,
        ..string.clone()
    };
    assert_eq!(same, plain.string());
    assert!(string.kelvin_voigt_bending.is_none() && string.relaxation_bending.is_none());
    assert_eq!(string.damping_ratio, 0.0);
    // Binding the span after the loss gives the same string, bit for bit.
    let reordered = ResolvedMeasuredStringSpecimen::try_new(96, &e2_source(), selectors())
        .unwrap()
        .with_loss_factors(loss("eta_F_e2", "eta_B_e2"))
        .unwrap()
        .on_span_from_speaking_length(0.68, 0.65 / 0.68)
        .unwrap();
    assert_eq!(reordered.string(), string);
}

#[test]
fn no_selection_or_the_wrong_rows_miss_the_printed_law() {
    // Planted negatives. Each wrong selection lands on its own worked value, and the
    // comparison that passes above would fail on it by the stated factor.
    let unselected = e2_on_span();
    assert!(unselected.loss_properties().is_none());
    for (mode, _, _) in E2_LAW {
        assert_eq!(zeta(&unselected.string(), mode), 0.0);
    }
    let first = E2_LAW[0].2;
    let swapped = e2_on_span().with_loss_factors(loss("eta_B_e2", "eta_F_e2")).unwrap();
    let rate = zeta(&swapped.string(), 1);
    assert!((rate / 9.999_830_253_095_942_56e-3 - 1.0).abs() < 1e-13, "{rate:e}");
    assert!(rate / first > 980.0, "swapped rows only moved {}", rate / first);
    let other_string = e2_on_span().with_loss_factors(loss("eta_F_d3", "eta_B_e2")).unwrap();
    let rate = zeta(&other_string.string(), 1);
    assert!((rate / 2.516_949_202_882_612_65e-5 - 1.0).abs() < 1e-13, "{rate:e}");
    assert!(rate / first > 2.4, "the D3 row only moved {}", rate / first);
}

#[test]
fn loss_rows_refuse_wrong_quantities_missing_rows_repeats_and_nonphysical_values() {
    let plain = e2_on_span;
    let refused = |what: &str, selection: MeasuredStringLossProperties| {
        match plain().with_loss_factors(selection) {
            Err(AcousticRealizeError::InvalidDescription { what: actual }) => assert_eq!(actual, what),
            other => panic!("wanted {what}, got {:?}", other.map(|_| "a bound specimen")),
        }
    };
    // eta_A's s^-1 row cannot enter either dimensionless role.
    let scalars = "measured string needs positive dimension-only SI scalars with validity-wide constant support";
    refused(scalars, loss("eta_A_e2", "eta_B_e2"));
    refused(scalars, loss("eta_F_e2", "eta_A_e2"));
    let missing = "measured string is missing an explicitly selected source property";
    refused(missing, loss("mean_Q", "eta_B_e2"));
    refused(missing, loss("eta_F_e2", "mean_Q"));
    refused("invalid measured string loss selectors", loss("eta_F_e2", "eta_F_e2"));
    assert!(loss("eta_F_e2", "eta_F_e2").requirements().is_err());
    assert_eq!(
        loss("eta_F_e2", "eta_B_e2")
            .requirements()
            .unwrap()
            .iter()
            .map(|requirement| (requirement.name().to_owned(), requirement.dims()))
            .collect::<Vec<_>>(),
        [("eta_F_e2".to_owned(), Dims::NONE), ("eta_B_e2".to_owned(), Dims::NONE)]
    );
    let once = plain().with_loss_factors(loss("eta_F_e2", "eta_B_e2")).unwrap();
    assert!(once.with_loss_factors(loss("eta_F_d3", "eta_B_e2")).is_err());
    for bad in [0.0, -2e-5] {
        let rows = [
            ("eta_F_e2", bad, Dims::NONE, None),
            ("eta_B_e2", 0.02, Dims::NONE, None),
        ];
        let state = source_with(E2, None, false, &rows);
        let specimen = ResolvedMeasuredStringSpecimen::try_new(8, &state, selectors()).unwrap();
        assert!(specimen.with_loss_factors(loss("eta_F_e2", "eta_B_e2")).is_err(), "accepted {bad}");
    }
}

#[test]
fn a_loss_row_with_a_stated_band_narrows_the_specimen_and_readmits_every_mode() {
    let f1 = fs_nlmodal::prestressed_beam_omega(E2[0], E2[2], E2[3], E2[4], 1);
    let wide = Some((0.0, 100.0 * f1));
    let rows = [
        ("eta_F_e2", 2e-5, Dims::NONE, Some((0.0, 3.5 * f1))),
        ("eta_B_e2", 0.02, Dims::NONE, Some((0.0, 1.5 * f1))),
    ];
    let state = source_with(E2, wide, false, &rows);
    let selection = loss("eta_F_e2", "eta_B_e2");
    let one = ResolvedMeasuredStringSpecimen::try_new(1, &state, selectors()).unwrap();
    assert_eq!(one.omega_band_rad_s(), wide);
    let bound = one.with_loss_factors(selection.clone()).unwrap();
    assert_eq!(bound.omega_band_rad_s(), Some((0.0, 1.5 * f1)));
    // The second mode fits the five roles' band but not the loss rows' own.
    let two = ResolvedMeasuredStringSpecimen::try_new(2, &state, selectors()).unwrap();
    assert!(two.with_loss_factors(selection).is_err());
}
