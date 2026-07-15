//! Conformance suite for fs-qty (plan §13.3): any reimplementation must pass
//! these cases. Each case logs a JSON-lines verdict so failures are
//! diagnosable from output alone (observability standard; fs-obs schema
//! adoption pending that crate's landing).

use fs_qty::parse::{ParseErrorKind, parse_qty};
use fs_qty::semantic::{
    AngleDomain, CompositionBasis, FormRequirement, PhasorAmplitude, PhasorQty, QuantityKind,
    SemanticError, SemanticQty, SemanticType, StrainBasis, StrainComponent, ValueForm,
    amount_concentration_to_mass_concentration, amount_to_mass,
    mass_concentration_to_amount_concentration, mass_to_amount,
};
use fs_qty::{Dims, DynViscosity, Length, Pressure, QtyAny, Time};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-qty/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

/// qty-001: the FrankenScript literal battery — exact values and dimensions
/// for every literal form the plan's example studies use (Appendix C).
#[test]
fn qty_001_appendix_c_literal_battery() {
    let cases: &[(&str, f64, [i8; 6])] = &[
        ("0.12Pa*s", 0.12, [-1, 1, -1, 0, 0, 0]),
        ("0.061N/m", 0.061, [0, 1, -2, 0, 0, 0]),
        ("0.5L/s", 5e-4, [3, 0, -1, 0, 0, 0]),
        ("12mm", 0.012, [1, 0, 0, 0, 0, 0]),
        (
            "65deg",
            65.0 * core::f64::consts::PI / 180.0,
            [0, 0, 0, 0, 0, 0],
        ),
        ("2h", 7200.0, [0, 0, 1, 0, 0, 0]),
        ("0.03m2/s3", 0.03, [2, 0, -3, 0, 0, 0]),
        ("15rad/s", 15.0, [0, 0, -1, 0, 0, 0]),
        ("2e-2", 0.02, [0, 0, 0, 0, 0, 0]),
        ("1mol", 1.0, [0, 0, 0, 0, 0, 1]),
    ];
    for (text, value, dims) in cases {
        let q = parse_qty(text).unwrap_or_else(|e| panic!("{text}: {e}"));
        let ok = (q.value - value).abs() <= 1e-12 * value.abs().max(1.0) && q.dims == Dims(*dims);
        verdict(
            &format!("qty-001/{text}"),
            ok,
            &format!("value={} dims={:?}", q.value, q.dims),
        );
    }
}

/// qty-002: typed and erased algebra agree bit-for-bit.
#[test]
fn qty_002_typed_erased_agreement() {
    let typed = (Length::new(0.37) / Time::new(1.61)).value();
    let erased = (Length::new(0.37).erase() / Time::new(1.61).erase()).value;
    verdict(
        "qty-002/bit-agreement",
        typed.to_bits() == erased.to_bits(),
        &format!("typed={typed:?} erased={erased:?}"),
    );
}

/// qty-003: JSON round-trip is bit-exact and shape-canonical.
#[test]
fn qty_003_json_round_trip() {
    let q = DynViscosity::new(0.12).erase();
    let text = fs_qty::json::to_json(q).expect("finite");
    let canonical = text == r#"{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}"#;
    let back = fs_qty::json::from_json(&text).expect("parses");
    verdict(
        "qty-003/round-trip",
        canonical && back.value.to_bits() == q.value.to_bits() && back.dims == q.dims,
        &text,
    );
}

/// qty-003b: legacy bytes survive unchanged and decode only with an immutable
/// five-to-six semantic-crosswalk receipt.
#[test]
fn qty_003b_legacy_json_crosswalk() {
    const OLD: &str = r#"{"value":0.12,"dims":[-1,1,-1,0,0]}"#;
    let decoded = fs_qty::json::decode_json(OLD).expect("legacy decode");
    let receipt = decoded.migration().expect("receipt required");
    let new = fs_qty::json::to_json(decoded.qty()).expect("canonical v2");
    verdict(
        "qty-003b/legacy-crosswalk",
        decoded.qty().dims == Dims([-1, 1, -1, 0, 0, 0])
            && fs_qty::json::to_legacy_json(decoded.qty()).as_deref() == Ok(OLD)
            && receipt.verifies(OLD.as_bytes(), new.as_bytes()),
        &format!(
            "source={:?} target={:?} rule={:?} old_hash={} new_hash={}",
            receipt.source_version(),
            receipt.target_version(),
            receipt.rule(),
            receipt.old_hash(),
            receipt.new_hash()
        ),
    );
}

/// qty-004: dimension safety — runtime mismatches produce structured,
/// teaching errors; downcasts check dimensions.
#[test]
fn qty_004_dimension_safety() {
    let e = Pressure::new(1.0)
        .erase()
        .try_add(Time::new(1.0).erase())
        .unwrap_err();
    verdict(
        "qty-004/mismatch-error",
        e.to_string().contains("dimension mismatch"),
        &e.to_string(),
    );
    let bad: Result<Pressure, _> = QtyAny::dimensionless(1.0).to_typed();
    verdict(
        "qty-004/downcast-checked",
        bad.is_err(),
        "dimensionless -> Pressure must fail",
    );
}

/// qty-005: the parser never panics (a compressed re-run of the hardening
/// battery at conformance level — reimplementations must hold this too).
#[test]
fn qty_005_parser_total_over_garbage() {
    let mut state: u64 = 0x00C0_FFEE;
    for _ in 0..2_000 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bytes: Vec<u8> = (0..(state % 16))
            .map(|i| ((state >> (i % 56)) & 0x7F) as u8)
            .collect();
        let s = String::from_utf8_lossy(&bytes);
        let _ = parse_qty(&s);
    }
    verdict("qty-005/no-panic", true, "2000 garbage inputs, no panic");
}

fn logged_parser_refusal(case: &str, input: &str, expected: &ParseErrorKind) {
    let outcome = std::panic::catch_unwind(|| parse_qty(input));
    match outcome {
        Ok(Err(error)) => verdict(
            case,
            error.input == input
                && &error.kind == expected
                && error.at <= input.len()
                && !error.help.is_empty(),
            &format!(
                "input={input}; expected={expected:?}; actual={:?}; at={}; help={}",
                error.kind, error.at, error.help
            ),
        ),
        Ok(Ok(quantity)) => verdict(
            case,
            false,
            &format!("input={input}; unexpectedly admitted={quantity:?}"),
        ),
        Err(_) => verdict(case, false, &format!("input={input}; parser panicked")),
    }
}

/// qty-006: exponent and floating-point authority boundaries refuse before
/// narrowing or retention. Each case logs its exact structured error.
#[test]
fn qty_006_exponent_and_nonfinite_boundaries_are_total() {
    let positive_cap = parse_qty("1m^60").expect("positive cap");
    let negative_cap = parse_qty("1m^-60").expect("negative cap");
    let derived_cap = parse_qty("1Pa^30").expect("derived cap");
    verdict(
        "qty-006/exact-cap",
        positive_cap.value.is_finite()
            && negative_cap.value.is_finite()
            && derived_cap.value.is_finite()
            && positive_cap.dims == Dims([60, 0, 0, 0, 0, 0])
            && negative_cap.dims == Dims([-60, 0, 0, 0, 0, 0])
            && derived_cap.dims == Dims([-30, 30, -60, 0, 0, 0]),
        &format!(
            "positive={:?}; negative={:?}; derived={:?}",
            positive_cap.dims, negative_cap.dims, derived_cap.dims
        ),
    );
    for input in ["0", "-0.0", "0e999", "-0e999"] {
        let zero = parse_qty(input).unwrap_or_else(|error| panic!("{input}: {error}"));
        verdict(
            &format!("qty-006/textual-zero/{input}"),
            zero.value == 0.0 && zero.dims.is_none(),
            &format!("input={input}; admitted={zero:?}"),
        );
    }

    for (case, input) in [
        ("qty-006/cap-plus-one", "1m^61"),
        ("qty-006/i8-min", "1m^-128"),
        ("qty-006/derived-over-cap", "1Pa^31"),
        ("qty-006/derived-hostile", "1Pa^127"),
    ] {
        logged_parser_refusal(case, input, &ParseErrorKind::BadExponent);
    }

    let repeated_positive = format!("1{}", vec!["m"; 61].join("*"));
    let repeated_negative = format!("1rad{}", "/m".repeat(61));
    logged_parser_refusal(
        "qty-006/repeated-positive",
        &repeated_positive,
        &ParseErrorKind::BadExponent,
    );
    logged_parser_refusal(
        "qty-006/repeated-negative",
        &repeated_negative,
        &ParseErrorKind::BadExponent,
    );

    for (case, input) in [
        ("qty-006/raw-overflow", "1e999"),
        ("qty-006/unit-raw-overflow", "1e999m"),
        ("qty-006/raw-positive-underflow", "1e-999"),
        ("qty-006/raw-negative-underflow", "-1e-999"),
        ("qty-006/prefix-product-overflow", "1e308TW"),
        ("qty-006/scale-power-overflow", "1Trad^26"),
        ("qty-006/scale-power-underflow", "1THz^-30"),
        ("qty-006/product-underflow", "1e-308pHz^2"),
        ("qty-006/division-underflow", "1e-308rad/THz^25"),
        ("qty-006/division-by-underflow", "1rad/THz^-30"),
    ] {
        logged_parser_refusal(case, input, &ParseErrorKind::NonFiniteValue);
    }
}

fn all_semantic_kinds() -> Vec<QuantityKind> {
    vec![
        QuantityKind::AbsoluteTemperature,
        QuantityKind::TemperatureDifference,
        QuantityKind::Angle(AngleDomain::Mechanical),
        QuantityKind::Angle(AngleDomain::Electrical),
        QuantityKind::AngularVelocity(AngleDomain::Mechanical),
        QuantityKind::AngularVelocity(AngleDomain::Electrical),
        QuantityKind::Torque,
        QuantityKind::Energy,
        QuantityKind::Pressure,
        QuantityKind::Stress,
        QuantityKind::Strain {
            basis: StrainBasis::Tensor,
            component: StrainComponent::Normal,
        },
        QuantityKind::Strain {
            basis: StrainBasis::Tensor,
            component: StrainComponent::Shear,
        },
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Normal,
        },
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Shear,
        },
        QuantityKind::Composition(CompositionBasis::MassFraction),
        QuantityKind::Composition(CompositionBasis::MoleFraction),
        QuantityKind::Composition(CompositionBasis::VolumeFraction),
        QuantityKind::Mass,
        QuantityKind::Amount,
        QuantityKind::MolarMass,
        QuantityKind::MassConcentration,
        QuantityKind::AmountConcentration,
        QuantityKind::Entropy,
        QuantityKind::HeatCapacity,
        QuantityKind::AcousticPressure,
        QuantityKind::AcousticPower,
    ]
}

const fn expected_waveform_kind(kind: QuantityKind) -> bool {
    match kind {
        QuantityKind::TemperatureDifference
        | QuantityKind::Angle(_)
        | QuantityKind::AngularVelocity(_)
        | QuantityKind::Torque
        | QuantityKind::Pressure
        | QuantityKind::Stress
        | QuantityKind::Strain { .. }
        | QuantityKind::AcousticPressure => true,
        QuantityKind::AbsoluteTemperature
        | QuantityKind::Energy
        | QuantityKind::Composition(_)
        | QuantityKind::Mass
        | QuantityKind::Amount
        | QuantityKind::MolarMass
        | QuantityKind::MassConcentration
        | QuantityKind::AmountConcentration
        | QuantityKind::Entropy
        | QuantityKind::HeatCapacity
        | QuantityKind::AcousticPower => false,
    }
}

/// qty-007: every parameterized kind crosses every scalar/phasor form through
/// the public carriers. Each decision is emitted as one deterministic JSONL
/// row, and failures accumulate so later matrix entries remain observable.
#[test]
#[allow(clippy::too_many_lines)] // The complete logged Cartesian policy matrix is one auditable G0 artifact.
fn qty_007_semantic_kind_form_matrix_is_sealed_and_logged() {
    let mut failures = Vec::new();
    for kind in all_semantic_kinds() {
        let waveform = expected_waveform_kind(kind);
        for form in [
            ValueForm::Static,
            ValueForm::Instantaneous,
            ValueForm::Peak,
            ValueForm::Rms,
        ] {
            let expected = form == ValueForm::Static || waveform;
            let semantic_type = SemanticType::new(kind, form);
            let outcome = SemanticQty::new(QtyAny::new(1.0, kind.expected_dims()), semantic_type);
            let actual = outcome.is_ok();
            let typed = match outcome {
                Ok(_) => expected,
                Err(SemanticError::UnsupportedForm {
                    source,
                    requirement: FormRequirement::StaticOnly,
                    ..
                }) => !expected && source == semantic_type,
                Err(_) => false,
            };
            let pass = actual == expected && typed && kind.admits_scalar_form(form) == expected;
            println!(
                "{{\"suite\":\"fs-qty/conformance\",\"case\":\"qty-007/scalar-form\",\"kind\":\"{kind:?}\",\"form\":\"{form:?}\",\"expected\":\"{}\",\"actual\":\"{}\",\"verdict\":\"{}\"}}",
                if expected { "admit" } else { "refuse" },
                if actual { "admit" } else { "refuse" },
                if pass { "pass" } else { "fail" },
            );
            if !pass {
                failures.push(format!("scalar {kind:?}/{form:?}"));
            }
        }

        for amplitude in [PhasorAmplitude::Peak, PhasorAmplitude::Rms] {
            let dims = kind.expected_dims();
            let outcome = PhasorQty::new(
                QtyAny::new(-1.0, dims),
                QtyAny::new(2.0, dims),
                kind,
                amplitude,
            );
            let actual = outcome.is_ok();
            let requested_form = match amplitude {
                PhasorAmplitude::Peak => ValueForm::Peak,
                PhasorAmplitude::Rms => ValueForm::Rms,
            };
            let typed = match outcome {
                Ok(phasor) => {
                    waveform
                        && phasor.real().value.to_bits() == (-1.0f64).to_bits()
                        && phasor.imaginary().value.to_bits() == 2.0f64.to_bits()
                }
                Err(SemanticError::UnsupportedForm {
                    source,
                    requirement: FormRequirement::StaticOnly,
                    ..
                }) => !waveform && source == SemanticType::new(kind, requested_form),
                Err(_) => false,
            };
            let pass = actual == waveform && typed && kind.admits_phasor() == waveform;
            println!(
                "{{\"suite\":\"fs-qty/conformance\",\"case\":\"qty-007/phasor-form\",\"kind\":\"{kind:?}\",\"form\":\"Phasor{amplitude:?}\",\"expected\":\"{}\",\"actual\":\"{}\",\"verdict\":\"{}\"}}",
                if waveform { "admit" } else { "refuse" },
                if actual { "admit" } else { "refuse" },
                if pass { "pass" } else { "fail" },
            );
            if !pass {
                failures.push(format!("phasor {kind:?}/{amplitude:?}"));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "semantic matrix failures: {failures:?}"
    );
}

fn semantic_from_literal(text: &str, kind: QuantityKind) -> SemanticQty {
    let quantity = parse_qty(text).unwrap_or_else(|error| panic!("{text}: {error}"));
    SemanticQty::new(quantity, SemanticType::new(kind, ValueForm::Static))
        .unwrap_or_else(|error| panic!("{text}: {error}"))
}

fn within_conversion_enclosure(actual: f64, expected: f64) -> bool {
    actual > 0.0
        && (actual - expected).abs()
            <= 16.0 * f64::EPSILON * actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE)
}

/// qty-008: equivalent unit spellings commute with mass/amount and
/// concentration-basis conversion through the public parser and semantic API.
#[test]
fn qty_008_mass_amount_unit_rescaling_is_metamorphic() {
    let cases = [
        ("18g", "18g/mol"),
        ("0.018kg", "0.018kg/mol"),
        ("18000mg", "18000mg/mol"),
    ];
    let mut amounts = Vec::new();
    for (mass_text, molar_mass_text) in cases {
        let mass = semantic_from_literal(mass_text, QuantityKind::Mass);
        let molar_mass = semantic_from_literal(molar_mass_text, QuantityKind::MolarMass);
        let amount = mass_to_amount(mass, molar_mass).expect("positive conversion");
        let recovered = amount_to_mass(amount, molar_mass).expect("inverse conversion");
        let pass = within_conversion_enclosure(amount.value(), 1.0)
            && within_conversion_enclosure(recovered.value(), mass.value());
        println!(
            "{{\"suite\":\"fs-qty/conformance\",\"case\":\"qty-008/unit-rescaling\",\"mass\":\"{mass_text}\",\"molar_mass\":\"{molar_mass_text}\",\"amount\":{},\"recovered_mass\":{},\"verdict\":\"{}\"}}",
            amount.value(),
            recovered.value(),
            if pass { "pass" } else { "fail" },
        );
        assert!(pass, "unit-rescaling conversion failed for {mass_text}");
        amounts.push(amount.value());
    }
    assert!(
        amounts
            .iter()
            .all(|&amount| within_conversion_enclosure(amount, amounts[0]))
    );

    let mass_concentration = semantic_from_literal("18g/m3", QuantityKind::MassConcentration);
    let molar_mass = semantic_from_literal("18g/mol", QuantityKind::MolarMass);
    let amount_concentration =
        mass_concentration_to_amount_concentration(mass_concentration, molar_mass)
            .expect("concentration basis conversion");
    let recovered = amount_concentration_to_mass_concentration(amount_concentration, molar_mass)
        .expect("concentration inverse");
    verdict(
        "qty-008/concentration-rescaling",
        within_conversion_enclosure(amount_concentration.value(), 1.0)
            && within_conversion_enclosure(recovered.value(), mass_concentration.value()),
        &format!(
            "amount_concentration={} recovered_mass_concentration={}",
            amount_concentration.value(),
            recovered.value()
        ),
    );
}

// ---------------------------------------------------------------------------
// G0 property adoption (bead frankensim-4nh8): the dimension algebra laws,
// generated + shrunk via fs-propcheck. The fixed cases above REMAIN as
// regression pins; these properties cover the space between them and
// deliver minimal counterexamples on failure.
// ---------------------------------------------------------------------------

/// Generate a small Dims vector (exponents in [-3, 3] — the physically
/// meaningful range; overflow semantics are a separate documented bound).
fn gen_dims(s: &mut fs_propcheck::Stream) -> Vec<i64> {
    (0..6).map(|_| s.int_in(-3, 3)).collect()
}

fn to_dims(v: &[i64]) -> fs_qty::Dims {
    let mut a = [0i8; 6];
    for (slot, &x) in a.iter_mut().zip(v) {
        *slot = x as i8;
    }
    fs_qty::Dims(a)
}

#[test]
fn g0_dims_plus_commutes_and_minus_inverts() {
    fs_propcheck::check(
        "dims-plus-commutes",
        0x971_0001,
        400,
        |s| (gen_dims(s), gen_dims(s)),
        |(a, b)| {
            let (da, db) = (to_dims(a), to_dims(b));
            da.plus(db) == db.plus(da)
        },
    );
    fs_propcheck::check(
        "dims-minus-inverts-plus",
        0x971_0002,
        400,
        |s| (gen_dims(s), gen_dims(s)),
        |(a, b)| {
            let (da, db) = (to_dims(a), to_dims(b));
            da.plus(db).minus(db) == da
        },
    );
    println!(
        "{{\"suite\":\"fs-qty\",\"case\":\"g0-dims-laws\",\"verdict\":\"pass\",\"detail\":\"800 generated cases, shrink-armed\"}}"
    );
}

#[test]
fn g0_dims_times_distributes_over_plus() {
    fs_propcheck::check(
        "dims-times-distributes",
        0x971_0003,
        400,
        |s| (gen_dims(s), gen_dims(s), s.int_in(-3, 3)),
        |args| {
            let ((a, b), n) = ((&args.0, &args.1), args.2 as i8);
            let (da, db) = (to_dims(a), to_dims(b));
            da.plus(db).times(n) == da.times(n).plus(db.times(n))
        },
    );
    println!(
        "{{\"suite\":\"fs-qty\",\"case\":\"g0-times-distributes\",\"verdict\":\"pass\",\"detail\":\"400 generated cases\"}}"
    );
}
