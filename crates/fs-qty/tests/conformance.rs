//! Conformance suite for fs-qty (plan §13.3): any reimplementation must pass
//! these cases. Each case logs a JSON-lines verdict so failures are
//! diagnosable from output alone (observability standard; fs-obs schema
//! adoption pending that crate's landing).

use fs_qty::parse::{ParseBudget, ParseErrorKind, ParseResource, parse_qty, parse_qty_with_budget};
use fs_qty::semantic::{
    AngleDomain, CompositionBasis, FormRequirement, PhasorAmplitude, PhasorQty, QuantityKind,
    SemanticError, SemanticQty, SemanticType, StrainBasis, StrainComponent, ValueForm,
    amount_concentration_to_mass_concentration, amount_to_mass,
    mass_concentration_to_amount_concentration, mass_to_amount,
};
use fs_qty::{
    Dims, DynViscosity, Length, Pressure, QUANTITY_SPEC_ENCODED_LEN, QtyAny, QuantitySpec,
    QuantitySpecDecodeError, Time,
};

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
            error.verifies_source(input)
                && error.input_bytes == input.len()
                && error.preview.len() <= ParseBudget::DEFAULT.max_diagnostic_bytes()
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

/// qty-005b: hostile sources, factors, and tokens are bounded before retained
/// diagnostics can scale with the rejected input.
#[test]
fn qty_005b_parser_budget_and_diagnostic_bounds() {
    let budget = ParseBudget::new(96, 3, 8, 12);
    let at_factor_cap = "1rad*rad*rad";
    assert!(parse_qty_with_budget(at_factor_cap, budget).is_ok());

    let cases = [
        (format!("1m{}", " ".repeat(95)), ParseResource::InputBytes),
        ("1rad*rad*rad*rad".to_string(), ParseResource::Factors),
        ("1abcdefghijk".to_string(), ParseResource::TokenBytes),
        ("1m^000000000".to_string(), ParseResource::TokenBytes),
    ];
    for (input, resource) in cases {
        let first = parse_qty_with_budget(&input, budget).expect_err("budget must refuse");
        let second = parse_qty_with_budget(&input, budget).expect_err("repeat must refuse");
        let pass = first == second
            && first.preview.len() <= budget.max_diagnostic_bytes()
            && first.input_bytes == input.len()
            && matches!(
                &first.kind,
                ParseErrorKind::BudgetExceeded {
                    resource: actual,
                    ..
                } if *actual == resource
            );
        verdict(
            "qty-005b/bounded-refusal",
            pass,
            &format!(
                "resource={resource:?} input_bytes={} preview_bytes={} hash_present={}",
                input.len(),
                first.preview.len(),
                first.source_hash.is_some()
            ),
        );
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

fn quantity_spec_kinds() -> [QuantityKind; 26] {
    [
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

/// qty-009: one reusable, injective identity token covers every dimensional
/// and semantic scalar descriptor without laundering dimension-only aliases.
#[test]
fn qty_009_quantity_spec_codec_is_canonical_and_injective() {
    let dimensional = QuantitySpec::dimensional(Dims([1, -2, 3, -4, 5, -6]));
    let dimensional_bytes = dimensional.canonical_bytes();
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let strain = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Shear,
        },
        ValueForm::Static,
    ));
    let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let temperature_difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));

    let goldens = dimensional_bytes == [1, 1, 254, 3, 252, 5, 250, 0, 0, 0, 0, 0]
        && pressure.canonical_bytes() == [1, 255, 1, 254, 0, 0, 0, 1, 7, 0, 0, 1]
        && strain.canonical_bytes() == [1, 0, 0, 0, 0, 0, 0, 1, 9, 2, 2, 1];
    let temperatures_are_distinct = QuantitySpec::dimensional(Dims([0, 0, 0, 1, 0, 0]))
        != absolute_temperature
        && absolute_temperature != temperature_difference
        && absolute_temperature.canonical_bytes() != temperature_difference.canonical_bytes();

    let mut encodings = std::collections::BTreeSet::new();
    encodings.insert(dimensional_bytes);
    let mut all_round_trip = true;
    let mut semantic_cases = 0usize;
    for kind in quantity_spec_kinds() {
        for form in [
            ValueForm::Static,
            ValueForm::Instantaneous,
            ValueForm::Peak,
            ValueForm::Rms,
        ] {
            let spec = QuantitySpec::semantic(SemanticType::new(kind, form));
            let bytes = spec.canonical_bytes();
            all_round_trip &= QuantitySpec::from_canonical_bytes(&bytes) == Ok(spec);
            all_round_trip &= encodings.insert(bytes);
            semantic_cases += 1;
        }
    }
    let extreme_dims = QuantitySpec::dimensional(Dims([i8::MIN, -1, 0, 1, i8::MAX, 42]));
    let extreme_round_trip =
        QuantitySpec::from_canonical_bytes(&extreme_dims.canonical_bytes()) == Ok(extreme_dims);

    verdict(
        "qty-009/canonical-injective",
        QUANTITY_SPEC_ENCODED_LEN == 12
            && goldens
            && temperatures_are_distinct
            && semantic_cases == 104
            && encodings.len() == 105
            && all_round_trip
            && extreme_round_trip,
        &format!(
            "semantic_cases={semantic_cases} unique_tokens={} dimensional={dimensional_bytes:?}",
            encodings.len()
        ),
    );
}

/// qty-009b: the decoder accepts only canonical schema tokens. A descriptor
/// can name an unsupported kind/form for diagnostics, but value admission
/// remains independently fail-closed through `SemanticQty`.
#[test]
#[allow(clippy::too_many_lines)] // One mutation matrix exercises every codec refusal class.
fn qty_009b_quantity_spec_decode_refuses_aliases() {
    let dimensional = QuantitySpec::dimensional(Dims([1, -2, 3, -4, 5, -6]));
    let dimensional_bytes = dimensional.canonical_bytes();
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let pressure_bytes = pressure.canonical_bytes();
    let strain = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Shear,
        },
        ValueForm::Static,
    ));

    let short = matches!(
        QuantitySpec::from_canonical_bytes(&dimensional_bytes[..11]),
        Err(QuantitySpecDecodeError::Length { actual: 11 })
    );
    let long = matches!(
        QuantitySpec::from_canonical_bytes(&[0; 13]),
        Err(QuantitySpecDecodeError::Length { actual: 13 })
    );
    let mut wrong_version = dimensional_bytes;
    wrong_version[0] = 2;
    let version = matches!(
        QuantitySpec::from_canonical_bytes(&wrong_version),
        Err(QuantitySpecDecodeError::Version { actual: 2 })
    );
    let mut wrong_variant = dimensional_bytes;
    wrong_variant[7] = 2;
    let variant = matches!(
        QuantitySpec::from_canonical_bytes(&wrong_variant),
        Err(QuantitySpecDecodeError::Variant { actual: 2 })
    );
    let mut dimensional_payload = dimensional_bytes;
    dimensional_payload[10] = 1;
    let payload = matches!(
        QuantitySpec::from_canonical_bytes(&dimensional_payload),
        Err(QuantitySpecDecodeError::DimensionalPayload {
            index: 10,
            value: 1
        })
    );
    let mut unknown_kind = pressure_bytes;
    unknown_kind[8] = 255;
    let kind = matches!(
        QuantitySpec::from_canonical_bytes(&unknown_kind),
        Err(QuantitySpecDecodeError::SemanticKind { kind: 255, .. })
    );
    let mut noncanonical_parameters = pressure_bytes;
    noncanonical_parameters[9] = 1;
    let parameters = matches!(
        QuantitySpec::from_canonical_bytes(&noncanonical_parameters),
        Err(QuantitySpecDecodeError::SemanticKind {
            kind: 7,
            parameter_a: 1,
            parameter_b: 0
        })
    );
    let mut unknown_form = pressure_bytes;
    unknown_form[11] = 0;
    let form = matches!(
        QuantitySpec::from_canonical_bytes(&unknown_form),
        Err(QuantitySpecDecodeError::ValueForm { actual: 0 })
    );
    let mut wrong_semantic_dims = pressure_bytes;
    wrong_semantic_dims[1] = 0;
    let dimensions = matches!(
        QuantitySpec::from_canonical_bytes(&wrong_semantic_dims),
        Err(QuantitySpecDecodeError::SemanticDimensions {
            actual: Dims([0, 1, -2, 0, 0, 0]),
            expected: Dims([-1, 1, -2, 0, 0, 0]),
            ..
        })
    );

    let mut mutation_law = true;
    for canonical in [dimensional_bytes, pressure_bytes, strain.canonical_bytes()] {
        for index in 0..canonical.len() {
            for value in u8::MIN..=u8::MAX {
                let mut mutated = canonical;
                mutated[index] = value;
                if let Ok(decoded) = QuantitySpec::from_canonical_bytes(&mutated) {
                    mutation_law &= decoded.canonical_bytes() == mutated;
                }
            }
        }
    }

    let unsupported_type = SemanticType::new(QuantityKind::Energy, ValueForm::Rms);
    let unsupported_spec = QuantitySpec::semantic(unsupported_type);
    let descriptor_round_trip =
        QuantitySpec::from_canonical_bytes(&unsupported_spec.canonical_bytes())
            == Ok(unsupported_spec);
    let value_refused = matches!(
        SemanticQty::new(QtyAny::new(1.0, unsupported_type.expected_dims()), unsupported_type),
        Err(SemanticError::UnsupportedForm { source, .. }) if source == unsupported_type
    );

    verdict(
        "qty-009/noncanonical-refusal",
        short
            && long
            && version
            && variant
            && payload
            && kind
            && parameters
            && form
            && dimensions
            && mutation_law
            && descriptor_round_trip
            && value_refused,
        &format!(
            "lengths={short}/{long} tags={version}/{variant}/{kind}/{form} canonical_mutations={mutation_law} value_refused={value_refused}"
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
    // The algebra laws hold EXACTLY on the checked authoritative ops
    // (bead sj31i.11) — saturation broke plus/minus inversion, which
    // is precisely how a clamped exponent aliased false physics.
    fs_propcheck::check(
        "dims-plus-commutes",
        0x971_0001,
        400,
        |s| (gen_dims(s), gen_dims(s)),
        |(a, b)| {
            let (da, db) = (to_dims(a), to_dims(b));
            da.checked_plus(db) == db.checked_plus(da)
        },
    );
    fs_propcheck::check(
        "dims-minus-inverts-plus",
        0x971_0002,
        400,
        |s| (gen_dims(s), gen_dims(s)),
        |(a, b)| {
            let (da, db) = (to_dims(a), to_dims(b));
            da.checked_plus(db)
                .is_none_or(|sum| sum.checked_minus(db) == Some(da))
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
            let lhs = da.checked_plus(db).and_then(|sum| sum.checked_times(n));
            let rhs = da
                .checked_times(n)
                .zip(db.checked_times(n))
                .and_then(|(x, y)| x.checked_plus(y));
            match (lhs, rhs) {
                (Some(x), Some(y)) => x == y,
                _ => true, // distributivity is claimed only where both sides are defined
            }
        },
    );
    println!(
        "{{\"suite\":\"fs-qty\",\"case\":\"g0-times-distributes\",\"verdict\":\"pass\",\"detail\":\"400 generated cases\"}}"
    );
}

/// Bead sj31i.11 boundary battery: every i8 edge refuses typed on the
/// authoritative ops, saturate-then-cancel aliasing is impossible, the
/// diagnostic ops keep their explicitly non-authoritative clamping, and
/// runtime refusal matches the compile-time const-eval refusal class.
#[test]
fn g0_checked_dimension_authority_refuses_every_i8_boundary() {
    use fs_qty::{Dims, QtyAny};

    let hi = Dims([127, 0, 0, 0, 0, 0]);
    let lo = Dims([-128, 0, 0, 0, 0, 0]);
    let one = Dims([1, 0, 0, 0, 0, 0]);

    // +127 + 1 and -128 - 1 refuse on every component position.
    for position in 0..6 {
        let mut top = [0i8; 6];
        top[position] = 127;
        let mut bottom = [0i8; 6];
        bottom[position] = -128;
        let mut unit = [0i8; 6];
        unit[position] = 1;
        assert_eq!(Dims(top).checked_plus(Dims(unit)), None);
        assert_eq!(Dims(bottom).checked_minus(Dims(unit)), None);
        assert_eq!(Dims(bottom).checked_times(-1), None); // -128 * -1
        assert_eq!(Dims(top).checked_times(2), None);
    }
    // Exact boundaries ADMIT.
    assert_eq!(hi.checked_plus(Dims::NONE), Some(hi));
    assert_eq!(lo.checked_minus(Dims::NONE), Some(lo));
    assert_eq!(lo.checked_times(1), Some(lo));

    // THE ALIASING KILL: saturate-then-cancel re-enters the valid range
    // with wrong physics; the checked chain refuses at the overflow
    // instead of laundering 97 for the true 120.
    let a = Dims([100, 0, 0, 0, 0, 0]);
    let b = Dims([50, 0, 0, 0, 0, 0]);
    let c = Dims([30, 0, 0, 0, 0, 0]);
    let laundered = a.saturating_plus(b).saturating_minus(c);
    assert_eq!(laundered, Dims([97, 0, 0, 0, 0, 0]));
    assert_ne!(laundered.0[0] as i32, 100 + 50 - 30);
    assert_eq!(a.checked_plus(b), None);
    assert_eq!(a.checked_plus(b).and_then(|s| s.checked_minus(c)), None);

    // Long multiply/divide trees through QtyAny: overflow refuses typed
    // mid-chain; a dividing tail cannot resurrect the expression.
    let big = QtyAny::new(2.0, Dims([100, 0, 0, 0, 0, 0]));
    let boost = QtyAny::new(3.0, Dims([40, 0, 0, 0, 0, 0]));
    let refusal = big.try_mul(boost).expect_err("140 exceeds i8");
    assert_eq!(refusal.op, "mul");
    assert_eq!(refusal.rhs, Some(Dims([40, 0, 0, 0, 0, 0])));
    let chain = big
        .try_mul(QtyAny::new(1.0, Dims([27, 0, 0, 0, 0, 0])))
        .expect("127 is the exact ceiling")
        .try_mul(QtyAny::new(1.0, one))
        .expect_err("128 leaves i8");
    assert_eq!(chain.op, "mul");
    let divides = QtyAny::new(1.0, lo)
        .try_div(QtyAny::new(1.0, one))
        .expect_err("-129 leaves i8");
    assert_eq!(divides.op, "div");

    // Runtime/const parity: the same ceiling that refuses here refuses
    // at COMPILE TIME on the const-generic path (`Qty<{127}, …> *
    // Qty<{1}, …>` is a const-eval overflow error — see the
    // compile-fail coverage in the crate doctests).
    assert_eq!(hi.checked_plus(one), None);
}
