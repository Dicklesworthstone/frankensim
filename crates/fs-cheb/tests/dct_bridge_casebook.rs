//! Structured BEDROCK cross-crate evidence for the public FFT-to-Chebyshev
//! seam (bead 6ys.18.8).
//!
//! This cheap PR case synthesizes first-kind values through fs-cheb's public
//! DCT-III bridge, analyzes them through fs-fft's public DCT-II, applies the
//! documented `2/n` scale, and reconstructs a `Cheb1` for point evaluation.
//! It closes only that explicit integration seam. FrankenScipy oracles,
//! Dual-through-fs-la support, the full BEDROCK G5 sweep, and fresh dual-ISA
//! execution evidence remain separate proof work.

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_cheb::{Cheb1, values_from_coeffs};
use fs_fft::{dct2, dct3};
use fs_math::det;

const SUITE: &str = "fs-bedrock/fft-cheb-bridge-v1";
const SAMPLE_COUNT: usize = 8;
const DOMAIN: [f64; 2] = [-1.0, 1.0];
const COEFFICIENTS: [f64; 4] = [2.0, 3.0, 0.0, 5.0];
const EVALUATION_POINTS: [f64; 4] = [-1.0, 0.0, 0.5, 1.0];
const EXPECTED_EVALUATIONS: [f64; 4] = [-7.0, 1.0, -2.5, 9.0];
const ABSOLUTE_TOLERANCE: f64 = 1.0e-11;
const DCT_II_CONVENTION: &str = "X[k]=sum_j(x[j]*cos(pi*k*(2*j+1)/(2*n)));analysis-scale=2/n:v1";
const CHEB_CONVENTION: &str = "first-kind-roots;stored-c0-unhalved;clenshaw-half-c0:v1";
const SAMPLE_GRID_CONVENTION: &str = "t[j]=cos(pi*(j+0.5)/n);j-ascending;roots-descending:v1";

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("conformance fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, value.to_bits());
    }
}

fn push_u64s(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
    }
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, frame: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, frame.len());
    bytes.extend_from_slice(frame);
}

fn bridge_inputs() -> Vec<u8> {
    let mut bytes = b"fs-bedrock:fft-cheb-dct-bridge:v1".to_vec();
    push_text(
        &mut bytes,
        "fs_cheb::values_from_coeffs+fs_fft::dct3+fs_fft::dct2+fs_math::det::cos+Cheb1::from_coeffs+Cheb1::eval",
    );
    push_text(&mut bytes, DCT_II_CONVENTION);
    push_text(&mut bytes, CHEB_CONVENTION);
    push_text(&mut bytes, SAMPLE_GRID_CONVENTION);
    push_len(&mut bytes, SAMPLE_COUNT);
    push_text(&mut bytes, "domain");
    push_f64s(&mut bytes, &DOMAIN);
    push_text(&mut bytes, "stored-c0-unhalved-coefficients");
    push_f64s(&mut bytes, &COEFFICIENTS);
    push_text(&mut bytes, "evaluation-points");
    push_f64s(&mut bytes, &EVALUATION_POINTS);
    push_text(&mut bytes, "expected-source-evaluations");
    push_f64s(&mut bytes, &EXPECTED_EVALUATIONS);
    push_text(&mut bytes, "absolute-tolerance");
    push_u64(&mut bytes, ABSOLUTE_TOLERANCE.to_bits());
    bytes
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn bitwise_equal(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

#[derive(Debug)]
struct BridgeMeasurement {
    synthesis: Vec<f64>,
    direct_dct3: Vec<f64>,
    independent_samples: Vec<f64>,
    recovered: Vec<f64>,
    evaluations: [f64; 4],
    resynthesis: Vec<f64>,
}

fn measure_bridge() -> BridgeMeasurement {
    let synthesis = values_from_coeffs(&COEFFICIENTS, SAMPLE_COUNT);
    let mut padded = COEFFICIENTS.to_vec();
    padded.resize(SAMPLE_COUNT, 0.0);
    let direct_dct3 = dct3(&padded);
    let source = Cheb1::from_coeffs(DOMAIN[0], DOMAIN[1], COEFFICIENTS.to_vec());
    let independent_samples = (0..SAMPLE_COUNT)
        .map(|sample| {
            let angle = std::f64::consts::PI * (sample as f64 + 0.5) / SAMPLE_COUNT as f64;
            source.eval(det::cos(angle))
        })
        .collect();
    let mut recovered = dct2(&synthesis);
    let scale = 2.0 / SAMPLE_COUNT as f64;
    for coefficient in &mut recovered {
        *coefficient *= scale;
    }
    let reconstructed = Cheb1::from_coeffs(DOMAIN[0], DOMAIN[1], recovered.clone());
    let evaluations = EVALUATION_POINTS.map(|point| reconstructed.eval(point));
    let resynthesis = values_from_coeffs(&recovered, SAMPLE_COUNT);
    BridgeMeasurement {
        synthesis,
        direct_dct3,
        independent_samples,
        recovered,
        evaluations,
        resynthesis,
    }
}

fn measurements_replay_bitwise(left: &BridgeMeasurement, right: &BridgeMeasurement) -> bool {
    bitwise_equal(&left.synthesis, &right.synthesis)
        && bitwise_equal(&left.direct_dct3, &right.direct_dct3)
        && bitwise_equal(&left.independent_samples, &right.independent_samples)
        && bitwise_equal(&left.recovered, &right.recovered)
        && bitwise_equal(&left.evaluations, &right.evaluations)
        && bitwise_equal(&left.resynthesis, &right.resynthesis)
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    component: usize,
    bit: u32,
}

fn context(corruption: Option<Corruption>) -> String {
    match corruption {
        Some(corruption) => format!(
            "seed=0x{:016x}; component={}; bit={}",
            corruption.seed, corruption.component, corruption.bit,
        ),
        None => "mode=canonical".to_owned(),
    }
}

#[allow(clippy::too_many_lines)] // One auditable seam keeps every cross-crate stage in log order.
fn bridge_outcome(reference_coefficients: [f64; 4], corruption: Option<Corruption>) -> CaseOutcome {
    let run = measure_bridge();
    let replay = measure_bridge();
    let context = context(corruption);

    if !measurements_replay_bitwise(&run, &replay) {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; synthesis_a={:016x?}; synthesis_b={:016x?}; independent_samples_a={:016x?}; independent_samples_b={:016x?}; recovered_a={:016x?}; recovered_b={:016x?}; evaluations_a={:016x?}; evaluations_b={:016x?}",
            bits(&run.synthesis),
            bits(&replay.synthesis),
            bits(&run.independent_samples),
            bits(&replay.independent_samples),
            bits(&run.recovered),
            bits(&replay.recovered),
            bits(&run.evaluations),
            bits(&replay.evaluations),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#determinism-class")
        .with_evidence("crates/fs-cheb/CONTRACT.md#determinism-class");
    }

    if !bitwise_equal(&run.synthesis, &run.direct_dct3) {
        return CaseOutcome::fail(format!(
            "{context}; stage=public-dct3-synthesis-identity; values_from_coeffs_bits={:016x?}; direct_dct3_bits={:016x?}",
            bits(&run.synthesis),
            bits(&run.direct_dct3),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }

    if run.synthesis.len() != run.independent_samples.len() {
        return CaseOutcome::fail(format!(
            "{context}; stage=first-kind-sample-shape; synthesis_count={}; independent_count={}; n={SAMPLE_COUNT}",
            run.synthesis.len(),
            run.independent_samples.len(),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }
    for (sample, (&computed, &reference)) in run
        .synthesis
        .iter()
        .zip(&run.independent_samples)
        .enumerate()
    {
        let absolute_error = (computed - reference).abs();
        if !absolute_error.is_finite() || absolute_error > ABSOLUTE_TOLERANCE {
            return CaseOutcome::fail(format!(
                "{context}; stage=first-kind-sample-semantics; convention={SAMPLE_GRID_CONVENTION}; sample={sample}; computed={computed}; independent={reference}; computed_bits=0x{:016x}; independent_bits=0x{:016x}; absolute_error={absolute_error}; tolerance={ABSOLUTE_TOLERANCE}; synthesis_bits={:016x?}; independent_bits={:016x?}",
                computed.to_bits(),
                reference.to_bits(),
                bits(&run.synthesis),
                bits(&run.independent_samples),
            ))
            .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
        }
    }

    let source = Cheb1::from_coeffs(DOMAIN[0], DOMAIN[1], COEFFICIENTS.to_vec());
    let source_evaluations = EVALUATION_POINTS.map(|point| source.eval(point));
    if !bitwise_equal(&source_evaluations, &EXPECTED_EVALUATIONS) {
        return CaseOutcome::fail(format!(
            "{context}; stage=source-chebyshev-known-answer; computed_bits={:016x?}; reference_bits={:016x?}; computed={source_evaluations:?}; reference={EXPECTED_EVALUATIONS:?}",
            bits(&source_evaluations),
            bits(&EXPECTED_EVALUATIONS),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }

    let mut padded_reference = reference_coefficients.to_vec();
    padded_reference.resize(SAMPLE_COUNT, 0.0);
    if run.recovered.len() != padded_reference.len() {
        return CaseOutcome::fail(format!(
            "{context}; stage=dct2-analysis-shape; computed_count={}; reference_count={}; n={SAMPLE_COUNT}",
            run.recovered.len(),
            padded_reference.len(),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }
    for (component, (&computed, &reference)) in
        run.recovered.iter().zip(&padded_reference).enumerate()
    {
        let absolute_error = (computed - reference).abs();
        if !absolute_error.is_finite() || absolute_error > ABSOLUTE_TOLERANCE {
            return CaseOutcome::fail(format!(
                "{context}; stage=dct2-analysis; convention={DCT_II_CONVENTION}; component={component}; computed={computed}; reference={reference}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; absolute_error={absolute_error}; tolerance={ABSOLUTE_TOLERANCE}; recovered_bits={:016x?}; synthesis_bits={:016x?}",
                computed.to_bits(),
                reference.to_bits(),
                bits(&run.recovered),
                bits(&run.synthesis),
            ))
            .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
        }
    }

    for (component, (&computed, &reference)) in run
        .evaluations
        .iter()
        .zip(&EXPECTED_EVALUATIONS)
        .enumerate()
    {
        let absolute_error = (computed - reference).abs();
        if !absolute_error.is_finite() || absolute_error > ABSOLUTE_TOLERANCE {
            return CaseOutcome::fail(format!(
                "{context}; stage=reconstructed-chebyshev-evaluation; convention={CHEB_CONVENTION}; point={}; component={component}; computed={computed}; reference={reference}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; absolute_error={absolute_error}; tolerance={ABSOLUTE_TOLERANCE}; evaluation_bits={:016x?}; recovered_bits={:016x?}",
                EVALUATION_POINTS[component],
                computed.to_bits(),
                reference.to_bits(),
                bits(&run.evaluations),
                bits(&run.recovered),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
        }
    }

    if run.resynthesis.len() != run.synthesis.len() {
        return CaseOutcome::fail(format!(
            "{context}; stage=dct-sample-roundtrip-shape; resynthesis_count={}; synthesis_count={}; n={SAMPLE_COUNT}",
            run.resynthesis.len(),
            run.synthesis.len(),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
        .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
    }
    for (sample, (&computed, &reference)) in run.resynthesis.iter().zip(&run.synthesis).enumerate()
    {
        let absolute_error = (computed - reference).abs();
        if !absolute_error.is_finite() || absolute_error > ABSOLUTE_TOLERANCE {
            return CaseOutcome::fail(format!(
                "{context}; stage=dct-sample-roundtrip; sample={sample}; computed={computed}; reference={reference}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; absolute_error={absolute_error}; tolerance={ABSOLUTE_TOLERANCE}; resynthesis_bits={:016x?}; synthesis_bits={:016x?}",
                computed.to_bits(),
                reference.to_bits(),
                bits(&run.resynthesis),
                bits(&run.synthesis),
            ))
            .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
        }
    }

    CaseOutcome::pass(
        "n=8; synthesis=values_from_coeffs==dct3-bitwise; first_kind_samples=source_eval<=1e-11; analysis_scale=2/n; coefficient_error<=1e-11; evaluation_error<=1e-11; sample_roundtrip_error<=1e-11; same_run_bits=identical",
    )
    .with_evidence("crates/fs-fft/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn fft_to_chebyshev_bridge_emits_a_replay_complete_green_record() {
    let inputs_digest = fnv1a64(&bridge_inputs());
    assert_eq!(inputs_digest, 0x135b_0e52_a5ed_3d5f);

    let report = Suite::new(SUITE)
        .case(
            "dct-analysis-chebyshev-reconstruction",
            inputs_digest,
            ToleranceSpec::AbsoluteLe(ABSOLUTE_TOLERANCE),
            || bridge_outcome(COEFFICIENTS, None),
        )
        .run();

    report.assert_green();
    assert_eq!(report.records.len(), 1);
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-bedrock/fft-cheb-bridge-v1\",",
                "\"case\":\"dct-analysis-chebyshev-reconstruction\",\"inputs_digest\":\"135b0e52a5ed3d5f\",",
                "\"tolerance\":\"abs<=1e-11\",\"pass\":true,",
                "\"details\":\"n=8; synthesis=values_from_coeffs==dct3-bitwise; first_kind_samples=source_eval<=1e-11; analysis_scale=2/n; coefficient_error<=1e-11; evaluation_error<=1e-11; sample_roundtrip_error<=1e-11; same_run_bits=identical\",",
                "\"evidence\":[\"crates/fs-fft/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-cheb/CONTRACT.md#public-types-and-semantics\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the cross-crate bridge record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_is_replay_identical_and_turns_the_bridge_red() {
    const CORRUPTION_SEED: u64 = 0xF5BC_0001;
    let component = (CORRUPTION_SEED & 0x3) as usize;
    let bit = 40 + ((CORRUPTION_SEED >> 2) & 0x7) as u32;
    assert_eq!(component, 1);
    assert_eq!(bit, 40);
    let canonical_bits = COEFFICIENTS.map(f64::to_bits);
    let mut corrupted_bits = canonical_bits;
    corrupted_bits[component] ^= 1_u64 << bit;
    assert_eq!(corrupted_bits[component], 0x4008_0100_0000_0000);
    let corrupted = corrupted_bits.map(f64::from_bits);
    assert!((corrupted[component] - COEFFICIENTS[component]).abs() > ABSOLUTE_TOLERANCE);
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        component,
        bit,
    };

    let bridge = bridge_inputs();
    let mut inputs = b"fs-bedrock:seeded-fft-cheb-reference-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, component);
    push_u64(&mut inputs, u64::from(bit));
    push_nested(&mut inputs, "nested-fft-cheb-bridge-frame", &bridge);
    push_text(&mut inputs, "canonical-coefficient-bits");
    push_u64s(&mut inputs, &canonical_bits);
    push_text(&mut inputs, "corrupted-coefficient-bits");
    push_u64s(&mut inputs, &corrupted_bits);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x30d3_f3f7_8100_831c);

    let make_report = || {
        Suite::new(SUITE)
            .case(
                "seeded-dct-reference-corruption",
                inputs_digest,
                ToleranceSpec::AbsoluteLe(ABSOLUTE_TOLERANCE),
                move || bridge_outcome(corrupted, Some(corruption)),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();

    assert!(!first.all_passed());
    assert!(!replay.all_passed());
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one structured failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed corruption must produce exactly one structured failure");
    };
    assert_eq!(first_failure.case, "seeded-dct-reference-corruption");
    assert_eq!(first_failure.inputs_digest, "30d3f3f78100831c");
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains(&format!("component={component}; bit={bit}"))
    );
    assert!(first_failure.details.contains("stage=dct2-analysis"));
    assert!(first_failure.details.contains("recovered_bits=["));
    assert!(first_failure.details.contains("synthesis_bits=["));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"abs<=1e-11\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| first.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-dct-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
