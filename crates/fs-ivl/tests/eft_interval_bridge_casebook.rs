//! Structured evidence for two public BEDROCK arithmetic seams (bead
//! 6ys.18.9).
//!
//! The EFT case compares the public double-double components consumed by
//! fs-la's mixed-precision ladder with fs-ivl's public expansion components;
//! fs-la's private residual accumulator itself is outside this fixture. The
//! interval case compares point evaluations with fs-math's deterministic
//! elementary results nudged by the documented outward ULP budgets. This
//! cheap G0 slice is not a FrankenScipy oracle battery, the full per-crate
//! fs-ivl suite, a performance result, a full G5 sweep, or fresh dual-ISA
//! execution evidence.

use core::fmt::Write as _;
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_ivl::{
    Interval,
    expansion::{fast_expansion_sum_zeroelim, scale_expansion_zeroelim},
};
use fs_math::{
    dd::Dd,
    det,
    eft::{two_prod, two_sum},
    next_down, next_up,
};

const SUITE: &str = "fs-bedrock/eft-interval-bridge-v1";

const SUM_A: u64 = 0x4341_c379_37e0_8000;
const SUM_B: u64 = 0x3ff0_0000_0000_0000;
const SUM_COMPONENTS: [u64; 2] = [SUM_B, SUM_A];

const PRODUCT_A: u64 = 0x3ff0_0000_0000_0001;
const PRODUCT_B: u64 = 0x3ff0_0000_0000_0001;
const PRODUCT_COMPONENTS: [u64; 2] = [0x3970_0000_0000_0000, 0x3ff0_0000_0000_0002];

#[derive(Debug, Clone, Copy)]
enum Elementary {
    Exp,
    Ln,
    Tanh,
}

impl Elementary {
    const fn name(self) -> &'static str {
        match self {
            Self::Exp => "exp",
            Self::Ln => "ln",
            Self::Tanh => "tanh",
        }
    }

    const fn budget(self) -> u32 {
        match self {
            Self::Exp | Self::Ln => 3,
            Self::Tanh => 5,
        }
    }

    fn deterministic_value(self, x: f64) -> f64 {
        match self {
            Self::Exp => det::exp(x),
            Self::Ln => det::ln(x),
            Self::Tanh => det::tanh(x),
        }
    }

    fn interval_value(self, x: f64) -> Interval {
        match self {
            Self::Exp => Interval::point(x).exp(),
            Self::Ln => Interval::point(x).ln(),
            Self::Tanh => Interval::point(x).tanh(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct ElementaryVector {
    function: Elementary,
    input: u64,
    deterministic: u64,
    endpoints: [u64; 2],
}

const ELEMENTARY_VECTORS: [ElementaryVector; 3] = [
    ElementaryVector {
        function: Elementary::Exp,
        input: 0x0000_0000_0000_0000,
        deterministic: 0x3ff0_0000_0000_0000,
        endpoints: [0x3fef_ffff_ffff_fffd, 0x3ff0_0000_0000_0003],
    },
    ElementaryVector {
        function: Elementary::Ln,
        input: 0x3ff0_0000_0000_0000,
        deterministic: 0x0000_0000_0000_0000,
        endpoints: [0x8000_0000_0000_0003, 0x0000_0000_0000_0003],
    },
    ElementaryVector {
        function: Elementary::Tanh,
        input: 0x0000_0000_0000_0000,
        deterministic: 0x0000_0000_0000_0000,
        endpoints: [0x8000_0000_0000_0005, 0x0000_0000_0000_0005],
    },
];

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

fn push_words(bytes: &mut Vec<u8>, words: &[u64]) {
    push_len(bytes, words.len());
    for &word in words {
        push_u64(bytes, word);
    }
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, frame: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, frame.len());
    bytes.extend_from_slice(frame);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn hex_words(words: &[u64]) -> String {
    let mut encoded = String::from("[");
    for (index, word) in words.iter().enumerate() {
        if index > 0 {
            encoded.push(',');
        }
        write!(&mut encoded, "0x{word:016x}").expect("writing to String cannot fail");
    }
    encoded.push(']');
    encoded
}

fn eft_inputs() -> Vec<u8> {
    let mut bytes = b"fs-bedrock:shared-eft-expansion-components:v1".to_vec();
    push_text(
        &mut bytes,
        "fs_math::eft::{two_sum,two_prod}+fs_math::dd::Dd+fs_ivl::expansion::{fast_expansion_sum_zeroelim,scale_expansion_zeroelim}",
    );
    push_text(
        &mut bytes,
        "component-order=least-significant-to-leading:v1",
    );
    push_text(&mut bytes, "absorbed-sum");
    push_words(&mut bytes, &[SUM_A, SUM_B]);
    push_words(&mut bytes, &SUM_COMPONENTS);
    push_text(&mut bytes, "fma-product-tail");
    push_words(&mut bytes, &[PRODUCT_A, PRODUCT_B]);
    push_words(&mut bytes, &PRODUCT_COMPONENTS);
    bytes
}

fn nudge_down(mut value: f64, steps: u32) -> f64 {
    for _ in 0..steps {
        value = next_down(value);
    }
    value
}

fn nudge_up(mut value: f64, steps: u32) -> f64 {
    for _ in 0..steps {
        value = next_up(value);
    }
    value
}

fn derived_interval_bits(vector: ElementaryVector) -> [u64; 2] {
    let value = f64::from_bits(vector.deterministic);
    [
        nudge_down(value, vector.function.budget()).to_bits(),
        nudge_up(value, vector.function.budget()).to_bits(),
    ]
}

fn interval_inputs() -> Vec<u8> {
    let mut bytes = b"fs-bedrock:deterministic-elementary-interval-budgets:v1".to_vec();
    push_text(
        &mut bytes,
        "fs_math::{det,next_down,next_up}+fs_ivl::Interval::{point,exp,ln,tanh}",
    );
    push_text(&mut bytes, "outward-loop=exactly-k-next-float-steps:v1");
    push_len(&mut bytes, ELEMENTARY_VECTORS.len());
    for vector in ELEMENTARY_VECTORS {
        push_text(&mut bytes, vector.function.name());
        push_u64(&mut bytes, vector.input);
        push_u64(&mut bytes, u64::from(vector.function.budget()));
        push_u64(&mut bytes, vector.deterministic);
        push_words(&mut bytes, &vector.endpoints);
    }
    bytes
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct EftMeasurement {
    sum_eft: Vec<u64>,
    sum_dd: Vec<u64>,
    sum_expansion: Vec<u64>,
    product_eft: Vec<u64>,
    product_dd: Vec<u64>,
    product_expansion: Vec<u64>,
}

fn pair_components(pair: (f64, f64)) -> Vec<u64> {
    vec![pair.1.to_bits(), pair.0.to_bits()]
}

fn dd_components(value: Dd) -> Vec<u64> {
    vec![value.lo.to_bits(), value.hi.to_bits()]
}

fn expansion_components(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn measure_eft() -> EftMeasurement {
    let (sum_a, sum_b) = (f64::from_bits(SUM_A), f64::from_bits(SUM_B));
    let (product_a, product_b) = (f64::from_bits(PRODUCT_A), f64::from_bits(PRODUCT_B));
    EftMeasurement {
        sum_eft: pair_components(two_sum(sum_a, sum_b)),
        sum_dd: dd_components(Dd::from_f64(sum_a) + Dd::from_f64(sum_b)),
        sum_expansion: expansion_components(&fast_expansion_sum_zeroelim(&[sum_a], &[sum_b])),
        product_eft: pair_components(two_prod(product_a, product_b)),
        product_dd: dd_components(Dd::from_f64(product_a) * Dd::from_f64(product_b)),
        product_expansion: expansion_components(&scale_expansion_zeroelim(&[product_a], product_b)),
    }
}

fn eft_outcome() -> CaseOutcome {
    let run = measure_eft();
    let replay = measure_eft();
    let inputs_hex = hex_bytes(&eft_inputs());
    if run != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class")
        .with_evidence("crates/fs-ivl/CONTRACT.md#determinism-class");
    }

    let expected_sum = SUM_COMPONENTS.to_vec();
    if run.sum_eft != expected_sum
        || run.sum_dd != expected_sum
        || run.sum_expansion != expected_sum
    {
        return CaseOutcome::fail(format!(
            "stage=absorbed-sum-components; expected={}; eft={}; dd={}; expansion={}; inputs_hex={inputs_hex}",
            hex_words(&expected_sum),
            hex_words(&run.sum_eft),
            hex_words(&run.sum_dd),
            hex_words(&run.sum_expansion),
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }

    let expected_product = PRODUCT_COMPONENTS.to_vec();
    if run.product_eft != expected_product
        || run.product_dd != expected_product
        || run.product_expansion != expected_product
    {
        return CaseOutcome::fail(format!(
            "stage=fma-product-components; expected={}; eft={}; dd={}; expansion={}; inputs_hex={inputs_hex}",
            hex_words(&expected_product),
            hex_words(&run.product_eft),
            hex_words(&run.product_dd),
            hex_words(&run.product_expansion),
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }

    CaseOutcome::pass(
        "sum_components=[tail,rounded]; product_components=[tail,rounded]; direct_eft=dd=expansion; same_run_bits=identical",
    )
    .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics")
}

fn measure_intervals() -> [[u64; 2]; 3] {
    ELEMENTARY_VECTORS.map(|vector| {
        let interval = vector.function.interval_value(f64::from_bits(vector.input));
        [interval.lo().to_bits(), interval.hi().to_bits()]
    })
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    function: usize,
    endpoint: usize,
    bit: u32,
}

fn interval_outcome(
    reference: [[u64; 2]; 3],
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let run = measure_intervals();
    let replay = measure_intervals();
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |corruption| {
            format!(
                "seed=0x{:016x}; function={}; endpoint={}; bit={}",
                corruption.seed, corruption.function, corruption.endpoint, corruption.bit,
            )
        },
    );
    if run != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={run:016x?}; second={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class")
        .with_evidence("crates/fs-ivl/CONTRACT.md#determinism-class");
    }

    for (index, vector) in ELEMENTARY_VECTORS.into_iter().enumerate() {
        let deterministic = vector
            .function
            .deterministic_value(f64::from_bits(vector.input));
        if deterministic.to_bits() != vector.deterministic {
            return CaseOutcome::fail(format!(
                "{context}; stage=deterministic-elementary-known-answer; function={}; input_bits=0x{:016x}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; all_computed={run:016x?}; inputs_hex={inputs_hex}",
                vector.function.name(),
                vector.input,
                deterministic.to_bits(),
                vector.deterministic,
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#invariants");
        }
        let derived = derived_interval_bits(vector);
        if derived != vector.endpoints {
            return CaseOutcome::fail(format!(
                "{context}; stage=outward-reference-derivation; function={}; deterministic_bits=0x{:016x}; budget={}; derived=[0x{:016x},0x{:016x}]; literal=[0x{:016x},0x{:016x}]; inputs_hex={inputs_hex}",
                vector.function.name(),
                vector.deterministic,
                vector.function.budget(),
                derived[0],
                derived[1],
                vector.endpoints[0],
                vector.endpoints[1],
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#invariants")
            .with_evidence("crates/fs-ivl/CONTRACT.md#invariants");
        }
        if run[index] != reference[index] {
            return CaseOutcome::fail(format!(
                "{context}; stage=elementary-interval-budget; function={}; input_bits=0x{:016x}; deterministic_bits=0x{:016x}; budget={}; computed=[0x{:016x},0x{:016x}]; reference=[0x{:016x},0x{:016x}]; all_computed={run:016x?}; all_reference={reference:016x?}; inputs_hex={inputs_hex}",
                vector.function.name(),
                vector.input,
                deterministic.to_bits(),
                vector.function.budget(),
                run[index][0],
                run[index][1],
                reference[index][0],
                reference[index][1],
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-ivl/CONTRACT.md#invariants");
        }
    }

    CaseOutcome::pass(
        "functions=exp,ln,tanh; point_inputs=3; outward_budgets=3,3,5; endpoint_bits=exact; same_run_bits=identical",
    )
    .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-ivl/CONTRACT.md#invariants")
}

#[test]
fn eft_and_interval_bridges_emit_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let eft_frame = eft_inputs();
    let interval_frame = interval_inputs();
    let eft_digest = fnv1a64(&eft_frame);
    let interval_digest = fnv1a64(&interval_frame);
    assert_eq!(eft_digest, 0xde07_2355_a326_fc99);
    assert_eq!(interval_digest, 0xcccd_cab1_6b4d_6a7a);
    let references = ELEMENTARY_VECTORS.map(|vector| vector.endpoints);

    let report = Suite::new(SUITE)
        .case(
            "shared-eft-expansion-components",
            eft_digest,
            ToleranceSpec::Exact,
            eft_outcome,
        )
        .case(
            "deterministic-elementary-interval-budgets",
            interval_digest,
            ToleranceSpec::Exact,
            move || interval_outcome(references, None, &interval_frame),
        )
        .run();

    report.assert_green();
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "shared-eft-expansion-components",
            "deterministic-elementary-interval-budgets",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-bedrock/eft-interval-bridge-v1\",",
                "\"case\":\"shared-eft-expansion-components\",\"inputs_digest\":\"de072355a326fc99\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"sum_components=[tail,rounded]; product_components=[tail,rounded]; direct_eft=dd=expansion; same_run_bits=identical\",",
                "\"evidence\":[\"crates/fs-math/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-ivl/CONTRACT.md#public-types-and-semantics\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the cross-crate bridge record schema and field order are contract"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Keep the corruption frame and every refusal assertion together.
fn disclosed_seeded_corruption_is_replay_identical_and_turns_the_bridge_red() {
    const CORRUPTION_SEED: u64 = 0xF51C_0001;
    let function = (CORRUPTION_SEED & 0x3) as usize;
    let endpoint = ((CORRUPTION_SEED >> 2) & 0x1) as usize;
    let bit = ((CORRUPTION_SEED >> 32) & 0x3f) as u32;
    assert_eq!((function, endpoint, bit), (1, 0, 0));

    let canonical = ELEMENTARY_VECTORS.map(|vector| vector.endpoints);
    let mut corrupted = canonical;
    corrupted[function][endpoint] ^= 1_u64 << bit;
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        function,
        endpoint,
        bit,
    };

    let interval_frame = interval_inputs();
    let mut inputs = b"fs-bedrock:seeded-interval-reference-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, function);
    push_len(&mut inputs, endpoint);
    push_u64(&mut inputs, u64::from(bit));
    push_nested(
        &mut inputs,
        "nested-elementary-interval-frame",
        &interval_frame,
    );
    push_text(&mut inputs, "canonical-endpoint-bits");
    for endpoints in canonical {
        push_words(&mut inputs, &endpoints);
    }
    push_text(&mut inputs, "corrupted-endpoint-bits");
    for endpoints in corrupted {
        push_words(&mut inputs, &endpoints);
    }
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x9fa2_ea2b_44d2_9082);

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-interval-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || interval_outcome(corrupted, Some(corruption), &input_frame),
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
    assert_eq!(first_failure.case, "seeded-interval-reference-corruption");
    assert_eq!(first_failure.inputs_digest, "9fa2ea2b44d29082");
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains("function=1; endpoint=0; bit=0")
    );
    assert!(first_failure.details.contains("function=ln"));
    assert!(
        first_failure
            .details
            .contains("stage=elementary-interval-budget")
    );
    assert!(first_failure.details.contains("all_computed="));
    assert!(first_failure.details.contains("all_reference="));
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| first.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-interval-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
