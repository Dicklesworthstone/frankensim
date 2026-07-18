//! G0/G3/G5 full-study replay for production Monte Carlo hypervolume.
//!
//! Two dyadic two-box fixtures isolate both sampling regimes. The 10-objective
//! fixture makes only Sobol dimensions 8 and 9 decision-relevant; the
//! 12-objective fixture makes the ten-dimensional Sobol head
//! dominance-neutral and isolates the Philox tail in dimensions 10 and 11.
//! Each fixture has an exact inclusion-exclusion volume, while the retained
//! public result binds the estimate bits and hit count. Independent runs must
//! reproduce the complete canonical frame. A disclosed evidence-generator
//! stream flips one returned estimate mantissa bit; stale, resealed-reference,
//! semantic-accounting, stable-red, and merge gates all refuse it. Separate
//! supplied-digest-integrity and retained-fixture-identity tripwires prove
//! those admission causes remain distinguishable.
//!
//! This target does not claim an exact or certified Monte Carlo volume, an IID
//! standard error or formal confidence guarantee, an independent Sobol/Philox
//! implementation, sample-trace or checkpoint replay, arbitrary
//! inputs/seeds/configurations, cross-process or cross-ISA equality,
//! cancellation, persistence, authentication, performance, optimizer quality,
//! or archive correctness.

#![deny(unsafe_code)]

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::mc_hypervolume;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/mc-hypervolume-study-replay";
const CASE: &str = "sobol-head-and-philox-tail-dyadic-union";
const RED_CASE: &str = "seeded-returned-estimate-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-dfo-mc-hypervolume-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-dfo-mc-hypervolume-result-v1";
const FIXTURE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.mc-hypervolume-fixture.v1";
const RESULT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.mc-hypervolume-result.v1";
const EVENT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.mc-hypervolume-event.v1";
const SOURCE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.mc-hypervolume-source.v1";
const SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN: &str =
    "frankensim.fs-dfo.mc-hypervolume-supplied-fixture-digest-tripwire.v1";
const MC_HYPERVOLUME_SOURCE: &[u8] = include_bytes!("../src/moo.rs");

const SOBOL_DIMENSIONS: usize = 10;
const HYBRID_DIMENSIONS: usize = 12;
const SAMPLES: usize = 65_536;
const REFERENCE_COORDINATE: f64 = 1.0;
const NEUTRAL_COORDINATE: f64 = 0.5;
const LOW_COORDINATE: f64 = 0.25;
const HIGH_COORDINATE: f64 = 0.75;

// With the first d-2 coordinates fixed at 1/2 and the final two boxes
// [1/4, 3/4] and [3/4, 1/4], inclusion-exclusion gives
// (3/16 + 3/16 - 1/16) * (1/2)^(d-2).
const SOBOL_BOX_VOLUME: f64 = 9.0 / 4096.0;
const SOBOL_TRUE_VOLUME: f64 = 5.0 / 4096.0;
const HYBRID_BOX_VOLUME: f64 = 9.0 / 16_384.0;
const HYBRID_TRUE_VOLUME: f64 = 5.0 / 16_384.0;
const HEURISTIC_ANALYTIC_BAND_MULTIPLIER: f64 = 8.0;

// Mirrored from the production loop in fs_dfo::mc_hypervolume. The fixture
// identity also binds the exact source bytes so this declaration never floats
// free of the implementation used to produce the evidence.
const PRODUCTION_SOBOL_FIRST_POINT_INDEX: u32 = 1;
const PRODUCTION_PHILOX_TAIL_KERNEL: u32 = 0x0871;
const PRODUCTION_PHILOX_TAIL_TILE: u32 = 0;

const SOBOL_SEED: u64 = 0x4D43_4856_0000_0048;
const HYBRID_SEED: u64 = 0x4D43_4856_0000_1048;
const MUTATION_SEED: u64 = 0x4D43_FA11_0000_0048;
const MUTATION_KERNEL: u32 = 0x4848;
const MUTATION_TILE: u32 = 0;
const MUTATION_BIT_COUNT: u64 = 8;
const WRONG_FIXTURE_SOBOL_SEED: u64 = SOBOL_SEED ^ 1;

const _: () = assert!(SOBOL_DIMENSIONS == 10);
const _: () = assert!(HYBRID_DIMENSIONS == SOBOL_DIMENSIONS + 2);
const _: () = assert!(SAMPLES.is_power_of_two());
const _: () = assert!(SOBOL_TRUE_VOLUME < SOBOL_BOX_VOLUME);
const _: () = assert!(HYBRID_TRUE_VOLUME < HYBRID_BOX_VOLUME);

#[derive(Debug, Clone, Copy, PartialEq)]
struct FixtureSpec {
    label: &'static str,
    dimensions: usize,
    seed: u64,
    box_volume: f64,
    true_volume: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct EstimateBits {
    label: &'static str,
    dimensions: usize,
    estimate: u64,
    hits: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    estimates: Vec<EstimateBits>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    fixture_digest: ContentHash,
    record: StudyRecord,
    result: ReplayIdentity,
    result_digest: ContentHash,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    SuppliedFixtureDigestMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    RetainedFixtureIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
    ResultPayloadIdentityMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    ReferenceIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
    SemanticInconsistency(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    estimate_index: usize,
    mantissa_bit: u32,
    selector_draws: u64,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_error: AdmissionError,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FixtureTripwires {
    supplied_digest_run: StudyRun,
    supplied_digest_error: AdmissionError,
    wrong_fixture_run: StudyRun,
    wrong_fixture_error: AdmissionError,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixed MC-hypervolume fixture cardinality fits u64")
}

fn digest_bytes(digest: ContentHash) -> [u8; 32] {
    *digest.as_bytes()
}

fn fixture_specs() -> [FixtureSpec; 2] {
    assert_eq!(
        fs_rand::qmc::MAX_SOBOL_DIM,
        SOBOL_DIMENSIONS,
        "fixture partition must follow the embedded Sobol table boundary"
    );
    [
        FixtureSpec {
            label: "sobol-head-10d",
            dimensions: SOBOL_DIMENSIONS,
            seed: SOBOL_SEED,
            box_volume: SOBOL_BOX_VOLUME,
            true_volume: SOBOL_TRUE_VOLUME,
        },
        FixtureSpec {
            label: "philox-tail-12d",
            dimensions: HYBRID_DIMENSIONS,
            seed: HYBRID_SEED,
            box_volume: HYBRID_BOX_VOLUME,
            true_volume: HYBRID_TRUE_VOLUME,
        },
    ]
}

fn two_box_front(dimensions: usize) -> Vec<Vec<f64>> {
    assert!(dimensions >= 2);
    let mut first = vec![NEUTRAL_COORDINATE; dimensions];
    let mut second = first.clone();
    first[dimensions - 2] = LOW_COORDINATE;
    first[dimensions - 1] = HIGH_COORDINATE;
    second[dimensions - 2] = HIGH_COORDINATE;
    second[dimensions - 1] = LOW_COORDINATE;
    vec![first, second]
}

fn reference(dimensions: usize) -> Vec<f64> {
    vec![REFERENCE_COORDINATE; dimensions]
}

fn fixture_identity() -> ReplayIdentity {
    fixture_identity_for(fixture_specs())
}

fn production_source_digest() -> ContentHash {
    hash_domain(SOURCE_DIGEST_DOMAIN, MC_HYPERVOLUME_SOURCE)
}

#[allow(clippy::too_many_lines)] // One canonical frame keeps every bound input ordered together.
fn fixture_identity_for(specs: [FixtureSpec; 2]) -> ReplayIdentity {
    let source_digest = production_source_digest();
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_dfo::mc_hypervolume")
        .str("algorithm-randomness", "seeded-scrambled-Sobol-plus-Philox-tail")
        .str("objective-semantics", "minimization-dominated-volume")
        .str("coordinate-units", "dimensionless")
        .str("reference-semantics", "strictly-worse-componentwise")
        .str("sampling-box", "componentwise-front-minimum-to-reference")
        .str("returned-fields", "estimate;hit-count")
        .str("estimate-accounting", "box-volume*hits/sample-count")
        .str("oracle", "dyadic-two-box-inclusion-exclusion")
        .u64("fixture-count", usize_u64(specs.len()))
        .u64("samples-per-fixture", usize_u64(SAMPLES))
        .u64(
            "maximum-sobol-dimensions",
            usize_u64(fs_rand::qmc::MAX_SOBOL_DIM),
        )
        .str(
            "production-sobol-point-schedule",
            "zero-based-sample-s;Sobol::point(s+1);point-zero-skipped",
        )
        .u64(
            "production-sobol-first-point-index",
            u64::from(PRODUCTION_SOBOL_FIRST_POINT_INDEX),
        )
        .str("production-sobol-seed-source", "fixture-algorithm-seed")
        .str(
            "production-philox-tail-schedule",
            "sample-major;increasing-objective-index;one-next_f64-per-tail-coordinate",
        )
        .str(
            "production-philox-tail-constructor",
            "fs_rand::StreamKey{seed,kernel,tile}.stream()",
        )
        .str(
            "production-philox-tail-seed-source",
            "fixture-algorithm-seed",
        )
        .u64(
            "production-philox-tail-kernel",
            u64::from(PRODUCTION_PHILOX_TAIL_KERNEL),
        )
        .u64(
            "production-philox-tail-tile",
            u64::from(PRODUCTION_PHILOX_TAIL_TILE),
        )
        .u64(
            "production-hybrid-tail-draws-per-sample",
            usize_u64(HYBRID_DIMENSIONS - SOBOL_DIMENSIONS),
        )
        .str(
            "analytic-sanity-band-semantics",
            "heuristic-only;binomial-scale-shape-plus-one-hit-resolution;randomized-QMC-not-IID;not-a-standard-error-or-confidence-bound",
        )
        .f64_bits(
            "heuristic-analytic-band-multiplier",
            HEURISTIC_ANALYTIC_BAND_MULTIPLIER,
        )
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "fs-rand-stream-position-domain",
            fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
        )
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .str(
            "mutation-selector",
            "fs_rand::StreamKey::next_below;estimate-index-then-low-mantissa-bit",
        )
        .u64("mutation-mantissa-bit-count", MUTATION_BIT_COUNT)
        .u64("wrong-fixture-sobol-seed", WRONG_FIXTURE_SOBOL_SEED)
        .bytes("mc-hypervolume-source-blake3", source_digest.as_bytes())
        .str("mc-hypervolume-source-digest-domain", SOURCE_DIGEST_DOMAIN)
        .str(
            "supplied-fixture-digest-tripwire-domain",
            SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN,
        )
        .str("fixture-digest-domain", FIXTURE_DIGEST_DOMAIN)
        .str("result-digest-domain", RESULT_DIGEST_DOMAIN)
        .str("event-digest-domain", EVENT_DIGEST_DOMAIN)
        .str(
            "no-claims",
            "exact-MC;certified-MC;IID-standard-error;QMC-error-bound;formal-confidence;independent-sampler;independent-schedule-verifier;sample-trace;checkpoint;arbitrary-inputs;cross-process;cross-ISA;Cx;persistence;authentication;performance;optimizer-quality;archive",
        );
    for (fixture_index, spec) in specs.into_iter().enumerate() {
        builder = builder
            .u64("fixture-index", usize_u64(fixture_index))
            .str("fixture-label", spec.label)
            .u64("dimensions", usize_u64(spec.dimensions))
            .u64("algorithm-seed", spec.seed)
            .f64_bits("analytic-box-volume", spec.box_volume)
            .f64_bits("analytic-true-volume", spec.true_volume);
        let front = two_box_front(spec.dimensions);
        builder = builder.u64("front-points", usize_u64(front.len()));
        for (point_index, point) in front.iter().enumerate() {
            builder = builder
                .u64("front-point-index", usize_u64(point_index))
                .u64("front-point-dimensions", usize_u64(point.len()));
            for (coordinate_index, &coordinate) in point.iter().enumerate() {
                builder = builder
                    .u64("front-coordinate-index", usize_u64(coordinate_index))
                    .f64_bits("front-coordinate", coordinate);
            }
        }
        let reference = reference(spec.dimensions);
        builder = builder.u64("reference-dimensions", usize_u64(reference.len()));
        for (coordinate_index, coordinate) in reference.into_iter().enumerate() {
            builder = builder
                .u64("reference-coordinate-index", usize_u64(coordinate_index))
                .f64_bits("reference-coordinate", coordinate);
        }
    }
    builder.finish()
}

fn fixture_digest(fixture: &ReplayIdentity) -> ContentHash {
    hash_domain(FIXTURE_DIGEST_DOMAIN, fixture.canonical_bytes())
}

fn result_identity(
    fixture: &ReplayIdentity,
    strong_fixture: ContentHash,
    record: &StudyRecord,
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .child("fixture-compatibility-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .bytes("fixture-blake3", strong_fixture.as_bytes())
        .u64("estimate-count", usize_u64(record.estimates.len()));
    for (index, estimate) in record.estimates.iter().enumerate() {
        builder = builder
            .u64("estimate-index", usize_u64(index))
            .str("estimate-label", estimate.label)
            .u64("estimate-dimensions", usize_u64(estimate.dimensions))
            .f64_bits("returned-estimate", f64::from_bits(estimate.estimate))
            .u64("returned-hits", usize_u64(estimate.hits));
    }
    builder.finish()
}

fn result_digest(result: &ReplayIdentity) -> ContentHash {
    hash_domain(RESULT_DIGEST_DOMAIN, result.canonical_bytes())
}

fn event_digest(event: &Event) -> ContentHash {
    hash_domain(
        EVENT_DIGEST_DOMAIN,
        event.content_identity().canonical_bytes(),
    )
}

fn run_study() -> StudyRun {
    let estimates = fixture_specs()
        .into_iter()
        .map(|spec| {
            let front = two_box_front(spec.dimensions);
            let reference = reference(spec.dimensions);
            let (estimate, hits) = mc_hypervolume(&front, &reference, SAMPLES, spec.seed);
            EstimateBits {
                label: spec.label,
                dimensions: spec.dimensions,
                estimate: estimate.to_bits(),
                hits,
            }
        })
        .collect();
    let record = StudyRecord { estimates };
    let fixture = fixture_identity();
    let fixture_digest = fixture_digest(&fixture);
    let result = result_identity(&fixture, fixture_digest, &record);
    let result_digest = result_digest(&result);
    StudyRun {
        fixture,
        fixture_digest,
        record,
        result,
        result_digest,
    }
}

// Diagnostic-only randomized-QMC sanity scale. Its binomial-shaped term is a
// convenient magnitude check, not an IID standard error, QMC error bound, or
// confidence interval.
fn heuristic_analytic_sanity_band(spec: FixtureSpec) -> f64 {
    let probability = spec.true_volume / spec.box_volume;
    let binomial_scale_shape =
        spec.box_volume * fs_math::det::sqrt(probability * (1.0 - probability) / SAMPLES as f64);
    HEURISTIC_ANALYTIC_BAND_MULTIPLIER * binomial_scale_shape + spec.box_volume / SAMPLES as f64
}

fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    let specs = fixture_specs();
    if record.estimates.len() != specs.len() {
        return Some(format!(
            "estimate-count:{}!={}",
            record.estimates.len(),
            specs.len()
        ));
    }
    for (index, (retained, spec)) in record.estimates.iter().zip(specs).enumerate() {
        if retained.label != spec.label || retained.dimensions != spec.dimensions {
            return Some(format!(
                "estimate[{index}]-schema:label={};dimensions={};expected-label={};expected-dimensions={}",
                retained.label, retained.dimensions, spec.label, spec.dimensions
            ));
        }
        if retained.hits > SAMPLES {
            return Some(format!(
                "estimate[{index}]-hits:{}>{SAMPLES}",
                retained.hits
            ));
        }
        let estimate = f64::from_bits(retained.estimate);
        if !estimate.is_finite() || estimate < 0.0 || estimate > spec.box_volume {
            return Some(format!(
                "estimate[{index}]-range:{estimate:.17e};box={:.17e}",
                spec.box_volume
            ));
        }
        let accounted = spec.box_volume * retained.hits as f64 / SAMPLES as f64;
        if accounted.to_bits() != retained.estimate {
            return Some(format!(
                "estimate[{index}]-hit-accounting:reported=0x{:016x};accounted=0x{:016x};hits={}",
                retained.estimate,
                accounted.to_bits(),
                retained.hits
            ));
        }
        let sanity_band = heuristic_analytic_sanity_band(spec);
        if (estimate - spec.true_volume).abs() > sanity_band {
            return Some(format!(
                "estimate[{index}]-heuristic-analytic-sanity-band:reported={estimate:.17e};truth={:.17e};band={sanity_band:.17e};not-IID-standard-error;not-confidence-bound",
                spec.true_volume
            ));
        }
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let expected_fixture = fixture_identity();
    let expected_fixture_digest = fixture_digest(&expected_fixture);
    let computed_fixture_digest = fixture_digest(&run.fixture);
    if computed_fixture_digest != run.fixture_digest {
        return Err(AdmissionError::SuppliedFixtureDigestMismatch {
            declared: digest_bytes(run.fixture_digest),
            computed: digest_bytes(computed_fixture_digest),
        });
    }
    if run.fixture.canonical_bytes() != expected_fixture.canonical_bytes() {
        return Err(AdmissionError::RetainedFixtureIdentityMismatch {
            expected: digest_bytes(expected_fixture_digest),
            found: digest_bytes(computed_fixture_digest),
        });
    }
    let computed_result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    let computed_result_digest = result_digest(&computed_result);
    if run.result.canonical_bytes() != computed_result.canonical_bytes()
        || run.result_digest != computed_result_digest
    {
        return Err(AdmissionError::ResultPayloadIdentityMismatch {
            declared: digest_bytes(run.result_digest),
            computed: digest_bytes(computed_result_digest),
        });
    }
    Ok(())
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    match semantic_mismatch(&run.record) {
        Some(mismatch) => Err(AdmissionError::SemanticInconsistency(mismatch)),
        None => Ok(()),
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.result.canonical_bytes()
        && run.result_digest == reference.result_digest
    {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: digest_bytes(reference.result_digest),
            found: digest_bytes(run.result_digest),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    run.result_digest = result_digest(&run.result);
}

const fn admission_error_name(error: &AdmissionError) -> &'static str {
    match error {
        AdmissionError::SuppliedFixtureDigestMismatch { .. } => "SuppliedFixtureDigestMismatch",
        AdmissionError::RetainedFixtureIdentityMismatch { .. } => "RetainedFixtureIdentityMismatch",
        AdmissionError::ResultPayloadIdentityMismatch { .. } => "ResultPayloadIdentityMismatch",
        AdmissionError::ReferenceIdentityMismatch { .. } => "ReferenceIdentityMismatch",
        AdmissionError::SemanticInconsistency(_) => "SemanticInconsistency",
    }
}

fn fixture_tripwires(reference: &StudyRun) -> FixtureTripwires {
    let mut supplied_digest_run = reference.clone();
    supplied_digest_run.fixture_digest = hash_domain(
        SUPPLIED_FIXTURE_DIGEST_TRIPWIRE_DOMAIN,
        reference.fixture.canonical_bytes(),
    );
    assert_ne!(
        supplied_digest_run.fixture_digest, reference.fixture_digest,
        "fixture-digest tripwire must supply a distinct digest"
    );
    reseal(&mut supplied_digest_run);
    assert_eq!(
        supplied_digest_run.result,
        result_identity(
            &supplied_digest_run.fixture,
            supplied_digest_run.fixture_digest,
            &supplied_digest_run.record,
        ),
        "fixture-digest tripwire result frame must be self-consistent"
    );
    assert_eq!(
        supplied_digest_run.result_digest,
        result_digest(&supplied_digest_run.result),
        "fixture-digest tripwire result digest must be resealed"
    );
    let supplied_digest_error = validate_payload(&supplied_digest_run)
        .expect_err("self-consistent result with a false supplied fixture digest must refuse");

    let mut wrong_specs = fixture_specs();
    wrong_specs[0].seed = WRONG_FIXTURE_SOBOL_SEED;
    let mut wrong_fixture_run = reference.clone();
    wrong_fixture_run.fixture = fixture_identity_for(wrong_specs);
    wrong_fixture_run.fixture_digest = fixture_digest(&wrong_fixture_run.fixture);
    reseal(&mut wrong_fixture_run);
    assert_eq!(
        wrong_fixture_run.fixture_digest,
        fixture_digest(&wrong_fixture_run.fixture),
        "wrong fixture digest must be internally valid"
    );
    assert_eq!(
        wrong_fixture_run.result,
        result_identity(
            &wrong_fixture_run.fixture,
            wrong_fixture_run.fixture_digest,
            &wrong_fixture_run.record,
        ),
        "wrong fixture result frame must be self-consistent"
    );
    assert_eq!(
        wrong_fixture_run.result_digest,
        result_digest(&wrong_fixture_run.result),
        "wrong fixture result digest must be resealed"
    );
    let wrong_fixture_error = validate_payload(&wrong_fixture_run)
        .expect_err("self-consistent but unretained fixture identity must refuse");

    FixtureTripwires {
        supplied_digest_run,
        supplied_digest_error,
        wrong_fixture_run,
        wrong_fixture_error,
    }
}

fn exact_estimate_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(reference_row) = reference.record.estimates.get(mutation.estimate_index) else {
        return false;
    };
    let Some(mutant_row) = mutant.record.estimates.get(mutation.estimate_index) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference.fixture_digest != mutant.fixture_digest
        || reference_row.estimate != mutation.before
        || mutant_row.estimate != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.estimates[mutation.estimate_index].estimate = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let estimate_index =
        usize::try_from(selector.next_below(usize_u64(reference.record.estimates.len())))
            .expect("selected estimate index fits usize");
    let mantissa_bit = u32::try_from(selector.next_below(MUTATION_BIT_COUNT))
        .expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.estimates[estimate_index].estimate;
    let after = before ^ (1u64 << mantissa_bit);
    run.record.estimates[estimate_index].estimate = after;
    let stale_error = validate_payload(&run).expect_err("unsealed MC-HV mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed MC-HV mutation must miss retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed MC-HV estimate mutation must fail hit accounting");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            estimate_index,
            mantissa_bit,
            selector_draws,
            before,
            after,
        },
        stale_error,
        reference_error,
        semantic_error,
    }
}

fn green_receipt(run: &StudyRun, tripwires: &FixtureTripwires) -> Event {
    let sobol = &run.record.estimates[0];
    let hybrid = &run.record.estimates[1];
    let source_digest = production_source_digest();
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "mc-hypervolume-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"fixture_blake3\":\"{}\",",
                    "\"result_identity\":\"{}\",\"result_blake3\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::mc_hypervolume\",",
                    "\"production_source_blake3\":\"{}\",",
                    "\"samples_per_fixture\":{},\"max_sobol_dimensions\":{},",
                    "\"production_schedule\":{{\"sobol_point_index\":\"s+1\",",
                    "\"sobol_first_index\":{},\"philox_tail_kernel\":\"0x{:04x}\",",
                    "\"philox_tail_tile\":{},\"philox_tail_constructor\":\"fs_rand::StreamKey\",",
                    "\"hybrid_tail_draws_per_sample\":{}}},",
                    "\"analytic_sanity_band\":{{\"classification\":",
                    "\"heuristic-only;randomized-QMC-not-IID;not-standard-error-or-confidence-bound\",",
                    "\"multiplier_bits\":\"0x{:016x}\",\"sobol_band_bits\":\"0x{:016x}\",",
                    "\"hybrid_band_bits\":\"0x{:016x}\"}},",
                    "\"fixture_tripwires\":{{\"supplied_digest_refusal\":\"{}\",",
                    "\"supplied_digest\":\"{}\",\"resealed_wrong_fixture_refusal\":\"{}\",",
                    "\"wrong_fixture_sobol_seed\":{},\"wrong_fixture_blake3\":\"{}\"}},",
                    "\"sobol\":{{\"dimensions\":{},\"seed\":{},",
                    "\"estimate_bits\":\"0x{:016x}\",\"hits\":{},",
                    "\"true_volume_bits\":\"0x{:016x}\"}},",
                    "\"hybrid\":{{\"dimensions\":{},\"seed\":{},",
                    "\"estimate_bits\":\"0x{:016x}\",\"hits\":{},",
                    "\"true_volume_bits\":\"0x{:016x}\"}},",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\",",
                    "\"stream_semantics\":{}}},",
                    "\"no_claims\":[\"exact-MC\",\"certified-MC\",",
                    "\"IID-standard-error\",\"QMC-error-bound\",",
                    "\"formal-confidence\",\"independent-sampler\",",
                    "\"independent-schedule-verifier\",",
                    "\"sample-trace\",\"arbitrary-inputs\",\"cross-process\",",
                    "\"cross-ISA\",\"cancellation\",\"persistence\",",
                    "\"authentication\",\"performance\",\"optimizer-quality\",",
                    "\"archive\"]}}"
                ),
                run.fixture.hex(),
                run.fixture_digest.to_hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                source_digest.to_hex(),
                SAMPLES,
                fs_rand::qmc::MAX_SOBOL_DIM,
                PRODUCTION_SOBOL_FIRST_POINT_INDEX,
                PRODUCTION_PHILOX_TAIL_KERNEL,
                PRODUCTION_PHILOX_TAIL_TILE,
                HYBRID_DIMENSIONS - SOBOL_DIMENSIONS,
                HEURISTIC_ANALYTIC_BAND_MULTIPLIER.to_bits(),
                heuristic_analytic_sanity_band(fixture_specs()[0]).to_bits(),
                heuristic_analytic_sanity_band(fixture_specs()[1]).to_bits(),
                admission_error_name(&tripwires.supplied_digest_error),
                tripwires.supplied_digest_run.fixture_digest.to_hex(),
                admission_error_name(&tripwires.wrong_fixture_error),
                WRONG_FIXTURE_SOBOL_SEED,
                tripwires.wrong_fixture_run.fixture_digest.to_hex(),
                sobol.dimensions,
                SOBOL_SEED,
                sobol.estimate,
                sobol.hits,
                SOBOL_TRUE_VOLUME.to_bits(),
                hybrid.dimensions,
                HYBRID_SEED,
                hybrid.estimate,
                hybrid.hits,
                HYBRID_TRUE_VOLUME.to_bits(),
                fs_dfo::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
                fs_rand::STREAM_SEMANTICS_VERSION,
            ),
        },
        None,
    )
}

fn green_verdict(run: &StudyRun, tripwires: &FixtureTripwires) -> Event {
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: format!(
                "fixture={}; result={}; blake3={}; source-blake3={}; sobol-seed=0x{SOBOL_SEED:016x}; hybrid-seed=0x{HYBRID_SEED:016x}; schedule=Sobol::point(s+1),StreamKey(seed,0x{PRODUCTION_PHILOX_TAIL_KERNEL:04x},{PRODUCTION_PHILOX_TAIL_TILE}).next_f64; analytic-band=heuristic-QMC-sanity-only-not-IID-standard-error-or-confidence; supplied-digest-tripwire={}; resealed-wrong-fixture-tripwire={}@sobol-seed-0x{WRONG_FIXTURE_SOBOL_SEED:016x}; estimates={}; composite aggregate seed zero",
                run.fixture.hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                production_source_digest().to_hex(),
                admission_error_name(&tripwires.supplied_digest_error),
                admission_error_name(&tripwires.wrong_fixture_error),
                run.record.estimates.len(),
            ),
            seed: 0,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let mutation = corruption.mutation;
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=estimates[{}].estimate; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result_digest.to_hex(),
        corruption.run.result_digest.to_hex(),
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.estimate_index,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
        },
        None,
    )
}

fn assert_mergeable(event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn assert_event_pair(first: &Event, second: &Event, label: &str) {
    assert_eq!(
        first.content_identity().canonical_bytes(),
        second.content_identity().canonical_bytes(),
        "{label} content must replay byte-for-byte"
    );
    assert_eq!(event_digest(first), event_digest(second));
    for event in [first, second] {
        fs_obs::lint_failure_record(event).expect("MC-HV evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("MC-HV evidence is fs-obs wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("MC-HV evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay plus every refusal gate.
fn hybrid_mc_hypervolume_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();

    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(validate_semantics(&original), Ok(()));
    assert_eq!(validate_semantics(&replay), Ok(()));
    assert_eq!(admit_reference(&original, &replay), Ok(()));
    assert_eq!(admit_reference(&replay, &original), Ok(()));
    assert_eq!(original.record, replay.record);
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.fixture_digest, replay.fixture_digest);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.result_digest, replay.result_digest);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete MC-HV result frames must replay byte-for-byte"
    );

    let first_tripwires = fixture_tripwires(&original);
    let second_tripwires = fixture_tripwires(&replay);
    assert_eq!(
        first_tripwires, second_tripwires,
        "fixture-integrity and retained-identity refusals must replay exactly"
    );
    assert_eq!(
        validate_semantics(&first_tripwires.supplied_digest_run),
        Ok(()),
        "supplied-digest tripwire changes no scientific result field"
    );
    assert_eq!(
        first_tripwires.supplied_digest_run.fixture, original.fixture,
        "supplied-digest tripwire retains the exact expected fixture"
    );
    assert!(matches!(
        &first_tripwires.supplied_digest_error,
        AdmissionError::SuppliedFixtureDigestMismatch { declared, computed }
            if declared == first_tripwires.supplied_digest_run.fixture_digest.as_bytes()
                && computed == original.fixture_digest.as_bytes()
    ));
    assert_eq!(
        validate_semantics(&first_tripwires.wrong_fixture_run),
        Ok(()),
        "wrong-fixture tripwire retains a numerically self-consistent record"
    );
    assert_ne!(
        first_tripwires.wrong_fixture_run.fixture, original.fixture,
        "wrong-fixture tripwire must change the canonical fixture"
    );
    assert!(matches!(
        &first_tripwires.wrong_fixture_error,
        AdmissionError::RetainedFixtureIdentityMismatch { expected, found }
            if expected == original.fixture_digest.as_bytes()
                && found == first_tripwires.wrong_fixture_run.fixture_digest.as_bytes()
    ));

    let first_receipt = green_receipt(&original, &first_tripwires);
    let second_receipt = green_receipt(&replay, &second_tripwires);
    assert_event_pair(&first_receipt, &second_receipt, "green MC-HV receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original, &first_tripwires);
    let second_green = green_verdict(&replay, &second_tripwires);
    assert_event_pair(&first_green, &second_green, "green MC-HV verdict");
    assert_mergeable(&first_green);
    assert_mergeable(&second_green);
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded MC-HV corruption must replay exactly");
    assert!(
        exact_estimate_bit_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one returned estimate bit"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed MC-HV mutation must be internally self-consistent"
    );
    let after = f64::from_bits(first.mutation.after);
    assert!(after.is_finite() && after > 0.0);
    assert!(matches!(
        &first.stale_error,
        AdmissionError::ResultPayloadIdentityMismatch { declared, computed }
            if declared == original.result_digest.as_bytes()
                && computed == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result_digest.as_bytes()
                && found == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.semantic_error,
        AdmissionError::SemanticInconsistency(mismatch)
            if mismatch.contains("hit-accounting")
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_pair(&first_red, &second_red, "red MC-HV evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded MC-HV corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!(
        "estimates[{}].estimate",
        first.mutation.estimate_index
    )));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("ResultPayloadIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
