//! Canonical heterogeneous-bound WFG adapter matrix (7tv.24.7).
//!
//! The production WFG family admits both normalized coordinates and the
//! standard canonical domain `z_i in [0, 2(i + 1)]`. This crate-external G0/G3
//! battery binds their equivalence, indexed endpoint admission, structured
//! refusal, and deterministic same-process repeatability through public APIs.
//!
//! This does not claim external executable parity, optimizer convergence,
//! cancellation, cross-ISA bit stability, or performance.

#![deny(unsafe_code)]

use fs_dfo::wfg::{Wfg1, Wfg2, Wfg3, Wfg4, Wfg5, Wfg6, Wfg7, Wfg8, Wfg9, WfgError, WfgEvaluation};

const OBJECTIVES: usize = 4;
const POSITION_PARAMETERS: usize = 6;
const DISTANCE_PARAMETERS: usize = 6;
const DIMENSION: usize = POSITION_PARAMETERS + DISTANCE_PARAMETERS;

const ASYMMETRIC_NORMALIZED: [f64; DIMENSION] = [
    0.031_25, 0.093_75, 0.156_25, 0.218_75, 0.281_25, 0.343_75, 0.406_25, 0.468_75, 0.531_25,
    0.593_75, 0.718_75, 0.906_25,
];

#[derive(Debug, Clone, Copy)]
#[allow(clippy::enum_variant_names)] // Exact public names make matrix failures actionable.
enum PublicWfg {
    Wfg1(Wfg1),
    Wfg2(Wfg2),
    Wfg3(Wfg3),
    Wfg4(Wfg4),
    Wfg5(Wfg5),
    Wfg6(Wfg6),
    Wfg7(Wfg7),
    Wfg8(Wfg8),
    Wfg9(Wfg9),
}

impl PublicWfg {
    fn family() -> Result<[Self; 9], WfgError> {
        Ok([
            Self::Wfg1(Wfg1::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg2(Wfg2::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg3(Wfg3::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg4(Wfg4::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg5(Wfg5::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg6(Wfg6::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg7(Wfg7::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg8(Wfg8::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
            Self::Wfg9(Wfg9::new(
                OBJECTIVES,
                POSITION_PARAMETERS,
                DISTANCE_PARAMETERS,
            )?),
        ])
    }

    const fn name(self) -> &'static str {
        match self {
            Self::Wfg1(_) => "WFG1",
            Self::Wfg2(_) => "WFG2",
            Self::Wfg3(_) => "WFG3",
            Self::Wfg4(_) => "WFG4",
            Self::Wfg5(_) => "WFG5",
            Self::Wfg6(_) => "WFG6",
            Self::Wfg7(_) => "WFG7",
            Self::Wfg8(_) => "WFG8",
            Self::Wfg9(_) => "WFG9",
        }
    }

    fn evaluate_normalized(self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        match self {
            Self::Wfg1(problem) => problem.evaluate_normalized(input),
            Self::Wfg2(problem) => problem.evaluate_normalized(input),
            Self::Wfg3(problem) => problem.evaluate_normalized(input),
            Self::Wfg4(problem) => problem.evaluate_normalized(input),
            Self::Wfg5(problem) => problem.evaluate_normalized(input),
            Self::Wfg6(problem) => problem.evaluate_normalized(input),
            Self::Wfg7(problem) => problem.evaluate_normalized(input),
            Self::Wfg8(problem) => problem.evaluate_normalized(input),
            Self::Wfg9(problem) => problem.evaluate_normalized(input),
        }
    }

    fn evaluate_canonical(self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        match self {
            Self::Wfg1(problem) => problem.evaluate_canonical(input),
            Self::Wfg2(problem) => problem.evaluate_canonical(input),
            Self::Wfg3(problem) => problem.evaluate_canonical(input),
            Self::Wfg4(problem) => problem.evaluate_canonical(input),
            Self::Wfg5(problem) => problem.evaluate_canonical(input),
            Self::Wfg6(problem) => problem.evaluate_canonical(input),
            Self::Wfg7(problem) => problem.evaluate_canonical(input),
            Self::Wfg8(problem) => problem.evaluate_canonical(input),
            Self::Wfg9(problem) => problem.evaluate_canonical(input),
        }
    }
}

fn canonical_upper_bound(component: usize) -> f64 {
    2.0 * (component as f64 + 1.0)
}

fn canonical_coordinates(normalized: &[f64; DIMENSION]) -> [f64; DIMENSION] {
    let mut canonical = [0.0; DIMENSION];
    for (component, (&value, output)) in normalized.iter().zip(&mut canonical).enumerate() {
        *output = value * canonical_upper_bound(component);
    }
    canonical
}

fn assert_slice_bits_eq(actual: &[f64], expected: &[f64], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: length");
    for (component, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}[{component}]: actual={actual:.17e}, expected={expected:.17e}"
        );
    }
}

fn assert_evaluation_bits_eq(actual: &WfgEvaluation, expected: &WfgEvaluation, context: &str) {
    assert_slice_bits_eq(actual.transformed(), expected.transformed(), context);
    assert_slice_bits_eq(actual.reduced(), expected.reduced(), context);
    assert_slice_bits_eq(actual.positioned(), expected.positioned(), context);
    assert_slice_bits_eq(actual.shape(), expected.shape(), context);
    assert_slice_bits_eq(actual.objectives(), expected.objectives(), context);
}

#[test]
fn canonical_family_matches_normalized_kernels_and_repeats_bitwise() {
    let probes = [
        [0.0; DIMENSION],
        [0.25; DIMENSION],
        ASYMMETRIC_NORMALIZED,
        [1.0; DIMENSION],
    ];

    for problem in PublicWfg::family().unwrap() {
        for normalized in &probes {
            let name = problem.name();
            let canonical = canonical_coordinates(normalized);
            let expected = problem.evaluate_normalized(normalized).unwrap();
            let first = problem.evaluate_canonical(&canonical).unwrap();
            let second = problem.evaluate_canonical(&canonical).unwrap();

            assert_evaluation_bits_eq(&first, &expected, name);
            assert_evaluation_bits_eq(&first, &second, name);
        }
    }
}

#[test]
fn canonical_bounds_are_inclusive_and_index_specific() {
    let problem = Wfg9::new(OBJECTIVES, POSITION_PARAMETERS, DISTANCE_PARAMETERS).unwrap();

    for component in 0..DIMENSION {
        let upper_bound = canonical_upper_bound(component);
        let mut canonical = [0.0; DIMENSION];
        canonical[component] = upper_bound;
        let actual = problem.evaluate_canonical(&canonical).unwrap();

        let mut normalized = [0.0; DIMENSION];
        normalized[component] = 1.0;
        let expected = problem.evaluate_normalized(&normalized).unwrap();
        assert_evaluation_bits_eq(&actual, &expected, "inclusive indexed upper bound");

        let above = f64::from_bits(upper_bound.to_bits() + 1);
        canonical[component] = above;
        assert_eq!(
            problem.evaluate_canonical(&canonical).unwrap_err(),
            WfgError::CanonicalInputOutOfRange {
                component,
                bits: above.to_bits(),
                upper_bound_bits: upper_bound.to_bits(),
            }
        );

        canonical[component] = -f64::EPSILON;
        assert_eq!(
            problem.evaluate_canonical(&canonical).unwrap_err(),
            WfgError::CanonicalInputOutOfRange {
                component,
                bits: (-f64::EPSILON).to_bits(),
                upper_bound_bits: upper_bound.to_bits(),
            }
        );
    }

    let mut canonical = [0.0; DIMENSION];
    canonical[5] = -0.0;
    let actual = problem.evaluate_canonical(&canonical).unwrap();
    let mut normalized = [0.0; DIMENSION];
    normalized[5] = -0.0;
    let expected = problem.evaluate_normalized(&normalized).unwrap();
    assert_evaluation_bits_eq(&actual, &expected, "canonical negative zero");
}

#[test]
fn canonical_family_refuses_shape_and_nonfinite_inputs_uniformly() {
    for problem in PublicWfg::family().unwrap() {
        let name = problem.name();
        let mut malformed_short = [0.0; DIMENSION - 1];
        malformed_short[0] = f64::NAN;
        assert_eq!(
            problem.evaluate_canonical(&malformed_short).unwrap_err(),
            WfgError::WrongInputLength {
                expected: DIMENSION,
                actual: DIMENSION - 1,
            },
            "{name} short canonical input"
        );
        assert_eq!(
            problem
                .evaluate_canonical(&[0.0; DIMENSION + 1])
                .unwrap_err(),
            WfgError::WrongInputLength {
                expected: DIMENSION,
                actual: DIMENSION + 1,
            },
            "{name} long canonical input"
        );

        let payload_nan = f64::from_bits(0x7ff8_1234_5678_9abc);
        for value in [payload_nan, f64::INFINITY, f64::NEG_INFINITY] {
            let mut canonical = [0.0; DIMENSION];
            canonical[5] = value;
            assert_eq!(
                problem.evaluate_canonical(&canonical).unwrap_err(),
                WfgError::NonFiniteCanonicalInput {
                    component: 5,
                    bits: value.to_bits(),
                },
                "{name} non-finite canonical input"
            );
        }
    }
}
