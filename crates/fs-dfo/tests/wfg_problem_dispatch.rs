//! Typed public WFG family-dispatch conformance (7tv.24.9).
//!
//! This crate-external G0/G3 battery binds the stable WFG1-WFG9 catalog and
//! proves that dynamic study selection delegates bitwise to every concrete
//! typed evaluator in both normalized and canonical input domains.
//!
//! It does not claim optimizer quality, adaptive selection, cancellation,
//! external executable parity, cross-ISA bit stability, or performance.

#![deny(unsafe_code)]

use fs_dfo::wfg::{
    Wfg1, Wfg2, Wfg3, Wfg4, Wfg5, Wfg6, Wfg7, Wfg8, Wfg9, WfgError, WfgEvaluation, WfgProblem,
    WfgVariant,
};

#[derive(Debug, Clone, Copy)]
struct DispatchCase {
    objectives: usize,
    position_parameters: usize,
    distance_parameters: usize,
}

impl DispatchCase {
    const fn dimension(self) -> usize {
        self.position_parameters + self.distance_parameters
    }
}

const CASES: [DispatchCase; 3] = [
    DispatchCase {
        objectives: 2,
        position_parameters: 1,
        distance_parameters: 2,
    },
    DispatchCase {
        objectives: 4,
        position_parameters: 6,
        distance_parameters: 6,
    },
    DispatchCase {
        objectives: 5,
        position_parameters: 12,
        distance_parameters: 8,
    },
];

fn normalized_probe(dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|component| ((component * 11 + 5) % 31 + 1) as f64 / 32.0)
        .collect()
}

fn canonical_coordinates(normalized: &[f64]) -> Vec<f64> {
    normalized
        .iter()
        .enumerate()
        .map(|(component, &value)| value * 2.0 * (component as f64 + 1.0))
        .collect()
}

fn direct_evaluations(
    variant: WfgVariant,
    case: DispatchCase,
    normalized: &[f64],
    canonical: &[f64],
) -> Result<(WfgEvaluation, WfgEvaluation), WfgError> {
    macro_rules! evaluate {
        ($problem:expr) => {{
            let problem = $problem?;
            Ok((
                problem.evaluate_normalized(normalized)?,
                problem.evaluate_canonical(canonical)?,
            ))
        }};
    }

    let DispatchCase {
        objectives,
        position_parameters,
        distance_parameters,
    } = case;
    match variant {
        WfgVariant::Wfg1 => evaluate!(Wfg1::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg2 => evaluate!(Wfg2::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg3 => evaluate!(Wfg3::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg4 => evaluate!(Wfg4::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg5 => evaluate!(Wfg5::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg6 => evaluate!(Wfg6::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg7 => evaluate!(Wfg7::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg8 => evaluate!(Wfg8::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
        WfgVariant::Wfg9 => evaluate!(Wfg9::new(
            objectives,
            position_parameters,
            distance_parameters
        )),
    }
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
fn public_variant_catalog_has_stable_complete_family_order() {
    let expected = [
        "WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9",
    ];
    assert_eq!(WfgVariant::ALL.map(WfgVariant::name), expected);
    for (variant, expected) in WfgVariant::ALL.into_iter().zip(expected) {
        assert_eq!(variant.to_string(), expected);
    }
}

#[test]
fn public_dispatch_matches_every_typed_evaluator_in_both_domains() {
    for case in CASES {
        let normalized = normalized_probe(case.dimension());
        let canonical = canonical_coordinates(&normalized);

        for variant in WfgVariant::ALL {
            let problem = WfgProblem::new(
                variant,
                case.objectives,
                case.position_parameters,
                case.distance_parameters,
            )
            .unwrap();
            assert_eq!(problem.variant(), variant);
            assert_eq!(problem.objectives(), case.objectives);
            assert_eq!(problem.position_parameters(), case.position_parameters);
            assert_eq!(problem.distance_parameters(), case.distance_parameters);
            assert_eq!(problem.dimension(), case.dimension());

            let dispatched_normalized = problem.evaluate_normalized(&normalized).unwrap();
            let dispatched_canonical = problem.evaluate_canonical(&canonical).unwrap();
            let (direct_normalized, direct_canonical) =
                direct_evaluations(variant, case, &normalized, &canonical).unwrap();
            assert_evaluation_bits_eq(&dispatched_normalized, &direct_normalized, variant.name());
            assert_evaluation_bits_eq(&dispatched_canonical, &direct_canonical, variant.name());
            assert_evaluation_bits_eq(
                &dispatched_canonical,
                &dispatched_normalized,
                variant.name(),
            );
        }
    }
}

#[test]
fn public_dispatch_preserves_selected_variant_refusals() {
    for variant in WfgVariant::ALL {
        assert_eq!(
            WfgProblem::new(variant, 1, 1, 2).unwrap_err(),
            WfgError::TooFewObjectives { objectives: 1 },
            "{variant} constructor refusal"
        );

        let problem = WfgProblem::new(variant, 4, 6, 6).unwrap();
        let mut normalized = [0.0; 12];
        normalized[2] = f64::from_bits(0x7ff8_1234_5678_9abc);
        assert_eq!(
            problem.evaluate_normalized(&normalized).unwrap_err(),
            WfgError::NonFiniteInput {
                component: 2,
                bits: normalized[2].to_bits(),
            },
            "{variant} normalized refusal"
        );

        let mut canonical = [0.0; 12];
        let upper_bound = 6.0_f64;
        canonical[2] = f64::from_bits(upper_bound.to_bits() + 1);
        assert_eq!(
            problem.evaluate_canonical(&canonical).unwrap_err(),
            WfgError::CanonicalInputOutOfRange {
                component: 2,
                bits: canonical[2].to_bits(),
                upper_bound_bits: upper_bound.to_bits(),
            },
            "{variant} canonical refusal"
        );
    }

    for variant in [WfgVariant::Wfg2, WfgVariant::Wfg3] {
        assert_eq!(
            WfgProblem::new(variant, 4, 6, 3).unwrap_err(),
            WfgError::DistanceParametersNotEven {
                distance_parameters: 3,
            },
            "{variant} odd-distance refusal"
        );
    }
}
