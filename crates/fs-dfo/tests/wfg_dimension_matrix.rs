//! Cross-objective WFG1-WFG9 dimension matrix (7tv.24.8).
//!
//! This crate-external G0/G3 battery exercises the public normalized family at
//! minimum and wider valid partitions from two through eight objectives. It
//! binds trace dimensions, fixed objective scaling, WFG3 degeneracy semantics,
//! variant distinction, and deterministic same-process repeatability.
//!
//! It does not claim optimizer convergence, external executable parity,
//! cancellation, cross-ISA bit stability, or performance.

#![deny(unsafe_code)]

use fs_dfo::wfg::{Wfg1, Wfg2, Wfg3, Wfg4, Wfg5, Wfg6, Wfg7, Wfg8, Wfg9, WfgError, WfgEvaluation};

#[derive(Debug, Clone, Copy)]
struct MatrixCase {
    objectives: usize,
    position_parameters: usize,
    distance_parameters: usize,
}

impl MatrixCase {
    const fn dimension(self) -> usize {
        self.position_parameters + self.distance_parameters
    }
}

const CASES: [MatrixCase; 8] = [
    MatrixCase {
        objectives: 2,
        position_parameters: 1,
        distance_parameters: 2,
    },
    MatrixCase {
        objectives: 2,
        position_parameters: 4,
        distance_parameters: 8,
    },
    MatrixCase {
        objectives: 3,
        position_parameters: 2,
        distance_parameters: 2,
    },
    MatrixCase {
        objectives: 3,
        position_parameters: 6,
        distance_parameters: 8,
    },
    MatrixCase {
        objectives: 5,
        position_parameters: 4,
        distance_parameters: 2,
    },
    MatrixCase {
        objectives: 5,
        position_parameters: 12,
        distance_parameters: 8,
    },
    MatrixCase {
        objectives: 8,
        position_parameters: 7,
        distance_parameters: 2,
    },
    MatrixCase {
        objectives: 8,
        position_parameters: 21,
        distance_parameters: 8,
    },
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
    fn family(case: MatrixCase) -> Result<[Self; 9], WfgError> {
        let MatrixCase {
            objectives,
            position_parameters,
            distance_parameters,
        } = case;
        Ok([
            Self::Wfg1(Wfg1::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg2(Wfg2::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg3(Wfg3::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg4(Wfg4::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg5(Wfg5::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg6(Wfg6::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg7(Wfg7::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg8(Wfg8::new(
                objectives,
                position_parameters,
                distance_parameters,
            )?),
            Self::Wfg9(Wfg9::new(
                objectives,
                position_parameters,
                distance_parameters,
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

    const fn objectives(self) -> usize {
        match self {
            Self::Wfg1(problem) => problem.objectives(),
            Self::Wfg2(problem) => problem.objectives(),
            Self::Wfg3(problem) => problem.objectives(),
            Self::Wfg4(problem) => problem.objectives(),
            Self::Wfg5(problem) => problem.objectives(),
            Self::Wfg6(problem) => problem.objectives(),
            Self::Wfg7(problem) => problem.objectives(),
            Self::Wfg8(problem) => problem.objectives(),
            Self::Wfg9(problem) => problem.objectives(),
        }
    }

    const fn position_parameters(self) -> usize {
        match self {
            Self::Wfg1(problem) => problem.position_parameters(),
            Self::Wfg2(problem) => problem.position_parameters(),
            Self::Wfg3(problem) => problem.position_parameters(),
            Self::Wfg4(problem) => problem.position_parameters(),
            Self::Wfg5(problem) => problem.position_parameters(),
            Self::Wfg6(problem) => problem.position_parameters(),
            Self::Wfg7(problem) => problem.position_parameters(),
            Self::Wfg8(problem) => problem.position_parameters(),
            Self::Wfg9(problem) => problem.position_parameters(),
        }
    }

    const fn distance_parameters(self) -> usize {
        match self {
            Self::Wfg1(problem) => problem.distance_parameters(),
            Self::Wfg2(problem) => problem.distance_parameters(),
            Self::Wfg3(problem) => problem.distance_parameters(),
            Self::Wfg4(problem) => problem.distance_parameters(),
            Self::Wfg5(problem) => problem.distance_parameters(),
            Self::Wfg6(problem) => problem.distance_parameters(),
            Self::Wfg7(problem) => problem.distance_parameters(),
            Self::Wfg8(problem) => problem.distance_parameters(),
            Self::Wfg9(problem) => problem.distance_parameters(),
        }
    }

    const fn dimension(self) -> usize {
        match self {
            Self::Wfg1(problem) => problem.dimension(),
            Self::Wfg2(problem) => problem.dimension(),
            Self::Wfg3(problem) => problem.dimension(),
            Self::Wfg4(problem) => problem.dimension(),
            Self::Wfg5(problem) => problem.dimension(),
            Self::Wfg6(problem) => problem.dimension(),
            Self::Wfg7(problem) => problem.dimension(),
            Self::Wfg8(problem) => problem.dimension(),
            Self::Wfg9(problem) => problem.dimension(),
        }
    }

    const fn transformed_dimension(self) -> usize {
        match self {
            Self::Wfg2(problem) => {
                problem.position_parameters() + problem.distance_parameters() / 2
            }
            Self::Wfg3(problem) => {
                problem.position_parameters() + problem.distance_parameters() / 2
            }
            _ => self.dimension(),
        }
    }

    const fn is_wfg3(self) -> bool {
        matches!(self, Self::Wfg3(_))
    }

    fn evaluate(self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
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
}

fn normalized_probe(dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|component| ((component * 7 + 3) % 31 + 1) as f64 / 32.0)
        .collect()
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

fn assert_repeats(first: &WfgEvaluation, second: &WfgEvaluation, context: &str) {
    assert_slice_bits_eq(first.transformed(), second.transformed(), context);
    assert_slice_bits_eq(first.reduced(), second.reduced(), context);
    assert_slice_bits_eq(first.positioned(), second.positioned(), context);
    assert_slice_bits_eq(first.shape(), second.shape(), context);
    assert_slice_bits_eq(first.objectives(), second.objectives(), context);
}

fn assert_unit_trace(values: &[f64], name: &str, trace: &str) {
    for (component, &value) in values.iter().enumerate() {
        assert!(
            value.is_finite() && (0.0..=1.0).contains(&value),
            "{name} {trace}[{component}] escaped [0,1]: {value:.17e}"
        );
    }
}

fn assert_positioning(problem: PublicWfg, evaluation: &WfgEvaluation) {
    let name = problem.name();
    if !problem.is_wfg3() {
        assert_slice_bits_eq(evaluation.positioned(), evaluation.reduced(), name);
        return;
    }

    let last = problem.objectives() - 1;
    let distance = evaluation.reduced()[last];
    let mut changed = false;
    for (axis, (&reduced, &positioned)) in evaluation
        .reduced()
        .iter()
        .zip(evaluation.positioned())
        .enumerate()
    {
        let expected = if axis == 0 || axis == last {
            reduced
        } else {
            distance.mul_add(reduced - 0.5, 0.5)
        };
        assert_eq!(
            positioned.to_bits(),
            expected.to_bits(),
            "WFG3 positioned axis {axis}"
        );
        changed |= positioned.to_bits() != reduced.to_bits();
    }
    assert_eq!(changed, problem.objectives() > 2, "WFG3 degenerate axes");
}

fn exercise_problem(problem: PublicWfg, case: MatrixCase, input: &[f64]) -> Vec<f64> {
    let name = problem.name();
    assert_eq!(problem.objectives(), case.objectives, "{name} objectives");
    assert_eq!(
        problem.position_parameters(),
        case.position_parameters,
        "{name} position parameters"
    );
    assert_eq!(
        problem.distance_parameters(),
        case.distance_parameters,
        "{name} distance parameters"
    );
    assert_eq!(problem.dimension(), case.dimension(), "{name} dimension");

    let first = problem.evaluate(input).unwrap();
    let second = problem.evaluate(input).unwrap();
    assert_eq!(
        first.transformed().len(),
        problem.transformed_dimension(),
        "{name} transformed dimension"
    );
    assert_eq!(first.reduced().len(), case.objectives, "{name} reductions");
    assert_eq!(
        first.positioned().len(),
        case.objectives,
        "{name} positions"
    );
    assert_eq!(first.shape().len(), case.objectives, "{name} shapes");
    assert_eq!(first.objectives().len(), case.objectives, "{name} outputs");

    for (trace, values) in [
        ("transformed", first.transformed()),
        ("reduced", first.reduced()),
        ("positioned", first.positioned()),
        ("shape", first.shape()),
    ] {
        assert_unit_trace(values, name, trace);
    }
    assert_repeats(&first, &second, name);

    let distance = first.positioned()[case.objectives - 1];
    for (objective, (&shape, &actual)) in first.shape().iter().zip(first.objectives()).enumerate() {
        let scale = 2.0 * (objective + 1) as f64;
        assert_eq!(
            actual.to_bits(),
            scale.mul_add(shape, distance).to_bits(),
            "{name} objective {objective} scale reconstruction"
        );
    }

    assert_positioning(problem, &first);
    first.into_objectives()
}

#[test]
fn public_family_spans_minimum_and_wide_objective_dimension_cases() {
    for case in CASES {
        let input = normalized_probe(case.dimension());
        let mut observed = Vec::with_capacity(9);

        for problem in PublicWfg::family(case).unwrap() {
            observed.push((problem.name(), exercise_problem(problem, case, &input)));
        }

        for (left_index, (left_name, left)) in observed.iter().enumerate() {
            for (right_name, right) in observed.iter().skip(left_index + 1) {
                assert!(
                    left.iter()
                        .zip(right)
                        .any(|(&left, &right)| left.to_bits() != right.to_bits()),
                    "{left_name} and {right_name} aliased for {case:?}"
                );
            }
        }
    }
}
