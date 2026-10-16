//! Public WFG1-WFG9 family conformance matrix (7tv.24.6).
//!
//! This crate-external battery exercises the typed normalized evaluators only
//! through `fs_dfo::wfg`. It binds their shared admission and trace contract,
//! a frozen asymmetric family probe, boundary behavior, and deterministic
//! replay without duplicating the retained WFG4 optimizer campaign. The
//! frozen values come from an independent direct f64 port of the corrected
//! WFG equations at jMetal revision
//! `ea7e882f6b8f94b99535921674e62cda7986f20e`.
//!
//! No claim is made here about heterogeneous canonical-bound adaptation,
//! optimizer convergence, external executable parity, cancellation,
//! cross-ISA bit stability, or performance.

#![deny(unsafe_code)]

use fs_dfo::wfg::{Wfg1, Wfg2, Wfg3, Wfg4, Wfg5, Wfg6, Wfg7, Wfg8, Wfg9, WfgError, WfgEvaluation};

const OBJECTIVES: usize = 4;
const POSITION_PARAMETERS: usize = 6;
const DISTANCE_PARAMETERS: usize = 6;
const DIMENSION: usize = POSITION_PARAMETERS + DISTANCE_PARAMETERS;
const COMPRESSED_DIMENSION: usize = POSITION_PARAMETERS + DISTANCE_PARAMETERS / 2;
const TOLERANCE: f64 = 1.0e-10;

const MATRIX_INPUT: [f64; DIMENSION] = [
    0.07, 0.16, 0.29, 0.41, 0.58, 0.73, 0.09, 0.24, 0.39, 0.57, 0.76, 0.93,
];

const EXPECTED_OBJECTIVES: [[f64; OBJECTIVES]; 9] = [
    [
        2.768_600_224_231_849,
        0.982_622_086_858_670_3,
        0.985_242_546_246_867_8,
        1.067_266_014_389_427_3,
    ],
    [
        0.548_353_791_428_739_2,
        0.547_405_966_338_882,
        0.592_649_364_215_561_7,
        8.495_894_747_158_395,
    ],
    [
        0.602_251_353_590_325,
        0.625_916_340_438_397_5,
        0.947_546_031_746_031_8,
        7.626_031_746_031_746,
    ],
    [
        0.436_597_272_751_891_1,
        0.705_029_059_823_277_9,
        4.869_552_170_575_179,
        5.606_144_039_584_852,
    ],
    [
        1.340_066_301_510_183_5,
        1.892_466_234_493_387_6,
        1.153_408_558_093_649_3,
        7.263_161_918_167_952,
    ],
    [
        0.838_132_085_900_593_7,
        0.955_676_195_184_115_6,
        1.814_273_919_577_827_2,
        8.504_172_026_146_684,
    ],
    [
        0.503_492_830_836_990_5,
        0.513_972_095_359_592,
        1.635_336_019_191_497_8,
        8.351_205_718_756_75,
    ],
    [
        0.667_775_849_823_971,
        0.700_603_860_343_086_6,
        1.426_047_269_151_465_8,
        8.376_760_479_496_88,
    ],
    [
        0.439_037_889_780_964_75,
        0.603_439_658_694_161_6,
        3.846_968_585_159_4,
        6.882_147_066_421_952,
    ],
];

#[derive(Debug, Clone, Copy)]
#[allow(clippy::enum_variant_names)] // Exact public type names keep failure output unambiguous.
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
            Self::Wfg2(_) | Self::Wfg3(_) => COMPRESSED_DIMENSION,
            _ => DIMENSION,
        }
    }

    const fn has_degenerate_positioning(self) -> bool {
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

fn assert_close(actual: f64, expected: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= TOLERANCE,
        "{context}: actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn assert_slice_close(actual: &[f64], expected: &[f64], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: length");
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert_close(actual, expected, &format!("{context}[{index}]"));
    }
}

fn assert_slice_bits_eq(actual: &[f64], expected: &[f64], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: length");
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}[{index}]: actual={actual:.17e}, expected={expected:.17e}"
        );
    }
}

fn assert_replay(first: &WfgEvaluation, second: &WfgEvaluation, name: &str) {
    assert_slice_bits_eq(first.transformed(), second.transformed(), name);
    assert_slice_bits_eq(first.reduced(), second.reduced(), name);
    assert_slice_bits_eq(first.positioned(), second.positioned(), name);
    assert_slice_bits_eq(first.shape(), second.shape(), name);
    assert_slice_bits_eq(first.objectives(), second.objectives(), name);
}

fn assert_unit_trace(values: &[f64], name: &str, trace: &str) {
    for (index, &value) in values.iter().enumerate() {
        assert!(
            value.is_finite() && (0.0..=1.0).contains(&value),
            "{name} {trace}[{index}] escaped [0,1]: {value:.17e}"
        );
    }
}

#[test]
fn public_family_matrix_matches_the_frozen_common_probe() {
    let mut observed = Vec::with_capacity(9);

    for (problem, expected) in PublicWfg::family()
        .unwrap()
        .into_iter()
        .zip(EXPECTED_OBJECTIVES)
    {
        let name = problem.name();
        assert_eq!(problem.objectives(), OBJECTIVES, "{name} objectives");
        assert_eq!(
            problem.position_parameters(),
            POSITION_PARAMETERS,
            "{name} position parameters"
        );
        assert_eq!(
            problem.distance_parameters(),
            DISTANCE_PARAMETERS,
            "{name} distance parameters"
        );
        assert_eq!(problem.dimension(), DIMENSION, "{name} dimension");

        let first = problem.evaluate(&MATRIX_INPUT).unwrap();
        let second = problem.evaluate(&MATRIX_INPUT).unwrap();
        assert_eq!(
            first.transformed().len(),
            problem.transformed_dimension(),
            "{name} transformed dimension"
        );
        assert_eq!(first.reduced().len(), OBJECTIVES, "{name} reductions");
        assert_eq!(first.positioned().len(), OBJECTIVES, "{name} positions");
        assert_eq!(first.shape().len(), OBJECTIVES, "{name} shapes");
        assert_eq!(first.objectives().len(), OBJECTIVES, "{name} outputs");

        for (trace, values) in [
            ("transformed", first.transformed()),
            ("reduced", first.reduced()),
            ("positioned", first.positioned()),
            ("shape", first.shape()),
        ] {
            assert_unit_trace(values, name, trace);
        }
        assert!(
            first
                .objectives()
                .iter()
                .all(|value| value.is_finite() && *value >= 0.0),
            "{name} objective range"
        );

        if problem.has_degenerate_positioning() {
            assert!(
                first
                    .reduced()
                    .iter()
                    .zip(first.positioned())
                    .any(|(&reduced, &positioned)| reduced.to_bits() != positioned.to_bits()),
                "WFG3 must expose its degenerate reconstruction"
            );
        } else {
            assert_slice_bits_eq(first.positioned(), first.reduced(), name);
        }

        assert_slice_close(first.objectives(), &expected, name);
        assert_replay(&first, &second, name);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives(), name);
        observed.push((name, first.into_objectives()));
    }

    for (index, (left_name, left)) in observed.iter().enumerate() {
        for (right_name, right) in observed.iter().skip(index + 1) {
            assert!(
                left.iter()
                    .zip(right)
                    .any(|(&left, &right)| left.to_bits() != right.to_bits()),
                "{left_name} and {right_name} unexpectedly aliased"
            );
        }
    }
}

#[test]
fn public_family_admission_is_structured_and_uniform() {
    let family = PublicWfg::family().unwrap();

    for problem in family {
        let name = problem.name();
        assert_eq!(
            problem.evaluate(&[0.0; DIMENSION - 1]).unwrap_err(),
            WfgError::WrongInputLength {
                expected: DIMENSION,
                actual: DIMENSION - 1,
            },
            "{name} short input"
        );
        assert_eq!(
            problem.evaluate(&[0.0; DIMENSION + 1]).unwrap_err(),
            WfgError::WrongInputLength {
                expected: DIMENSION,
                actual: DIMENSION + 1,
            },
            "{name} long input"
        );

        for (value, expected) in [
            (
                f64::NAN,
                WfgError::NonFiniteInput {
                    component: 5,
                    bits: f64::NAN.to_bits(),
                },
            ),
            (
                f64::INFINITY,
                WfgError::NonFiniteInput {
                    component: 5,
                    bits: f64::INFINITY.to_bits(),
                },
            ),
            (
                f64::NEG_INFINITY,
                WfgError::NonFiniteInput {
                    component: 5,
                    bits: f64::NEG_INFINITY.to_bits(),
                },
            ),
            (
                -f64::EPSILON,
                WfgError::InputOutOfRange {
                    component: 5,
                    bits: (-f64::EPSILON).to_bits(),
                },
            ),
            (
                1.0 + f64::EPSILON,
                WfgError::InputOutOfRange {
                    component: 5,
                    bits: (1.0 + f64::EPSILON).to_bits(),
                },
            ),
        ] {
            let mut input = [0.0; DIMENSION];
            input[5] = value;
            assert_eq!(problem.evaluate(&input).unwrap_err(), expected, "{name}");
        }
    }
}

#[test]
fn public_family_constructor_rules_preserve_variant_constraints() {
    assert_eq!(
        Wfg1::new(1, 1, 2).unwrap_err(),
        WfgError::TooFewObjectives { objectives: 1 }
    );
    assert_eq!(
        Wfg8::new(4, 4, 2).unwrap_err(),
        WfgError::PositionParametersNotDivisible {
            position_parameters: 4,
            groups: 3,
        }
    );
    assert_eq!(
        Wfg9::new(3, 4, 0).unwrap_err(),
        WfgError::NoDistanceParameters
    );
    for error in [
        Wfg2::new(4, 6, 3).unwrap_err(),
        Wfg3::new(4, 6, 3).unwrap_err(),
    ] {
        assert_eq!(
            error,
            WfgError::DistanceParametersNotEven {
                distance_parameters: 3,
            }
        );
    }
    assert!(Wfg1::new(4, 6, 3).is_ok());
    assert!(Wfg4::new(4, 6, 3).is_ok());
    assert!(Wfg5::new(4, 6, 3).is_ok());
    assert!(Wfg6::new(4, 6, 3).is_ok());
    assert!(Wfg7::new(4, 6, 3).is_ok());
    assert!(Wfg8::new(4, 6, 3).is_ok());
    assert!(Wfg9::new(4, 6, 3).is_ok());
}

#[test]
fn public_family_boundaries_remain_finite_and_bitwise_replayable() {
    for fill in [0.0, 0.35, 1.0] {
        let input = [fill; DIMENSION];
        for problem in PublicWfg::family().unwrap() {
            let name = problem.name();
            let first = problem.evaluate(&input).unwrap();
            let second = problem.evaluate(&input).unwrap();
            assert_replay(&first, &second, name);
            assert!(
                first.objectives().iter().all(|value| value.is_finite()),
                "{name} boundary {fill}"
            );
        }
    }
}
