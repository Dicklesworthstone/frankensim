//! G0/G1/G3/G5 battery for block composition, solver admission, FGMRES replay, and
//! Newton--Krylov globalization. Each case emits one deterministic JSON-line
//! summary; assertion messages retain the exact failing iteration/decision.

use fs_solver::{
    BlockError, BlockOperator2, BlockOperator3, BlockSchur2, DefinitenessEvidence, FgmresState,
    FlexiblePreconditioner, Globalization, GlobalizationDecision, LineSearchConfig, LinearOp,
    LinearSolverKind, LinearSystemFinding, LinearSystemVerifier, LinearVerificationError,
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonStallDiagnosis, NonlinearProblem,
    NullspaceEvidence, PreconditionerClass, RealEquivalentComplexOp, RectLinearOp, SchurSolveSign,
    SolverAdmissionError, SourceCompatibility, SquareBlock, SymmetryEvidence, TrustRegionConfig,
    ZeroBlock, admit_linear_solver, verify_linear_system,
};
use fs_sparse::precond::Precond;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-solver/nonlinear-block\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

#[derive(Debug, Clone)]
struct DenseRect {
    rows: usize,
    cols: usize,
    values: Vec<f64>,
}

impl DenseRect {
    fn new(rows: usize, cols: usize, values: &[f64]) -> Self {
        assert_eq!(values.len(), rows * cols);
        Self {
            rows,
            cols,
            values: values.to_vec(),
        }
    }
}

impl RectLinearOp for DenseRect {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.cols);
        assert_eq!(y.len(), self.rows);
        y.fill(0.0);
        for (row, output) in y.iter_mut().enumerate() {
            for (column, input) in x.iter().copied().enumerate() {
                *output = self.values[row * self.cols + column].mul_add(input, *output);
            }
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.rows);
        assert_eq!(y.len(), self.cols);
        y.fill(0.0);
        for (row, input) in x.iter().copied().enumerate() {
            for (column, output) in y.iter_mut().enumerate() {
                *output = self.values[row * self.cols + column].mul_add(input, *output);
            }
        }
    }
}

#[derive(Debug, Clone)]
struct DenseSquare(DenseRect);

impl DenseSquare {
    fn new(n: usize, values: &[f64]) -> Self {
        Self(DenseRect::new(n, n, values))
    }
}

impl LinearOp for DenseSquare {
    fn n(&self) -> usize {
        self.0.rows
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply(x, y);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply_transpose(x, y);
    }
}

#[test]
fn block_operator_refuses_zero_aggregate_dimension() {
    let zero = ZeroBlock::new(0, 0);
    let blocks: [[&dyn RectLinearOp; 2]; 2] = [[&zero, &zero], [&zero, &zero]];

    assert!(matches!(
        BlockOperator2::new(blocks),
        Err(BlockError::Empty)
    ));
    verdict(
        "block-zero-dimension-refusal",
        "a populated block table cannot admit a vacuous zero-dimensional solver operator",
    );
}

#[test]
fn block_2x2_3x3_and_real_equivalent_match_dense_actions() {
    let a00 = DenseSquare::new(2, &[2.0, 1.0, -1.0, 3.0]);
    let a01 = DenseRect::new(2, 1, &[4.0, 5.0]);
    let a10 = DenseRect::new(1, 2, &[6.0, -2.0]);
    let a11 = DenseSquare::new(1, &[7.0]);
    let a00_block = SquareBlock::new(&a00);
    let a11_block = SquareBlock::new(&a11);
    let blocks: [[&dyn RectLinearOp; 2]; 2] = [[&a00_block, &a01], [&a10, &a11_block]];
    let operator = BlockOperator2::new(blocks).expect("valid 2x2 partition");
    let x = [1.0, 2.0, -1.0];
    let mut y = [0.0; 3];
    LinearOp::apply(&operator, &x, &mut y);
    assert_eq!(y, [0.0, 0.0, -5.0]);
    let mut yt = [0.0; 3];
    LinearOp::apply_transpose(&operator, &x, &mut yt);
    assert_eq!(yt, [-6.0, 9.0, 7.0]);

    let one = DenseSquare::new(1, &[1.0]);
    let two = DenseSquare::new(1, &[2.0]);
    let three = DenseSquare::new(1, &[3.0]);
    let four = DenseSquare::new(1, &[4.0]);
    let five = DenseSquare::new(1, &[5.0]);
    let six = DenseSquare::new(1, &[6.0]);
    let seven = DenseSquare::new(1, &[7.0]);
    let eight = DenseSquare::new(1, &[8.0]);
    let nine = DenseSquare::new(1, &[9.0]);
    let one_block = SquareBlock::new(&one);
    let two_block = SquareBlock::new(&two);
    let three_block = SquareBlock::new(&three);
    let four_block = SquareBlock::new(&four);
    let five_block = SquareBlock::new(&five);
    let six_block = SquareBlock::new(&six);
    let seven_block = SquareBlock::new(&seven);
    let eight_block = SquareBlock::new(&eight);
    let nine_block = SquareBlock::new(&nine);
    let blocks3: [[&dyn RectLinearOp; 3]; 3] = [
        [&one_block, &two_block, &three_block],
        [&four_block, &five_block, &six_block],
        [&seven_block, &eight_block, &nine_block],
    ];
    let operator3 = BlockOperator3::new(blocks3).expect("valid 3x3 partition");
    let mut y3 = [0.0; 3];
    LinearOp::apply(&operator3, &[1.0, -1.0, 2.0], &mut y3);
    assert_eq!(y3, [5.0, 11.0, 17.0]);
    let mut y3_transpose = [0.0; 3];
    LinearOp::apply_transpose(&operator3, &[1.0, -1.0, 2.0], &mut y3_transpose);
    assert_eq!(y3_transpose, [11.0, 13.0, 15.0]);

    let real = DenseSquare::new(2, &[2.0, 0.0, 0.0, 3.0]);
    let imaginary = DenseSquare::new(2, &[0.0, 1.0, 2.0, 0.0]);
    let complex = RealEquivalentComplexOp::new(&real, &imaginary).expect("same dimension");
    let mut yc = [0.0; 4];
    LinearOp::apply(&complex, &[1.0, 2.0, -1.0, 4.0], &mut yc);
    assert_eq!(yc, [-2.0, 8.0, 0.0, 14.0]);
    let mut yc_transpose = [0.0; 4];
    LinearOp::apply_transpose(&complex, &[1.0, 2.0, -1.0, 4.0], &mut yc_transpose);
    assert_eq!(yc_transpose, [10.0, 5.0, -6.0, 11.0]);
    verdict(
        "block-actions",
        "2x2/3x3 rectangular composition, transpose, and real-equivalent complex action match dense references",
    );
}

#[derive(Debug, Clone)]
struct DiagonalInverse(Vec<f64>);

impl Precond for DiagonalInverse {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        assert_eq!(residual.len(), self.0.len());
        assert_eq!(output.len(), self.0.len());
        for ((out, right), inverse) in output.iter_mut().zip(residual).zip(&self.0) {
            *out = right * inverse;
        }
    }
}

#[test]
fn exact_block_schur_matches_monolithic_saddle_reference() {
    // K = [diag(2,3)  [1;1]; [1,1]  0]. Its positive complement is
    // S = B A^-1 B^T = 1/2 + 1/3 = 5/6.
    let a = DenseSquare::new(2, &[2.0, 0.0, 0.0, 3.0]);
    let bt = DenseRect::new(2, 1, &[1.0, 1.0]);
    let b = DenseRect::new(1, 2, &[1.0, 1.0]);
    let zero = DenseSquare::new(1, &[0.0]);
    let a_block = SquareBlock::new(&a);
    let zero_block = SquareBlock::new(&zero);
    let blocks: [[&dyn RectLinearOp; 2]; 2] = [[&a_block, &bt], [&b, &zero_block]];
    let saddle = BlockOperator2::new(blocks).expect("saddle partition");
    let a_inverse = DiagonalInverse(vec![0.5, 1.0 / 3.0]);
    let schur_inverse = DiagonalInverse(vec![6.0 / 5.0]);
    let inverse = BlockSchur2::new(
        2,
        1,
        &a_inverse,
        &schur_inverse,
        &b,
        &bt,
        SchurSolveSign::Negative,
    )
    .expect("exact block LDU");
    for rhs in [[1.0, 2.0, 3.0], [-4.0, 0.5, 2.0], [0.0, 0.0, 1.0]] {
        let mut solution = [0.0; 3];
        Precond::apply(&inverse, &rhs, &mut solution);
        let mut reconstructed = [0.0; 3];
        LinearOp::apply(&saddle, &solution, &mut reconstructed);
        for (index, (actual, expected)) in reconstructed.iter().zip(rhs).enumerate() {
            assert!(
                (actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs().max(1.0),
                "Schur reconstruction entry {index}: actual={actual:.17e}, expected={expected:.17e}, solution={solution:?}"
            );
        }
    }
    verdict(
        "schur-reference",
        "exact injected A and positive-complement inverses reconstruct three manufactured saddle right-hand sides",
    );
}

#[derive(Debug, Clone, Copy)]
struct FindingVerifier {
    finding: LinearSystemFinding,
}

impl LinearSystemVerifier for FindingVerifier {
    fn verifier_id(&self) -> &str {
        "test/exact-dense-structure-verifier/v1"
    }

    fn verify(
        &self,
        _operator: &dyn LinearOp,
        _rhs: &[f64],
    ) -> Result<LinearSystemFinding, LinearVerificationError> {
        Ok(self.finding)
    }
}

#[test]
fn admission_refuses_physics_name_shortcuts_and_variable_plain_gmres() {
    let nonsymmetric = DenseSquare::new(2, &[2.0, 1.0, 0.0, 1.0]);
    let general = verify_linear_system(
        &nonsymmetric,
        &[1.0, 2.0],
        &FindingVerifier {
            finding: LinearSystemFinding {
                symmetry: SymmetryEvidence::Nonsymmetric,
                definiteness: DefinitenessEvidence::Unknown,
                nullspace: NullspaceEvidence::Trivial,
                source: SourceCompatibility::Compatible,
                preconditioner: PreconditionerClass::Variable,
            },
        },
    )
    .expect("coherent nonsymmetric finding");
    assert_eq!(
        admit_linear_solver(LinearSolverKind::Cg, general.clone()),
        Err(SolverAdmissionError::SymmetryRequired)
    );
    assert_eq!(
        admit_linear_solver(LinearSolverKind::Gmres, general.clone()),
        Err(SolverAdmissionError::PreconditionerIncompatible {
            solver: LinearSolverKind::Gmres,
            preconditioner: PreconditionerClass::Variable,
        })
    );
    assert_eq!(
        admit_linear_solver(LinearSolverKind::Fgmres, general)
            .expect("FGMRES admits general variable-preconditioned system")
            .kind(),
        LinearSolverKind::Fgmres
    );
    let unresolved = verify_linear_system(
        &nonsymmetric,
        &[1.0, 2.0],
        &FindingVerifier {
            finding: LinearSystemFinding {
                symmetry: SymmetryEvidence::Nonsymmetric,
                definiteness: DefinitenessEvidence::Unknown,
                nullspace: NullspaceEvidence::Unresolved,
                source: SourceCompatibility::Compatible,
                preconditioner: PreconditionerClass::Variable,
            },
        },
    )
    .expect("unresolved nullspace is coherent evidence, but not admissible");
    assert_eq!(
        admit_linear_solver(LinearSolverKind::Fgmres, unresolved),
        Err(SolverAdmissionError::NullspaceUnresolved)
    );
    assert_eq!(
        verify_linear_system(
            &nonsymmetric,
            &[1.0, 2.0],
            &FindingVerifier {
                finding: LinearSystemFinding {
                    symmetry: SymmetryEvidence::Nonsymmetric,
                    definiteness: DefinitenessEvidence::Unknown,
                    nullspace: NullspaceEvidence::Projected { dimension: 3 },
                    source: SourceCompatibility::Compatible,
                    preconditioner: PreconditionerClass::Variable,
                },
            },
        ),
        Err(LinearVerificationError::ProjectedNullspaceTooLarge {
            dimension: 3,
            operator: 2,
        })
    );
    assert_eq!(
        verify_linear_system(
            &nonsymmetric,
            &[1.0, 2.0],
            &FindingVerifier {
                finding: LinearSystemFinding {
                    symmetry: SymmetryEvidence::Nonsymmetric,
                    definiteness: DefinitenessEvidence::Indefinite,
                    nullspace: NullspaceEvidence::Trivial,
                    source: SourceCompatibility::Compatible,
                    preconditioner: PreconditionerClass::FixedSpd,
                },
            },
        ),
        Err(LinearVerificationError::ContradictoryIndefiniteFinding)
    );

    let spd = DenseSquare::new(2, &[2.0, -1.0, -1.0, 2.0]);
    let verified_spd = verify_linear_system(
        &spd,
        &[1.0, 0.0],
        &FindingVerifier {
            finding: LinearSystemFinding {
                symmetry: SymmetryEvidence::Symmetric,
                definiteness: DefinitenessEvidence::PositiveDefinite,
                nullspace: NullspaceEvidence::Trivial,
                source: SourceCompatibility::Compatible,
                preconditioner: PreconditionerClass::FixedSpd,
            },
        },
    )
    .expect("coherent SPD finding");
    assert_eq!(
        admit_linear_solver(LinearSolverKind::Minres, verified_spd.clone()),
        Err(SolverAdmissionError::IndefiniteRequired)
    );
    assert!(admit_linear_solver(LinearSolverKind::Cg, verified_spd).is_ok());
    let indefinite = DenseSquare::new(2, &[1.0, 0.0, 0.0, -1.0]);
    let verified_indefinite = verify_linear_system(
        &indefinite,
        &[1.0, 1.0],
        &FindingVerifier {
            finding: LinearSystemFinding {
                symmetry: SymmetryEvidence::Symmetric,
                definiteness: DefinitenessEvidence::Indefinite,
                nullspace: NullspaceEvidence::Trivial,
                source: SourceCompatibility::Compatible,
                preconditioner: PreconditionerClass::FixedSpd,
            },
        },
    )
    .expect("coherent symmetric-indefinite finding");
    assert!(admit_linear_solver(LinearSolverKind::Minres, verified_indefinite).is_ok());
    verdict(
        "solver-admission",
        "nonsymmetric-to-CG, SPD-to-MINRES, and variable-to-GMRES refuse; verified SPD/indefinite/general findings admit CG/MINRES/FGMRES",
    );
}

#[derive(Debug, Clone, Copy)]
struct CyclingDiagonal;

impl FlexiblePreconditioner for CyclingDiagonal {
    fn apply(&self, logical_iteration: usize, residual: &[f64], output: &mut [f64]) {
        let diagonal = match logical_iteration % 3 {
            0 => [0.5, 1.0, 1.5],
            1 => [1.5, 0.75, 1.0],
            _ => [1.0, 1.25, 0.5],
        };
        assert_eq!(residual.len(), diagonal.len());
        assert_eq!(output.len(), diagonal.len());
        for ((out, value), scale) in output.iter_mut().zip(residual).zip(diagonal) {
            *out = scale * value;
        }
    }
}

fn bit_pattern(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn fgmres_variable_preconditioner_resume_is_bitwise() {
    let operator = DenseSquare::new(3, &[4.0, 1.0, 0.0, -1.0, 3.0, 1.0, 0.0, 2.0, 5.0]);
    let rhs = [1.0, -2.0, 3.0];
    let mut straight = FgmresState::new(&rhs, 2);
    let straight_report = straight.run(&operator, &CyclingDiagonal, &rhs, 1.0e-12, 20);
    assert!(
        straight_report.converged,
        "straight FGMRES failed: {straight_report:?}"
    );

    let mut split = FgmresState::new(&rhs, 2);
    let prefix = split.run(&operator, &CyclingDiagonal, &rhs, 1.0e-12, 2);
    assert!(
        !prefix.converged,
        "fixture must exercise a real resume boundary"
    );
    let split_report = split.run(&operator, &CyclingDiagonal, &rhs, 1.0e-12, 18);
    assert!(
        split_report.converged,
        "resumed FGMRES failed: {split_report:?}"
    );
    assert_eq!(bit_pattern(&straight.x), bit_pattern(&split.x));
    assert_eq!(
        bit_pattern(&straight.history),
        bit_pattern(&split.history),
        "cycle-end residual history moved across split run"
    );
    verdict(
        "fgmres-resume",
        "logical-iteration-keyed non-collinear diagonal preconditioners exercise stored z_j directions and split cycles remain bitwise equal",
    );
}

#[derive(Debug, Clone, Copy)]
struct SquareRootTwo;

impl NonlinearProblem for SquareRootTwo {
    fn dimension(&self) -> usize {
        1
    }

    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        residual[0] = x[0].mul_add(x[0], -2.0);
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output[0] = 2.0 * x[0] * direction[0];
    }
}

#[derive(Debug, Clone, Copy)]
struct NoRealRoot;

impl NonlinearProblem for NoRealRoot {
    fn dimension(&self) -> usize {
        1
    }

    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        residual[0] = x[0].mul_add(x[0], 1.0);
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output[0] = 2.0 * x[0] * direction[0];
    }
}

#[derive(Debug, Clone, Copy)]
struct HugeResidual;

impl NonlinearProblem for HugeResidual {
    fn dimension(&self) -> usize {
        2
    }

    fn residual(&self, _x: &[f64], residual: &mut [f64]) {
        residual.fill(f64::MAX);
    }

    fn jacobian_apply(&self, _x: &[f64], _direction: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
}

#[derive(Debug, Clone, Copy)]
struct BratuTwo {
    lambda: f64,
}

impl NonlinearProblem for BratuTwo {
    fn dimension(&self) -> usize {
        2
    }

    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        residual[0] = 2.0 * x[0] - x[1] - self.lambda * fs_math::det::exp(x[0]);
        residual[1] = 2.0 * x[1] - x[0] - self.lambda * fs_math::det::exp(x[1]);
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output[0] = (2.0 - self.lambda * fs_math::det::exp(x[0])) * direction[0] - direction[1];
        output[1] = (2.0 - self.lambda * fs_math::det::exp(x[1])) * direction[1] - direction[0];
    }
}

#[derive(Debug, Clone, Copy)]
struct SaturatingDiffusion;

impl NonlinearProblem for SaturatingDiffusion {
    fn dimension(&self) -> usize {
        2
    }

    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        residual[0] = 2.0 * x[0] - x[1] + fs_math::det::tanh(x[0]) - 1.0;
        residual[1] = 2.0 * x[1] - x[0] + fs_math::det::tanh(x[1]) + 0.5;
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        let tanh0 = fs_math::det::tanh(x[0]);
        let tanh1 = fs_math::det::tanh(x[1]);
        output[0] = (3.0 - tanh0 * tanh0) * direction[0] - direction[1];
        output[1] = (3.0 - tanh1 * tanh1) * direction[1] - direction[0];
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One G1 globalization/resume/convergence receipt.
fn newton_krylov_globalizes_logs_and_resumes_bitwise() {
    let line_config = NewtonKrylovConfig {
        absolute_tolerance: 1.0e-14,
        relative_tolerance: 1.0e-13,
        linear_restart: 4,
        max_linear_cycles: 4,
        forcing_minimum: 1.0e-14,
        forcing_maximum: 0.1,
        forcing_gamma: 0.5,
        forcing_exponent: 1.5,
        globalization: Globalization::LineSearch(LineSearchConfig::default()),
    };
    let initial = NewtonKrylovState::new(&SquareRootTwo, vec![1.5], line_config)
        .expect("sqrt-two checkpoint");
    let mut straight = initial.clone();
    let straight_report = straight.run(&SquareRootTwo, 12);
    assert!(
        straight_report.converged,
        "sqrt-two Newton: {straight_report:?}"
    );
    assert!((straight.x[0] - fs_math::det::sqrt(2.0)).abs() <= 4.0 * f64::EPSILON);
    for entry in straight_report
        .history
        .iter()
        .filter(|entry| entry.residual_before > 1.0e-8)
        .skip(1)
    {
        assert!(
            entry.residual_after <= 2.0 * entry.residual_before * entry.residual_before,
            "quadratic-rate miss at iteration {}: before={:.17e}, after={:.17e}",
            entry.iteration,
            entry.residual_before,
            entry.residual_after
        );
    }
    let mut split = initial;
    let prefix = split.run(&SquareRootTwo, 2);
    assert!(!prefix.converged, "split fixture needs a resumable prefix");
    let resumed = split.run(&SquareRootTwo, 10);
    assert!(resumed.converged);
    assert_eq!(bit_pattern(&straight.x), bit_pattern(&split.x));
    assert_eq!(straight_report.history, resumed.history);

    let rejection_config = NewtonKrylovConfig {
        globalization: Globalization::LineSearch(LineSearchConfig {
            minimum_step: 0.25,
            ..LineSearchConfig::default()
        }),
        ..NewtonKrylovConfig::default()
    };
    let mut rejected = NewtonKrylovState::new(&NoRealRoot, vec![0.1], rejection_config)
        .expect("line-search rejection checkpoint");
    let rejection_report = rejected.run(&NoRealRoot, 1);
    assert_eq!(
        rejection_report.diagnosis,
        Some(NewtonStallDiagnosis::LineSearchRejected)
    );
    assert_eq!(rejection_report.history.len(), 1);
    assert_eq!(
        rejection_report.history[0].decision,
        GlobalizationDecision::LineSearchRejected
    );
    assert_eq!(rejection_report.history[0].step_length, 0.25);
    assert!(matches!(
        NewtonKrylovState::new(&HugeResidual, vec![0.0, 0.0], NewtonKrylovConfig::default()),
        Err(NewtonError::NonFiniteResidualNorm { iteration: 0 })
    ));

    let bratu_config = NewtonKrylovConfig {
        globalization: Globalization::LineSearch(LineSearchConfig::default()),
        ..NewtonKrylovConfig::default()
    };
    let mut bratu =
        NewtonKrylovState::new(&BratuTwo { lambda: 0.25 }, vec![0.0, 0.0], bratu_config)
            .expect("Bratu checkpoint");
    let bratu_report = bratu.run(&BratuTwo { lambda: 0.25 }, 24);
    assert!(
        bratu_report.converged,
        "Bratu-type Newton: {bratu_report:?}"
    );

    let trust_config = NewtonKrylovConfig {
        absolute_tolerance: 1.0e-12,
        relative_tolerance: 1.0e-10,
        globalization: Globalization::TrustRegion(TrustRegionConfig {
            initial_radius: 0.05,
            minimum_radius: 1.0e-12,
            maximum_radius: 10.0,
            acceptance_ratio: 0.05,
            expansion_ratio: 0.7,
            shrink: 0.25,
            expansion: 2.0,
        }),
        ..NewtonKrylovConfig::default()
    };
    let mut saturating =
        NewtonKrylovState::new(&SaturatingDiffusion, vec![4.0, -3.0], trust_config)
            .expect("saturating-diffusion checkpoint");
    let saturation_report = saturating.run(&SaturatingDiffusion, 64);
    assert!(
        saturation_report.converged,
        "trust-globalized saturation solve: {saturation_report:?}"
    );
    assert!(
        saturation_report.history.iter().any(|entry| {
            entry.decision == GlobalizationDecision::TrustRegionAccepted && entry.step_length < 1.0
        }),
        "small initial trust radius must be visible in telemetry: {:?}",
        saturation_report.history
    );
    assert!(saturation_report.history.iter().all(|entry| {
        entry.forcing.is_finite()
            && entry.linear_relative_residual.is_finite()
            && entry.newton_step_norm.is_finite()
            && entry.step_length.is_finite()
    }));
    verdict(
        "newton-globalization",
        &format!(
            "sqrt2={} iterations with quadratic tail and bitwise split replay; Bratu={} iterations; \
             saturating diffusion={} trust attempts with radius-limited decisions",
            straight_report.iterations, bratu_report.iterations, saturation_report.iterations
        ),
    );
}
