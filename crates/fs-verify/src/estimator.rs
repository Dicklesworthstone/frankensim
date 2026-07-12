//! The VERIFIER: equilibrated-flux a-posteriori bounds (Prager–Synge,
//! 1D elliptic class), interval-evaluated to VERIFIED color.
//!
//! The rigor structure: ANY σ with `σ′ = −f` yields the guaranteed
//! bound `‖(u − u_h)′‖ ≤ ‖σ − u_h′‖` — so the free constant in
//! `σ = c − F` is optimized in plain f64 for TIGHTNESS while the bound
//! itself is evaluated with outward-rounded intervals over exact Gauss
//! quadrature (polynomial data ⇒ the quadrature identity is exact;
//! only rounding needs enclosing). Malformed inputs and unusable
//! enclosures FAIL CLOSED as structured refusals: no color, ever.

use crate::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS, MmsProblem, gauss5,
    require_converged, true_energy_error, try_zeroed, validate_candidate, validate_problem,
};
use crate::interval::Iv;
use fs_evidence::Color;
use std::fmt::Write as _;

/// Largest mesh admitted by the synchronous v0 verifier.
pub const MAX_VERIFIER_MESH_NODES: usize = MAX_FEM1D_MESH_NODES;
/// Exactness envelope for the manufactured solution: degree at most five.
pub const MAX_VERIFIER_POLY_COEFFICIENTS: usize = MAX_FEM1D_POLY_COEFFICIENTS;

/// Estimator families (Proposal D's independence escalation needs at
/// least two registered per class).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimatorFamily {
    /// Equilibrated flux (guaranteed, constant-free — the verifier).
    EquilibratedFlux,
    /// Hierarchical (refined-mesh comparison — independent, NOT
    /// guaranteed; the falsifier's cross-check).
    Hierarchical,
}

/// Polynomial role carried by a structured verifier refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierPolynomial {
    /// Manufactured exact-solution metadata (`u`).
    ExactSolution,
    /// Canonical forcing (`f = -u''`).
    Forcing,
    /// Canonical zero-constant antiderivative of the forcing (`big_f`).
    ForcingAntiderivative,
}

/// Stable, structured reason why no verifier authority was issued.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierRefusal {
    /// Mesh length is outside the bounded `2..=MAX_VERIFIER_MESH_NODES` class.
    MeshNodeCount,
    /// One polynomial is empty or exceeds the degree-five exactness class.
    PolynomialCoefficientCount {
        /// Polynomial whose resource envelope was violated.
        polynomial: VerifierPolynomial,
    },
    /// Candidate and mesh lengths differ.
    CandidateLength,
    /// Tolerance is non-finite or non-positive.
    InvalidTolerance,
    /// Mesh endpoints are not canonical `+0.0` and `1.0`.
    MeshDomain,
    /// A mesh coordinate is non-finite or the mesh is not strictly increasing.
    MeshCoordinates,
    /// A candidate value is non-finite.
    CandidateNonFinite,
    /// Candidate endpoints are not canonical homogeneous `+0.0` values.
    CandidateBoundary,
    /// One polynomial contains a non-finite coefficient.
    PolynomialNonFinite {
        /// Polynomial containing the non-finite coefficient.
        polynomial: VerifierPolynomial,
    },
    /// The exact-solution polynomial does not vanish canonically at both ends.
    ExactSolutionBoundary,
    /// A public derived polynomial differs from the canonical value recomputed from `u`.
    DerivedPolynomialMismatch {
        /// Public derived polynomial that did not match its canonical value.
        polynomial: VerifierPolynomial,
    },
    /// The optional tightness constant could not be computed finitely.
    NonFiniteTightness,
    /// Interval construction produced a non-finite, reversed, or unusable enclosure.
    InvalidEnclosure,
}

impl VerifierRefusal {
    /// Stable identifier for diagnostics and ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            Self::MeshNodeCount => "mesh-node-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-coefficient-count",
            Self::CandidateLength => "candidate-length",
            Self::InvalidTolerance => "invalid-tolerance",
            Self::MeshDomain => "mesh-domain",
            Self::MeshCoordinates => "mesh-coordinates",
            Self::CandidateNonFinite => "candidate-non-finite",
            Self::CandidateBoundary => "candidate-boundary",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-non-finite",
            Self::ExactSolutionBoundary => "exact-solution-boundary",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::Forcing,
            } => "derived-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "derived-big-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "derived-u-mismatch",
            Self::NonFiniteTightness => "non-finite-tightness",
            Self::InvalidEnclosure => "invalid-enclosure",
        }
    }
}

impl EstimatorFamily {
    /// Stable id for ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            EstimatorFamily::EquilibratedFlux => "equilibrated-flux-1d",
            EstimatorFamily::Hierarchical => "hierarchical-h2",
        }
    }
}

/// The verifier's verdict on one candidate.
#[derive(Debug, Clone)]
pub struct VerifierReport {
    /// The certified error-bound enclosure (energy norm).
    pub bound: Iv,
    /// Accept ⟺ `bound.hi ≤ tolerance` for an admitted finite report.
    pub accept: bool,
    /// The verified color carried by an ACCEPT (`None` on reject or refusal).
    pub color: Option<Color>,
    /// The tolerance tested against (feeds the planner).
    pub tolerance: f64,
    /// Estimator family id.
    pub family: &'static str,
    /// FNV hash of the reconstructed flux (ledger identity).
    pub flux_hash: u64,
    /// Structured refusal (`None` only when a finite bound was produced).
    pub refusal: Option<VerifierRefusal>,
}

impl VerifierReport {
    /// The review-round-3 ledger row (structured, never stdout).
    #[must_use]
    pub fn to_row(&self, problem: &str, oracle_error: f64) -> String {
        let problem = json_escape(problem);
        let family = json_escape(self.family);
        let bound_lo = finite_scientific(self.bound.lo);
        let bound_hi = finite_scientific(self.bound.hi);
        let oracle = finite_scientific(oracle_error);
        let tolerance = finite_scientific(self.tolerance);
        let effectivity = if self.refusal.is_none()
            && oracle_error.is_finite()
            && oracle_error > 0.0
            && self.bound.hi.is_finite()
        {
            finite_fixed(self.bound.hi / oracle_error)
        } else if self.refusal.is_none() && oracle_error == 0.0 {
            "1.0000".to_string()
        } else {
            "null".to_string()
        };
        let refusal = self.refusal.map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", reason.id()),
        );
        let verdict = if self.refusal.is_some() {
            "refused"
        } else if self.accept {
            "accept"
        } else {
            "reject"
        };
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"problem\":\"{problem}\",\"estimator_family_id\":\"{}\",\
             \"flux_hash\":\"{:016X}\",\"bound_lo\":{bound_lo},\"bound_hi\":{bound_hi},\
             \"oracle_true_error\":{oracle},\"effectivity\":{effectivity},\
             \"verdict\":\"{verdict}\",\"tolerance\":{tolerance},\"refusal\":{refusal}}}",
            family, self.flux_hash,
        );
        s
    }
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

fn finite_scientific(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

fn finite_fixed(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.4}")
    } else {
        "null".to_string()
    }
}

fn fnv_extend(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn hash_polynomial(mut hash: u64, polynomial: &crate::fem1d::Poly) -> u64 {
    hash = fnv_extend(
        hash,
        &(polynomial.coefficients().len() as u64).to_le_bytes(),
    );
    for coefficient in polynomial.coefficients() {
        hash = fnv_extend(hash, &coefficient.to_bits().to_le_bytes());
    }
    hash
}

fn flux_hash(c_star: f64, f: &crate::fem1d::Poly, big_f: &crate::fem1d::Poly) -> u64 {
    let mut hash = fnv_extend(0xcbf2_9ce4_8422_2325, &c_star.to_bits().to_le_bytes());
    hash = hash_polynomial(hash, f);
    hash_polynomial(hash, big_f)
}

fn poly_bits_equal(left: &crate::fem1d::Poly, right: &crate::fem1d::Poly) -> bool {
    left.coefficients().len() == right.coefficients().len()
        && left
            .coefficients()
            .iter()
            .zip(right.coefficients())
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

fn validate_inputs(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
) -> Result<(crate::fem1d::Poly, crate::fem1d::Poly), VerifierRefusal> {
    if !(2..=MAX_VERIFIER_MESH_NODES).contains(&problem.mesh().len()) {
        return Err(VerifierRefusal::MeshNodeCount);
    }
    if candidate.len() != problem.mesh().len() {
        return Err(VerifierRefusal::CandidateLength);
    }
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        if !(1..=MAX_VERIFIER_POLY_COEFFICIENTS).contains(&polynomial.coefficients().len()) {
            return Err(VerifierRefusal::PolynomialCoefficientCount { polynomial: role });
        }
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(VerifierRefusal::InvalidTolerance);
    }
    if problem.mesh().first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || problem.mesh().last().map(|value| value.to_bits()) != Some(1.0_f64.to_bits())
    {
        return Err(VerifierRefusal::MeshDomain);
    }
    if !problem.mesh().iter().all(|value| value.is_finite())
        || !problem.mesh().windows(2).all(|pair| pair[0] < pair[1])
    {
        return Err(VerifierRefusal::MeshCoordinates);
    }
    if !candidate.iter().all(|value| value.is_finite()) {
        return Err(VerifierRefusal::CandidateNonFinite);
    }
    if candidate.first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || candidate.last().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
    {
        return Err(VerifierRefusal::CandidateBoundary);
    }
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        if !polynomial
            .coefficients()
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(VerifierRefusal::PolynomialNonFinite { polynomial: role });
        }
    }
    if problem
        .exact_solution()
        .coefficients()
        .first()
        .map(|value| value.to_bits())
        != Some(0.0_f64.to_bits())
        || !problem.exact_solution().is_exactly_zero_at_one()
    {
        return Err(VerifierRefusal::ExactSolutionBoundary);
    }

    let expected_f = problem
        .exact_solution()
        .derive()
        .and_then(|derivative| derivative.derive())
        .map(|second_derivative| second_derivative.neg())
        .map_err(|_| VerifierRefusal::PolynomialNonFinite {
            polynomial: VerifierPolynomial::Forcing,
        })?;
    if !poly_bits_equal(problem.forcing(), &expected_f) {
        return Err(VerifierRefusal::DerivedPolynomialMismatch {
            polynomial: VerifierPolynomial::Forcing,
        });
    }
    let expected_big_f =
        expected_f
            .antiderive()
            .map_err(|_| VerifierRefusal::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            })?;
    if !poly_bits_equal(problem.rounded_forcing_antiderivative(), &expected_big_f) {
        return Err(VerifierRefusal::DerivedPolynomialMismatch {
            polynomial: VerifierPolynomial::ForcingAntiderivative,
        });
    }
    Ok((expected_f, expected_big_f))
}

fn tightness_constant(
    problem: &MmsProblem,
    candidate: &[f64],
    big_f: &crate::fem1d::Poly,
) -> Result<f64, VerifierRefusal> {
    let mut mean = 0.0;
    for element in 0..problem.mesh().len() - 1 {
        let (x0, x1) = (problem.mesh()[element], problem.mesh()[element + 1]);
        let h = x1 - x0;
        let slope = (candidate[element + 1] - candidate[element]) / h;
        if !h.is_finite() || h <= 0.0 || !slope.is_finite() {
            return Err(VerifierRefusal::NonFiniteTightness);
        }
        for (point, weight) in gauss5(x0, x1) {
            let value = big_f.eval(point) + slope;
            let contribution = weight * value;
            if !point.is_finite()
                || !weight.is_finite()
                || !value.is_finite()
                || !contribution.is_finite()
            {
                return Err(VerifierRefusal::NonFiniteTightness);
            }
            mean += contribution;
            if !mean.is_finite() {
                return Err(VerifierRefusal::NonFiniteTightness);
            }
        }
    }
    Ok(mean)
}

fn finite_interval(interval: Iv) -> Result<Iv, VerifierRefusal> {
    if interval.lo.is_finite() && interval.hi.is_finite() && interval.lo <= interval.hi {
        Ok(interval)
    } else {
        Err(VerifierRefusal::InvalidEnclosure)
    }
}

fn interval_element_geometry(x0: f64, x1: f64) -> Result<(Iv, Iv, Iv), VerifierRefusal> {
    let x0 = Iv::point(x0);
    let x1 = Iv::point(x1);
    let h = finite_interval(x1.sub(x0))?;
    if h.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    let midpoint = finite_interval(x0.add(x1).scale_pos(0.5))?;
    let half = finite_interval(h.scale_pos(0.5))?;
    if half.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((h, midpoint, half))
}

fn interval_candidate_slope(first: f64, second: f64, h: Iv) -> Result<Iv, VerifierRefusal> {
    let difference = finite_interval(Iv::point(second).sub(Iv::point(first)))?;
    finite_interval(difference.div_pos(h))
}

fn interval_quadrature_geometry(
    midpoint: Iv,
    half: Iv,
    node_constant: f64,
    weight_constant: f64,
) -> Result<(Iv, Iv), VerifierRefusal> {
    let node = finite_interval(midpoint.add(half.mul(iv_c(node_constant))))?;
    let weight = finite_interval(half.mul(iv_c(weight_constant)))?;
    if weight.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((node, weight))
}

fn interval_antiderivative_coefficient(
    coefficient: f64,
    exponent: usize,
) -> Result<Iv, VerifierRefusal> {
    finite_interval(Iv::point(coefficient).div_pos(Iv::point(exponent as f64)))
}

fn interval_forcing_antiderivative(
    forcing: &crate::fem1d::Poly,
    x: Iv,
) -> Result<Iv, VerifierRefusal> {
    // F(x) = x * Horner(f_k / (k + 1)). Coefficient division is itself
    // intervalized: the rounded coefficients in `big_f` are replay metadata,
    // not point enclosures of the exact antiderivative of the authoritative f.
    let mut accumulated = Iv::zero();
    for (degree, coefficient) in forcing.coefficients().iter().copied().enumerate().rev() {
        let antiderivative_coefficient =
            interval_antiderivative_coefficient(coefficient, degree + 1)?;
        accumulated = finite_interval(accumulated.mul(x).add(antiderivative_coefficient))?;
    }
    finite_interval(x.mul(accumulated))
}

fn equilibrated_bound(
    problem: &MmsProblem,
    candidate: &[f64],
    forcing: &crate::fem1d::Poly,
    c_star: f64,
) -> Result<Iv, VerifierRefusal> {
    let mut eta_sq = Iv::zero();
    for element in 0..problem.mesh().len() - 1 {
        let (h, midpoint, half) =
            interval_element_geometry(problem.mesh()[element], problem.mesh()[element + 1])?;
        let slope = interval_candidate_slope(candidate[element], candidate[element + 1], h)?;
        for (node_constant, weight_constant) in GAUSS5_REF {
            let (node, weight) =
                interval_quadrature_geometry(midpoint, half, node_constant, weight_constant)?;
            let antiderivative = interval_forcing_antiderivative(forcing, node)?;
            let residual = finite_interval(Iv::point(c_star).sub(antiderivative).sub(slope))?;
            let contribution = finite_interval(weight.mul(residual.sq()))?;
            eta_sq = finite_interval(eta_sq.add(contribution))?;
        }
    }
    let bound = finite_interval(eta_sq.sqrt())?;
    if bound.lo < 0.0 {
        Err(VerifierRefusal::InvalidEnclosure)
    } else {
        Ok(bound)
    }
}

fn refused(tolerance: f64, reason: VerifierRefusal) -> VerifierReport {
    VerifierReport {
        bound: Iv {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        },
        accept: false,
        color: None,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: 0,
        refusal: Some(reason),
    }
}

/// The equilibrated-flux VERIFIER: certify (or reject) a candidate's
/// nodal values against `tolerance`. The returned bound is a TRUE
/// upper bound on `‖(u − u_h)′‖` whenever the candidate satisfies the
/// boundary conditions; the enclosure is rigorous by outward rounding.
#[must_use]
pub fn verify(problem: &MmsProblem, candidate: &[f64], tolerance: f64) -> VerifierReport {
    let (canonical_f, canonical_big_f) = match validate_inputs(problem, candidate, tolerance) {
        Ok(polynomials) => polynomials,
        Err(reason) => return refused(tolerance, reason),
    };
    // Any finite c is sound. This rounded optimizer affects tightness only.
    let c_star = match tightness_constant(problem, candidate, &canonical_big_f) {
        Ok(value) => value,
        Err(reason) => return refused(tolerance, reason),
    };
    let bound = match equilibrated_bound(problem, candidate, &canonical_f, c_star) {
        Ok(bound) => bound,
        Err(reason) => return refused(tolerance, reason),
    };
    let accept = bound.hi <= tolerance;
    let color = if accept {
        Some(Color::Verified {
            lo: 0.0,
            hi: bound.hi,
        })
    } else {
        None
    };
    VerifierReport {
        bound,
        accept,
        color,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: flux_hash(c_star, &canonical_f, &canonical_big_f),
        refusal: None,
    }
}

const GAUSS5_REF: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_08),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.0, 0.568_888_888_888_888_9),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.906_179_845_938_664, 0.236_926_885_056_189_08),
];

/// One-ulp-widened constant (the tabulated Gauss data carries ~1 ulp
/// of transcription error; widening keeps enclosures honest).
fn iv_c(v: f64) -> Iv {
    Iv {
        lo: crate::interval::down(v),
        hi: crate::interval::up(v),
    }
}

/// The INDEPENDENT second family: hierarchical estimate from a
/// uniformly refined solve (`h/2`). Not guaranteed — the falsifier's
/// cross-check, never a color source.
///
/// # Errors
/// Returns [`Fem1dError`] for malformed inputs, refinement overflow/resource
/// excess, allocation failure, or a non-finite estimate.
pub fn hierarchical_estimate(problem: &MmsProblem, candidate: &[f64]) -> Result<f64, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let fine_nodes = problem
        .mesh()
        .len()
        .checked_mul(2)
        .and_then(|nodes| nodes.checked_sub(1))
        .ok_or(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: usize::MAX,
            limit: MAX_FEM1D_MESH_NODES,
        })?;
    if fine_nodes > MAX_FEM1D_MESH_NODES {
        return Err(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: fine_nodes,
            limit: MAX_FEM1D_MESH_NODES,
        });
    }
    let mut fine_mesh = Vec::new();
    fine_mesh
        .try_reserve_exact(fine_nodes)
        .map_err(|_| Fem1dError::AllocationFailed {
            stage: "hierarchical mesh",
            requested: fine_nodes,
        })?;
    for w in problem.mesh().windows(2) {
        fine_mesh.push(w[0]);
        fine_mesh.push(f64::midpoint(w[0], w[1]));
    }
    fine_mesh.push(problem.mesh()[problem.mesh().len() - 1]);
    let fine = problem.with_mesh(fine_mesh)?;
    let fine_u = crate::fem1d::solve_p1(&fine)?;
    // ‖u_{h/2}′ − u_h′‖ over the fine mesh.
    let mut acc = 0.0;
    for e in 0..fine.mesh().len() - 1 {
        let (x0, x1) = (fine.mesh()[e], fine.mesh()[e + 1]);
        let h = x1 - x0;
        let fine_slope = (fine_u[e + 1] - fine_u[e]) / h;
        // The coarse element containing this fine element.
        let coarse_e = e / 2;
        let ch = problem.mesh()[coarse_e + 1] - problem.mesh()[coarse_e];
        let coarse_slope = (candidate[coarse_e + 1] - candidate[coarse_e]) / ch;
        let d = fine_slope - coarse_slope;
        let updated = (h * d).mul_add(d, acc);
        if !fine_slope.is_finite()
            || !coarse_slope.is_finite()
            || !d.is_finite()
            || !updated.is_finite()
        {
            return Err(Fem1dError::NonFiniteIntermediate {
                stage: "hierarchical estimate",
                index: Some(e),
            });
        }
        acc = updated;
    }
    let estimate = acc.sqrt();
    if estimate.is_finite() {
        Ok(estimate)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "hierarchical estimate",
            index: None,
        })
    }
}

/// The nonlinear WARM-START fallback: the candidate is accepted only
/// as a starting point; the measured value is iteration savings and
/// the color is ESTIMATED, never verified (the honest R1 boundary).
#[derive(Debug, Clone)]
pub struct WarmStartReport {
    /// Newton iterations from a cold start (zero).
    pub cold_iterations: u32,
    /// Newton iterations from the candidate.
    pub warm_iterations: u32,
    /// The color of the claim (always `Estimated`).
    pub color: Color,
}

/// Measure warm-start savings on the toy nonlinear class.
///
/// # Errors
/// Returns [`Fem1dError`] when either run is malformed, unusable, or does not
/// converge within the admitted budget. Nonconvergence never becomes savings.
pub fn warm_start(
    problem: &MmsProblem,
    candidate: &[f64],
    max_iter: u32,
) -> Result<WarmStartReport, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let zero = try_zeroed("cold nonlinear start", problem.mesh().len())?;
    let cold = crate::fem1d::solve_nonlinear(problem, &zero, max_iter)?;
    require_converged(&cold, "cold nonlinear solve")?;
    let warm = crate::fem1d::solve_nonlinear(problem, candidate, max_iter)?;
    require_converged(&warm, "warm nonlinear solve")?;
    Ok(WarmStartReport {
        cold_iterations: cold.iterations,
        warm_iterations: warm.iterations,
        color: Color::Estimated {
            estimator: "warm-start-iteration-savings".to_string(),
            dispersion: f64::INFINITY,
        },
    })
}

/// Convenience for the batteries: effectivity of a report against the
/// oracle.
///
/// # Errors
/// Returns [`Fem1dError`] when the independent oracle or report bound is not a
/// usable finite value. Oracle failure is never mapped to effectivity `1.0`.
pub fn effectivity(
    problem: &MmsProblem,
    candidate: &[f64],
    report: &VerifierReport,
) -> Result<f64, Fem1dError> {
    if report.refusal.is_some() {
        return Err(Fem1dError::InvalidScalar {
            field: "verifier report",
            reason: "refused reports have no defined effectivity",
        });
    }
    let truth = true_energy_error(problem, candidate)?;
    if !report.bound.hi.is_finite() || report.bound.hi < 0.0 {
        return Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity report bound",
            index: None,
        });
    }
    if truth == 0.0 {
        return Err(Fem1dError::InvalidScalar {
            field: "oracle true error",
            reason: "effectivity is undefined for a zero denominator",
        });
    }
    let value = report.bound.hi / truth;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity ratio",
            index: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::dd::Dd;

    fn contains_dd(interval: Iv, exact: Dd) -> bool {
        !exact.lt(Dd::from_f64(interval.lo)) && !Dd::from_f64(interval.hi).lt(exact)
    }

    #[test]
    fn intervalized_element_inputs_cover_independent_dd_oracle() {
        // These decimal-looking f64 inputs make every legacy point computation
        // below round away a nonzero residual. The double-double oracle therefore
        // detects removal of intervalized mesh, slope, node, or weight arithmetic.
        let (x0, x1) = (0.1, 0.4);
        let (candidate0, candidate1) = (0.1, 0.2);
        let (node_constant, weight_constant) = GAUSS5_REF[0];
        let (dx0, dx1) = (Dd::from_f64(x0), Dd::from_f64(x1));
        let half_constant = Dd::from_f64(0.5);

        let (h, midpoint, half) = interval_element_geometry(x0, x1).unwrap();
        let exact_h = dx1 - dx0;
        let exact_midpoint = (dx0 + dx1) * half_constant;
        let exact_half = exact_h * half_constant;
        assert_ne!(exact_h, Dd::from_f64(x1 - x0));
        assert_ne!(exact_midpoint, Dd::from_f64(f64::midpoint(x0, x1)));
        assert_ne!(exact_half, Dd::from_f64((x1 - x0) * 0.5));
        assert!(contains_dd(h, exact_h));
        assert!(contains_dd(midpoint, exact_midpoint));
        assert!(contains_dd(half, exact_half));

        let slope = interval_candidate_slope(candidate0, candidate1, h).unwrap();
        let exact_difference = Dd::from_f64(candidate1) - Dd::from_f64(candidate0);
        let rounded_slope = (candidate1 - candidate0) / (x1 - x0);
        assert_ne!(Dd::from_f64(rounded_slope) * exact_h, exact_difference);
        assert!(!(exact_difference).lt(Dd::from_f64(slope.lo) * exact_h));
        assert!(!(Dd::from_f64(slope.hi) * exact_h).lt(exact_difference));

        let (node, weight) =
            interval_quadrature_geometry(midpoint, half, node_constant, weight_constant).unwrap();
        let exact_node = exact_midpoint + exact_half * Dd::from_f64(node_constant);
        let exact_weight = exact_half * Dd::from_f64(weight_constant);
        let rounded_node = f64::midpoint(x0, x1) + (x1 - x0) * 0.5 * node_constant;
        let rounded_weight = (x1 - x0) * 0.5 * weight_constant;
        assert_ne!(exact_node, Dd::from_f64(rounded_node));
        assert_ne!(exact_weight, Dd::from_f64(rounded_weight));
        assert!(contains_dd(node, exact_node));
        assert!(contains_dd(weight, exact_weight));

        // `1/3` is not representable. The coefficient interval must reach the
        // side of the rounded quotient selected by the exact FMA residual;
        // treating the rounded antiderivative coefficient as a point fails it.
        let coefficient = interval_antiderivative_coefficient(1.0, 3).unwrap();
        let rounded = 1.0_f64 / 3.0;
        let residual = rounded.mul_add(3.0, -1.0);
        if residual > 0.0 {
            assert!(coefficient.lo <= crate::interval::down(rounded));
        } else if residual < 0.0 {
            assert!(coefficient.hi >= crate::interval::up(rounded));
        } else {
            assert!(coefficient.lo <= rounded && rounded <= coefficient.hi);
        }
    }

    #[test]
    fn gauss_constants_enclose_independent_truth_brackets() {
        // Each bit pair is the adjacent-f64 bracket around the corresponding
        // high-precision Gauss-Legendre constant, derived independently from
        // the decimal reference values. Fifteen-digit literals miss some
        // weights by up to eight ulps, so this locks the certified quadrature
        // inputs rather than merely checking that `iv_c` widens its input.
        let positive_constants = [
            (
                GAUSS5_REF[4].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[3].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
            (
                GAUSS5_REF[0].1,
                0x3fce_539e_c36e_038c,
                0x3fce_539e_c36e_038d,
            ),
            (
                GAUSS5_REF[1].1,
                0x3fde_a1da_25ae_415a,
                0x3fde_a1da_25ae_415b,
            ),
            (
                GAUSS5_REF[2].1,
                0x3fe2_3456_789a_bcdf,
                0x3fe2_3456_789a_bce0,
            ),
        ];
        for (constant, lower_bits, upper_bits) in positive_constants {
            let interval = iv_c(constant);
            assert!(interval.lo <= f64::from_bits(lower_bits));
            assert!(interval.hi >= f64::from_bits(upper_bits));
        }

        for (constant, positive_lower_bits, positive_upper_bits) in [
            (
                GAUSS5_REF[0].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[1].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
        ] {
            let interval = iv_c(constant);
            assert!(interval.lo <= -f64::from_bits(positive_upper_bits));
            assert!(interval.hi >= -f64::from_bits(positive_lower_bits));
        }
    }
}
