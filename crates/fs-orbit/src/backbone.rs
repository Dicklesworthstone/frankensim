//! Conservative backbone charts in the complete original physical state.
//!
//! Fixed-amplitude entries permute the full problem into the existing solver's
//! first-coordinate chart without changing its equations or refinement owner.
//! Family continuation retains the conservation-dependent balance row and
//! uses a rectangular phase/arclength corrector over all states and frequency.

pub(crate) mod continuation;
pub use continuation::{
    AdaptiveBackboneOrbit, BackboneBranch, BackboneContinuationBudget, BackboneContinuationError,
    BackboneContinuationPoint, BackboneRefinementControls, continue_backbone_branch,
    solve_backbone_adaptive,
};

use core::cell::RefCell;

use crate::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbRefinementBudget, HbScales, IslandStateJacobian,
    OrbitError, OrbitProblem, ScaledHbOrbit, scaling::AdmittedScales,
    solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded,
};
use fs_math::c64::C64;

/// Resolve a backbone with the first-harmonic peak of a chosen state pinned.
///
/// `coordinate` is an index in the ORIGINAL complete state. `amplitude` is
/// `2*abs(X_1[coordinate])` in that coordinate's units, not its turning-point
/// displacement or full-waveform peak. The physical packed seed includes every
/// original state and harmonic followed by the unknown angular frequency.
/// Output coefficients keep that same original state ordering.
///
/// This changes only the chart used by the existing scaled backbone solver.
/// All island equations, frequency-dependent port rows, declared exact state
/// tangents, physical scales and independent sampling checks remain active.
/// The phase remains free, as in [`HbAnchor::Backbone`]. A static or nodal first
/// state is admissible when the selected coordinate oscillates. No damping,
/// coupling or material-memory state is removed to obtain a conservative orbit.
///
/// The source must declare itself autonomous. A nonzero amplitude does not
/// establish a periodic solution, stability, a primitive period or convergence
/// of harmonic truncation. Use the adaptive entry for measured order refinement.
///
/// Apart from the solver's existing work, preparation retains one state buffer
/// and permuted scales; island calls allocate nothing, and each port matrix is
/// permuted in place. Selecting coordinate zero calls the original owner directly.
///
/// # Errors
/// An invalid coordinate, nonautonomous source, invalid seed/scales/budgets, or
/// a refusal from the original fixed-order solver. Admission precedes source
/// evaluation; a declared exact tangent's refusal never becomes finite differences.
pub fn solve_hb_backbone_coordinate_scaled_seeded<P: OrbitProblem>(
    problem: &P,
    coordinate: usize,
    amplitude: f64,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, OrbitError> {
    let anchor = admit(problem, coordinate, amplitude, &seed, budget, scales)?;
    if coordinate == 0 {
        return solve_hb_scaled_seeded(problem, anchor, seed, budget, scales);
    }
    let mapped = SwappedProblem::new(problem, coordinate);
    let mut result = solve_hb_scaled_seeded(
        &mapped,
        anchor,
        swap_seed(seed, problem.dim(), coordinate),
        budget,
        &swap_scales(scales, coordinate),
    )?;
    result.orbit.coeffs.swap(0, coordinate);
    Ok(result)
}

/// The same physical-coordinate backbone with automatic harmonic refinement.
///
/// Every larger solve pins the SAME original coordinate and amplitude. State
/// and equation tolerances are supplied in original order and units. The shared
/// refinement owner still requires independent sampling, retained/omitted
/// equation agreement, waveform agreement between two solved orders, and
/// frequency agreement. Its free-phase waveform alignment uses the selected
/// coordinate rather than a possibly static first state. Returned coefficients
/// are restored to original order; the scalar refinement evidence is retained.
/// These measurements are not an infinite-tail or forward-error certificate.
///
/// # Errors
/// The fixed-order entry's refusals and the existing adaptive solver's physical
/// resolution and resource admission failures.
pub fn solve_hb_backbone_coordinate_adaptive_scaled_seeded<P: OrbitProblem>(
    problem: &P,
    coordinate: usize,
    amplitude: f64,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
) -> Result<AdaptiveHbOrbit, OrbitError> {
    let anchor = admit(problem, coordinate, amplitude, &seed, budget, scales)?;
    if refinement.residual_tolerances.len() != problem.dim()
        || refinement.state_tolerances.len() != problem.dim()
    {
        return Err(invalid(
            "backbone refinement needs original state and equation tolerances",
        ));
    }
    if coordinate == 0 {
        return solve_hb_adaptive_scaled_seeded(problem, anchor, seed, budget, refinement, scales);
    }
    let mapped = SwappedProblem::new(problem, coordinate);
    let mut limits = refinement.clone();
    limits.residual_tolerances.swap(0, coordinate);
    limits.state_tolerances.swap(0, coordinate);
    let mut result = solve_hb_adaptive_scaled_seeded(
        &mapped,
        anchor,
        swap_seed(seed, problem.dim(), coordinate),
        budget,
        &limits,
        &swap_scales(scales, coordinate),
    )?;
    result.orbit.coeffs.swap(0, coordinate);
    Ok(result)
}

fn invalid(what: &'static str) -> OrbitError {
    OrbitError::BadParameter { what }
}

fn admit<P: OrbitProblem>(
    problem: &P,
    coordinate: usize,
    amplitude: f64,
    seed: &[f64],
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<HbAnchor, OrbitError> {
    if coordinate >= problem.dim() || !problem.autonomous() {
        return Err(invalid(
            "backbone needs an original coordinate of an autonomous problem",
        ));
    }
    let anchor = HbAnchor::Backbone {
        amplitude,
        omega_guess: seed.last().copied().unwrap_or(f64::NAN),
    };
    AdmittedScales::new(problem.dim(), anchor, seed, budget, scales)?;
    Ok(anchor)
}

fn swap_seed(mut seed: Vec<f64>, dimension: usize, coordinate: usize) -> Vec<f64> {
    // Every DC/real/imaginary block contains one complete state. The final
    // frequency is a scalar, not a state, and must retain its position/value.
    let coefficient_count = seed.len() - 1;
    for block in seed[..coefficient_count].chunks_exact_mut(dimension) {
        block.swap(0, coordinate);
    }
    seed
}

fn swap_scales(scales: &HbScales, coordinate: usize) -> HbScales {
    let mut result = scales.clone();
    result.states.swap(0, coordinate);
    result.equations.swap(0, coordinate);
    result
}

fn swap_matrix<T>(matrix: &mut [T], dimension: usize, coordinate: usize) {
    for j in 0..dimension {
        matrix.swap(j, coordinate * dimension + j);
    }
    for i in 0..dimension {
        matrix.swap(i * dimension, i * dimension + coordinate);
    }
}

struct SwappedProblem<'a, P> {
    original: &'a P,
    coordinate: usize,
    state: RefCell<Vec<f64>>,
    jacobian: Option<IslandStateJacobian<P>>,
}

impl<'a, P: OrbitProblem> SwappedProblem<'a, P> {
    fn new(original: &'a P, coordinate: usize) -> Self {
        Self {
            original,
            coordinate,
            state: RefCell::new(vec![0.0; original.dim()]),
            jacobian: original.island_state_jacobian(),
        }
    }

    fn state_jacobian(&self, t: f64, x: &[f64], out: &mut [f64]) -> bool {
        let mut state = self.state.borrow_mut();
        state.copy_from_slice(x);
        state.swap(0, self.coordinate);
        let callback = self.jacobian.expect("declared original tangent");
        if !callback(self.original, t, &state, out) {
            return false;
        }
        swap_matrix(out, self.dim(), self.coordinate);
        true
    }
}

impl<P: OrbitProblem> OrbitProblem for SwappedProblem<'_, P> {
    fn dim(&self) -> usize {
        self.original.dim()
    }

    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let mut state = self.state.borrow_mut();
        state.copy_from_slice(x);
        state.swap(0, self.coordinate);
        self.original.island(t, &state, out);
        out.swap(0, self.coordinate);
    }

    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.jacobian
            .map(|_| Self::state_jacobian as IslandStateJacobian<Self>)
    }

    fn port(&self, s: C64) -> Vec<C64> {
        let mut matrix = self.original.port(s);
        let dimension = self.dim();
        if matrix.len() != dimension * dimension {
            return vec![C64::new(f64::NAN, f64::NAN); dimension * dimension];
        }
        swap_matrix(&mut matrix, dimension, self.coordinate);
        matrix
    }

    fn autonomous(&self) -> bool {
        self.original.autonomous()
    }
}
