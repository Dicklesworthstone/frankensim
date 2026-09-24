//! Primal/dual continuum bounds for cell-weighted integrals, not point samples.
//! For a(z,w)=J(w), J(u)-J(v)=R_v(z_h)+a(u-v,z-z_h).
//! The residual correction is evaluated outward; neither candidate is assumed
//! to solve its discrete equations. Cauchy-Schwarz bounds the last term by the
//! product of the independently equilibrated primal and dual energy bounds.
use super::*;

#[derive(Debug, Clone)]
pub struct GoalBound {
    /// Continuum enclosure of integral w*u over the declared domain.
    pub enclosure: Iv,
    /// Integral of the supplied P1 primal field, evaluated outward.
    pub candidate_value: Iv,
    /// l(z_h)-a(v,z_h), including Neumann and Robin loads.
    pub residual_correction: Iv,
    /// Upper bound on the absolute primal/dual error product.
    pub remainder_upper: f64,
    pub domain_volume: Iv,
    pub primal: EnergyBound,
    pub dual: EnergyBound,
}

#[derive(Debug, Clone)]
pub struct MeanBound {
    /// Whole-domain volume average, divided by the OUTWARD geometric volume.
    pub enclosure: Iv,
    pub candidate_mean: Iv,
    /// Underlying unit-weight integral, with both flux reconstructions.
    pub integral: GoalBound,
}

fn linear_integral(values: &[Iv], measure: Iv) -> Iv {
    measure.mul(values.iter().copied().fold(Iv::zero(), Iv::add))
        .div_pos(Iv::point(values.len() as f64))
}
fn product_integral(a: &[Iv], b: &[Iv], measure: Iv) -> Iv {
    let sa = a.iter().copied().fold(Iv::zero(), Iv::add);
    let sb = b.iter().copied().fold(Iv::zero(), Iv::add);
    let diagonal = a.iter().zip(b).fold(Iv::zero(), |s, (a,b)| s.add(a.mul(*b)));
    measure.mul(diagonal.add(sa.mul(sb)))
        .div_pos(Iv::point((a.len()*(a.len()+1)) as f64))
}

/// Bound a cell-weighted volume integral with a conforming homogeneous-Dirichlet
/// dual candidate. `weights` declares one exact constant per tetrahedron. The
/// dual equation is constructed here from that same functional and primal
/// operator; callers cannot substitute a different PDE or Robin coefficient.
pub fn goal_bound(
    problem: &TetProblem<'_>, candidate: &[f64], dual_candidate: &[f64],
    weights: &[f64], budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
) -> Result<GoalBound, TetError> {
    goal_bound_impl(&Problem::from(problem), candidate, dual_candidate, weights, budget, &mut keep_going)
}

pub(super) fn goal_bound_impl(
    problem: &Problem<'_>, candidate: &[f64], dual_candidate: &[f64],
    weights: &[f64], budget: FluxBudget, keep_going: &mut impl FnMut() -> bool,
) -> Result<GoalBound, TetError> {
    poll(keep_going)?;
    if problem.tets.len() > budget.max_cells { return Err(TetError::Budget); }
    if weights.len() != problem.tets.len() || weights.iter().any(|w| !w.is_finite()) {
        return Err(TetError::Invalid("finite cell weights required for the exact functional"));
    }
    let primal = energy_bound_impl(problem, candidate, budget, keep_going)?;
    let mut dual_boundary = Vec::with_capacity(problem.boundary.len());
    for face in problem.boundary {
        poll(keep_going)?;
        dual_boundary.push(BoundaryFace { vertices: face.vertices, condition: match face.condition {
            BoundaryCondition::Dirichlet(_) => BoundaryCondition::Dirichlet([0.0; 3]),
            BoundaryCondition::Neumann(_) => BoundaryCondition::Neumann(0.0),
            BoundaryCondition::Robin { h, .. } => BoundaryCondition::Robin { h, reference: [0.0; 3] },
        }});
    }
    let dual_problem = Problem { source: weights, boundary: &dual_boundary, ..*problem };
    let dual = energy_bound_impl(&dual_problem, dual_candidate, budget, keep_going)?;
    let (cells, faces) = build(problem, candidate, budget, keep_going)?;
    let (dual_cells, _) = build(&dual_problem, dual_candidate, budget, keep_going)?;
    let mut value = Iv::zero();
    let mut residual = Iv::zero();
    let mut volume = Iv::zero();
    for (e, cell) in cells.iter().enumerate() {
        poll(keep_going)?;
        let u = problem.tets[e].map(|i| Iv::point(candidate[i]));
        let z = problem.tets[e].map(|i| Iv::point(dual_candidate[i]));
        volume = volume.add(cell.volume);
        value = value.add(linear_integral(&u, cell.volume).mul(Iv::point(weights[e])));
        residual = residual.add(linear_integral(&z, cell.volume).mul(Iv::point(problem.source[e])))
            .sub(problem.conductivity.bilinear_integral(e, cell.gradient, dual_cells[e].gradient, cell.volume));
    }
    for face in &faces {
        poll(keep_going)?;
        let z = face.vertices.map(|i| Iv::point(dual_candidate[i]));
        match face.condition {
            Some(BoundaryCondition::Neumann(q)) => {
                residual = residual.sub(linear_integral(&z, face.area).mul(Iv::point(q)));
            }
            Some(BoundaryCondition::Robin { h, reference }) => {
                let departure = std::array::from_fn::<_, 3, _>(|i| {
                    Iv::point(candidate[face.vertices[i]]).sub(Iv::point(reference[i]))
                });
                residual = residual.sub(product_integral(&departure, &z, face.area).mul(Iv::point(h)));
            }
            _ => {}
        }
    }
    let remainder = Iv::point(primal.energy_error_upper).mul(Iv::point(dual.energy_error_upper));
    let enclosure = value.add(residual).add(Iv { lo: -remainder.hi, hi: remainder.hi });
    poll(keep_going)?;
    if enclosure.is_unbounded() || volume.is_unbounded() || volume.lo <= 0.0 {
        return Err(TetError::Unbounded);
    }
    Ok(GoalBound { enclosure, candidate_value: value, residual_correction: residual,
        remainder_upper: remainder.hi, domain_volume: volume, primal, dual })
}

/// Enclose the whole-domain mean temperature. The dual candidate belongs to
/// the unit volumetric-source equation, not to a nodal maximum or a rounded
/// reciprocal-volume source. Division by exact-domain volume happens last.
pub fn mean_bound(
    problem: &TetProblem<'_>, candidate: &[f64], dual_candidate: &[f64],
    budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
) -> Result<MeanBound, TetError> {
    mean_bound_impl(&Problem::from(problem), candidate, dual_candidate, budget, &mut keep_going)
}

pub(super) fn mean_bound_impl(
    problem: &Problem<'_>, candidate: &[f64], dual_candidate: &[f64],
    budget: FluxBudget, keep_going: &mut impl FnMut() -> bool,
) -> Result<MeanBound, TetError> {
    poll(keep_going)?;
    if problem.tets.len() > budget.max_cells { return Err(TetError::Budget); }
    let weights = vec![1.0; problem.tets.len()];
    let integral = goal_bound_impl(problem, candidate, dual_candidate, &weights, budget, keep_going)?;
    let enclosure = integral.enclosure.div_pos(integral.domain_volume);
    let candidate_mean = integral.candidate_value.div_pos(integral.domain_volume);
    if enclosure.is_unbounded() || candidate_mean.is_unbounded() { return Err(TetError::Unbounded); }
    Ok(MeanBound { enclosure, candidate_mean, integral })
}
