//! Outward-evaluated equilibrated RT0 majorants on a conforming tetrahedral domain.
//!
//! The admitted PDE is -div(k grad u)=f, with element-constant scalar k>0
//! and f, P1 Dirichlet data, constant outward Neumann flux, and constant
//! positive Robin h with P1 reference data. Coordinates and coefficients are
//! interpreted as exact real values of their binary64 encodings.
//!
//! A graph solve only proposes face fluxes. A rooted forest then defines an
//! EXACT conservative flux by real-arithmetic elimination, enclosed by `Iv`:
//! every cell's outgoing flux sums to f*volume and shared faces use one flux
//! with opposite signs. No small floating-point residual is called equilibrium.
//! For q in H(div), div q=f and q.n=g_N, integration by parts gives
//!   ||u-v||_a <= (||q+k grad v||^2_(1/k)
//!                  + ||q.n-h(v-u_ref)||^2_(1/h,Robin))^(1/2).
//! All integrands are quadratic polynomials and simplex moments are evaluated
//! outward, including geometry. This includes algebraic error in ANY admitted
//! P1 candidate, not just a converged Galerkin solution.
//!
//! Domain-conditional: the caller must supply a conforming, non-overlapping
//! tetrahedralization of its declared polyhedral domain. Local incidence,
//! face orientation and nondegeneracy are checked; global intersection freedom
//! and agreement with an original CAD surface are NOT certified here. No
//! nonlinear, contact, uncertain-coefficient, or point-maximum bound is claimed.
//! The construction does not claim local efficiency or contrast robustness.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::interval::Iv;

mod geometry;
mod goal;
pub use goal::{GoalBound, MeanBound, goal_bound, mean_bound};
use geometry::{Cell, Face, build, dot, integral_square, scale, sub};

/// Data on one exterior face, ordered as `BoundaryFace::vertices`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Piecewise affine prescribed trace. The candidate must match exactly.
    Dirichlet([f64; 3]),
    /// Constant outward flux: positive means heat leaving the domain.
    Neumann(f64),
    /// q.n = h*(u-reference), with constant h>0 and affine reference.
    Robin { h: f64, reference: [f64; 3] },
}

/// Exactly one declaration is required for every exterior face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryFace {
    pub vertices: [usize; 3],
    pub condition: BoundaryCondition,
}

/// Borrowed linear scalar problem. Source and conductivity have one value per tet.
#[derive(Debug, Clone, Copy)]
pub struct TetProblem<'a> {
    pub vertices: &'a [[f64; 3]],
    pub tets: &'a [[usize; 4]],
    pub conductivity: &'a [f64],
    pub source: &'a [f64],
    pub boundary: &'a [BoundaryFace],
}

/// Admission and work limits. CG is a bounded proposal, not a proof obligation.
#[derive(Debug, Clone, Copy)]
pub struct FluxBudget {
    pub max_cells: usize,
    pub max_iterations: usize,
}

impl Default for FluxBudget {
    fn default() -> Self {
        Self { max_cells: 65_536, max_iterations: 128 }
    }
}

/// No bound is returned on a refusal or cancellation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TetError {
    Invalid(&'static str),
    Unsupported(&'static str),
    Unbounded,
    Budget,
    Cancelled,
}
impl std::fmt::Display for TetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "tetrahedral flux bound refused: {self:?}")
    }
}
impl std::error::Error for TetError {}

/// Bound on sqrt(integral k|grad(u-v)|^2 + integral_R h(u-v)^2).
/// Does not include geometry/model/input uncertainty and is not a maximum bound.
#[derive(Debug, Clone)]
pub struct EnergyBound {
    pub energy_error_upper: f64,
    /// Enclosure of the computable majorant squared, NOT of the true error squared.
    pub majorant_squared: Iv,
    /// Outward upper contributions, including a cell's exterior Robin faces.
    pub cell_majorant_squared_upper: Vec<f64>,
    /// Correlated conservative face-flux enclosures, outward from each cell.
    /// Entry i is opposite the tet's vertex i. Independent choices from these
    /// boxes need not be conservative; the forest construction establishes existence.
    pub outward_flux_integrals: Vec<[Iv; 4]>,
    pub proposal_iterations: usize,
}

fn poll(keep_going: &mut impl FnMut() -> bool) -> Result<(), TetError> {
    if keep_going() { Ok(()) } else { Err(TetError::Cancelled) }
}
fn midpoint(x: Iv) -> Result<f64, TetError> {
    if x.is_unbounded() { return Err(TetError::Unbounded); }
    let m = 0.5 * x.lo + 0.5 * x.hi;
    if m.is_finite() { Ok(m) } else { Err(TetError::Unbounded) }
}
fn signed(x: Iv, sign: f64) -> Iv { x.mul(Iv::point(sign)) }

/// Construct and integrate an equilibrated flux for an arbitrary conforming P1 field.
/// Cancellation is polled during geometry, proposal iterations, elimination and integration.
pub fn energy_bound(
    problem: &TetProblem<'_>,
    candidate: &[f64],
    budget: FluxBudget,
    mut keep_going: impl FnMut() -> bool,
) -> Result<EnergyBound, TetError> {
    let (cells, faces) = build(problem, candidate, budget, &mut keep_going)?;
    let (mut flux, proposal_iterations) = propose(problem, &cells, &faces, budget, &mut keep_going)?;
    equilibrate(problem, &cells, &faces, &mut flux, &mut keep_going)?;
    let mut total = Iv::zero();
    let mut contributions = Vec::with_capacity(cells.len());
    let mut outward = Vec::with_capacity(cells.len());
    for (e, cell) in cells.iter().enumerate() {
        poll(&mut keep_going)?;
        let local = std::array::from_fn(|i| {
            let face = &faces[cell.faces[i]];
            signed(flux[cell.faces[i]], face.sign(e))
        });
        let mut defect = [[Iv::zero(); 3]; 4];
        for (j, value) in defect.iter_mut().enumerate() {
            *value = scale(cell.gradient, Iv::point(problem.conductivity[e]));
            for (i, integral) in local.iter().enumerate() {
                let basis = scale(sub(cell.points[j], cell.points[i]),
                    integral.div_pos(cell.volume.scale_pos(3.0)));
                for d in 0..3 { value[d] = value[d].add(basis[d]); }
            }
        }
        let squared: Vec<Iv> = (0..3).map(|d| {
            integral_square(&defect.map(|v| v[d]), cell.volume, 20.0)
        }).collect();
        let mut eta = squared.into_iter().fold(Iv::zero(), Iv::add)
            .div_pos(Iv::point(problem.conductivity[e]));
        for (i, &f) in cell.faces.iter().enumerate() {
            if let Some(BoundaryCondition::Robin { h, reference }) = faces[f].condition {
                let qn = local[i].div_pos(faces[f].area);
                let residual = std::array::from_fn::<_, 3, _>(|j| {
                    qn.sub(Iv::point(h).mul(Iv::point(candidate[faces[f].vertices[j]])
                        .sub(Iv::point(reference[j]))))
                });
                eta = eta.add(integral_square(&residual, faces[f].area, 12.0)
                    .div_pos(Iv::point(h)));
            }
        }
        if eta.is_unbounded() || eta.hi < 0.0 { return Err(TetError::Unbounded); }
        contributions.push(eta.hi);
        total = total.add(eta);
        outward.push(local);
    }
    poll(&mut keep_going)?;
    let root = total.sqrt();
    if root.is_unbounded() { return Err(TetError::Unbounded); }
    Ok(EnergyBound {
        energy_error_upper: root.hi,
        majorant_squared: total,
        cell_majorant_squared_upper: contributions,
        outward_flux_integrals: outward,
        proposal_iterations,
    })
}

/// Weighted graph projection improves the proposal. It is intentionally not
/// trusted: early CG termination is safe because exact forest repair follows.
fn propose(
    problem: &TetProblem<'_>, cells: &[Cell], faces: &[Face], budget: FluxBudget,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<(Vec<Iv>, usize), TetError> {
    let n = cells.len();
    let mut flux = Vec::with_capacity(faces.len());
    let mut weights = Vec::with_capacity(faces.len());
    let mut rhs = problem.source.iter().zip(cells)
        .map(|(&f, c)| midpoint(c.volume.mul(Iv::point(f)))).collect::<Result<Vec<_>, _>>()?;
    let mut diag = vec![0.0; n];
    for face in faces {
        poll(keep_going)?;
        let a = face.sides[0].0;
        let fixed = matches!(face.condition, Some(BoundaryCondition::Neumann(_)));
        let proposed = if let Some(BoundaryCondition::Neumann(q)) = face.condition {
            face.area.mul(Iv::point(q))
        } else {
            let mut q = -problem.conductivity[a] * midpoint(dot(cells[a].gradient, face.normal_area))?;
            if let Some(&(b, _)) = face.sides.get(1) {
                q = 0.5 * q - 0.5 * problem.conductivity[b] * midpoint(dot(cells[b].gradient, face.normal_area))?;
            }
            Iv::point(q)
        };
        let f = midpoint(proposed)?;
        rhs[a] -= f;
        if let Some(&(b, _)) = face.sides.get(1) { rhs[b] += f; }
        let resistance = face.sides.iter().map(|&(e, _)| {
            midpoint(cells[e].volume).map(|v| v / problem.conductivity[e])
        }).collect::<Result<Vec<_>, _>>()?.into_iter().sum::<f64>();
        let area = midpoint(face.area)?;
        let weight = if fixed { 0.0 } else { area * area / resistance };
        if !weight.is_finite() || (!fixed && weight <= 0.0) || !f.is_finite() {
            return Err(TetError::Unbounded);
        }
        diag[a] += weight;
        if let Some(&(b, _)) = face.sides.get(1) { diag[b] += weight; }
        flux.push(proposed);
        weights.push(weight);
    }
    if diag.iter().any(|d| !d.is_finite() || *d <= 0.0) {
        return Err(TetError::Unsupported("each component needs a Dirichlet or positive Robin boundary"));
    }
    let apply = |x: &[f64]| {
        let mut y = vec![0.0; n];
        for (face, &w) in faces.iter().zip(&weights) {
            let a = face.sides[0].0;
            let b = face.sides.get(1).map(|s| s.0);
            let value = w * (x[a] - b.map_or(0.0, |b| x[b]));
            y[a] += value;
            if let Some(b) = b { y[b] -= value; }
        }
        y
    };
    let scalar = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(a, b)| a*b).sum::<f64>();
    let mut x = vec![0.0; n];
    let mut r = rhs;
    let mut p: Vec<_> = r.iter().zip(&diag).map(|(r, d)| r / d).collect();
    let mut rz = scalar(&r, &p);
    let initial = rz;
    let mut iterations = 0;
    for _ in 0..budget.max_iterations {
        poll(keep_going)?;
        if !rz.is_finite() || rz <= 0.0 || rz <= initial * 1e-24 { break; }
        let ap = apply(&p);
        let pap = scalar(&p, &ap);
        if !pap.is_finite() || pap <= 0.0 { break; }
        let alpha = rz / pap;
        if !alpha.is_finite() { break; }
        for i in 0..n { x[i] += alpha*p[i]; r[i] -= alpha*ap[i]; }
        let z: Vec<_> = r.iter().zip(&diag).map(|(r, d)| r/d).collect();
        let next = scalar(&r, &z);
        let beta = next / rz;
        for i in 0..n { p[i] = z[i] + beta*p[i]; }
        rz = next;
        iterations += 1;
    }
    for ((face, &w), value) in faces.iter().zip(&weights).zip(&mut flux) {
        poll(keep_going)?;
        if w == 0.0 { continue; }
        let a = face.sides[0].0;
        let delta = w * (x[a] - face.sides.get(1).map_or(0.0, |s| x[s.0]));
        let proposed = midpoint(*value)? + delta;
        if !proposed.is_finite() { return Err(TetError::Unbounded); }
        *value = Iv::point(proposed);
    }
    Ok((flux, iterations))
}

fn equilibrate(
    problem: &TetProblem<'_>, cells: &[Cell], faces: &[Face], flux: &mut [Iv],
    keep_going: &mut impl FnMut() -> bool,
) -> Result<(), TetError> {
    let mut parent = vec![None; cells.len()];
    let mut order = Vec::with_capacity(cells.len());
    let mut queue = VecDeque::new();
    // Multi-root forest: each cell adjacent to a free exterior face can drain
    // directly. This shortens correction paths without trusting a PDE solve.
    for (f, face) in faces.iter().enumerate() {
        poll(keep_going)?;
        if face.sides.len() == 1 && !matches!(face.condition, Some(BoundaryCondition::Neumann(_))) {
            let e = face.sides[0].0;
            if parent[e].is_none() { parent[e] = Some(f); queue.push_back(e); }
        }
    }
    while let Some(e) = queue.pop_front() {
        poll(keep_going)?;
        order.push(e);
        for &f in &cells[e].faces {
            for &(other, _) in &faces[f].sides {
                if parent[other].is_none() {
                    parent[other] = Some(f);
                    queue.push_back(other);
                }
            }
        }
    }
    if order.len() != cells.len() {
        return Err(TetError::Unsupported("unanchored pure-Neumann component"));
    }
    for &e in order.iter().rev() {
        poll(keep_going)?;
        let pf = parent[e].ok_or(TetError::Invalid("missing forest parent"))?;
        let mut balance = cells[e].volume.mul(Iv::point(problem.source[e]));
        for &f in &cells[e].faces {
            if f != pf { balance = balance.sub(signed(flux[f], faces[f].sign(e))); }
        }
        flux[pf] = signed(balance, faces[pf].sign(e));
        if flux[pf].is_unbounded() { return Err(TetError::Unbounded); }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
