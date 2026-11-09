//! Deterministic ring-Galerkin solution of the diffuse Clausing equation.
//!
//! With x=z/D, E(x)=(sqrt(1+x*x)-x)^2/(2*sqrt(1+x*x)) and
//! K(x)=-E'(x), w(x)=E(h-x)+integral K(|x-y|)w(y)dy.
//! Exact ring averages preserve reciprocity; CG solves the antisymmetric half
//! system for w-1/2 with O(N) storage. Refinement estimates spatial error but
//! is not a certified enclosure. Roundoff is excluded from the residual bound.
//! Equations and independent benchmarks: Matuda, DOI 10.3131/jvsj2.58.311,
//! and Gomez-Goni/Lobo, DOI 10.1116/1.1568746.

use super::super::super::poll;
use super::{CylindricalMolecularGeometry, GasEnergyError, det, invalid, positive};

pub(super) mod accommodation;

/// Explicit accuracy and work budgets. At least three nested grids are used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClausingControls {
    pub initial_rings: usize,
    pub maximum_rings: usize,
    /// Maximum CG steps per grid; cancellation is checked every sixteen rows.
    pub maximum_iterations: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
}

impl Default for ClausingControls {
    fn default() -> Self {
        Self {
            initial_rings: 32,
            maximum_rings: 2048,
            maximum_iterations: 512,
            absolute_tolerance: 1e-6,
            relative_tolerance: 1e-4,
        }
    }
}

impl ClausingControls {
    fn validate(self) -> Result<(), GasEnergyError> {
        if self.initial_rings < 4
            || self.initial_rings % 2 != 0
            || self.maximum_rings > 2048
            || self.initial_rings > self.maximum_rings / 4
            || self.maximum_iterations == 0
            || self.maximum_iterations > 4096
            || ![self.absolute_tolerance, self.relative_tolerance]
                .iter()
                .all(|x| x.is_finite() && *x >= 0.0)
            || !positive(self.absolute_tolerance + self.relative_tolerance)
        {
            return Err(invalid("invalid Clausing accuracy or work budgets"));
        }
        Ok(())
    }
}

/// Evidence for the returned probabilities, not an error certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClausingResolution {
    pub controls: ClausingControls,
    /// Zero for the exact zero-length aperture.
    pub rings: usize,
    /// Number of ring-count doublings.
    pub refinements: usize,
    /// Total CG steps over all grids.
    pub iterations: usize,
    /// Full-system L1 residual, explicitly recomputed on the final grid.
    pub residual_l1: f64,
    /// Latest absolute change in the wall-transmitted probability.
    pub discretization_estimate: f64,
    /// Discrete transmission residual bound, excluding floating-point rounding.
    pub linear_solve_estimate: f64,
    pub previous_linear_solve_estimate: f64,
    /// Grid change plus both linear-solve estimates; not a continuum enclosure.
    pub estimated_absolute_error: f64,
}

pub(super) fn resolve(
    mut geometry: CylindricalMolecularGeometry,
    controls: ClausingControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CylindricalMolecularGeometry, GasEnergyError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if geometry.length_m == 0.0 {
        geometry.resolution = Some(ClausingResolution {
            controls,
            rings: 0,
            refinements: 0,
            iterations: 0,
            residual_l1: 0.0,
            discretization_estimate: 0.0,
            linear_solve_estimate: 0.0,
            previous_linear_solve_estimate: 0.0,
            estimated_absolute_error: 0.0,
        });
        return Ok(geometry);
    }
    let h = geometry.length_m / geometry.diameter_m;
    let linear_target =
        0.01 * (controls.absolute_tolerance + controls.relative_tolerance * geometry.transmission);
    if !positive(linear_target) {
        return Err(invalid("Clausing solve tolerance is not representable"));
    }
    let mut rings = controls.initial_rings;
    let mut refinements = 0;
    let mut iterations = 0;
    let mut previous: Option<GridSolution> = None;
    loop {
        let grid = Grid::new(h, rings, &mut cancelled)?;
        let current = grid.solve(
            controls.maximum_iterations,
            linear_target,
            geometry.wall_hit,
            &mut cancelled,
        )?;
        iterations += current.iterations;
        if let Some(old) = previous {
            let change = (current.wall_transmitted - old.wall_transmitted).abs();
            let estimate = change + current.linear_estimate + old.linear_estimate;
            let transmission = geometry.direct + current.wall_transmitted;
            let tolerance =
                controls.absolute_tolerance + controls.relative_tolerance * transmission;
            if refinements >= 2 && estimate <= tolerance {
                if !positive(transmission) || transmission > 1.0 {
                    return Err(invalid(
                        "resolved Clausing transmission is not representable",
                    ));
                }
                poll(&mut cancelled)?;
                geometry.transmission = transmission;
                geometry.wall_transmitted = current.wall_transmitted;
                geometry.reflected = current.reflected;
                geometry.resolution = Some(ClausingResolution {
                    controls,
                    rings,
                    refinements,
                    iterations,
                    residual_l1: current.residual_l1,
                    discretization_estimate: change,
                    linear_solve_estimate: current.linear_estimate,
                    previous_linear_solve_estimate: old.linear_estimate,
                    estimated_absolute_error: estimate,
                });
                return Ok(geometry);
            }
        }
        previous = Some(current);
        if rings > controls.maximum_rings / 2 {
            return Err(invalid("Clausing spatial refinement budget exhausted"));
        }
        rings *= 2;
        refinements += 1;
    }
}

struct Grid {
    delta: f64,
    diagonal: f64,
    wall: Vec<f64>,
    end: Vec<f64>,
}

#[derive(Clone)]
struct GridSolution {
    wall_transmitted: f64,
    reflected: f64,
    residual_l1: f64,
    linear_estimate: f64,
    iterations: usize,
    half_escape: Vec<f64>,
}

impl Grid {
    fn new(
        h: f64,
        rings: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        let delta = h / rings as f64;
        if !positive(delta) {
            return Err(invalid("Clausing ring width is not representable"));
        }
        let mut roots = Vec::with_capacity(rings + 1);
        let mut q = Vec::with_capacity(rings + 1);
        for i in 0..=rings {
            if i % 16 == 0 {
                poll(cancelled)?;
            }
            let x = i as f64 * delta;
            let (root, reciprocal) = if x <= 1.0 {
                let root = det::sqrt(1.0 + x * x);
                (root, 1.0 / (root + x))
            } else {
                let inverse = 1.0 / x;
                let scaled = det::sqrt(1.0 + inverse * inverse);
                (x * scaled, inverse / (scaled + 1.0))
            };
            roots.push(root);
            q.push(reciprocal);
        }
        let mut end = Vec::with_capacity(rings);
        let mut wall = vec![0.0; rings];
        wall[0] = if delta <= 1.0 {
            delta * (1.0 - delta / (roots[1] + 1.0))
        } else {
            1.0 - q[1]
        };
        for i in 0..rings {
            if i % 16 == 0 {
                poll(cancelled)?;
            }
            let sum = q[i] + q[i + 1];
            end.push(sum * sum / (4.0 * (roots[i] + roots[i + 1])));
            if i > 0 {
                // Positive exact divided difference of the disc view factor;
                // avoids cancellation for both tiny rings and distant rings.
                let b = roots[i - 1] + roots[i];
                let c = roots[i] + roots[i + 1];
                let outer = q[i - 1] + q[i + 1];
                let right = q[i] + q[i + 1];
                wall[i] = delta / (2.0 * (roots[i - 1] + roots[i + 1]))
                    * (outer * (outer + 2.0 * q[i]) / b
                        + (2.0 * (i as f64 * delta) / b) * right * right / c);
            }
        }
        if !wall
            .iter()
            .chain(&end)
            .all(|x| x.is_finite() && *x >= 0.0 && *x <= 1.0)
        {
            return Err(invalid("Clausing view factors are not representable"));
        }
        Ok(Self {
            delta,
            diagonal: q[1],
            wall,
            end,
        })
    }

    fn apply(
        &self,
        input: &[f64],
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        let rings = self.wall.len();
        for i in 0..input.len() {
            if i % 16 == 0 {
                poll(cancelled)?;
            }
            output[i] = compensated(input.iter().enumerate().map(|(j, value)| {
                let entry = if i == j {
                    self.diagonal
                } else {
                    -self.wall[i.abs_diff(j)]
                };
                (entry + self.wall[rings - 1 - i - j]) * value
            }));
        }
        Ok(())
    }

    fn solve(
        &self,
        maximum_iterations: usize,
        target: f64,
        wall_hit: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GridSolution, GasEnergyError> {
        let rings = self.wall.len();
        let half = rings / 2;
        let rhs: Vec<_> = (0..half)
            .map(|i| 0.5 * (self.end[rings - 1 - i] - self.end[i]))
            .collect();
        let solved = conjugate_gradient(
            &rhs,
            maximum_iterations,
            target / (8.0 * self.delta),
            |input, output| self.apply(input, output, cancelled),
        )?;
        let solution = solved.values;
        let residual_l1 = 2.0 * solved.residual_l1;
        let linear_estimate = 4.0 * self.delta * residual_l1;
        if solution.iter().any(|x| !x.is_finite() || x.abs() > 0.5) {
            return Err(invalid(
                "Clausing exit probabilities left their physical bounds",
            ));
        }
        // Integrate both positive paths: 1-transmission loses tiny wall heat.
        let transmitted = 4.0
            * self.delta
            * compensated((0..half).map(|i| {
                self.end[i] * (0.5 + solution[i]) + self.end[rings - 1 - i] * (0.5 - solution[i])
            }));
        let reflected = 4.0
            * self.delta
            * compensated((0..half).map(|i| {
                self.end[i] * (0.5 - solution[i]) + self.end[rings - 1 - i] * (0.5 + solution[i])
            }));
        let normalization = wall_hit / (transmitted + reflected);
        if !positive(normalization)
            || (normalization - 1.0).abs() > 64.0 * rings as f64 * f64::EPSILON
        {
            return Err(invalid(
                "Clausing integrated wall probability balance failed",
            ));
        }
        let wall_transmitted = transmitted * normalization;
        let reflected = reflected * normalization;
        if ![wall_transmitted, reflected]
            .iter()
            .all(|x| positive(*x) && *x <= wall_hit)
        {
            return Err(invalid(
                "Clausing wall path probabilities are not representable",
            ));
        }
        Ok(GridSolution {
            wall_transmitted,
            reflected,
            residual_l1,
            linear_estimate,
            iterations: solved.iterations,
            half_escape: solution,
        })
    }
}

struct LinearSolution {
    values: Vec<f64>,
    residual_l1: f64,
    iterations: usize,
}

fn conjugate_gradient(
    rhs: &[f64],
    maximum_iterations: usize,
    target_l1: f64,
    mut apply: impl FnMut(&[f64], &mut [f64]) -> Result<(), GasEnergyError>,
) -> Result<LinearSolution, GasEnergyError> {
    let half = rhs.len();
    let mut solution = vec![0.0; half];
    let mut residual = rhs.to_vec();
    let mut direction = residual.clone();
    let mut product = vec![0.0; half];
    let mut rho = dot(&residual, &residual);
    for iteration in 0..=maximum_iterations {
        let estimated = compensated(residual.iter().map(|x| x.abs()));
        if estimated <= target_l1 || iteration == maximum_iterations {
            // Never accept a recursively updated CG residual as final proof.
            apply(&solution, &mut product)?;
            for i in 0..half {
                residual[i] = rhs[i] - product[i];
            }
            let residual_l1 = compensated(residual.iter().map(|x| x.abs()));
            if residual_l1 <= target_l1 {
                return Ok(LinearSolution {
                    values: solution,
                    residual_l1,
                    iterations: iteration,
                });
            }
            if iteration == maximum_iterations {
                return Err(invalid("Clausing linear solve budget exhausted"));
            }
            direction.copy_from_slice(&residual);
            rho = dot(&residual, &residual);
        }
        apply(&direction, &mut product)?;
        let denominator = dot(&direction, &product);
        let alpha = rho / denominator;
        if !positive(rho) || !positive(denominator) || !positive(alpha) {
            return Err(invalid("Clausing conjugate-gradient arithmetic failed"));
        }
        for i in 0..half {
            solution[i] += alpha * direction[i];
            residual[i] -= alpha * product[i];
        }
        let next_rho = dot(&residual, &residual);
        let beta = next_rho / rho;
        if !beta.is_finite() || beta < 0.0 {
            return Err(invalid(
                "Clausing conjugate-gradient residual is not representable",
            ));
        }
        for i in 0..half {
            direction[i] = residual[i] + beta * direction[i];
        }
        rho = next_rho;
    }
    Err(invalid("Clausing linear solve budget exhausted"))
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    compensated(left.iter().zip(right).map(|(a, b)| a * b))
}

fn compensated(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
mod tests;
