//! Repeated, diffuse wall encounters with partial total-energy accommodation.
//!
//! A collision replaces the mean energy excess over the wall by (1-alpha)
//! times that excess. Thus an N-hit trajectory thermalizes 1-(1-alpha)^N.
//! The single-collision definition is Chiggiato, arXiv:1501.07162, Sec. 2.2.2.
//! Here it is composed with the Clausing Markov kernel, not with a specular
//! fraction: angles remain diffuse and independent of energy at every hit.
//!
//! For A=I-W, B=I-(1-alpha)W, w=A^-1 eR and y=B^-1 eL, the incident-normalized
//! transmitted heat weight is 4*delta*alpha*y'w. Replacing w by 1-w gives the
//! returned heat weight. These resolvents account for arbitrarily many hits.

use super::{
    ClausingControls, ClausingResolution, CylindricalMolecularGeometry, GasEnergyError, Grid,
    GridSolution, compensated, conjugate_gradient, invalid, poll, positive,
};

/// Geometry-bound first-energy-moment response. A material consumer must qualify
/// the coefficient for its original species and temperature/pressure support.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularAccommodation {
    coefficient: f64,
    transmitted_heat_probability: f64,
    reflected_heat_probability: f64,
    resolution: Option<AccommodationResolution>,
}

impl MolecularAccommodation {
    pub const fn coefficient(&self) -> f64 {
        self.coefficient
    }
    /// Mean fraction thermalized, weighted by eventual transmission, per entry.
    pub const fn transmitted_heat_probability(&self) -> f64 {
        self.transmitted_heat_probability
    }
    /// Mean fraction thermalized, weighted by return to the source, per entry.
    pub const fn reflected_heat_probability(&self) -> f64 {
        self.reflected_heat_probability
    }
    /// None denotes an exact algebraic endpoint law relative to the supplied
    /// geometry, whose transmission may still carry a numerical estimate.
    pub const fn resolution(&self) -> Option<&AccommodationResolution> {
        self.resolution.as_ref()
    }
}

/// Numerical evidence for both energy weights; estimates exclude roundoff and
/// are not a certified continuum enclosure. The returned geometry also retains
/// its independent mass-transmission report on these same rings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AccommodationResolution {
    pub controls: ClausingControls,
    pub rings: usize,
    pub refinements: usize,
    /// All mass and thermal CG steps over every attempted grid.
    pub iterations: usize,
    pub symmetric_residual_l1: f64,
    pub antisymmetric_residual_l1: f64,
    pub transmitted_grid_change: f64,
    pub reflected_grid_change: f64,
    pub linear_solve_estimate: f64,
    pub previous_linear_solve_estimate: f64,
    pub transmitted_error_estimate: f64,
    pub reflected_error_estimate: f64,
}

pub(in super::super) fn resolve(
    mut geometry: CylindricalMolecularGeometry,
    coefficient: f64,
    controls: ClausingControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(CylindricalMolecularGeometry, MolecularAccommodation), GasEnergyError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !coefficient.is_finite() || !(0.0..=1.0).contains(&coefficient) {
        return Err(invalid(
            "thermal accommodation coefficient must be finite and in [0,1]",
        ));
    }
    let mut response = MolecularAccommodation {
        coefficient,
        transmitted_heat_probability: 0.0,
        reflected_heat_probability: 0.0,
        resolution: None,
    };
    if coefficient == 0.0 || coefficient == 1.0 || geometry.length_m == 0.0 {
        if coefficient == 1.0 {
            response.transmitted_heat_probability = geometry.wall_transmitted;
            response.reflected_heat_probability = geometry.reflected;
        }
        return Ok((geometry, response));
    }
    let target =
        0.01 * (controls.absolute_tolerance + controls.relative_tolerance * geometry.transmission);
    if !positive(target) {
        return Err(invalid(
            "accommodation solve tolerance is not representable",
        ));
    }
    let h = geometry.length_m / geometry.diameter_m;
    let mut rings = controls.initial_rings;
    let mut refinements = 0;
    let mut iterations = 0;
    let mut mass_iterations = 0;
    let mut previous: Option<(GridSolution, HeatSolution)> = None;
    loop {
        let grid = Grid::new(h, rings, &mut cancelled)?;
        let mass = grid.solve(
            controls.maximum_iterations,
            target,
            geometry.wall_hit,
            &mut cancelled,
        )?;
        let heat = solve_heat(
            &grid,
            &mass,
            coefficient,
            controls.maximum_iterations,
            target,
            &mut cancelled,
        )?;
        mass_iterations += mass.iterations;
        iterations += mass.iterations + heat.iterations;
        if let Some((old_mass, old_heat)) = previous {
            let mass_change = (mass.wall_transmitted - old_mass.wall_transmitted).abs();
            let mass_error = mass_change + mass.linear_estimate + old_mass.linear_estimate;
            let transmitted_change = (heat.transmitted - old_heat.transmitted).abs();
            let reflected_change = (heat.reflected - old_heat.reflected).abs();
            let transmitted_error =
                transmitted_change + heat.linear_estimate + old_heat.linear_estimate;
            let reflected_error =
                reflected_change + heat.linear_estimate + old_heat.linear_estimate;
            let transmission = geometry.direct + mass.wall_transmitted;
            let tolerance =
                |value: f64| controls.absolute_tolerance + controls.relative_tolerance * value;
            if refinements >= 2
                && mass_error <= tolerance(transmission)
                && transmitted_error <= tolerance(heat.transmitted)
                && reflected_error <= tolerance(heat.reflected)
            {
                poll(&mut cancelled)?;
                geometry.transmission = transmission;
                geometry.wall_transmitted = mass.wall_transmitted;
                geometry.reflected = mass.reflected;
                geometry.resolution = Some(ClausingResolution {
                    controls,
                    rings,
                    refinements,
                    iterations: mass_iterations,
                    residual_l1: mass.residual_l1,
                    discretization_estimate: mass_change,
                    linear_solve_estimate: mass.linear_estimate,
                    previous_linear_solve_estimate: old_mass.linear_estimate,
                    estimated_absolute_error: mass_error,
                });
                response.transmitted_heat_probability = heat.transmitted;
                response.reflected_heat_probability = heat.reflected;
                response.resolution = Some(AccommodationResolution {
                    controls,
                    rings,
                    refinements,
                    iterations,
                    symmetric_residual_l1: heat.residuals[0],
                    antisymmetric_residual_l1: heat.residuals[1],
                    transmitted_grid_change: transmitted_change,
                    reflected_grid_change: reflected_change,
                    linear_solve_estimate: heat.linear_estimate,
                    previous_linear_solve_estimate: old_heat.linear_estimate,
                    transmitted_error_estimate: transmitted_error,
                    reflected_error_estimate: reflected_error,
                });
                return Ok((geometry, response));
            }
        }
        previous = Some((mass, heat));
        if rings > controls.maximum_rings / 2 {
            return Err(invalid("accommodation spatial refinement budget exhausted"));
        }
        rings *= 2;
        refinements += 1;
    }
}

struct HeatSolution {
    transmitted: f64,
    reflected: f64,
    linear_estimate: f64,
    residuals: [f64; 2],
    iterations: usize,
}

fn solve_heat(
    grid: &Grid,
    mass: &GridSolution,
    coefficient: f64,
    maximum_iterations: usize,
    target: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HeatSolution, GasEnergyError> {
    let rings = grid.wall.len();
    let half = rings / 2;
    let mut fields = Vec::with_capacity(2);
    for symmetric in [true, false] {
        let rhs: Vec<_> = (0..half)
            .map(|i| {
                0.5 * (grid.end[i]
                    + if symmetric {
                        grid.end[rings - 1 - i]
                    } else {
                        -grid.end[rings - 1 - i]
                    })
            })
            .collect();
        // Resolve normalized fields even when alpha is tiny; otherwise an
        // absolute heat tolerance could erase a physically nonzero response.
        let target_l1 = (target / (8.0 * grid.delta))
            .min(1e-11 * compensated(rhs.iter().map(|value| value.abs())));
        fields.push(conjugate_gradient(
            &rhs,
            maximum_iterations,
            target_l1,
            |input, output| {
                apply_discounted(grid, input, output, coefficient, symmetric, cancelled)
            },
        )?);
    }
    let even = &fields[0];
    let odd = &fields[1];
    let pairs: Vec<_> = even
        .values
        .iter()
        .zip(&odd.values)
        .map(|(a, b)| [a + b, a - b])
        .collect();
    if pairs
        .iter()
        .flatten()
        .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
    {
        return Err(invalid(
            "discounted Clausing probabilities left their physical bounds",
        ));
    }
    let norm = (mass.wall_transmitted + mass.reflected)
        / (4.0 * grid.delta * compensated(grid.end.iter().copied()));
    let retention = 1.0 - coefficient;
    let residual_l1 = 2.0 * (even.residual_l1 + odd.residual_l1);
    let (transmitted, reflected, linear_estimate) = if coefficient < 0.5 {
        let weight = 4.0 * grid.delta * coefficient * norm;
        let transmitted = weight
            * compensated((0..half).map(|i| {
                pairs[i][0] * (0.5 + mass.half_escape[i])
                    + pairs[i][1] * (0.5 - mass.half_escape[i])
            }));
        let reflected = weight
            * compensated((0..half).map(|i| {
                pairs[i][0] * (0.5 - mass.half_escape[i])
                    + pairs[i][1] * (0.5 + mass.half_escape[i])
            }));
        let escape_floor = (0..half)
            .map(|i| grid.end[i] + grid.end[rings - 1 - i])
            .fold(1.0, f64::min);
        let factor = coefficient / (coefficient + retention * escape_floor);
        // B has row deficit >= alpha+(1-alpha)*escape_floor. The resolvent
        // identity alpha*A^-1*B^-1=A^-1-(1-alpha)*B^-1 propagates mass residual.
        let estimate = factor * (mass.linear_estimate + 4.0 * grid.delta * residual_l1);
        (transmitted, reflected, estimate)
    } else {
        // Retain the tiny unaccommodated complement near alpha=1 instead of
        // subtracting two almost equal heat weights after the solve.
        let weight = 4.0 * grid.delta * retention * norm;
        let retained_transmitted = weight
            * compensated(
                (0..half)
                    .map(|i| pairs[i][0] * grid.end[rings - 1 - i] + pairs[i][1] * grid.end[i]),
            );
        let retained_reflected = weight
            * compensated(
                (0..half)
                    .map(|i| pairs[i][0] * grid.end[i] + pairs[i][1] * grid.end[rings - 1 - i]),
            );
        (
            mass.wall_transmitted - retained_transmitted,
            mass.reflected - retained_reflected,
            mass.linear_estimate + 4.0 * grid.delta * retention * residual_l1,
        )
    };
    if ![
        (transmitted, mass.wall_transmitted),
        (reflected, mass.reflected),
    ]
    .iter()
    .all(|(value, bound)| positive(*value) && value <= bound)
        || !linear_estimate.is_finite()
    {
        return Err(invalid(
            "accommodation heat probabilities are not representable",
        ));
    }
    Ok(HeatSolution {
        transmitted,
        reflected,
        linear_estimate,
        residuals: [even.residual_l1, odd.residual_l1],
        iterations: even.iterations + odd.iterations,
    })
}

fn apply_discounted(
    grid: &Grid,
    input: &[f64],
    output: &mut [f64],
    coefficient: f64,
    symmetric: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasEnergyError> {
    let rings = grid.wall.len();
    for i in 0..input.len() {
        if i % 16 == 0 {
            poll(cancelled)?;
        }
        let ordinary = compensated(input.iter().enumerate().map(|(j, value)| {
            let entry = if i == j {
                grid.diagonal
            } else {
                -grid.wall[i.abs_diff(j)]
            };
            let mirror = grid.wall[rings - 1 - i - j];
            (entry + if symmetric { -mirror } else { mirror }) * value
        }));
        output[i] = coefficient * input[i] + (1.0 - coefficient) * ordinary;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
