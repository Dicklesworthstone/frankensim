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
//! Species may supply different coefficients. They share the geometric kernel
//! and mass solution, with joint refinement and one active thermal solve at a time.

use super::super::super::super::MAX_EFFUSIVE_SPECIES;
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
    /// Authored coefficient count, including duplicates and exact endpoints.
    pub batch_size: usize,
    /// Distinct coefficients strictly between zero and one solved on each grid.
    pub thermal_coefficients: usize,
    /// All shared mass and distinct thermal CG steps over every attempted grid.
    /// This is the same batch total in each resolved coefficient's report.
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
    geometry: CylindricalMolecularGeometry,
    coefficient: f64,
    controls: ClausingControls,
    cancelled: impl FnMut() -> bool,
) -> Result<(CylindricalMolecularGeometry, MolecularAccommodation), GasEnergyError> {
    let (geometry, mut responses) = resolve_batch(geometry, &[coefficient], controls, cancelled)?;
    let response = responses
        .pop()
        .ok_or_else(|| invalid("accommodation scalar response was not resolved"))?;
    Ok((geometry, response))
}

pub(in super::super) fn resolve_batch(
    mut geometry: CylindricalMolecularGeometry,
    coefficients: &[f64],
    controls: ClausingControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(CylindricalMolecularGeometry, Vec<MolecularAccommodation>), GasEnergyError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if coefficients.is_empty() || coefficients.len() > MAX_EFFUSIVE_SPECIES {
        return Err(invalid(
            "accommodation batch exceeds effusive species bounds",
        ));
    }
    if coefficients
        .iter()
        .any(|coefficient| !coefficient.is_finite() || !(0.0..=1.0).contains(coefficient))
    {
        return Err(invalid(
            "thermal accommodation coefficient must be finite and in [0,1]",
        ));
    }
    // Canonical solve order makes work and results independent of species order.
    let mut distinct: Vec<_> = coefficients
        .iter()
        .copied()
        .filter(|coefficient| *coefficient > 0.0 && *coefficient < 1.0)
        .collect();
    distinct.sort_by(f64::total_cmp);
    distinct.dedup();
    if distinct.is_empty() || geometry.length_m == 0.0 {
        let responses = coefficients
            .iter()
            .map(|&coefficient| exact_response(&geometry, coefficient))
            .collect();
        return Ok((geometry, responses));
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
    let mut previous: Option<(GridSolution, Vec<HeatSolution>)> = None;
    loop {
        let grid = Grid::new(h, rings, &mut cancelled)?;
        let mass = grid.solve(
            controls.maximum_iterations,
            target,
            geometry.wall_hit,
            &mut cancelled,
        )?;
        mass_iterations += mass.iterations;
        iterations += mass.iterations;
        let mut heat = Vec::with_capacity(distinct.len());
        for &coefficient in &distinct {
            poll(&mut cancelled)?;
            let response = solve_heat(
                &grid,
                &mass,
                coefficient,
                controls.maximum_iterations,
                target,
                &mut cancelled,
            )?;
            iterations += response.iterations;
            heat.push(response);
        }
        if let Some((old_mass, old_heat)) = &previous {
            let mass_change = (mass.wall_transmitted - old_mass.wall_transmitted).abs();
            let mass_error = mass_change + mass.linear_estimate + old_mass.linear_estimate;
            let transmission = geometry.direct + mass.wall_transmitted;
            let tolerance =
                |value: f64| controls.absolute_tolerance + controls.relative_tolerance * value;
            if refinements >= 2
                && mass_error <= tolerance(transmission)
                && heat.iter().zip(old_heat).all(|(current, old)| {
                    (current.transmitted - old.transmitted).abs()
                        + current.linear_estimate
                        + old.linear_estimate
                        <= tolerance(current.transmitted)
                        && (current.reflected - old.reflected).abs()
                            + current.linear_estimate
                            + old.linear_estimate
                            <= tolerance(current.reflected)
                })
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
                let resolved: Vec<_> = distinct
                    .iter()
                    .zip(heat.iter().zip(old_heat))
                    .map(|(&coefficient, (current, old))| {
                        let transmitted_change = (current.transmitted - old.transmitted).abs();
                        let reflected_change = (current.reflected - old.reflected).abs();
                        MolecularAccommodation {
                            coefficient,
                            transmitted_heat_probability: current.transmitted,
                            reflected_heat_probability: current.reflected,
                            resolution: Some(AccommodationResolution {
                                controls,
                                rings,
                                refinements,
                                batch_size: coefficients.len(),
                                thermal_coefficients: distinct.len(),
                                iterations,
                                symmetric_residual_l1: current.residuals[0],
                                antisymmetric_residual_l1: current.residuals[1],
                                transmitted_grid_change: transmitted_change,
                                reflected_grid_change: reflected_change,
                                linear_solve_estimate: current.linear_estimate,
                                previous_linear_solve_estimate: old.linear_estimate,
                                transmitted_error_estimate: transmitted_change
                                    + current.linear_estimate
                                    + old.linear_estimate,
                                reflected_error_estimate: reflected_change
                                    + current.linear_estimate
                                    + old.linear_estimate,
                            }),
                        }
                    })
                    .collect();
                let responses = coefficients
                    .iter()
                    .map(|&coefficient| {
                        if coefficient == 0.0 || coefficient == 1.0 {
                            Ok(exact_response(&geometry, coefficient))
                        } else {
                            distinct
                                .binary_search_by(|value| value.total_cmp(&coefficient))
                                .ok()
                                .and_then(|index| resolved.get(index).cloned())
                                .ok_or_else(|| {
                                    invalid("accommodation batch response was not resolved")
                                })
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                return Ok((geometry, responses));
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

fn exact_response(
    geometry: &CylindricalMolecularGeometry,
    coefficient: f64,
) -> MolecularAccommodation {
    MolecularAccommodation {
        coefficient,
        transmitted_heat_probability: if coefficient == 1.0 {
            geometry.wall_transmitted
        } else {
            0.0
        },
        reflected_heat_probability: if coefficient == 1.0 {
            geometry.reflected
        } else {
            0.0
        },
        resolution: None,
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
