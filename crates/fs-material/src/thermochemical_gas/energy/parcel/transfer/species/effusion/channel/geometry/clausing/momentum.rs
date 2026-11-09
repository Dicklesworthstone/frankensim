//! Axial angular moments of the original circular diffuse Clausing transport.
//!
//! Cosine-law incidence/reflection and angular beaming are described by Radke
//! and Jacobson, NASA/TM-2016-219118, sections 2.2, 3 and 4.1:
//! <https://ntrs.nasa.gov/api/citations/20180007747/downloads/20180007747.pdf>.
//! The following integrals are derived from that same geometric law. They are
//! not mass-transmission multipliers or empirical force factors.
//!
//! Let a be the distance from a wall ring to an end, in diameters, and let psi
//! be the azimuth about the tube axis measured toward the inward wall normal.
//! A diffuse direction reaches the end when mu > a/sqrt(a*a+cos(psi)^2).
//! Integrating mu times the cosine-law solid-angle density gives
//! M(a) = 1/(3*pi) integral cos(psi)^4/(a*a+cos(psi)^2)^(3/2) dpsi,
//! over -pi/2..pi/2. Its ring average is analytic in a; only psi is quadrature.
//! The direct moment is 4 integral_h^infinity M(a) da, evaluated as a positive
//! integral to retain strongly collimated small direct paths. Existing Clausing
//! escape probabilities give the wall visitation density by reciprocity.
//!
//! No speed distribution, molecular energy, moving-wall work or partial-speed
//! accommodation is inferred. A material consumer binds the returned geometry
//! and moments to its original physically qualified channel and speed laws.

use super::{
    ClausingControls, ClausingResolution, CylindricalMolecularGeometry, GasEnergyError, Grid,
    GridSolution, compensated, det, invalid, poll, positive,
};

/// Incident-normalized angular expectations. Returned and transmitted values
/// use each exit's OUTWARD axial direction; both are nonnegative.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularAxialMoments {
    direct: f64,
    wall_transmitted: f64,
    reflected: f64,
    resolution: AxialMomentumResolution,
}
impl MolecularAxialMoments {
    pub const fn direct_axial_cosine(&self) -> f64 {
        self.direct
    }
    pub const fn wall_transmitted_axial_cosine(&self) -> f64 {
        self.wall_transmitted
    }
    pub const fn reflected_axial_cosine(&self) -> f64 {
        self.reflected
    }
    pub fn transmitted_axial_cosine(&self) -> f64 {
        self.direct + self.wall_transmitted
    }
    pub const fn resolution(&self) -> &AxialMomentumResolution {
        &self.resolution
    }
}

/// Joint spatial/angular numerical evidence, excluding floating-point rounding.
/// Grid and quadrature differences are estimates, not certified enclosures.
/// Angular integration uses eight positive Gauss points per panel, with at most
/// `controls.maximum_rings/2` panels. Work and scratch are explicitly bounded.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialMomentumResolution {
    pub controls: ClausingControls,
    pub rings: usize,
    pub angular_panels: usize,
    pub refinements: usize,
    /// Shared mass CG iterations over every attempted spatial grid.
    pub iterations: usize,
    /// Angular-kernel evaluations over all coarse/fine quadrature attempts.
    pub quadrature_evaluations: usize,
    pub direct_quadrature_estimate: f64,
    pub wall_transmitted_quadrature_estimate: f64,
    pub reflected_quadrature_estimate: f64,
    pub wall_transmitted_grid_change: f64,
    pub reflected_grid_change: f64,
    pub linear_solve_estimate: f64,
    pub previous_linear_solve_estimate: f64,
    pub direct_error_estimate: f64,
    pub wall_transmitted_error_estimate: f64,
    pub reflected_error_estimate: f64,
    /// Direct + wall-transmitted + returned moments minus 2/3. Equal diffuse
    /// reservoirs must reproduce the complete equilibrium cosine distribution.
    pub reciprocity_residual: f64,
}

#[derive(Clone, Copy)]
struct Moments {
    values: [f64; 3],
    quadrature: [f64; 3],
}

pub(in super::super) fn resolve(
    mut geometry: CylindricalMolecularGeometry,
    controls: ClausingControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(CylindricalMolecularGeometry, MolecularAxialMoments), GasEnergyError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if geometry.length_m == 0.0 {
        let geometry = super::resolve(geometry, controls, &mut cancelled)?;
        return Ok((
            geometry,
            MolecularAxialMoments {
                direct: 2.0 / 3.0,
                wall_transmitted: 0.0,
                reflected: 0.0,
                resolution: AxialMomentumResolution {
                    controls,
                    rings: 0,
                    angular_panels: 0,
                    refinements: 0,
                    iterations: 0,
                    quadrature_evaluations: 0,
                    direct_quadrature_estimate: 0.0,
                    wall_transmitted_quadrature_estimate: 0.0,
                    reflected_quadrature_estimate: 0.0,
                    wall_transmitted_grid_change: 0.0,
                    reflected_grid_change: 0.0,
                    linear_solve_estimate: 0.0,
                    previous_linear_solve_estimate: 0.0,
                    direct_error_estimate: 0.0,
                    wall_transmitted_error_estimate: 0.0,
                    reflected_error_estimate: 0.0,
                    reciprocity_residual: 0.0,
                },
            },
        ));
    }
    let h = geometry.length_m / geometry.diameter_m;
    let tolerance = |value: f64| controls.absolute_tolerance + controls.relative_tolerance * value;
    let target = 0.01 * tolerance(geometry.transmission);
    if !positive(target) {
        return Err(invalid("axial moment solve tolerance is not representable"));
    }
    let mut rings = controls.initial_rings;
    let (mut refinements, mut iterations, mut evaluations) = (0, 0, 0);
    let mut panels = 4;
    let mut previous: Option<(GridSolution, Moments)> = None;
    loop {
        let grid = Grid::new(h, rings, &mut cancelled)?;
        let mass = grid.solve(
            controls.maximum_iterations,
            target,
            geometry.wall_hit,
            &mut cancelled,
        )?;
        iterations += mass.iterations;
        let mut coarse = angular_moments(&grid, &mass, h, panels / 2, &mut cancelled)?;
        evaluations += (rings + 1) * 8 * (panels / 2);
        let moments = loop {
            let fine = angular_moments(&grid, &mass, h, panels, &mut cancelled)?;
            evaluations += (rings + 1) * 8 * panels;
            let quadrature = core::array::from_fn(|i| (fine[i] - coarse[i]).abs());
            if fine
                .iter()
                .zip(quadrature)
                .all(|(&value, error)| error <= 0.1 * tolerance(value))
            {
                break Moments {
                    values: fine,
                    quadrature,
                };
            }
            if panels > controls.maximum_rings / 4 {
                return Err(invalid("axial moment angular quadrature budget exhausted"));
            }
            panels *= 2;
            coarse = fine;
        };
        if let Some((old_mass, old)) = &previous {
            let mass_change = (mass.wall_transmitted - old_mass.wall_transmitted).abs();
            let mass_error = mass_change + mass.linear_estimate + old_mass.linear_estimate;
            let transmission = geometry.direct + mass.wall_transmitted;
            let changes: [f64; 3] =
                core::array::from_fn(|i| (moments.values[i] - old.values[i]).abs());
            let errors = [
                moments.quadrature[0],
                changes[1]
                    + moments.quadrature[1]
                    + old.quadrature[1]
                    + mass.linear_estimate
                    + old_mass.linear_estimate,
                changes[2]
                    + moments.quadrature[2]
                    + old.quadrature[2]
                    + mass.linear_estimate
                    + old_mass.linear_estimate,
            ];
            if refinements >= 2
                && mass_error <= tolerance(transmission)
                && moments
                    .values
                    .iter()
                    .zip(errors)
                    .all(|(&value, error)| error <= tolerance(value))
            {
                poll(&mut cancelled)?;
                geometry.transmission = transmission;
                geometry.wall_transmitted = mass.wall_transmitted;
                geometry.reflected = mass.reflected;
                geometry.resolution = Some(ClausingResolution {
                    controls,
                    rings,
                    refinements,
                    iterations,
                    residual_l1: mass.residual_l1,
                    discretization_estimate: mass_change,
                    linear_solve_estimate: mass.linear_estimate,
                    previous_linear_solve_estimate: old_mass.linear_estimate,
                    estimated_absolute_error: mass_error,
                });
                let [direct, wall_transmitted, reflected] = moments.values;
                if ![
                    (direct, geometry.direct),
                    (wall_transmitted, geometry.wall_transmitted),
                    (reflected, geometry.reflected),
                ]
                .iter()
                .all(|(value, bound)| positive(*value) && value <= bound)
                {
                    return Err(invalid("axial channel moments left their physical bounds"));
                }
                let reciprocity_residual = compensated(moments.values.into_iter()) - 2.0 / 3.0;
                if reciprocity_residual.abs()
                    > errors.into_iter().sum::<f64>() + 64.0 * rings as f64 * f64::EPSILON
                {
                    return Err(invalid("axial channel moment reciprocity failed"));
                }
                let response = MolecularAxialMoments {
                    direct,
                    wall_transmitted,
                    reflected,
                    resolution: AxialMomentumResolution {
                        controls,
                        rings,
                        angular_panels: panels,
                        refinements,
                        iterations,
                        quadrature_evaluations: evaluations,
                        direct_quadrature_estimate: moments.quadrature[0],
                        wall_transmitted_quadrature_estimate: moments.quadrature[1],
                        reflected_quadrature_estimate: moments.quadrature[2],
                        wall_transmitted_grid_change: changes[1],
                        reflected_grid_change: changes[2],
                        linear_solve_estimate: mass.linear_estimate,
                        previous_linear_solve_estimate: old_mass.linear_estimate,
                        direct_error_estimate: errors[0],
                        wall_transmitted_error_estimate: errors[1],
                        reflected_error_estimate: errors[2],
                        reciprocity_residual,
                    },
                };
                return Ok((geometry, response));
            }
        }
        previous = Some((mass, moments));
        if rings > controls.maximum_rings / 2 {
            return Err(invalid("axial moment spatial refinement budget exhausted"));
        }
        rings *= 2;
        refinements += 1;
    }
}

// Eight-point Gauss-Legendre nodes/weights, paired about each panel center.
// Positive weights preserve nonnegative path moments; panel refinement supplies
// the explicitly retained angular quadrature estimate.
const GAUSS: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];

fn angular_moments(
    grid: &Grid,
    mass: &GridSolution,
    h: f64,
    panels: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], GasEnergyError> {
    let width = core::f64::consts::FRAC_PI_2 / panels as f64;
    let mut nodes = Vec::with_capacity(8 * panels);
    for panel in 0..panels {
        if panel % 16 == 0 {
            poll(cancelled)?;
        }
        let middle = (panel as f64 + 0.5) * width;
        for (abscissa, weight) in GAUSS {
            for sign in [-1.0, 1.0] {
                let cosine = det::cos(middle + sign * 0.5 * width * abscissa);
                nodes.push((cosine, weight * width / (3.0 * core::f64::consts::PI)));
            }
        }
    }
    let direct = 4.0
        * angular_integral(&nodes, cancelled, |c| {
            let c2 = c * c;
            let integrand = if h <= 1.0 {
                let root = det::sqrt(h * h + c2);
                c2 * (c / root) * (c / (root + h))
            } else {
                let inverse = 1.0 / h;
                let q = c * inverse;
                let root = det::sqrt(1.0 + q * q);
                c2 * q * q / (root * (root + 1.0))
            };
            integrand
        })?;
    let rings = grid.end.len();
    let mut end = Vec::with_capacity(rings);
    for i in 0..rings {
        poll(cancelled)?;
        let (a, b) = (i as f64 * grid.delta, (i + 1) as f64 * grid.delta);
        end.push(angular_integral(&nodes, cancelled, |c| {
            ring_integrand(a, b, c)
        })?);
    }
    let norm = (mass.wall_transmitted + mass.reflected)
        / (4.0 * grid.delta * compensated(grid.end.iter().copied()));
    let transmitted = 4.0
        * grid.delta
        * norm
        * compensated((0..rings / 2).map(|i| {
            end[i] * (0.5 + mass.half_escape[i])
                + end[rings - 1 - i] * (0.5 - mass.half_escape[i])
        }));
    let reflected = 4.0
        * grid.delta
        * norm
        * compensated((0..rings / 2).map(|i| {
            end[i] * (0.5 - mass.half_escape[i])
                + end[rings - 1 - i] * (0.5 + mass.half_escape[i])
        }));
    let values = [direct, transmitted, reflected];
    if !values.into_iter().all(positive) {
        return Err(invalid("axial angular integral is not representable"));
    }
    Ok(values)
}

fn angular_integral(
    nodes: &[(f64, f64)],
    cancelled: &mut impl FnMut() -> bool,
    value: impl Fn(f64) -> f64,
) -> Result<f64, GasEnergyError> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for (i, &(cosine, weight)) in nodes.iter().enumerate() {
        if i % 16 == 0 {
            poll(cancelled)?;
        }
        let term = weight * value(cosine);
        let next = sum + term;
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
    }
    Ok(sum + correction)
}

// cos^2(psi) * [b/sqrt(b^2+c^2)-a/sqrt(a^2+c^2)]/(b-a).
// Positive divided differences avoid subtracting nearby directions; scaled
// forms also retain very short rings and strongly collimated distant rings.
fn ring_integrand(a: f64, b: f64, c: f64) -> f64 {
    if b <= c {
        let (u, v) = (a / c, b / c);
        let (r0, r1) = (det::sqrt(1.0 + u * u), det::sqrt(1.0 + v * v));
        c * ((u + v) / (v * r0 + u * r1)) / (r0 * r1)
    } else {
        let q = c / b;
        let r1 = det::sqrt(1.0 + q * q);
        if a == 0.0 {
            return (c * c / b) / r1;
        }
        let u = a / b;
        let r0 = det::sqrt(u * u + q * q);
        (c * c * q * q / b) * ((u + 1.0) / (r0 + u * r1)) / (r0 * r1)
    }
}

#[cfg(test)]
mod tests;
