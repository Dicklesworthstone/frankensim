//! Joint exit angle and collision history for Bernoulli wall-speed resets.
//!
//! The cosine transport law and its angular kernels belong to the parent.
//! With q=1-alpha, A=I-W, B=I-qW, y=B^-1 eL and u=A^-1 mL, reciprocity gives
//! retained wall moments C*q*y'm_exit and reset moments C*alpha*y'u_exit,
//! C=4*delta times the common incident normalization. The identity
//! A^-1-q*B^-1=alpha*A^-1*B^-1 includes every encounter, including the first.
//! The reset event is independent of speed and direction, and wall temperature
//! is constant along the bore: the last reset has that same speed distribution.
//!
//! Positive formulas retain both tiny-alpha reset weights and the tiny retained
//! complement near alpha=1. After checking the resolvent partition, only the
//! larger component is formed as total minus smaller. Mean energy accommodation
//! by itself does not establish this speed-mixture law.

use super::super::accommodation::{
    self, AccommodationResolution, HeatSolution, MolecularAccommodation, heat_from_field,
    solve_discounted,
};
use super::{
    AngularIntegrals, ClausingControls, CylindricalMolecularGeometry, GasEnergyError, Grid,
    GridSolution, MolecularAxialMoments, Moments, compensated, invalid, poll, positive,
};

/// Geometric weights for original-speed and wall-reset molecules, per original
/// incident molecule. Retained transmission INCLUDES the direct flight.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularAxialAccommodation {
    accommodation: MolecularAccommodation,
    retained: [f64; 2],
    reset: [f64; 2],
    resolution: Option<AxialAccommodationResolution>,
}
impl MolecularAxialAccommodation {
    pub const fn coefficient(&self) -> f64 {
        self.accommodation.coefficient()
    }
    /// Heat probabilities from the same final mass/angular grid and coefficient.
    pub const fn accommodation(&self) -> &MolecularAccommodation {
        &self.accommodation
    }
    pub const fn retained_transmitted_axial_cosine(&self) -> f64 {
        self.retained[0]
    }
    pub const fn reset_transmitted_axial_cosine(&self) -> f64 {
        self.reset[0]
    }
    pub const fn retained_reflected_axial_cosine(&self) -> f64 {
        self.retained[1]
    }
    pub const fn reset_reflected_axial_cosine(&self) -> f64 {
        self.reset[1]
    }
    /// Endpoints and zero-length apertures have an exact algebraic split,
    /// relative to the common geometry and total-angular resolution.
    pub const fn resolution(&self) -> Option<&AxialAccommodationResolution> {
        self.resolution.as_ref()
    }
}

/// Every two-element array is ordered [transmitted, reflected]. Estimates
/// exclude roundoff and are not certified continuum enclosures. The retained
/// transmission error includes the common direct-flight quadrature error;
/// its grid change and linear error concern only wall-encountered molecules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialAccommodationResolution {
    pub controls: ClausingControls,
    pub rings: usize,
    pub angular_panels: usize,
    pub refinements: usize,
    pub batch_size: usize,
    pub thermal_coefficients: usize,
    /// All shared mass/angular and distinct discounted CG steps on all grids.
    pub iterations: usize,
    /// The one coefficient-independent angular parity pair on each grid.
    pub angular_escape_iterations: usize,
    pub symmetric_residual_l1: f64,
    pub antisymmetric_residual_l1: f64,
    pub retained_grid_change: [f64; 2],
    pub reset_grid_change: [f64; 2],
    pub retained_quadrature_estimate: [f64; 2],
    pub reset_quadrature_estimate: [f64; 2],
    pub retained_linear_solve_estimate: [f64; 2],
    pub reset_linear_solve_estimate: [f64; 2],
    pub previous_retained_linear_solve_estimate: [f64; 2],
    pub previous_reset_linear_solve_estimate: [f64; 2],
    pub retained_error_estimate: [f64; 2],
    pub reset_error_estimate: [f64; 2],
}

pub(super) fn exact_response(
    geometry: &CylindricalMolecularGeometry,
    moments: &MolecularAxialMoments,
    coefficient: f64,
) -> MolecularAxialAccommodation {
    let full = coefficient == 1.0;
    MolecularAxialAccommodation {
        accommodation: accommodation::exact_response(geometry, coefficient),
        retained: if full {
            [moments.direct, 0.0]
        } else {
            [moments.transmitted_axial_cosine(), moments.reflected]
        },
        reset: if full {
            [moments.wall_transmitted, moments.reflected]
        } else {
            [0.0; 2]
        },
        resolution: None,
    }
}

// [retained-transmitted-wall, reset-transmitted, retained-returned, reset-returned].
struct JointSolution {
    heat: HeatSolution,
    values: [f64; 4],
    linear: [f64; 4],
    quadrature: [f64; 4],
}
impl JointSolution {
    fn changes(&self, old: &Self) -> [f64; 4] {
        core::array::from_fn(|i| (self.values[i] - old.values[i]).abs())
    }
    fn errors(&self, old: &Self) -> [f64; 4] {
        let changes = self.changes(old);
        core::array::from_fn(|i| {
            changes[i] + self.linear[i] + old.linear[i] + self.quadrature[i] + old.quadrature[i]
        })
    }
}

#[derive(Default)]
pub(super) struct GridResponses {
    rows: Vec<JointSolution>,
    pub(super) iterations: usize,
    pub(super) angular_escape_iterations: usize,
}
impl GridResponses {
    pub(super) fn converged(
        &self,
        old: &Self,
        moments: &Moments,
        tolerance: impl Fn(f64) -> f64,
    ) -> bool {
        self.rows.iter().zip(&old.rows).all(|(current, old)| {
            let heat_error = current.heat.linear_estimate + old.heat.linear_estimate;
            let errors = current.errors(old);
            (current.heat.transmitted - old.heat.transmitted).abs() + heat_error
                <= tolerance(current.heat.transmitted)
                && (current.heat.reflected - old.heat.reflected).abs() + heat_error
                    <= tolerance(current.heat.reflected)
                && current
                    .values
                    .iter()
                    .zip(errors)
                    .all(|(&value, error)| error <= tolerance(value))
                && errors[0] + moments.quadrature[0]
                    <= tolerance(current.values[0] + moments.values[0])
        })
    }
}

pub(super) fn solve_grid(
    grid: &Grid,
    mass: &GridSolution,
    angular: &AngularIntegrals,
    moments: &Moments,
    kernel_relative_error: f64,
    coefficients: &[f64],
    controls: ClausingControls,
    target: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GridResponses, GasEnergyError> {
    let mut result = GridResponses::default();
    if coefficients.is_empty() {
        return Ok(result);
    }
    let rings = grid.end.len();
    // This one angular escape field serves every coefficient. Mirroring gives
    // the opposite exit, so neither a second mass solve nor species-sized
    // collections of ring fields are necessary.
    let escape = solve_discounted(
        grid,
        &angular.end,
        0.0,
        controls.maximum_iterations,
        target,
        cancelled,
    )?;
    result.angular_escape_iterations = escape.iterations;
    result.iterations = escape.iterations;
    let escape_residual = 2.0 * compensated(escape.residuals.into_iter());
    let weight = (mass.wall_transmitted + mass.reflected) / compensated(grid.end.iter().copied());
    let escape_floor = (0..rings / 2)
        .map(|i| grid.end[i] + grid.end[rings - 1 - i])
        .fold(1.0, f64::min);
    for &coefficient in coefficients {
        poll(cancelled)?;
        let field = solve_discounted(
            grid,
            &grid.end,
            coefficient,
            controls.maximum_iterations,
            target,
            cancelled,
        )?;
        result.iterations += field.iterations;
        let heat = heat_from_field(grid, mass, coefficient, &field)?;
        let retention = 1.0 - coefficient;
        let factor = coefficient / (coefficient + retention * escape_floor);
        let field_residual = 2.0 * compensated(field.residuals.into_iter());
        let retained_error = weight * retention * field_residual;
        // The discounted inverse has row sum <= 1/(alpha+q*escape_floor).
        // Also alpha*A^-1*B^-1*eL <= factor. These two positive resolvent
        // bounds propagate both residuals without erasing tiny-alpha accuracy.
        let reset_error = weight * factor * (field_residual + escape_residual);
        let mut values = [0.0; 4];
        let mut linear = [0.0; 4];
        let mut quadrature = [0.0; 4];
        for side in 0..2 {
            let retained = weight
                * retention
                * compensated((0..rings / 2).map(|i| {
                    let pair = field.pairs[i];
                    if side == 0 {
                        pair[0] * angular.end[rings - 1 - i] + pair[1] * angular.end[i]
                    } else {
                        pair[0] * angular.end[i] + pair[1] * angular.end[rings - 1 - i]
                    }
                }));
            let reset = weight
                * coefficient
                * compensated((0..rings / 2).map(|i| {
                    let pair = field.pairs[i];
                    let escape = escape.pairs[i];
                    if side == 0 {
                        pair[0] * escape[1] + pair[1] * escape[0]
                    } else {
                        pair[0] * escape[0] + pair[1] * escape[1]
                    }
                }));
            let total = moments.values[side + 1];
            let discrepancy = compensated([retained, reset, -total].into_iter()).abs();
            if ![retained, reset].into_iter().all(positive)
                || discrepancy
                    > retained_error
                        + reset_error
                        + mass.linear_estimate
                        + 64.0 * rings as f64 * f64::EPSILON * total
            {
                return Err(invalid(
                    "joint angular/collision moments failed their resolvent partition",
                ));
            }
            let (small, error, index) = if retained < reset {
                (retained, retained_error, 0)
            } else {
                (reset, reset_error, 1)
            };
            let large = total - small;
            if !positive(large) || small > total {
                return Err(invalid(
                    "joint angular/collision moments left their physical bounds",
                ));
            }
            let small_index = 2 * side + index;
            let large_index = 2 * side + 1 - index;
            values[small_index] = small;
            values[large_index] = large;
            linear[small_index] = error;
            linear[large_index] = mass.linear_estimate + error;
            // Each exact split is a positive linear functional of the exit
            // kernel, so a componentwise relative kernel difference bounds
            // its quadrature change, including arbitrarily small weights.
            quadrature[small_index] = kernel_relative_error * (small + error);
            quadrature[large_index] = moments.quadrature[side + 1] + quadrature[small_index];
        }
        result.rows.push(JointSolution {
            heat,
            values,
            linear,
            quadrature,
        });
    }
    Ok(result)
}

pub(super) fn responses(
    geometry: &CylindricalMolecularGeometry,
    moments: &MolecularAxialMoments,
    coefficients: &[f64],
    distinct: &[f64],
    grids: (&GridResponses, &GridResponses),
    iterations: usize,
    angular_escape_iterations: usize,
) -> Result<Vec<MolecularAxialAccommodation>, GasEnergyError> {
    let report = moments.resolution();
    coefficients
        .iter()
        .map(|&coefficient| {
            if coefficient == 0.0 || coefficient == 1.0 {
                return Ok(exact_response(geometry, moments, coefficient));
            }
            let index = distinct
                .binary_search_by(|value| value.total_cmp(&coefficient))
                .map_err(|_| invalid("joint angular accommodation response was not resolved"))?;
            let (current, old) = (&grids.0.rows[index], &grids.1.rows[index]);
            let changes = current.changes(old);
            let errors = current.errors(old);
            let heat = &current.heat;
            let old_heat = &old.heat;
            let heat_change = [
                (heat.transmitted - old_heat.transmitted).abs(),
                (heat.reflected - old_heat.reflected).abs(),
            ];
            let accommodation = MolecularAccommodation {
                coefficient,
                transmitted_heat_probability: heat.transmitted,
                reflected_heat_probability: heat.reflected,
                resolution: Some(AccommodationResolution {
                    controls: report.controls,
                    rings: report.rings,
                    refinements: report.refinements,
                    batch_size: coefficients.len(),
                    thermal_coefficients: distinct.len(),
                    iterations: iterations - angular_escape_iterations,
                    symmetric_residual_l1: heat.residuals[0],
                    antisymmetric_residual_l1: heat.residuals[1],
                    transmitted_grid_change: heat_change[0],
                    reflected_grid_change: heat_change[1],
                    linear_solve_estimate: heat.linear_estimate,
                    previous_linear_solve_estimate: old_heat.linear_estimate,
                    transmitted_error_estimate: heat_change[0]
                        + heat.linear_estimate
                        + old_heat.linear_estimate,
                    reflected_error_estimate: heat_change[1]
                        + heat.linear_estimate
                        + old_heat.linear_estimate,
                }),
            };
            Ok(MolecularAxialAccommodation {
                accommodation,
                retained: [moments.direct + current.values[0], current.values[2]],
                reset: [current.values[1], current.values[3]],
                resolution: Some(AxialAccommodationResolution {
                    controls: report.controls,
                    rings: report.rings,
                    angular_panels: report.angular_panels,
                    refinements: report.refinements,
                    batch_size: coefficients.len(),
                    thermal_coefficients: distinct.len(),
                    iterations,
                    angular_escape_iterations,
                    symmetric_residual_l1: heat.residuals[0],
                    antisymmetric_residual_l1: heat.residuals[1],
                    retained_grid_change: [changes[0], changes[2]],
                    reset_grid_change: [changes[1], changes[3]],
                    retained_quadrature_estimate: [
                        current.quadrature[0] + report.direct_quadrature_estimate,
                        current.quadrature[2],
                    ],
                    reset_quadrature_estimate: [current.quadrature[1], current.quadrature[3]],
                    retained_linear_solve_estimate: [current.linear[0], current.linear[2]],
                    reset_linear_solve_estimate: [current.linear[1], current.linear[3]],
                    previous_retained_linear_solve_estimate: [old.linear[0], old.linear[2]],
                    previous_reset_linear_solve_estimate: [old.linear[1], old.linear[3]],
                    retained_error_estimate: [errors[0] + report.direct_error_estimate, errors[2]],
                    reset_error_estimate: [errors[1], errors[3]],
                }),
            })
        })
        .collect()
}
