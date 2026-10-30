//! Source-owned Fourier conduction and full three-dimensional Newtonian stress.
//!
//! Every face evaluates the original transport law at its actual temperature,
//! pressure and frozen composition. Stress work and Fourier heat share one
//! conservative total-energy flux; kinetic energy lost to viscosity becomes
//! internal energy without an additional heat source. These operators do not
//! resolve a physical shock thickness without an adequate grid and source data.

use super::*;
use crate::compressible::molecular::parabolic;
pub use crate::compressible::molecular::{MolecularControl, MolecularError, MolecularStepControl};
use fs_gas::multicomponent::dissipation::TransportLaw;
use fs_qty::{Pressure, Temperature};

mod coupled;
pub use coupled::{MolecularEuler3d, MolecularStepReport3d};

/// One complete molecular interval, without advancing the enclosing flow clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport3d {
    /// Complete represented duration [s].
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps.
    pub accepted_steps: usize,
    /// Rejected internal attempts.
    pub rejected_steps: usize,
    /// Largest accepted-stage parabolic estimate.
    pub maximum_parabolic_number: f64,
    /// Original extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Complete inward molecular boundary transfer.
    pub boundary_change: Totals,
    /// `after - before - boundary_change`.
    pub conservation_residual: Totals,
}

struct MolecularWorkspace {
    current: Vec<State>,
    stage: Vec<State>,
    candidate: Vec<State>,
    primitive: Vec<Primitive>,
    // gradient[velocity component][derivative axis], prepared once per cell.
    gradient: Vec<[[f64; 3]; 3]>,
    flux: [Vec<SpatialFaceFlux>; 3],
    binary: Vec<f64>,
    boundary_step: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
}

/// Prepared three-dimensional molecular operator bound to the exact gas source.
///
/// Periodic faces share one flux. Reflecting faces are stationary adiabatic
/// slip walls: normal viscous stress loads the support, while tangential stress,
/// heat and work vanish at the wall. Transmissive ends have zero molecular
/// flux. These conditions do not represent no-slip or prescribed-temperature
/// surfaces. All state, gradient and face workspaces are preallocated.
pub struct MolecularTransport3d<L> {
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    law: L,
    boundary_impulses: Vec<BoundaryImpulse>,
    work: MolecularWorkspace,
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Bind complete original source models, fixed composition and geometry.
    pub fn new(
        flow: &Euler3d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        if law.source().models() != flow.gas.source().components() {
            return Err(MolecularError::Invalid(
                "transport and Euler caloric sources differ",
            ));
        }
        let n = flow.states.len();
        let species = law.source().models().len();
        let binary = species.checked_mul(species).ok_or(MolecularError::Invalid(
            "molecular species workspace overflow",
        ))?;
        let faces = flow.grid.boundary_len();
        let result = Self {
            gas: flow.gas.clone(),
            grid: flow.grid,
            boundaries: flow.boundaries,
            law,
            boundary_impulses: allocated(faces, Totals::default())?,
            work: MolecularWorkspace {
                current: allocated(n, State::default())?,
                stage: allocated(n, State::default())?,
                candidate: allocated(n, State::default())?,
                primitive: allocated(n, Primitive::default())?,
                gradient: allocated(n, [[0.0; 3]; 3])?,
                flux: [
                    allocated(flow.grid.face_len(0), SpatialFaceFlux::default())?,
                    allocated(flow.grid.face_len(1), SpatialFaceFlux::default())?,
                    allocated(flow.grid.face_len(2), SpatialFaceFlux::default())?,
                ],
                binary: allocated(binary, 0.0)?,
                boundary_step: allocated(faces, Totals::default())?,
                boundary_trial: allocated(faces, Totals::default())?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original immutable material/collision source and its reusable workspace.
    pub const fn law(&self) -> &L {
        &self.law
    }

    /// Outward molecular transfer during the last complete accepted interval.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.grid.boundary_range(side)]
    }

    /// Advance the whole molecular duration atomically, leaving time unchanged.
    ///
    /// The composing owner owns physical time. Existing Euler step metadata
    /// continues to describe its original advective interval; the returned
    /// molecular report and this operator's loads describe this update.
    pub fn advance(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport3d, MolecularError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid(
                "molecular operator and flow differ",
            ));
        }
        if !positive(duration.value()) {
            return Err(MolecularError::Invalid("molecular duration"));
        }
        let before = checked_totals(&flow.states, self.grid.cell_volume(), &mut cancelled)?;
        self.work.current.copy_from_slice(&flow.states);
        self.work.boundary_trial.fill(Totals::default());
        let mut elapsed = 0.0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut maximum_parabolic: f64 = 0.0;
        while elapsed < duration.value() {
            poll(&mut cancelled)?;
            if accepted + rejected >= controls.maximum_steps {
                return Err(MolecularError::WorkLimit {
                    elapsed_s: elapsed,
                    accepted,
                    rejected,
                });
            }
            let first_rate = spatial(
                &self.gas,
                &self.grid,
                self.boundaries,
                &self.work.current,
                &mut self.law,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.flux,
                &mut self.work.binary,
                &mut cancelled,
            )?;
            let remaining = duration.value() - elapsed;
            let mut dt = remaining.min(controls.maximum_step_s);
            if first_rate > 0.0 {
                dt = dt.min(controls.parabolic_cfl / first_rate);
            }
            for halving in 0..=controls.maximum_halvings {
                poll(&mut cancelled)?;
                if accepted + rejected >= controls.maximum_steps {
                    return Err(MolecularError::WorkLimit {
                        elapsed_s: elapsed,
                        accepted,
                        rejected,
                    });
                }
                if !positive(dt)
                    || elapsed + dt <= elapsed
                    || (dt < controls.minimum_step_s && dt != remaining)
                {
                    return Err(MolecularError::Invalid("unrepresentable molecular substep"));
                }
                let trial = (|| {
                    // First-trial faces are already prepared for the CFL. Only
                    // a retry needs to reconstruct them after stage two reused
                    // the buffers; this avoids one whole source/face pass.
                    let r1 = if halving == 0 {
                        first_rate
                    } else {
                        spatial(
                            &self.gas,
                            &self.grid,
                            self.boundaries,
                            &self.work.current,
                            &mut self.law,
                            &mut self.work.primitive,
                            &mut self.work.gradient,
                            &mut self.work.flux,
                            &mut self.work.binary,
                            &mut cancelled,
                        )?
                    };
                    parabolic(r1 * dt, controls.parabolic_cfl)?;
                    evolve::update(
                        &self.gas,
                        &self.grid,
                        &self.work.current,
                        &self.work.current,
                        &self.work.flux,
                        &mut self.work.stage,
                        dt,
                        false,
                        &mut cancelled,
                    )?;
                    evolve::capture_boundaries(
                        &self.grid,
                        &self.work.flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        false,
                        &mut cancelled,
                    )?;
                    let r2 = spatial(
                        &self.gas,
                        &self.grid,
                        self.boundaries,
                        &self.work.stage,
                        &mut self.law,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
                        &mut self.work.flux,
                        &mut self.work.binary,
                        &mut cancelled,
                    )?;
                    parabolic(r2 * dt, controls.parabolic_cfl)?;
                    evolve::update(
                        &self.gas,
                        &self.grid,
                        &self.work.current,
                        &self.work.stage,
                        &self.work.flux,
                        &mut self.work.candidate,
                        dt,
                        true,
                        &mut cancelled,
                    )?;
                    evolve::capture_boundaries(
                        &self.grid,
                        &self.work.flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        true,
                        &mut cancelled,
                    )?;
                    Ok::<f64, MolecularError>(r1.max(r2) * dt)
                })();
                match trial {
                    Ok(number) => {
                        maximum_parabolic = maximum_parabolic.max(number);
                        for (i, (total, step)) in self
                            .work
                            .boundary_trial
                            .iter_mut()
                            .zip(&self.work.boundary_step)
                            .enumerate()
                        {
                            if i % 64 == 0 {
                                poll(&mut cancelled)?;
                            }
                            *total = add_totals(*total, *step);
                        }
                        std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                        elapsed = if dt == remaining {
                            duration.value()
                        } else {
                            elapsed + dt
                        };
                        accepted += 1;
                        break;
                    }
                    Err(MolecularError::Cancelled) => return Err(MolecularError::Cancelled),
                    Err(error) => {
                        rejected += 1;
                        if halving == controls.maximum_halvings
                            || dt * 0.5 < controls.minimum_step_s
                        {
                            return Err(MolecularError::StepFailure {
                                elapsed_s: elapsed,
                                step_s: dt,
                                last: Box::new(error),
                            });
                        }
                        dt *= 0.5;
                    }
                }
            }
        }
        let after = checked_totals(&self.work.current, self.grid.cell_volume(), &mut cancelled)?;
        let outward = sum_arrays(
            self.work.boundary_trial.iter().map(|q| q.array()),
            &mut cancelled,
        )?;
        let boundary_change = Totals::from_array(outward.map(|v| -v));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(&mut cancelled)?;
        flow.states.copy_from_slice(&self.work.current);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        Ok(MolecularReport3d {
            elapsed: duration,
            accepted_steps: accepted,
            rejected_steps: rejected,
            maximum_parabolic_number: maximum_parabolic,
            before,
            after,
            boundary_change,
            conservation_residual,
        })
    }
}

fn add_totals(a: Totals, b: Totals) -> Totals {
    let a = a.array();
    let b = b.array();
    Totals::from_array(std::array::from_fn(|k| a[k] + b[k]))
}

fn balance(before: Totals, after: Totals, boundary: Totals) -> Result<Totals, MolecularError> {
    let a = after.array();
    let b = before.array();
    let c = boundary.array();
    let result = Totals::from_array(std::array::from_fn(|k| a[k] - b[k] - c[k]));
    if !a
        .into_iter()
        .chain(b)
        .chain(c)
        .chain(result.array())
        .all(f64::is_finite)
    {
        return Err(MolecularError::Invalid(
            "unrepresentable molecular extensive balance",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn spatial<L: TransportLaw>(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    states: &[State],
    law: &mut L,
    primitive: &mut [Primitive],
    gradient: &mut [[[f64; 3]; 3]],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    binary: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MolecularError> {
    let mut minimum_density = f64::INFINITY;
    let mut minimum_cv = f64::INFINITY;
    for (q, p) in states.iter().zip(primitive.iter_mut()) {
        *p = q.primitive(gas, &mut *cancelled)?;
        let props = gas
            .source()
            .evaluate_properties(Temperature::new(p.temperature), Pressure::new(p.pressure))
            .map_err(Error::from)?;
        let volumetric_cv = p.density * props.mass_specific_cv().value();
        if !positive(volumetric_cv) {
            return Err(MolecularError::Invalid(
                "molecular volumetric heat capacity",
            ));
        }
        minimum_density = minimum_density.min(p.density);
        minimum_cv = minimum_cv.min(volumetric_cv);
    }
    // Cell central gradients are reused by the six adjacent face traces. This
    // retains the same centered stencil without redoing every neighbor lookup.
    for (cell, derivatives) in gradient.iter_mut().enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        let c = grid.coordinates(cell);
        for axis in 0..3 {
            let lower = evolve::neighbour(primitive, grid, boundaries, c, axis, false);
            let upper = evolve::neighbour(primitive, grid, boundaries, c, axis, true);
            for component in 0..3 {
                derivatives[component][axis] =
                    0.5 * (upper[component + 1] - lower[component + 1]) / grid.spacing[axis];
            }
        }
    }
    let mut maximum_diffusivity: f64 = 0.0;
    for axis in 0..3 {
        let count = grid.shape[axis];
        let pair = boundaries.axis(axis);
        let tangential = tangents(axis);
        let flux = &mut fluxes[axis];
        let mut shape = grid.shape;
        shape[axis] += 1;
        for index in 0..flux.len() {
            poll(cancelled)?;
            let mut c = [
                index % shape[0],
                (index / shape[0]) % shape[1],
                index / (shape[0] * shape[1]),
            ];
            let face = c[axis];
            if face == count && pair == AxisBoundaries::Periodic {
                c[axis] = 0;
                flux[index] = flux[grid.cell_face(axis, c)];
                continue;
            }
            let exterior = face == 0 || face == count;
            let wall = if let AxisBoundaries::Ends { lower, upper } = pair {
                exterior && (if face == 0 { lower } else { upper }) == BoundaryCondition::Reflecting
            } else {
                false
            };
            if exterior && pair != AxisBoundaries::Periodic && !wall {
                flux[index] = SpatialFaceFlux::default();
                continue;
            }
            let mut left = c;
            let mut right = c;
            left[axis] = if face == 0 {
                if wall { 0 } else { count - 1 }
            } else {
                face - 1
            };
            right[axis] = if face == count { count - 1 } else { face };
            let left = grid.index(left);
            let right = grid.index(right);
            let pl = primitive[left];
            let pr = primitive[right];
            let temperature = 0.5 * pl.temperature + 0.5 * pr.temperature;
            let pressure = 0.5 * pl.pressure + 0.5 * pr.pressure;
            let coefficients = law.evaluate(
                Temperature::new(temperature),
                Pressure::new(pressure),
                gas.source().mole_composition().fractions(),
                binary,
                cancelled,
            )?;
            let mu = coefficients.dynamic_viscosity.value();
            let bulk = coefficients.bulk_viscosity.value();
            let conductivity = coefficients.thermal_conductivity.value();
            if ![mu, bulk, conductivity]
                .into_iter()
                .all(|v| v.is_finite() && v >= 0.0)
            {
                return Err(MolecularError::Invalid(
                    "nonnegative finite molecular coefficients required",
                ));
            }
            let longitudinal = (4.0 / 3.0) * mu + bulk;
            let cross = bulk - (2.0 / 3.0) * mu;
            // Cover both tangential cross-derivative stencils, variable face
            // coefficients and the doubled half-cell reflecting-wall gradient.
            let diffusivity = ((longitudinal + 2.0 * cross.abs() + 4.0 * mu) / minimum_density)
                .max(conductivity / minimum_cv);
            if !diffusivity.is_finite() {
                return Err(MolecularError::Invalid(
                    "molecular diffusion scale overflow",
                ));
            }
            maximum_diffusivity = maximum_diffusivity.max(diffusivity);
            let normal_derivative = if wall {
                let sign = if face == 0 { 1.0 } else { -1.0 };
                sign * 2.0 * pl.velocity[axis] / grid.spacing[axis]
            } else {
                (pr.velocity[axis] - pl.velocity[axis]) / grid.spacing[axis]
            };
            let transverse_divergence: f64 = tangential
                .into_iter()
                .map(|t| 0.5 * gradient[left][t][t] + 0.5 * gradient[right][t][t])
                .sum();
            let normal_stress = longitudinal * normal_derivative + cross * transverse_divergence;
            let shear_stress: [f64; 2] = std::array::from_fn(|i| {
                let t = tangential[i];
                if wall {
                    0.0
                } else {
                    mu * ((pr.velocity[t] - pl.velocity[t]) / grid.spacing[axis]
                        + 0.5 * gradient[left][axis][t]
                        + 0.5 * gradient[right][axis][t])
                }
            });
            let energy = if wall {
                0.0
            } else {
                -normal_stress * (0.5 * pl.velocity[axis] + 0.5 * pr.velocity[axis])
                    - shear_stress[0]
                        * (0.5 * pl.velocity[tangential[0]] + 0.5 * pr.velocity[tangential[0]])
                    - shear_stress[1]
                        * (0.5 * pl.velocity[tangential[1]] + 0.5 * pr.velocity[tangential[1]])
                    - conductivity * (pr.temperature - pl.temperature) / grid.spacing[axis]
            };
            if ![normal_stress, shear_stress[0], shear_stress[1], energy]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(MolecularError::Invalid("nonfinite molecular face flux"));
            }
            flux[index] = SpatialFaceFlux {
                mass: 0.0,
                normal_momentum: -normal_stress,
                transverse_momentum: shear_stress.map(|v| -v),
                energy,
                signal_speed: 0.0,
                used_hll_fallback: false,
            };
        }
    }
    // A conservative step estimate, not a fully discrete entropy theorem.
    // Original source admission is checked again after every conservative stage.
    let rate = if maximum_diffusivity == 0.0 {
        0.0
    } else {
        4.0 * maximum_diffusivity
            * grid
                .spacing
                .into_iter()
                .map(|h| h.recip().powi(2))
                .sum::<f64>()
    };
    if !rate.is_finite() {
        return Err(MolecularError::Invalid("parabolic scale overflow"));
    }
    Ok(rate)
}

#[cfg(test)]
#[path = "molecular_tests.rs"]
mod tests;
