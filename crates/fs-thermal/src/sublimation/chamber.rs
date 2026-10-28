//! Simultaneous sublimation/deposition into one finite, thermostatted chamber.
//!
//! Every surface sees the SAME endpoint vapor partial pressure. A backward-Euler
//! mass closure couples that pressure to the original conservative per-surface
//! enthalpy solves. Warm surfaces can supply cold surfaces during one interval;
//! no surface receives an independent copy of the initial vapor inventory.
//!
//! Work is O(number of surfaces) per pressure trial, with a bounded scalar outer
//! solve and no dense matrix. An initial frozen-temperature pressure predictor
//! and safeguarded interpolation avoid a scan over the full solid inventory.
//! Allocation is bounded by max_surfaces. Refine time for first-order accuracy.
//! Gas temperature is prescribed; this module does not solve gas sensible heat,
//! molecular heat accommodation, changing geometry or mixed-phase inventories.

use super::{
    MOLAR_GAS_CONSTANT, SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw,
    SurfaceKinetics, SurfaceStep, VaporReservoir, advance_surface, coefficients, finite,
    nonnegative, positive, properties, validate_surface,
};

/// One pure-solid inventory, exposed area, and its independently selected kinetics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChamberSurface {
    pub solid: SolidSurfaceState,
    pub kinetics: SurfaceKinetics,
}

/// Numerical admission and explicit work/allocation limits, not model-error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChamberSolverConfig {
    pub surface_solver: SolverConfig,
    pub relative_mass_tolerance: f64,
    pub absolute_mass_tolerance_kg: f64,
    /// Total pressure trials, including bracket construction and refused trials.
    pub max_iterations: usize,
    pub max_surfaces: usize,
}
impl Default for ChamberSolverConfig {
    fn default() -> Self {
        Self {
            surface_solver: SolverConfig::default(),
            relative_mass_tolerance: 1.0e-10,
            absolute_mass_tolerance_kg: 1.0e-15,
            max_iterations: 96,
            max_surfaces: 1024,
        }
    }
}

/// One accepted joint exchange. Surface order matches input order.
#[derive(Clone, Debug, PartialEq)]
pub struct ChamberStep {
    /// Each surface contains the SAME final reservoir, not an independent bath.
    pub surfaces: Vec<SurfaceStep>,
    pub reservoir: VaporReservoir,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Includes sensible enthalpy on the solid laws' common original datum.
    /// The thermostatted chamber's external thermal owner receives this energy.
    pub vapor_energy_j: f64,
    /// Actual final vapor mass minus original vapor mass minus accepted transfer.
    pub mass_residual_kg: f64,
    /// Trial-pressure implied vapor mass minus original mass minus transfer.
    /// Final mass is reconstructed from accepted transfers, so this nonlinear
    /// residual cannot silently manufacture or destroy vapor inventory.
    pub nonlinear_mass_residual_kg: f64,
    /// Sum of original solid stored-energy plus vapor-energy residuals.
    pub energy_residual_j: f64,
    pub pressure_iterations: usize,
}

fn compensated_add(total: &mut f64, correction: &mut f64, value: f64) {
    let next = *total + value;
    *correction += if total.abs() >= value.abs() {
        (*total - next) + value
    } else {
        (value - next) + *total
    };
    *total = next;
}

struct Trial {
    pressure_pa: f64,
    surfaces: Vec<SurfaceStep>,
    transfer_kg: f64,
    transfer_magnitude_kg: f64,
    remaining_vapor_kg: f64,
    residual_kg: f64,
}

fn admitted(trial: &Trial, initial_mass: f64, config: ChamberSolverConfig) -> bool {
    let scale = initial_mass
        .max(trial.remaining_vapor_kg.abs())
        .max(trial.transfer_magnitude_kg);
    let tolerance = config.absolute_mass_tolerance_kg + config.relative_mass_tolerance * scale;
    tolerance.is_finite() && trial.remaining_vapor_kg >= 0.0 && trial.residual_kg.abs() <= tolerance
}

/// Advance pure-solid surfaces sharing one volatile species and one finite bath.
///
/// All caloric callbacks must retain their original enthalpy datum; surfaces of
/// the same volatile must use a common datum. The law's molar mass must be
/// temperature independent. Surface temperature comes from temperature(index,h).
/// The callback must describe a continuous nondecreasing pure-solid caloric
/// branch. No liquid fraction is inferred from temperature.
///
/// Each local solve uses the common ENDPOINT pressure, not a copied finite bath.
/// An accepted trial reconstructs final vapor mass from precisely the accepted
/// surface transfers. The separately exposed nonlinear mass residual bounds the
/// trial-pressure closure. A zero solid mass is extinction for a topology owner
/// to handle/refuse; this scalar module does not erase cells.
///
/// No input is mutated. Cancellation is checked at every nonlinear evaluation.
/// Unsupported caloric trial bounds are backed away from during bracketing.
/// If an admissible bracket cannot be established within the budget, the step
/// refuses; a smaller time interval may be needed even when inventory is finite.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber<L: SublimationLaw + ?Sized>(
    law: &L,
    surfaces: &[ChamberSurface],
    reservoir: VaporReservoir,
    duration_s: f64,
    config: ChamberSolverConfig,
    temperature: impl Fn(usize, f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ChamberStep, SublimationError> {
    checkpoint()?;
    if !nonnegative(duration_s)
        || !positive(config.relative_mass_tolerance)
        || config.relative_mass_tolerance >= 1.0
        || !nonnegative(config.absolute_mass_tolerance_kg)
        || config.max_iterations == 0
        || config.max_iterations > 4096
        || surfaces.len() > config.max_surfaces
        || !positive(config.surface_solver.relative_tolerance)
        || config.surface_solver.relative_tolerance >= 1.0
        || !nonnegative(config.surface_solver.absolute_energy_tolerance_j)
        || config.surface_solver.max_iterations == 0
        || config.surface_solver.max_iterations > 4096
    {
        return Err(SublimationError::Invalid(
            "invalid chamber solve or surface budget",
        ));
    }
    let VaporReservoir::Finite {
        vapor_mass_kg: initial_mass,
        volume_m3,
        temperature_k: gas_temperature,
    } = reservoir
    else {
        return Err(SublimationError::Invalid(
            "a shared chamber requires a finite reservoir",
        ));
    };
    if !nonnegative(initial_mass) || !positive(volume_m3) || !positive(gas_temperature) {
        return Err(SublimationError::Invalid("invalid finite chamber state"));
    }
    if surfaces.is_empty() {
        checkpoint()?;
        return Ok(ChamberStep {
            surfaces: Vec::new(),
            reservoir,
            transferred_mass_kg: 0.0,
            latent_energy_j: 0.0,
            vapor_energy_j: 0.0,
            mass_residual_kg: 0.0,
            nonlinear_mass_residual_kg: 0.0,
            energy_residual_j: 0.0,
            pressure_iterations: 0,
        });
    }
    let mut molar_mass = 0.0;
    let mut source = 0.0;
    let mut source_correction = 0.0;
    let mut conductance = 0.0;
    let mut conductance_correction = 0.0;
    let mut maximum_equilibrium_pressure: f64 = 0.0;
    for (index, surface) in surfaces.iter().enumerate() {
        checkpoint()?;
        validate_surface(surface.solid, duration_s)?;
        let t = temperature(index, surface.solid.specific_enthalpy_j_kg)?;
        let p = properties(law, t)?;
        if index == 0 {
            molar_mass = p.molar_mass_kg_mol;
        } else if p.molar_mass_kg_mol != molar_mass {
            return Err(SublimationError::Invalid(
                "shared vapor requires one constant molar mass",
            ));
        }
        let (a, b) = coefficients(p, t, gas_temperature, surface.kinetics)?;
        if surface.solid.area_m2 > 0.0 && a > 0.0 {
            if b == 0.0 {
                return Err(SublimationError::Invalid(
                    "underflowed chamber pressure conductance",
                ));
            }
            compensated_add(
                &mut source,
                &mut source_correction,
                finite(surface.solid.area_m2 * a * p.saturation_pressure_pa)?,
            );
            compensated_add(
                &mut conductance,
                &mut conductance_correction,
                finite(surface.solid.area_m2 * b)?,
            );
            maximum_equilibrium_pressure = maximum_equilibrium_pressure.max(finite(
                p.saturation_pressure_pa * (gas_temperature / t).sqrt(),
            )?);
        }
    }
    source = finite(source + source_correction)?;
    conductance = finite(conductance + conductance_correction)?;
    let beta = finite((MOLAR_GAS_CONSTANT / molar_mass) * gas_temperature / volume_m3)?;
    if !positive(beta) {
        return Err(SublimationError::Invalid(
            "unrepresentable chamber compressibility",
        ));
    }
    let initial_pressure = finite(initial_mass * beta)?;
    let upper_pressure = maximum_equilibrium_pressure.max(initial_pressure);
    let prediction = if conductance > 0.0 && duration_s > 0.0 {
        let scale = conductance * duration_s * beta;
        if scale.is_nan() {
            return Err(SublimationError::Invalid("invalid chamber time scale"));
        }
        let weight = 1.0 / (1.0 + scale);
        finite(weight * initial_pressure + (1.0 - weight) * (source / conductance))?
    } else {
        initial_pressure
    };
    let evaluations = std::cell::Cell::new(0usize);
    let mut evaluate = |pressure_pa: f64| -> Result<Trial, SublimationError> {
        checkpoint()?;
        if evaluations.get() >= config.max_iterations {
            return Err(SublimationError::NonConverged);
        }
        evaluations.set(evaluations.get() + 1);
        if !nonnegative(pressure_pa) {
            return Err(SublimationError::Invalid("invalid trial chamber pressure"));
        }
        let mut values = Vec::with_capacity(surfaces.len());
        let mut transfer = 0.0;
        let mut correction = 0.0;
        let mut magnitude = 0.0;
        for (index, surface) in surfaces.iter().enumerate() {
            let value = advance_surface(
                law,
                surface.kinetics,
                surface.solid,
                VaporReservoir::Prescribed {
                    partial_pressure_pa: pressure_pa,
                    temperature_k: gas_temperature,
                },
                duration_s,
                config.surface_solver,
                |h| temperature(index, h),
                &mut checkpoint,
            )?;
            if properties(law, value.surface_temperature_k)?.molar_mass_kg_mol != molar_mass {
                return Err(SublimationError::Invalid(
                    "volatile molar mass changed during phase exchange",
                ));
            }
            compensated_add(&mut transfer, &mut correction, value.transferred_mass_kg);
            magnitude = finite(magnitude + value.transferred_mass_kg.abs())?;
            values.push(value);
        }
        let transfer = finite(transfer + correction)?;
        let remaining = finite(initial_mass + transfer)?;
        // Difference form avoids subtracting large pressure-derived inventories.
        let residual = finite((pressure_pa - initial_pressure) / beta - transfer)?;
        Ok(Trial {
            pressure_pa,
            surfaces: values,
            transfer_kg: transfer,
            transfer_magnitude_kg: magnitude,
            remaining_vapor_kg: remaining,
            residual_kg: residual,
        })
    };
    let initial = evaluate(prediction)?;
    let endpoint = if admitted(&initial, initial_mass, config) {
        initial
    } else {
        let sign = initial.residual_kg.signum();
        let mut near = initial;
        let mut candidate_pressure = if sign > 0.0 { 0.0 } else { upper_pressure };
        let mut invalid_bound = None;
        let far = loop {
            if candidate_pressure == near.pressure_pa {
                return Err(SublimationError::NonConverged);
            }
            match evaluate(candidate_pressure) {
                Ok(value) => {
                    if admitted(&value, initial_mass, config) || value.residual_kg.signum() != sign
                    {
                        break value;
                    }
                    near = value;
                    // The initial-temperature equilibrium envelope normally
                    // brackets the cooperative pure-solid response. Expand
                    // only if a caller's admissible constitutive law needs it.
                    candidate_pressure = if let Some(bound) = invalid_bound {
                        near.pressure_pa + (bound - near.pressure_pa) * 0.5
                    } else if sign < 0.0 {
                        finite(near.pressure_pa * 2.0 + initial_pressure.max(1.0))?
                    } else {
                        near.pressure_pa * 0.5
                    };
                }
                Err(SublimationError::Cancelled) => return Err(SublimationError::Cancelled),
                Err(error) => {
                    if evaluations.get() >= config.max_iterations {
                        return Err(error);
                    }
                    invalid_bound = Some(candidate_pressure);
                    candidate_pressure =
                        near.pressure_pa + (candidate_pressure - near.pressure_pa) * 0.5;
                }
            }
        };
        if admitted(&far, initial_mass, config) {
            far
        } else {
            let (mut left, mut right) = if near.pressure_pa < far.pressure_pa {
                (near, far)
            } else {
                (far, near)
            };
            let mut count = 0;
            loop {
                let span = right.pressure_pa - left.pressure_pa;
                let fraction = -left.residual_kg / (right.residual_kg - left.residual_kg);
                // Every third trial bisects, bounding interval contraction even
                // when the secant is repeatedly attracted to one endpoint.
                let fraction = if count % 3 == 2 || !(0.05..=0.95).contains(&fraction) {
                    0.5
                } else {
                    fraction
                };
                let p = left.pressure_pa + span * fraction;
                if p == left.pressure_pa || p == right.pressure_pa {
                    return Err(SublimationError::NonConverged);
                }
                let middle = evaluate(p)?;
                if admitted(&middle, initial_mass, config) {
                    break middle;
                }
                if middle.residual_kg.signum() == left.residual_kg.signum() {
                    left = middle;
                } else {
                    right = middle;
                }
                count += 1;
            }
        }
    };
    drop(evaluate);
    checkpoint()?;
    let next_reservoir = VaporReservoir::Finite {
        vapor_mass_kg: endpoint.remaining_vapor_kg,
        volume_m3,
        temperature_k: gas_temperature,
    };
    let mut latent = 0.0;
    let mut latent_correction = 0.0;
    let mut energy = 0.0;
    let mut energy_correction = 0.0;
    let mut residual = 0.0;
    let mut residual_correction = 0.0;
    let mut values = endpoint.surfaces;
    for value in &mut values {
        value.reservoir = next_reservoir;
        compensated_add(&mut latent, &mut latent_correction, value.latent_energy_j);
        compensated_add(&mut energy, &mut energy_correction, value.vapor_energy_j);
        compensated_add(
            &mut residual,
            &mut residual_correction,
            value.energy_residual_j,
        );
    }
    Ok(ChamberStep {
        surfaces: values,
        reservoir: next_reservoir,
        transferred_mass_kg: endpoint.transfer_kg,
        latent_energy_j: finite(latent + latent_correction)?,
        vapor_energy_j: finite(energy + energy_correction)?,
        mass_residual_kg: finite(
            endpoint.remaining_vapor_kg - initial_mass - endpoint.transfer_kg,
        )?,
        nonlinear_mass_residual_kg: endpoint.residual_kg,
        energy_residual_j: finite(residual + residual_correction)?,
        pressure_iterations: evaluations.get(),
    })
}

#[cfg(test)]
mod tests;
