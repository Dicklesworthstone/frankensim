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
    SublimationProperties, SurfaceKinetics, SurfaceStep, VaporReservoir, advance_surface_with_heat,
    coefficients, finite, nonnegative, positive, properties, validate_surface,
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
    /// Signed heat supplied to all solids, counted once in their local balances.
    pub supplied_heat_j: f64,
    /// Sum of local storage plus carried vapor energy minus supplied heat.
    pub energy_residual_j: f64,
    pub pressure_iterations: usize,
}

pub(crate) fn compensated_add(total: &mut f64, correction: &mut f64, value: f64) {
    let next = *total + value;
    *correction += if total.abs() >= value.abs() {
        (*total - next) + value
    } else {
        (value - next) + *total
    };
    *total = next;
}

pub(crate) struct InitialSurface {
    pub temperature_k: f64,
    pub properties: SublimationProperties,
    pub area_m2: f64,
    pub kinetics: SurfaceKinetics,
    /// The predictor must not invent emission from an initially bare substrate.
    pub has_deposit: bool,
}

/// Only the physical phase account is shared. A supported deposit also retains
/// its independent substrate storage in the original local result.
pub(crate) trait ChamberExchange {
    fn phase(&self) -> &SurfaceStep;
    fn phase_mut(&mut self) -> &mut SurfaceStep;
}

impl ChamberExchange for SurfaceStep {
    fn phase(&self) -> &SurfaceStep {
        self
    }
    fn phase_mut(&mut self) -> &mut SurfaceStep {
        self
    }
}

pub(crate) struct ChamberSolution<T> {
    pub surfaces: Vec<T>,
    pub reservoir: VaporReservoir,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    pub mass_residual_kg: f64,
    pub nonlinear_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub pressure_iterations: usize,
}

struct Trial<T> {
    pressure_pa: f64,
    surfaces: Vec<T>,
    transfer_kg: f64,
    transfer_magnitude_kg: f64,
    remaining_vapor_kg: f64,
    residual_kg: f64,
}

fn admitted<T>(trial: &Trial<T>, initial_mass: f64, config: ChamberSolverConfig) -> bool {
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
    checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ChamberStep, SublimationError> {
    advance_chamber_impl(
        law,
        surfaces,
        reservoir,
        duration_s,
        |_| 0.0,
        config,
        temperature,
        checkpoint,
    )
}

/// Advance one finite chamber with signed heat supplied to each solid during
/// the interval. Heat entries follow the INPUT surface order; positive values
/// heat the solid and negative values cool it. Every entry must be finite.
///
/// Heat, latent transfer and shared endpoint pressure are solved together.
/// Each original caloric chart remains authoritative, and each local residual
/// is `m0*(h1-h0) + dm*L - Q`. The complete carried vapor energy includes
/// `dm*h1`; latent energy is not an additional heat source. A zero-duration
/// interval admits caloric impulses but no mass exchange. No input is mutated.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_with_heat<L: SublimationLaw + ?Sized>(
    law: &L,
    surfaces: &[ChamberSurface],
    reservoir: VaporReservoir,
    duration_s: f64,
    supplied_heat_j: &[f64],
    config: ChamberSolverConfig,
    temperature: impl Fn(usize, f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ChamberStep, SublimationError> {
    checkpoint()?;
    if supplied_heat_j.len() != surfaces.len() {
        return Err(SublimationError::Invalid(
            "chamber heat must match the input surface count",
        ));
    }
    advance_chamber_impl(
        law,
        surfaces,
        reservoir,
        duration_s,
        |index| supplied_heat_j[index],
        config,
        temperature,
        checkpoint,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_chamber_impl<L: SublimationLaw + ?Sized>(
    law: &L,
    surfaces: &[ChamberSurface],
    reservoir: VaporReservoir,
    duration_s: f64,
    supplied_heat_j: impl Fn(usize) -> f64,
    config: ChamberSolverConfig,
    temperature: impl Fn(usize, f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ChamberStep, SublimationError> {
    let solution = solve_chamber(
        surfaces.len(),
        reservoir,
        duration_s,
        config,
        |index| {
            let surface = surfaces[index];
            if !supplied_heat_j(index).is_finite() {
                return Err(SublimationError::Invalid("nonfinite supplied chamber heat"));
            }
            validate_surface(surface.solid, duration_s)?;
            let t = temperature(index, surface.solid.specific_enthalpy_j_kg)?;
            Ok(InitialSurface {
                temperature_k: t,
                properties: properties(law, t)?,
                area_m2: surface.solid.area_m2,
                kinetics: surface.kinetics,
                has_deposit: true,
            })
        },
        |index, reservoir, checkpoint| {
            let surface = surfaces[index];
            let value = advance_surface_with_heat(
                law,
                surface.kinetics,
                surface.solid,
                reservoir,
                duration_s,
                supplied_heat_j(index),
                config.surface_solver,
                |h| temperature(index, h),
                checkpoint,
            )?;
            let molar_mass = properties(law, value.surface_temperature_k)?.molar_mass_kg_mol;
            Ok((value, molar_mass))
        },
        &mut checkpoint,
    )?;
    let mut supplied = 0.0;
    let mut supplied_correction = 0.0;
    for (index, surface) in solution.surfaces.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        compensated_add(
            &mut supplied,
            &mut supplied_correction,
            surface.supplied_heat_j,
        );
    }
    Ok(ChamberStep {
        surfaces: solution.surfaces,
        reservoir: solution.reservoir,
        transferred_mass_kg: solution.transferred_mass_kg,
        latent_energy_j: solution.latent_energy_j,
        vapor_energy_j: solution.vapor_energy_j,
        mass_residual_kg: solution.mass_residual_kg,
        nonlinear_mass_residual_kg: solution.nonlinear_mass_residual_kg,
        supplied_heat_j: finite(supplied + supplied_correction)?,
        energy_residual_j: solution.energy_residual_j,
        pressure_iterations: solution.pressure_iterations,
    })
}

/// One bounded pressure closure for free solids and supported deposits. The
/// local callback owns its original storage law and returns the endpoint molar
/// mass so a changing species identity cannot pass the shared mass balance.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_chamber<T: ChamberExchange, C: FnMut() -> Result<(), SublimationError>>(
    surface_count: usize,
    reservoir: VaporReservoir,
    duration_s: f64,
    config: ChamberSolverConfig,
    initial_surface: impl Fn(usize) -> Result<InitialSurface, SublimationError>,
    mut evaluate_surface: impl FnMut(
        usize,
        VaporReservoir,
        &mut C,
    ) -> Result<(T, f64), SublimationError>,
    mut checkpoint: C,
) -> Result<ChamberSolution<T>, SublimationError> {
    checkpoint()?;
    if !nonnegative(duration_s)
        || !positive(config.relative_mass_tolerance)
        || config.relative_mass_tolerance >= 1.0
        || !nonnegative(config.absolute_mass_tolerance_kg)
        || config.max_iterations == 0
        || config.max_iterations > 4096
        || surface_count > config.max_surfaces
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
    if surface_count == 0 {
        checkpoint()?;
        return Ok(ChamberSolution {
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
    for index in 0..surface_count {
        checkpoint()?;
        let surface = initial_surface(index)?;
        let t = surface.temperature_k;
        let p = surface.properties;
        if index == 0 {
            molar_mass = p.molar_mass_kg_mol;
        } else if p.molar_mass_kg_mol != molar_mass {
            return Err(SublimationError::Invalid(
                "shared vapor requires one constant molar mass",
            ));
        }
        let (a, b) = coefficients(p, t, gas_temperature, surface.kinetics)?;
        if surface.area_m2 > 0.0 && a > 0.0 {
            if b == 0.0 {
                return Err(SublimationError::Invalid(
                    "underflowed chamber pressure conductance",
                ));
            }
            if surface.has_deposit {
                compensated_add(
                    &mut source,
                    &mut source_correction,
                    finite(surface.area_m2 * a * p.saturation_pressure_pa)?,
                );
            }
            compensated_add(
                &mut conductance,
                &mut conductance_correction,
                finite(surface.area_m2 * b)?,
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
    let mut evaluate = |pressure_pa: f64| -> Result<Trial<T>, SublimationError> {
        checkpoint()?;
        if evaluations.get() >= config.max_iterations {
            return Err(SublimationError::NonConverged);
        }
        evaluations.set(evaluations.get() + 1);
        if !nonnegative(pressure_pa) {
            return Err(SublimationError::Invalid("invalid trial chamber pressure"));
        }
        let mut values = Vec::with_capacity(surface_count);
        let mut transfer = 0.0;
        let mut correction = 0.0;
        let mut magnitude = 0.0;
        for index in 0..surface_count {
            let (value, endpoint_molar_mass) = evaluate_surface(
                index,
                VaporReservoir::Prescribed {
                    partial_pressure_pa: pressure_pa,
                    temperature_k: gas_temperature,
                },
                &mut checkpoint,
            )?;
            if endpoint_molar_mass != molar_mass {
                return Err(SublimationError::Invalid(
                    "volatile molar mass changed during phase exchange",
                ));
            }
            compensated_add(
                &mut transfer,
                &mut correction,
                value.phase().transferred_mass_kg,
            );
            magnitude = finite(magnitude + value.phase().transferred_mass_kg.abs())?;
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
    let initial = match evaluate(prediction) {
        Ok(value) => value,
        Err(SublimationError::Cancelled) => return Err(SublimationError::Cancelled),
        Err(initial_error) => {
            // This frozen-temperature predictor omits supplied heat and the
            // latent buffering of a newly forming deposit. Its local caloric
            // endpoint can therefore leave the material domain even when the
            // joint pressure/temperature solution is admissible. Find one
            // physical pressure anchor before the ordinary bracket solve;
            // none of these trials changes the finite inventory or charts.
            // Existing max_iterations charges every refused pressure trial.
            let pressures = [initial_pressure, 0.0, upper_pressure];
            let mut anchor = None;
            for (index, pressure) in pressures.into_iter().enumerate() {
                if evaluations.get() >= config.max_iterations {
                    break;
                }
                if pressure == prediction || pressures[..index].contains(&pressure) {
                    continue;
                }
                match evaluate(pressure) {
                    Ok(value) => {
                        anchor = Some(value);
                        break;
                    }
                    Err(SublimationError::Cancelled) => return Err(SublimationError::Cancelled),
                    Err(_) => {}
                }
            }
            anchor.ok_or(initial_error)?
        }
    };
    let endpoint = if admitted(&initial, initial_mass, config) {
        initial
    } else {
        let sign = initial.residual_kg.signum();
        let mut near = initial;
        let mut candidate_pressure = if sign > 0.0 {
            0.0
        } else if upper_pressure > near.pressure_pa {
            upper_pressure
        } else {
            // Supplied heat can raise the endpoint above every INITIAL
            // equilibrium pressure, including an initially equilibrated bath.
            finite(near.pressure_pa * 2.0 + initial_pressure.max(1.0))?
        };
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
    for exchange in &mut values {
        let value = exchange.phase_mut();
        value.reservoir = next_reservoir;
        compensated_add(&mut latent, &mut latent_correction, value.latent_energy_j);
        compensated_add(&mut energy, &mut energy_correction, value.vapor_energy_j);
        compensated_add(
            &mut residual,
            &mut residual_correction,
            value.energy_residual_j,
        );
    }
    Ok(ChamberSolution {
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
