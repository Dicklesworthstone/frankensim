//! Optional mean-flow inertia on the original sealed-gas capillary owner.
//!
//! The original gas supplies every midpoint pressure, endpoint pressure and
//! pressure sensitivity. Velocity-zero events consume their actual elapsed
//! heat. The existing adaptive owner resolves the remaining physical interval;
//! no composed gas receipt or substitute equation of state is invented.

use super::{
    SealedCapillaryCandidate, SealedCapillaryError, SealedCapillaryGeometry, SealedCapillaryInput,
    SealedCapillarySolve, SealedCapillaryStep, SealedGasCapillary,
};
use crate::capillary::inertia::{
    self as mechanics, InertialCapillaryInput, InertialCapillarySolve, InertialContactLine,
};
use crate::capillary::material::inertia::CapillaryInertiaParameters;
use crate::capillary::material::{MaterialCapillaryError, MaterialMeniscusInput};
use crate::capillary::{
    ContactLineMotion, InclinedTwoFluidCircularMeniscusInput, TwoFluidCircularMeniscusInput,
};
use core::f64::consts::PI;
use fs_material::thermochemical_gas::energy::GasEnergyError;

/// Independent velocity and event-time accuracy for the optional inertial state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryInertiaAccuracy {
    /// Used by both the scalar velocity solve and adaptive velocity comparison.
    pub velocity_absolute_tolerance_m_s: f64,
    /// Accuracy of a resolved numerical velocity-zero event.
    pub stopping_time_absolute_tolerance_s: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct InertiaState {
    pub parameters: CapillaryInertiaParameters,
    pub velocity_m_s: f64,
    pub accuracy: SealedCapillaryInertiaAccuracy,
}

impl<'a> SealedGasCapillary<'a> {
    /// Add explicit liquid/end inertia to this exact original gas and wetting owner.
    ///
    /// Geometry must match the sealed tube. Displaced-fluid end terms must be
    /// zero because this spatially uniform gas model carries no gas bulk inertia.
    /// The two reservoir coefficients describe the declared liquid reservoir,
    /// including its kinetic-energy exchange; no additional default loss is added.
    pub fn with_inertia(
        mut self,
        parameters: CapillaryInertiaParameters,
        initial_velocity_m_s: f64,
        accuracy: SealedCapillaryInertiaAccuracy,
    ) -> Result<Self, SealedCapillaryError> {
        let g = self.geometry;
        if parameters.radius_m != g.radius_m
            || parameters.tube_length_m != g.tube_length_m
            || !initial_velocity_m_s.is_finite()
            || initial_velocity_m_s.abs() > parameters.maximum_velocity_m_s
            || ![
                accuracy.velocity_absolute_tolerance_m_s,
                accuracy.stopping_time_absolute_tolerance_s,
            ]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
        {
            return Err(SealedCapillaryError::Invalid(
                "matching sealed geometry, admitted velocity and positive inertial accuracy required",
            ));
        }
        let rho = self.tube.liquid.density_kg_m3();
        let (added, _, _) = parameters
            .checked_coefficients(rho, 0.0)
            .map_err(parameter_error)?;
        let mass = PI * g.radius_m * g.radius_m * rho * g.wet_length_m + added;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(SealedCapillaryError::Invalid(
                "positive initial effective mass required; dry startup needs explicit added mass",
            ));
        }
        if let Some(line) = self.tube.contact_line {
            if initial_velocity_m_s.abs() > line.limits.maximum_speed_m_s {
                return Err(SealedCapillaryError::Invalid(
                    "initial inertial velocity exceeds the original contact-line source domain",
                ));
            }
        }
        self.inertia = Some(InertiaState {
            parameters,
            velocity_m_s: initial_velocity_m_s,
            accuracy,
        });
        Ok(self)
    }

    /// Physical column velocity, present only after explicit inertia admission.
    #[must_use]
    pub fn velocity_m_s(&self) -> Option<f64> {
        self.inertia.map(|state| state.velocity_m_s)
    }

    #[must_use]
    pub fn inertia_parameters(&self) -> Option<&CapillaryInertiaParameters> {
        self.inertia.as_ref().map(|state| &state.parameters)
    }

    #[must_use]
    pub fn inertia_accuracy(&self) -> Option<SealedCapillaryInertiaAccuracy> {
        self.inertia.map(|state| state.accuracy)
    }
}

fn parameter_error(error: MaterialCapillaryError) -> SealedCapillaryError {
    match error {
        MaterialCapillaryError::Unsupported(reason) => SealedCapillaryError::Invalid(reason),
        _ => SealedCapillaryError::Invalid("invalid explicit inertial end parameters"),
    }
}

fn mechanical_error(error: mechanics::InertialCapillaryError) -> SealedCapillaryError {
    use mechanics::InertialCapillaryError as Error;
    match error {
        Error::Invalid(reason) => SealedCapillaryError::Invalid(reason),
        Error::Cancelled => SealedCapillaryError::Cancelled,
        Error::StepTooLarge(reason) => SealedCapillaryError::StepTooLarge(reason),
        Error::NotConverged {
            iterations,
            residual_pa,
        } => SealedCapillaryError::NotConverged {
            iterations,
            residual_pa,
        },
        Error::EnergyBalance {
            residual_j,
            tolerance_j,
        } => SealedCapillaryError::EnergyBalance {
            residual_j,
            tolerance_j,
        },
    }
}

/// One authentic gas operation and one conservative inertial segment.
pub(super) fn step<'a>(
    owner: &SealedGasCapillary<'a>,
    input: SealedCapillaryInput,
    controls: SealedCapillarySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SealedCapillaryCandidate<'a>, SealedCapillaryError> {
    if !input.inlet_pressure_pa.is_finite()
        || input.inlet_pressure_pa < 0.0
        || !input.duration_s.is_finite()
        || input.duration_s <= 0.0
        || !controls.length_absolute_tolerance_m.is_finite()
        || controls.length_absolute_tolerance_m <= 0.0
    {
        return Err(SealedCapillaryError::Invalid(
            "nonnegative absolute inlet pressure, positive time and position accuracy required",
        ));
    }
    let inertia = owner
        .inertia
        .expect("inertial dispatch requires an admitted state");
    let g = owner.geometry;
    let rho = owner.tube.liquid.density_kg_m3();
    let area = PI * g.radius_m * g.radius_m;
    let (added, advancing, receding) = inertia
        .parameters
        .checked_coefficients(rho, 0.0)
        .map_err(parameter_error)?;
    let line = owner.tube.contact_line;
    let maximum_velocity = line.map_or(inertia.parameters.maximum_velocity_m_s, |source| {
        inertia
            .parameters
            .maximum_velocity_m_s
            .min(source.limits.maximum_speed_m_s)
    });
    // A correction in velocity also changes x by dt/2. Respect the existing
    // position tolerance as well as the independently declared velocity budget.
    let velocity_tolerance = inertia
        .accuracy
        .velocity_absolute_tolerance_m_s
        .min(2.0 * controls.length_absolute_tolerance_m / input.duration_s);
    let raw = InertialCapillaryInput {
        meniscus: InclinedTwoFluidCircularMeniscusInput {
            meniscus: TwoFluidCircularMeniscusInput {
                meniscus: owner.tube.primitive_input(MaterialMeniscusInput {
                    radius_m: g.radius_m,
                    tube_length_m: g.tube_length_m,
                    wet_length_m: g.wet_length_m,
                    applied_pressure_pa: input.inlet_pressure_pa,
                    duration_s: input.duration_s,
                }),
                displaced_dynamic_viscosity_pa_s: 0.0,
            },
            density_kg_m3: rho,
            displaced_density_kg_m3: 0.0,
            opposing_gravity_m_s2: input.opposing_gravity_m_s2,
        },
        initial_velocity_m_s: inertia.velocity_m_s,
        added_mass_kg: added,
        advancing_reservoir_loss_kg_m: advancing,
        receding_reservoir_loss_kg_m: receding,
        maximum_velocity_m_s: maximum_velocity,
        contact_line: line.map(|source| InertialContactLine {
            friction_pa_s: source.source.friction_pa_s(),
            advancing_angle_radians: source.source.advancing_angle_radians(),
            receding_angle_radians: source.source.receding_angle_radians(),
            maximum_angle_departure_radians: source.limits.maximum_angle_departure_radians,
        }),
    };
    let result = mechanics::solve_with_pressure(
        raw,
        InertialCapillarySolve {
            maximum_iterations: controls.maximum_iterations,
            pressure_absolute_tolerance_pa: controls.pressure_absolute_tolerance_pa,
            velocity_absolute_tolerance_m_s: velocity_tolerance,
            time_absolute_tolerance_s: inertia.accuracy.stopping_time_absolute_tolerance_s,
            energy_absolute_tolerance_j: controls.energy_absolute_tolerance_j,
            energy_relative_tolerance: controls.energy_relative_tolerance,
        },
        |length, elapsed, cancel| {
            let volume = owner.gas.volume_m3() - area * (length - g.wet_length_m);
            let heat = input.gas_heat_power_w * elapsed;
            if !heat.is_finite() || (input.gas_heat_power_w != 0.0 && elapsed > 0.0 && heat == 0.0)
            {
                return Err(GasEnergyError::Invalid(
                    "actual inertial segment heat is unrepresentable",
                ));
            }
            let candidate = owner
                .gas
                .midpoint_heat_volume(volume, heat, controls.gas, cancel)?;
            Ok(mechanics::Pressure {
                downstream_pressure_pa: candidate.midpoint().state().pressure,
                endpoint_pressure_pa: candidate.endpoint().state().state().pressure,
                derivative_wet_length_pa_m: -area
                    * candidate.step().pressure_final_volume_derivative_pa_m3,
                derivative_time_pa_s: input.gas_heat_power_w
                    * candidate.step().pressure_heat_derivative_pa_j,
                payload: candidate,
            })
        },
        cancelled,
    )
    .map_err(|error| match error {
        mechanics::InertialSolveError::Physics(error) => mechanical_error(error),
        mechanics::InertialSolveError::Pressure(GasEnergyError::Cancelled) => {
            SealedCapillaryError::Cancelled
        }
        mechanics::InertialSolveError::Pressure(error) => SealedCapillaryError::Gas(error),
    })?;
    let mechanical = result.step;
    // Heating can move the midpoint drive into the static band or reverse it
    // while the original pressure requires immediate depinning. A rest state
    // must begin on its original active branch; refine to resolve a later stop.
    // A moving column is deliberately free to coast against the drive.
    if inertia.velocity_m_s == 0.0 {
        let capillary_scale = 2.0 * owner.tube.surface_tension_n_m / g.radius_m;
        let initial_pressure = owner.gas.state().state().pressure;
        let initial_drive = input.inlet_pressure_pa - initial_pressure
            + capillary_scale * owner.tube.contact_angle_cosine
            - rho * input.opposing_gravity_m_s2 * g.wet_length_m;
        let (advancing_threshold, receding_threshold) = line.map_or((0.0, 0.0), |source| {
            (
                capillary_scale
                    * (owner.tube.contact_angle_cosine
                        - source.source.advancing_angle_radians().cos()),
                capillary_scale
                    * (source.source.receding_angle_radians().cos()
                        - owner.tube.contact_angle_cosine),
            )
        });
        if (initial_drive > advancing_threshold + controls.pressure_absolute_tolerance_pa
            && mechanical.motion != ContactLineMotion::Advancing)
            || (initial_drive < -receding_threshold - controls.pressure_absolute_tolerance_pa
                && mechanical.motion != ContactLineMotion::Receding)
        {
            return Err(SealedCapillaryError::StepTooLarge(
                "initial gas pressure requires a different depinning branch; refine",
            ));
        }
        if let Some(source) = line.filter(|_| mechanical.motion == ContactLineMotion::Pinned) {
            let cosine = (initial_pressure + rho * input.opposing_gravity_m_s2 * g.wet_length_m
                - input.inlet_pressure_pa)
                / capillary_scale;
            let departure =
                (cosine.clamp(-1.0, 1.0).acos() - owner.tube.contact_angle_cosine.acos()).abs();
            if departure > source.limits.maximum_angle_departure_radians {
                return Err(SealedCapillaryError::StepTooLarge(
                    "initial pinned angle exceeds the original contact-line source domain",
                ));
            }
        }
    }
    let midpoint_pressure = result.pressure.downstream_pressure_pa;
    let (gas, gas_step) = result.pressure.payload.into_parts();
    let geometry = SealedCapillaryGeometry {
        wet_length_m: mechanical.wet_length_m,
        ..g
    };
    geometry.require_volume(gas.volume_m3())?;
    let mass = rho * mechanical.transferred_volume_m3;
    let total_loss = mechanical.liquid_viscous_dissipation_j
        + mechanical.contact_line_viscous_dissipation_j
        + mechanical.pinning_dissipation_j
        + mechanical.reservoir_dissipation_j;
    let residual = (mechanical.inlet_pressure_work_j + gas_step.supplied_heat_j
        - gas_step.stored_energy_change_j)
        - mechanical.kinetic_energy_change_j
        - mechanical.wall_surface_energy_change_j
        - mechanical.gravitational_storage_change_j
        - total_loss;
    let scale = mechanical.inlet_pressure_work_j.abs()
        + gas_step.supplied_heat_j.abs()
        + gas_step.stored_energy_change_j.abs()
        + mechanical.kinetic_energy_change_j.abs()
        + mechanical.wall_surface_energy_change_j.abs()
        + mechanical.gravitational_storage_change_j.abs()
        + total_loss;
    let tolerance =
        controls.energy_absolute_tolerance_j + controls.energy_relative_tolerance * scale;
    if !mass.is_finite()
        || (mechanical.transferred_volume_m3 != 0.0 && mass == 0.0)
        || !residual.is_finite()
        || !tolerance.is_finite()
    {
        return Err(SealedCapillaryError::Invalid(
            "inertial mass or whole-system energy is unrepresentable",
        ));
    }
    if residual.abs() > tolerance {
        return Err(SealedCapillaryError::EnergyBalance {
            residual_j: residual,
            tolerance_j: tolerance,
        });
    }
    let apparent = line.map(|source| {
        let mean_velocity = 0.5 * mechanical.initial_velocity_m_s + 0.5 * mechanical.velocity_m_s;
        let cosine = match mechanical.motion {
            ContactLineMotion::Advancing => {
                source.source.advancing_angle_radians().cos()
                    - source.source.friction_pa_s() * mean_velocity / owner.tube.surface_tension_n_m
            }
            ContactLineMotion::Receding => {
                source.source.receding_angle_radians().cos()
                    - source.source.friction_pa_s() * mean_velocity / owner.tube.surface_tension_n_m
            }
            ContactLineMotion::Pinned => {
                let head = rho * input.opposing_gravity_m_s2 * g.wet_length_m;
                (midpoint_pressure + head - input.inlet_pressure_pa) * g.radius_m
                    / (2.0 * owner.tube.surface_tension_n_m)
            }
        };
        cosine.clamp(-1.0, 1.0).acos()
    });
    if cancelled() {
        return Err(SealedCapillaryError::Cancelled);
    }
    Ok(SealedCapillaryCandidate {
        endpoint: SealedGasCapillary {
            tube: owner.tube,
            geometry,
            gas,
            inertia: Some(InertiaState {
                velocity_m_s: mechanical.velocity_m_s,
                ..inertia
            }),
        },
        step: SealedCapillaryStep {
            initial_wet_length_m: g.wet_length_m,
            wet_length_m: mechanical.wet_length_m,
            elapsed_time_s: mechanical.elapsed_time_s,
            transferred_liquid_volume_m3: mechanical.transferred_volume_m3,
            transferred_liquid_mass_kg: mass,
            inlet_pressure_work_j: mechanical.inlet_pressure_work_j,
            wall_surface_energy_change_j: mechanical.wall_surface_energy_change_j,
            gravitational_energy_change_j: mechanical.gravitational_storage_change_j,
            liquid_viscous_dissipation_j: mechanical.liquid_viscous_dissipation_j,
            contact_line_viscous_dissipation_j: mechanical.contact_line_viscous_dissipation_j,
            pinning_dissipation_j: mechanical.pinning_dissipation_j,
            reservoir_dissipation_j: mechanical.reservoir_dissipation_j,
            kinetic_energy_change_j: mechanical.kinetic_energy_change_j,
            initial_velocity_m_s: Some(mechanical.initial_velocity_m_s),
            velocity_m_s: Some(mechanical.velocity_m_s),
            inertial_end: Some(mechanical.end),
            contact_line_motion: line.map(|_| mechanical.motion),
            apparent_contact_angle_radians: apparent,
            gas: gas_step,
            midpoint_gas_pressure_pa: midpoint_pressure,
            pressure_residual_pa: mechanical.pressure_residual_pa,
            energy_balance_residual_j: residual,
            iterations: mechanical.iterations,
        },
    })
}
