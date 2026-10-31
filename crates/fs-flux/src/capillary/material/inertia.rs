//! Source-owned inertial capillary motion through an open circular tube.
//!
//! Material coefficients remain at their original admitted state. Added end
//! inertance and directional reservoir exchange are explicit geometry/model
//! inputs, never material-name defaults. The original displaced-fluid owner is
//! borrowed and its transport is evaluated once at binding.
//!
//! Complete steps resolve discrete stopping events and restart from rest.
//! They report kinetic storage and each dissipative channel separately. These
//! are lumped inertance dynamics with Poiseuille drag, not resolved unsteady
//! velocity profiles, deposited films or meniscus end-cap topology.

pub mod advance;

use super::{MaterialCapillaryError, MaterialCapillaryTube};
use crate::capillary::inertia as mechanics;
use crate::capillary::{
    CircularMeniscusInput, ContactLineMotion, InclinedTwoFluidCircularMeniscusInput,
    TwoFluidCircularMeniscusInput,
};
use core::f64::consts::PI;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Explicit tube geometry and reservoir/end models, with no inferred defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryInertiaParameters {
    pub radius_m: f64,
    pub tube_length_m: f64,
    /// Geometry-equivalent added inertance length of the invading liquid.
    pub liquid_added_length_m: f64,
    /// Geometry-equivalent added length of the displaced phase.
    pub displaced_added_length_m: f64,
    pub liquid_entry_loss_coefficient: f64,
    pub liquid_exit_loss_coefficient: f64,
    pub displaced_entry_loss_coefficient: f64,
    pub displaced_exit_loss_coefficient: f64,
    /// Explicit operating speed limit for this reduced inertance model.
    pub maximum_velocity_m_s: f64,
}

impl CapillaryInertiaParameters {
    /// Geometry and original densities determine mass and force coefficients.
    /// A zero displaced density explicitly denotes the massless vented phase.
    pub(super) fn checked_coefficients(
        self,
        density: f64,
        displaced_density: f64,
    ) -> Result<(f64, f64, f64), MaterialCapillaryError> {
        let invalid = MaterialCapillaryError::Unsupported;
        if ![
            self.radius_m,
            self.tube_length_m,
            self.maximum_velocity_m_s,
            density,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
            || ![
                self.liquid_added_length_m,
                self.displaced_added_length_m,
                self.liquid_entry_loss_coefficient,
                self.liquid_exit_loss_coefficient,
                self.displaced_entry_loss_coefficient,
                self.displaced_exit_loss_coefficient,
                displaced_density,
            ]
            .into_iter()
            .all(|v| v.is_finite() && v >= 0.0)
        {
            return Err(invalid(
                "finite positive inertial geometry/speed/density and nonnegative end coefficients required",
            ));
        }
        if displaced_density == 0.0
            && (self.displaced_added_length_m != 0.0
                || self.displaced_entry_loss_coefficient != 0.0
                || self.displaced_exit_loss_coefficient != 0.0)
        {
            return Err(invalid(
                "a massless vented phase cannot own displaced-fluid end inertance or loss",
            ));
        }
        let area = PI * self.radius_m * self.radius_m;
        let mass = area
            * (density * self.liquid_added_length_m
                + displaced_density * self.displaced_added_length_m);
        let advancing = 0.5
            * area
            * (density * self.liquid_entry_loss_coefficient
                + displaced_density * self.displaced_exit_loss_coefficient);
        let receding = 0.5
            * area
            * (density * self.liquid_exit_loss_coefficient
                + displaced_density * self.displaced_entry_loss_coefficient);
        if !area.is_finite()
            || area <= 0.0
            || ![mass, advancing, receding]
                .into_iter()
                .all(|v| v.is_finite() && v >= 0.0)
            || (mass == 0.0
                && (self.liquid_added_length_m > 0.0 || self.displaced_added_length_m > 0.0))
            || (advancing == 0.0
                && (self.liquid_entry_loss_coefficient > 0.0
                    || self.displaced_exit_loss_coefficient > 0.0))
            || (receding == 0.0
                && (self.liquid_exit_loss_coefficient > 0.0
                    || self.displaced_entry_loss_coefficient > 0.0))
        {
            return Err(invalid(
                "source-derived end inertance or reservoir force is unrepresentable",
            ));
        }
        Ok((mass, advancing, receding))
    }
}

/// Physical position and signed velocity, not an integrator auxiliary velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialMeniscusState {
    pub wet_length_m: f64,
    pub velocity_m_s: f64,
}

/// Constant boundary values over one requested interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialCapillaryBoundary {
    /// Inlet minus downstream reservoir pressure, before hydrostatic head.
    pub applied_pressure_pa: f64,
    pub opposing_gravity_m_s2: f64,
    pub duration_s: f64,
}

/// Bounded event resolution around the existing inertial scalar solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialCapillaryIntervalControls {
    pub mechanics: mechanics::InertialCapillarySolve,
    /// Maximum discrete velocity-zero events within one complete step.
    pub maximum_stopping_events: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MaterialInertiaError {
    Material(MaterialCapillaryError),
    Mechanics(mechanics::InertialCapillaryError),
    Invalid(&'static str),
    Cancelled,
    StoppingBudget { events: usize },
    TimeStagnation,
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl core::fmt::Display for MaterialInertiaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source inertial capillary: {self:?}")
    }
}
impl std::error::Error for MaterialInertiaError {}
impl From<MaterialCapillaryError> for MaterialInertiaError {
    fn from(value: MaterialCapillaryError) -> Self {
        Self::Material(value)
    }
}
impl From<mechanics::InertialCapillaryError> for MaterialInertiaError {
    fn from(value: mechanics::InertialCapillaryError) -> Self {
        match value {
            mechanics::InertialCapillaryError::Cancelled => Self::Cancelled,
            other => Self::Mechanics(other),
        }
    }
}

/// Immutable original material owners and one physical meniscus state.
#[derive(Clone, Copy)]
pub struct SourceInertialCapillary<'a, 'd> {
    tube: MaterialCapillaryTube<'a>,
    displaced_fluid: Option<&'d dyn FluidTransport>,
    displaced_transport: Option<FluidTransportState>,
    parameters: CapillaryInertiaParameters,
    state: InertialMeniscusState,
    added_mass_kg: f64,
    advancing_reservoir_loss_kg_m: f64,
    receding_reservoir_loss_kg_m: f64,
}
impl core::fmt::Debug for SourceInertialCapillary<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SourceInertialCapillary")
            .field("tube", &self.tube)
            .field("displaced_transport", &self.displaced_transport)
            .field("parameters", &self.parameters)
            .field("state", &self.state)
            .finish_non_exhaustive()
    }
}

/// Accepted full-interval accounts, reconstructed from physical endpoint states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialInertialCapillaryStep {
    pub initial_state: InertialMeniscusState,
    pub state: InertialMeniscusState,
    pub elapsed_time_s: f64,
    pub transferred_volume_m3: f64,
    pub transferred_liquid_mass_kg: f64,
    pub displaced_mass_kg: f64,
    pub initial_kinetic_energy_j: f64,
    pub kinetic_energy_j: f64,
    pub kinetic_energy_change_j: f64,
    pub pressure_work_j: f64,
    pub wall_surface_energy_change_j: f64,
    pub gravitational_storage_change_j: f64,
    pub outlet_potential_energy_j: f64,
    pub liquid_viscous_dissipation_j: f64,
    pub displaced_viscous_dissipation_j: f64,
    pub contact_line_viscous_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub reservoir_dissipation_j: f64,
    pub total_dissipation_j: f64,
    pub energy_balance_residual_j: f64,
    pub initial_apparent_angle_radians: f64,
    pub final_apparent_angle_radians: f64,
    pub final_motion: ContactLineMotion,
    pub stopping_events: usize,
    /// First accepted discrete stop, measured from this interval's start.
    /// Adaptivity uses this to avoid hiding the same early stop in both trials.
    pub first_stopping_time_s: Option<f64>,
    pub segments: usize,
    pub iterations: usize,
    pub maximum_pressure_residual_pa: f64,
}

#[derive(Clone, Debug)]
pub struct MaterialInertialCapillaryCandidate<'a, 'd> {
    endpoint: SourceInertialCapillary<'a, 'd>,
    step: MaterialInertialCapillaryStep,
}
impl<'a, 'd> MaterialInertialCapillaryCandidate<'a, 'd> {
    #[must_use]
    pub const fn endpoint(&self) -> &SourceInertialCapillary<'a, 'd> {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &MaterialInertialCapillaryStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        SourceInertialCapillary<'a, 'd>,
        MaterialInertialCapillaryStep,
    ) {
        (self.endpoint, self.step)
    }
}

impl<'a> MaterialCapillaryTube<'a> {
    /// Add physical velocity and explicitly declared liquid/reservoir inertance.
    /// The exterior phase has negligible inertia and hydraulic resistance.
    pub fn with_inertia(
        self,
        parameters: CapillaryInertiaParameters,
        state: InertialMeniscusState,
    ) -> Result<SourceInertialCapillary<'a, 'static>, MaterialInertiaError> {
        SourceInertialCapillary::bind(self, parameters, state, None)
    }

    /// Bind both original fluids at their frozen transport/source states.
    /// The caller assigns the actual fluid pair admitted by tension and wetting.
    pub fn with_inertia_two_fluid<'d>(
        self,
        parameters: CapillaryInertiaParameters,
        state: InertialMeniscusState,
        displaced_fluid: &'d dyn FluidTransport,
    ) -> Result<SourceInertialCapillary<'a, 'd>, MaterialInertiaError> {
        SourceInertialCapillary::bind(self, parameters, state, Some(displaced_fluid))
    }
}

impl<'a, 'd> SourceInertialCapillary<'a, 'd> {
    fn bind(
        tube: MaterialCapillaryTube<'a>,
        parameters: CapillaryInertiaParameters,
        state: InertialMeniscusState,
        displaced_fluid: Option<&'d dyn FluidTransport>,
    ) -> Result<Self, MaterialInertiaError> {
        let displaced_transport = displaced_fluid
            .map(|source| source.transport_state())
            .transpose()
            .map_err(MaterialCapillaryError::Transport)?;
        let density = tube.liquid().density_kg_m3();
        let density_d = displaced_transport.map_or(0.0, |s| s.density_kg_m3());
        let (mass, advancing, receding) = parameters.checked_coefficients(density, density_d)?;
        let speed_limit =
            tube.contact_line_limits()
                .map_or(parameters.maximum_velocity_m_s, |limits| {
                    parameters
                        .maximum_velocity_m_s
                        .min(limits.maximum_speed_m_s)
                });
        let area = PI * parameters.radius_m * parameters.radius_m;
        let moving_mass = mass
            + area
                * (density * state.wet_length_m
                    + density_d * (parameters.tube_length_m - state.wet_length_m));
        if !state.wet_length_m.is_finite()
            || state.wet_length_m < 0.0
            || state.wet_length_m > parameters.tube_length_m
            || !state.velocity_m_s.is_finite()
            || state.velocity_m_s.abs() > speed_limit
            || !moving_mass.is_finite()
            || moving_mass <= 0.0
        {
            return Err(MaterialInertiaError::Invalid(
                "in-tube state, admitted speed and positive physical moving mass required",
            ));
        }
        Ok(Self {
            tube,
            displaced_fluid,
            displaced_transport,
            parameters,
            state,
            added_mass_kg: mass,
            advancing_reservoir_loss_kg_m: advancing,
            receding_reservoir_loss_kg_m: receding,
        })
    }

    #[must_use]
    pub const fn tube(&self) -> &MaterialCapillaryTube<'a> {
        &self.tube
    }
    #[must_use]
    pub const fn parameters(&self) -> CapillaryInertiaParameters {
        self.parameters
    }
    #[must_use]
    pub const fn state(&self) -> InertialMeniscusState {
        self.state
    }
    #[must_use]
    pub const fn displaced_fluid(&self) -> Option<&'d dyn FluidTransport> {
        self.displaced_fluid
    }
    #[must_use]
    pub const fn displaced_transport(&self) -> Option<&FluidTransportState> {
        self.displaced_transport.as_ref()
    }

    fn raw_input(&self, boundary: InertialCapillaryBoundary) -> mechanics::InertialCapillaryInput {
        let line = self.tube.contact_line_source();
        let line_limits = self.tube.contact_line_limits();
        mechanics::InertialCapillaryInput {
            meniscus: InclinedTwoFluidCircularMeniscusInput {
                meniscus: TwoFluidCircularMeniscusInput {
                    meniscus: CircularMeniscusInput {
                        dynamic_viscosity_pa_s: self.tube.liquid().dynamic_viscosity_pa_s(),
                        radius_m: self.parameters.radius_m,
                        tube_length_m: self.parameters.tube_length_m,
                        wet_length_m: self.state.wet_length_m,
                        surface_tension_n_m: self.tube.surface_tension_n_m,
                        contact_angle_cosine: self.tube.contact_angle_cosine,
                        applied_pressure_pa: boundary.applied_pressure_pa,
                        duration_s: boundary.duration_s,
                    },
                    displaced_dynamic_viscosity_pa_s: self
                        .displaced_transport
                        .map_or(0.0, |s| s.dynamic_viscosity_pa_s()),
                },
                density_kg_m3: self.tube.liquid().density_kg_m3(),
                displaced_density_kg_m3: self
                    .displaced_transport
                    .map_or(0.0, |s| s.density_kg_m3()),
                opposing_gravity_m_s2: boundary.opposing_gravity_m_s2,
            },
            initial_velocity_m_s: self.state.velocity_m_s,
            added_mass_kg: self.added_mass_kg,
            advancing_reservoir_loss_kg_m: self.advancing_reservoir_loss_kg_m,
            receding_reservoir_loss_kg_m: self.receding_reservoir_loss_kg_m,
            maximum_velocity_m_s: line_limits.map_or(
                self.parameters.maximum_velocity_m_s,
                |limits| {
                    self.parameters
                        .maximum_velocity_m_s
                        .min(limits.maximum_speed_m_s)
                },
            ),
            contact_line: line.zip(line_limits).map(|(source, limits)| {
                mechanics::InertialContactLine {
                    friction_pa_s: source.friction_pa_s(),
                    advancing_angle_radians: source.advancing_angle_radians(),
                    receding_angle_radians: source.receding_angle_radians(),
                    maximum_angle_departure_radians: limits.maximum_angle_departure_radians,
                }
            }),
        }
    }

    /// Resolve all admitted velocity-zero events and consume the complete interval.
    /// Tube ends, cancellation and exhausted work refuse the entire candidate.
    /// The event times belong to the midpoint discretization; temporal accuracy
    /// requires refinement or the adaptive interval method.
    pub fn step_candidate(
        &self,
        boundary: InertialCapillaryBoundary,
        controls: InertialCapillaryIntervalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialInertialCapillaryCandidate<'a, 'd>, MaterialInertiaError> {
        poll(&mut cancelled)?;
        validate_boundary(boundary, controls)?;
        let mut endpoint = *self;
        let mut elapsed = 0.0;
        let mut sums = [Sum::default(); 10];
        let mut stops = 0usize;
        let mut first_stop = None;
        let mut segments = 0usize;
        let mut iterations = 0usize;
        let mut max_residual = 0.0_f64;
        let mut initial_angle = None;
        let last = loop {
            poll(&mut cancelled)?;
            let remaining = boundary.duration_s - elapsed;
            let step = mechanics::step_inertial_circular_meniscus(
                endpoint.raw_input(InertialCapillaryBoundary {
                    duration_s: remaining,
                    ..boundary
                }),
                controls.mechanics,
                &mut cancelled,
            )?;
            if step.elapsed_time_s > remaining
                || step.elapsed_time_s < 0.0
                || (remaining > 0.0 && step.elapsed_time_s == 0.0)
            {
                return Err(MaterialInertiaError::TimeStagnation);
            }
            let next = if step.elapsed_time_s == remaining {
                boundary.duration_s
            } else {
                elapsed + step.elapsed_time_s
            };
            if remaining > 0.0 && (next <= elapsed || !next.is_finite()) {
                return Err(MaterialInertiaError::TimeStagnation);
            }
            initial_angle.get_or_insert(step.initial_apparent_angle_radians);
            accumulate_raw(&mut sums, &step);
            iterations = iterations
                .checked_add(step.iterations)
                .ok_or(MaterialInertiaError::Invalid("iteration count overflow"))?;
            max_residual = max_residual.max(step.pressure_residual_pa.abs());
            segments += 1;
            endpoint.state = InertialMeniscusState {
                wet_length_m: step.wet_length_m,
                velocity_m_s: step.velocity_m_s,
            };
            elapsed = next;
            if step.end == mechanics::InertialCapillaryEnd::Stopped {
                stops += 1;
                first_stop.get_or_insert(elapsed);
                if stops > controls.maximum_stopping_events {
                    return Err(MaterialInertiaError::StoppingBudget { events: stops });
                }
            } else if elapsed != boundary.duration_s {
                return Err(MaterialInertiaError::Invalid(
                    "incomplete segment without a stopping event",
                ));
            }
            if elapsed == boundary.duration_s {
                break step;
            }
        };
        let step = self.finish_step(
            endpoint,
            boundary,
            sums,
            controls.mechanics,
            initial_angle.expect("at least one raw segment is evaluated"),
            last.final_apparent_angle_radians,
            last.motion,
            stops,
            first_stop,
            segments,
            iterations,
            max_residual,
        )?;
        poll(&mut cancelled)?;
        Ok(MaterialInertialCapillaryCandidate { endpoint, step })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_step(
        &self,
        endpoint: Self,
        boundary: InertialCapillaryBoundary,
        sums: [Sum; 10],
        controls: mechanics::InertialCapillarySolve,
        initial_angle: f64,
        final_angle: f64,
        final_motion: ContactLineMotion,
        stops: usize,
        first_stop: Option<f64>,
        segments: usize,
        iterations: usize,
        max_residual: f64,
    ) -> Result<MaterialInertialCapillaryStep, MaterialInertiaError> {
        let v = sums.map(Sum::value);
        let rho = self.tube.liquid().density_kg_m3();
        let rho_d = self.displaced_transport.map_or(0.0, |s| s.density_kg_m3());
        let area = PI * self.parameters.radius_m * self.parameters.radius_m;
        let x0 = self.state.wet_length_m;
        let x1 = endpoint.state.wet_length_m;
        let u0 = self.state.velocity_m_s;
        let u1 = endpoint.state.velocity_m_s;
        let mass0 =
            self.added_mass_kg + area * (rho * x0 + rho_d * (self.parameters.tube_length_m - x0));
        let mass_change = area * (rho - rho_d) * (x1 - x0);
        let kinetic0 = 0.5 * mass0 * u0 * u0;
        let mass1 =
            self.added_mass_kg + area * (rho * x1 + rho_d * (self.parameters.tube_length_m - x1));
        let kinetic1 = 0.5 * mass1 * u1 * u1;
        let mean_velocity = 0.5 * u0 + 0.5 * u1;
        let kinetic_change = if mass_change >= 0.0 {
            mass0 * mean_velocity * (u1 - u0) + 0.5 * mass_change * u1 * u1
        } else {
            mass1 * mean_velocity * (u1 - u0) + 0.5 * mass_change * u0 * u0
        };
        // Reconstruct net material inventory and conservative work from the
        // actual endpoints, independently of the accepted segment reductions.
        let volume = area * (x1 - x0);
        let pressure_work = boundary.applied_pressure_pa * volume;
        let wall_energy = -(2.0 * self.tube.surface_tension_n_m * self.tube.contact_angle_cosine
            / self.parameters.radius_m)
            * volume;
        let gravity_energy =
            (rho - rho_d) * boundary.opposing_gravity_m_s2 * (0.5 * x0 + 0.5 * x1) * volume;
        let outlet_potential =
            rho_d * boundary.opposing_gravity_m_s2 * self.parameters.tube_length_m * volume;
        let total_loss = v[5] + v[6] + v[7] + v[8] + v[9];
        let residual = (pressure_work - wall_energy - gravity_energy - outlet_potential)
            - kinetic_change
            - total_loss;
        let scale = pressure_work.abs()
            + wall_energy.abs()
            + gravity_energy.abs()
            + outlet_potential.abs()
            + kinetic_change.abs()
            + total_loss;
        let tolerance =
            controls.energy_absolute_tolerance_j + controls.energy_relative_tolerance * scale;
        let mass_l = rho * volume;
        let mass_d = rho_d * volume;
        if v.into_iter().any(|x| !x.is_finite())
            || ![
                volume,
                pressure_work,
                wall_energy,
                gravity_energy,
                outlet_potential,
                kinetic0,
                kinetic1,
                kinetic_change,
                total_loss,
                residual,
                tolerance,
                mass_l,
                mass_d,
            ]
            .into_iter()
            .all(f64::is_finite)
            || (volume != 0.0 && (mass_l == 0.0 || (rho_d > 0.0 && mass_d == 0.0)))
        {
            return Err(MaterialInertiaError::Invalid(
                "complete-interval physical accounting is unrepresentable",
            ));
        }
        if residual.abs() > tolerance {
            return Err(MaterialInertiaError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        Ok(MaterialInertialCapillaryStep {
            initial_state: self.state,
            state: endpoint.state,
            elapsed_time_s: boundary.duration_s,
            transferred_volume_m3: volume,
            transferred_liquid_mass_kg: mass_l,
            displaced_mass_kg: mass_d,
            initial_kinetic_energy_j: kinetic0,
            kinetic_energy_j: kinetic1,
            kinetic_energy_change_j: kinetic_change,
            pressure_work_j: pressure_work,
            wall_surface_energy_change_j: wall_energy,
            gravitational_storage_change_j: gravity_energy,
            outlet_potential_energy_j: outlet_potential,
            liquid_viscous_dissipation_j: v[5],
            displaced_viscous_dissipation_j: v[6],
            contact_line_viscous_dissipation_j: v[7],
            pinning_dissipation_j: v[8],
            reservoir_dissipation_j: v[9],
            total_dissipation_j: total_loss,
            energy_balance_residual_j: residual,
            initial_apparent_angle_radians: initial_angle,
            final_apparent_angle_radians: final_angle,
            final_motion,
            stopping_events: stops,
            first_stopping_time_s: first_stop,
            segments,
            iterations,
            maximum_pressure_residual_pa: max_residual,
        })
    }
}

fn validate_boundary(
    boundary: InertialCapillaryBoundary,
    controls: InertialCapillaryIntervalControls,
) -> Result<(), MaterialInertiaError> {
    let c = controls.mechanics;
    if ![
        boundary.applied_pressure_pa,
        boundary.opposing_gravity_m_s2,
        boundary.duration_s,
    ]
    .into_iter()
    .all(f64::is_finite)
        || boundary.duration_s < 0.0
        || !(1..=128).contains(&c.maximum_iterations)
        || ![
            c.pressure_absolute_tolerance_pa,
            c.velocity_absolute_tolerance_m_s,
            c.time_absolute_tolerance_s,
            c.energy_absolute_tolerance_j,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
        || !c.energy_relative_tolerance.is_finite()
        || c.energy_relative_tolerance < 0.0
        || controls.maximum_stopping_events >= usize::MAX / 2
    {
        return Err(MaterialInertiaError::Invalid(
            "finite forward boundaries and bounded positive solve controls required",
        ));
    }
    Ok(())
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialInertiaError> {
    if cancelled() {
        Err(MaterialInertiaError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }
    fn value(self) -> f64 {
        self.total + self.correction
    }
}
fn accumulate_raw(sums: &mut [Sum; 10], step: &mechanics::InertialCapillaryStep) {
    let values = [
        step.transferred_volume_m3,
        step.pressure_work_j,
        step.wall_surface_energy_change_j,
        step.gravitational_storage_change_j,
        step.outlet_potential_energy_j,
        step.liquid_viscous_dissipation_j,
        step.displaced_viscous_dissipation_j,
        step.contact_line_viscous_dissipation_j,
        step.pinning_dissipation_j,
        step.reservoir_dissipation_j,
    ];
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value);
    }
}
fn accumulate_interval(sums: &mut [Sum; 10], step: &MaterialInertialCapillaryStep) {
    let values = [
        step.transferred_volume_m3,
        step.pressure_work_j,
        step.wall_surface_energy_change_j,
        step.gravitational_storage_change_j,
        step.outlet_potential_energy_j,
        step.liquid_viscous_dissipation_j,
        step.displaced_viscous_dissipation_j,
        step.contact_line_viscous_dissipation_j,
        step.pinning_dissipation_j,
        step.reservoir_dissipation_j,
    ];
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value);
    }
}
