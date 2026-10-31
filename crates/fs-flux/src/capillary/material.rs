//! Source-bound, quasisteady capillary transport through a circular tube.
//!
//! The original liquid and ordered substrate/liquid interface own density,
//! viscosity, surface tension and equilibrium wetting. This adapter supplies
//! no material defaults and evolves only the meniscus position. Coefficients
//! remain frozen at the complete admitted source point; neither heat nor a
//! dynamic contact-line law is inferred from an equilibrium angle.

use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::liquid::surface::{ResolvedLiquidSurfaceState, ResolvedLiquidWettingState};
use fs_material::state_point::MaterialStatePointError;

use super::{
    CircularCapillaryError, CircularMeniscusInput, CircularMeniscusStep,
    InclinedCircularMeniscusInput, InclinedCircularMeniscusStep,
    InclinedTwoFluidCircularMeniscusInput, InclinedTwoFluidCircularMeniscusStep,
    TwoFluidCircularMeniscusInput, TwoFluidCircularMeniscusStep,
    step_newtonian_inclined_circular_meniscus, step_newtonian_inclined_two_fluid_circular_meniscus,
};

/// Geometry, current meniscus position and prescribed loading, all in SI.
/// Material coefficients come exclusively from the retained source owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialMeniscusInput {
    /// Radius of the rigid, circular, no-slip tube [m].
    pub radius_m: f64,
    /// Available tube length [m].
    pub tube_length_m: f64,
    /// Current liquid-column length [m], between zero and the tube length.
    pub wet_length_m: f64,
    /// Inlet liquid pressure minus prescribed downstream pressure [Pa].
    pub applied_pressure_pa: f64,
    /// Requested nonnegative advance [s]; an end event may consume less.
    pub duration_s: f64,
}

/// Refusal before publishing an accepted meniscus or physical transfers.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialCapillaryError {
    /// The original wetting and liquid source assignments do not agree.
    Source(MaterialStatePointError),
    /// The displaced-fluid owner refused its transport evaluation.
    Transport(MaterialError),
    /// An evolution obligation or required source is outside this reduced law.
    Unsupported(&'static str),
    /// The existing capillary operator refused its numerical inputs or result.
    Flow(CircularCapillaryError),
    /// A signed mass transfer cannot be represented from the accepted volume.
    NonFiniteMass,
}

impl core::fmt::Display for MaterialCapillaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "material capillary source: {error}"),
            Self::Transport(error) => write!(f, "displaced capillary fluid: {error}"),
            Self::Unsupported(reason) => write!(f, "material capillary: {reason}"),
            Self::Flow(error) => write!(f, "material capillary flow: {error:?}"),
            Self::NonFiniteMass => {
                write!(f, "material capillary mass transfer is not representable")
            }
        }
    }
}

impl std::error::Error for MaterialCapillaryError {}

impl From<MaterialStatePointError> for MaterialCapillaryError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}

impl From<CircularCapillaryError> for MaterialCapillaryError {
    fn from(error: CircularCapillaryError) -> Self {
        Self::Flow(error)
    }
}

/// A frozen liquid/interface pair feeding the existing circular-tube laws.
///
/// The references preserve every original claim and query condition. Binding
/// checks the liquid identity and complete point carried by the ordered wetting
/// card, including an adhesion-derived angle's original surface-tension source.
/// It does not equate a named interface environment with an arbitrary gas EOS.
#[derive(Clone, Copy, Debug)]
pub struct MaterialCapillaryTube<'a> {
    liquid: &'a ResolvedLiquidSurfaceState,
    wetting: &'a ResolvedLiquidWettingState,
    surface_tension_n_m: f64,
    contact_angle_cosine: f64,
}

impl<'a> MaterialCapillaryTube<'a> {
    /// Bind one actively sourced liquid surface to its equilibrium wall law.
    ///
    /// # Errors
    /// Refuses mismatched source assignments, neglected capillarity and liquid
    /// bundles requiring shear-rate/temperature evolution. Resolving a source
    /// once does not grant this Newtonian law an evolving rheology capability.
    pub fn new(
        liquid: &'a ResolvedLiquidSurfaceState,
        wetting: &'a ResolvedLiquidWettingState,
    ) -> Result<Self, MaterialCapillaryError> {
        wetting.require_liquid_state(liquid)?;
        if liquid.rheology_axes().is_some() {
            return Err(MaterialCapillaryError::Unsupported(
                "quasisteady tube transport requires a frozen Newtonian liquid",
            ));
        }
        let surface_tension_n_m =
            liquid
                .surface_tension_n_m()
                .ok_or(MaterialCapillaryError::Unsupported(
                    "material capillary transport requires an active sourced surface tension",
                ))?;
        Ok(Self {
            liquid,
            wetting,
            surface_tension_n_m,
            contact_angle_cosine: wetting.contact_angle_radians().cos(),
        })
    }

    /// Exact original liquid state, including its viscosity basis and receipts.
    #[must_use]
    pub const fn liquid(&self) -> &'a ResolvedLiquidSurfaceState {
        self.liquid
    }

    /// Exact original ordered substrate/liquid interface and ambient context.
    #[must_use]
    pub const fn wetting(&self) -> &'a ResolvedLiquidWettingState {
        self.wetting
    }

    /// Horizontal Washburn transport with negligible displaced-gas resistance.
    /// This retains the primitive's exact endpoint events and wall-energy work.
    pub fn step_horizontal(
        &self,
        input: MaterialMeniscusInput,
    ) -> Result<MaterialMeniscusStep<'a>, MaterialCapillaryError> {
        self.step_inclined(input, 0.0)
    }

    /// Rise or drainage with signed gravity opposing advance [m/s²].
    ///
    /// The inlet is at the reservoir free-surface elevation. Exterior gas
    /// pressure is uniform and prescribed, and its inertia/resistance are
    /// neglected. Density, hydrostatic feedback and signed transferred mass
    /// all use the same source state. Source coefficients do not evolve as the
    /// meniscus moves; dynamic wetting, entrance losses and inertia are excluded.
    pub fn step_inclined(
        &self,
        input: MaterialMeniscusInput,
        opposing_gravity_m_s2: f64,
    ) -> Result<MaterialMeniscusStep<'a>, MaterialCapillaryError> {
        let discrete = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
            meniscus: self.primitive_input(input),
            density_kg_m3: self.liquid.density_kg_m3(),
            opposing_gravity_m_s2,
        })?;
        let transferred_mass_kg = signed_mass(
            self.liquid.density_kg_m3(),
            discrete.meniscus.transferred_volume_m3,
        )?;
        Ok(MaterialMeniscusStep {
            tube: *self,
            input,
            opposing_gravity_m_s2,
            discrete,
            transferred_mass_kg,
        })
    }

    /// Horizontal displacement with both fluids' finite Poiseuille resistance.
    ///
    /// The displaced phase is vented to its prescribed outlet pressure and is
    /// homogeneous, Newtonian and frozen at the supplied transport state. A gas
    /// owner supplies its viscosity and density, not compressible gas storage:
    /// this is a low-pressure-variation approximation, not a sealed-bubble model.
    /// A second liquid is admissible only when the retained tension and wetting
    /// sources describe that actual fluid pair. The caller owns that scene
    /// assignment; the generic transport view cannot establish it.
    ///
    /// The original displaced-fluid owner and its single checked evaluation
    /// remain in the result. A raw transport fixture gains no sourced authority.
    /// No heat capacity or conductivity defaults are supplied by this adapter.
    pub fn step_two_fluid_horizontal<'d, F: FluidTransport + ?Sized>(
        &self,
        input: MaterialMeniscusInput,
        displaced_fluid: &'d F,
    ) -> Result<MaterialTwoFluidMeniscusStep<'a, 'd, F>, MaterialCapillaryError> {
        self.step_two_fluid_inclined(input, displaced_fluid, 0.0)
    }

    /// Two-fluid displacement with both columns' hydrostatic and viscous loads.
    ///
    /// Both fluids remain homogeneous and incompressible at their retained
    /// source states. Applied pressure is inlet minus outlet pressure, before
    /// either hydrostatic column is subtracted. The result separates the tube's
    /// gravitational storage change from potential energy transported through
    /// the outlet; equal densities cancel only the storage change.
    pub fn step_two_fluid_inclined<'d, F: FluidTransport + ?Sized>(
        &self,
        input: MaterialMeniscusInput,
        displaced_fluid: &'d F,
        opposing_gravity_m_s2: f64,
    ) -> Result<MaterialTwoFluidMeniscusStep<'a, 'd, F>, MaterialCapillaryError> {
        let displaced_transport = displaced_fluid
            .transport_state()
            .map_err(MaterialCapillaryError::Transport)?;
        let discrete = step_newtonian_inclined_two_fluid_circular_meniscus(
            InclinedTwoFluidCircularMeniscusInput {
                meniscus: TwoFluidCircularMeniscusInput {
                    meniscus: self.primitive_input(input),
                    displaced_dynamic_viscosity_pa_s: displaced_transport.dynamic_viscosity_pa_s(),
                },
                density_kg_m3: self.liquid.density_kg_m3(),
                displaced_density_kg_m3: displaced_transport.density_kg_m3(),
                opposing_gravity_m_s2,
            },
        )?;
        let volume = discrete.meniscus.meniscus.transferred_volume_m3;
        let transferred_liquid_mass_kg = signed_mass(self.liquid.density_kg_m3(), volume)?;
        let displaced_mass_kg = signed_mass(displaced_transport.density_kg_m3(), volume)?;
        Ok(MaterialTwoFluidMeniscusStep {
            tube: *self,
            input,
            displaced_fluid,
            displaced_transport,
            opposing_gravity_m_s2,
            discrete,
            transferred_liquid_mass_kg,
            displaced_mass_kg,
        })
    }

    fn primitive_input(&self, input: MaterialMeniscusInput) -> CircularMeniscusInput {
        CircularMeniscusInput {
            dynamic_viscosity_pa_s: self.liquid.dynamic_viscosity_pa_s(),
            radius_m: input.radius_m,
            tube_length_m: input.tube_length_m,
            wet_length_m: input.wet_length_m,
            surface_tension_n_m: self.surface_tension_n_m,
            contact_angle_cosine: self.contact_angle_cosine,
            applied_pressure_pa: input.applied_pressure_pa,
            duration_s: input.duration_s,
        }
    }
}

/// Accepted horizontal/inclined transport with its actual material owners.
#[derive(Clone, Copy, Debug)]
pub struct MaterialMeniscusStep<'a> {
    tube: MaterialCapillaryTube<'a>,
    input: MaterialMeniscusInput,
    opposing_gravity_m_s2: f64,
    discrete: InclinedCircularMeniscusStep,
    transferred_mass_kg: f64,
}

impl<'a> MaterialMeniscusStep<'a> {
    /// Original liquid and wetting states, unchanged by the advance.
    #[must_use]
    pub const fn tube(&self) -> &MaterialCapillaryTube<'a> {
        &self.tube
    }

    /// Geometry, starting position and prescribed pressure/time actually used.
    #[must_use]
    pub const fn input(&self) -> MaterialMeniscusInput {
        self.input
    }

    /// Signed acceleration opposing advance [m/s²]; zero for horizontal flow.
    #[must_use]
    pub const fn opposing_gravity_m_s2(&self) -> f64 {
        self.opposing_gravity_m_s2
    }

    /// Accepted position, endpoint event, signed volume and wall/pressure work.
    #[must_use]
    pub const fn meniscus(&self) -> &CircularMeniscusStep {
        &self.discrete.meniscus
    }

    /// Complete physical result, including hydrostatic pressure and gravity work.
    #[must_use]
    pub const fn discrete(&self) -> &InclinedCircularMeniscusStep {
        &self.discrete
    }

    /// Signed liquid mass entering the tube [kg], using the original density.
    #[must_use]
    pub const fn transferred_mass_kg(&self) -> f64 {
        self.transferred_mass_kg
    }
}

/// Accepted two-fluid transport retaining both original constitutive owners.
pub struct MaterialTwoFluidMeniscusStep<'a, 'd, F: FluidTransport + ?Sized> {
    tube: MaterialCapillaryTube<'a>,
    input: MaterialMeniscusInput,
    displaced_fluid: &'d F,
    displaced_transport: FluidTransportState,
    opposing_gravity_m_s2: f64,
    discrete: InclinedTwoFluidCircularMeniscusStep,
    transferred_liquid_mass_kg: f64,
    displaced_mass_kg: f64,
}

impl<'a, 'd, F: FluidTransport + ?Sized> MaterialTwoFluidMeniscusStep<'a, 'd, F> {
    /// Original liquid and equilibrium-wetting sources.
    #[must_use]
    pub const fn tube(&self) -> &MaterialCapillaryTube<'a> {
        &self.tube
    }

    /// Geometry, starting position and prescribed pressure/time actually used.
    #[must_use]
    pub const fn input(&self) -> MaterialMeniscusInput {
        self.input
    }

    /// Exact displaced-fluid owner, including any source and query context.
    #[must_use]
    pub const fn displaced_fluid(&self) -> &'d F {
        self.displaced_fluid
    }

    /// The one checked, frozen evaluation consumed by the accepted step.
    #[must_use]
    pub const fn displaced_transport(&self) -> &FluidTransportState {
        &self.displaced_transport
    }

    /// Endpoint, pressure/surface work and separately integrated fluid losses.
    #[must_use]
    pub const fn discrete(&self) -> &TwoFluidCircularMeniscusStep {
        &self.discrete.meniscus
    }

    /// Full inclined result, including stored and advected gravitational energy.
    #[must_use]
    pub const fn inclined_discrete(&self) -> &InclinedTwoFluidCircularMeniscusStep {
        &self.discrete
    }

    /// Signed gravity opposing advance; zero for the horizontal entry point.
    #[must_use]
    pub const fn opposing_gravity_m_s2(&self) -> f64 {
        self.opposing_gravity_m_s2
    }

    /// Signed liquid mass entering the tube [kg].
    #[must_use]
    pub const fn transferred_liquid_mass_kg(&self) -> f64 {
        self.transferred_liquid_mass_kg
    }

    /// Signed displaced-fluid mass leaving the outlet [kg].
    /// Both fluids exchange the same volume, not generally the same mass.
    #[must_use]
    pub const fn displaced_mass_kg(&self) -> f64 {
        self.displaced_mass_kg
    }
}

fn signed_mass(density: f64, volume: f64) -> Result<f64, MaterialCapillaryError> {
    let mass = density * volume;
    if !mass.is_finite() || (volume != 0.0 && mass == 0.0) {
        return Err(MaterialCapillaryError::NonFiniteMass);
    }
    Ok(mass)
}
