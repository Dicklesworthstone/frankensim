//! Optional sourced contact-line friction and hysteresis for the same tube owner.
//!
//! This low-speed, long-tube closure retains equilibrium wall free energy.
//! Moving apparent angles supply separate line-friction and pinning losses;
//! changing meniscus end-cap area/volume and liquid inertia are not resolved.

use super::{MaterialCapillaryError, MaterialCapillaryTube, MaterialMeniscusInput};
use crate::capillary::{
    ContactLineCircularMeniscusInput, InclinedTwoFluidCircularMeniscusInput,
    TwoFluidCircularMeniscusInput,
};
use fs_material::liquid::surface::contact_line::FrozenContactLineWetting;

/// Explicit operating domain of the sourced linear contact-line closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineTubeLimits {
    /// Maximum admitted contact-line speed [m/s]; never inferred from viscosity.
    pub maximum_speed_m_s: f64,
    /// Maximum apparent-angle departure from the original Young angle [rad].
    pub maximum_angle_departure_radians: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct BoundContactLine<'a> {
    pub source: FrozenContactLineWetting<'a>,
    pub limits: ContactLineTubeLimits,
}

impl<'a> MaterialCapillaryTube<'a> {
    /// Enable friction and advancing/receding thresholds from this exact
    /// original wetting bundle. Missing coefficients remain a source refusal.
    /// The same optional law is consumed by open flow and sealed gas coupling.
    pub fn with_contact_line(
        mut self,
        limits: ContactLineTubeLimits,
    ) -> Result<Self, MaterialCapillaryError> {
        if !limits.maximum_speed_m_s.is_finite()
            || limits.maximum_speed_m_s <= 0.0
            || !limits.maximum_angle_departure_radians.is_finite()
            || limits.maximum_angle_departure_radians <= 0.0
            || limits.maximum_angle_departure_radians > core::f64::consts::PI
        {
            return Err(MaterialCapillaryError::Unsupported(
                "finite positive contact-line speed/angle bounds are required",
            ));
        }
        let source = FrozenContactLineWetting::new(self.liquid, self.wetting)?;
        self.contact_line = Some(BoundContactLine { source, limits });
        Ok(self)
    }

    /// Original checked contact-line source, absent for equilibrium-only flow.
    #[must_use]
    pub fn contact_line_source(&self) -> Option<&FrozenContactLineWetting<'a>> {
        self.contact_line.as_ref().map(|line| &line.source)
    }

    #[must_use]
    pub fn contact_line_limits(&self) -> Option<ContactLineTubeLimits> {
        self.contact_line.map(|line| line.limits)
    }

    pub(super) fn contact_line_input(
        &self,
        input: MaterialMeniscusInput,
        displaced_viscosity: f64,
        displaced_density: f64,
        gravity: f64,
    ) -> Option<ContactLineCircularMeniscusInput> {
        self.contact_line
            .map(|line| ContactLineCircularMeniscusInput {
                meniscus: InclinedTwoFluidCircularMeniscusInput {
                    meniscus: TwoFluidCircularMeniscusInput {
                        meniscus: self.primitive_input(input),
                        displaced_dynamic_viscosity_pa_s: displaced_viscosity,
                    },
                    density_kg_m3: self.liquid.density_kg_m3(),
                    displaced_density_kg_m3: displaced_density,
                    opposing_gravity_m_s2: gravity,
                },
                contact_line_friction_pa_s: line.source.friction_pa_s(),
                advancing_angle_radians: line.source.advancing_angle_radians(),
                receding_angle_radians: line.source.receding_angle_radians(),
                max_speed_m_s: line.limits.maximum_speed_m_s,
                max_angle_departure_radians: line.limits.maximum_angle_departure_radians,
            })
    }
}
