//! Set-valued junction shear on a frozen, already resolved normal patch.
use super::*;
use crate::{DryTangentialState, TangentialLaw, TangentialStep, TangentialStepInput};

/// One explicitly frozen real area and repulsive load. The normal solver owns
/// their evolution; constructing this value adds no tangential elastic storage.
/// Changing the patch requires a new value, not relabeling signed net force as
/// repulsion. Shear-induced area change and peeling are outside this model.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveTangentialContact {
    law: AdhesiveFrictionLaw,
    contact: AdhesiveContactState,
}

impl AdhesiveTangentialContact {
    /// Freeze the checked area/load and the original junction/pressure law.
    pub fn new(
        law: AdhesiveFrictionLaw,
        contact: AdhesiveContactState,
    ) -> Result<Self, TriboError> {
        law.static_capacity_n(contact)?;
        Ok(Self { law, contact })
    }
    /// Original complete constitutive law.
    #[must_use]
    pub const fn law(&self) -> &AdhesiveFrictionLaw {
        &self.law
    }
    /// Actual real area and repulsive force supplied by the normal owner.
    #[must_use]
    pub const fn contact(&self) -> AdhesiveContactState {
        self.contact
    }

    fn kinetic_capacity(&self, speed: f64) -> Result<f64, TriboError> {
        let mu = if speed == 0.0 {
            self.law.pressure_law.zero_speed_kinetic_coefficient()?
        } else {
            self.law.pressure_law.kinetic_coefficient(speed)?
        };
        capacity(self.law.kinetic_shear_strength_pa, mu, self.contact)
    }

    // Solve the real offset-plus-pressure load line directly. A clamped linear
    // pressure coefficient can vanish while junction shear remains nonzero.
    fn affine_capacity(
        &self,
        distance: f64,
        compliance: f64,
        dt: f64,
    ) -> Result<Option<f64>, TriboError> {
        let junction = checked_mul(
            self.law.kinetic_shear_strength_pa,
            self.contact.area_m2,
            "junction shear capacity",
        )?;
        let normal = self.contact.repulsive_normal_force_n;
        if normal == 0.0 {
            return Ok(Some(junction));
        }
        if junction == 0.0 {
            return self
                .law
                .pressure_law
                .affine_sliding_traction_n(normal, distance, compliance, dt);
        }
        match self.law.pressure_law {
            FrictionLaw::Coulomb { kinetic_mu, .. } => Ok(Some(capacity(
                self.law.kinetic_shear_strength_pa,
                kinetic_mu,
                self.contact,
            )?)),
            FrictionLaw::VelocityDependent {
                mu_zero,
                slope_per_speed,
                ..
            } => {
                let slip_at_junction = distance - compliance * junction;
                if slip_at_junction > 0.0
                    && self
                        .law
                        .pressure_law
                        .kinetic_coefficient(slip_at_junction / dt)?
                        == 0.0
                {
                    return Ok(Some(junction));
                }
                let speed = distance / dt;
                finite(speed, "slip_speed_mps")?;
                let released = junction + normal * (mu_zero + slope_per_speed * speed);
                finite(released, "affine released traction")?;
                let gain = normal * slope_per_speed * (compliance / dt);
                let denominator = 1.0 + gain;
                if !gain.is_finite() || denominator <= 64.0 * f64::EPSILON * (1.0 + gain.abs()) {
                    return Ok(None);
                }
                let force = released / denominator;
                Ok((force.is_finite() && force >= junction).then_some(force))
            }
            FrictionLaw::Stribeck { .. } | FrictionLaw::Hyperbolic { .. } => Ok(None),
        }
    }
}

impl TangentialLaw for AdhesiveTangentialContact {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        "adhesive-junction-held-effort"
    }
    fn initial_state(&self) -> Self::State {
        self.law.pressure_law.initial_state()
    }
    fn stored_energy_j(&self, _s: &Self::State) -> f64 {
        0.0
    }
    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        // A redundant load is checked, never silently ignored or interpreted as
        // the signed adhesive net reaction. Zero repulsion can still carry shear.
        if x.normal_n != self.contact.repulsive_normal_force_n {
            return Err(TriboError::InvalidInput {
                field: "held adhesive input must retain the frozen repulsive normal load",
            });
        }
        crate::tangential::advance_with_capacity(
            s,
            x,
            self.contact.area_m2 == 0.0,
            self.law.static_capacity_n(self.contact)?,
            |speed| self.kinetic_capacity(speed),
            |distance, compliance, dt| self.affine_capacity(distance, compliance, dt),
        )
    }
}
