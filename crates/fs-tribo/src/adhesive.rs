//! Adhesive junction shear on a caller-resolved real contact area.
//!
//! The capacity is `tau * A + mu * N_rep`, where `A` is the actual attached
//! area and `N_rep` is the nonnegative repulsive contribution to normal load.
//! It is NOT the signed net normal reaction: an adhesive junction can transmit
//! shear at zero or tensile net load. Setting the pressure law's coefficients
//! to zero gives the constant-shear-strength `F = tau A` model observed in
//! Carpick, Ogletree and Salmeron, Langmuir 12 (1996), DOI 10.1021/la9509007.
//!
//! Both shear strengths are explicit interface data, not bulk yield strength
//! and not inferred from work of adhesion. This is an uncoupled junction model:
//! the normal owner supplies area; shear-induced peeling, partial slip, static
//! reaction selection, ageing and wear require additional models. It does not
//! upgrade the caller's provenance. Evaluation is allocation-free and O(1).

use crate::{
    FrictionLaw, FrictionRegime, InputProvenance, InterfaceSystemRef, TangentialSlip, TriboError,
    checked_add, checked_mul, finite, nonnegative_finite, positive_finite,
};

/// Geometry and repulsive load of one already-resolved adhesive contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveContactState {
    area_m2: f64,
    repulsive_normal_force_n: f64,
}

impl AdhesiveContactState {
    /// Supply actual contact area [m2] and repulsive normal load [N]. A detached
    /// junction has both zero. A signed net force must never be supplied here.
    /// A normal-contact solver (for example JKR) owns these two quantities.
    pub fn new(area_m2: f64, repulsive_normal_force_n: f64) -> Result<Self, TriboError> {
        nonnegative_finite(area_m2, "adhesive contact area")?;
        nonnegative_finite(repulsive_normal_force_n, "adhesive repulsive load")?;
        if area_m2 == 0.0 && repulsive_normal_force_n != 0.0 {
            return Err(TriboError::InvalidInput {
                field: "nonzero repulsive load requires nonzero contact area",
            });
        }
        Ok(Self {
            area_m2,
            repulsive_normal_force_n,
        })
    }

    /// Real attached area, m2; not a nominal geometric overlap area.
    #[must_use]
    pub const fn area_m2(self) -> f64 {
        self.area_m2
    }

    /// Repulsive contribution, N, before subtracting the adhesive attraction.
    #[must_use]
    pub const fn repulsive_normal_force_n(self) -> f64 {
        self.repulsive_normal_force_n
    }
}

/// Constant junction shear plus the existing speed-dependent pressure law.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveFrictionLaw {
    pressure_law: FrictionLaw,
    static_shear_strength_pa: f64,
    kinetic_shear_strength_pa: f64,
}

impl AdhesiveFrictionLaw {
    /// All parameters must be measured or explicitly authored for the same
    /// interface. The caller must avoid double-counting adhesion already fitted
    /// into the pressure law. Zero shear strengths recover the pressure law.
    pub fn new(
        pressure_law: FrictionLaw,
        static_shear_strength_pa: f64,
        kinetic_shear_strength_pa: f64,
    ) -> Result<Self, TriboError> {
        pressure_law.validate()?;
        nonnegative_finite(static_shear_strength_pa, "static junction shear strength")?;
        nonnegative_finite(kinetic_shear_strength_pa, "kinetic junction shear strength")?;
        if kinetic_shear_strength_pa > static_shear_strength_pa {
            return Err(TriboError::InvalidInput {
                field: "kinetic junction shear strength exceeds static strength",
            });
        }
        Ok(Self {
            pressure_law,
            static_shear_strength_pa,
            kinetic_shear_strength_pa,
        })
    }

    /// The explicitly selected pressure contribution; no material-name dispatch.
    #[must_use]
    pub const fn pressure_law(&self) -> &FrictionLaw {
        &self.pressure_law
    }

    /// Set-valued static capacity [N]; the mechanical solver chooses the reaction.
    pub fn static_capacity_n(&self, contact: AdhesiveContactState) -> Result<f64, TriboError> {
        capacity(
            self.static_shear_strength_pa,
            self.pressure_law.static_coefficient()?,
            contact,
        )
    }

    /// Positive-speed sliding capacity [N]. Unlike an effective Coulomb
    /// coefficient, this remains well-defined at zero net normal load.
    pub fn sliding_capacity_n(
        &self,
        contact: AdhesiveContactState,
        speed_m_per_s: f64,
    ) -> Result<f64, TriboError> {
        capacity(
            self.kinetic_shear_strength_pa,
            self.pressure_law.kinetic_coefficient(speed_m_per_s)?,
            contact,
        )
    }

    /// Traction opposes checked tangential slip; lost mechanical power is
    /// returned once for routing through the existing thermal partition owner.
    /// At rest this returns capacity and zero force, not an invented reaction.
    pub fn evaluate<'a>(
        &self,
        interface: &'a InterfaceSystemRef,
        contact: AdhesiveContactState,
        slip: TangentialSlip,
    ) -> Result<AdhesiveFrictionResponse<'a>, TriboError> {
        interface.validate()?;
        let static_capacity_n = self.static_capacity_n(contact)?;
        let speed = slip.speed()?;
        let (regime, sliding_capacity_n, traction_n, dissipated_power_w) = if speed == 0.0 {
            (FrictionRegime::Sticking, None, [0.0; 3], 0.0)
        } else {
            let magnitude = self.sliding_capacity_n(contact, speed)?;
            // Normalize before multiplying: magnitude/speed can overflow even
            // when the physically bounded traction is representable.
            let traction = slip.velocity_mps().map(|v| -(v / speed) * magnitude);
            for x in traction {
                finite(x, "adhesive friction traction")?;
            }
            let power = checked_mul(magnitude, speed, "adhesive friction power")?;
            (FrictionRegime::Sliding, Some(magnitude), traction, power)
        };
        Ok(AdhesiveFrictionResponse {
            regime,
            static_capacity_n,
            sliding_capacity_n,
            traction_n,
            dissipated_power_w,
            provenance: &interface.provenance,
        })
    }

    /// Explicit-ODE scalar adapter. Below the caller's numerical speed scale,
    /// linearly join zero to the complete kinetic capacity at that scale.
    /// This is continuous and dissipative; it does not solve static friction.
    pub fn regularized_traction_1d(
        &self,
        contact: AdhesiveContactState,
        velocity_m_per_s: f64,
        regularization_m_per_s: f64,
    ) -> Result<f64, TriboError> {
        finite(velocity_m_per_s, "adhesive relative velocity")?;
        positive_finite(regularization_m_per_s, "adhesive regularization speed")?;
        // Validate the actual area/load capacity even when velocity is zero.
        self.static_capacity_n(contact)?;
        let speed = velocity_m_per_s.abs();
        if speed == 0.0 {
            return Ok(0.0);
        }
        let limit = self.sliding_capacity_n(contact, speed.max(regularization_m_per_s))?;
        let ratio = (speed / regularization_m_per_s).min(1.0);
        let force = -velocity_m_per_s.signum() * (limit * ratio);
        finite(force, "adhesive regularized traction")?;
        Ok(force)
    }
}

/// One immutable, allocation-free adhesive friction query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveFrictionResponse<'a> {
    /// Rest/sliding branch; a static capacity does not select a static reaction.
    pub regime: FrictionRegime,
    /// Total static force capacity, N, including the junction shear term.
    pub static_capacity_n: f64,
    /// Total positive-speed sliding capacity, N; absent at exact rest.
    pub sliding_capacity_n: Option<f64>,
    /// Force on the body whose relative tangent velocity was supplied, N.
    pub traction_n: [f64; 3],
    /// Nonnegative `-force dot relative_velocity`, W, for one heat allocation.
    pub dissipated_power_w: f64,
    /// Original caller provenance, borrowed without promotion or allocation.
    pub provenance: &'a InputProvenance,
}

fn capacity(strength: f64, mu: f64, contact: AdhesiveContactState) -> Result<f64, TriboError> {
    checked_add(
        checked_mul(strength, contact.area_m2, "junction shear capacity")?,
        checked_mul(
            mu,
            contact.repulsive_normal_force_n,
            "repulsive friction capacity",
        )?,
        "adhesive total friction capacity",
    )
}
