//! Original interface sources for a frozen low-speed contact-line law.
//!
//! Line friction is independent of bulk liquid viscosity. Advancing/receding
//! thresholds are independent of the equilibrium surface-energy potential.
//! This binding preserves those distinctions and the complete ordered source
//! bundle; geometry, motion, rate limits and dissipation belong to the consumer.

use super::adhesion::{EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin};
use super::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, ResolvedLiquidSurfaceState, ResolvedLiquidWettingState,
    SURFACE_TENSION_DIMS, require_frozen_coefficient,
};
use crate::liquid::DYNAMIC_VISCOSITY_DIMS;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedScalarProperty,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use core::f64::consts::PI;
use fs_matdb::InterfaceSystemCard;
use fs_qty::{Dims, Time};

pub const CONTACT_LINE_FRICTION_PROPERTY: &str = "contact-line-friction";
pub const ADVANCING_CONTACT_ANGLE_PROPERTY: &str = "advancing-contact-angle";
pub const RECEDING_CONTACT_ANGLE_PROPERTY: &str = "receding-contact-angle";

/// Resolve the original equilibrium and contact-line claims at one liquid point.
///
/// The result retains every original receipt. Resolution alone does not admit
/// the law for frozen use: [`FrozenContactLineWetting::new`] checks quantities,
/// coefficient support and angle ordering before an actual consumer uses it.
pub fn resolve_contact_line_wetting(
    card: &InterfaceSystemCard,
    liquid: &ResolvedLiquidSurfaceState,
    origin: EquilibriumAngleOrigin,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidWettingState, MaterialStatePointError> {
    let (name, dims) = match origin {
        EquilibriumAngleOrigin::SourcedAngle => (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        EquilibriumAngleOrigin::YoungDupre => {
            (EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, SURFACE_TENSION_DIMS)
        }
    };
    let requirements = [
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)?,
        ScalarPropertyRequirement::try_new(
            CONTACT_LINE_FRICTION_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            ADVANCING_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            RECEDING_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let resolved =
        resolve_interface_state_point(card, &liquid.material().query(), &requirements, selection)?;
    let wetting = match origin {
        EquilibriumAngleOrigin::SourcedAngle => {
            ResolvedLiquidWettingState::try_from_resolved(card, resolved)
        }
        EquilibriumAngleOrigin::YoungDupre => {
            ResolvedLiquidWettingState::try_from_adhesion_work(card, resolved, liquid)
        }
    }?;
    wetting.require_liquid_state(liquid)?;
    Ok(wetting)
}

/// Checked, allocation-free view of one original contact-line source bundle.
///
/// Coefficients may be frozen at their admitted temperature and fixed external
/// conditions. Evolving speed, rate or numeric history support is not silently
/// converted into this constant-coefficient law. Named preparation/history in
/// the interface context remains part of the original source identity.
#[derive(Clone, Copy, Debug)]
pub struct FrozenContactLineWetting<'a> {
    wetting: &'a ResolvedLiquidWettingState,
    friction_pa_s: f64,
    advancing_angle_radians: f64,
    receding_angle_radians: f64,
}

impl<'a> FrozenContactLineWetting<'a> {
    pub fn new(
        liquid: &ResolvedLiquidSurfaceState,
        wetting: &'a ResolvedLiquidWettingState,
    ) -> Result<Self, MaterialStatePointError> {
        wetting.require_liquid_state(liquid)?;
        let equilibrium_name = match wetting.angle_origin() {
            EquilibriumAngleOrigin::SourcedAngle => EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
            EquilibriumAngleOrigin::YoungDupre => EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
        };
        let equilibrium = wetting.resolved().property(equilibrium_name).ok_or(
            MaterialStatePointError::InvalidDerived {
                quantity: equilibrium_name,
            },
        )?;
        require_no_numeric_history(equilibrium)?;
        let mut values = [0.0; 3];
        for (slot, (name, dims)) in [
            (CONTACT_LINE_FRICTION_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
            (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
            (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        ]
        .into_iter()
        .enumerate()
        {
            let property = wetting
                .resolved()
                .property(name)
                .ok_or(MaterialStatePointError::InvalidDerived { quantity: name })?;
            require_frozen_coefficient(property, name, dims)?;
            require_no_numeric_history(property)?;
            let value = property.value_si();
            if !value.is_finite() || value <= 0.0 {
                return Err(MaterialStatePointError::OutsideConsumerDomain {
                    property: name.to_owned(),
                    value,
                    admissibility: ScalarAdmissibility::StrictlyPositive,
                });
            }
            values[slot] = value;
        }
        let equilibrium = wetting.contact_angle_radians();
        if !(values[2] <= equilibrium && equilibrium <= values[1] && values[1] < PI) {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "ordered receding/equilibrium/advancing contact angles in (0, pi)",
            });
        }
        Ok(Self {
            wetting,
            friction_pa_s: values[0],
            advancing_angle_radians: values[1],
            receding_angle_radians: values[2],
        })
    }

    /// Original ordered interface, all selected claims and their complete query.
    #[must_use]
    pub const fn wetting(&self) -> &'a ResolvedLiquidWettingState {
        self.wetting
    }
    /// Positive contact-line friction [Pa s], independent of bulk viscosity.
    #[must_use]
    pub const fn friction_pa_s(&self) -> f64 {
        self.friction_pa_s
    }
    #[must_use]
    pub const fn equilibrium_angle_radians(&self) -> f64 {
        self.wetting.contact_angle_radians()
    }
    #[must_use]
    pub const fn advancing_angle_radians(&self) -> f64 {
        self.advancing_angle_radians
    }
    #[must_use]
    pub const fn receding_angle_radians(&self) -> f64 {
        self.receding_angle_radians
    }
}

fn require_no_numeric_history(
    property: &ResolvedScalarProperty,
) -> Result<(), MaterialStatePointError> {
    for (axis, _) in &property.answer().receipt.query_point {
        if property
            .answer()
            .evidence
            .model
            .validity
            .bound(axis)
            .is_none()
        {
            continue;
        }
        let temporal = property
            .answer()
            .receipt
            .axis_quantities
            .get(axis)
            .is_some_and(|quantity| quantity.dims() == Time::DIMS);
        if temporal
            || [
                "time",
                "elapsed_time",
                "age",
                "history",
                "wetting_history",
                "contact_line_history",
                "cycle_count",
                "cycles",
                "contact_line_displacement",
            ]
            .contains(&axis.as_str())
        {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: property.requirement().name().to_owned(),
                reason: "contact-line history-bounded coefficients need an evolving source owner",
            });
        }
    }
    Ok(())
}
