//! Source-defined, velocity-regularized kinetic friction on damaged bond area.
//!
//! The original cohesive material supplies `P = integral d(s) Kc <-gn(s)> ds`
//! along the actual straight jump path. With constant path slip speed `v`, the
//! mean resisting shear is `mu P v/max(|v|,v_reg)`. The residual therefore uses
//! the same damage, compressive contact and work-conjugate trace map as the bond.
//! This Huber regularization reproduces kinetic Coulomb above an explicit speed
//! and continuous viscous-like creep below it. It is not ideal static friction,
//! an elastic friction spring, a second damage state or a calibrated rough-
//! surface mixture. The numerical speed is never inferred from a material name.

use fs_exec::Cx;
use fs_material::cohesive::source::RESIDUAL_KINETIC_FRICTION_PROPERTY;
use fs_qty::{Dims, QuantitySpec};

use super::{CohesiveMotion, CohesiveSystem, DynamicsError, Error, finite, positive};
use crate::cohesive::CohesivePatch;
use std::sync::Arc;

/// Kinetic coefficient from the patch's original ordered interface source and
/// an explicit numerical regularization speed. Source receipts remain in that
/// same patch's complete `ResolvedBilinearCohesive` bundle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DamagedAreaFriction {
    kinetic_coefficient: f64,
    regularization_speed_m_s: f64,
}

impl DamagedAreaFriction {
    /// Source-defined residual coefficient; independent of bond strength/Gc.
    #[must_use]
    pub const fn kinetic_coefficient(self) -> f64 {
        self.kinetic_coefficient
    }
    /// Explicit velocity below which friction is linearly regularized.
    #[must_use]
    pub const fn regularization_speed_m_s(self) -> f64 {
        self.regularization_speed_m_s
    }

    fn from_patch(patch: &CohesivePatch<'_>, speed: f64) -> Result<Self, Error> {
        positive(speed)?;
        let source = patch.source.ok_or(Error::Invalid(
            "damaged-area friction requires the original sourced bond",
        ))?;
        let property = source
            .resolved()
            .property(RESIDUAL_KINETIC_FRICTION_PROPERTY)
            .ok_or(Error::Invalid(
                "original cohesive bundle lacks residual kinetic friction",
            ))?;
        if property.requirement().quantity() != QuantitySpec::dimensional(Dims::NONE)
            || !property.is_static_at_fixed_external_conditions()
        {
            return Err(Error::Invalid(
                "residual kinetic friction must be dimensionless and fixed at the source query",
            ));
        }
        let mu = finite(property.value_si())?;
        if mu < 0.0 {
            return Err(Error::Invalid(
                "residual kinetic friction must be nonnegative",
            ));
        }
        Ok(Self {
            kinetic_coefficient: mu,
            regularization_speed_m_s: speed,
        })
    }

    pub(super) fn active(self) -> bool {
        self.kinetic_coefficient != 0.0
    }

    /// Mean traction in Pa and exact work per original area in J/m2 for this
    /// declared straight path. Pressure is the material's actual path integral,
    /// never total normal traction, endpoint damage or initial compression.
    pub(super) fn response(
        self,
        pressure: f64,
        delta: [f64; 2],
        dt: f64,
    ) -> Result<([f64; 2], f64), Error> {
        let slip = finite(delta[0].hypot(delta[1]))?;
        let regularization = finite(self.regularization_speed_m_s * dt)?;
        positive(regularization)?;
        finite(pressure)?;
        if pressure < 0.0 {
            return Err(Error::Invalid("negative damaged compression measure"));
        }
        if slip == 0.0 || pressure == 0.0 {
            return Ok(([0.0; 2], 0.0));
        }
        let cap = finite(self.kinetic_coefficient * pressure)?;
        let magnitude = finite(cap * (slip / slip.max(regularization)))?;
        let traction = delta.map(|d| magnitude * (d / slip));
        let work = finite(magnitude * slip)?;
        Ok((traction, work))
    }
}

impl<'m, 'a> CohesiveMotion<'m, 'a> {
    /// Bind residual friction once at trajectory adoption, retaining the original
    /// displacement, damage, velocity and clock. One entry per original patch:
    /// `None` leaves it frictionless; `Some(v_reg)` selects the sourced kinetic
    /// property already retained in that cohesive bundle. Missing data refuse.
    /// All subsequent mechanical steps retain this choice; an explicit source
    /// temperature transfer refreshes the same selected property and no others.
    pub fn with_damaged_area_friction(
        self,
        speeds_m_s: &[Option<f64>],
        cx: &Cx<'_>,
    ) -> Result<Self, DynamicsError> {
        self.with_damaged_area_friction_controlled(speeds_m_s, &mut || cx.checkpoint().is_err())
    }

    /// Predicate counterpart; no partial configuration is returned on refusal.
    pub fn with_damaged_area_friction_controlled<C: FnMut() -> bool>(
        mut self,
        speeds_m_s: &[Option<f64>],
        cancelled: &mut C,
    ) -> Result<Self, DynamicsError> {
        poll(cancelled)?;
        if self.friction.is_some()
            || self.coulomb.is_some()
            || speeds_m_s.len() != self.state.points.len()
        {
            return Err(
                Error::Invalid("bind friction once with one choice per original patch").into(),
            );
        }
        let mut models = Vec::with_capacity(speeds_m_s.len());
        for (patch, speed) in self.state.system.patches.iter().zip(speeds_m_s) {
            poll(cancelled)?;
            models.push(
                speed
                    .map(|v| DamagedAreaFriction::from_patch(patch, v))
                    .transpose()?,
            );
        }
        let models = models.into();
        poll(cancelled)?;
        self.friction = Some(models);
        Ok(self)
    }

    /// Bound coefficients and regularization speeds, in original patch order.
    #[must_use]
    pub fn damaged_area_friction(&self) -> Option<&[Option<DamagedAreaFriction>]> {
        self.friction.as_deref()
    }
}

// The enclosing temperature transfer has already checked the SAME original
// card, selected claims, ordered surfaces, nonthermal query and geometry. Since
// friction is read from that complete bond bundle, no parallel source identity
// can drift, and changing mu at fixed jump supplies no friction work or heat.
pub(super) fn temperature_successors(
    original: Option<&[Option<DamagedAreaFriction>]>,
    target: &CohesiveSystem<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Arc<[Option<DamagedAreaFriction>]>>, Error> {
    let Some(original) = original else {
        return Ok(None);
    };
    let mut models = Vec::with_capacity(original.len());
    for (patch, old) in target.patches.iter().zip(original) {
        poll(cancelled)?;
        models.push(
            old.map(|m| DamagedAreaFriction::from_patch(patch, m.regularization_speed_m_s))
                .transpose()?,
        );
    }
    let models = models.into();
    poll(cancelled)?;
    Ok(Some(models))
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
