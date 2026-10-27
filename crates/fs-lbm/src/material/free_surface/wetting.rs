//! Source-pinned hysteresis drives the existing planar capillary boundary.
//!
//! This explicit geometric closure estimates contact-line speed from successive
//! accepted wall-adjacent VOF fields. It does not identify wall-fluid velocity
//! with a no-slip contact-line velocity. The estimate is first order in space
//! at the wall and backward in time; its angle drives the NEXT native step.
//! In the deadband, the measured angle is projected onto the sourced interval.
//! That is a bounded hysteresis closure, not a proof of exact pinning, calibrated
//! dynamic wetting, or closed surface/mechanical/thermal energy conservation.

use std::collections::BTreeMap;

use super::*;
use crate::d3q19::freesurface3::wetting::WallFill3;
use crate::freesurface::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, ContactLineRegime2, DYNAMIC_WETTING_SPEED_AXIS,
    DYNAMIC_WETTING_SPEED_DIMS, RECEDING_CONTACT_ANGLE_PROPERTY,
};
use fs_matdb::{InterfaceSystemCard, MatDbError};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::{Dims, QuantitySpec};

/// Explicit model/numerical inputs; the deadband has no invented source receipt.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceWettingControls {
    /// Symmetric physical speed deadband [m/s], finite and nonnegative.
    pub pinning_speed_m_s: f64,
    /// Maximum active wall-adjacent cells, including gas history.
    pub maximum_wall_cells: usize,
}

/// One accepted contact estimate and the two ORIGINAL source endpoint claims.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceWettingSample {
    /// Positive is advancing liquid along its outward in-wall normal [m/s].
    pub signed_speed_m_s: f64,
    /// Direction/deadband selection; `Pinned` does not certify zero line motion.
    pub regime: ContactLineRegime2,
    /// Geometric angle through liquid from the accepted first-cell fill [rad].
    pub geometric_angle_radians: f64,
    /// Angle used by the next native curvature reconstruction [rad].
    pub applied_angle_radians: f64,
    /// Both pinned endpoint claims, ordered surfaces, source point and receipts.
    pub material: ResolvedInterfaceStatePoint,
}

/// A source, geometry, numerical or cancellation refusal; no state is published.
#[derive(Debug)]
pub enum FreeSurfaceWettingError {
    Material(MaterialStatePointError),
    Query(MatDbError),
    Flow(MaterialLbmError),
    Invalid(&'static str),
}

impl core::fmt::Display for FreeSurfaceWettingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "free-surface wetting material: {e}"),
            Self::Query(e) => write!(f, "free-surface wetting query: {e}"),
            Self::Flow(e) => write!(f, "free-surface wetting: {e}"),
            Self::Invalid(message) => write!(f, "free-surface wetting: {message}"),
        }
    }
}
impl std::error::Error for FreeSurfaceWettingError {}
impl From<MaterialStatePointError> for FreeSurfaceWettingError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<MatDbError> for FreeSurfaceWettingError {
    fn from(e: MatDbError) -> Self {
        Self::Query(e)
    }
}
impl From<crate::d3q19::freesurface3::FreeSurfaceError3> for FreeSurfaceWettingError {
    fn from(e: crate::d3q19::freesurface3::FreeSurfaceError3) -> Self {
        Self::Flow(e.into())
    }
}

type Error = FreeSurfaceWettingError;

/// One frozen isothermal liquid and its evolving planar contact-line history.
///
/// Bulk density, viscosity and sigma remain the original liquid's coefficients.
/// Surface A is the assigned stationary substrate, surface B must be this exact
/// liquid phase/process/revision. The original full liquid query is retained;
/// only the explicit physical contact-line speed coordinate evolves. Named gas
/// and surface/history context belongs to the interface card, not a guessed EOS.
pub struct MaterialFreeSurfaceWetting3<'s, 'm> {
    flow: MaterialFreeSurface3<'s>,
    card: &'m InterfaceSystemCard,
    source: ResolvedInterfaceStatePoint,
    wall: PlanarWall3,
    controls: FreeSurfaceWettingControls,
    geometry: BTreeMap<[usize; 3], WallFill3>,
    samples: BTreeMap<[usize; 3], FreeSurfaceWettingSample>,
}

impl<'s> MaterialFreeSurface3<'s> {
    /// Bind a rate-dependent advancing/receding interface to one planar wall.
    /// Initial history is stationary: both source claims must support speed zero.
    /// The caller assigns substrate and named environment to the scene. Existing
    /// static wetting and the legacy Wetting ghost cannot be combined with this
    /// law. Other domain walls and internal solid tiles remain neutral.
    pub fn with_dynamic_planar_wetting<'m>(
        mut self,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        selection: MaterialPropertySelection,
        controls: FreeSurfaceWettingControls,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWetting3<'s, 'm>, Error> {
        checkpoint(gate)?;
        if !controls.pinning_speed_m_s.is_finite()
            || controls.pinning_speed_m_s < 0.0
            || controls.maximum_wall_cells == 0
            || controls.maximum_wall_cells > 262_144
        {
            return Err(Error::Invalid(
                "invalid speed deadband or bounded wall-cell budget",
            ));
        }
        if self.scaling.planar_wetting.is_some() || !self.lattice.admits_dynamic_planar_wetting() {
            return Err(Error::Invalid(
                "dynamic planar wetting requires otherwise neutral walls",
            ));
        }
        let liquid = self.scaling.material.material();
        if &card.surface_b().material != liquid.material()
            || liquid
                .query_point()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
        {
            return Err(Error::Invalid(
                "wetting requires this liquid identity and a separate evolving speed axis",
            ));
        }
        let point = liquid.query().with_quantity(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
            0.0,
        )?;
        let requirements = [
            RECEDING_CONTACT_ANGLE_PROPERTY,
            ADVANCING_CONTACT_ANGLE_PROPERTY,
        ]
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: 0.0,
                    upper: core::f64::consts::PI,
                },
            )
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
        let source = resolve_interface_state_point(card, &point, &requirements, selection)?;
        require_owned_axes(&source)?;
        bracket(&source)?;
        let geometry = self
            .lattice
            .planar_wall_fill(wall, controls.maximum_wall_cells, gate)?;
        let mut candidate = MaterialFreeSurfaceWetting3 {
            flow: self,
            card,
            source,
            wall,
            controls,
            geometry,
            samples: BTreeMap::new(),
        };
        let samples = candidate.resolve_samples(&candidate.geometry, None, gate)?;
        candidate.flow.lattice.set_planar_contact_field(
            wall,
            samples
                .iter()
                .map(|(&p, s)| (p, s.applied_angle_radians))
                .collect(),
            gate,
        )?;
        checkpoint(gate)?;
        candidate.samples = samples;
        Ok(candidate)
    }
}

impl MaterialFreeSurfaceWetting3<'_, '_> {
    /// Read-only actual flow, native populations, mass ledger and physical units.
    #[must_use]
    pub const fn flow(&self) -> &MaterialFreeSurface3<'_> {
        &self.flow
    }

    /// Original two-claim binding at stationary initialization.
    #[must_use]
    pub const fn source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }

    /// Accepted endpoint queries and next-step angles by wall-cell coordinate.
    #[must_use]
    pub const fn samples(&self) -> &BTreeMap<[usize; 3], FreeSurfaceWettingSample> {
        &self.samples
    }

    /// Advance native momentum/mass and atomically refresh the contact history.
    /// Endpoint source-domain, degeneracy, budget and cancellation refusals keep
    /// populations, fill, carry, conversions, angle field and receipts unchanged.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), Error> {
        checkpoint(gate)?;
        let mut lattice = self.flow.lattice.clone();
        lattice.step_with_gate(gate)?;
        let geometry =
            lattice.planar_wall_fill(self.wall, self.controls.maximum_wall_cells, gate)?;
        let samples = self.resolve_samples(&geometry, Some(&self.geometry), gate)?;
        lattice.set_planar_contact_field(
            self.wall,
            samples
                .iter()
                .map(|(&p, s)| (p, s.applied_angle_radians))
                .collect(),
            gate,
        )?;
        checkpoint(gate)?;
        self.flow.lattice = lattice;
        self.geometry = geometry;
        self.samples = samples;
        Ok(())
    }

    fn resolve_samples(
        &self,
        geometry: &BTreeMap<[usize; 3], WallFill3>,
        previous: Option<&BTreeMap<[usize; 3], WallFill3>>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 3], FreeSurfaceWettingSample>, Error> {
        let mut result = BTreeMap::new();
        for (&position, current) in geometry {
            checkpoint(gate)?;
            let old_fill = previous.map_or(current.fill, |map| {
                map.get(&position).map_or(0.0, |g| g.fill)
            });
            let speed = contact_speed(
                current.fill - old_fill,
                current.tangent_gradient,
                self.flow.scaling.cell_size_m(),
                self.flow.scaling.timestep_s(),
            )?;
            let point = self
                .source
                .query()
                .with(DYNAMIC_WETTING_SPEED_AXIS, speed.abs())?;
            let material = self.source.at_query_point(self.card, &point)?;
            let (receding, advancing) = bracket(&material)?;
            let (regime, applied_angle_radians) = select_angle(
                speed,
                self.controls.pinning_speed_m_s,
                current.angle_radians,
                receding,
                advancing,
            );
            result.insert(
                position,
                FreeSurfaceWettingSample {
                    signed_speed_m_s: speed,
                    regime,
                    geometric_angle_radians: current.angle_radians,
                    applied_angle_radians,
                    material,
                },
            );
        }
        Ok(result)
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}

// At a fixed point of a translating fill profile, d(phi)/dt = V |grad_t(phi)|
// for outward liquid normal -grad_t(phi)/|grad_t(phi)|. Gradient here is per
// LATTICE CELL, so dx/dt converts its displacement to physical metres/second.
fn contact_speed(
    delta_fill: f64,
    gradient_per_cell: f64,
    dx_m: f64,
    dt_s: f64,
) -> Result<f64, Error> {
    if gradient_per_cell == 0.0 {
        return if delta_fill == 0.0 {
            Ok(0.0)
        } else {
            Err(Error::Invalid(
                "moving contact geometry has no resolved tangential gradient",
            ))
        };
    }
    let speed = (delta_fill / gradient_per_cell) * (dx_m / dt_s);
    if !speed.is_finite() || delta_fill != 0.0 && speed == 0.0 {
        return Err(Error::Invalid("contact-line speed is unrepresentable"));
    }
    Ok(speed)
}

fn select_angle(
    speed: f64,
    deadband: f64,
    geometric: f64,
    receding: f64,
    advancing: f64,
) -> (ContactLineRegime2, f64) {
    if speed > deadband {
        (ContactLineRegime2::Advancing, advancing)
    } else if speed < -deadband {
        (ContactLineRegime2::Receding, receding)
    } else {
        (
            ContactLineRegime2::Pinned,
            geometric.clamp(receding, advancing),
        )
    }
}

fn bracket(material: &ResolvedInterfaceStatePoint) -> Result<(f64, f64), Error> {
    let angle = |name| {
        material
            .property(name)
            .map(|p| p.value_si())
            .ok_or(Error::Invalid(
                "both original contact-angle endpoint claims are required",
            ))
    };
    let receding = angle(RECEDING_CONTACT_ANGLE_PROPERTY)?;
    let advancing = angle(ADVANCING_CONTACT_ANGLE_PROPERTY)?;
    if receding > advancing {
        return Err(Error::Invalid("receding angle exceeds advancing angle"));
    }
    Ok((receding, advancing))
}

// The contact history owns this speed only. Temperature and fixed environment
// remain the frozen liquid context; frequency, strain, stress, or another rate
// cannot silently become constants merely because an initial sample exists.
fn require_owned_axes(source: &ResolvedInterfaceStatePoint) -> Result<(), Error> {
    for property in source.properties() {
        if !property.has_only_owned_and_fixed_external_conditions(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
        ) {
            return Err(Error::Invalid(
                "wetting source needs its owned speed and canonical fixed external conditions",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_translating_fill_has_physical_speed_and_refinement_invariance() {
        // Independent physical plane phi(x,t)=1/2 - a*(x - V*t), a [1/m].
        // Its lattice derivative is a*dx and its fixed-cell change is a*V*dt.
        for speed in [-0.021, 0.0, 0.037] {
            for dx in [0.001, 0.0005] {
                for dt in [0.0001, 0.000025] {
                    let slope_per_m = 31.0;
                    let actual =
                        contact_speed(slope_per_m * speed * dt, slope_per_m * dx, dx, dt).unwrap();
                    assert!((actual - speed).abs() < 1e-16);
                }
            }
        }
        assert_eq!(contact_speed(0.0, 0.0, 0.001, 0.0001).unwrap(), 0.0);
        assert!(contact_speed(0.001, 0.0, 0.001, 0.0001).is_err());
        assert!(contact_speed(1.0, f64::MIN_POSITIVE, 1.0, 1e-4).is_err());
    }

    #[test]
    fn g1_direction_and_closed_deadband_use_the_sourced_hysteresis_bracket() {
        let receding = 0.9;
        let advancing = 2.1;
        assert_eq!(
            select_angle(0.02, 0.01, 1.5, receding, advancing),
            (ContactLineRegime2::Advancing, advancing)
        );
        assert_eq!(
            select_angle(-0.02, 0.01, 1.5, receding, advancing),
            (ContactLineRegime2::Receding, receding)
        );
        for speed in [-0.01, 0.0, 0.01] {
            assert_eq!(
                select_angle(speed, 0.01, 1.5, receding, advancing),
                (ContactLineRegime2::Pinned, 1.5)
            );
            assert_eq!(
                select_angle(speed, 0.01, 0.3, receding, advancing),
                (ContactLineRegime2::Pinned, receding)
            );
            assert_eq!(
                select_angle(speed, 0.01, 2.9, receding, advancing),
                (ContactLineRegime2::Pinned, advancing)
            );
        }
    }
}
