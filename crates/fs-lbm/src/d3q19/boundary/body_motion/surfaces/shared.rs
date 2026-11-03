//! Explicit rigid carriers shared by several moving missing-film patches.

use super::*;
use crate::d3q19::boundary::partial_saturation::{
    collect_carrier_surfaces, prepare_carrier_surfaces,
};
use crate::d3q19::{
    CarrierFaceExchange3, CarrierLubricationExchange3, LubricationPlanePatch3,
    SharedLubricationSurface3,
};

/// A prescribed reservoir or a body-fixed patch on an indexed rigid carrier.
///
/// Several attached patches may name the same carrier. Its complete inertia,
/// authored load, free rotation and received momentum are owned once. Equal
/// body values at different indices still identify independent physical bodies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MovingCarrierSurface3 {
    /// Original prescribed translation and axial spin, advanced on acceptance.
    Prescribed(LubricationPlane3),
    /// Local planar missing-film geometry attached to one finite carrier.
    Attached {
        /// Stable index into the separately supplied finite carrier registry.
        carrier: usize,
        /// Point relative to the COM and normal in the carrier's body axes.
        patch_body: LubricationPlanePatch3,
    },
}

impl MovingCarrierSurface3 {
    /// Actual world patch geometry at the supplied accepted carrier state.
    ///
    /// # Errors
    /// Refuses an invalid registry index or unrepresentable transformed geometry.
    pub fn geometry(
        &self,
        carriers: &[MovingPartialSaturationBody3],
    ) -> Result<LubricationPlanePatch3> {
        match *self {
            Self::Prescribed(plane) => Ok(LubricationPlanePatch3::from_plane(plane)),
            Self::Attached {
                carrier,
                patch_body,
            } => {
                let body = carrier_at(carriers, carrier)?;
                let arm = array(native(body.orientation.rotate(vector(patch_body.point())))?);
                let normal = array(native(
                    body.orientation.rotate(vector(patch_body.normal())),
                )?);
                patch_body.at(core::array::from_fn(|i| body.center[i] + arm[i]), normal)
            }
        }
    }

    /// Actual rigid point velocity, in lattice units, without advancing state.
    ///
    /// # Errors
    /// Refuses an invalid index, nonfinite point or inadmissible local speed.
    pub fn velocity_at(
        &self,
        carriers: &[MovingPartialSaturationBody3],
        point: [f64; 3],
    ) -> Result<[f64; 3]> {
        if !point.into_iter().all(f64::is_finite) {
            return Err(invalid("nonfinite shared surface observation point"));
        }
        match *self {
            Self::Prescribed(plane) => plane.velocity_at(point),
            Self::Attached { carrier, .. } => {
                let body = carrier_at(carriers, carrier)?;
                let velocity = point_velocity(
                    body.velocity,
                    body.angular_velocity_world()?,
                    core::array::from_fn(|i| point[i] - body.center[i]),
                );
                if !low_mach(velocity) {
                    return Err(invalid("shared surface point speed is inadmissible"));
                }
                Ok(velocity)
            }
        }
    }

    fn owner(self) -> Option<usize> {
        match self {
            Self::Prescribed(_) => None,
            Self::Attached { carrier, .. } => Some(carrier),
        }
    }

    fn translated(self, duration: f64) -> Result<Self> {
        match self {
            Self::Prescribed(plane) => Ok(Self::Prescribed(plane.translated(duration)?)),
            Self::Attached { .. } => Ok(self),
        }
    }
}

/// One physical carrier's aggregate film transfer and actual moving trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingCarrierExchange3 {
    /// All attached patch reactions combined once at the native midpoint.
    pub film: CarrierLubricationExchange3,
    /// Actual pose evolution, fixed-origin momentum, storage and free-drift defect.
    pub motion: MovingFinitePlaneExchange3,
}

/// One shared moving sphere/fluid step with distinct patch and carrier accounts.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MovingLubricatedCarriersExchange3<const N: usize> {
    /// Original moving sphere/fluid exchange, including every film mode.
    pub motion: MovingPartialSaturationBodiesExchange3<N>,
    /// Per-patch reaction and work in the supplied surface order.
    pub surfaces: Vec<CarrierFaceExchange3>,
    /// One aggregate trajectory account per physical registry body.
    pub carriers: Vec<MovingCarrierExchange3>,
}

/// One physical dry-body/fluid step with shared finite carriers.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct DryLubricatedCarriersExchange3<const N: usize> {
    /// Existing dry sphere/exterior-fluid subsystem, retaining all film reactions.
    pub motion: DryBodiesExchange3<N>,
    /// Per-patch received wrench/work; no repeated carrier kinetic storage.
    pub surfaces: Vec<CarrierFaceExchange3>,
    /// One actual trajectory and aggregate film response per finite carrier.
    pub carriers: Vec<MovingCarrierExchange3>,
}

impl BoundaryGrid3 {
    /// Advance spheres, shared rigid film carriers and fluid in one native solve.
    ///
    /// Each carrier takes one free half-drift, one aggregate momentum update,
    /// and one final half-drift. Every attached patch contributes to the SAME
    /// full-inertia mobility, including cross-patch resistance entries. Loads
    /// contain one world COM impulse/couple per carrier, separate from film.
    /// Prescribed patches carry their original actuator work and no load slot.
    ///
    /// At most 64 carriers and 64 patches are admitted, with at least one patch
    /// per carrier. Inactive films still retain a carrier's free/authored motion.
    /// Original sphere, film, speed and conservative sweep limits apply. These
    /// patches do not create resolved fluid walls, finite extents, contact, a
    /// GCL, displaced volume, inferred gravity/buoyancy or a thermal destination.
    ///
    /// # Errors
    /// Refuses invalid/unused carrier associations, load count, unrepresentable
    /// motion, native solve/sweep failures or cancellation before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_spheres_and_carriers<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        carriers: &mut [MovingPartialSaturationBody3],
        carrier_external: &[[f64; 6]],
        surfaces: &mut [MovingCarrierSurface3],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingLubricatedCarriersExchange3<N>> {
        model.admit(radii)?;
        let prepared = PreparedSurfaces::new_shared(carriers, carrier_external, surfaces, gate)?;
        let shapes = radii.map(|r| BodyFrameEllipsoid3::sphere(r).expect("admitted radius"));
        let (report, candidate) = self.moving_partial_saturation_resistance_transition(
            samplers,
            core::array::from_fn(|a| &shapes[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            gate,
            Some((model, &prepared.initial_planes)),
            true,
            Some(&prepared),
        )?;
        let candidate = candidate
            .expect("shared surfaces prepare an accepted candidate")
            .shared
            .expect("shared surface candidate retains its registry");
        carriers.copy_from_slice(&candidate.carrier_states);
        surfaces.copy_from_slice(&candidate.surface_states);
        Ok(MovingLubricatedCarriersExchange3 {
            motion: report.motion,
            surfaces: candidate.surface_reports,
            carriers: candidate.carrier_reports,
        })
    }

    /// Advance physical dry spheres and shared carriers on the original fluid.
    ///
    /// Every retained-fluid Newton probe starts from the same accepted registry;
    /// carrier recoil and all cross-patch reactions are solved inside that probe.
    /// Only the accepted fluid, spheres, carriers and prescribed patch positions
    /// publish. Carrier momentum/kinetic storage appears once in the registry
    /// report; individual patches retain only their received wrench and work.
    ///
    /// # Errors
    /// Inherits the shared moving-surface bounds and original dry-body budgets.
    /// Any refusal preserves every accepted physical state and surface position.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres_and_carriers<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        carriers: &mut [MovingPartialSaturationBody3],
        carrier_external: &[[f64; 6]],
        surfaces: &mut [MovingCarrierSurface3],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedCarriersExchange3<N>> {
        model.admit(radii)?;
        let prepared = PreparedSurfaces::new_shared(carriers, carrier_external, surfaces, gate)?;
        let shapes = radii.map(|r| BodyFrameEllipsoid3::sphere(r).expect("admitted radius"));
        let (report, candidate) = self.step_with_dry_partial_saturation_resistance(
            samplers,
            core::array::from_fn(|a| &shapes[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            config,
            gate,
            Some((model, &prepared.initial_planes)),
            Some(&prepared),
        )?;
        let candidate = candidate
            .expect("shared surfaces prepare an accepted candidate")
            .shared
            .expect("shared surface candidate retains its registry");
        carriers.copy_from_slice(&candidate.carrier_states);
        surfaces.copy_from_slice(&candidate.surface_states);
        Ok(DryLubricatedCarriersExchange3 {
            motion: report.motion,
            surfaces: candidate.surface_reports,
            carriers: candidate.carrier_reports,
        })
    }
}

pub(super) struct PreparedSharedSurfaces {
    initial_carriers: Vec<MovingPartialSaturationBody3>,
    midpoint_carriers: Vec<MovingPartialSaturationBody3>,
    external: Vec<[f64; 6]>,
    initial_surfaces: Vec<MovingCarrierSurface3>,
    midpoint_surfaces: Vec<SharedLubricationSurface3>,
    midpoint_world_bodies: Vec<PartialSaturationBody3>,
}

pub(super) struct SharedSurfaceCandidate {
    carrier_states: Vec<MovingPartialSaturationBody3>,
    surface_states: Vec<MovingCarrierSurface3>,
    carrier_reports: Vec<MovingCarrierExchange3>,
    surface_reports: Vec<CarrierFaceExchange3>,
}

impl SharedSurfaceCandidate {
    pub(super) fn empty() -> Self {
        Self {
            carrier_states: Vec::new(),
            surface_states: Vec::new(),
            carrier_reports: Vec::new(),
            surface_reports: Vec::new(),
        }
    }
}

impl PreparedSurfaces {
    fn new_shared(
        carriers: &[MovingPartialSaturationBody3],
        external: &[[f64; 6]],
        surfaces: &[MovingCarrierSurface3],
        gate: &CancelGate,
    ) -> Result<Self> {
        checkpoint(gate)?;
        if carriers.len() > 64 || surfaces.len() > 64 || carriers.len() != external.len() {
            return Err(invalid(
                "invalid shared moving carrier, patch or load count",
            ));
        }
        let mut used = [false; 64];
        for &surface in surfaces {
            if let Some(index) = surface.owner() {
                carrier_at(carriers, index)?;
                used[index] = true;
            }
        }
        if used[..carriers.len()].iter().any(|&present| !present) {
            return Err(invalid(
                "every shared moving carrier requires an attached patch",
            ));
        }
        let mut initial_carriers = reserved(carriers.len())?;
        let mut midpoint_carriers = reserved(carriers.len())?;
        let mut initial_world_bodies = reserved(carriers.len())?;
        let mut midpoint_world_bodies = reserved(carriers.len())?;
        let mut loads = reserved(carriers.len())?;
        for (&body, &load) in carriers.iter().zip(external) {
            checkpoint(gate)?;
            if !load.into_iter().all(f64::is_finite) {
                return Err(invalid("nonfinite shared carrier external impulse"));
            }
            initial_carriers.push(body);
            let midpoint_body = drift(body, 0.5)?;
            midpoint_carriers.push(midpoint_body);
            initial_world_bodies.push(world_body(body)?);
            midpoint_world_bodies.push(world_body(midpoint_body)?);
            loads.push(load);
        }
        let mut initial_surfaces = reserved(surfaces.len())?;
        let mut initial = reserved(surfaces.len())?;
        let mut midpoint = reserved(surfaces.len())?;
        let mut initial_descriptors = reserved(surfaces.len())?;
        let mut midpoint_surfaces = reserved(surfaces.len())?;
        for &surface in surfaces {
            checkpoint(gate)?;
            let mid = surface.translated(0.5)?;
            initial_surfaces.push(surface);
            initial.push(flatten(surface, &initial_carriers)?);
            midpoint.push(flatten(mid, &midpoint_carriers)?);
            initial_descriptors.push(frozen(surface, &initial_carriers)?);
            midpoint_surfaces.push(frozen(mid, &midpoint_carriers)?);
        }
        let initial_planes =
            prepare_carrier_surfaces(&initial_descriptors, &initial_world_bodies, &loads, gate)?;
        let midpoint_planes =
            prepare_carrier_surfaces(&midpoint_surfaces, &midpoint_world_bodies, &loads, gate)?;
        admit_planes(&initial_planes)?;
        admit_planes(&midpoint_planes)?;
        admit_distinct_sweep(&initial, &midpoint, gate)?;
        checkpoint(gate)?;
        Ok(Self {
            initial,
            midpoint,
            external: Vec::new(),
            initial_planes,
            midpoint_planes,
            shared: Some(PreparedSharedSurfaces {
                initial_carriers,
                midpoint_carriers,
                external: loads,
                initial_surfaces,
                midpoint_surfaces,
                midpoint_world_bodies,
            }),
        })
    }
}

impl PreparedSharedSurfaces {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish<const N: usize>(
        &self,
        exchanges: &[PlaneLubricationExchange3],
        before: [[f64; 3]; N],
        after: [[f64; 3]; N],
        radii: [f64; N],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<SurfaceCandidate> {
        if exchanges.len() != self.initial_surfaces.len() {
            return Err(invalid("shared surface/carrier response count mismatch"));
        }
        let (surface_reports, carrier_exchanges) = collect_carrier_surfaces(
            &self.midpoint_surfaces,
            &self.midpoint_world_bodies,
            exchanges,
            gate,
        )?;
        let mut carrier_states = reserved(carrier_exchanges.len())?;
        let mut kicked_carriers = reserved(carrier_exchanges.len())?;
        let mut carrier_reports = reserved(carrier_exchanges.len())?;
        for (index, &film) in carrier_exchanges.iter().enumerate() {
            checkpoint(gate)?;
            let aggregate = PlaneLubricationExchange3 {
                impulse: film.impulse,
                angular_impulse: film.angular_impulse,
                work: film.work,
                finite_body: Some(film.motion),
                ..PlaneLubricationExchange3::default()
            };
            let (kicked, accepted, motion) = finish_carrier(
                self.initial_carriers[index],
                self.midpoint_carriers[index],
                aggregate,
                self.external[index],
                origin,
            )?;
            kicked_carriers.push(kicked);
            carrier_states.push(accepted);
            carrier_reports.push(MovingCarrierExchange3 { film, motion });
        }
        let mut surface_states = reserved(exchanges.len())?;
        let mut exchanged = reserved(exchanges.len())?;
        let mut states = reserved(exchanges.len())?;
        for &surface in &self.initial_surfaces {
            checkpoint(gate)?;
            let end = surface.translated(1.0)?;
            exchanged.push(flatten(surface.translated(0.5)?, &kicked_carriers)?);
            states.push(flatten(end, &carrier_states)?);
            surface_states.push(end);
        }
        admit_spheres(&exchanged, &states, before, after, radii, gate)?;
        admit_distinct_sweep(&exchanged, &states, gate)?;
        checkpoint(gate)?;
        Ok(SurfaceCandidate {
            states,
            reports: Vec::new(),
            shared: Some(SharedSurfaceCandidate {
                carrier_states,
                surface_states,
                carrier_reports,
                surface_reports,
            }),
        })
    }
}

fn carrier_at(
    carriers: &[MovingPartialSaturationBody3],
    index: usize,
) -> Result<MovingPartialSaturationBody3> {
    carriers
        .get(index)
        .copied()
        .ok_or_else(|| invalid("moving film carrier index out of range"))
}

// The existing surface geometry and sweep owner receives a view of each face.
// These views do not advance or own a second independent rigid state.
fn flatten(
    surface: MovingCarrierSurface3,
    carriers: &[MovingPartialSaturationBody3],
) -> Result<MovingLubricationSurface3> {
    match surface {
        MovingCarrierSurface3::Prescribed(plane) => {
            Ok(MovingLubricationSurface3::Prescribed(plane))
        }
        MovingCarrierSurface3::Attached {
            carrier,
            patch_body,
        } => {
            let mut plane = MovingFiniteLubricationPlane3::new(
                carrier_at(carriers, carrier)?,
                patch_body.point(),
                patch_body.normal(),
                patch_body.cutoff_gap(),
                patch_body.roughness_gap(),
            )?;
            if let Some((cutoff, roughness)) = patch_body.tangential_gaps() {
                plane = plane.with_tangential(cutoff, roughness)?;
            }
            Ok(MovingLubricationSurface3::Finite(plane))
        }
    }
}

fn frozen(
    surface: MovingCarrierSurface3,
    carriers: &[MovingPartialSaturationBody3],
) -> Result<SharedLubricationSurface3> {
    match surface {
        MovingCarrierSurface3::Prescribed(plane) => {
            Ok(SharedLubricationSurface3::Prescribed(plane))
        }
        MovingCarrierSurface3::Attached { carrier, .. } => {
            Ok(surface.geometry(carriers)?.carried(carrier))
        }
    }
}
