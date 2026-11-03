//! Dense contact-line history using the same source and Cox law as sparse VOF.

use super::super::wetting::{
    CoxVoinovLaw, FreeSurfaceWettingControls, FreeSurfaceWettingError as Error,
    FreeSurfaceWettingSample, WettingPhysics, checkpoint, contact_speed,
};
use super::{
    FreeSurface, MaterialFreeSurface2, MaterialFreeSurfaceMomentumError2,
    MaterialFreeSurfaceMomentumReport2, MaterialGasPressureLink2, PlanarWall2,
};
use crate::freesurface::FreeSurfaceMomentumError2;
use crate::freesurface::wetting::WallFill2;
use fs_exec::CancelGate;
use fs_matdb::InterfaceSystemCard;
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use std::collections::BTreeMap;

/// Source-pinned advancing/receding adhesion on the actual dense D2Q9 wall.
///
/// Accepted wall-adjacent fill changes determine signed speed in m/s. The
/// first-cell/backward-time estimate drives the next step's capillary ghosts.
/// Pinned means projection into the sourced bracket, not exact immobilization.
/// The existing single-liquid, isothermal, unit-depth flow remains the owner
/// of populations, conversions and mass; no surface-energy/heat law is inferred.
pub struct MaterialFreeSurfaceWetting2<'s, 'm> {
    flow: MaterialFreeSurface2<'s>,
    physics: WettingPhysics<'m>,
    wall: PlanarWall2,
    geometry: BTreeMap<[usize; 2], WallFill2>,
    samples: BTreeMap<[usize; 2], FreeSurfaceWettingSample>,
}

impl<'s> MaterialFreeSurface2<'s> {
    /// Attach one named solid/liquid interface and its original endpoint claims
    /// to a complete stationary wall. Existing equilibrium/legacy wetting refuses.
    /// The substrate is surface A; surface B must be this exact liquid identity.
    /// Other source coordinates stay at the original liquid query point.
    pub fn with_dynamic_planar_wetting<'m>(
        mut self,
        wall: PlanarWall2,
        card: &'m InterfaceSystemCard,
        selection: MaterialPropertySelection,
        controls: FreeSurfaceWettingControls,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWetting2<'s, 'm>, Error> {
        checkpoint(gate)?;
        if self.scaling.planar_wetting.is_some() || !self.lattice.admits_dynamic_planar_wetting() {
            return Err(Error::Invalid(
                "dynamic planar wetting requires otherwise neutral walls",
            ));
        }
        let physics = WettingPhysics::try_new(&self.scaling, card, selection, controls, gate)?;
        let geometry = self
            .lattice
            .planar_wall_fill(wall, controls.maximum_wall_cells, gate)?;
        let mut samples = BTreeMap::new();
        for (&p, current) in &geometry {
            checkpoint(gate)?;
            samples.insert(p, physics.sample(current.angle_radians, 0.0, None, gate)?);
        }
        install(&mut self.lattice, wall, &samples, gate)?;
        checkpoint(gate)?;
        Ok(MaterialFreeSurfaceWetting2 {
            flow: self,
            physics,
            wall,
            geometry,
            samples,
        })
    }
}

impl MaterialFreeSurfaceWetting2<'_, '_> {
    /// Apply the existing full-angle Cox hydrodynamic prediction to STATIC
    /// source endpoints. Explicit microscopic/observation lengths, phase
    /// viscosity ratio, low-Ca and angular bounds belong to `law`.
    /// A measured speed-dependent endpoint curve refuses double rate counting.
    pub fn with_cox_voinov(mut self, law: CoxVoinovLaw, gate: &CancelGate) -> Result<Self, Error> {
        let physics = self
            .physics
            .with_cox_voinov(law, &self.flow.scaling, gate)?;
        let mut samples = BTreeMap::new();
        for (&p, current) in &self.geometry {
            checkpoint(gate)?;
            samples.insert(
                p,
                physics.sample(
                    current.angle_radians,
                    self.samples[&p].signed_speed_m_s,
                    None,
                    gate,
                )?,
            );
        }
        install(&mut self.flow.lattice, self.wall, &samples, gate)?;
        checkpoint(gate)?;
        self.physics = physics;
        self.samples = samples;
        Ok(self)
    }

    /// Actual physical-unit flow and immutable native state.
    #[must_use]
    pub const fn flow(&self) -> &MaterialFreeSurface2<'_> {
        &self.flow
    }

    /// Original stationary endpoint binding; selected claims never reselect.
    #[must_use]
    pub const fn source(&self) -> &ResolvedInterfaceStatePoint {
        &self.physics.source
    }

    /// Accepted speeds, next-step angles and original query receipts.
    #[must_use]
    pub const fn samples(&self) -> &BTreeMap<[usize; 2], FreeSurfaceWettingSample> {
        &self.samples
    }

    /// Optional explicit hydrodynamic law with its cached angular integral.
    #[must_use]
    pub fn cox_voinov_law(&self) -> Option<&CoxVoinovLaw> {
        self.physics.cox_voinov_law()
    }

    /// Transactional native momentum/mass step and contact-history refresh.
    /// Material-domain, degenerate-geometry, Cox-limit or cancellation refusals
    /// preserve all accepted populations, mass/carry, conversions and receipts.
    /// Uses the same conservative conversion and SI admission as
    /// [`Self::step_with_momentum`], retaining this entry's error type.
    /// Cancellation is checked during native streaming and at every wall query;
    /// the dense collision sweep itself remains noninterruptible.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), Error> {
        self.step_with_momentum(gate)
            .map(|_| ())
            .map_err(|error| match error {
                MaterialFreeSurfaceMomentumError2::Native(FreeSurfaceMomentumError2::Cancelled) => {
                    crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into()
                }
                MaterialFreeSurfaceMomentumError2::Native(FreeSurfaceMomentumError2::Invalid(
                    reason,
                )) => Error::Invalid(reason),
                MaterialFreeSurfaceMomentumError2::Units(error) => Error::Flow(error),
                MaterialFreeSurfaceMomentumError2::Pressure(error) => Error::Flow(error.into()),
                MaterialFreeSurfaceMomentumError2::Wetting(error) => error,
            })
    }

    /// Advance the original dense solver with signed, capacity-bounded local
    /// momentum conversion, then refresh this same sourced wetting history.
    /// Native conversion, SI overflow, source-domain and cancellation refusals
    /// leave populations, mass, geometry and all angle receipts unchanged.
    /// One native candidate is shared by both operations. The returned kinetic
    /// change is numerical remap energy per depth, never physical wetting heat.
    pub fn step_with_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (lattice, native) = self.flow.lattice.momentum_candidate(gate)?;
        let report = self.flow.momentum_report_si(native)?;
        self.finish_step(lattice, gate)?;
        Ok(report)
    }

    /// Advance prescribed absolute gas pressure [Pa], conservative conversion
    /// and this same sourced contact-line history in one atomic proposal.
    /// Coordinates supplied to the pressure field are liquid cell centres [m].
    /// The flow must have an explicit fixed pressure reference. Pressure,
    /// conversion, SI-range, wetting-source and cancellation refusals leave
    /// both the accepted liquid and the complete contact history unchanged.
    pub fn step_with_gas_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 2]) -> f64,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (lattice, report) = self.flow.gas_pressure_momentum_candidate(
            gate, maximum_density_excursion, pressure_pa,
        )?;
        self.finish_step(lattice, gate)?;
        Ok(report)
    }

    /// Use independent absolute pressures [Pa] at exposed liquid/gas links,
    /// conservative conversion and sourced wetting in the same proposal.
    /// Geometry is supplied in metres. The conversion receipt does not include
    /// pressure work, gas storage or a physical contact-line heat allocation.
    pub fn step_with_gas_link_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink2) -> f64,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (lattice, report) = self.flow.gas_link_pressure_momentum_candidate(
            gate, maximum_density_excursion, pressure_pa,
        )?;
        self.finish_step(lattice, gate)?;
        Ok(report)
    }

    /// Preserve cell-centre angular momentum in each local conversion component
    /// while advancing the same sourced contact-line history. The native
    /// nonperiodic, chosen-frame and geometric-support restrictions apply.
    /// Flow, SI account and angle history are admitted as one transaction.
    pub fn step_with_angular_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (lattice, native) = self.flow.lattice.angular_momentum_candidate(gate)?;
        let report = self.flow.momentum_report_si(native)?;
        self.finish_step(lattice, gate)?;
        Ok(report)
    }

    fn finish_step(&mut self, mut lattice: FreeSurface, gate: &CancelGate) -> Result<(), Error> {
        checkpoint(gate)?;
        let geometry =
            lattice.planar_wall_fill(self.wall, self.physics.controls.maximum_wall_cells, gate)?;
        let mut samples = BTreeMap::new();
        for (&p, current) in &geometry {
            checkpoint(gate)?;
            let previous = self.geometry.get(&p).ok_or(Error::Invalid(
                "dense contact-line history topology changed",
            ))?;
            let speed = contact_speed(
                current.fill - previous.fill,
                current.tangent_gradient,
                self.flow.scaling.cell_size_m(),
                self.flow.scaling.timestep_s(),
            )?;
            samples.insert(
                p,
                self.physics
                    .sample(current.angle_radians, speed, None, gate)?,
            );
        }
        install(&mut lattice, self.wall, &samples, gate)?;
        checkpoint(gate)?;
        self.flow.lattice = lattice;
        self.geometry = geometry;
        self.samples = samples;
        Ok(())
    }
}

fn install(
    lattice: &mut FreeSurface,
    wall: PlanarWall2,
    samples: &BTreeMap<[usize; 2], FreeSurfaceWettingSample>,
    gate: &CancelGate,
) -> Result<(), Error> {
    lattice.set_planar_contact_field(
        wall,
        samples
            .iter()
            .map(|(&p, s)| (p, s.applied_angle_radians))
            .collect(),
        gate,
    )?;
    Ok(())
}
