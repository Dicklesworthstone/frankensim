//! One source/SGS constitutive root in the original moving-wall transaction.

use super::*;
use crate::d3q19::freesurface3::PlanarWallFieldControls3;
use crate::material::free_surface::{
    MaterialPlanarWallFieldReport3, MaterialPlanarWallVelocityField3,
};

impl<'s, 'm> MaterialFreeSurfaceRheology3<'s, 'm> {
    /// Bind source rheology and SGS using the prescribed wall's current SI
    /// velocity. This supplies the actual halfway boundary to WALE/Vreman,
    /// including a two-layer liquid film whose liquid-only stencil refuses.
    /// Gas and empty interface cells never supply a velocity or material rate.
    ///
    /// The original joint constitutive root still uses force-corrected kinetic
    /// stress for molecular strain; spatial rotation enters only the selected
    /// full-gradient eddy closure. No populations or physical time advance.
    /// Supply the actual held wall field again on each wall step; it is not
    /// silently retained as a future wall actuation. Native flooded-plane and
    /// wall-exchange admission remain part of the subsequent stepping path.
    pub fn with_sgs_and_planar_wall_velocity_field_si(
        self,
        config: MaterialFreeSurfaceSgsConfig3,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<MaterialFreeSurfaceRheologySgs3<'s, 'm>, Error> {
        poll(gate)?;
        let (velocities, origin) = self
            .flow
            .prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        self.bind_sgs(
            config,
            gate,
            kinematics,
            workspace,
            Some(PlanarWallVelocityField3 {
                wall: motion.wall,
                moment_origin: origin,
                tangential_velocity: &velocities,
            }),
        )
    }
}

impl MaterialFreeSurfaceRheologySgs3<'_, '_> {
    /// Advance the actual source-resolved liquid with a moving planar wall
    /// and the already bound SGS closure. SI wall exchange is returned first;
    /// the second report contains the input state's physical viscous powers.
    ///
    /// The current prescribed wall velocity enters the same WALE/Vreman
    /// gradient as the native reflection. A fresh joint molecular/eddy root
    /// supplies total shear relaxation AND the original source's independent
    /// bulk response before that collision. Changing the wall field therefore
    /// cannot leave stale molecular or bulk viscosity in the applied operator.
    /// Smagorinsky keeps its stress-based joint root without a gradient cache.
    ///
    /// Actual endpoint pressure and Guo-corrected stress re-query the same
    /// source before publication, with this tick's held wall field supplying
    /// the endpoint gradient. The reusable relaxation workspace then contains
    /// endpoint total tau; material samples retain molecular tau separately.
    /// Zero coefficient reproduces the original molecular wall trajectory.
    ///
    /// One native candidate carries collision, gas reconstruction, mass
    /// transport and wall impulse/work. No thermal route log or second fluid
    /// clone is allocated. Invalid wall data, either source epoch, SI output,
    /// convergence or cancellation leave flow and source samples unchanged.
    /// Existing fixed, fully flooded, tangent-only wall limits apply. This is
    /// explicit isothermal LES exploration, without heating or a wall law.
    pub fn step_with_planar_wall_velocity_field_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        mut kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<
        (
            MaterialPlanarWallFieldReport3,
            MaterialFreeSurfaceSgsReport3,
        ),
        Error,
    > {
        poll(gate)?;
        let (velocities, origin) = self
            .flow
            .flow
            .prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let wall = PlanarWallVelocityField3 {
            wall: motion.wall,
            moment_origin: origin,
            tangential_velocity: &velocities,
        };
        let units = SgsUnits::new(&self.flow.flow.scaling)?;
        let disabled = self.native_config.model.coefficient() == 0.0;
        let (samples, report) = self.prepare_with_wall(
            self.lattice(),
            gate,
            kinematics.as_deref_mut(),
            workspace,
            Some(wall),
        )?;
        let report = units.report(report)?;
        let bulk = if disabled {
            None // Keep the exact accepted molecular source epoch.
        } else {
            self.flow.flow.scaling.source_bulk_relaxations(
                samples
                    .iter()
                    .map(|(&location, cell)| (location, cell.material.material())),
                gate,
            )?
        };
        drop(samples);
        let (mut lattice, exchange) = self
            .lattice()
            .planar_wall_field_candidate_with_relaxation_and_bulk(
                wall,
                controls,
                Some(&workspace.tau),
                bulk,
                gate,
            )?;
        let samples = if disabled {
            // Preserve the original source root and warm-start history exactly.
            self.flow.sample_lattice(&lattice, gate)?
        } else {
            self.prepare_with_wall(&lattice, gate, kinematics, workspace, Some(wall))?
                .0
        };
        self.rebind_bulk(&mut lattice, &samples, gate)?;
        let exchange =
            self.flow
                .flow
                .planar_wall_field_report_si(exchange, motion.moment_origin_m, gate)?;
        if disabled {
            workspace
                .tau
                .retain(|&(key, lane), _| samples.contains_key(&MassLocation3::Cell(key, lane)));
            for (i, (&location, cell)) in samples.iter().enumerate() {
                if i % 64 == 0 {
                    poll(gate)?;
                }
                if let MassLocation3::Cell(key, lane) = location {
                    workspace.tau.insert((key, lane), cell.sample.collision_tau);
                }
            }
        }
        poll(gate)?;
        self.flow.flow.lattice = lattice;
        self.flow.samples = samples;
        Ok((exchange, report))
    }
}
