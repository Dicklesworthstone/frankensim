//! Original link-resolved rigid wall work on the thermal liquid trajectory.

use super::*;

mod pair;
pub use pair::{FreeSurfaceRigidWallPairEnthalpyReport, FreeSurfaceRigidWallPairViscousReport};

/// A complete unpublished thermal tick and its actual rigid wall response.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceRigidWallEnthalpyReport {
    pub flow: FreeSurfaceEnthalpyReport,
    /// One native reflection supplies impulse, moment, work and the full
    /// translation/spin response about the supplied common origin.
    pub wall: MaterialPlanarRigidWallReport3,
}

/// The same tick with the optionally selected initial-state molecular heat samples.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceRigidWallViscousReport {
    pub flow: FreeSurfaceRigidWallEnthalpyReport,
    /// Empty when the source descriptor disables molecular heating.
    pub applied: Vec<viscous::FreeSurfaceViscousCell>,
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Propose a thermal tick with two tangent translations and normal spin.
    ///
    /// The original rigid operator samples velocity at each actual reflected
    /// link, preserving its reciprocal angular impulse and rotational work.
    /// Heat refreshes the original cellwise material relaxation and force
    /// frame before collision; the SAME tick supplies the native mass routes
    /// used for extensive enthalpy transport. Endpoint material, rheology,
    /// wetting, source speed, wall mass and cancellation gates precede return.
    ///
    /// This borrows the original thermal state and clones its native fluid
    /// exactly once. The original source bindings remain attached, and no
    /// prepared conduction/electrical operator is copied: those coupled
    /// boundaries retain their ordinary owning-timeline refusal.
    ///
    /// Every candidate begins at the same accepted state and advances one
    /// physical tick only if the outer owner accepts it. A finite-body owner
    /// can therefore measure a pilot response and replay its solved midpoint
    /// velocity without advancing heat or time twice. Work is mechanical, not
    /// an additional caloric source. Fixed flooded-plane and single-phase
    /// assumptions remain; normal coverage, tilt and finite edges are absent.
    pub fn planar_rigid_wall_motion_candidate_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallEnthalpyReport), Error> {
        let (candidate, report) = self.planar_rigid_wall_motion_candidate_with_sources_si(
            motion,
            maximum_relative_mass_defect,
            FreeSurfaceWallSources {
                heat,
                ..FreeSurfaceWallSources::default()
            },
            gate,
        )?;
        Ok((candidate, report.flow))
    }

    /// Use the original force-corrected strain and sourced molecular viscosity.
    ///
    /// Held initial-state molecular heat is applied once in this candidate,
    /// separately from authored heat and wall work. This retains the existing
    /// first-order source split and requires timestep refinement for strong
    /// thermal feedback. Numerical flow loss is never turned into heat.
    pub fn planar_rigid_wall_motion_candidate_with_viscous_heating_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallViscousReport), Error> {
        self.planar_rigid_wall_motion_candidate_with_sources_si(
            motion,
            maximum_relative_mass_defect,
            FreeSurfaceWallSources {
                heat,
                molecular_heating: true,
                gas_pressure: None,
            },
            gate,
        )
    }

    /// Propose one thermal rigid-wall tick with all selected sources.
    ///
    /// The immutable pressure field and authored heat are reused unchanged for
    /// any finite-body pilot and replay. Each candidate samples molecular heat
    /// from the same accepted initial state when selected, then advances one
    /// original collision and extensive enthalpy transport. Gas pressure and
    /// wall work remain mechanical loading; neither is added to heat. The
    /// returned molecular sample list is empty when that source is disabled.
    pub fn planar_rigid_wall_motion_candidate_with_sources_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        sources: FreeSurfaceWallSources<'_>,
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallViscousReport), Error> {
        let native = self
            .flow
            .prepare_planar_rigid_wall_motion_si(motion, gate)?;
        let prepared = self.prepare_wall_sources(sources, gate)?;
        let (candidate, flow) =
            sources.with_lattice_gas_pressure(self.flow.scaling, |gas_pressure| {
                self.rigid_wall_enthalpy_candidate(
                    ThermalWallBoundary::Rigid {
                        native,
                        motion_si: motion,
                        maximum_relative_mass_defect,
                    },
                    prepared.heated,
                    prepared.external,
                    prepared.viscous,
                    prepared.absolute_heat,
                    gas_pressure,
                    gate,
                )
            })?;
        Ok((
            candidate,
            FreeSurfaceRigidWallViscousReport {
                flow,
                applied: prepared.applied,
            },
        ))
    }

    fn rigid_wall_enthalpy_candidate(
        &self,
        boundary: ThermalWallBoundary<'_>,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        absolute_heat: f64,
        gas_pressure: Option<(f64, GasPressureField<'_>)>,
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallEnthalpyReport), Error> {
        let (candidate, flow, wall) = self.wall_enthalpy_candidate(
            boundary,
            heated,
            external,
            viscous,
            absolute_heat,
            gas_pressure,
            gate,
        )?;
        let wall = match wall {
            ThermalWallReport::Rigid(report) => report,
            _ => unreachable!("rigid boundary retains its original wall receipt"),
        };
        Ok((candidate, FreeSurfaceRigidWallEnthalpyReport { flow, wall }))
    }

    // One construction seam retains the same source and thermal history for
    // either rigid boundary owner; no population or prepared operator is cloned
    // after the original thermal pipeline has produced its single candidate.
    fn wall_enthalpy_candidate(
        &self,
        boundary: ThermalWallBoundary<'_>,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        absolute_heat: f64,
        gas_pressure: Option<(f64, GasPressureField<'_>)>,
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceEnthalpyReport, ThermalWallReport), Error> {
        // Both entry points use prepare_heat, whose original electrical
        // ownership gate prevents dropping or cloning a prepared operator.
        let tick = self.prepare_enthalpy_tick_with_wall(
            heated,
            external,
            viscous,
            0.0,
            absolute_heat,
            None,
            None,
            None,
            gas_pressure,
            Some(boundary),
            gate,
        )?;
        let wall = tick
            .wall
            .expect("thermal boundary retains its original wall receipt");
        let candidate = Self {
            flow: MaterialFreeSurface3 {
                lattice: tick.lattice,
                scaling: self.flow.scaling,
            },
            law: self.law,
            caloric: self.caloric,
            controls: self.controls,
            time_origin_s: self.time_origin_s,
            accepted_steps: tick.accepted_steps,
            density_accuracy_kg_m3: self.density_accuracy_kg_m3,
            volume_slab: None,
            h: tick.h,
            states: tick.states,
            phase_conduction: self.phase_conduction,
            buoyancy: self.buoyancy,
            rheology: self.rheology.clone(),
            rheology_samples: tick.rheology_samples,
            wetting: tick.wetting,
            thermocapillary: self.thermocapillary,
            surface_temperature_forces: tick.surface_temperature_forces,
            #[cfg(feature = "enthalpy-coupling")]
            electrical: None,
        };
        poll(gate)?;
        Ok((candidate, tick.report, wall))
    }
}
