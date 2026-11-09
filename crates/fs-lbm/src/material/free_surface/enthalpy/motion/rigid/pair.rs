//! Two rigid wall exchanges on one original material and enthalpy trajectory.

use super::*;

/// One unpublished native thermal tick and both opposing rigid-wall exchanges.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceRigidWallPairEnthalpyReport {
    /// Original extensive enthalpy transport, sources, mass and endpoint report.
    pub flow: FreeSurfaceEnthalpyReport,
    /// Both reciprocal physical wrenches and their common liquid inventory.
    /// Mechanical boundary work is not deposited as caloric heat in `flow`.
    pub walls: MaterialPlanarRigidWallPairReport3,
}

/// The same shared tick with optionally selected original molecular heating.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceRigidWallPairViscousReport {
    pub flow: FreeSurfaceRigidWallPairEnthalpyReport,
    /// Held initial-state source viscosity and force-corrected strain samples.
    /// Empty when the source descriptor disables molecular heating.
    pub applied: Vec<viscous::FreeSurfaceViscousCell>,
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Propose two opposing flooded rigid walls on ONE thermal liquid tick.
    ///
    /// Authored heat refreshes the original local viscosity/rheology and
    /// optional Boussinesq force before the same collision and both wall
    /// reflections. Their single native mass-route trace transports extensive
    /// enthalpy; material, force-frame and wetting refresh then admit the actual
    /// endpoint before either physical wrench or the thermal state can escape.
    ///
    /// Both walls use the same original incoming fluid. Their separate affine
    /// responses and common-origin impulse/moment/work accounts retain the native
    /// pair's caller ordering. The shared thermal owner clones its fluid once;
    /// each candidate starts from this unchanged accepted source history. A
    /// finite-body pilot and replay can therefore accept heat and time once.
    ///
    /// Wall work remains mechanical exchange. The original fixed fully flooded
    /// plane, low-Mach and single-phase assumptions remain. Normal motion, tilt,
    /// changing coverage and finite edges are unsupported. Prepared conduction,
    /// electrical and temperature-dependent surface boundaries retain their
    /// existing ownership refusals; no incomplete heat operator is cloned.
    pub fn planar_rigid_wall_pair_candidate_si(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallPairEnthalpyReport), Error> {
        let (candidate, report) = self.planar_rigid_wall_pair_candidate_with_sources_si(
            motions,
            maximum_relative_mass_defect,
            FreeSurfaceWallSources {
                heat,
                ..FreeSurfaceWallSources::default()
            },
            gate,
        )?;
        Ok((candidate, report.flow))
    }

    /// Include the original held initial-state molecular stress-power source.
    ///
    /// Original force-corrected strain and source viscosity determine molecular
    /// heating once per candidate, independently of the two wall works. Authored
    /// and molecular heat share the same source refresh and native mass routes.
    /// The original first-order source split requires timestep refinement for
    /// strong feedback. Numerical flow loss is never inferred to be heat.
    pub fn planar_rigid_wall_pair_candidate_with_viscous_heating_si(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallPairViscousReport), Error> {
        self.planar_rigid_wall_pair_candidate_with_sources_si(
            motions,
            maximum_relative_mass_defect,
            FreeSurfaceWallSources {
                heat,
                molecular_heating: true,
                gas_pressure: None,
            },
            gate,
        )
    }

    /// Propose both rigid walls, prescribed gas load and heat on one thermal tick.
    ///
    /// One frozen source descriptor drives the original collision and the
    /// resulting mass/enthalpy trace. Each gas link and each of the disjoint
    /// flooded walls retains its original native population owner. Finite-body
    /// pilots and final replay must reuse the same pressure field and heat;
    /// only the accepted candidate advances the liquid clock. Molecular samples
    /// are empty when disabled. Pressure and both wall works remain mechanical
    /// exchange, without inferred compressive or numerical-dissipation heat.
    pub fn planar_rigid_wall_pair_candidate_with_sources_si(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        sources: FreeSurfaceWallSources<'_>,
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallPairViscousReport), Error> {
        let boundary = self.thermal_rigid_wall_pair(motions, maximum_relative_mass_defect, gate)?;
        let prepared = self.prepare_wall_sources(sources, gate)?;
        let (candidate, flow) =
            sources.with_lattice_gas_pressure(self.flow.scaling, |gas_pressure| {
                self.rigid_wall_pair_enthalpy_candidate(
                    boundary,
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
            FreeSurfaceRigidWallPairViscousReport {
                flow,
                applied: prepared.applied,
            },
        ))
    }

    fn thermal_rigid_wall_pair(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<ThermalWallBoundary<'static>, Error> {
        let native = [
            self.flow
                .prepare_planar_rigid_wall_motion_si(motions[0], gate)?,
            self.flow
                .prepare_planar_rigid_wall_motion_si(motions[1], gate)?,
        ];
        Ok(ThermalWallBoundary::RigidPair {
            native,
            motions_si: motions,
            maximum_relative_mass_defect,
        })
    }

    fn rigid_wall_pair_enthalpy_candidate(
        &self,
        boundary: ThermalWallBoundary<'_>,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        absolute_heat: f64,
        gas_pressure: Option<(f64, GasPressureField<'_>)>,
        gate: &CancelGate,
    ) -> Result<(Self, FreeSurfaceRigidWallPairEnthalpyReport), Error> {
        let (candidate, flow, walls) = self.wall_enthalpy_candidate(
            boundary,
            heated,
            external,
            viscous,
            absolute_heat,
            gas_pressure,
            gate,
        )?;
        let walls = match walls {
            ThermalWallReport::RigidPair(report) => report,
            _ => unreachable!("rigid pair retains both original wall receipts"),
        };
        Ok((
            candidate,
            FreeSurfaceRigidWallPairEnthalpyReport { flow, walls },
        ))
    }
}
