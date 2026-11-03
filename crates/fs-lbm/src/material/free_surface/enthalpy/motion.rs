//! Wall momentum on the existing source-driven thermal free-surface owner.

use super::*;
use crate::d3q19::freesurface3::transport::MassTransport3;
use crate::d3q19::freesurface3::{
    PlanarRigidWallMotion3, PlanarRigidWallReport3, PlanarWallFieldControls3,
    PlanarWallFieldReport3, PlanarWallVelocityField3,
};
use crate::material::free_surface::{
    MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallReport3, MaterialPlanarWallFieldReport3,
    MaterialPlanarWallVelocityField3,
};

mod rigid;
pub use rigid::{FreeSurfaceRigidWallEnthalpyReport, FreeSurfaceRigidWallViscousReport};

/// One accepted thermal flow tick and its actual distributed wall exchange.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfacePlanarWallEnthalpyReport {
    /// Extensive enthalpy transported by the exact native mass routes.
    pub flow: FreeSurfaceEnthalpyReport,
    /// Reciprocal SI impulse, torque, face couples, work and pressure references.
    /// Mechanical wall work is not an additional caloric source in `flow`.
    pub wall: MaterialPlanarWallFieldReport3,
}

/// The same wall-driven tick with explicitly selected molecular stress heating.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfacePlanarWallViscousReport {
    pub flow: FreeSurfacePlanarWallEnthalpyReport,
    /// Held initial-state molecular power from the original viscosity/strain owner.
    pub applied: Vec<viscous::FreeSurfaceViscousCell>,
}

#[derive(Clone, Copy)]
pub(super) enum ThermalWallBoundary<'a> {
    Field {
        native: PlanarWallVelocityField3<'a>,
        controls: PlanarWallFieldControls3,
        moment_origin_m: [f64; 3],
    },
    Rigid {
        native: PlanarRigidWallMotion3,
        motion_si: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
    },
}

pub(super) enum NativeThermalWallReport {
    Field(PlanarWallFieldReport3),
    Rigid(PlanarRigidWallReport3),
}

pub(super) enum ThermalWallReport {
    Field(MaterialPlanarWallFieldReport3),
    Rigid(MaterialPlanarRigidWallReport3),
}

impl ThermalWallBoundary<'_> {
    pub(super) fn step(
        self,
        lattice: &mut FreeSurface3,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        maximum_routes: usize,
        gate: &CancelGate,
    ) -> Result<(MassTransport3, NativeThermalWallReport), Error> {
        match self {
            Self::Field {
                native, controls, ..
            } => {
                let (trace, report) = lattice.step_with_planar_wall_field_transport(
                    native,
                    controls,
                    tau,
                    forces,
                    maximum_routes,
                    gate,
                )?;
                Ok((trace, NativeThermalWallReport::Field(report)))
            }
            Self::Rigid {
                native,
                maximum_relative_mass_defect,
                ..
            } => {
                let (trace, report) = lattice.step_with_planar_rigid_wall_transport(
                    native,
                    maximum_relative_mass_defect,
                    tau,
                    forces,
                    maximum_routes,
                    gate,
                )?;
                Ok((trace, NativeThermalWallReport::Rigid(report)))
            }
        }
    }

    pub(super) fn finish(
        self,
        flow: &MaterialFreeSurface3<'_>,
        lattice: &FreeSurface3,
        report: NativeThermalWallReport,
        gate: &CancelGate,
    ) -> Result<ThermalWallReport, Error> {
        match (self, report) {
            (
                Self::Field {
                    native,
                    controls,
                    moment_origin_m,
                },
                NativeThermalWallReport::Field(mut report),
            ) => {
                lattice.refresh_planar_wall_field_endpoint(native, controls, &mut report, gate)?;
                Ok(ThermalWallReport::Field(flow.planar_wall_field_report_si(
                    report,
                    moment_origin_m,
                    gate,
                )?))
            }
            (
                Self::Rigid {
                    motion_si,
                    maximum_relative_mass_defect,
                    ..
                },
                NativeThermalWallReport::Rigid(mut report),
            ) => {
                lattice.refresh_planar_rigid_wall_endpoint(
                    &mut report,
                    maximum_relative_mass_defect,
                    gate,
                )?;
                Ok(ThermalWallReport::Rigid(
                    flow.planar_rigid_wall_report_si(motion_si, report, gate)?,
                ))
            }
            _ => Err(Error::Invalid(
                "native thermal wall receipt changed representation",
            )),
        }
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Advance a nonuniform tangential wall field with evolving liquid enthalpy.
    ///
    /// Authored heat refreshes the original temperature/rheology source before
    /// native collision. Its cellwise viscosity and optional Boussinesq force
    /// act in the SAME native step as wall reflection. Actual native mass routes
    /// then carry extensive enthalpy, including conversion and carry, before the
    /// accepted source/Guo frame is refreshed at the advected endpoint.
    ///
    /// Every local wall receipt uses those actual reflected populations. Wall
    /// work remains mechanical exchange; it is not deposited as enthalpy. Use
    /// the molecular-heating method below when that constitutive source is
    /// wanted. The original pressure-work and single-phase limitations remain.
    ///
    /// Geometry is a fixed, completely flooded plane. Face velocities are held
    /// for one physical tick. No normal cover/vacate, finite plate edge, moving
    /// contact line or gas-pressure actuator is added. A previously selected
    /// thermal surface-force boundary still requires its owning heat solve.
    ///
    /// Count/field budgets precede the sole native candidate clone. Native
    /// mass/positivity/Mach, thermal source/transport, actual final speed, SI
    /// representability and cancellation must all pass before populations,
    /// enthalpy, rheology, wetting and the thermal clock accept together.
    pub fn step_with_planar_wall_velocity_field_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfacePlanarWallEnthalpyReport, Error> {
        let (velocities, origin) = self
            .flow
            .prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let (heated, external, absolute_heat) = self.prepare_heat(heat, gate)?;
        self.advance_wall_enthalpy(
            ThermalWallBoundary::Field {
                native: PlanarWallVelocityField3 {
                    wall: motion.wall,
                    moment_origin: origin,
                    tangential_velocity: &velocities,
                },
                controls,
                moment_origin_m: motion.moment_origin_m,
            },
            heated,
            external,
            0.0,
            absolute_heat,
            gate,
        )
    }

    /// Drive the thermal liquid and retain its existing molecular heat source.
    ///
    /// The original force-corrected strain estimator and sourced viscosity
    /// provide held initial-state `2 mu D:D`. Authored and molecular heat then
    /// refresh the source before the single wall/flow/transport step. This is
    /// the existing first-order source split, including any selected memoryless
    /// rheology; strong thermal feedback requires timestep refinement.
    ///
    /// Molecular heat, actuator work and numerical dissipation remain separate.
    /// This does not assert exact discrete kinetic-plus-thermal balance or
    /// create interfacial friction heat. A late refusal accepts no component.
    pub fn step_with_planar_wall_velocity_field_and_viscous_heating_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfacePlanarWallViscousReport, Error> {
        let (velocities, origin) = self
            .flow
            .prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let prepared = self.prepare_viscous_heat(heat, gate)?;
        let flow = self.advance_wall_enthalpy(
            ThermalWallBoundary::Field {
                native: PlanarWallVelocityField3 {
                    wall: motion.wall,
                    moment_origin: origin,
                    tangential_velocity: &velocities,
                },
                controls,
                moment_origin_m: motion.moment_origin_m,
            },
            prepared.heated,
            prepared.external,
            prepared.viscous,
            prepared.absolute_heat,
            gate,
        )?;
        Ok(FreeSurfacePlanarWallViscousReport {
            flow,
            applied: prepared.applied,
        })
    }

    fn advance_wall_enthalpy(
        &mut self,
        boundary: ThermalWallBoundary<'_>,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        absolute_heat: f64,
        gate: &CancelGate,
    ) -> Result<FreeSurfacePlanarWallEnthalpyReport, Error> {
        let (flow, wall) = self.advance_prepared_enthalpy_with_wall(
            heated,
            external,
            viscous,
            0.0,
            absolute_heat,
            None,
            None,
            None,
            None,
            Some(boundary),
            gate,
        )?;
        Ok(FreeSurfacePlanarWallEnthalpyReport {
            flow,
            wall: match wall {
                Some(ThermalWallReport::Field(report)) => report,
                _ => unreachable!("field boundary retains its original wall receipt"),
            },
        })
    }
}
