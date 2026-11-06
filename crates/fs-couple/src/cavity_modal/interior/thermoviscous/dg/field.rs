//! Native prepared-field dispatch for the one reciprocal FE/DG solver.
//! The sealed implementations preserve each original field owner's state,
//! mean-port map, source-work integration and publication rules.

use super::*;
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousError, ThermoviscousFe,
    time::{
        PreparedThermoviscousVolumeTime, TangentialWallConductivity, ThermoviscousCaloricReservoir,
    },
};

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::PreparedThermoviscousTime<'_, '_> {}
    impl Sealed for super::PreparedThermoviscousVolumeTime<'_, '_> {}
}

/// Original declarations used only to admit the unchanged physical connection.
#[doc(hidden)]
#[derive(Clone, Copy)]
pub struct ThermoviscousDgFieldBinding<'s, 'c> {
    pub model: &'s ThermoviscousFe<'c>,
    pub duration_s: f64,
    pub frozen_cp_j_kg_k: f64,
    pub reservoirs: &'s [ThermoviscousCaloricReservoir],
    pub conductive_wall_count: usize,
    pub maximum_wall_frequency_hz: Option<f64>,
    pub tangential_wall_conductivity: Option<&'s [TangentialWallConductivity]>,
    pub body_force_count: usize,
    pub heat_count: usize,
}

/// Dispatch restricted to the two original complete prepared field owners.
/// This is not an extension point for replacement field or source laws.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub trait ThermoviscousDgFields<'c>: sealed::Sealed {
    fn binding(&self) -> ThermoviscousDgFieldBinding<'_, 'c>;
    fn zero_state(&self) -> ThermoviscousTimeState;
    fn stored_energy_j(
        &self,
        state: &ThermoviscousTimeState,
    ) -> Result<f64, ThermoviscousTimeError>;
    fn average_relation(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        port: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ThermoviscousTimeError>;
    fn average_relations(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        ports: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError>;
    fn step_fields(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            ThermoviscousTimeStep,
            Option<ThermoviscousDgVolumeObservations>,
        ),
        ThermoviscousTimeError,
    >;
}

/// Original interval observations of the volume inputs on the last accepted
/// coupled frame. Work is signed and included once in that frame's
/// authored_external_work_j. The source shapes are never normalized here.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoviscousDgVolumeObservations {
    /// Mean integral(shape_s . v dV) [m/s], including moving-wall velocity.
    pub mean_body_force_velocity_m_s: Vec<f64>,
    /// Signed mechanical source work [J], in the original source order.
    pub body_force_work_by_source_j: Vec<f64>,
    /// Mean integral(shape_s * T'/T0 dV), dimensionless.
    pub mean_heat_source_temperature_ratio: Vec<f64>,
    /// Signed acoustic available-energy work from heating [J].
    pub heat_source_work_by_source_j: Vec<f64>,
}

fn passive_sources(force: &[f64], heat: &[f64]) -> Result<(), ThermoviscousTimeError> {
    if force.is_empty() && heat.is_empty() {
        Ok(())
    } else {
        Err(ThermoviscousError::InvalidInput(
            "passive prepared fields have no authored volume-source columns",
        )
        .into())
    }
}

#[allow(clippy::too_many_arguments)]
impl<'c> ThermoviscousDgFields<'c> for PreparedThermoviscousTime<'_, 'c> {
    fn binding(&self) -> ThermoviscousDgFieldBinding<'_, 'c> {
        ThermoviscousDgFieldBinding {
            model: self.inner(),
            duration_s: self.duration_s(),
            frozen_cp_j_kg_k: self.frozen_cp_j_kg_k(),
            reservoirs: self.caloric_reservoirs(),
            conductive_wall_count: self.conductive_wall_layouts().len(),
            maximum_wall_frequency_hz: self.maximum_wall_frequency_hz(),
            tangential_wall_conductivity: self.tangential_wall_conductivity(),
            body_force_count: 0,
            heat_count: 0,
        }
    }

    fn zero_state(&self) -> ThermoviscousTimeState {
        PreparedThermoviscousTime::zero_state(self)
    }

    fn stored_energy_j(
        &self,
        state: &ThermoviscousTimeState,
    ) -> Result<f64, ThermoviscousTimeError> {
        PreparedThermoviscousTime::stored_energy_j(self, state)
    }

    fn average_relation(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        port: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ThermoviscousTimeError> {
        passive_sources(force, heat)?;
        PreparedThermoviscousTime::average_relation(self, initial, traction, wall, port, cancelled)
    }

    fn average_relations(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        ports: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError> {
        passive_sources(force, heat)?;
        PreparedThermoviscousTime::average_relations(
            self, initial, traction, wall, ports, cancelled,
        )
    }

    fn step_fields(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            ThermoviscousTimeStep,
            Option<ThermoviscousDgVolumeObservations>,
        ),
        ThermoviscousTimeError,
    > {
        passive_sources(force, heat)?;
        Ok((self.step(initial, traction, wall, cancelled)?, None))
    }
}

#[allow(clippy::too_many_arguments)]
impl<'c> ThermoviscousDgFields<'c> for PreparedThermoviscousVolumeTime<'_, 'c> {
    fn binding(&self) -> ThermoviscousDgFieldBinding<'_, 'c> {
        ThermoviscousDgFieldBinding {
            model: self.inner(),
            duration_s: self.duration_s(),
            frozen_cp_j_kg_k: self.frozen_cp_j_kg_k(),
            reservoirs: self.caloric_reservoirs(),
            conductive_wall_count: self.conductive_wall_layouts().len(),
            maximum_wall_frequency_hz: self.maximum_wall_frequency_hz(),
            tangential_wall_conductivity: self.tangential_wall_conductivity(),
            body_force_count: self.body_force_source_count(),
            heat_count: self.heat_source_count(),
        }
    }

    fn zero_state(&self) -> ThermoviscousTimeState {
        PreparedThermoviscousVolumeTime::zero_state(self)
    }

    fn stored_energy_j(
        &self,
        state: &ThermoviscousTimeState,
    ) -> Result<f64, ThermoviscousTimeError> {
        PreparedThermoviscousVolumeTime::stored_energy_j(self, state)
    }

    fn average_relation(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        port: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ThermoviscousTimeError> {
        PreparedThermoviscousVolumeTime::average_relation(
            self, initial, traction, wall, force, heat, port, cancelled,
        )
    }

    fn average_relations(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        ports: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError> {
        PreparedThermoviscousVolumeTime::average_relations(
            self, initial, traction, wall, force, heat, ports, cancelled,
        )
    }

    fn step_fields(
        &self,
        initial: &ThermoviscousTimeState,
        traction: &[f64],
        wall: &[f64],
        force: &[f64],
        heat: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            ThermoviscousTimeStep,
            Option<ThermoviscousDgVolumeObservations>,
        ),
        ThermoviscousTimeError,
    > {
        let candidate = self.step(initial, traction, wall, force, heat, cancelled)?;
        Ok((
            candidate.field,
            Some(ThermoviscousDgVolumeObservations {
                mean_body_force_velocity_m_s: candidate.mean_body_force_velocity_m_s,
                body_force_work_by_source_j: candidate.body_force_work_by_source_j,
                mean_heat_source_temperature_ratio: candidate.mean_heat_source_temperature_ratio,
                heat_source_work_by_source_j: candidate.heat_source_work_by_source_j,
            }),
        ))
    }
}
