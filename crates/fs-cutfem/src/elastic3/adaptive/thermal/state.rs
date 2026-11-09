//! Instantaneous thermal states through the original hanging-node map.
use std::ops::ControlFlow;

use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::eigenstrain::PreparedEigenstrain3;
use crate::elastic3::mass::MaterialMass3;
use crate::elastic3::thermal::equilibrium::EquilibriumOwner3;
use crate::elastic3::thermal::state::recover_state;
use crate::elastic3::thermal::{
    ThermalCoordinateMap3, ThermoelasticState3, ThermoelasticStateSettings3,
};
use crate::elastic3::ElasticityError3;

impl AdaptiveElasticity3 {
    /// Admit a physical instantaneous displacement on original independent
    /// coordinates. Reconstructs cell strains and scatters full resistance with
    /// the SAME hanging map and transpose; no static force gate is imposed.
    pub fn prepare_thermoelastic_state(
        &self,
        u: &[f64],
        settings: ThermoelasticStateSettings3,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticState3, ElasticityError3> {
        recover_state(EquilibriumOwner3::Adaptive(self),
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            u, None, settings,
            |node, component, value, output|
                self.scatter_full_force(node, component, value, output),
            poll)
    }

    /// Consume original preparations while recovering a new physical endpoint.
    /// Source/geometry/material bindings and retained mass limits are checked.
    /// No geometry, hanging-map, mass or free-strain force reassembly occurs.
    pub fn recover_thermoelastic_state(
        &self,
        u: &[f64],
        eigenstrain: PreparedEigenstrain3,
        mass: MaterialMass3,
        settings: ThermoelasticStateSettings3,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticState3, ElasticityError3> {
        recover_state(EquilibriumOwner3::Adaptive(self),
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            u, Some((eigenstrain, mass)), settings,
            |node, component, value, output|
                self.scatter_full_force(node, component, value, output),
            poll)
    }
}
