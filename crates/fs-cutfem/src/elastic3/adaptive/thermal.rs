//! Native thermal source recovery through the original adaptive constraint map.

use std::ops::ControlFlow;

use fs_material::elastic::thermal::IsothermalThermoelasticState;
use fs_solver::op::LinearOp;

use super::AdaptiveElasticity3;
use crate::elastic3::ElasticityError3;
use crate::elastic3::eigenstrain::PreparedEigenstrain3;
use crate::elastic3::thermal::ThermoelasticResponse3;

mod entropy;
mod equilibrium;
mod exchange;
mod jacobian;
mod sensitivity;
mod stress;
mod state;

impl AdaptiveElasticity3 {
    /// Complete cell-constant thermal source in the original physical leaf.
    #[must_use]
    pub fn cell_thermoelastic_state(&self, cell: usize) -> Option<&IsothermalThermoelasticState> {
        self.raw.cell_thermoelastic_state(cell)
    }

    /// Prepare original integrated free strain on the retained cut geometry.
    /// Existing adaptive loads apply the original hanging-node transpose.
    pub fn prepare_thermoelastic(
        &self,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedEigenstrain3, ElasticityError3> {
        self.raw.prepare_thermoelastic(max_points, poll)
    }

    /// Recover the native potential and mechanical response from independent
    /// coordinates, reconstructing only one cell at a time through the SAME
    /// constraint map and homogeneous strong-support mask as the operator.
    pub fn thermoelastic_response(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticResponse3>, ElasticityError3> {
        self.raw.thermoelastic_response_mapped::<true>(
            u,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |cell, field| {
                std::array::from_fn(|i| {
                    self.physical_coefficient(field, cell.nodes[i / 3], i % 3)
                })
            },
            poll,
        )
    }

    /// Allocation-free native total/elastic domain and potential admission.
    pub fn admit_thermoelastic(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.raw
            .thermoelastic_response_mapped::<false>(
                u,
                self.n(),
                maximum_kinematic_strain,
                max_points,
                |cell, field| {
                    std::array::from_fn(|i| {
                        self.physical_coefficient(field, cell.nodes[i / 3], i % 3)
                    })
                },
                poll,
            )
            .map(|_| ())
    }
}
