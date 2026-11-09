//! Finite physical exchange through the original adaptive constraint map.

use std::ops::ControlFlow;

use fs_material::elastic::thermal::AffineThermoelastic;
use fs_solver::op::LinearOp;

use crate::HexCell;
use crate::elastic3::ElasticityError3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::thermal::{ThermoelasticCellExchange3, ThermoelasticExchangeLimits3};

impl AdaptiveElasticity3 {
    /// Original physical leaf bounds, with no new geometry evaluation.
    #[must_use]
    pub fn reference_cell_bounds(&self, cell: usize) -> Option<HexCell> {
        self.raw.reference_cell_bounds(cell)
    }

    /// Physical cut volume on the retained original reference rule.
    #[must_use]
    pub fn reference_cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.raw.reference_cell_volume_m3(cell)
    }

    #[must_use]
    pub fn has_embedded_supports(&self) -> bool {
        self.raw.has_embedded_supports()
    }

    /// Source-owned finite exchange on the SAME adaptive independent space.
    ///
    /// Reconstructs only two cell-local fields through the existing hanging
    /// rows and homogeneous strong-support mask. No global physical field or
    /// second constraint owner is allocated. The original laws, geometry and
    /// cell scales must remain fixed over the interval. Current native sources
    /// are the AFTER endpoint; all bulk/domain and budget semantics are those
    /// of CutElasticity3::thermoelastic_exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn thermoelastic_exchange(
        &self,
        laws: &[AffineThermoelastic],
        before_temperatures_k: &[f64],
        before_u: &[f64],
        after_u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermoelasticExchangeLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticCellExchange3>, ElasticityError3> {
        self.raw.thermoelastic_exchange_mapped(
            laws, before_temperatures_k, before_u, after_u, self.n(),
            maximum_kinematic_strain, limits,
            |cell, field| std::array::from_fn(|i| {
                self.physical_coefficient(field, cell.nodes[i / 3], i % 3)
            }),
            checkpoint,
        )
    }
}
