//! Native stress-objective derivatives at fixed reference geometry.
//!
//! The source potential supplies the mixed stress/temperature partial, including
//! both C'(T) and the derivative of integrated, possibly anisotropic expansion.
//! The original reference stress transpose owns the displacement derivative.

use std::ops::ControlFlow;

use fs_material::tensor::from_mandel;
use fs_solver::op::LinearOp;

use super::checkpoint;
use crate::elastic3::eigenstrain::physical_response;
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3};

/// Explicit partial derivatives of a scalar objective of native thermal stress.
/// These hold geometry fixed and do not differentiate the solved equilibrium.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoelasticStressPullback3 {
    /// `dJ/du` in the operator's independent displacement coordinates.
    pub displacement: Vec<f64>,
    /// `dJ/dscale_cell` at fixed displacement and source temperature.
    pub scales: Vec<f64>,
    /// `dJ/dT_cell` at fixed displacement and stiffness scale, per kelvin.
    /// Entries follow the same original active-cell order as `scales`.
    pub temperatures: Vec<f64>,
}

/// Only compact cotangents survive the streaming physical-response pass.
pub(in crate::elastic3) struct StressPartials3 {
    pub(in crate::elastic3) scaled_cotangents: Vec<[f64; 6]>,
    pub(in crate::elastic3) scales: Vec<f64>,
    pub(in crate::elastic3) temperatures: Vec<f64>,
}

impl CutElasticity3 {
    /// Exact stress-objective partials with respect to u, scales and temperature.
    ///
    /// Supply one cotangent `[xx, yy, zz, xy, yz, zx]` per ORIGINAL bulk point,
    /// pairing directly with physical stress components. Include quadrature
    /// weights in the cotangents when differentiating a volume integral; no
    /// weight or tensor-shear multiplicity is inserted here. Their exact count
    /// bounds all point storage before allocation. Cotangents are algebraic,
    /// finite quantities and do not consume a physical material-strain budget.
    ///
    /// The physical displacement passes the same native total/elastic strain
    /// admission as `thermoelastic_response`, including the original B-bar map.
    /// Each source supplies C(T), integrated expansion and its mixed partial at
    /// fixed constitutive TOTAL strain. The temperature derivative therefore
    /// includes C'(T) and anisotropic, temperature-dependent expansion, with the
    /// stiffness scale applied exactly once. Density and geometry stay fixed.
    ///
    /// For equilibrium `R(u,T)=0`, solve `K^T z = displacement` and subtract
    /// `thermal_residual_temperature_pullback(u,z,...)` from `temperatures`.
    /// Temperature-dependent external loads retain their own derivative terms.
    /// No derivative matrix or finite-difference material action is constructed.
    pub fn thermoelastic_stress_pullback(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        derivatives: &[[f64; 6]],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticStressPullback3, ElasticityError3> {
        let direct = self.thermoelastic_stress_partials_mapped(
            u,
            self.n(),
            maximum_kinematic_strain,
            derivatives,
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            &mut poll,
        )?;
        // Integrated free strain shifts stress but does not change dσ/du.
        // This owner retains the original tensor, B-bar and clamp transpose.
        let displacement =
            self.reference_bulk_stress_pullback(&direct.scaled_cotangents, &mut poll)?;
        checkpoint(&mut poll)?;
        Ok(ThermoelasticStressPullback3 {
            displacement,
            scales: direct.scales,
            temperatures: direct.temperatures,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::elastic3) fn thermoelastic_stress_partials_mapped(
        &self,
        u: &[f64],
        n: usize,
        maximum_kinematic_strain: f64,
        derivatives: &[[f64; 6]],
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<StressPartials3, ElasticityError3> {
        checkpoint(&mut poll)?;
        if u.len() != n
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(ElasticityError3::Invalid(
                "invalid thermoelastic displacement or kinematic bound",
            ));
        }
        for chunk in u.chunks(256) {
            checkpoint(&mut poll)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "nonfinite thermoelastic displacement",
                ));
            }
        }
        let mut count = 0usize;
        for cell in &self.cells {
            checkpoint(&mut poll)?;
            if cell.material.thermal_state().is_none() {
                return Err(ElasticityError3::Invalid(
                    "thermoelastic stress pullback requires a native thermal source in every cell",
                ));
            }
            count = count
                .checked_add(cell.rules.bulk().len())
                .ok_or(ElasticityError3::Invalid(
                    "thermoelastic point count overflow",
                ))?;
        }
        if count != derivatives.len() {
            return Err(ElasticityError3::Invalid(
                "one thermoelastic stress cotangent per bulk point required",
            ));
        }
        for derivative in derivatives {
            checkpoint(&mut poll)?;
            if derivative.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "nonfinite thermoelastic stress cotangent",
                ));
            }
        }
        let mut direct = StressPartials3 {
            scaled_cotangents: Vec::with_capacity(count),
            scales: vec![0.0; self.cells.len()],
            temperatures: vec![0.0; self.cells.len()],
        };
        for (id, cell) in self.cells.iter().enumerate() {
            checkpoint(&mut poll)?;
            let source = cell
                .material
                .thermal_state()
                .expect("thermal sources checked before allocation")
                .source();
            let free_mandel = source.free_strain_mandel();
            let free_strain = from_mandel(&free_mandel);
            let mean_trace = free_mandel[..3].iter().sum::<f64>();
            let local = local(cell, u);
            if local.iter().any(|value| !value.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "thermoelastic reconstructed displacement overflow",
                ));
            }
            let scale = self.scales[id];
            for &(position, weight) in cell.rules.bulk() {
                checkpoint(&mut poll)?;
                let mechanical = physical_response(
                    cell,
                    &local,
                    position,
                    free_strain,
                    mean_trace,
                    scale,
                    maximum_kinematic_strain,
                )?;
                let potential = source
                    .state(mechanical.constitutive_total_mandel, weight)
                    .map_err(|_| {
                        ElasticityError3::Invalid(
                            "native thermoelastic potential is outside its admitted domain",
                        )
                    })?;
                let beta = from_mandel(&potential.stress_temperature_derivative_mandel_pa_k)
                    .map(|value| scale * value);
                if [
                    potential.free_energy_j,
                    potential.entropy_j_k,
                    potential.internal_energy_j,
                    potential.heat_capacity_j_k,
                ]
                .iter()
                .any(|value| !(scale * value).is_finite())
                    || beta.iter().any(|value| !value.is_finite())
                {
                    return Err(ElasticityError3::Invalid(
                        "scaled thermoelastic potential overflow",
                    ));
                }
                let derivative = derivatives[direct.scaled_cotangents.len()];
                direct.scales[id] += derivative
                    .iter()
                    .zip(mechanical.reference_stress)
                    .map(|(d, stress)| d * stress)
                    .sum::<f64>();
                direct.temperatures[id] += derivative
                    .iter()
                    .zip(beta)
                    .map(|(d, slope)| d * slope)
                    .sum::<f64>();
                let scaled = derivative.map(|value| scale * value);
                if !direct.scales[id].is_finite()
                    || !direct.temperatures[id].is_finite()
                    || scaled.iter().any(|value| !value.is_finite())
                {
                    return Err(ElasticityError3::Invalid(
                        "thermoelastic stress pullback overflow",
                    ));
                }
                direct.scaled_cotangents.push(scaled);
            }
        }
        checkpoint(&mut poll)?;
        Ok(direct)
    }
}
