//! Fixed noncondensable inventory on the shared liquid-pocket trajectory.
//!
//! The original BubbleGas owns p_g(V). This coupling uses its barotropic work
//! potential, never labels it latent heat or general gas caloric energy.
//! A positive inventory is always an open gas/vapor pocket: a declared positive
//! minimum volume is a refusal boundary, not a bounce or an artificial closure.
//!
//! At every gas pocket p_abs = p_sat + p_g(V) and dV/dt = Q. Differentiating
//! this nonlinear pressure constraint adds -dp_g/dV to the original coupled
//! liquid Gram matrix. The one liquid/gas trajectory exchanges work through
//! the original pressure/flow ports. Total mechanical storage includes
//! (p_reference-p_sat)V - integral(V0..V) p_g(v)dv.
//!
//! The existing VaporCavity mode denotes an open pocket, including this gas
//! and saturated-vapor mixture. Only zero-inventory pure vapor can collapse.
//! Gas mass transport, dissolution, heat transport and phase kinetics are
//! outside this fixed-inventory barotropic law.

use super::*;
use fs_flux::cavitation::BubbleGas;
use fs_flux::cavitation::gas_volume::GasVolumeResponse;
use fs_la::LuWorkspace;

/// Optional gas inventory attached to one already declared vapor-pocket port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrappedGasPocket {
    pub port: PortId,
    pub gas: BubbleGas,
    /// Positive lower volume bound for POSITIVE gas inventory [m³].
    /// An exactly zero original gas inventory retains the pure-vapor path.
    pub minimum_volume_m3: f64,
}
impl TrappedGasPocket {
    pub(super) fn has_inventory(self) -> bool {
        self.gas.reference_pressure_pa() > 0.0
    }
}

/// One bounded allocation at interval entry; every dynamic factor/solve reuses
/// these buffers. Native cancellation is polled around each RHS evaluation.
pub(super) struct GasFlowWorkspace {
    solver: LuWorkspace,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    solution: Vec<f64>,
}
impl GasFlowWorkspace {
    pub(super) fn new(count: usize) -> Result<Self, DuctCavityError> {
        Ok(Self {
            solver: LuWorkspace::new(count).map_err(|_| {
                DuctCavityError::ActiveConstraints("gas pressure workspace allocation")
            })?,
            matrix: vec![0.0; count * count],
            rhs: vec![0.0; count],
            solution: vec![0.0; count],
        })
    }
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuctPockets<F> {
    /// Original fixed gas law, including an explicitly authored zero inventory.
    pub fn trapped_gas(&self, port: PortId) -> Option<TrappedGasPocket> {
        self.sites
            .iter()
            .find(|s| s.specification.port == port)
            .and_then(|s| s.gas)
    }

    pub(super) fn has_inventory(&self, index: usize) -> bool {
        self.sites[index].caloric.is_some()
            || self.sites[index]
                .gas
                .is_some_and(TrappedGasPocket::has_inventory)
    }

    pub(super) fn gas_response(
        &self,
        index: usize,
        volume: f64,
    ) -> Result<GasVolumeResponse, DuctCavityError> {
        let Some(gas) = self.sites[index].gas.filter(|g| g.has_inventory()) else {
            return Ok(GasVolumeResponse {
                pressure_pa: 0.0,
                stiffness_pa_per_m3: 0.0,
                mechanical_potential_j: 0.0,
            });
        };
        if !volume.is_finite() || volume < gas.minimum_volume_m3 {
            return Err(DuctCavityError::OutsideDomain(
                "positive gas inventory reaches its minimum admitted volume",
            ));
        }
        gas.gas.volume_response(volume).map_err(|_| {
            DuctCavityError::OutsideDomain("original gas volume law is unrepresentable")
        })
    }

    pub(super) fn pocket_potential(
        &self,
        index: usize,
        volume: f64,
    ) -> Result<f64, DuctCavityError> {
        let value = (self.limits.reference_pressure_pa - self.material.saturation_pressure_pa())
            * volume
            + self.gas_response(index, volume)?.mechanical_potential_j;
        if !value.is_finite() {
            return Err(DuctCavityError::OutsideDomain(
                "combined cavity mechanical potential overflows",
            ));
        }
        Ok(value)
    }

    pub(super) fn gas_stiffness(
        &self,
        volumes: &[f64],
    ) -> Result<[f64; MAX_DUCT_POCKETS], DuctCavityError> {
        let mut result = [0.0; MAX_DUCT_POCKETS];
        let bulk = self.duct.acoustic_material().adiabatic_bulk_modulus_pa();
        for (i, value) in result[..self.sites.len()].iter_mut().enumerate() {
            if self.sites[i].caloric.is_some() || !self.has_inventory(i) {
                continue;
            }
            let response = self.gas_response(i, volumes[i])?;
            let absolute = self.material.saturation_pressure_pa() + response.pressure_pa;
            if !absolute.is_finite()
                || ((absolute - self.limits.reference_pressure_pa) / bulk).abs()
                    > self.limits.maximum_liquid_strain
            {
                return Err(DuctCavityError::OutsideDomain(
                    "gas-cushion pressure exceeds the native liquid strain domain",
                ));
            }
            *value = response.stiffness_pa_per_m3;
        }
        Ok(result)
    }

    pub(super) fn gas_flows(
        &self,
        active_mask: u16,
        x: &[f64],
        scale: f64,
        volumes: &[f64],
        thermal: Option<&CaloricStage>,
        work: &mut GasFlowWorkspace,
    ) -> Result<[f64; MAX_DUCT_POCKETS], DuctCavityError> {
        let count = self.sites.len();
        let mut stiffness = self.gas_stiffness(volumes)?;
        if let Some(stage) = thermal {
            for (i, value) in stiffness[..count].iter_mut().enumerate() {
                if self.sites[i].caloric.is_some() {
                    *value = stage.stiffness[i];
                }
            }
        }
        work.matrix.fill(0.0);
        work.rhs.fill(0.0);
        for i in 0..count {
            if !active(active_mask, i) {
                if self.has_inventory(i) {
                    return Err(DuctCavityError::ActiveConstraints(
                        "positive gas inventory cannot close",
                    ));
                }
                work.matrix[i * count + i] = 1.0;
                continue;
            }
            for j in 0..count {
                if active(active_mask, j) {
                    work.matrix[i * count + j] =
                        self.gram[i * count + j] / self.column_norms[i] / self.column_norms[j];
                }
            }
            work.matrix[i * count + i] +=
                stiffness[i] / self.column_norms[i] / self.column_norms[i];
            work.rhs[i] = (-scale * dot(&self.pressure_rates[i * x.len()..(i + 1) * x.len()], x)
                + thermal.map_or(0.0, |stage| stage.forcing[i]))
                / self.column_norms[i];
        }
        work.solver
            .solve_into(&work.matrix, &work.rhs, &mut work.solution)
            .map_err(|_| {
                DuctCavityError::ActiveConstraints(
                    "coupled gas/liquid pressure solve is unrepresentable",
                )
            })?;
        let mut result = [0.0; MAX_DUCT_POCKETS];
        for i in 0..count {
            if active(active_mask, i) {
                result[i] = work.solution[i] / self.column_norms[i];
            }
            if !result[i].is_finite() {
                return Err(DuctCavityError::ActiveConstraints(
                    "gas/liquid flow overflows",
                ));
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
