//! Finite-temperature noncondensable inventory on the original DG trajectory.
//!
//! Each original ClosedGasParcel retains its mass, composition, NASA-9 calorics
//! and source axes. One additional conservative coordinate is Delta U from its
//! accepted state. The shared liquid/pocket solve uses p_g=m R T(U)/V and
//! U'=heat_power-p_g Q, with the original temperature-dependent Cv. Consequently
//! K_g=gamma(T) p_g/V and heat adds R/(Cv V)*heat_power to the pressure-rate
//! constraint. No endpoint pressure projection or frozen-gamma refresh is used.
//!
//! Liquid material and its saturation pressure remain frozen at their original
//! temperature. The gas may have a different uniform temperature. Vapor mass,
//! latent heat, phase kinetics, dissolution and gas transport are not supplied.
//! Prescribed heat is an external power port, including a caller-owned thermal
//! resistance law; its exact trajectory integral is included in the balance.

use super::*;
use fs_material::thermochemical_gas::energy::increment::GasCaloricBranch;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// One source-backed finite gas, replacing the polytropic law at this port.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricGasPocket {
    pub port: PortId,
    pub gas: ClosedGasParcel,
    pub minimum_volume_m3: f64,
    /// Original continuous temperature/source branch and bounded caloric solve.
    pub solve: GasEnergySolve,
}

/// Instantaneous original-caloric state supplied to a deterministic heat law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricGasObservation {
    pub port: PortId,
    pub volume_m3: f64,
    pub mass_kg: f64,
    pub temperature_k: f64,
    pub gas_pressure_pa: f64,
    pub cv_j_kg_k: f64,
}

/// Independent accepted caloric, heat-port and native pressure-work accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricGasStep {
    pub port: PortId,
    pub mass_kg: f64,
    pub temperature_k: [f64; 2],
    pub gas_pressure_pa: [f64; 2],
    pub internal_energy_change_j: f64,
    pub supplied_heat_j: f64,
    /// Work from ORIGINAL DG port work plus (p_reference-p_sat) Delta V.
    pub work_by_gas_j: f64,
    pub energy_balance_residual_j: f64,
    /// Actual original-law Delta U minus the integrated conservative coordinate.
    pub source_energy_residual_j: f64,
}

pub(super) type HeatLaw<'a> =
    dyn Fn(f64, &CaloricGasObservation) -> Result<f64, DuctCavityError> + 'a;

#[derive(Clone, Copy, Default)]
pub(super) struct CaloricStage {
    pub pressure: [f64; MAX_DUCT_POCKETS],
    pub stiffness: [f64; MAX_DUCT_POCKETS],
    pub forcing: [f64; MAX_DUCT_POCKETS],
    pub heat: [f64; MAX_DUCT_POCKETS],
}

pub(super) struct CaloricEndpoint {
    pub index: usize,
    pub gas: ClosedGasParcel,
    pub step: CaloricGasStep,
}

pub(super) fn error(error: GasEnergyError) -> DuctCavityError {
    if error == GasEnergyError::Cancelled {
        DuctCavityError::Cancelled
    } else {
        DuctCavityError::Gas(error)
    }
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuctPockets<F> {
    /// The original accepted finite gas, including source model and query axes.
    pub fn caloric_gas(&self, port: PortId) -> Option<&ClosedGasParcel> {
        self.sites
            .iter()
            .find(|s| s.specification.port == port)
            .and_then(|s| s.caloric.as_ref())
            .map(|p| &p.gas)
    }

    /// Advance the SAME acoustic/pocket/gas trajectory with signed heat power.
    /// The callback returns watts INTO one gas and must be deterministic in
    /// time and the supplied state; rejected RK trials can revisit times.
    /// Resolve drive frequencies with maximum_step_s. Callback side effects
    /// are not transactional; all accepted physical states are transactional.
    pub fn advance_with_gas_heat<Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError>
    where
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, DuctCavityError>,
        Cancel: FnMut() -> bool,
    {
        self.advance_kernel(end_time_s, numerics, heat, cancelled)
    }

    pub(super) fn prepare_caloric(
        &self,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<[Option<GasCaloricBranch<'_>>; MAX_DUCT_POCKETS], DuctCavityError> {
        let mut result = [None; MAX_DUCT_POCKETS];
        for (i, site) in self.sites.iter().enumerate() {
            if let Some(pocket) = &site.caloric {
                result[i] = Some(
                    pocket
                        .gas
                        .model()
                        .prepare_caloric_branch(pocket.solve, &mut || cancelled())
                        .map_err(error)?,
                );
            }
        }
        Ok(result)
    }

    pub(super) fn caloric_stage(
        &self,
        branches: &[Option<GasCaloricBranch<'_>>; MAX_DUCT_POCKETS],
        time_s: f64,
        volumes: &[f64],
        energy_changes: &[f64],
        heat: &HeatLaw<'_>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CaloricStage, DuctCavityError> {
        let mut stage = CaloricStage::default();
        let bulk = self.duct.acoustic_material().adiabatic_bulk_modulus_pa();
        for (i, site) in self.sites.iter().enumerate() {
            let Some(pocket) = &site.caloric else {
                continue;
            };
            if cancelled() {
                return Err(DuctCavityError::Cancelled);
            }
            let volume = volumes[i];
            if !volume.is_finite()
                || volume < pocket.minimum_volume_m3
                || volume > site.specification.maximum_volume_m3
            {
                return Err(DuctCavityError::OutsideDomain(
                    "caloric gas reaches its admitted volume boundary",
                ));
            }
            let density = pocket.gas.mass_kg() / volume;
            let increment = energy_changes[i] / pocket.gas.mass_kg();
            if !positive(density)
                || !increment.is_finite()
                || (energy_changes[i] != 0.0 && increment == 0.0)
            {
                return Err(DuctCavityError::OutsideDomain(
                    "caloric gas density or energy coordinate is unrepresentable",
                ));
            }
            let branch = branches[i].as_ref().ok_or(DuctCavityError::InvalidInput(
                "missing original caloric branch",
            ))?;
            let sample = branch
                .sample_density_internal_energy_increment(
                    pocket.gas.point(),
                    density,
                    increment,
                    &mut || cancelled(),
                )
                .map_err(error)?;
            let gas_constant = pocket.gas.state().state().specific_gas_constant;
            let ratio = gas_constant / sample.cv_j_kg_k;
            let gamma = 1.0 + ratio;
            let pressure = sample.pressure_pa;
            let stiffness = gamma * (pressure / volume);
            let absolute = self.material.saturation_pressure_pa() + pressure;
            if !positive(ratio)
                || gamma <= 1.0
                || gamma > 5.0 / 3.0 + 1e-9
                || !positive(stiffness)
                || !absolute.is_finite()
                || ((absolute - self.limits.reference_pressure_pa) / bulk).abs()
                    > self.limits.maximum_liquid_strain
            {
                return Err(DuctCavityError::OutsideDomain(
                    "caloric gas pressure or native gamma/strain domain exceeded",
                ));
            }
            let observation = CaloricGasObservation {
                port: site.specification.port,
                volume_m3: volume,
                mass_kg: pocket.gas.mass_kg(),
                temperature_k: sample.temperature_k,
                gas_pressure_pa: pressure,
                cv_j_kg_k: sample.cv_j_kg_k,
            };
            let power = heat(time_s, &observation)?;
            let forcing = (ratio / volume) * power;
            if !power.is_finite() || !forcing.is_finite() || (power != 0.0 && forcing == 0.0) {
                return Err(DuctCavityError::OutsideDomain(
                    "caloric gas heat power or pressure-rate forcing is unrepresentable",
                ));
            }
            stage.pressure[i] = pressure;
            stage.stiffness[i] = stiffness;
            stage.forcing[i] = forcing;
            stage.heat[i] = power;
        }
        Ok(stage)
    }

    pub(super) fn pocket_gas_pressure(
        &self,
        index: usize,
        volume: f64,
        stage: Option<&CaloricStage>,
    ) -> Result<f64, DuctCavityError> {
        if self.sites[index].caloric.is_some() {
            stage
                .map(|s| s.pressure[index])
                .ok_or(DuctCavityError::InvalidInput(
                    "missing caloric pressure state",
                ))
        } else {
            Ok(self.gas_response(index, volume)?.pressure_pa)
        }
    }

    pub(super) fn recover_caloric(
        &self,
        volumes: &[f64],
        energies: &[f64],
        heats: &[f64],
        receipts: &[DgTrajectoryPort],
        numerics: &DuctCavityNumerics,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<CaloricEndpoint>, DuctCavityError> {
        let mut endpoints = Vec::new();
        for (i, site) in self.sites.iter().enumerate() {
            let Some(pocket) = &site.caloric else {
                continue;
            };
            let (gas, actual) = pocket
                .gas
                .at_volume_internal_energy_increment(
                    volumes[i],
                    energies[i],
                    pocket.solve,
                    &mut || cancelled(),
                )
                .map_err(error)?;
            let work = receipts[i].work_into_duct_j
                + (self.limits.reference_pressure_pa - self.material.saturation_pressure_pa())
                    * (volumes[i] - pocket.gas.volume_m3());
            let residual = actual + work - heats[i];
            let tolerance = numerics.absolute_energy_tolerance_j
                + numerics.relative_tolerance * (actual.abs() + work.abs() + heats[i].abs());
            if !work.is_finite() || !residual.is_finite() || residual.abs() > tolerance {
                return Err(DuctCavityError::EnergyBalance);
            }
            let step = CaloricGasStep {
                port: site.specification.port,
                mass_kg: gas.mass_kg(),
                temperature_k: [
                    pocket.gas.state().state().temperature,
                    gas.state().state().temperature,
                ],
                gas_pressure_pa: [
                    pocket.gas.state().state().pressure,
                    gas.state().state().pressure,
                ],
                internal_energy_change_j: actual,
                supplied_heat_j: heats[i],
                work_by_gas_j: work,
                energy_balance_residual_j: residual,
                source_energy_residual_j: actual - energies[i],
            };
            endpoints.push(CaloricEndpoint {
                index: i,
                gas,
                step,
            });
        }
        Ok(endpoints)
    }
}

#[cfg(test)]
mod tests;
