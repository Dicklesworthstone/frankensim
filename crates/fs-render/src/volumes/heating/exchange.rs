//! Passive volumetric emission, reabsorption and black-bath exchange on the
//! original thermal clock. The immutable sampled optical operator is reciprocal;
//! signed Planck differences create one equal/opposite native pair transfer.
//! Uniform material/bath temperature is therefore exactly at equilibrium for
//! every sample, without subtracting independently noisy positive powers.
//!
//! This owner retains fixed optical coefficients over explicitly supplied
//! temperature domains. Thermal temperatures and latent storage remain native.
//! Spectral weights integrate only the declared source window. A sampled
//! operator's uncertainty is correlated throughout its entire trajectory;
//! conditional watt diagnostics are not final thermal-state error bars.

use super::{CheckedSum, OpticalHeatingError};
use crate::thermal_surface::{ThermalSurfaceError, planck_radiative_exchange_w_and_w_k};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::exchange::{
    ExchangeCoefficientEstimate, RadiativeExchangeOperator,
};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// The open index-one boundary sees either vacuum or one isotropic black bath.
/// Material temperatures remain strictly positive; vacuum does not extend a
/// material's caloric or coefficient domain to zero kelvin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RadiativeBath {
    Vacuum,
    Blackbody { temperature_k: f64 },
}

/// Explicit validity of the fixed optical snapshot for one original material
/// cell. Every positive-density cell requires a binding, including transparent
/// cells and pure scatterers that another thermal owner can heat.
/// Multiple optical cells can map to one native cell; their valid-temperature
/// intervals are intersected. The binding does not infer a caloric law or mass.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeThermalCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RadiativeExchangeBudget {
    pub max_bindings: usize,
    pub max_thermal_pairs: usize,
    /// Counts all source preparations over one complete native advance.
    pub max_source_evaluations: usize,
}

/// Explicit admission of finite-sample coefficient and heat diagnostics.
/// Errors combine conservatively over potentially correlated cells/pairs/lanes;
/// none is a confidence interval, a rare-path guarantee or optional-stopping
/// bound. Unresolved exchange is never assigned to an invented heat receiver.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeExchangeTolerance {
    pub coefficient_roundoff_m2_sr: f64,
    pub power_roundoff_w: f64,
    pub max_cell_standard_error_w: f64,
    /// Sum of receiving-cell error allowances; internal pairs count twice.
    pub max_total_standard_error_w: f64,
    /// Estimated L1 scale of omitted heat over every possible tail endpoint.
    pub max_unresolved_exchange_w: f64,
    pub max_unresolved_standard_error_w: f64,
}

/// Current signed exchange at one accepted start state. Positive bath power
/// heats the material; positive pair magnitude counts each unordered pair once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeExchangePower<const N: usize> {
    pub medium_identity: ContentHash,
    pub net_bath_to_material_w: [f64; N],
    pub internal_transfer_magnitude_w: [f64; N],
    pub thermal_external_w: f64,
    pub maximum_cell_standard_error_w: f64,
    pub total_standard_error_w: f64,
    pub unresolved_exchange_w: f64,
    pub unresolved_standard_error_w: f64,
}

/// Actual accepted exchange, not a second thermal store. Internal exchange is
/// carried by native HeatTransfer and contributes zero external heat. Only the
/// signed bath account contributes to the original external-heat audit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeExchangeEnergy<const N: usize> {
    pub duration_s: f64,
    pub net_bath_to_material_j: [f64; N],
    pub internal_transfer_magnitude_j: [f64; N],
    pub thermal_external_j: f64,
    /// Integrated estimated L1 tail scale, without claiming a true error bound.
    pub unresolved_exchange_scale_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeExchangeError {
    Invalid(&'static str),
    Budget,
    TemperatureDomain {
        thermal_cell: CellId,
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    Transport(VolumeTransportError),
    Planck(ThermalSurfaceError),
    Heating(OpticalHeatingError),
    Thermal(DepositionError),
}
impl From<VolumeTransportError> for RadiativeExchangeError {
    fn from(e: VolumeTransportError) -> Self {
        Self::Transport(e)
    }
}
impl From<ThermalSurfaceError> for RadiativeExchangeError {
    fn from(e: ThermalSurfaceError) -> Self {
        Self::Planck(e)
    }
}
impl From<OpticalHeatingError> for RadiativeExchangeError {
    fn from(e: OpticalHeatingError) -> Self {
        Self::Heating(e)
    }
}
impl From<DepositionError> for RadiativeExchangeError {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl core::fmt::Display for RadiativeExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "volumetric radiative exchange refused: {self:?}")
    }
}
impl core::error::Error for RadiativeExchangeError {}

#[derive(Clone, Debug)]
struct BoundCell<const N: usize> {
    id: CellId,
    minimum_k: f64,
    maximum_k: f64,
    can_absorb: bool,
    escape: ExchangeCoefficientEstimate<N>,
    unresolved: ExchangeCoefficientEstimate<N>,
}
#[derive(Clone, Debug)]
struct BoundPair<const N: usize> {
    a: usize,
    b: usize,
    coefficient: ExchangeCoefficientEstimate<N>,
}
#[derive(Debug)]
struct Model<const N: usize> {
    operator: Arc<RadiativeExchangeOperator<N>>,
    weights_nm: [f64; N],
    bath: RadiativeBath,
    cells: Vec<BoundCell<N>>,
    pairs: Vec<BoundPair<N>>,
    budget: RadiativeExchangeBudget,
    tolerance: RadiativeExchangeTolerance,
}
#[derive(Clone, Debug)]
struct Prepared<const N: usize> {
    at_s: f64,
    temperatures_k: Vec<f64>,
    sources_w: BTreeMap<CellId, f64>,
    transfers: Vec<HeatTransfer>,
    conductance_bounds_w_k: BTreeMap<CellId, f64>,
    cell_errors_w: BTreeMap<CellId, f64>,
    power: RadiativeExchangePower<N>,
}

/// A fixed reciprocal optical operator with evolving native thermal state.
/// Immutable geometry/material/sampling data are shared in Arc; every prepared
/// rate and energy account is candidate-local and rolls back with the thermal
/// advance. Combine with incident heating through the existing CoupledObservers.
#[derive(Clone, Debug)]
pub struct VolumeRadiativeExchange<const N: usize> {
    model: Arc<Model<N>>,
    prepared: Prepared<N>,
    time_s: f64,
    evaluations: usize,
    evaluation_limit: usize,
    sums: EnergySums<N>,
    energy: RadiativeExchangeEnergy<N>,
}

impl<const N: usize> VolumeRadiativeExchange<N> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        operator: RadiativeExchangeOperator<N>,
        wavelength_weights_nm: [f64; N],
        bath: RadiativeBath,
        bindings: &[RadiativeThermalCell],
        thermal: &DepositionSimulation<'_>,
        budget: RadiativeExchangeBudget,
        tolerance: RadiativeExchangeTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, RadiativeExchangeError> {
        checkpoint()?;
        if thermal.cells().is_empty() {
            return Err(invalid(
                "radiative exchange requires an active native thermal state",
            ));
        }
        if N == 0
            || budget.max_source_evaluations == 0
            || bindings.len() > budget.max_bindings
            || operator.material_cells.len() > budget.max_bindings
        {
            return Err(RadiativeExchangeError::Budget);
        }
        if !wavelength_weights_nm.iter().copied().all(positive)
            || ![
                tolerance.coefficient_roundoff_m2_sr,
                tolerance.power_roundoff_w,
                tolerance.max_cell_standard_error_w,
                tolerance.max_total_standard_error_w,
                tolerance.max_unresolved_exchange_w,
                tolerance.max_unresolved_standard_error_w,
            ]
            .into_iter()
            .all(nonnegative)
            || matches!(bath, RadiativeBath::Blackbody { temperature_k } if !positive(temperature_k))
        {
            return Err(invalid(
                "positive quadrature/bath temperature and finite nonnegative error allowances required",
            ));
        }
        transport_checkpoint(&mut checkpoint, |check| {
            operator.validate_with_checkpoint(check, tolerance.coefficient_roundoff_m2_sr)
        })?;
        let mut selected = BTreeMap::new();
        for binding in bindings {
            checkpoint()?;
            if !positive(binding.minimum_temperature_k)
                || !positive(binding.maximum_temperature_k)
                || binding.minimum_temperature_k > binding.maximum_temperature_k
            {
                return Err(invalid(
                    "fixed optical snapshot requires an explicit positive temperature domain",
                ));
            }
            if selected.insert(binding.grid_cell, *binding).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell).into());
            }
            if !thermal.cells().contains_key(&binding.thermal_cell) {
                return Err(OpticalHeatingError::InactiveCell(binding.thermal_cell).into());
            }
        }
        let mut grouped = BTreeMap::<CellId, CellSums<N>>::new();
        let mut receivers = BTreeMap::new();
        for &grid_cell in &operator.material_cells {
            checkpoint()?;
            let binding = selected
                .remove(&grid_cell)
                .ok_or(OpticalHeatingError::UnmappedCell(grid_cell))?;
            let sums = grouped.entry(binding.thermal_cell).or_insert_with(|| {
                CellSums::new(binding.minimum_temperature_k, binding.maximum_temperature_k)
            });
            sums.minimum_k = sums.minimum_k.max(binding.minimum_temperature_k);
            sums.maximum_k = sums.maximum_k.min(binding.maximum_temperature_k);
            if sums.minimum_k > sums.maximum_k {
                return Err(invalid(
                    "one thermal receiver has disjoint optical validity domains",
                ));
            }
            receivers.insert(grid_cell, binding.thermal_cell);
        }
        // Scattering/transparent materials still constrain the fixed snapshot's
        // validity, even though they have no sparse emissive row or heat rate.
        for cell in &operator.cells {
            checkpoint()?;
            let receiver = receivers
                .get(&cell.cell_index)
                .ok_or_else(|| invalid("emitter is absent from the original material-cell set"))?;
            let sums = grouped
                .get_mut(receiver)
                .expect("validated material receiver");
            sums.can_absorb = true;
            sums.escape.add(&cell.escape)?;
            sums.unresolved.add(&cell.unresolved)?;
        }
        if !selected.is_empty() {
            return Err(invalid(
                "radiative binding names a cell outside the positive-density material field",
            ));
        }
        let mut pair_sums = BTreeMap::<(CellId, CellId), EstimateSums<N>>::new();
        for pair in &operator.pairs {
            checkpoint()?;
            let a = *receivers
                .get(&pair.cell_a)
                .ok_or_else(|| invalid("missing pair receiver"))?;
            let b = *receivers
                .get(&pair.cell_b)
                .ok_or_else(|| invalid("missing pair receiver"))?;
            if a == b {
                continue;
            }
            let key = if a < b { (a, b) } else { (b, a) };
            if !pair_sums.contains_key(&key) && pair_sums.len() >= budget.max_thermal_pairs {
                return Err(RadiativeExchangeError::Budget);
            }
            pair_sums.entry(key).or_default().add(&pair.exchange)?;
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(grouped.len())
            .map_err(|_| RadiativeExchangeError::Budget)?;
        let mut indices = BTreeMap::new();
        for (id, sums) in grouped {
            checkpoint()?;
            indices.insert(id, cells.len());
            cells.push(BoundCell {
                id,
                minimum_k: sums.minimum_k,
                maximum_k: sums.maximum_k,
                can_absorb: sums.can_absorb,
                escape: sums.escape.value()?,
                unresolved: sums.unresolved.value()?,
            });
        }
        let mut pairs = Vec::new();
        pairs
            .try_reserve_exact(pair_sums.len())
            .map_err(|_| RadiativeExchangeError::Budget)?;
        for ((a, b), sums) in pair_sums {
            checkpoint()?;
            pairs.push(BoundPair {
                a: indices[&a],
                b: indices[&b],
                coefficient: sums.value()?,
            });
        }
        let model = Arc::new(Model {
            operator: Arc::new(operator),
            weights_nm: wavelength_weights_nm,
            bath,
            cells,
            pairs,
            budget,
            tolerance,
        });
        let prepared = model.prepare(thermal, &mut checkpoint)?;
        let sums = EnergySums::default();
        let energy = sums.value()?;
        Ok(Self {
            model,
            prepared,
            time_s: thermal.time_s(),
            evaluations: 0,
            evaluation_limit: budget.max_source_evaluations,
            sums,
            energy,
        })
    }

    pub fn operator(&self) -> &RadiativeExchangeOperator<N> {
        &self.model.operator
    }
    pub fn wavelength_weights_nm(&self) -> &[f64; N] {
        &self.model.weights_nm
    }
    pub fn bath(&self) -> RadiativeBath {
        self.model.bath
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn budget(&self) -> RadiativeExchangeBudget {
        RadiativeExchangeBudget {
            max_source_evaluations: self.evaluation_limit,
            ..self.model.budget
        }
    }
    /// Tighten the next complete native advance without cloning immutable data.
    pub fn with_work_limit(&self, evaluations: usize) -> Self {
        let mut result = self.clone();
        result.evaluation_limit = evaluations.min(self.evaluation_limit);
        result
    }
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.prepared.sources_w
    }
    pub fn transfers(&self) -> &[HeatTransfer] {
        &self.prepared.transfers
    }
    pub fn conductance_bounds_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.prepared.conductance_bounds_w_k
    }
    pub fn conditional_cell_standard_errors_w(&self) -> &BTreeMap<CellId, f64> {
        &self.prepared.cell_errors_w
    }
    /// Rates used over the last interval, or initial rates before an advance.
    pub fn power(&self) -> &RadiativeExchangePower<N> {
        &self.prepared.power
    }
    pub fn energy(&self) -> &RadiativeExchangeEnergy<N> {
        &self.energy
    }

    fn check_clock(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), RadiativeExchangeError> {
        if self.time_s != thermal.time_s() {
            return Err(invalid("radiation and thermal clocks differ"));
        }
        Ok(())
    }
    fn check_prepared(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), RadiativeExchangeError> {
        self.check_clock(thermal)?;
        if self.prepared.at_s != thermal.time_s() {
            return Err(invalid(
                "radiative source must be prepared before transfer or CFL queries",
            ));
        }
        for (cell, &temperature) in self.model.cells.iter().zip(&self.prepared.temperatures_k) {
            check()?;
            if thermal.temperature_k(cell.id)? != temperature {
                return Err(invalid(
                    "radiative hooks require the unchanged prepared thermal state",
                ));
            }
        }
        Ok(())
    }
}

impl<const N: usize> Model<N> {
    fn temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<f64>, RadiativeExchangeError> {
        let mut temperatures = Vec::new();
        temperatures
            .try_reserve_exact(self.cells.len())
            .map_err(|_| RadiativeExchangeError::Budget)?;
        for cell in &self.cells {
            check()?;
            let t = thermal.temperature_k(cell.id)?;
            if !t.is_finite() || !(cell.minimum_k..=cell.maximum_k).contains(&t) {
                return Err(RadiativeExchangeError::TemperatureDomain {
                    thermal_cell: cell.id,
                    temperature_k: t,
                    minimum_k: cell.minimum_k,
                    maximum_k: cell.maximum_k,
                });
            }
            temperatures.push(t);
        }
        Ok(temperatures)
    }

    fn prepare(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Prepared<N>, RadiativeExchangeError> {
        let temperatures = self.temperatures(thermal, check)?;
        let count = self.cells.len();
        let mut errors = allocate(count, CheckedSum::default())?;
        let mut conductances = allocate(count, CheckedSum::default())?;
        let mut bath_sources = BTreeMap::new();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(self.pairs.len())
            .map_err(|_| RadiativeExchangeError::Budget)?;
        let mut bath_by_lane = [CheckedSum::default(); N];
        let mut internal_by_lane = [CheckedSum::default(); N];
        for pair in &self.pairs {
            check()?;
            let ta = temperatures[pair.a];
            let tb = temperatures[pair.b];
            let mut watts = CheckedSum::default();
            for lane in 0..N {
                let lambda = self.operator.wavelengths_nm[lane];
                let weight = self.weights_nm[lane];
                let (value, g) = planck_radiative_exchange_w_and_w_k(
                    ta,
                    Some(tb),
                    lambda,
                    pair.coefficient.mean_m2_sr[lane],
                    weight,
                )?;
                watts.add(value)?;
                internal_by_lane[lane].add(value.abs())?;
                let error = planck_radiative_exchange_w_and_w_k(
                    ta,
                    Some(tb),
                    lambda,
                    pair.coefficient.standard_error_m2_sr[lane],
                    weight,
                )?
                .0
                .abs();
                errors[pair.a].add(error)?;
                errors[pair.b].add(error)?;
                conductances[pair.a].add(g)?;
                conductances[pair.b].add(g)?;
            }
            let power_w = watts.value()?;
            if power_w != 0.0 {
                transfers.push(HeatTransfer {
                    a: self.cells[pair.a].id,
                    b: self.cells[pair.b].id,
                    power_w,
                });
            }
        }
        let absorbing_temperatures = self
            .cells
            .iter()
            .zip(&temperatures)
            .filter(|(cell, _)| cell.can_absorb)
            .map(|(_, t)| *t);
        let mut minimum_k = absorbing_temperatures.clone().fold(f64::INFINITY, f64::min);
        let mut maximum_k = absorbing_temperatures.fold(0.0_f64, f64::max);
        if let RadiativeBath::Blackbody { temperature_k } = self.bath {
            minimum_k = minimum_k.min(temperature_k);
            maximum_k = maximum_k.max(temperature_k);
        }
        let mut unresolved = CheckedSum::default();
        let mut unresolved_error = CheckedSum::default();
        for (index, cell) in self.cells.iter().enumerate() {
            check()?;
            let t = temperatures[index];
            let mut bath_w = CheckedSum::default();
            for lane in 0..N {
                let lambda = self.operator.wavelengths_nm[lane];
                let weight = self.weights_nm[lane];
                let reference = match self.bath {
                    RadiativeBath::Vacuum => None,
                    RadiativeBath::Blackbody { temperature_k } => Some(temperature_k),
                };
                let (outward_w, bath_g) = planck_radiative_exchange_w_and_w_k(
                    t,
                    reference,
                    lambda,
                    cell.escape.mean_m2_sr[lane],
                    weight,
                )?;
                let incoming = -outward_w;
                bath_w.add(incoming)?;
                bath_by_lane[lane].add(incoming)?;
                errors[index].add(
                    planck_radiative_exchange_w_and_w_k(
                        t,
                        reference,
                        lambda,
                        cell.escape.standard_error_m2_sr[lane],
                        weight,
                    )?
                    .0
                    .abs(),
                )?;
                conductances[index].add(bath_g)?;
                // A truncated directed path can become one half of a symmetric
                // pair (two receiving signs) or one full bath escape. Both are
                // bounded in heat-vector L1 by the maximum endpoint difference.
                // Scale BEFORE exponentiating, just as for actual heat: neither
                // an unrepresentable K*weight nor an underflowing raw Planck
                // difference may erase a representable final watt allowance.
                let cold = if matches!(self.bath, RadiativeBath::Vacuum) {
                    None
                } else {
                    Some(minimum_k)
                };
                unresolved.add(tail_scale(
                    t,
                    cold,
                    maximum_k,
                    lambda,
                    cell.unresolved.mean_m2_sr[lane],
                    weight,
                )?)?;
                unresolved_error.add(tail_scale(
                    t,
                    cold,
                    maximum_k,
                    lambda,
                    cell.unresolved.standard_error_m2_sr[lane],
                    weight,
                )?)?;
            }
            let power = bath_w.value()?;
            if power != 0.0 {
                bath_sources.insert(cell.id, power);
            }
        }
        let mut power = RadiativeExchangePower {
            medium_identity: self.operator.medium_identity,
            net_bath_to_material_w: [0.0; N],
            internal_transfer_magnitude_w: [0.0; N],
            thermal_external_w: 0.0,
            maximum_cell_standard_error_w: 0.0,
            total_standard_error_w: 0.0,
            unresolved_exchange_w: unresolved.value()?,
            unresolved_standard_error_w: unresolved_error.value()?,
        };
        let mut conductance_bounds_w_k = BTreeMap::new();
        let mut cell_errors_w = BTreeMap::new();
        let mut total_error = CheckedSum::default();
        for (index, cell) in self.cells.iter().enumerate() {
            check()?;
            let error = errors[index].value()?;
            total_error.add(error)?;
            power.maximum_cell_standard_error_w = power.maximum_cell_standard_error_w.max(error);
            cell_errors_w.insert(cell.id, error);
            let g = conductances[index].value()?;
            if g > 0.0 {
                conductance_bounds_w_k.insert(cell.id, g);
            }
        }
        power.total_standard_error_w = total_error.value()?;
        let mut spectral_external = CheckedSum::default();
        for lane in 0..N {
            power.net_bath_to_material_w[lane] = bath_by_lane[lane].value()?;
            power.internal_transfer_magnitude_w[lane] = internal_by_lane[lane].value()?;
            spectral_external.add(power.net_bath_to_material_w[lane])?;
        }
        let mut external = CheckedSum::default();
        for value in bath_sources.values() {
            external.add(*value)?;
        }
        power.thermal_external_w = external.value()?;
        if (power.thermal_external_w - spectral_external.value()?).abs()
            > self.tolerance.power_roundoff_w
        {
            return Err(invalid(
                "mapped bath power disagrees with spectral exchange",
            ));
        }
        if power.maximum_cell_standard_error_w > self.tolerance.max_cell_standard_error_w
            || power.total_standard_error_w > self.tolerance.max_total_standard_error_w
            || power.unresolved_exchange_w > self.tolerance.max_unresolved_exchange_w
            || power.unresolved_standard_error_w > self.tolerance.max_unresolved_standard_error_w
        {
            return Err(invalid(
                "conditional exchange sampling error or unresolved heat exceeds allowance",
            ));
        }
        check()?;
        Ok(Prepared {
            at_s: thermal.time_s(),
            temperatures_k: temperatures,
            sources_w: bath_sources,
            transfers,
            conductance_bounds_w_k,
            cell_errors_w,
            power,
        })
    }
}

impl<const N: usize> DepositionObserver for VolumeRadiativeExchange<N> {
    type Error = RadiativeExchangeError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        Ok(())
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        if self.evaluations >= self.evaluation_limit {
            return Err(RadiativeExchangeError::Budget);
        }
        let prepared = self.model.prepare(thermal, check)?;
        let sources = prepared.sources_w.clone();
        self.prepared = prepared;
        self.evaluations += 1;
        Ok(sources)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.check_prepared(thermal, check)?;
        Ok(self.prepared.transfers.clone())
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.check_prepared(thermal, check)?;
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.prepared.conductance_bounds_w_k,
            check,
        )?)
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        self.model.temperatures(thermal, &mut || Ok(()))?;
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || self.prepared.at_s != start_s
            || !positive(duration_s)
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "radiative energy requires the matching accepted native interval",
            ));
        }
        self.model.temperatures(thermal, check)?;
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        sums.external
            .add(finite(self.prepared.power.thermal_external_w * duration_s)?)?;
        sums.unresolved.add(finite(
            self.prepared.power.unresolved_exchange_w * duration_s,
        )?)?;
        for lane in 0..N {
            sums.bath[lane].add(finite(
                self.prepared.power.net_bath_to_material_w[lane] * duration_s,
            )?)?;
            sums.internal[lane].add(finite(
                self.prepared.power.internal_transfer_magnitude_w[lane] * duration_s,
            )?)?;
        }
        let energy = sums.value()?;
        check()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct EstimateSums<const N: usize> {
    mean: [CheckedSum; N],
    error: [CheckedSum; N],
}
impl<const N: usize> Default for EstimateSums<N> {
    fn default() -> Self {
        Self {
            mean: [CheckedSum::default(); N],
            error: [CheckedSum::default(); N],
        }
    }
}
impl<const N: usize> EstimateSums<N> {
    fn add(
        &mut self,
        value: &ExchangeCoefficientEstimate<N>,
    ) -> Result<(), RadiativeExchangeError> {
        for lane in 0..N {
            self.mean[lane].add(value.mean_m2_sr[lane])?;
            self.error[lane].add(value.standard_error_m2_sr[lane])?;
        }
        Ok(())
    }
    fn value(&self) -> Result<ExchangeCoefficientEstimate<N>, RadiativeExchangeError> {
        let mut value = ExchangeCoefficientEstimate {
            mean_m2_sr: [0.0; N],
            standard_error_m2_sr: [0.0; N],
        };
        for lane in 0..N {
            value.mean_m2_sr[lane] = self.mean[lane].value()?;
            value.standard_error_m2_sr[lane] = self.error[lane].value()?;
        }
        Ok(value)
    }
}
struct CellSums<const N: usize> {
    minimum_k: f64,
    maximum_k: f64,
    can_absorb: bool,
    escape: EstimateSums<N>,
    unresolved: EstimateSums<N>,
}
impl<const N: usize> CellSums<N> {
    fn new(minimum_k: f64, maximum_k: f64) -> Self {
        Self {
            minimum_k,
            maximum_k,
            can_absorb: false,
            escape: EstimateSums::default(),
            unresolved: EstimateSums::default(),
        }
    }
}
#[derive(Clone, Debug)]
struct EnergySums<const N: usize> {
    duration: CheckedSum,
    bath: [CheckedSum; N],
    internal: [CheckedSum; N],
    external: CheckedSum,
    unresolved: CheckedSum,
}
impl<const N: usize> Default for EnergySums<N> {
    fn default() -> Self {
        Self {
            duration: CheckedSum::default(),
            bath: [CheckedSum::default(); N],
            internal: [CheckedSum::default(); N],
            external: CheckedSum::default(),
            unresolved: CheckedSum::default(),
        }
    }
}
impl<const N: usize> EnergySums<N> {
    fn value(&self) -> Result<RadiativeExchangeEnergy<N>, RadiativeExchangeError> {
        let mut result = RadiativeExchangeEnergy {
            duration_s: self.duration.value()?,
            net_bath_to_material_j: [0.0; N],
            internal_transfer_magnitude_j: [0.0; N],
            thermal_external_j: self.external.value()?,
            unresolved_exchange_scale_j: self.unresolved.value()?,
        };
        for lane in 0..N {
            result.net_bath_to_material_j[lane] = self.bath[lane].value()?;
            result.internal_transfer_magnitude_j[lane] = self.internal[lane].value()?;
        }
        Ok(result)
    }
}
fn allocate<T: Clone>(count: usize, value: T) -> Result<Vec<T>, RadiativeExchangeError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| RadiativeExchangeError::Budget)?;
    result.resize(count, value);
    Ok(result)
}
fn transport_checkpoint<T>(
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(
        &mut dyn FnMut() -> Result<(), VolumeTransportError>,
    ) -> Result<T, VolumeTransportError>,
) -> Result<T, RadiativeExchangeError> {
    let mut refusal = None;
    let result = operation(&mut || match check() {
        Ok(()) => Ok(()),
        Err(e) => {
            refusal = Some(e);
            Err(VolumeTransportError::Cancelled)
        }
    });
    if let Some(e) = refusal {
        return Err(e.into());
    }
    result.map_err(Into::into)
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn nonnegative(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}
fn invalid(s: &'static str) -> RadiativeExchangeError {
    RadiativeExchangeError::Invalid(s)
}
fn finite(v: f64) -> Result<f64, RadiativeExchangeError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("unrepresentable radiative exchange"))
    }
}
fn tail_scale(
    t: f64,
    cold: Option<f64>,
    hot: f64,
    wavelength: f64,
    coefficient: f64,
    weight: f64,
) -> Result<f64, RadiativeExchangeError> {
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    let toward_cold =
        planck_radiative_exchange_w_and_w_k(t, cold, wavelength, coefficient, weight)?
            .0
            .abs();
    let toward_hot =
        planck_radiative_exchange_w_and_w_k(t, Some(hot), wavelength, coefficient, weight)?
            .0
            .abs();
    Ok(toward_cold.max(toward_hot))
}
