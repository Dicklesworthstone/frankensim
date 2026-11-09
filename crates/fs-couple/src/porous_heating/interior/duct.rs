//! Original molecular duct viscosity on the native fixed-volume thermal clock.
//!
//! The native cylinder/cone forms supply positive, already integrated and
//! cycle-averaged viscous watts from one loaded source bank. This adapter adds
//! them once to native internal-energy cells. Fourier/molecular availability
//! and mechanical-wall work are observed separately, never deposited.
//!
//! Acoustic material and coherent source statistics are held during an advance.
//! The original caloric object, its Cv/datum and its correspondence to the
//! source gas are caller-supplied physics: equal rho/T and an identity do not
//! establish shared provenance. No storage law, pressure work, streaming,
//! acoustic startup, material refresh or within-Magnus-panel field is invented.

use super::{CellBinding, InteriorHeatCell, InteriorHeatingError, finite, same_measure};
use fs_duct::DuctError;
use fs_duct::network::field::MAX_FIELD_PROBES;
use fs_duct::network::molecular::{
    MolecularNetworkCylinderAbsorption, PreparedMolecularNetworkConeAbsorption,
    PreparedMolecularNetworkCylinderAbsorption, PreparedMolecularNetworkFields,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig, ThermalStorageBasis,
};
use std::collections::BTreeMap;

/// Geometry identity in the original authored graph, without remeshing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DuctHeatGeometry {
    Cylinder,
    Cone { slice: usize },
}

/// One original fluid interval and its explicitly authored caloric destination.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctCellHeating {
    pub cell: CellId,
    pub edge_id: String,
    pub segment: usize,
    pub geometry: DuctHeatGeometry,
    /// Segment-local interval; different cone slices share this coordinate.
    pub interval_m: [f64; 2],
    pub volume_m3: f64,
    pub mass_kg: f64,
    /// Native integrated cycle-mean powers. Only viscous_power_w is deposited.
    pub power: MolecularNetworkCylinderAbsorption,
}

/// Original acoustic, thermal or retained physical-binding refusal.
#[derive(Debug)]
pub enum DuctHeatingError {
    Acoustic(DuctError),
    Thermal(DepositionError),
    Binding(InteriorHeatingError),
    TemperatureExcursion {
        cell: CellId,
        actual_k: f64,
        maximum_k: f64,
    },
}
impl core::fmt::Display for DuctHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular duct viscous heating refused: {self:?}")
    }
}
impl std::error::Error for DuctHeatingError {}
impl From<DuctError> for DuctHeatingError {
    fn from(value: DuctError) -> Self {
        Self::Acoustic(value)
    }
}
impl From<DepositionError> for DuctHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<InteriorHeatingError> for DuctHeatingError {
    fn from(value: InteriorHeatingError) -> Self {
        Self::Binding(value)
    }
}
fn invalid(message: &'static str) -> DuctHeatingError {
    InteriorHeatingError::Invalid(message).into()
}

struct BoundCell<'t> {
    heat: DuctCellHeating,
    binding: CellBinding<'t>,
    reference_temperature_k: f64,
}

/// Source-retaining held heat from one original loaded molecular field bank.
///
/// Both optional native forms and the explicit mean/covariance factors stay
/// borrowed. Cold preparation evaluates coherence once per bank; native
/// thermal substeps only use cached watts and original binding checks.
pub struct MolecularDuctViscousHeating<'p, 'f, 'n, 't> {
    cylinders: Option<&'p PreparedMolecularNetworkCylinderAbsorption<'f, 'n>>,
    cones: Option<&'p PreparedMolecularNetworkConeAbsorption<'f, 'n>>,
    fields: &'f PreparedMolecularNetworkFields<'n>,
    mean: Option<&'p [C64]>,
    factors: &'p [&'p [C64]],
    cells: Vec<BoundCell<'t>>,
    sources_w: BTreeMap<CellId, f64>,
    maximum_temperature_change_k: f64,
}

/// Separate mechanisms integrated over accepted native thermal durations.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DuctHeatingEnergy {
    pub viscous_heat_j: f64,
    /// Fourier availability destruction; not a local caloric source.
    pub thermal_availability_loss_j: f64,
    /// Molecular availability destruction; not a local caloric source.
    pub molecular_availability_loss_j: f64,
    /// Mechanical-wall work, possibly including radiation; not deposited.
    pub wall_work_j: f64,
    /// Viscous heat only, already included once in native external-heat work.
    pub deposited_j: f64,
}

pub struct DuctHeatingStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// All original thermal sources/exchanges plus viscous heat exactly once.
    pub heat_report: EnergyAudit,
    pub acoustic: DuctHeatingEnergy,
    /// Viscous joules, in cylinder-then-cone input order, never CellId order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

impl<'p, 'f, 'n, 't> MolecularDuctViscousHeating<'p, 'f, 'n, 't> {
    /// Bind every selected cylinder/cone interval to a distinct native cell.
    /// Destinations follow cylinder cells, then cone cells, each in its original
    /// order. Selected intervals may cover a subset of the graph, but must have
    /// positive volume and must not overlap on an authored edge/segment.
    ///
    /// Both forms must borrow the SAME prepared field object. Mean and factors
    /// use its canonical source coordinates and retain interference within each
    /// vector. None plus no factors means zero excitation, not default sources.
    /// No cross-frequency or finite-shutter source interpretation is inferred.
    ///
    /// Original interval length, volume, density*volume and initial temperature
    /// must match each destination to roundoff. Its exact caloric object must
    /// identify InternalEnergy storage; enthalpy/Cp is not fixed-volume gas
    /// storage. World placement/orientation and original Cv/source provenance
    /// remain caller-owned. Exposed convection area is left unchanged.
    ///
    /// Every native event and substep checks the declared temperature excursion
    /// from the frozen gas state. It is an approximation admission, not a
    /// certificate for unchanged constitutive data, phase or composition.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cylinders: Option<&'p PreparedMolecularNetworkCylinderAbsorption<'f, 'n>>,
        cones: Option<&'p PreparedMolecularNetworkConeAbsorption<'f, 'n>>,
        mean: Option<&'p [C64]>,
        factors: &'p [&'p [C64]],
        thermal: &DepositionSimulation<'_>,
        destinations: &[InteriorHeatCell<'t>],
        maximum_cells: usize,
        maximum_temperature_change_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DuctHeatingError> {
        checkpoint()?;
        let fields = cylinders
            .map(|value| value.fields())
            .or_else(|| cones.map(|value| value.fields()))
            .ok_or_else(|| invalid("duct heating needs an original native loss bank"))?;
        if let (Some(a), Some(b)) = (cylinders, cones) {
            if !core::ptr::eq(a.fields(), b.fields()) {
                return Err(invalid(
                    "duct heat banks must share the exact original loaded fields",
                ));
            }
        }
        let count = cylinders.map_or(0, |v| v.cells().len()) + cones.map_or(0, |v| v.cells().len());
        if maximum_cells == 0
            || maximum_cells > MAX_FIELD_PROBES
            || count == 0
            || count > maximum_cells
        {
            return Err(DepositionError::Budget("molecular duct heat cells").into());
        }
        if destinations.len() != count
            || !maximum_temperature_change_k.is_finite()
            || maximum_temperature_change_k <= 0.0
        {
            return Err(invalid(
                "duct heat needs one destination per cell and a finite positive temperature excursion",
            ));
        }
        let cylindrical = cylinders
            .into_iter()
            .flat_map(|bank| bank.cells())
            .map(|cell| {
                (
                    cell.probe_indices[0],
                    cell.segment,
                    [cell.start_m, cell.end_m],
                    cell.volume_m3,
                    DuctHeatGeometry::Cylinder,
                )
            });
        let conical = cones.into_iter().flat_map(|bank| bank.cells()).map(|cell| {
            (
                cell.probe_indices[0],
                cell.segment,
                [cell.start_m, cell.end_m],
                cell.volume_m3,
                DuctHeatGeometry::Cone { slice: cell.slice },
            )
        });
        let mut cells = Vec::with_capacity(count);
        let mut sources_w = BTreeMap::new();
        let mut intervals: BTreeMap<(String, usize), Vec<[f64; 2]>> = BTreeMap::new();
        for ((probe, segment, interval, volume, geometry), route) in
            cylindrical.chain(conical).zip(destinations)
        {
            checkpoint()?;
            if sources_w.insert(route.cell, 0.0).is_some() {
                return Err(invalid("duct heat destinations must be distinct"));
            }
            if interval[1] <= interval[0] || !volume.is_finite() || volume <= 0.0 {
                return Err(invalid(
                    "duct thermal intervals must have positive length and volume",
                ));
            }
            let location = &fields.fields().locations()[probe];
            let edge_id = location.probe.edge_id.clone();
            let gas = fields
                .network()
                .edge(&edge_id)
                .expect("loss bank retains its original edge")
                .state()
                .gas();
            let binding = CellBinding::try_new(thermal, route.cell, route.caloric)?;
            let delta = std::array::from_fn::<_, 3, _>(|i| {
                binding.geometry.end_m[i] - binding.geometry.start_m[i]
            });
            let length = det::hypot(det::hypot(delta[0], delta[1]), delta[2]);
            let mass = finite(gas.density * volume)?;
            if !same_measure(length, interval[1] - interval[0])
                || !same_measure(binding.geometry.volume_m3, volume)
                || !same_measure(binding.mass_kg, mass)
                || !same_measure(thermal.temperature_k(route.cell)?, gas.temperature)
                || binding.storage_basis != ThermalStorageBasis::InternalEnergy
                || binding.caloric_identity.is_none()
            {
                return Err(invalid(
                    "duct heat needs original length, volume, fluid mass, temperature and identified internal-energy storage",
                ));
            }
            intervals
                .entry((edge_id.clone(), segment))
                .or_default()
                .push(interval);
            cells.push(BoundCell {
                heat: DuctCellHeating {
                    cell: route.cell,
                    edge_id,
                    segment,
                    geometry,
                    interval_m: interval,
                    volume_m3: volume,
                    mass_kg: mass,
                    power: MolecularNetworkCylinderAbsorption::default(),
                },
                binding,
                reference_temperature_k: gas.temperature,
            });
        }
        for entries in intervals.values_mut() {
            checkpoint()?;
            entries.sort_by(|a, b| a[0].total_cmp(&b[0]));
            if entries.windows(2).any(|pair| pair[0][1] > pair[1][0]) {
                return Err(invalid(
                    "duct heat intervals overlap on the same original edge and segment",
                ));
            }
        }
        let mut refusal = None;
        let computed = (|| {
            let mut cancelled = || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            };
            let mut powers = Vec::with_capacity(count);
            if let Some(bank) = cylinders {
                powers.extend(bank.observe_coherence(mean, factors, &mut cancelled)?);
            }
            if let Some(bank) = cones {
                powers.extend(bank.observe_coherence(mean, factors, &mut cancelled)?);
            }
            Ok::<_, DuctError>(powers)
        })();
        if let Some(error) = refusal {
            return Err(error.into());
        }
        for (cell, power) in cells.iter_mut().zip(computed?) {
            checkpoint()?;
            if channels(power)
                .iter()
                .chain([power.absorbed_power_w].iter())
                .any(|v| !v.is_finite() || *v < 0.0)
            {
                return Err(invalid(
                    "duct heat needs finite passive original loss forms",
                ));
            }
            *sources_w
                .get_mut(&cell.heat.cell)
                .expect("destination admitted") = power.viscous_power_w;
            cell.heat.power = power;
        }
        checkpoint()?;
        Ok(Self {
            cylinders,
            cones,
            fields,
            mean,
            factors,
            cells,
            sources_w,
            maximum_temperature_change_k,
        })
    }

    pub fn fields(&self) -> &'f PreparedMolecularNetworkFields<'n> {
        self.fields
    }
    pub fn cylinders(&self) -> Option<&'p PreparedMolecularNetworkCylinderAbsorption<'f, 'n>> {
        self.cylinders
    }
    pub fn cones(&self) -> Option<&'p PreparedMolecularNetworkConeAbsorption<'f, 'n>> {
        self.cones
    }
    pub fn mean(&self) -> Option<&'p [C64]> {
        self.mean
    }
    pub fn factors(&self) -> &'p [&'p [C64]] {
        self.factors
    }
    pub fn cells(&self) -> impl ExactSizeIterator<Item = &DuctCellHeating> {
        self.cells.iter().map(|cell| &cell.heat)
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DuctHeatingError> {
        for cell in &self.cells {
            checkpoint()?;
            cell.binding.validate(thermal)?;
            let change =
                (thermal.temperature_k(cell.heat.cell)? - cell.reference_temperature_k).abs();
            if change > self.maximum_temperature_change_k {
                return Err(DuctHeatingError::TemperatureExcursion {
                    cell: cell.heat.cell,
                    actual_k: change,
                    maximum_k: self.maximum_temperature_change_k,
                });
            }
        }
        Ok(())
    }

    /// Start a fresh account. Continue with the accepted observer returned by
    /// native advance_observed when accumulating across more than one call.
    pub fn observer(&self) -> DuctHeatingObserver<'_, 'p, 'f, 'n, 't> {
        DuctHeatingObserver {
            owner: self,
            cell_energy: vec![Dd::ZERO; self.cells.len()],
            mechanisms: [Dd::ZERO; 4],
        }
    }

    /// Held cycle-mean viscosity joins the original authored thermal sources,
    /// links, reservoirs and events. Do not author these same watts separately.
    /// Native candidate state and account publish together; refusal leaves
    /// the supplied thermal state and source bank unchanged.
    pub fn advance_to<'s>(
        &self,
        initial: &DepositionSimulation<'s>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<DuctHeatingStep<'s>, DuctHeatingError> {
        self.validate(initial, &mut checkpoint)?;
        let (thermal, observer, heat_report) =
            initial.advance_observed(&self.observer(), target_s, config, &mut checkpoint)?;
        self.validate(&thermal, &mut checkpoint)?;
        let acoustic = observer.energy()?;
        let cell_energy_j = observer.cell_energy_j()?;
        checkpoint()?;
        Ok(DuctHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

fn channels(power: MolecularNetworkCylinderAbsorption) -> [f64; 4] {
    [
        power.viscous_power_w,
        power.thermal_power_w,
        power.molecular_power_w,
        power.wall_work_w,
    ]
}

/// Cloneable native observer: immutable source/bindings stay borrowed while
/// candidate accepted-work accounts are copied, never shared mutably.
#[derive(Clone)]
pub struct DuctHeatingObserver<'h, 'p, 'f, 'n, 't> {
    owner: &'h MolecularDuctViscousHeating<'p, 'f, 'n, 't>,
    cell_energy: Vec<Dd>,
    mechanisms: [Dd; 4],
}
impl DuctHeatingObserver<'_, '_, '_, '_, '_> {
    pub fn energy(&self) -> Result<DuctHeatingEnergy, DuctHeatingError> {
        let deposited = self
            .cell_energy
            .iter()
            .copied()
            .fold(Dd::ZERO, |a, b| a + b);
        Ok(DuctHeatingEnergy {
            viscous_heat_j: finite(self.mechanisms[0].to_f64())?,
            thermal_availability_loss_j: finite(self.mechanisms[1].to_f64())?,
            molecular_availability_loss_j: finite(self.mechanisms[2].to_f64())?,
            wall_work_j: finite(self.mechanisms[3].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
        })
    }
    pub fn cell_energy_j(&self) -> Result<Vec<(CellId, f64)>, DuctHeatingError> {
        self.owner
            .cells
            .iter()
            .zip(&self.cell_energy)
            .map(|(cell, energy)| Ok((cell.heat.cell, finite(energy.to_f64())?)))
            .collect()
    }
}
impl DepositionObserver for DuctHeatingObserver<'_, '_, '_, '_, '_> {
    type Error = DuctHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        Ok(self.owner.sources_w.clone())
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.owner.validate(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        for (cell, energy) in self.owner.cells.iter().zip(&mut self.cell_energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(cell.heat.power.viscous_power_w * dt)?);
            for (sum, watts) in self.mechanisms.iter_mut().zip(channels(cell.heat.power)) {
                *sum = *sum + Dd::from_f64(finite(watts * dt)?);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
