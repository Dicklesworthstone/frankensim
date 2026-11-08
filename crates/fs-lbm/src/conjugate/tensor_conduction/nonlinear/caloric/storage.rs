//! Original caloric admission and centered cell-energy coordinates.

use fs_exec::CancelGate;

use super::super::super::StorageRate;
use super::super::coordinate::CellCoordinate;
use super::{TensorCaloricStorage, source_error};
use crate::conjugate::tensor_conduction::invalid;
use crate::conjugate::{ChtError, ThermalSetup, VoxelDomain, finite, finite_positive, poll};

pub(super) struct PreparedCaloric<'a> {
    pub input: TensorCaloricStorage<'a>,
    coefficient_kg_s: Vec<f64>,
    bounds_j_kg: Vec<[f64; 2]>,
    previous_temperature_k: Vec<f64>,
}

pub(super) struct CaloricStage<'s, 'a> {
    pub storage: &'s PreparedCaloric<'a>,
    prescribed: Vec<Option<f64>>,
}

fn field<T>(cells: usize) -> Result<Vec<T>, ChtError> {
    cells
        .checked_mul(std::mem::size_of::<T>())
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| invalid("caloric tensor field exceeds addressable storage"))?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(cells)
        .map_err(|_| invalid("caloric tensor field allocation refused"))?;
    Ok(values)
}

impl<'a> PreparedCaloric<'a> {
    pub fn new(
        domain: &VoxelDomain,
        conductivity_entries: usize,
        input: TensorCaloricStorage<'a>,
        stage_interval_s: f64,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        finite_positive("tensor.caloric.stage_interval_s", stage_interval_s)?;
        let cells = domain.cell_count();
        if input.laws.len() != cells
            || input.reference_mass_kg.len() != cells
            || input.previous_specific_energy_j_kg.len() != cells
            || conductivity_entries != cells
        {
            return Err(invalid(
                "caloric tensor conduction needs one original law, mass, energy and tensor law per cell",
            ));
        }
        let mut coefficient_kg_s = field(cells)?;
        let mut bounds_j_kg = field(cells)?;
        let mut previous_temperature_k = field(cells)?;
        for cell in 0..cells {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            let mass = input.reference_mass_kg[cell];
            finite_positive("tensor.caloric.reference_mass_kg", mass)?;
            let volume = domain.volume(cell);
            finite_positive("tensor.caloric.reference_volume_m3", volume)?;
            input.laws[cell]
                .validate_reference_cell(mass, volume)
                .map_err(source_error)?;
            let bounds = input.laws[cell].specific_energy_bounds_j_kg();
            if bounds.iter().any(|e| !e.is_finite()) || bounds[0] >= bounds[1] {
                return Err(invalid(
                    "original caloric bounds must be finite and strictly ordered",
                ));
            }
            let previous = input.previous_specific_energy_j_kg[cell];
            if !previous.is_finite() || !(bounds[0]..=bounds[1]).contains(&previous) {
                return Err(invalid(
                    "previous specific energy is outside the original caloric bounds",
                ));
            }
            let [temperature, slope] = input.laws[cell]
                .temperature_and_slope_at_specific_energy(previous)
                .map_err(source_error)?;
            admit_state(temperature, slope)?;
            let coefficient = mass / stage_interval_s;
            finite_positive("tensor.caloric.mass_rate_kg_s", coefficient)?;
            coefficient_kg_s.push(coefficient);
            bounds_j_kg.push(bounds);
            previous_temperature_k.push(temperature);
        }
        poll(gate)?;
        Ok(Self {
            input,
            coefficient_kg_s,
            bounds_j_kg,
            previous_temperature_k,
        })
    }

    fn admit_energy(&self, cell: usize, energy: f64) -> Result<f64, ChtError> {
        let bounds = self.bounds_j_kg[cell];
        if energy.is_finite() && (bounds[0]..=bounds[1]).contains(&energy) {
            Ok(energy)
        } else {
            Err(invalid(
                "trial specific energy is outside the original caloric bounds",
            ))
        }
    }

    pub fn stage(
        &self,
        setup: &ThermalSetup,
        prescribed_specific_energy_j_kg: &[(usize, f64)],
        gate: &CancelGate,
    ) -> Result<CaloricStage<'_, 'a>, ChtError> {
        poll(gate)?;
        let cells = self.input.laws.len();
        if prescribed_specific_energy_j_kg.len() != setup.fixed_temperature.len() {
            return Err(invalid(
                "explicit energies must exactly cover all prescribed-temperature cells",
            ));
        }
        let mut prescribed = field(cells)?;
        prescribed.resize(cells, None);
        for (index, &(cell, energy)) in prescribed_specific_energy_j_kg.iter().enumerate() {
            if index.is_multiple_of(4096) {
                poll(gate)?;
            }
            if cell >= cells || prescribed[cell].is_some() {
                return Err(invalid(
                    "prescribed specific energy needs unique in-range cells",
                ));
            }
            prescribed[cell] = Some(self.admit_energy(cell, energy)?);
        }
        // Geometry preparation independently admits unique fixed-temperature
        // cells. Equal counts and coverage then establish exact correspondence.
        for (index, &(cell, temperature)) in setup.fixed_temperature.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            finite_positive("thermal.fixed_temperature", temperature)?;
            let energy = prescribed.get(cell).copied().flatten().ok_or_else(|| {
                invalid("every prescribed temperature needs an explicit specific energy")
            })?;
            let [actual, slope] = self.input.laws[cell]
                .temperature_and_slope_at_specific_energy(energy)
                .map_err(source_error)?;
            admit_state(actual, slope)?;
            if (actual - temperature).abs() > 64.0 * f64::EPSILON * actual.max(temperature) {
                return Err(invalid(
                    "prescribed specific energy does not give the authored cell temperature",
                ));
            }
            finite(
                "tensor.caloric.prescribed_energy_change_j_kg",
                energy - self.input.previous_specific_energy_j_kg[cell],
            )?;
        }
        poll(gate)?;
        Ok(CaloricStage {
            storage: self,
            prescribed,
        })
    }

    pub fn stored_j(&self, cell: usize, specific_energy_j_kg: f64) -> f64 {
        self.input.reference_mass_kg[cell]
            * (specific_energy_j_kg - self.input.previous_specific_energy_j_kg[cell])
    }

    pub fn storage_rate_w(&self, cell: usize, specific_energy_j_kg: f64) -> f64 {
        self.coefficient_kg_s[cell]
            * (specific_energy_j_kg - self.input.previous_specific_energy_j_kg[cell])
    }
}

fn admit_state(temperature: f64, slope: f64) -> Result<(), ChtError> {
    finite_positive("tensor.caloric.temperature_k", temperature)?;
    finite("tensor.caloric.temperature_slope_kg_k_j", slope)?;
    if slope < 0.0 {
        return Err(invalid(
            "original caloric temperature derivative must be nonnegative",
        ));
    }
    Ok(())
}

impl CaloricStage<'_, '_> {
    pub fn specific_energy(&self, cell: usize, value: f64) -> Result<f64, ChtError> {
        let energy = self.prescribed[cell]
            .unwrap_or_else(|| self.storage.input.previous_specific_energy_j_kg[cell] + value);
        self.storage.admit_energy(cell, energy)
    }
}

impl CellCoordinate for CaloricStage<'_, '_> {
    const IS_TEMPERATURE: bool = false;

    fn reference_temperatures(&self) -> Option<&[f64]> {
        Some(&self.storage.previous_temperature_k)
    }

    fn normalize(&self, cell: usize, value: f64) -> Result<f64, ChtError> {
        finite("tensor.caloric.energy_increment_j_kg", value)?;
        let energy = self.specific_energy(cell, value)?;
        let change = energy - self.storage.input.previous_specific_energy_j_kg[cell];
        finite("tensor.caloric.energy_change_j_kg", change)?;
        Ok(change)
    }

    fn temperature_and_slope(&self, cell: usize, value: f64) -> Result<[f64; 2], ChtError> {
        let energy = self.specific_energy(cell, value)?;
        let state = self.storage.input.laws[cell]
            .temperature_and_slope_at_specific_energy(energy)
            .map_err(source_error)?;
        admit_state(state[0], state[1])?;
        Ok(state)
    }

    fn prescribed_value(&self, cell: usize, _temperature_k: f64) -> Result<f64, ChtError> {
        let energy = self.prescribed[cell]
            .ok_or_else(|| invalid("prescribed cell lacks an explicit specific energy"))?;
        Ok(energy - self.storage.input.previous_specific_energy_j_kg[cell])
    }

    fn storage_rate(&self, cell: usize, value: f64) -> Result<StorageRate, ChtError> {
        // Audit the same representable energy returned to the caller. Do not
        // accept a heat increment hidden by rounding the absolute energy datum.
        let energy = self.specific_energy(cell, value)?;
        let power_w = self.storage.storage_rate_w(cell, energy);
        finite("tensor.caloric.storage_power_w", power_w)?;
        Ok(StorageRate {
            power_w,
            arithmetic_w: power_w.abs(),
        })
    }

    fn storage_tangent(&self, cell: usize) -> f64 {
        self.storage.coefficient_kg_s[cell]
    }
}
