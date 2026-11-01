//! Optical constitutive data on the original native specific-enthalpy axis.
//!
//! A sharp latent interval has many states at one temperature. These tables
//! retain those states without a second caloric solver or an inferred phase
//! mixing rule. Each optical source names the exact immutable caloric identity
//! and therefore its energy datum; cells additionally name their native material.
//! Geometry, density and angular phase functions remain explicitly fixed.

use super::*;
use crate::volumes::transport::MassAbsorptionSpectrum;
use fs_thermal::enthalpy::ThermalStorageBasis;

#[derive(Clone, Debug)]
struct EnthalpyTable<S> {
    source_identity: ContentHash,
    caloric_identity: ContentHash,
    temperature_bounds_k: (f64, f64),
    wavelength_bounds_nm: (f64, f64),
    rows: Vec<(f64, S)>,
}

impl<S: Clone + TemperatureSpectrum> EnthalpyTable<S> {
    fn try_new(
        source_identity: ContentHash,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, S)],
        max_enthalpy_knots: usize,
        spectral_bounds: fn(&S) -> (f64, f64),
    ) -> Result<Self, OpticalFeedbackError> {
        if source_identity == ContentHash([0; 32])
            || caloric_identity == ContentHash([0; 32])
            || rows.len() < 2
            || rows.len() > max_enthalpy_knots
            || !temperature_bounds_k.0.is_finite()
            || !temperature_bounds_k.1.is_finite()
            || temperature_bounds_k.0 <= 0.0
            || temperature_bounds_k.0 > temperature_bounds_k.1
            || rows.iter().any(|(h, _)| !h.is_finite())
            || rows.windows(2).any(|pair| {
                let span = pair[1].0 - pair[0].0;
                !span.is_finite() || span <= 0.0
            })
        {
            return Err(invalid(
                "bounded ordered finite enthalpy rows, positive kelvin domain and optical/caloric identities required",
            ));
        }
        let wavelength_bounds_nm = spectral_bounds(&rows[0].1);
        if rows
            .iter()
            .any(|(_, spectrum)| spectral_bounds(spectrum) != wavelength_bounds_nm)
        {
            return Err(invalid("enthalpy rows must share one wavelength domain"));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        retained.extend_from_slice(rows);
        Ok(Self {
            source_identity,
            caloric_identity,
            temperature_bounds_k,
            wavelength_bounds_nm,
            rows: retained,
        })
    }

    fn compile<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
        hash: &mut DomainHasher,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CompiledOpticalTable<N>, OpticalFeedbackError> {
        hash.update(self.caloric_identity.as_bytes());
        hash.update(&self.temperature_bounds_k.0.to_bits().to_le_bytes());
        hash.update(&self.temperature_bounds_k.1.to_bits().to_le_bytes());
        Ok(CompiledOpticalTable {
            table: CompiledTable::compile(
                self.source_identity,
                self.wavelength_bounds_nm,
                &self.rows,
                wavelengths_nm,
                hash,
                checkpoint,
            )?,
            temperature_bounds_k: self.temperature_bounds_k,
            enthalpy: Some(EnthalpyTableDomain {
                caloric_identity: self.caloric_identity,
                minimum_j_kg: self.rows[0].0,
                maximum_j_kg: self.rows[self.rows.len() - 1].0,
            }),
        })
    }
}

/// Source-bound mass absorption [m²/kg] versus the original caloric chart's
/// specific enthalpy [J/kg]. Interpolation is linear in h and in the original
/// spectrum's wavelength coordinate. Neither coordinate extrapolates.
#[derive(Clone, Debug)]
pub struct EnthalpyAbsorptionTable(EnthalpyTable<MassAbsorptionSpectrum>);
impl EnthalpyAbsorptionTable {
    pub fn try_new(
        source_identity: ContentHash,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, MassAbsorptionSpectrum)],
        max_enthalpy_knots: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        EnthalpyTable::try_new(
            source_identity,
            caloric_identity,
            temperature_bounds_k,
            rows,
            max_enthalpy_knots,
            MassAbsorptionSpectrum::wavelength_bounds_nm,
        )
        .map(Self)
    }
    pub fn source_identity(&self) -> ContentHash {
        self.0.source_identity
    }
    pub fn caloric_identity(&self) -> ContentHash {
        self.0.caloric_identity
    }
    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        self.0.temperature_bounds_k
    }
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.0.wavelength_bounds_nm
    }
    pub fn specific_enthalpy_bounds_j_kg(&self) -> (f64, f64) {
        (self.0.rows[0].0, self.0.rows[self.0.rows.len() - 1].0)
    }
}

/// Source-bound mass scattering [m²/kg] versus native specific enthalpy [J/kg].
/// The angular scattering phase law remains the field's explicit fixed input.
#[derive(Clone, Debug)]
pub struct EnthalpyScatteringTable(EnthalpyTable<MassScatteringSpectrum>);
impl EnthalpyScatteringTable {
    pub fn try_new(
        source_identity: ContentHash,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, MassScatteringSpectrum)],
        max_enthalpy_knots: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        EnthalpyTable::try_new(
            source_identity,
            caloric_identity,
            temperature_bounds_k,
            rows,
            max_enthalpy_knots,
            MassScatteringSpectrum::wavelength_bounds_nm,
        )
        .map(Self)
    }
    pub fn source_identity(&self) -> ContentHash {
        self.0.source_identity
    }
    pub fn caloric_identity(&self) -> ContentHash {
        self.0.caloric_identity
    }
    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        self.0.temperature_bounds_k
    }
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.0.wavelength_bounds_nm
    }
    pub fn specific_enthalpy_bounds_j_kg(&self) -> (f64, f64) {
        (self.0.rows[0].0, self.0.rows[self.0.rows.len() - 1].0)
    }
}

/// One optical cell bound to the exact original native material and enthalpy
/// receiver. Optical tables must name that registered caloric law's identity.
#[derive(Clone, Copy, Debug)]
pub struct EnthalpyScatteringCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub material_id: MaterialId,
    pub absorption_table: usize,
    pub scattering_table: usize,
    pub phase: usize,
}

/// Immutable field storage bounds; photon and trajectory work remain owned by
/// the existing incident/radiative consumers.
#[derive(Clone, Copy, Debug)]
pub struct EnthalpyOpticalFieldBudget {
    pub max_grid_cells: usize,
    pub max_tables: usize,
    pub max_total_enthalpy_knots: usize,
}

impl<'g, const N: usize> TemperatureOpticalField<'g, N> {
    /// Compile one field whose absorption and scattering read native h. Every
    /// cell must use enthalpy storage and the declared original caloric owner.
    /// The supplied kelvin domain is an additional source restriction, not an
    /// independently integrated temperature or inferred h-to-phase conversion.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_enthalpy(
        geometry: AbsorbingGrid<'g>,
        wavelengths_nm: [f64; N],
        absorption: &[EnthalpyAbsorptionTable],
        scattering: &[EnthalpyScatteringTable],
        phases: &[PhaseFunction],
        bindings: &[EnthalpyScatteringCell],
        thermal: &DepositionSimulation<'_>,
        budget: EnthalpyOpticalFieldBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        let count = geometry.cell_count();
        let tables = absorption
            .len()
            .checked_add(scattering.len())
            .ok_or(OpticalFeedbackError::Budget)?;
        if N == 0
            || count > budget.max_grid_cells
            || tables > budget.max_tables
            || bindings.len() > count
            || phases.len() > count
        {
            return Err(OpticalFeedbackError::Budget);
        }
        if absorption.is_empty()
            || scattering.is_empty()
            || phases.is_empty()
            || thermal.cells().is_empty()
        {
            return Err(invalid(
                "nonempty optical laws and active thermal cells required",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.enthalpy-optical-field.v1");
        hash.update(geometry.content_identity().as_bytes());
        hash.update(&(N as u64).to_le_bytes());
        for wavelength in wavelengths_nm {
            hash.update(&wavelength.to_bits().to_le_bytes());
        }
        let mut bounds = (0.0_f64, f64::MAX);
        let mut knots = 0usize;
        let absorption = compile_tables(
            absorption.iter().map(|table| &table.0),
            &wavelengths_nm,
            &mut bounds,
            &mut knots,
            budget.max_total_enthalpy_knots,
            &mut hash,
            &mut checkpoint,
        )?;
        let scattering = compile_tables(
            scattering.iter().map(|table| &table.0),
            &wavelengths_nm,
            &mut bounds,
            &mut knots,
            budget.max_total_enthalpy_knots,
            &mut hash,
            &mut checkpoint,
        )?;
        Self::compile_prepared(
            geometry,
            wavelengths_nm,
            PreparedTables {
                absorption,
                scattering,
                bounds_nm: bounds,
                knots,
            },
            phases,
            bindings.iter().map(|binding| PreparedBinding {
                cell: ThermalScatteringCell {
                    grid_cell: binding.grid_cell,
                    thermal_cell: binding.thermal_cell,
                    absorption_table: binding.absorption_table,
                    scattering_table: binding.scattering_table,
                    phase: binding.phase,
                },
                material_id: Some(binding.material_id),
            }),
            thermal,
            hash,
            &mut checkpoint,
        )
        .map(|(field, _)| field)
    }
}

fn compile_tables<'a, S: Clone + TemperatureSpectrum + 'a, const N: usize>(
    tables: impl ExactSizeIterator<Item = &'a EnthalpyTable<S>>,
    wavelengths_nm: &[f64; N],
    bounds: &mut (f64, f64),
    knots: &mut usize,
    max_knots: usize,
    hash: &mut DomainHasher,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<CompiledOpticalTable<N>>, OpticalFeedbackError> {
    let mut compiled = Vec::new();
    compiled
        .try_reserve_exact(tables.len())
        .map_err(|_| OpticalFeedbackError::Budget)?;
    hash.update(&(tables.len() as u64).to_le_bytes());
    for table in tables {
        checkpoint()?;
        *knots = knots
            .checked_add(table.rows.len())
            .ok_or(OpticalFeedbackError::Budget)?;
        if *knots > max_knots {
            return Err(OpticalFeedbackError::Budget);
        }
        bounds.0 = bounds.0.max(table.wavelength_bounds_nm.0);
        bounds.1 = bounds.1.min(table.wavelength_bounds_nm.1);
        compiled.push(table.compile(wavelengths_nm, hash, checkpoint)?);
    }
    Ok(compiled)
}

#[derive(Clone, Copy)]
pub(super) struct EnthalpyTableDomain {
    caloric_identity: ContentHash,
    minimum_j_kg: f64,
    maximum_j_kg: f64,
}
#[derive(Clone, Copy)]
pub(super) struct BoundEnthalpy {
    material_id: MaterialId,
    domain: EnthalpyTableDomain,
}
impl BoundEnthalpy {
    pub(super) fn from_tables(
        material_id: Option<MaterialId>,
        absorption: Option<EnthalpyTableDomain>,
        scattering: Option<EnthalpyTableDomain>,
    ) -> Result<Option<Self>, OpticalFeedbackError> {
        match (material_id, absorption, scattering) {
            (None, None, None) => Ok(None),
            (Some(material_id), Some(a), Some(s)) => Self::intersection(
                Some(Self {
                    material_id,
                    domain: a,
                }),
                Some(Self {
                    material_id,
                    domain: s,
                }),
            ),
            _ => Err(invalid(
                "optical coordinate and native material bindings disagree",
            )),
        }
    }
    pub(super) fn intersection(
        a: Option<Self>,
        b: Option<Self>,
    ) -> Result<Option<Self>, OpticalFeedbackError> {
        match (a, b) {
            (None, None) => Ok(None),
            (Some(mut a), Some(b))
                if a.material_id == b.material_id
                    && a.domain.caloric_identity == b.domain.caloric_identity =>
            {
                a.domain.minimum_j_kg = a.domain.minimum_j_kg.max(b.domain.minimum_j_kg);
                a.domain.maximum_j_kg = a.domain.maximum_j_kg.min(b.domain.maximum_j_kg);
                if a.domain.minimum_j_kg > a.domain.maximum_j_kg {
                    return Err(invalid("optical enthalpy domains do not overlap"));
                }
                Ok(Some(a))
            }
            _ => Err(invalid(
                "optical enthalpy tables must bind the same original caloric material and datum",
            )),
        }
    }
    pub(super) fn native_value(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<f64, OpticalFeedbackError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(OpticalHeatingError::InactiveCell(id))?;
        if cell.material_id != self.material_id
            || thermal.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy
            || thermal.material_caloric_identity(cell.material_id)?
                != Some(*self.domain.caloric_identity.as_bytes())
        {
            return Err(invalid(
                "native enthalpy optical law requires its exact material, caloric identity and storage basis",
            ));
        }
        let h = cell.specific_h_j_kg;
        if !h.is_finite() || !(self.domain.minimum_j_kg..=self.domain.maximum_j_kg).contains(&h) {
            return Err(OpticalFeedbackError::EnthalpyDomain {
                thermal_cell: id,
                specific_h_j_kg: h,
                minimum_j_kg: self.domain.minimum_j_kg,
                maximum_j_kg: self.domain.maximum_j_kg,
            });
        }
        Ok(h)
    }
}

#[cfg(test)]
mod tests;
