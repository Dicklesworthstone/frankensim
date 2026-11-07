//! Temperature-dependent two-photon absorption on the original thermal clock.
//!
//! The monochromatic path shares the original geometry cache, linear optical
//! tables, positive analytic attenuation kernel, native enthalpy observer and
//! complete-trajectory refinement. Only beta changes here: it is an explicit
//! bulk coefficient in m/W interpolated in kelvin, without density scaling.
//! Each source declares its temperature/wavelength support and local incident
//! irradiance ceiling. A changing upstream opacity can change downstream
//! admission, so both applied rates and final/event endpoint states are checked.
//!
//! Geometry, density, beam area and incident power remain prescribed. The
//! channel promptly thermalizes absorbed energy; carrier storage, saturation,
//! nonlinear refraction and cross-frequency absorption are not inferred. The
//! explicit native step freezes coefficients over one step. Refinement measures
//! that time error; the exact spatial cell solution does not remove it.

use super::*;
use crate::volumes::transport::deposition::nonlinear::absorb_two_photon;

#[derive(Debug)]
struct Source {
    identity: ContentHash,
    source_identity: ContentHash,
    rows: Vec<(f64, f64)>,
    wavelengths_nm: (f64, f64),
    maximum_irradiance_w_m2: f64,
}

/// Source-declared bulk beta(T), with linear interpolation and no extrapolation.
/// Clones share the immutable original samples. Beta is constant in wavelength
/// over the explicitly declared band and already includes both photons' energy.
#[derive(Clone, Debug)]
pub struct TemperatureTwoPhotonTable(Arc<Source>);

impl TemperatureTwoPhotonTable {
    /// Admit `(temperature_k, beta_m_w)` knots and their common source domain.
    /// The source must cover the original density, polarization and material
    /// state throughout this temperature band; no such coverage is inferred.
    pub fn try_new(
        source_identity: ContentHash,
        rows: &[(f64, f64)],
        wavelengths_nm: (f64, f64),
        maximum_irradiance_w_m2: f64,
        max_temperature_knots: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        if rows.len() > max_temperature_knots {
            return Err(OpticalFeedbackError::Budget);
        }
        let (lo, hi) = wavelengths_nm;
        if source_identity == ContentHash([0; 32])
            || rows.len() < 2
            || rows
                .iter()
                .any(|&(t, beta)| !t.is_finite() || t <= 0.0 || !nonnegative(beta))
            || rows.windows(2).any(|pair| pair[0].0 >= pair[1].0)
            || !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || hi <= lo
            || !maximum_irradiance_w_m2.is_finite()
            || maximum_irradiance_w_m2 <= 0.0
        {
            return Err(invalid(
                "ordered positive temperatures, nonnegative bulk beta and explicit nonlinear source support required",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        retained.extend_from_slice(rows);
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-two-photon-table.v1");
        hash.update(source_identity.as_bytes());
        hash.update(&(rows.len() as u64).to_le_bytes());
        for value in [lo, hi, maximum_irradiance_w_m2] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for &(temperature, beta) in rows {
            hash.update(&temperature.to_bits().to_le_bytes());
            hash.update(&beta.to_bits().to_le_bytes());
        }
        Ok(Self(Arc::new(Source {
            identity: hash.finalize(),
            source_identity,
            rows: retained,
            wavelengths_nm,
            maximum_irradiance_w_m2,
        })))
    }

    /// Original source identity, distinct from the numerical table identity.
    pub fn source_identity(&self) -> ContentHash {
        self.0.source_identity
    }

    /// Inclusive supported temperature interval, kelvin.
    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        (self.0.rows[0].0, self.0.rows[self.0.rows.len() - 1].0)
    }

    fn at(&self, cell: CellId, temperature_k: f64) -> Result<f64, OpticalFeedbackError> {
        let (minimum_k, maximum_k) = self.temperature_bounds_k();
        if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
            return Err(OpticalFeedbackError::TemperatureDomain {
                thermal_cell: cell,
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        let rows = &self.0.rows;
        let upper = rows
            .partition_point(|row| row.0 < temperature_k)
            .clamp(1, rows.len() - 1);
        let (ta, a) = rows[upper - 1];
        let (tb, b) = rows[upper];
        let fraction = (temperature_k - ta) / (tb - ta);
        // The convex form avoids an overflowing subtraction of coefficients.
        let beta = if a == b {
            a
        } else {
            (1.0 - fraction) * a + fraction * b
        };
        if !nonnegative(beta) {
            return Err(invalid("unrepresentable temperature-dependent beta"));
        }
        Ok(beta)
    }
}

/// Independent nonlinear material assignment on an original non-vacuum path cell.
/// The receiving thermal cell is inherited from `ThermalOpticalCell`.
#[derive(Clone, Copy, Debug)]
pub struct ThermalTwoPhotonCell {
    /// Original DDA cell index.
    pub grid_cell: usize,
    /// Index in the separately supplied nonlinear table array.
    pub two_photon_table: usize,
}

/// Explicit monochromatic nonlinear data, prepared together with linear optics.
#[derive(Clone, Copy, Debug)]
pub struct TwoPhotonFeedback<'a> {
    /// Constant physical cross-sectional area of the uniform beam, m².
    pub beam_area_m2: f64,
    /// Independently sourced beta(T) laws.
    pub tables: &'a [TemperatureTwoPhotonTable],
    /// Exactly one nonlinear assignment for each non-vacuum original path cell.
    pub cells: &'a [ThermalTwoPhotonCell],
}

#[derive(Debug)]
pub(super) struct CompiledTwoPhoton {
    area_m2: f64,
    tables: Vec<TemperatureTwoPhotonTable>,
    path_tables: Vec<usize>,
}

impl TemperatureOpticalHeating<1> {
    /// Prepare temperature-dependent linear and degenerate two-photon heating.
    ///
    /// This specialization accepts exactly one physical wavelength. Separate
    /// packets or overlapping beams cannot be added to synthesize a nonlinear
    /// interaction. Powers are already inside the medium. The original thermal
    /// owner and `advance_refined` integrate all heat and changing shielding.
    /// Combined linear/nonlinear table and knot counts obey the original budget;
    /// endpoint irradiance checks also consume its source-evaluation allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_two_photon(
        cx: &Cx<'_>,
        grid: &AbsorbingGrid<'_>,
        beam: IncidentBeam<1>,
        linear_tables: &[TemperatureAbsorptionTable],
        bindings: &[ThermalOpticalCell],
        nonlinear: TwoPhotonFeedback<'_>,
        thermal: &DepositionSimulation<'_>,
        budget: OpticalFeedbackBudget,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        admit(
            linear_tables.len(),
            linear_tables
                .iter()
                .map(TemperatureAbsorptionTable::knot_count),
            nonlinear,
            budget,
        )?;
        let result = Self::try_new(
            cx,
            grid,
            beam,
            linear_tables,
            bindings,
            thermal,
            budget,
            power_tolerance_w,
            &mut checkpoint,
        )?;
        Self::attach_two_photon(result, nonlinear, thermal, &mut checkpoint)
    }

    pub(super) fn attach_two_photon(
        mut result: Self,
        nonlinear: TwoPhotonFeedback<'_>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        // This is a newly constructed, unshared model: no accepted state or
        // candidate is mutated while adding the nonlinear constitutive owner.
        let model = Arc::get_mut(&mut result.model)
            .ok_or_else(|| invalid("unshared original optical preparation required"))?;
        let mut selected = BTreeMap::new();
        for cell in nonlinear.cells {
            checkpoint()?;
            if cell.two_photon_table >= nonlinear.tables.len() {
                return Err(invalid("unknown nonlinear absorption table"));
            }
            if selected
                .insert(cell.grid_cell, cell.two_photon_table)
                .is_some()
            {
                return Err(OpticalHeatingError::DuplicateCell(cell.grid_cell).into());
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-two-photon-beam.v1");
        hash.update(model.identity.as_bytes());
        hash.update(&nonlinear.beam_area_m2.to_bits().to_le_bytes());
        hash.update(&(nonlinear.tables.len() as u64).to_le_bytes());
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(nonlinear.tables.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        for table in nonlinear.tables {
            checkpoint()?;
            hash.update(table.0.identity.as_bytes());
            model.wavelength_bounds_nm.0 =
                model.wavelength_bounds_nm.0.max(table.0.wavelengths_nm.0);
            model.wavelength_bounds_nm.1 =
                model.wavelength_bounds_nm.1.min(table.0.wavelengths_nm.1);
            retained.push(table.clone());
        }
        let (lo, hi) = model.wavelength_bounds_nm;
        if lo >= hi || !(lo..=hi).contains(&model.beam.wavelengths_nm[0]) {
            return Err(VolumeTransportError::Wavelength.into());
        }
        let mut path_tables = Vec::new();
        path_tables
            .try_reserve_exact(model.path.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        for segment in &model.path {
            checkpoint()?;
            let table = selected
                .remove(&segment.path.cell_index)
                .ok_or(OpticalHeatingError::UnmappedCell(segment.path.cell_index))?;
            hash.update(&(segment.path.cell_index as u64).to_le_bytes());
            hash.update(&(table as u64).to_le_bytes());
            path_tables.push(table);
        }
        if !selected.is_empty() {
            return Err(invalid(
                "nonlinear binding lies outside the original non-vacuum path",
            ));
        }
        model.identity = hash.finalize();
        model.two_photon = Some(CompiledTwoPhoton {
            area_m2: nonlinear.beam_area_m2,
            tables: retained,
            path_tables,
        });
        let (power, sources) = model.evaluate(thermal, checkpoint)?;
        result.power = power;
        result.sources_w = sources;
        Ok(result)
    }
}

// Shared by the temperature and enthalpy linear-coordinate constructors.
pub(super) fn admit(
    linear_tables: usize,
    linear_knots: impl Iterator<Item = usize>,
    nonlinear: TwoPhotonFeedback<'_>,
    budget: OpticalFeedbackBudget,
) -> Result<(), OpticalFeedbackError> {
    if !nonlinear.beam_area_m2.is_finite()
        || nonlinear.beam_area_m2 <= 0.0
        || nonlinear.tables.is_empty()
    {
        return Err(invalid(
            "positive finite beam area and explicit nonlinear source tables required",
        ));
    }
    let total_tables = linear_tables
        .checked_add(nonlinear.tables.len())
        .ok_or(OpticalFeedbackError::Budget)?;
    let total_knots = linear_knots
        .chain(nonlinear.tables.iter().map(|t| t.0.rows.len()))
        .try_fold(0usize, |sum, n| {
            sum.checked_add(n).ok_or(OpticalFeedbackError::Budget)
        })?;
    if total_tables > budget.max_tables
        || total_knots > budget.max_total_temperature_knots
        || nonlinear.cells.len() > budget.max_path_cells
    {
        return Err(OpticalFeedbackError::Budget);
    }
    Ok(())
}

impl CompiledTwoPhoton {
    pub(super) fn absorb(
        &self,
        path_index: usize,
        segment: &BoundSegment,
        temperature_k: f64,
        power_w: f64,
        kappa_m2_kg: f64,
    ) -> Result<(f64, f64), OpticalFeedbackError> {
        let table = &self.tables[self.path_tables[path_index]];
        let beta = table.at(segment.thermal_cell, temperature_k)?;
        let irradiance = power_w / self.area_m2;
        if irradiance > table.0.maximum_irradiance_w_m2 {
            return Err(OpticalFeedbackError::IrradianceDomain {
                grid_cell: segment.path.cell_index,
                irradiance_w_m2: irradiance,
                maximum_w_m2: table.0.maximum_irradiance_w_m2,
            });
        }
        let (absorbed, transmitted) = absorb_two_photon(
            power_w,
            self.area_m2,
            segment.path.density_kg_m3,
            segment.path.length_m,
            kappa_m2_kg,
            beta,
        );
        if !nonnegative(absorbed) || !nonnegative(transmitted) {
            return Err(VolumeTransportError::Unrepresentable.into());
        }
        Ok((absorbed, transmitted))
    }
}
