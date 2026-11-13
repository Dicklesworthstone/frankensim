//! Explicit temperature feedback for coherent finite-footprint cell heating.
//!
//! Every original cell refreshes its complete Debye, Lorentz or Ohmic bundle
//! at its accepted native temperature. The original mean and real/proper factor
//! banks then rebuild one common two-face field across ALL those cells at each
//! authored frequency. Spatial/frequency interference is never replaced by
//! independent line powers or a layer-average optical temperature.
//!
//! One interval holds these refreshed coefficients and consumes the existing
//! spatial heat step exactly once. Its endpoint must remain within the declared
//! per-cell temperature excursion and admit every original endpoint source and
//! frequency stack before publication. The gate controls an explicit split; it
//! is not a coupled residual or a temporal-error certificate. The unchanged DD
//! clock advances only by the actual accepted native thermal duration.
//! The explicit field-heat advance preserves the original coherent reservoir
//! fields through that interval's volume integral, with separate query limits.
//!
//! Geometry, exterior source states and native thermal laws stay with their
//! original owners. The footprint remains a through-thickness thermal column,
//! with no resolved lateral temperature or patch-edge scattering. This does
//! not infer a causal reservoir model from visible/spectral n+k tables, refresh
//! transient polarization history, transfer parameter work, or add phase change.

use super::*;
use crate::layered_wave_heating::{
    LayeredWaveHeatingError,
    interference::{
        RefreshedFilmLaw, ThermalInterferenceDrive, ThermalInterferenceError,
        ThermalInterferenceLimits, ThermalInterferenceLine, Workspace, admit_drive,
        admit_refresh_source, bytes, checked_sum, reserve as reserve_sources, with_drive_ensemble,
    },
};
use fs_em::plane_wave::{
    DielectricResponse,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumFactorStatistics, AngularSpectrumGrid,
        PlanarAngularSpectrum, PreparedAngularSpectrumIllumination,
    },
    layers::PlanarDielectricLayers,
};
use fs_matdb::MaterialCard;

/// Complete original field bounds and a separate native-cell count bound.
/// fields.maximum_preparation_workspace_bytes bounds all coexisting numerical
/// preparation, including retained cell sources. Source payloads and allocator
/// overhead are excluded, as in the original interference preparation owner.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceFeedbackLimits {
    pub maximum_cells: usize,
    pub fields: ThermalInterferenceLimits,
}

#[derive(Debug)]
pub enum SpatialInterferenceFeedbackError {
    Invalid(&'static str),
    Sources(ThermalInterferenceError),
    Binding(LayeredWaveHeatingError),
    Spatial(SpatialInterferenceHeatingError),
    Thermal(DepositionError),
    TemperatureExcursion {
        cell: CellId,
        change_k: f64,
        maximum_k: f64,
    },
}
impl core::fmt::Display for SpatialInterferenceFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial interference feedback: {self:?}")
    }
}
impl std::error::Error for SpatialInterferenceFeedbackError {}
impl From<ThermalInterferenceError> for SpatialInterferenceFeedbackError {
    fn from(error: ThermalInterferenceError) -> Self {
        Self::Sources(error)
    }
}
impl From<LayeredWaveHeatingError> for SpatialInterferenceFeedbackError {
    fn from(error: LayeredWaveHeatingError) -> Self {
        Self::Binding(error)
    }
}
impl From<SpatialInterferenceHeatingError> for SpatialInterferenceFeedbackError {
    fn from(error: SpatialInterferenceHeatingError) -> Self {
        Self::Spatial(error)
    }
}
impl From<DepositionError> for SpatialInterferenceFeedbackError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Fixed geometry, source selections and coherent drive. Cards are in ORIGINAL
/// material-layer order; their original responses already belong to the
/// LayeredWaveHeating owner. Each native subdivision refreshes independently.
pub struct QuasiSteadySpatialInterferenceHeating<'h, 'b, 'a, 'm> {
    original: &'h LayeredWaveHeating<'b, 'a, 'm>,
    cards: &'h [&'h MaterialCard],
    grid: AngularSpectrumGrid,
    carrier_rad_m: [f64; 2],
    drive: ThermalInterferenceDrive<'h>,
    quadrature: &'h [SpatialInterferenceCellQuadrature<'h>],
    footprint: AngularSpectrumAperture,
    limits: SpatialInterferenceFeedbackLimits,
    frequencies: Vec<f64>,
    anchor_index: usize,
    line_count: usize,
    metadata_bytes: usize,
}

/// A jointly admitted native thermal endpoint and original material snapshots.
/// Applied snapshots generated this interval's heat; endpoint snapshots are
/// admission for subsequent rebuilding, never a second heat credit.
pub struct SpatialInterferenceFeedbackStep<'c, 't> {
    pub step: SpatialInterferenceStep<'c, 't>,
    pub applied_temperatures_k: Vec<f64>,
    pub endpoint_temperatures_k: Vec<f64>,
    pub maximum_temperature_excursion_k: f64,
    pub preparation_workspace_bytes: usize,
    applied_sources: Vec<RefreshedFilmLaw>,
    endpoint_sources: Vec<RefreshedFilmLaw>,
}
impl SpatialInterferenceFeedbackStep<'_, '_> {
    /// Original complete source bundle, indexed by native cell order.
    pub fn applied_response(&self, cell_index: usize) -> Option<DielectricResponse<'_>> {
        self.applied_sources
            .get(cell_index)
            .map(RefreshedFilmLaw::response)
    }
    pub fn endpoint_response(&self, cell_index: usize) -> Option<DielectricResponse<'_>> {
        self.endpoint_sources
            .get(cell_index)
            .map(RefreshedFilmLaw::response)
    }
}

impl<'h, 'b, 'a, 'm> QuasiSteadySpatialInterferenceHeating<'h, 'b, 'a, 'm> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        original: &'h LayeredWaveHeating<'b, 'a, 'm>,
        cards: &'h [&'h MaterialCard],
        grid: AngularSpectrumGrid,
        drive: ThermalInterferenceDrive<'h>,
        quadrature: &'h [SpatialInterferenceCellQuadrature<'h>],
        footprint: AngularSpectrumAperture,
        limits: SpatialInterferenceFeedbackLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SpatialInterferenceFeedbackError> {
        Self::try_new_with_carrier(
            original, cards, grid, [0.0; 2], drive, quadrature, footprint, limits, cancelled,
        )
    }

    /// Heat the original cells with a common physical transverse carrier
    /// [rad/m]. Both face arrays and every frequency/factor are envelopes
    /// with exp(-i k_carrier.x) removed. Their exact carrier phase reaches the
    /// finite footprint's spatial interference and irreversible heat.
    ///
    /// The authored wavevector is fixed across frequencies and per-cell
    /// temperature rebuilds; incidence angles follow each original exterior
    /// wavelength. The native two-port light-cone admission remains in force.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_carrier(
        original: &'h LayeredWaveHeating<'b, 'a, 'm>,
        cards: &'h [&'h MaterialCard],
        grid: AngularSpectrumGrid,
        carrier_rad_m: [f64; 2],
        drive: ThermalInterferenceDrive<'h>,
        quadrature: &'h [SpatialInterferenceCellQuadrature<'h>],
        footprint: AngularSpectrumAperture,
        limits: SpatialInterferenceFeedbackLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SpatialInterferenceFeedbackError> {
        poll(&mut cancelled)?;
        if carrier_rad_m.iter().any(|k| !k.is_finite()) {
            return Err(invalid("the transverse carrier must be finite in rad/m"));
        }
        let body = original.thermal_body().prepared_body();
        let n = body.cell_count();
        if n == 0
            || n > limits.maximum_cells
            || quadrature.len() != n
            || cards.len() != body.conduction().layers().len()
        {
            return Err(invalid(
                "feedback requires original layer cards and one quadrature per bounded native cell",
            ));
        }
        let area = footprint
            .extent_m()
            .into_iter()
            .fold(Dd::ONE, |area, length| area * Dd::from_f64(length));
        match_measure(area, Dd::from_f64(original.area_m2()))?;
        let line_count = admit_drive(drive, limits.fields, &mut cancelled)?;
        let points = grid
            .width
            .checked_mul(grid.height)
            .ok_or(invalid("original field grid size overflows"))?;
        if points == 0
            || points > limits.fields.angular.maximum_grid_points
            || grid
                .spacing_m
                .iter()
                .any(|dx| !dx.is_finite() || *dx <= 0.0)
        {
            return Err(invalid(
                "a bounded original field grid with positive finite spacings is required",
            ));
        }
        let mut cell = 0;
        for (layer_index, (layer, &card)) in
            body.conduction().layers().iter().zip(cards).enumerate()
        {
            for _ in 0..body.layer_cell_counts()[layer_index] {
                poll(&mut cancelled)?;
                let film = original
                    .optical_layers()
                    .layer(cell)
                    .ok_or(invalid("original cell film is absent"))?;
                admit_refresh_source(film.response, card, layer.temperature_axis())?;
                cell += 1;
            }
        }
        let metadata_bytes = checked_sum(&[
            bytes::<&ThermalInterferenceLine<'_>>(line_count)?,
            bytes::<f64>(line_count)?,
            bytes::<&PlanarDielectricLayers<'_>>(line_count)?,
            bytes::<LayeredWaveHeating<'_, '_, '_>>(line_count)?,
            bytes::<PlanarAngularSpectrum<'_, '_>>(line_count)?,
            bytes::<PreparedAngularSpectrumIllumination<'_, '_, '_>>(line_count)?,
            bytes::<AngularSpectrumExposure<'_, '_, '_>>(drive.factors.len() + 1)?,
            bytes::<AngularSpectrumFactorStatistics>(drive.factors.len())?,
            bytes::<DielectricResponse<'_>>(n)?,
            bytes::<RefreshedFilmLaw>(2 * n)?,
            bytes::<f64>(2 * n)?,
            bytes::<CellId>(
                line_count
                    .checked_mul(n)
                    .ok_or(invalid("cell source metadata overflows"))?,
            )?,
        ])?;
        Workspace::new(
            metadata_bytes,
            limits.fields.maximum_preparation_workspace_bytes,
        )?;
        let mut frequencies = reserve_sources(line_count)?;
        let mut anchor_frequency = None;
        for line in drive
            .mean
            .iter()
            .chain(drive.factors.iter().flat_map(|factor| factor.lines))
        {
            poll(&mut cancelled)?;
            if line.left_tangential_e_v_per_m.len() != points
                || line.right_tangential_e_v_per_m.len() != points
            {
                return Err(invalid(
                    "both physical faces of every source line must retain the original grid",
                ));
            }
            anchor_frequency.get_or_insert(line.angular_frequency_rad_s);
            frequencies.push(line.angular_frequency_rad_s);
        }
        frequencies.sort_unstable_by(f64::total_cmp);
        frequencies.dedup_by(|a, b| *a == *b);
        let anchor_frequency =
            anchor_frequency.ok_or(invalid("an original source line is required"))?;
        let anchor_index = frequencies
            .binary_search_by(|omega| omega.total_cmp(&anchor_frequency))
            .map_err(|_| invalid("original anchor frequency is absent"))?;
        poll(&mut cancelled)?;
        Ok(Self {
            original,
            cards,
            grid,
            carrier_rad_m,
            drive,
            quadrature,
            footprint,
            limits,
            frequencies,
            anchor_index,
            line_count,
            metadata_bytes,
        })
    }

    pub const fn original_owner(&self) -> &'h LayeredWaveHeating<'b, 'a, 'm> {
        self.original
    }
    pub const fn drive(&self) -> ThermalInterferenceDrive<'h> {
        self.drive
    }
    /// Authored common wavevector, unchanged across frequencies and refreshes.
    pub const fn carrier_transverse_wavenumber_rad_m(&self) -> [f64; 2] {
        self.carrier_rad_m
    }
    pub const fn footprint(&self) -> AngularSpectrumAperture {
        self.footprint
    }
    pub fn frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }

    fn workspace(&self) -> Result<Workspace, ThermalInterferenceError> {
        Workspace::new(
            self.metadata_bytes,
            self.limits.fields.maximum_preparation_workspace_bytes,
        )
    }

    fn refresh_sources(
        &self,
        temperatures: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<RefreshedFilmLaw>, SpatialInterferenceFeedbackError> {
        let body = self.original.thermal_body().prepared_body();
        if temperatures.len() != body.cell_count() {
            return Err(invalid(
                "one accepted temperature per original cell is required",
            ));
        }
        let mut sources = reserve_sources(temperatures.len())?;
        for (index, (layer, &card)) in body
            .conduction()
            .layers()
            .iter()
            .zip(self.cards)
            .enumerate()
        {
            for _ in 0..body.layer_cell_counts()[index] {
                poll(&mut *cancelled)?;
                let cell = sources.len();
                let original = self
                    .original
                    .optical_layers()
                    .layer(cell)
                    .ok_or(invalid("original cell film is absent"))?
                    .response;
                sources.push(RefreshedFilmLaw::at_temperature(
                    original,
                    card,
                    layer.temperature_axis(),
                    temperatures[cell],
                )?);
            }
        }
        poll(&mut *cancelled)?;
        Ok(sources)
    }

    fn frequency_owners<'s>(
        &'s self,
        sources: &'s [RefreshedFilmLaw],
        work: &mut Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<LayeredWaveHeating<'b, 'a, 's>>, SpatialInterferenceFeedbackError> {
        // Geometry/cell count are unchanged, so the original stack's complete
        // cache size bounds each rebuilt stack before constructing any of them.
        work.charge(
            self.original
                .optical_layers()
                .workspace_bytes()
                .checked_mul(self.frequencies.len())
                .ok_or(invalid("frequency cache storage overflows"))?,
        )?;
        let mut responses = reserve_sources(sources.len())?;
        responses.extend(sources.iter().map(RefreshedFilmLaw::response));
        let mut owners = reserve_sources(self.frequencies.len())?;
        for &omega in &self.frequencies {
            poll(&mut *cancelled)?;
            owners.push(LayeredWaveHeating::try_new_cellwise_with_incident_response(
                self.original.thermal_body(),
                self.original.optical_layers().incident_response(),
                self.original.optical_layers().exit_response(),
                omega,
                &responses,
                &mut *cancelled,
            )?);
        }
        Ok(owners)
    }

    /// Admit an existing native thermal state and declare its source epoch.
    /// Complete source refresh and every actual-frequency Maxwell stack are
    /// checked without depositing heat or changing the accepted timeline.
    pub fn initial_state<'t>(
        &self,
        thermal: DepositionSimulation<'t>,
        source_time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceState<'h, 't>, SpatialInterferenceFeedbackError> {
        poll(&mut cancelled)?;
        let source_time_s = checked_clock(source_time_s)?;
        let temperatures = self
            .original
            .thermal_body()
            .cell_temperatures(&thermal, || checkpoint(&mut cancelled))?;
        let mut work = self.workspace()?;
        let sources = self.refresh_sources(&temperatures, &mut cancelled)?;
        self.frequency_owners(&sources, &mut work, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(SpatialInterferenceState {
            thermal,
            source_clock: self.drive.source_clock,
            source_time_s,
        })
    }

    /// Rebuild at accepted per-cell temperatures, hold those laws over one
    /// original implicit thermal step, then admit all endpoint material laws.
    /// The excursion applies to EACH cell, not a mean-temperature change.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &SpatialInterferenceState<'_, 't>,
        target_s: f64,
        maximum_temperature_change_k: f64,
        temporal: SpatialInterferenceHeatingControls,
        conduction: LayerConductionConfig,
        solve: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceFeedbackStep<'h, 't>, SpatialInterferenceFeedbackError> {
        self.advance_to_impl(
            initial,
            target_s,
            maximum_temperature_change_k,
            temporal,
            None,
            conduction,
            solve,
            workspace,
            cancelled,
        )
    }

    /// Refresh every original cell's complete optical law at its accepted
    /// temperature, then integrate its coherent reservoir fields over the
    /// actual native thermal interval before reducing them to joules. The
    /// accepted DD source epoch, physical carrier and original continuous
    /// footprint enter that field integral unchanged.
    ///
    /// `field_limits` separately bound each native cell's explicit field query,
    /// including its original layer workspace; preparation remains subject to
    /// `SpatialInterferenceFeedbackLimits`. Both advance routes share the
    /// original ripple gate, single heat debit, temperature-excursion gate and
    /// complete endpoint material/frequency admission. Endpoint snapshots do
    /// not produce another heat credit. A refusal leaves the accepted state
    /// and source clock unchanged. Optics remain held within the interval;
    /// this does not transport reservoir transients or parameter work across
    /// the temperature refresh.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_field_heat<'t>(
        &self,
        initial: &SpatialInterferenceState<'_, 't>,
        target_s: f64,
        maximum_temperature_change_k: f64,
        temporal: SpatialInterferenceHeatingControls,
        field_limits: AngularSpectrumHeatWindowLimits,
        conduction: LayerConductionConfig,
        solve: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceFeedbackStep<'h, 't>, SpatialInterferenceFeedbackError> {
        self.advance_to_impl(
            initial,
            target_s,
            maximum_temperature_change_k,
            temporal,
            Some(field_limits),
            conduction,
            solve,
            workspace,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_to_impl<'t>(
        &self,
        initial: &SpatialInterferenceState<'_, 't>,
        target_s: f64,
        maximum_temperature_change_k: f64,
        temporal: SpatialInterferenceHeatingControls,
        field_limits: Option<AngularSpectrumHeatWindowLimits>,
        conduction: LayerConductionConfig,
        solve: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceFeedbackStep<'h, 't>, SpatialInterferenceFeedbackError> {
        poll(&mut cancelled)?;
        if initial.source_clock != self.drive.source_clock
            || !maximum_temperature_change_k.is_finite()
            || maximum_temperature_change_k < 0.0
        {
            return Err(invalid(
                "accepted source clock must match and temperature excursion must be finite and nonnegative",
            ));
        }
        let applied_temperatures_k = self
            .original
            .thermal_body()
            .cell_temperatures(&initial.thermal, || checkpoint(&mut cancelled))?;
        let mut work = self.workspace()?;
        let applied_sources = self.refresh_sources(&applied_temperatures_k, &mut cancelled)?;
        let step = {
            let owners = self.frequency_owners(&applied_sources, &mut work, &mut cancelled)?;
            let mut stacks = reserve_sources(owners.len())?;
            stacks.extend(owners.iter().map(LayeredWaveHeating::optical_layers));
            with_drive_ensemble(
                &stacks,
                &self.frequencies,
                self.grid,
                self.carrier_rad_m,
                self.drive,
                self.line_count,
                self.limits.fields,
                &mut work,
                &mut cancelled,
                |ensemble, work, cancelled| -> Result<_, SpatialInterferenceFeedbackError> {
                    let heating = SpatialInterferenceHeating::try_new(
                        &owners[self.anchor_index],
                        ensemble,
                        self.quadrature,
                        self.footprint,
                        self.drive.source_clock,
                        SpatialInterferenceHeatingLimits {
                            maximum_cells: self.limits.maximum_cells,
                            maximum_preparation_workspace_bytes: work.remaining(),
                            layer: self.limits.fields.layer,
                            cell_heat: self.limits.fields.heat,
                        },
                        &mut *cancelled,
                    )?;
                    work.charge(heating.preparation_workspace_bytes())?;
                    let step = if let Some(limits) = field_limits {
                        heating.advance_to_with_field_heat(
                            initial,
                            target_s,
                            temporal,
                            limits,
                            conduction,
                            solve,
                            workspace,
                            &mut *cancelled,
                        )?
                    } else {
                        heating.advance_to(
                            initial,
                            target_s,
                            temporal,
                            conduction,
                            solve,
                            workspace,
                            &mut *cancelled,
                        )?
                    };
                    // Only the clock label borrow changes owner. Its exact DD
                    // value and the original accepted native state move intact.
                    Ok(retain_clock(step, self.drive.source_clock))
                },
            )??
        };
        let endpoint_temperatures_k = self
            .original
            .thermal_body()
            .cell_temperatures(step.state.thermal(), || checkpoint(&mut cancelled))?;
        let maximum_temperature_excursion_k = temperature_excursion(
            self.original.cell_ids(),
            &applied_temperatures_k,
            &endpoint_temperatures_k,
            maximum_temperature_change_k,
            &mut cancelled,
        )?;
        let endpoint_sources = self.refresh_sources(&endpoint_temperatures_k, &mut cancelled)?;
        // Applied field caches have been dropped. Admit endpoint laws/frequency
        // support with a fresh bounded workspace, without a second source credit.
        let mut endpoint_work = self.workspace()?;
        self.frequency_owners(&endpoint_sources, &mut endpoint_work, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(SpatialInterferenceFeedbackStep {
            step,
            applied_temperatures_k,
            endpoint_temperatures_k,
            maximum_temperature_excursion_k,
            preparation_workspace_bytes: work.used().max(endpoint_work.used()),
            applied_sources,
            endpoint_sources,
        })
    }
}

fn retain_clock<'c, 't>(
    step: SpatialInterferenceStep<'_, 't>,
    clock: &'c str,
) -> SpatialInterferenceStep<'c, 't> {
    SpatialInterferenceStep {
        state: SpatialInterferenceState {
            thermal: step.state.thermal,
            source_clock: clock,
            source_time_s: step.state.source_time_s,
        },
        heat_report: step.heat_report,
        source_start_s: step.source_start_s,
        duration_s: step.duration_s,
        absorbed_energy_j: step.absorbed_energy_j,
        cell_energy_j: step.cell_energy_j,
        source_ripple_bound_j: step.source_ripple_bound_j,
    }
}

fn temperature_excursion(
    cells: &[CellId],
    initial: &[f64],
    endpoint: &[f64],
    maximum: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, SpatialInterferenceFeedbackError> {
    if cells.len() != initial.len()
        || initial.len() != endpoint.len()
        || !maximum.is_finite()
        || maximum < 0.0
    {
        return Err(invalid(
            "temperature excursion needs matching original cells and a finite nonnegative limit",
        ));
    }
    let mut largest = 0.0_f64;
    for ((&cell, &before), &after) in cells.iter().zip(initial).zip(endpoint) {
        poll(&mut *cancelled)?;
        let change = (after - before).abs();
        if !change.is_finite() || change > maximum {
            return Err(SpatialInterferenceFeedbackError::TemperatureExcursion {
                cell,
                change_k: change,
                maximum_k: maximum,
            });
        }
        largest = largest.max(change);
    }
    Ok(largest)
}

fn invalid(message: &'static str) -> SpatialInterferenceFeedbackError {
    SpatialInterferenceFeedbackError::Invalid(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn opposite_cell_changes_do_not_hide_temperature_excursion() {
        let cells = [CellId(41), CellId(42)];
        // The mean is unchanged, while each original cell moves by one kelvin.
        let start = [300.0, 320.0];
        let end = [301.0, 319.0];
        assert!(matches!(
            temperature_excursion(&cells, &start, &end, 0.5, &mut || false),
            Err(SpatialInterferenceFeedbackError::TemperatureExcursion {
                cell: CellId(41),
                change_k: 1.0,
                maximum_k: 0.5,
            })
        ));
        assert_eq!(
            temperature_excursion(&cells, &start, &end, 1.0, &mut || false).unwrap(),
            1.0
        );
        assert_eq!(
            temperature_excursion(&cells, &start, &start, 0.0, &mut || false).unwrap(),
            0.0
        );
    }

    #[test]
    fn temperature_gate_refuses_nonfinite_domains_and_polls_each_original_cell() {
        let cells = [CellId(41), CellId(42)];
        for maximum in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                temperature_excursion(&cells, &[300.0; 2], &[300.0; 2], maximum, &mut || false)
                    .is_err()
            );
        }
        assert!(
            temperature_excursion(&cells, &[300.0; 2], &[f64::NAN, 300.0], 1.0, &mut || false)
                .is_err()
        );
        assert!(temperature_excursion(&cells, &[300.0], &[300.0; 2], 1.0, &mut || false).is_err());
        let mut polls = 0;
        assert!(
            temperature_excursion(&cells, &[300.0; 2], &[300.0; 2], 1.0, &mut || {
                polls += 1;
                polls == 2
            })
            .is_err()
        );
        assert_eq!(polls, 2);
    }
}
