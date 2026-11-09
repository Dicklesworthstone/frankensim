//! Depth-dependent surface photoemission on the original thermal cells.
//!
//! Native coherent loss operators integrate the intersections of thermal cells
//! and source-authored escape bins. Only residual heat enters the native thermal
//! candidate; full optical absorption, removal work, escaping kinetic energy and
//! charge transfer remain separate accounts. Other films retain their original
//! absorbed heat. No original cell, material, source or boundary is replaced.
//!
//! [`LocalThermalization`] explicitly supplies the missing spatial closure:
//! residual energy thermalizes where the photon was absorbed. Conditional
//! escape probabilities and mean exit energies do not determine an electron
//! scattering/deposition kernel. This owner therefore does not model spatial
//! electron heat transport. Optical and emission properties are held over the
//! interval; charge is reported without advancing an electrical history,
//! collector, contact potential, space charge or field-dependent escape law.

use super::{
    LayeredWaveDrive, LayeredWaveEnergy, LayeredWaveHeatingError, PreparedLayeredWaveHeating,
    checkpoint, poll,
};
use fs_em::plane_wave::layers::{
    LayerIlluminationPower, MAX_LAYERS, PreparedLayerDepthInterval, photoelectric::EscapeFace,
};
use fs_material::{
    photoelectric::{
        ELEMENTARY_CHARGE_C, LIGHT_M_S, PhotoelectricError, PhotoelectricLaw,
        source::{
            PhotoelectricSourceError,
            depth::{DepthPhotoelectricMaterial, MAX_PHOTOELECTRIC_DEPTH_BINS},
        },
    },
    thermal::{
        conduction::LayerConductionConfig,
        transport::stack::transient::implicit::{
            LayeredThermalSolveControls, LayeredThermalSolveReport, LayeredThermalStep,
            LayeredThermalWorkspace,
        },
    },
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

/// Explicit spatial approximation, never inferred from a quantum yield.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocalThermalization {
    /// Energy remaining after removal work and escaping kinetic energy becomes
    /// heat at the photon's absorption position.
    AtAbsorptionPosition,
}

#[derive(Debug)]
pub enum LayeredPhotoelectricError {
    Invalid(&'static str),
    Heating(LayeredWaveHeatingError),
    Source(PhotoelectricSourceError),
    Emission(PhotoelectricError),
}
impl core::fmt::Display for LayeredPhotoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "layered photoelectric heat: {self:?}")
    }
}
impl std::error::Error for LayeredPhotoelectricError {}
impl From<LayeredWaveHeatingError> for LayeredPhotoelectricError {
    fn from(value: LayeredWaveHeatingError) -> Self {
        Self::Heating(value)
    }
}
impl From<PhotoelectricSourceError> for LayeredPhotoelectricError {
    fn from(value: PhotoelectricSourceError) -> Self {
        Self::Source(value)
    }
}
impl From<PhotoelectricError> for LayeredPhotoelectricError {
    fn from(value: PhotoelectricError) -> Self {
        Self::Emission(value)
    }
}

struct CellEscapeInterval<'q, 'p, 'm> {
    cell: usize,
    absorption: PreparedLayerDepthInterval<'q, 'p, 'm>,
    law: PhotoelectricLaw,
}

/// One exposed original film, with at most cells + bins - 1 cached intervals.
/// Optical cells are retained even when a source depth boundary cuts a cell.
pub struct PreparedLayeredPhotoelectricHeating<'q, 'p, 'b, 'a, 'm, 'd> {
    wave: &'q PreparedLayeredWaveHeating<'p, 'b, 'a, 'm>,
    source: &'d DepthPhotoelectricMaterial,
    layer_index: usize,
    first_cell: usize,
    cell_count: usize,
    face: EscapeFace,
    thermalization: LocalThermalization,
    wavelength_m: f64,
    intervals: Vec<CellEscapeInterval<'q, 'p, 'm>>,
}

/// Transfers over exactly the returned native thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoelectricEnergy {
    /// Original absorption in the emitting physical film, before partitioning.
    pub absorbed_emitter_j: f64,
    /// Additional heat deposited across ALL original cells, including other films.
    pub residual_heat_j: f64,
    pub removal_j: f64,
    pub kinetic_j: f64,
    pub emitted_electrons: f64,
    /// Positive emitter charging and equal negative escaping-electron charge.
    pub surface_charge_c: f64,
    pub emitted_charge_c: f64,
    /// Original total absorption minus heat, removal work and escaping kinetics.
    pub partition_residual_j: f64,
}

/// Joint candidate; a failed or cancelled solve publishes no history or transfer.
pub struct LayeredPhotoelectricStep<'t, 'd> {
    pub thermal: DepositionSimulation<'t>,
    pub heat_report: LayeredThermalSolveReport,
    /// Full optical absorption; only photoelectric.residual_heat_j enters heat.
    pub optical: LayeredWaveEnergy,
    pub optical_power: LayerIlluminationPower,
    pub photoelectric: PhotoelectricEnergy,
    pub source: &'d DepthPhotoelectricMaterial,
    pub emitting_layer_index: usize,
    pub escape_face: EscapeFace,
    pub thermalization: LocalThermalization,
    /// Residual optical heat once per original native CellId, in original order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

impl<'q, 'p, 'b, 'a, 'm, 'd> PreparedLayeredPhotoelectricHeating<'q, 'p, 'b, 'a, 'm, 'd> {
    /// Bind source depth measured from the explicitly declared vacuum face.
    /// The calibration thickness describes the ORIGINAL physical film, not an
    /// individual optical/thermal subcell. Every emitting cell must retain the
    /// complete same source snapshot; cellwise differing source states refuse.
    pub fn try_new(
        wave: &'q PreparedLayeredWaveHeating<'p, 'b, 'a, 'm>,
        layer_index: usize,
        source: &'d DepthPhotoelectricMaterial,
        face: EscapeFace,
        thermalization: LocalThermalization,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredPhotoelectricError> {
        poll(&mut cancelled)?;
        let invalid = LayeredPhotoelectricError::Invalid;
        let body = wave.owner.thermal.prepared_body();
        let layers = body.conduction().layers();
        let layer = layers
            .get(layer_index)
            .ok_or(invalid("emitting original film is out of range"))?;
        let exposed = match face {
            EscapeFace::EntranceVacuum => layer_index == 0,
            EscapeFace::ExitVacuum => layer_index + 1 == layers.len(),
        };
        if !exposed || layer.thickness_m() != source.thickness_m() {
            return Err(invalid(
                "emitting film must be exposed and retain its calibrated original thickness",
            ));
        }
        let first_cell: usize = body.layer_cell_counts()[..layer_index].iter().sum();
        let cell_count = body.layer_cell_counts()[layer_index];
        for cell in first_cell..first_cell + cell_count {
            let film = wave
                .owner
                .optical
                .layer(cell)
                .ok_or(invalid("original optical cell is missing"))?;
            if film.response.material() != source.resolved() {
                return Err(invalid(
                    "emission and every emitting optical cell require the same complete source",
                ));
            }
        }
        let wavelength_m =
            core::f64::consts::TAU * LIGHT_M_S / wave.owner.optical.angular_frequency_rad_s();
        let wavelength_nm = checked_scale(wavelength_m, 1e9)?;
        let capacity = cell_count + source.bins().len() - 1;
        if capacity > MAX_LAYERS + MAX_PHOTOELECTRIC_DEPTH_BINS - 1 {
            return Err(invalid(
                "cell/depth intersection count exceeds native bounds",
            ));
        }
        let mut intervals = Vec::with_capacity(capacity);
        for bin in 0..source.bins().len() {
            let law = source.at_wavelength(bin, wavelength_nm)?;
            let [a, b] = source.bin_bounds_fraction(bin)?;
            let [a, b] = match face {
                EscapeFace::EntranceVacuum => [a, b],
                EscapeFace::ExitVacuum => [1.0 - b, 1.0 - a],
            }
            .map(|fraction| fraction * cell_count as f64);
            if b <= a {
                return Err(invalid(
                    "positive source depth bin collapsed in native coordinates",
                ));
            }
            for local in 0..cell_count {
                let lo = a.max(local as f64);
                let hi = b.min((local + 1) as f64);
                if hi <= lo {
                    continue;
                }
                poll(&mut cancelled)?;
                let cell = first_cell + local;
                let width = wave.owner.optical.layer(cell).unwrap().thickness_m;
                let start = checked_scale(lo - local as f64, width)?;
                let end = checked_scale(hi - local as f64, width)?;
                if end <= start || intervals.len() == capacity {
                    return Err(invalid(
                        "positive cell/depth intersection is unrepresentable",
                    ));
                }
                let absorption = wave
                    .illumination
                    .prepare_depth_interval(cell, start, end)
                    .map_err(LayeredWaveHeatingError::from)?;
                intervals.push(CellEscapeInterval {
                    cell,
                    absorption,
                    law,
                });
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            wave,
            source,
            layer_index,
            first_cell,
            cell_count,
            face,
            thermalization,
            wavelength_m,
            intervals,
        })
    }

    /// Coherence work includes the original full-cell account AND these intervals.
    pub fn prepared_interval_count(&self) -> usize {
        self.intervals.len()
    }

    /// Advance one native thermal candidate with held illumination and local
    /// residual thermalization. For covariance factors, the supplied native
    /// budget must cover `(all_cells + prepared_interval_count) * factors.len()`.
    /// The cell/bin energy partition must meet controls.total_tolerance_j before
    /// the native nonlinear thermal owner applies its own acceptance checks.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        drive: LayeredWaveDrive<'_>,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayeredPhotoelectricStep<'t, 'd>, LayeredPhotoelectricError> {
        poll(&mut cancelled)?;
        let invalid = LayeredPhotoelectricError::Invalid;
        let n = self.wave.owner.cell_ids.len();
        if let LayeredWaveDrive::Coherence { factors, budget } = drive {
            if factors
                .len()
                .checked_mul(n + self.intervals.len())
                .is_none_or(|pairs| pairs > budget.maximum_layer_mode_pairs)
            {
                return Err(invalid(
                    "complete cell and escape-interval coherence work exceeds budget",
                ));
            }
        }
        let step = LayeredThermalStep::new(
            self.wave.owner.thermal,
            initial,
            target_s,
            conduction,
            || checkpoint(&mut cancelled),
        )
        .map_err(LayeredWaveHeatingError::from)?;
        let dt = step.duration_s();
        let mut heat = [0.0; MAX_LAYERS];
        let optical_power = self
            .wave
            .cell_power_into(drive, &mut heat[..n], &mut cancelled)?;
        let mut absorbed = Dd::ZERO;
        let mut emitter_absorbed = Dd::ZERO;
        for (cell, watts) in heat[..n].iter_mut().enumerate() {
            let energy = Dd::from_f64(checked_scale(*watts, dt)?);
            absorbed = absorbed + energy;
            if (self.first_cell..self.first_cell + self.cell_count).contains(&cell) {
                emitter_absorbed = emitter_absorbed + energy;
                *watts = 0.0;
            }
        }
        let mut local_heat = [Dd::ZERO; MAX_LAYERS];
        let mut removal = Dd::ZERO;
        let mut kinetic = Dd::ZERO;
        let mut electrons = Dd::ZERO;
        for interval in &self.intervals {
            poll(&mut cancelled)?;
            let absorption = match drive {
                LayeredWaveDrive::Coherent(input) => interval.absorption.evaluate(input),
                LayeredWaveDrive::Independent(input) => {
                    interval.absorption.evaluate_independent_sources(input)
                }
                LayeredWaveDrive::Coherence { factors, budget } => {
                    let mut interrupted = false;
                    let result =
                        interval
                            .absorption
                            .evaluate_coherence_modes(factors, budget, || {
                                let stop = cancelled();
                                interrupted |= stop;
                                stop
                            });
                    if interrupted {
                        return Err(LayeredWaveHeatingError::Cancelled.into());
                    }
                    result
                }
            }
            .map_err(LayeredWaveHeatingError::from)?;
            let watts = checked_scale(absorption, self.wave.owner.area_m2())?;
            let response = interval.law.evaluate(self.wavelength_m, watts)?;
            local_heat[interval.cell] =
                local_heat[interval.cell] + Dd::from_f64(response.heat_power_w);
            removal = removal + Dd::from_f64(response.removal_power_w);
            kinetic = kinetic + Dd::from_f64(response.kinetic_power_w);
            electrons = electrons + Dd::from_f64(response.emitted_electron_rate_s);
        }
        let mut total_heat = Dd::ZERO;
        for (cell, watts) in heat[..n].iter_mut().enumerate() {
            *watts = checked_scale(*watts + local_heat[cell].to_f64(), dt)?;
            total_heat = total_heat + Dd::from_f64(*watts);
        }
        let removal_j = checked_scale(removal.to_f64(), dt)?;
        let kinetic_j = checked_scale(kinetic.to_f64(), dt)?;
        let emitted_electrons = checked_scale(electrons.to_f64(), dt)?;
        let surface_charge_c = checked_scale(emitted_electrons, ELEMENTARY_CHARGE_C)?;
        let absorbed_j = checked_scale(absorbed.to_f64(), 1.0)?;
        let residual_heat_j = checked_scale(total_heat.to_f64(), 1.0)?;
        let partition_residual_j =
            (absorbed - total_heat - Dd::from_f64(removal_j) - Dd::from_f64(kinetic_j)).to_f64();
        if !partition_residual_j.is_finite()
            || !controls.total_tolerance_j.is_finite()
            || controls.total_tolerance_j < 0.0
            || partition_residual_j.abs() > controls.total_tolerance_j
        {
            return Err(invalid(
                "photoelectric partition exceeds the native total energy tolerance",
            ));
        }
        let mut incoming_j = [0.0; 2];
        let mut outgoing_j = [0.0; 2];
        for port in 0..2 {
            incoming_j[port] = checked_scale(
                checked_scale(optical_power.incoming_w_m2[port], self.wave.owner.area_m2())?,
                dt,
            )?;
            outgoing_j[port] = checked_scale(
                checked_scale(optical_power.outgoing_w_m2[port], self.wave.owner.area_m2())?,
                dt,
            )?;
        }
        let residual_j = (Dd::from_f64(incoming_j[0]) + Dd::from_f64(incoming_j[1])
            - Dd::from_f64(outgoing_j[0])
            - Dd::from_f64(outgoing_j[1])
            - absorbed)
            .to_f64();
        if !residual_j.is_finite() {
            return Err(invalid("optical interval account is unrepresentable"));
        }
        let cell_energy_j = step
            .cell_ids()
            .iter()
            .copied()
            .zip(heat[..n].iter().copied())
            .collect();
        let step = step
            .with_cell_energy_j(&heat[..n], || checkpoint(&mut cancelled))
            .map_err(LayeredWaveHeatingError::from)?;
        let (thermal, _, heat_report) = step
            .solve(controls, workspace, || checkpoint(&mut cancelled))
            .map_err(LayeredWaveHeatingError::from)?;
        poll(&mut cancelled)?;
        Ok(LayeredPhotoelectricStep {
            thermal,
            heat_report,
            optical: LayeredWaveEnergy {
                incoming_j,
                outgoing_j,
                absorbed_j,
                residual_j,
            },
            optical_power,
            photoelectric: PhotoelectricEnergy {
                absorbed_emitter_j: emitter_absorbed.to_f64(),
                residual_heat_j,
                removal_j,
                kinetic_j,
                emitted_electrons,
                surface_charge_c,
                emitted_charge_c: -surface_charge_c,
                partition_residual_j,
            },
            source: self.source,
            emitting_layer_index: self.layer_index,
            escape_face: self.face,
            thermalization: self.thermalization,
            cell_energy_j,
        })
    }
}

fn checked_scale(value: f64, scale: f64) -> Result<f64, LayeredPhotoelectricError> {
    let result = value * scale;
    if !result.is_finite() || value < 0.0 || scale <= 0.0 || (value > 0.0 && result == 0.0) {
        Err(LayeredPhotoelectricError::Invalid(
            "positive physical photoelectric scale is unrepresentable",
        ))
    } else {
        Ok(result)
    }
}
