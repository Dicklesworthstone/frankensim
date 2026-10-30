//! Authored signed heat power on the existing finite material-gas owner.
use super::{AperturePerformance, CoupledAperture, PlateValveInputError, Reader,
    ThermochemicalSources, GasRecord, CancelGate, MAX_CONTROLS, bad, checkpoint};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::finite_gas::{FiniteGasHeatProgram, GasHeatControl};
use crate::dg_port::finite_gas::{FiniteGasDgLoss, FiniteGasDgSolve, admitted_volume};
use fs_material::gas_species::input::{GasInputError, read_frozen_thermochemical_model};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_matdb::QueryPoint;
use std::sync::Arc;

pub(super) struct FiniteGasInput {
    model: ThermochemicalGasModel,
    point: QueryPoint,
    program: FiniteGasHeatProgram,
    line: usize,
}

impl FiniteGasInput {
    pub(super) fn read(reader: &mut Reader<'_>, sources: Option<&ThermochemicalSources<'_>>,
        gate: &CancelGate) -> Result<Self, PlateValveInputError>
    {
        let mut row = reader.row("finite_gas")?;
        let line = row.line;
        let sources = sources.ok_or_else(|| bad(line,
            "finite gas requires original thermochemical source packs and an explicit frozen-gas policy"))?;
        let (model, point) = read_frozen_thermochemical_model(&mut row.fields, sources)
            .map_err(|error| match error {
                GasInputError::Input(what) => bad(line, what),
                GasInputError::Material(source) => PlateValveInputError::Gas { line, source },
            })?;
        row.finish()?;
        let mut row = reader.row("gas_energy")?;
        let caloric = GasEnergySolve {
            minimum_temperature_k: row.scalar()?, maximum_temperature_k: row.scalar()?,
            energy_tolerance_j_kg: row.scalar()?, maximum_iterations: row.count(128)?,
        };
        row.finish()?;
        let mut row = reader.row("gas_feedback")?;
        let work_bracket_j = [row.scalar()?, row.scalar()?];
        let work_tolerance_j = row.scalar()?;
        let energy_tolerance_j = row.scalar()?;
        let maximum_iterations = row.count(128)?;
        let maximum_temperature_increment_k = row.scalar()?;
        let loss = match row.word()? {
            "lossless" => FiniteGasDgLoss::Lossless,
            "uniform-viscothermal" => FiniteGasDgLoss::UniformViscothermal { maximum_relative_error: row.scalar()? },
            "spatial-viscothermal" => FiniteGasDgLoss::SpatialViscothermal { maximum_relative_error: row.scalar()? },
            "molecular-diffusion" => FiniteGasDgLoss::UniformMolecularDiffusion { maximum_relative_error: row.scalar()? },
            _ => return Err(bad(row.line, "unknown original finite-gas loss owner")),
        };
        row.finish()?;
        let count = reader.count("gas_heat", MAX_CONTROLS)?;
        if count == 0 { return Err(bad(reader.line, "gas_heat requires an explicit initial power")); }
        let mut controls = Vec::new();
        controls.try_reserve_exact(count).map_err(|_| bad(reader.line, "cannot reserve heat controls"))?;
        for _ in 0..count {
            checkpoint(gate)?;
            let mut row = reader.row("heat")?;
            controls.push(GasHeatControl { sample: row.parse()?, power_w: row.scalar()? });
            row.finish()?;
        }
        Ok(Self { model, point, line, program: FiniteGasHeatProgram { controls, caloric,
            coupling: FiniteGasDgSolve { work_bracket_j, work_tolerance_j, energy_tolerance_j,
                maximum_iterations, maximum_temperature_increment_k, loss } } })
    }

    pub(super) fn bind(self, renderer: AperturePerformance, gate: &CancelGate,
        records: &mut Vec<GasRecord>) -> Result<AperturePerformance, PlateValveInputError>
    {
        checkpoint(gate)?;
        let CoupledAperture::Dg(system) = renderer.system() else {
            return Err(bad(self.line, "finite gas heat requires an original DG aperture"));
        };
        // Resolve the original source first, then ask the existing complete-port
        // owner for the volume of that SAME gas. The unit-volume temporary does
        // not advance or publish an inventory, and no geometry formula is copied.
        let probe = ClosedGasParcel::new(self.model, self.point, 1.0)
            .map_err(|source| PlateValveInputError::GasEnergy { line: self.line, source })?;
        let volume = admitted_volume::<GasEnergyError>(system.load(), &probe,
            &mut || gate.is_requested()).map_err(|error| PlateValveInputError::Physics(
                AcousticRealizeError::Nonlinear(error.to_string())))?;
        let gas = ClosedGasParcel::new(probe.model().clone(), probe.point().clone(), volume)
            .map_err(|source| PlateValveInputError::GasEnergy { line: self.line, source })?;
        let initial = GasRecord::Thermochemical(gas.state().clone());
        // The borrowed input gate belongs to cold construction. Playback owns a
        // separately accessible gate, so callers can retain it after moving the
        // renderer into an existing pressure/PCM consumer.
        let renderer = renderer.with_finite_gas_heat(gas, self.program,
            Arc::new(CancelGate::new_clock_free())).map_err(PlateValveInputError::Gesture)?;
        checkpoint(gate)?;
        records.push(initial);
        Ok(renderer)
    }
}
