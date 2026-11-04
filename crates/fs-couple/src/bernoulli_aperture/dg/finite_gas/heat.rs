//! Original gas-source heating on a solved moving-aperture trajectory.
use super::*;
use crate::bernoulli_aperture::dg::accept_acoustic_candidate_observed;
use crate::bernoulli_aperture::dynamic::{ApertureFrame, DynamicAperture};
use crate::bernoulli_aperture::wall::material::FiniteGasDgInternalHeat;
use crate::dg_port::finite_gas::{FiniteGasDgHeatBudget, check};
use fs_material::gas::GasState;

/// One complete, still-private mechanical/acoustic sample and its physical heat.
/// Thermal, chemical and transport trials reuse it without repeating either
/// the nonlinear source solve or the acoustic source-integral preparation.
#[derive(Clone)]
pub(crate) struct PreparedFiniteGasApertureHeat<A = LoadedAperture> {
    marched: ApertureDg<A>,
    frame: DgApertureFrame,
    initial_energy_j: f64,
    source_heat: SourceHeat,
}

impl<A> core::fmt::Debug for PreparedFiniteGasApertureHeat<A> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedFiniteGasApertureHeat")
            .field("frame", &self.frame)
            .field("source_heat", &self.source_heat)
            .finish_non_exhaustive()
    }
}

impl<A> PreparedFiniteGasApertureHeat<A> {
    pub(crate) fn heat(&self) -> FiniteGasDgInternalHeat {
        self.source_heat.gas_heat()
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare<A: Clone, E, C: FnMut() -> bool>(
    original: &ApertureDg<A>,
    gas: &ClosedGasParcel,
    drive: TubeDrive,
    initial_energy_j: f64,
    controls: FiniteGasDgSolve,
    budgets: FiniteGasDgHeatBudget,
    cancelled: &mut C,
    mechanics: impl FnOnce(
        &mut A,
        ApertureDrive,
    ) -> Result<(ApertureFrame, [MaterialPressureFrame; 2]), AcousticRealizeError>,
) -> Result<PreparedFiniteGasApertureHeat<A>, FiniteGasDgError<E>> {
    admit(original, gas, controls, cancelled)?;
    let mut marched = original.clone();
    let relation = marched.load.port_response().map_err(acoustic)?;
    let (aperture, spatial_pressure) = mechanics(
        &mut marched.aperture,
        ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: relation.history_pressure_pa / 2.0,
            body_flow_m3_s: drive.body_flow_m3_s,
        },
    )
    .map_err(acoustic)?;
    check(cancelled)?;
    let (frame, source_heat) = accept_acoustic_candidate_observed(
        &mut marched.load,
        relation,
        aperture,
        spatial_pressure,
        drive,
        // Keep the finite-gas refusal type intact. On an observer refusal the
        // native trial can only accept into this disposable PRIVATE clone.
        |trial| Ok(trial.finite_gas_heat(controls, budgets, cancelled)),
    )
    .map_err(acoustic)?;
    let source_heat = source_heat?;
    check(cancelled)?;
    Ok(PreparedFiniteGasApertureHeat {
        marched,
        frame,
        initial_energy_j,
        source_heat,
    })
}

#[allow(clippy::too_many_arguments)]
fn complete<A: Clone, S, E, C: FnMut() -> bool>(
    original: &ApertureDg<A>,
    prepared: &PreparedFiniteGasApertureHeat<A>,
    gas: &ClosedGasParcel,
    controls: FiniteGasDgSolve,
    evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
    cancelled: &mut C,
    rebind: impl FnOnce(&mut A, &GasState, f64) -> Result<(), AcousticRealizeError>,
) -> Result<(ApertureDg<A>, ClosedGasParcel, FiniteGasApertureDgStep<S>), FiniteGasDgError<E>> {
    let (load, next_gas, report) = finish(
        &original.load,
        &prepared.marched.load,
        gas,
        prepared.frame.clone(),
        prepared.initial_energy_j,
        controls,
        prepared.source_heat,
        evolve,
        cancelled,
    )?;
    // Only the mechanics need another private copy here: the existing caloric
    // owner has already produced the complete refreshed acoustic load.
    let mut aperture = prepared.marched.aperture.clone();
    rebind(
        &mut aperture,
        load.duct().inlet_gas(),
        load.impedance_pa_s_m3(),
    )
    .map_err(acoustic)?;
    check(cancelled)?;
    Ok((
        ApertureDg {
            aperture,
            load,
            upstream: None,
            supply_channel: None,
        },
        next_gas,
        report,
    ))
}

fn caloric_heat<A, E>(
    prepared: &PreparedFiniteGasApertureHeat<A>,
    external_heat_j: f64,
) -> Result<f64, FiniteGasDgError<E>> {
    let internal_heat_j = prepared.heat().total_j()?;
    let total = if internal_heat_j == 0.0 {
        external_heat_j
    } else {
        external_heat_j + internal_heat_j
    };
    if !external_heat_j.is_finite() || !total.is_finite() {
        return Err(FiniteGasDgError::Invalid(
            "external and original moving-load gas heat must have a finite sum",
        ));
    }
    Ok(total)
}

impl ApertureDg {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn prepare_finite_gas_source_heat<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        drive: TubeDrive,
        force_n: f64,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<PreparedFiniteGasApertureHeat<DynamicAperture>, FiniteGasDgError<E>> {
        prepare(
            self, gas, drive, self.stored_energy_j(), controls, budgets, cancelled,
            |aperture, drive| {
                let trial = aperture.preview_step_with_force(drive, force_n)?;
                let frame = trial.frame;
                aperture.accept_frame(trial);
                Ok((frame, [MaterialPressureFrame::default(); 2]))
            },
        )
    }

    /// Complete only the original caloric/material feedback around this sample.
    /// Include `prepared.heat()` in the actual caloric proposal; the returned
    /// external-heat account excludes that internal gas transfer.
    pub(crate) fn complete_finite_gas_source_heat<S, E, C: FnMut() -> bool>(
        &self,
        prepared: &PreparedFiniteGasApertureHeat<DynamicAperture>,
        gas: &ClosedGasParcel,
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<(Self, ClosedGasParcel, FiniteGasApertureDgStep<S>), FiniteGasDgError<E>> {
        complete(self, prepared, gas, controls, evolve, cancelled, DynamicAperture::rebind_exact_gas_load)
    }

    /// Return physical bore and main-inventory terminal viscosity/molecular heat
    /// to the original finite gas while advancing the moving source once.
    ///
    /// The original nonlinear solve supplies the actual held DG inlet flow;
    /// source observers also use the SAME solved outlet and branch inputs.
    /// The complete main-bore state and existing source budgets/128-DOF limit
    /// apply. Numerical, wall/solid, independent-branch and all aperture losses
    /// remain exported: no gas inventory or deposition partition is inferred
    /// for the moving gap, jet, contact or structural damping.
    ///
    /// External heat stays separate. The original fixed-volume, frozen-medium
    /// sample and first-order background split remain, including time/space
    /// refinement and explicit upstream/supply inventory refusals. Only the
    /// final successful proposal publishes mechanics, gas, load and histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_viscous_molecular_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        force_n: f64,
        external_heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        let prepared = self.prepare_finite_gas_source_heat(gas, drive, force_n, controls, budgets, cancelled)?;
        let total_heat = caloric_heat(&prepared, external_heat_j)?;
        let (next, next_gas, report) = self.complete_finite_gas_source_heat(
            &prepared, gas, controls,
            |old, _, work, cancel| {
                let mut background = heated_finite_gas(old, total_heat, work, gas_solve, cancel)?;
                background.supplied_heat_j = external_heat_j;
                Ok(background)
            },
            cancelled,
        )?;
        *self = next;
        *gas = next_gas;
        Ok(report)
    }
}

impl ApertureDg<LoadedAperture> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn prepare_finite_gas_source_heat<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<PreparedFiniteGasApertureHeat, FiniteGasDgError<E>> {
        prepare(
            self, gas, drive, self.stored_energy_j(), controls, budgets, cancelled,
            |aperture, drive| {
                let trial = aperture.preview(drive, forces_n, pressures_pa)?;
                let frame = trial.frame();
                let pressure = trial.spatial_pressure();
                aperture.accept(trial);
                Ok((frame, pressure))
            },
        )
    }

    /// Reuse original mechanics and acoustic history in each background trial.
    /// The callback includes prepared internal heat in caloric energy while
    /// retaining only external heat in `FiniteGasDgBackground::supplied_heat_j`.
    pub(crate) fn complete_finite_gas_source_heat<S, E, C: FnMut() -> bool>(
        &self,
        prepared: &PreparedFiniteGasApertureHeat,
        gas: &ClosedGasParcel,
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<(Self, ClosedGasParcel, FiniteGasApertureDgStep<S>), FiniteGasDgError<E>> {
        complete(self, prepared, gas, controls, evolve, cancelled, LoadedAperture::rebind_exact_gas_load)
    }

    /// The same original gas-loss heating with single, paired or material-modal
    /// mechanics and their independent force/pressure work ports. Physical heat
    /// ownership and the frozen-medium split are those of the scalar method.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_viscous_molecular_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        external_heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        let prepared = self.prepare_finite_gas_source_heat(
            gas, drive, forces_n, pressures_pa, controls, budgets, cancelled,
        )?;
        let total_heat = caloric_heat(&prepared, external_heat_j)?;
        let (next, next_gas, report) = self.complete_finite_gas_source_heat(
            &prepared, gas, controls,
            |old, _, work, cancel| {
                let mut background = heated_finite_gas(old, total_heat, work, gas_solve, cancel)?;
                background.supplied_heat_j = external_heat_j;
                Ok(background)
            },
            cancelled,
        )?;
        *self = next;
        *gas = next_gas;
        Ok(report)
    }
}
