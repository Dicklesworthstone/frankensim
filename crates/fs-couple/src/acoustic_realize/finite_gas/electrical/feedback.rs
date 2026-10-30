//! Endpoint electrical work on the original material/thermal/flow owners.
//!
//! For every trial deposited work J, restart the SAME accepted trajectory and
//! solve its thermal endpoint. The original conductor then supplies I_end and
//! R_end. A bounded residual closes J = dt I_end^2 R_end before publication.

use super::{
    FiniteGasElectrothermalStringFrame, FiniteGasStringError, FiniteGasStringHeat,
    FiniteGasStringSolve, ThermalMaterialStringRuntime, ThermalStringError,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::acoustic_realize::finite_gas::exchange::{
    FiniteGasStringExchangeError, FiniteGasStringExchangeFrame,
};
use crate::acoustic_realize::finite_gas::flow::FiniteGasStringFlowFrame;
use crate::cavity_phs::piston::thermal::{
    HeatExchangeSolve, feedback,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::conductor::OhmicDrive;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

/// Explicit source-work domain and bounded nonlinear work, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasElectricalFeedbackSolve {
    /// Upper bound of the nonnegative trial Joule-work bracket [J]. Every
    /// attempted endpoint must also satisfy the original material domains.
    pub maximum_joule_work_j: f64,
    /// Absolute tolerance on J - dt I_end^2 R_end [J].
    pub source_tolerance_j: f64,
    /// At most this many safeguarded iterations, plus two bracket evaluations.
    pub maximum_iterations: usize,
}

impl FiniteGasElectricalFeedbackSolve {
    pub fn validate(self) -> Result<(), FiniteGasElectricalFeedbackError> {
        if !self.maximum_joule_work_j.is_finite()
            || self.maximum_joule_work_j <= 0.0
            || !self.source_tolerance_j.is_finite()
            || self.source_tolerance_j <= 0.0
            || self.source_tolerance_j > self.maximum_joule_work_j
            || !(1..=128).contains(&self.maximum_iterations)
        {
            return Err(FiniteGasElectricalFeedbackError::Invalid(
                "electrical feedback requires positive finite work/tolerance and 1..128 iterations",
            ));
        }
        Ok(())
    }
}

/// No failure advances accepted material, gas, source or observer history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGasElectricalFeedbackError {
    Invalid(&'static str),
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
}
impl core::fmt::Display for FiniteGasElectricalFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas electrical feedback: {self:?}")
    }
}
impl std::error::Error for FiniteGasElectricalFeedbackError {}

/// Accepted endpoint electrical/thermal response and actual boundary-work audit.
#[derive(Debug)]
pub struct FiniteGasElectricalFeedbackFrame<R> {
    /// Original conductor and thermal/flow reports. `joule_heat_j` is the
    /// deposited trial J, not an additional heat credit beside source work.
    pub electrical: FiniteGasElectrothermalStringFrame<R>,
    /// Actual backward-Euler boundary work dt I_end^2 R_end [J].
    pub source_work_j: f64,
    /// Deposited work minus actual boundary work [J]. No state is renormalized
    /// to hide this admitted nonlinear residual.
    pub source_residual_j: f64,
    /// Number of complete private electrical/thermal proposals.
    pub trials: usize,
    /// Actual boundary work plus other inputs minus all storage/losses [J].
    /// The nested thermal account instead uses deposited heat J.
    pub energy_balance_residual_j: f64,
    /// Original coupled allowance, source residual allowance and roundoff [J].
    pub energy_tolerance_j: f64,
}

type ConductionFrame = FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>;
type ConductionError = FiniteGasStringExchangeError<LayerConductionError>;

struct Candidate {
    runtime: ThermalMaterialStringRuntime,
    gas: IsochoricGasParcel,
    frame: FiniteGasElectricalFeedbackFrame<ConductionFrame>,
}

fn source_error(error: FiniteGasElectricalFeedbackError) -> ConductionError {
    FiniteGasStringError::ElectricalFeedback(error).into()
}

fn solve_error(error: feedback::ScalarError<ConductionError>) -> ConductionError {
    use feedback::ScalarError;
    let error = match error {
        ScalarError::Evaluation(error) => return error,
        ScalarError::Cancelled => {
            return FiniteGasStringError::String(ThermalStringError::Acoustic(
                AcousticRealizeError::Cancelled,
            ))
            .into();
        }
        ScalarError::Invalid(message) => FiniteGasElectricalFeedbackError::Invalid(message),
        ScalarError::OutsideBracket { lower, upper } => {
            FiniteGasElectricalFeedbackError::OutsideBracket {
                lower_residual_j: lower,
                upper_residual_j: upper,
            }
        }
        ScalarError::NotConverged {
            iterations,
            residual,
        } => FiniteGasElectricalFeedbackError::NotConverged {
            iterations,
            residual_j: residual,
        },
    };
    source_error(error)
}

impl ThermalMaterialStringRuntime {
    /// Close Joule work and reciprocal material/gas conduction at one endpoint.
    ///
    /// All thermal and thermoelastic feedback uses the original owners. At
    /// each proposed J, the conductor is resolved from that SAME final solid
    /// state: held current gives I^2 R(T_end) dt and held voltage gives
    /// V^2/R(T_end) dt. No trial consumes heat, inventory or an audio sample.
    /// Source domains, temperature limits and inner thermal work caps apply
    /// at every trial, including the caller-authored bracket endpoint.
    ///
    /// The source residual and actual boundary-work balance are separately
    /// admitted before atomic publication. Mechanics remains frozen-initial
    /// over one sample; this is backward-Euler electrical/thermal feedback,
    /// not a circuit-storage, induction or second-order mechanical model.
    /// A bracketed solution is not a uniqueness or global-stability proof.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_implicit_drive_and_isochoric_gas_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<
        FiniteGasElectricalFeedbackFrame<FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>>,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        let Candidate {
            runtime,
            gas: next_gas,
            frame,
        } = self.propose_implicit_conduction(
            cx,
            expected_epoch,
            forces,
            drive,
            gas,
            heat,
            controls,
            exchange,
            layers,
            electrical,
        )?;
        super::super::check(cx)?;
        let _old_string = core::mem::replace(self, runtime);
        let _old_gas = core::mem::replace(gas, next_gas);
        Ok(frame)
    }

    /// The same converged electrical/thermal endpoint followed by native
    /// conservative nozzle transport over the same sample interval. Flow runs
    /// ONCE after convergence, not once per trial or as a second audio sample.
    /// It cannot alter the solid endpoint used by the accepted conductor law.
    ///
    /// The common flow owner audits actual electrical source work, including
    /// its explicitly bounded difference from deposited J, together with both
    /// extensive gas inventories. Final source/flow/observer refusal preserves
    /// all three accepted owners. Heat-then-flow remains first-order splitting.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<
        FiniteGasElectricalFeedbackFrame<
            FiniteGasStringFlowFrame<FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>>,
        >,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        let Candidate {
            runtime,
            gas: heated_gas,
            frame,
        } = self.propose_implicit_conduction(
            cx,
            expected_epoch,
            forces,
            drive,
            gas,
            heat,
            controls,
            exchange,
            layers,
            electrical,
        )?;
        let combined_heat = FiniteGasStringHeat {
            // This exact sum was already admitted by the shared source owner.
            string_heat_j: heat.string_heat_j + frame.electrical.joule_heat_j,
            gas_heat_j: heat.gas_heat_j,
        };
        let source_work_j = frame.source_work_j;
        let source_residual_j = frame.source_residual_j;
        let trials = frame.trials;
        let electrical_frame = frame.electrical.try_map_coupled(|thermal| {
            self.step_with_finite_gas_nozzle(
                cx,
                expected_epoch,
                forces,
                gas,
                reservoir,
                combined_heat,
                controls,
                nozzle,
                flow_controls,
                -source_residual_j,
                electrical.source_tolerance_j,
                |candidate, gas, _| {
                    *candidate = runtime;
                    *gas = heated_gas;
                    Ok::<_, ConductionError>(thermal)
                },
                |thermal| &mut thermal.coupled,
                |_, _| Ok(()),
            )
        })?;
        // The shared flow owner has already admitted and atomically published
        // the entire candidate. No fallible work follows that boundary.
        Ok(FiniteGasElectricalFeedbackFrame {
            source_work_j,
            source_residual_j,
            trials,
            energy_balance_residual_j: electrical_frame.coupled.energy_balance_residual_j,
            energy_tolerance_j: electrical_frame.coupled.energy_tolerance_j,
            electrical: electrical_frame,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn propose_implicit_conduction(
        &self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<Candidate, ConductionError> {
        electrical.validate().map_err(source_error)?;
        let (conductor_before, current_a, _) = self.prepare_drive(cx, expected_epoch, drive)?;
        let dt = self.mechanical.member.model.sample_period_s();
        let candidate = feedback::solve_scalar(
            [0.0, electrical.maximum_joule_work_j],
            electrical.source_tolerance_j,
            electrical.maximum_iterations,
            |joule_heat_j, trials, _| -> Result<_, ConductionError> {
                let (runtime, next_gas, report) = self.propose_joule_work_and_finite_gas(
                    drive,
                    gas,
                    heat,
                    conductor_before.clone(),
                    current_a,
                    joule_heat_j,
                    |candidate, gas, heat| {
                        candidate.step_with_isochoric_gas_conduction(
                            cx,
                            expected_epoch,
                            forces,
                            gas,
                            heat,
                            controls,
                            exchange,
                            layers,
                        )
                    },
                )?;
                let power = report
                    .conductor_after
                    .joule_power_w(report.final_current_a)
                    .map_err(ThermalStringError::Electrical)?;
                let source_work_j = dt * power;
                let source_residual_j = joule_heat_j - source_work_j;
                let coupled = &report.coupled.coupled;
                let energy_balance_residual_j =
                    coupled.energy_balance_residual_j - source_residual_j;
                let energy_tolerance_j = coupled.energy_tolerance_j
                    + electrical.source_tolerance_j
                    + 32.0 * f64::EPSILON * (joule_heat_j.abs() + source_work_j.abs());
                if !source_work_j.is_finite()
                    || (power > 0.0 && source_work_j == 0.0)
                    || !source_residual_j.is_finite()
                    || !energy_balance_residual_j.is_finite()
                    || !energy_tolerance_j.is_finite()
                {
                    return Err(source_error(FiniteGasElectricalFeedbackError::Invalid(
                        "endpoint source work or its complete energy allowance is unrepresentable",
                    )));
                }
                Ok((
                    Candidate {
                        runtime,
                        gas: next_gas,
                        frame: FiniteGasElectricalFeedbackFrame {
                            electrical: report,
                            source_work_j,
                            source_residual_j,
                            trials,
                            energy_balance_residual_j,
                            energy_tolerance_j,
                        },
                    },
                    source_residual_j,
                ))
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(solve_error)?;
        if candidate.frame.energy_balance_residual_j.abs() > candidate.frame.energy_tolerance_j {
            return Err(source_error(FiniteGasElectricalFeedbackError::Invalid(
                "actual endpoint electrical work does not close the combined storage balance",
            )));
        }
        Ok(candidate)
    }
}
