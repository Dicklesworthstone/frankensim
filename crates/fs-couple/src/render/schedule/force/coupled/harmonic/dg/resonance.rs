//! Actual driven peaks of the complete native modal/DG/terminal system. Every
//! query retains the same branches, source memories, finite gaps and closures.

pub use super::super::resonance::{
    DrivenResonance, HalfPowerBandwidth, ModalHarmonicObservation, ResonanceControls,
    ResonanceFrequency, ResonanceSearchError,
};
use super::*;
use fs_duct::dg::DgHarmonicDerivative;
use fs_phs::harmonic::resonance::search_resonance;

/// One fixed observation in the original source coordinates and physical units.
/// Port IDs refer to the declared source, never to its current vector position.
/// A main-bore side port and its finite branch's exterior are distinct locations.
#[derive(Clone, Copy, Debug)]
pub enum DgHarmonicObservation<'a> {
    /// Signed physical mechanical displacement, velocity or acceleration projection.
    Modal(ModalHarmonicObservation<'a>),
    /// Main-bore port pressure [Pa], including the driving inlet if selected.
    BorePressure(PortId),
    /// Flow INTO a main-bore port [m3/s]. At a branch junction its sign opposes
    /// flow into the branch; compressible branch storage is retained.
    BoreInflow(PortId),
    /// Pressure [Pa] at a local port of the named finite branch. `branch` is
    /// its main-bore attachment ID; `port` belongs to the branch duct itself.
    BranchPressure { branch: PortId, port: PortId },
    /// Flow INTO that local branch port [m3/s], with the same source ID rules.
    BranchInflow { branch: PortId, port: PortId },
    /// Exterior terminal pressure [Pa], after any finite passage pressure drop.
    TerminalPressure(PortId),
    /// Flow INTO the original rational terminal [m3/s].
    TerminalOutwardFlow(PortId),
    /// Duct-to-exterior pressure drop, or rigid closure reaction [Pa].
    ChannelPressureDrop(PortId),
    /// Flow of the original resolved finite channel field [m3/s]. Rigid
    /// closures have no oscillatory flow and cannot define a flow peak.
    ChannelFlow(PortId),
}

/// Search refusals retain the complete native modal/acoustic source error.
pub type DgDrivenResonanceError = ResonanceSearchError<DgHarmonicCouplingError>;

impl DgHarmonicObservation<'_> {
    fn validate(self, modes: usize) -> Result<(), DgHarmonicCouplingError> {
        if let Self::Modal(observation) = self {
            observation.validate(modes)?;
        }
        Ok(())
    }

    fn admit_loaded_flow(self, load: &DgPortLoad) -> Result<(), DgHarmonicCouplingError> {
        let closed = |port| {
            load.terminal_channels().iter().any(|channel| {
                channel.input().port == port
                    && matches!(channel.input().kind, DgTerminalChannelKind::Closed)
            })
        };
        let is_closed = match self {
            Self::TerminalOutwardFlow(port) | Self::ChannelFlow(port) => closed(port),
            Self::BoreInflow(port) => {
                closed(port) && !load.branches().iter().any(|branch| branch.port() == port)
            }
            Self::BranchInflow { branch, port } => {
                closed(branch)
                    && load.branches().iter().any(|source| {
                        source.port() == branch && source.duct().outlet() == DgBoundary::Port(port)
                    })
            }
            _ => false,
        };
        if is_closed {
            // The exact constraint is zero even when a dense solve retains tiny
            // roundoff. Never reinterpret that noise as a physical flow peak.
            return Err(invalid("a rigid closure has no oscillatory flow resonance").into());
        }
        Ok(())
    }

    fn observe(
        self,
        response: &DgLoadedModalHarmonicResponse,
        tangent: Option<&DgLoadedModalHarmonicFrequencyResponse>,
        gate: &CancelGate,
    ) -> Result<(C64, C64), DgHarmonicCouplingError> {
        poll(Some(gate))?;
        match self {
            Self::Modal(observation) => Ok(observation.observe(
                &response.modal,
                tangent.map(|value| &value.modal),
                gate,
            )?),
            Self::BorePressure(port) | Self::BoreInflow(port) => duct_observation(
                port,
                matches!(self, Self::BorePressure(_)),
                &response.duct,
                tangent.map(|value| &value.duct),
            ),
            Self::BranchPressure { branch, port } | Self::BranchInflow { branch, port } => {
                let source = response
                    .branches
                    .iter()
                    .find(|value| value.port == branch)
                    .ok_or_else(|| {
                        invalid("resonance observation requires the original finite branch")
                    })?;
                let derivative = tangent
                    .map(|value| {
                        value
                            .branches
                            .iter()
                            .find(|value| value.port == branch)
                            .map(|value| &value.duct)
                            .ok_or_else(|| {
                                invalid("resonance tangent lost the original finite branch")
                            })
                    })
                    .transpose()?;
                duct_observation(
                    port,
                    matches!(self, Self::BranchPressure { .. }),
                    &source.duct,
                    derivative,
                )
            }
            Self::TerminalPressure(port) | Self::TerminalOutwardFlow(port) => {
                let source = response.terminal.as_ref().ok_or_else(|| {
                    invalid("resonance observation requires the original rational terminal")
                })?;
                let index = source
                    .ports
                    .iter()
                    .position(|id| *id == port)
                    .ok_or_else(|| {
                        invalid("resonance observation names an absent terminal port")
                    })?;
                let pressure = matches!(self, Self::TerminalPressure(_));
                let value = if pressure {
                    source.pressure_pa[index]
                } else {
                    source.outward_flow_m3_s[index]
                };
                let derivative =
                    tangent
                        .map(|value| {
                            let source = value.terminal.as_ref().ok_or_else(|| {
                                invalid("resonance tangent lost the original terminal")
                            })?;
                            let index = source.ports.iter().position(|id| *id == port).ok_or_else(
                                || invalid("resonance tangent lost the original terminal port"),
                            )?;
                            Ok::<_, ModalCouplingError>(if pressure {
                                source.pressure_pa[index]
                            } else {
                                source.outward_flow_m3_s[index]
                            })
                        })
                        .transpose()?
                        .unwrap_or(C64::ZERO);
                Ok((value, derivative))
            }
            Self::ChannelPressureDrop(port) | Self::ChannelFlow(port) => {
                let source = response
                    .channels
                    .iter()
                    .find(|value| value.port == port)
                    .ok_or_else(|| {
                        invalid(
                            "resonance observation requires the original finite passage or closure",
                        )
                    })?;
                let pressure = matches!(self, Self::ChannelPressureDrop(_));
                let value = if pressure {
                    source.pressure_drop_pa
                } else {
                    source
                        .channel
                        .as_ref()
                        .ok_or_else(|| {
                            invalid("a rigid closure has no oscillatory flow resonance")
                        })?
                        .flow_m3_s
                };
                let derivative = tangent
                    .map(|value| {
                        let source = value
                            .channels
                            .iter()
                            .find(|value| value.port == port)
                            .ok_or_else(|| {
                                invalid("resonance tangent lost the original passage or closure")
                            })?;
                        Ok::<_, ModalCouplingError>(if pressure {
                            source.pressure_drop_pa
                        } else {
                            source
                                .channel
                                .as_ref()
                                .ok_or_else(|| {
                                    invalid("resonance tangent lost the original finite passage")
                                })?
                                .flow_m3_s
                        })
                    })
                    .transpose()?
                    .unwrap_or(C64::ZERO);
                Ok((value, derivative))
            }
        }
    }
}

fn duct_observation(
    port: PortId,
    pressure: bool,
    response: &DgHarmonicResponse,
    tangent: Option<&DgHarmonicDerivative>,
) -> Result<(C64, C64), DgHarmonicCouplingError> {
    let source = response
        .ports
        .iter()
        .find(|value| value.id == port)
        .ok_or_else(|| invalid("resonance observation names an absent original duct port"))?;
    let value = if pressure {
        source.pressure_pa
    } else {
        source.inflow_m3_s
    };
    let derivative = tangent
        .map(|value| {
            let source = value
                .ports
                .iter()
                .find(|value| value.id == port)
                .ok_or_else(|| invalid("resonance tangent lost the original duct port"))?;
            Ok::<_, ModalCouplingError>(if pressure {
                source.pressure_pa
            } else {
                source.inflow_m3_s
            })
        })
        .transpose()?
        .unwrap_or(C64::ZERO);
    Ok((value, derivative))
}

fn sample(
    observation: DgHarmonicObservation<'_>,
    response: &DgLoadedModalHarmonicResponse,
    tangent: Option<&DgLoadedModalHarmonicFrequencyResponse>,
    gate: &CancelGate,
) -> Result<(C64, C64, f64), DgHarmonicCouplingError> {
    let (value, derivative) = observation.observe(response, tangent, gate)?;
    let mut residual = response
        .solve_relative_residual
        .max(response.modal.solve_relative_residual)
        .max(response.duct.solve_relative_residual);
    for branch in &response.branches {
        residual = residual.max(branch.duct.solve_relative_residual);
    }
    for channel in &response.channels {
        if let Some(field) = &channel.channel {
            residual = residual.max(field.solve_relative_residual);
        }
    }
    if let Some(tangent) = tangent {
        residual = residual
            .max(tangent.solve_relative_residual)
            .max(tangent.modal.relative_derivative_residual)
            .max(tangent.duct.solve_relative_residual);
        for branch in &tangent.branches {
            residual = residual.max(branch.duct.solve_relative_residual);
        }
        for channel in &tangent.channels {
            if let Some(field) = &channel.channel {
                residual = residual.max(field.solve_relative_residual);
            }
        }
    }
    Ok((value, derivative, residual))
}

impl CoupledModalSystem {
    /// Resolve a driven peak and optional half-power width with the complete
    /// original bare DG duct attached. The selected response is under FIXED
    /// authored force, not an inferred pressure/flow impedance or free pole.
    /// Mechanical and duct observations can therefore have different peaks.
    ///
    /// The bracket, width and evaluation controls have the same meaning as
    /// [`Self::driven_resonance`]. Native source bands, complete-state work and
    /// residual admission remain active for every query; `maximum_work_per_frequency`
    /// bounds each solve and `controls.maximum_evaluations` bounds query count.
    /// Derivatives reuse the native joint LU and crossings need only a primary
    /// solve. No isolated duct impedance is inverted and no accepted state changes.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_resonance_with_duct(
        &self,
        force: &[C64],
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        observation: DgHarmonicObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        maximum_work_per_frequency: usize,
        gate: &CancelGate,
    ) -> Result<DrivenResonance, DgDrivenResonanceError> {
        poll(Some(gate)).map_err(|error| ResonanceSearchError::Source(error.into()))?;
        observation
            .validate(self.mode_count())
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(
            interval_rad_s,
            controls,
            || gate.is_requested(),
            |omega, derivative, _| {
                if derivative {
                    let result = self.harmonic_response_with_duct_frequency_derivative(
                        omega,
                        force,
                        None,
                        duct,
                        inlet,
                        flow_column,
                        maximum_work_per_frequency,
                        gate,
                    )?;
                    sample(observation, &result.response, Some(&result), gate)
                } else {
                    let result = self.harmonic_response_with_duct(
                        omega,
                        force,
                        duct,
                        inlet,
                        flow_column,
                        maximum_work_per_frequency,
                        gate,
                    )?;
                    sample(observation, &result, None, gate)
                }
            },
        )
    }

    /// Resolve an observed peak and linewidth of the COMPLETE source-bound
    /// load: modal and Maxwell states, reciprocal air, original DG gas and wall
    /// memory, every finite branch, rational MIMO state, gap cell and closure.
    /// Every frequency uses the existing joint physical response; no attachment,
    /// loss dispersion, mutual reaction or algebraic constraint is omitted.
    ///
    /// Excitation and observation stay fixed. Local pressure, modal motion and
    /// terminal flow peaks need not coincide. The optional width describes the
    /// selected amplitude squared, not material heat or a pole's decay quality.
    /// The caller must isolate the peak and both requested half-height crossings
    /// inside the original interval and admitted source bands. Numerical brackets
    /// do not certify continuum error or peak uniqueness.
    ///
    /// A known rigid closure's flow explicitly refuses: solve roundoff is not a
    /// physical flow resonance. Its pressure reaction remains observable, as does
    /// flow into a finite compressible branch attached to a closed exterior.
    /// Original Mach, source-band, work, residual and cancellation refusals are
    /// preserved. Driven mouths retain their native explicit scope refusal.
    /// No accepted state, source history or sample clock changes.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_resonance_with_load(
        &self,
        force: &[C64],
        load: &DgPortLoad,
        flow_column: &[f64],
        observation: DgHarmonicObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        maximum_work_per_frequency: usize,
        gate: &CancelGate,
    ) -> Result<DrivenResonance, DgDrivenResonanceError> {
        poll(Some(gate)).map_err(|error| ResonanceSearchError::Source(error.into()))?;
        observation
            .validate(self.mode_count())
            .map_err(ResonanceSearchError::Source)?;
        observation
            .admit_loaded_flow(load)
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(
            interval_rad_s,
            controls,
            || gate.is_requested(),
            |omega, derivative, _| {
                if derivative {
                    let result = self.harmonic_response_with_load_frequency_derivative(
                        omega,
                        force,
                        None,
                        load,
                        flow_column,
                        maximum_work_per_frequency,
                        gate,
                    )?;
                    sample(observation, &result.response, Some(&result), gate)
                } else {
                    let result = self.harmonic_response_with_load(
                        omega,
                        force,
                        load,
                        flow_column,
                        maximum_work_per_frequency,
                        gate,
                    )?;
                    sample(observation, &result, None, gate)
                }
            },
        )
    }
}
