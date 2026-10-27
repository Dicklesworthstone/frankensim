//! Terminal observation over the unchanged physical reed and pressure scheduler.
use super::*;
use crate::acoustic_realize::map_drive;
use crate::driving_point::{DrivingPointError, regional::RegionalDuct};
use crate::driving_point::regional::outlet::OutletPressure;
use crate::pcm_wav::baffled::CircularOutletReceiver;

pub(super) struct TerminalObservation {
    pressure: OutletPressure,
    holes: Vec<OutletPressure>,
    failure: Option<OutletError>,
}

pub(super) enum ObservationInput {
    Legacy(f64),
    Baffled { receiver: CircularOutletReceiver, exterior: Option<GasState> },
}
impl ObservationInput {
    pub(super) fn validate_policy(
        &self,
        policy: ObservationPolicy,
    ) -> Result<(), ReedPerformanceError> {
        if policy == ObservationPolicy::Acceptance && matches!(self, Self::Legacy(_)) {
            return Err(ReedPerformanceError::Observation(
                ObservationError::LegacyListenerInAcceptance,
            ));
        }
        Ok(())
    }
    pub(super) fn read(reader: &mut Reader<'_>) -> Result<Self, ReedPerformanceError> {
        Self::read_with_gas_sources(reader, None, &mut Vec::new())
    }

    pub(super) fn read_with_gas_sources(
        reader: &mut Reader<'_>,
        sources: Option<&ThermochemicalSources<'_>>,
        records: &mut Vec<GasRecord>,
    ) -> Result<Self, ReedPerformanceError> {
        let mut row = reader.next()?;
        let mut observation = match row.word()? {
            "listener_m" => Self::Legacy(row.positive()?),
            "baffled_outlet" => Self::Baffled { receiver: CircularOutletReceiver {
                position_m: [row.scalar()?, row.scalar()?, row.scalar()?],
                radial_rings: row.count(128)?, angular_points: row.count(512)?,
                maximum_frequency_hz: row.positive()?,
            }, exterior: None },
            _ => return Err(bad(row.line, "expected listener_m or baffled_outlet observation")),
        };
        row.finish()?;
        if reader.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("exterior") {
            let mut row = reader.row("exterior")?;
            let Self::Baffled { exterior, .. } = &mut observation else {
                return Err(bad(row.line, "an exterior gas requires the physical baffled_outlet observation"));
            };
            *exterior = Some(row.gas_with_sources(None, sources, records)?);
            row.finish()?;
        }
        Ok(observation)
    }

    /// Distance the observer law covers: the proxy's listener distance, or the
    /// baffled receiver's distance from the outlet centre.
    pub(super) fn range_m(&self) -> f64 {
        match self {
            Self::Legacy(distance) => *distance,
            Self::Baffled { receiver, .. } => {
                let [x, y, z] = receiver.position_m;
                x.hypot(y).hypot(z)
            }
        }
    }

    /// The receiver's declared band [Hz]; the legacy proxy declares none.
    pub(super) const fn declared_band_hz(&self) -> Option<f64> {
        match self {
            Self::Baffled { receiver, .. } => Some(receiver.maximum_frequency_hz),
            Self::Legacy(_) => None,
        }
    }

    /// The baffled receiver's position and its authored exterior gas.
    pub(super) fn baffled(&self) -> Option<([f64; 3], Option<GasState>)> {
        match self {
            Self::Legacy(_) => None,
            Self::Baffled { receiver, exterior } => Some((receiver.position_m, *exterior)),
        }
    }

    pub(super) fn validate_termination(&self, termination: Termination, line: usize)
        -> Result<(), ReedPerformanceError> {
        if matches!(self, Self::Baffled { .. }) && !matches!(termination, Termination::FlangedOpen) {
            return Err(bad(line, "baffled_outlet requires the same flanged terminal feedback load"));
        }
        Ok(())
    }

    /// `delay_speed_m_s` references the baffled receiver's propagation delays
    /// (the observer law's receipt names it); the proxy has no such delay.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn build(self, duct: &Duct, gases: &[GasState], reed: BeatingReed,
        termination: Termination, rate: u32, samples: usize, delay_speed_m_s: f64)
        -> Result<(ReedBoreVoice, Option<TerminalObservation>), ReedPerformanceError> {
        let regional = RegionalDuct::new(duct, gases)
            .map_err(|error| ReedPerformanceError::Model(map_drive(DrivingPointError::Duct(error))))?;
        self.build_regional(duct, regional, reed, termination, rate, samples, delay_speed_m_s)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn build_regional(self, duct: &Duct, regional: RegionalDuct, reed: BeatingReed,
        termination: Termination, rate: u32, samples: usize, delay_speed_m_s: f64)
        -> Result<(ReedBoreVoice, Option<TerminalObservation>), ReedPerformanceError> {
        let regional = regional.with_closing_pressure_pa(reed.closing_pressure_pa);
        if let Self::Legacy(distance) = self {
            let voice = ReedBoreVoice::with_line(duct, regional.inlet_gas(), reed,
                PlateBank::default(), distance, rate, |_| {
                    regional.characteristic_line_realized(termination, rate, samples, None, None)
                        .map_err(map_drive)
                }).map_err(ReedPerformanceError::Model)?;
            return Ok((voice, None));
        }
        let Self::Baffled { receiver, exterior } = self else { unreachable!("legacy returned above") };
        let regional = match exterior {
            Some(gas) => regional.with_exterior_gas(gas)
                .map_err(|error| ReedPerformanceError::Model(map_drive(DrivingPointError::Duct(error))))?,
            None => regional,
        };
        // Both feedback and observation see this ONE bound exterior. It never
        // changes internal gas, geometry or the source's Bernoulli density.
        // Use the same initializer, lay, source impedance and prime as every
        // other ReedBoreVoice. Its proxy is discarded only at the observation
        // boundary; no geometry, mechanics or pressure controls are rewritten.
        let [x, y, z] = receiver.position_m;
        let mut voice = ReedBoreVoice::with_line(duct, regional.inlet_gas(), reed,
            PlateBank::default(), x.hypot(y).hypot(z), rate, |_| {
                regional.characteristic_line_realized(termination, rate, samples, None, None)
                    .map_err(map_drive)
            }).map_err(ReedPerformanceError::Model)?;
        let grid = voice.line.history().len();
        let (mut pressure, mut holes, shunt) = regional
            .baffled_outlet_and_open_holes_and_shunt_delayed_at(
                rate,
                grid,
                receiver,
                fs_duct::WallModel::from_pin(None),
                delay_speed_m_s,
            )
            .map_err(ReedPerformanceError::Outlet)?;
        // Sample -1 exists in the original feedback history. Supply exactly
        // that same wave once, before sample zero, to the terminal transfer.
        pressure
            .step(voice.p_plus_prev)
            .map_err(ReedPerformanceError::Outlet)?;
        for hole in &mut holes {
            hole.step(voice.p_plus_prev)
                .map_err(ReedPerformanceError::Outlet)?;
        }
        if let Some(mut shunt) = shunt {
            shunt
                .prime(voice.p_plus_prev)
                .map_err(|error| ReedPerformanceError::Outlet(OutletError::Transfer(error)))?;
            voice.set_register_shunt(shunt);
        }
        Ok((
            voice,
            Some(TerminalObservation {
                pressure,
                holes,
                failure: None,
            }),
        ))
    }
}

#[cfg(test)]
#[path = "outlet/qualification.rs"]
mod qualification;

impl ReedPerformance {
    /// The selected physical receiver; None identifies the original bore/jet proxy.
    /// Axial and exterior travel delays remain on the mechanics clock, not the
    /// output-filter latency removed by ensemble alignment.
    #[must_use]
    pub fn outlet_receiver(&self) -> Option<OutletReceiverInfo> {
        self.terminal.as_ref().map(|terminal| terminal.pressure.info())
    }

    /// Original observer failure, retained even through the legacy PressureRenderer
    /// error envelope. Once present, all further stepping refuses as poisoned.
    #[must_use]
    pub fn last_outlet_error(&self) -> Option<&OutletError> {
        self.terminal.as_ref().and_then(|terminal| terminal.failure.as_ref())
    }

    /// Emit the selected pressure while retaining typed terminal-owner failures.
    ///
    /// The terminal path advances the existing scheduler one sample at a time
    /// to observe every accepted outgoing wave. The legacy path retains its
    /// original block call. This is an observation adapter, not a new mechanical
    /// or control stepper. No real-time performance claim is made for this path.
    ///
    /// # Errors
    /// Request shape/horizon refuses before output or state changes. Physics
    /// failures can leave a partial callback and poison the complete source;
    /// discard that callback rather than resuming either history separately.
    pub fn pressure_block(&mut self, output: &mut [f64]) -> Result<(), ReedPerformanceError> {
        self.validate_sample_count(output.len() as u64).map_err(ReedPerformanceError::Render)?;
        self.renderer.context().validate_block_len(output.len()).map_err(ReedPerformanceError::Render)?;
        // The observer law acts on emitted pressure only, after the receiver:
        // the mechanics and feedback never see it.
        let Some(terminal) = &mut self.terminal else {
            // The proxy's bore pressure travels no distance: the law acts on
            // its radiated part only. A lossless law has nothing to act on and
            // passes the voice's own sum, bit for bit.
            let lossless = self.observer.taps() == [1.0];
            for slot in output.iter_mut() {
                let mut proxy = [0.0];
                self.renderer.block(&mut proxy).map_err(ReedPerformanceError::Render)?;
                let RenderVoice::ReedBore(voice) = &self.renderer.context().voices[0] else {
                    unreachable!("reed performance owns exactly one reed voice")
                };
                let parts = voice.last_observer_parts();
                *slot = if lossless {
                    self.observer.process(proxy[0])
                } else {
                    parts.internal_pa + self.observer.process(parts.radiated_pa)
                };
            }
            return Ok(());
        };
        for slot in output {
            let mut proxy = [0.0];
            self.renderer.block(&mut proxy).map_err(ReedPerformanceError::Render)?;
            let RenderVoice::ReedBore(voice) = &self.renderer.context().voices[0] else {
                unreachable!("reed performance owns exactly one reed voice")
            };
            let p_plus = voice.p_plus_prev;
            let mut sum = match terminal.pressure.step(p_plus) {
                Ok(pressure) => pressure,
                Err(error) => {
                    terminal.failure = Some(error.clone());
                    return Err(ReedPerformanceError::Outlet(error));
                }
            };
            for hole in &mut terminal.holes {
                match hole.step(p_plus) {
                    Ok(hole_pressure) => sum += hole_pressure,
                    Err(error) => {
                        terminal.failure = Some(error.clone());
                        return Err(ReedPerformanceError::Outlet(error));
                    }
                }
            }
            *slot = self.observer.process(sum);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
