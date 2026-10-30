//! Terminal observation over the unchanged physical reed and pressure scheduler.
use super::*;
use crate::acoustic_realize::map_drive;
use crate::driving_point::{DrivingPointError, regional::RegionalDuct};
use crate::driving_point::molecular::MolecularDrivingPoint;
use crate::driving_point::regional::outlet::OutletPressure;
use crate::pcm_wav::baffled::CircularOutletReceiver;

pub(super) struct TerminalObservation {
    pressure: OutletPressure,
    holes: Vec<OutletPressure>,
    molecular_hole_observers: Option<Vec<MolecularHoleObserver>>,
    failure: Option<OutletError>,
}

struct MolecularHoleObserver {
    filter: ObserverFilter,
    omission: Result<ObserverOmission, AcousticRealizeError>,
}

impl TerminalObservation {
    pub(super) fn admit_molecular_hole_observers(
        &self,
        declared: bool,
    ) -> Result<(), AcousticRealizeError> {
        if let Some(observers) = &self.molecular_hole_observers {
            for observer in observers {
                let omission = observer.omission.as_ref().map_err(Clone::clone)?;
                if matches!(observer.filter.receipt().absorption.law, GasAbsorptionModel::StokesKirchhoff) {
                    omission.admit(true, declared)?;
                }
            }
        }
        Ok(())
    }
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

    /// Original molecular constitutive owner in both pressure feedback and
    /// physical terminal flow. There is no equilibrium reconstruction at the
    /// observation boundary and no independent reed state or material law.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn build_molecular(self, point: MolecularDrivingPoint, reed: BeatingReed,
        termination: Termination, rate: u32, samples: usize, delay_speed_m_s: f64,
        declared_law: Option<DeclaredObserverLaw>, records: &[GasRecord])
        -> Result<(ReedBoreVoice, Option<TerminalObservation>), ReedPerformanceError> {
        let delay = 2.0 * point.one_way_wavefront_s() * f64::from(rate);
        if !(delay >= 2.0 && delay < samples as f64 - 2.0) {
            return Err(ReedPerformanceError::Model(map_drive(DrivingPointError::Invalid {
                what: "round-trip delay does not fit the realized history",
            })));
        }
        match &self {
            Self::Legacy(_) if point.exterior_gas() != point.inlet_gas() => {
                return Err(ReedPerformanceError::Model(map_drive(DrivingPointError::Invalid {
                    what: "the compact reed proxy requires the same inlet and exterior gas",
                })));
            }
            Self::Baffled { exterior, .. } => {
                if termination != Termination::FlangedOpen
                    || exterior.is_some_and(|gas| gas != *point.exterior_gas()) {
                    return Err(ReedPerformanceError::Model(map_drive(DrivingPointError::Invalid {
                        what: "molecular outlet observation and feedback require the same flanged exterior load",
                    })));
                }
            }
            Self::Legacy(_) => {}
        }
        // Each mouth has its own distance and hence its own absorption,
        // dispersion and finite-FIR delay allocation. Filtering a sum at the
        // bell's range would silently assign that range to every side opening.
        let mut hole_observers = Vec::new();
        let mut hole_delay_speeds = Vec::new();
        if let Self::Baffled { receiver, .. } = &self {
            let length: f64 = point.geometry().segments.iter().map(Segment::axial_length).sum();
            for site in point.open_hole_sites() {
                let mut position = receiver.position_m;
                position[0] += length - site.axial_position_m;
                let range_m = position[0].hypot(position[1]).hypot(position[2]);
                let gas = point.exterior_gas();
                let filter = observer_filter(declared_law, gas, records, range_m, rate)?;
                admit_receiver_delay(&filter, Some((position, Some(*gas))),
                    site.hole_radius_m, gas.sound_speed)?;
                let omission = observer_omission(gas, records, range_m,
                    observer_band_top_hz(Some(receiver.maximum_frequency_hz), rate));
                hole_delay_speeds.push(filter.receipt().delay_speed_m_s.unwrap_or(gas.sound_speed));
                hole_observers.push(MolecularHoleObserver { filter, omission });
            }
        }
        let voice = ReedBoreVoice::new_with_molecular_driving_point(point, termination, reed,
            PlateBank::default(), self.range_m(), rate, &crate::driving_point::LineLadder::default())
            .map_err(ReedPerformanceError::Model)?;
        let Self::Baffled { receiver, .. } = self else {
            return Ok((voice, None));
        };
        let (mut pressure, mut holes) = voice.molecular_driving_point()
            .expect("molecular constructor retains its original source")
            .baffled_outlet_and_open_holes_with_delay_speeds(rate, voice.line.history().len(),
                receiver, delay_speed_m_s, &hole_delay_speeds)
            .map_err(ReedPerformanceError::Outlet)?;
        // Every opening receives the identical sample -1 already in feedback.
        pressure.step(voice.p_plus_prev).map_err(ReedPerformanceError::Outlet)?;
        for hole in &mut holes {
            hole.step(voice.p_plus_prev).map_err(ReedPerformanceError::Outlet)?;
        }
        Ok((voice, Some(TerminalObservation {
            pressure, holes, molecular_hole_observers: Some(hole_observers), failure: None,
        })))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn build_regional(self, duct: &Duct, regional: RegionalDuct, reed: BeatingReed,
        termination: Termination, rate: u32, samples: usize, delay_speed_m_s: f64)
        -> Result<(ReedBoreVoice, Option<TerminalObservation>), ReedPerformanceError> {
        if let Self::Legacy(distance) = self {
            let mut voice = ReedBoreVoice::with_line(duct, regional.inlet_gas(), reed,
                PlateBank::default(), distance, rate, |_| {
                    regional.characteristic_line_realized(termination, rate, samples, None, None)
                        .map_err(map_drive)
                }).map_err(ReedPerformanceError::Model)?;
            let grid = voice.line.history().len();
            if let Some(mut shunt) = regional
                .register_hole_shunt_line(rate, grid, fs_duct::WallModel::from_pin(None))
                .map_err(ReedPerformanceError::Outlet)?
            {
                shunt
                    .prime(voice.p_plus_prev)
                    .map_err(|error| ReedPerformanceError::Outlet(OutletError::Transfer(error)))?;
                voice.set_register_shunt(shunt);
            }
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
                molecular_hole_observers: None,
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

    /// Molecular side-mouth receivers and their individual exterior laws, in
    /// authored open/vented-hole order. Closed cavities have no radiating mouth.
    /// Each filter acts before summing pressure with the terminal contribution.
    pub fn molecular_open_hole_observers(
        &self,
    ) -> impl Iterator<Item = (OutletReceiverInfo, ReedObserverReceipt)> + '_ {
        self.terminal.iter().flat_map(|terminal| {
            terminal.molecular_hole_observers.iter().flat_map(|observers| {
                terminal.holes.iter().zip(observers).map(|(hole, observer)| {
                    let info = hole.info();
                    (info, ReedObserverReceipt {
                        filter: *observer.filter.receipt(),
                        range_spread_m: Some(range_spread(info.receiver.position_m, info.radius_m)),
                        consumer_delay_s: None,
                        internal_part_unpropagated: false,
                    })
                })
            })
        })
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
            if let Some(observers) = &mut terminal.molecular_hole_observers {
                sum = self.observer.process(sum);
                for (hole, observer) in terminal.holes.iter_mut().zip(observers) {
                    match hole.step(p_plus) {
                        Ok(pressure) => sum += observer.filter.process(pressure),
                        Err(error) => {
                            terminal.failure = Some(error.clone());
                            return Err(ReedPerformanceError::Outlet(error));
                        }
                    }
                }
            } else {
                for hole in &mut terminal.holes {
                    match hole.step(p_plus) {
                        Ok(pressure) => sum += pressure,
                        Err(error) => {
                            terminal.failure = Some(error.clone());
                            return Err(ReedPerformanceError::Outlet(error));
                        }
                    }
                }
                sum = self.observer.process(sum);
            }
            *slot = sum;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod molecular_tests {
    use super::*;
    use fs_duct::compose::ComposeOptions;
    use fs_duct::molecular::MolecularWallResolution;
    use fs_material::gas::GasSpec;
    use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
    use fs_thermochem::SpeciesId;

    fn gas(temperature_k: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap()
    }

    fn point(exterior: GasState) -> MolecularDrivingPoint {
        // Synthetic reservoir data exercise source coupling, not measured air.
        let data = "declared 1 synthetic 1 vibration 1 1300 1 power 2 293.15 0 200 500";
        let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
        let medium = RelaxingMedium::new(gas(293.15),
            vec![(SpeciesId::new("synthetic").unwrap(), 1.0)]).unwrap();
        MolecularDrivingPoint::new(&Duct { segments: vec![Segment::Cylinder {
            radius: 0.001, length: 0.15,
        }] }, law.evaluate(&medium).unwrap(), Some(MolecularWallResolution {
            poles_per_field: 16, band_rad_s: [1e-4, 30000.0], maximum_relative_error: 0.05,
        }), fs_duct::ConeSlices::Fixed(8), ComposeOptions::default(), exterior).unwrap()
    }

    fn reed() -> BeatingReed {
        BeatingReed { rest_opening_m: 4e-4, width_m: 0.003, closing_pressure_pa: 6000.0,
            blowing_pressure_pa: 2800.0, attack_s: 0.008,
            mass_kg: 0.0, stiffness_n_m: 0.0, damping_ratio: 0.0 }
    }

    fn receiver() -> CircularOutletReceiver {
        CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 2,
            angular_points: 8, maximum_frequency_hz: 800.0 }
    }

    #[test]
    fn molecular_observer_uses_the_same_prime_and_leaves_reed_feedback_unchanged() {
        let gas = gas(293.15);
        let source = point(gas);
        let (mut proxy, no_terminal) = ObservationInput::Legacy(1.0)
            .build_molecular(source.clone(), reed(), Termination::FlangedOpen,
                8000, 1024, gas.sound_speed, Some(DeclaredObserverLaw::Lossless), &[]).unwrap();
        assert!(no_terminal.is_none());
        let (mut physical, terminal) = ObservationInput::Baffled {
            receiver: receiver(), exterior: Some(gas),
        }.build_molecular(source, reed(), Termination::FlangedOpen,
            8000, 1024, gas.sound_speed, Some(DeclaredObserverLaw::Lossless), &[]).unwrap();
        let mut terminal = terminal.unwrap();
        assert!(terminal.holes.is_empty());
        let mut expected = physical.molecular_driving_point().unwrap()
            .baffled_outlet(8000, physical.line.history().len(), receiver()).unwrap();
        expected.step(proxy.p_plus_prev).unwrap();
        let mut moved = false;
        for _ in 0..512 {
            proxy.step_block(&mut [0.0]).unwrap();
            physical.step_block(&mut [0.0]).unwrap();
            let actual = terminal.pressure.step(physical.p_plus_prev).unwrap();
            let original = expected.step(proxy.p_plus_prev).unwrap();
            assert_eq!(actual.to_bits(), original.to_bits());
            assert_eq!(physical.p_plus_prev.to_bits(), proxy.p_plus_prev.to_bits());
            assert_eq!(physical.reed_y.to_bits(), proxy.reed_y.to_bits());
            assert_eq!(physical.reed_v.to_bits(), proxy.reed_v.to_bits());
            moved |= actual != 0.0;
        }
        assert!(moved);
        assert_eq!(physical.line.history(), proxy.line.history());
    }

    #[test]
    fn molecular_observer_refuses_a_different_load_or_undeclared_proxy_medium() {
        let source = point(gas(320.0));
        assert!(ObservationInput::Legacy(1.0).build_molecular(
            source.clone(), reed(), Termination::FlangedOpen,
            8000, 1024, gas(320.0).sound_speed, Some(DeclaredObserverLaw::Lossless), &[]).is_err());
        for (exterior, termination) in [
            (Some(gas(293.15)), Termination::FlangedOpen),
            (Some(gas(320.0)), Termination::UnflangedOpen),
        ] {
            assert!(ObservationInput::Baffled { receiver: receiver(), exterior }
                .build_molecular(source.clone(), reed(), termination,
                    8000, 1024, gas(320.0).sound_speed, Some(DeclaredObserverLaw::Lossless), &[]).is_err());
        }
    }
}
