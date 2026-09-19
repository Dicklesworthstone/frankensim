//! Persistent bowed-string stepping. The same exact-ZOH/friction loop serves
//! batch runs and interactive performance; stopping a block never resets a mode.

use super::{
    BowGesture, BowedRunConfig, BowedRunError, BowedStringCard, FrictionIsland,
    Termination, frame_total_energy, unit_shape, unit_shape_slope_at_bridge,
};
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use crate::stribeck_friction::StribeckFriction;
use crate::thin_plate::CompactBody;
use fs_exec::CancelGate;

/// Sample-clock binding of physical bow gesture tracks to this retained state.
pub mod schedule;

/// One completed audio sample. All mechanical channels describe the endpoint.
/// A rigid termination has no acoustic observer: its acoustic channels are None,
/// never velocity mislabeled as pressure.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BowedSample {
    /// Transverse velocity at the current bow station [m/s].
    pub bow_point_velocity_m_s: f64,
    /// Bow-minus-string velocity at that station [m/s].
    pub relative_velocity_m_s: f64,
    /// Transmitted bridge shear [N], including declared viscous bending.
    pub bridge_force_n: f64,
    /// Signed soundboard volume velocity [m^3/s], when a body is attached.
    pub body_volume_velocity_m3_s: Option<f64>,
    /// Compact-body observer pressure [Pa], when a body is attached.
    pub radiated_pressure_pa: Option<f64>,
    /// Endpoint string modal energy [J]; does not include soundboard energy.
    pub total_modal_energy_j: f64,
}

/// State of the original bowed-string physics across audio callbacks.
///
/// The initial configuration's `steps` is a batch-output length, not a limit on
/// this runtime. `max_block` bounds one callback and hence cancellation latency
/// in samples. Substeps per sample remain explicitly selected by the caller.
/// The string's candidate states reuse a workspace; a plate retains its own
/// allocation behavior. No blanket real-time guarantee is implied.
///
/// The existing first-order stiction correction and one-way rigid-bridge/body
/// approximation are unchanged. Control changes are physical inputs, not mode
/// resets, frequency assignments, or a new integration scheme.
pub struct BowedStringState {
    card: BowedStringCard,
    island: FrictionIsland,
    stribeck: Option<StribeckFriction>,
    gesture: BowGesture,
    model: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    shapes_at_bow: Vec<f64>,
    w_point: f64,
    bridge_force_factors: Vec<f64>,
    bridge_viscous_factors: Vec<f64>,
    plate: Option<(CompactBody, f64)>,
    listener_m: f64,
    forces: Vec<f64>,
    dt: f64,
    sub_dt: f64,
    subsamples: usize,
    capture_tol: f64,
    stuck: bool,
    prev_v_rel: f64,
    peak_energy_j: f64,
    samples: u64,
    max_block: usize,
    poisoned: bool,
}

fn invalid(what: &'static str) -> BowedRunError {
    BowedRunError::InvalidRequest { what }
}

impl BowedStringState {
    /// Admit once, beginning at rest. A zero-length batch still admits a usable
    /// voice. Gesture, interface, material-card and acoustic limits are the same
    /// as the original one-shot path. No sample is advanced here.
    pub fn new(config: &BowedRunConfig, max_block: usize) -> Result<Self, BowedRunError> {
        if max_block == 0 {
            return Err(invalid("bowed callback capacity must be positive"));
        }
        config.card.validate()?;
        BowGesture::admit(
            config.gesture.v_bow_m_s,
            config.gesture.normal_force_n,
            config.gesture.station_fraction,
        ).map_err(BowedRunError::Gesture)?;
        let stribeck = config.island.stribeck()?;
        config.island.check_state(config.gesture.v_bow_m_s, config.gesture.normal_force_n)?;
        match &config.island {
            FrictionIsland::Stribeck(_) | FrictionIsland::InterfaceStribeck { .. } => {
                let law = stribeck.ok_or(BowedRunError::Friction("missing Stribeck law"))?;
                config.island.sliding_traction(law, 0.0, config.gesture.normal_force_n)?;
                if !(law.mu_static * config.gesture.normal_force_n).is_finite() {
                    return Err(BowedRunError::Friction("static friction capacity must be finite"));
                }
            }
            FrictionIsland::ViscousOnly { viscous_n_s_per_m } => {
                if !(viscous_n_s_per_m.is_finite() && *viscous_n_s_per_m >= 0.0) {
                    return Err(BowedRunError::Friction("viscous friction must be finite and nonnegative"));
                }
            }
        }
        let plate = match &config.termination {
            Termination::PlateOnePort { body, ambient } => {
                if !(ambient.density.is_finite() && ambient.density > 0.0) {
                    return Err(BowedRunError::InvalidAmbientDensity { density_kg_m3: ambient.density });
                }
                if !(config.listener_m.is_finite() && config.listener_m > 0.0) {
                    return Err(BowedRunError::InvalidListenerDistance { distance_m: config.listener_m });
                }
                Some((body.as_ref().clone(), ambient.density))
            }
            Termination::Rigid => None,
        };
        let card = &config.card;
        let mu = card.linear_density_kg_m;
        let modes = (0..card.mode_count).map(|k| ModalAcousticMode {
            angular_frequency_rad_s: card.mode_omega_rad_s(k + 1),
            damping_ratio: card.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        }).collect();
        let model = ModalAcousticTimeModel::try_new(
            card.sample_rate_hz, modes, ModalAcousticTimeBudget::audible_reference(),
        ).map_err(BowedRunError::Model)?;
        let workspace = ModalAcousticWorkspace::new(&model);
        let shapes_at_bow: Vec<f64> = (0..card.mode_count)
            .map(|k| unit_shape(k, config.gesture.station_fraction, card.length_m, mu)).collect();
        let w_point = shapes_at_bow.iter().map(|phi| phi * phi).sum::<f64>();
        if !w_point.is_finite() || w_point <= 0.0 {
            return Err(invalid("bow station has no finite positive modal inverse mass"));
        }
        let bridge_force_factors = (0..card.mode_count).map(|k| {
            let kappa = (k + 1) as f64 * core::f64::consts::PI / card.length_m;
            (card.tension_n + card.bending_stiffness_n_m2 * kappa.powi(2))
                * unit_shape_slope_at_bridge(k, card.length_m, mu)
        }).collect();
        let bridge_viscous_factors = (0..card.mode_count).map(|k| {
            let kappa = (k + 1) as f64 * core::f64::consts::PI / card.length_m;
            card.viscous_bending_n_m2_s * kappa.powi(2)
                * unit_shape_slope_at_bridge(k, card.length_m, mu)
        }).collect();
        let dt = model.sample_period_s();
        let subsamples = config.subsamples.max(1); // preserve legacy zero-as-one
        let sub_dt = dt / subsamples as f64;
        if !sub_dt.is_finite() || sub_dt <= 0.0 {
            return Err(invalid("bow substep must have a positive finite duration"));
        }
        let peak_energy_j = frame_total_energy(&model);
        Ok(Self {
            card: card.clone(), island: config.island.clone(), stribeck,
            gesture: config.gesture, model, workspace, shapes_at_bow, w_point,
            bridge_force_factors, bridge_viscous_factors, plate,
            listener_m: config.listener_m, forces: vec![0.0; card.mode_count],
            dt, sub_dt, subsamples,
            capture_tol: stribeck.map_or(0.0, |law| law.stiction_m_s.max(0.02)),
            stuck: false, prev_v_rel: f64::NAN, peak_energy_j,
            samples: 0, max_block, poisoned: false,
        })
    }

    /// Number of complete audio samples; control changes never advance this clock.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.samples }

    /// Admitted audio clock [Hz].
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 { self.card.sample_rate_hz }

    /// Largest admitted output block.
    #[must_use]
    pub const fn max_block_len(&self) -> usize { self.max_block }

    /// Current string modal energy [J]. Soundboard energy is not included.
    #[must_use]
    pub fn total_modal_energy_j(&self) -> f64 { frame_total_energy(&self.model) }

    /// Largest string energy observed, including substep endpoints and initial rest.
    #[must_use]
    pub const fn peak_modal_energy_j(&self) -> f64 { self.peak_energy_j }

    /// Whether this voice has an actual pressure observer.
    #[must_use]
    pub fn has_radiation(&self) -> bool { self.plate.is_some() }

    /// Set physical bow inputs between samples. Signed speed permits reversals;
    /// zero normal force explicitly LIFTS the bow and bypasses friction, including
    /// the viscous falsifier. It does not silence or reset the string or body.
    ///
    /// Changes in speed or station invalidate the old contact's capture history,
    /// not the modal state. A load change alone retains stiction until its force
    /// cap is evaluated on the next substep. Refusals leave all state untouched.
    /// The positive-load, legacy [`BowGesture::admit`] contract is unchanged.
    pub fn set_bow(&mut self, velocity_m_s: f64, normal_force_n: f64, station: f64)
        -> Result<(), BowedRunError>
    {
        self.check_live()?;
        if !normal_force_n.is_finite() || normal_force_n < 0.0 {
            return Err(invalid("bow control normal force must be finite and nonnegative"));
        }
        BowGesture::admit(velocity_m_s, normal_force_n.max(f64::MIN_POSITIVE), station)
            .map_err(BowedRunError::Gesture)?;
        // Compute prospective shapes and the current velocity BEFORE changing
        // anything. Existing scalar laws own all interface applicability checks.
        let mut w_point = 0.0;
        let mut velocity = 0.0;
        for (k, state) in self.model.states().iter().enumerate() {
            let phi = unit_shape(k, station, self.card.length_m, self.card.linear_density_kg_m);
            w_point += phi * phi;
            velocity += phi * state.velocity_m_sqrt_kg_per_s;
        }
        if !w_point.is_finite() || w_point <= 0.0 || !velocity.is_finite() {
            return Err(invalid("bow control station has non-finite or zero modal response"));
        }
        if normal_force_n > 0.0 {
            self.island.check_state(velocity_m_s - velocity, normal_force_n)?;
            if let Some(law) = self.stribeck {
                if !(law.mu_static * normal_force_n).is_finite() {
                    return Err(BowedRunError::Friction("static friction capacity must be finite"));
                }
                self.island.sliding_traction(law, velocity_m_s - velocity, normal_force_n)?;
            }
        }
        if velocity_m_s.to_bits() != self.gesture.v_bow_m_s.to_bits()
            || station.to_bits() != self.gesture.station_fraction.to_bits()
            || normal_force_n == 0.0 || self.gesture.normal_force_n == 0.0
        {
            self.stuck = false;
            self.prev_v_rel = f64::NAN;
        }
        // Repeating an unchanged station preserves the original reduction bits.
        if station.to_bits() != self.gesture.station_fraction.to_bits() {
            for (k, phi) in self.shapes_at_bow.iter_mut().enumerate() {
                *phi = unit_shape(k, station, self.card.length_m, self.card.linear_density_kg_m);
            }
            self.w_point = w_point;
        }
        self.gesture = BowGesture { v_bow_m_s: velocity_m_s, normal_force_n, station_fraction: station };
        Ok(())
    }

    fn check_live(&self) -> Result<(), BowedRunError> {
        if self.poisoned { Err(BowedRunError::Poisoned) } else { Ok(()) }
    }

    fn validate_request(&self, len: usize) -> Result<(), BowedRunError> {
        self.check_live()?;
        if len == 0 || len > self.max_block {
            return Err(invalid("bowed block must be nonempty and within its admitted capacity"));
        }
        let count = u64::try_from(len).map_err(|_| invalid("bowed sample count exceeds u64"))?;
        self.samples.checked_add(count).ok_or_else(|| invalid("bowed sample clock overflow"))?;
        Ok(())
    }

    /// Advance one complete audio sample. A physical refusal permanently poisons
    /// this voice: substeps/body motion cannot be presented as a resumable sample.
    pub fn step(&mut self) -> Result<BowedSample, BowedRunError> {
        self.validate_request(1)?;
        match self.advance_sample() {
            Ok(frame) => { self.samples += 1; Ok(frame) }
            Err(error) => { self.poisoned = true; Err(error) }
        }
    }

    /// Fill one admitted block of endpoint diagnostics. On physics failure,
    /// discard this block. Cancellation belongs BETWEEN these calls, not inside
    /// a partially completed sample. No accumulated history is copied per block.
    pub fn block(&mut self, out: &mut [BowedSample]) -> Result<(), BowedRunError> {
        self.validate_request(out.len())?;
        for slot in out { *slot = self.step()?; }
        Ok(())
    }

    /// Render actual observer pressure without retaining per-sample diagnostics.
    /// A rigid termination refuses before any state/output change.
    pub fn pressure_block(&mut self, out: &mut [f64]) -> Result<(), BowedRunError> {
        self.validate_request(out.len())?;
        if !self.has_radiation() {
            return Err(invalid("rigid bowed string has no pressure observer; attach a physical body"));
        }
        for slot in out {
            // The termination cannot change after admission.
            *slot = self.step()?.radiated_pressure_pa.expect("admitted body observer");
        }
        Ok(())
    }

    fn advance_sample(&mut self) -> Result<BowedSample, BowedRunError> {
        const CAPTURE_BASIN_M_S: f64 = 0.15;
        for _ in 0..self.subsamples {
            let v_str: f64 = self.model.states().iter().zip(&self.shapes_at_bow)
                .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s).sum();
            let v_rel = self.gesture.v_bow_m_s - v_str;
            // No contact exists while the bow is lifted. In particular the
            // interface's positive-load domain is not queried at zero load.
            let traction = if self.gesture.normal_force_n == 0.0 {
                self.stuck = false;
                self.prev_v_rel = f64::NAN;
                0.0
            } else {
                self.island.check_state(v_rel, self.gesture.normal_force_n)?;
                let flip_speed = self.prev_v_rel.abs().max(v_rel.abs());
                let flipped = self.prev_v_rel.is_finite() && !self.stuck
                    && self.prev_v_rel.signum() != v_rel.signum()
                    && flip_speed <= CAPTURE_BASIN_M_S;
                self.prev_v_rel = v_rel;
                match &self.island {
                    FrictionIsland::Stribeck(_) | FrictionIsland::InterfaceStribeck { .. } => {
                        let law = self.stribeck.ok_or(BowedRunError::Friction("missing Stribeck law"))?;
                        let hold_cap = law.mu_static * self.gesture.normal_force_n;
                        // Preserve the legacy first-order pinning correction;
                        // the exact-ZOH step below is not an exact no-slip solve.
                        let pin = || -> f64 {
                            let accel_hold: f64 = self.model.modes().iter().zip(self.model.states())
                                .zip(&self.shapes_at_bow).map(|((m, s), phi)| {
                                    phi * (2.0 * m.damping_ratio * m.angular_frequency_rad_s
                                        * s.velocity_m_sqrt_kg_per_s
                                        + m.angular_frequency_rad_s * m.angular_frequency_rad_s
                                            * s.displacement_m_sqrt_kg)
                                }).sum::<f64>() / self.w_point;
                            accel_hold + v_rel / (self.w_point * self.sub_dt)
                        };
                        if self.stuck {
                            let p = pin();
                            if p.abs() <= hold_cap { p } else {
                                self.stuck = false;
                                self.island.sliding_traction(law, v_rel, self.gesture.normal_force_n)?
                            }
                        } else if flipped || v_rel.abs() <= self.capture_tol {
                            let p = pin();
                            if p.abs() <= hold_cap { self.stuck = true; p } else {
                                self.island.sliding_traction(law, v_rel, self.gesture.normal_force_n)?
                            }
                        } else {
                            self.island.sliding_traction(law, v_rel, self.gesture.normal_force_n)?
                        }
                    }
                    FrictionIsland::ViscousOnly { viscous_n_s_per_m } => viscous_n_s_per_m * v_rel,
                }
            };
            for (q, phi) in self.forces.iter_mut().zip(&self.shapes_at_bow) { *q = traction * phi; }
            let stepped = self.model.step_duration_into(&self.forces, self.sub_dt, &mut self.workspace)
                .map_err(BowedRunError::LimitExceeded)?;
            let energy = stepped.total_modal_energy_j;
            if self.gesture.normal_force_n > 0.0
                && matches!(&self.island, FrictionIsland::InterfaceStribeck { .. })
            {
                let endpoint: f64 = self.model.states().iter().zip(&self.shapes_at_bow)
                    .map(|(state, phi)| phi * state.velocity_m_sqrt_kg_per_s).sum();
                self.island.check_state(self.gesture.v_bow_m_s - endpoint, self.gesture.normal_force_n)?;
            }
            self.peak_energy_j = self.peak_energy_j.max(energy);
        }
        let bridge_force: f64 = self.model.states().iter().zip(&self.bridge_force_factors)
            .zip(&self.bridge_viscous_factors).map(|((s, factor), viscous)| {
                factor * s.displacement_m_sqrt_kg + viscous * s.velocity_m_sqrt_kg_per_s
            }).sum();
        let endpoint_velocity: f64 = self.model.states().iter().zip(&self.shapes_at_bow)
            .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s).sum();
        if !bridge_force.is_finite() || !endpoint_velocity.is_finite() {
            return Err(invalid("bowed endpoint observation is non-finite"));
        }
        let (body_volume_velocity_m3_s, radiated_pressure_pa) = if let Some((body, density)) = self.plate.as_mut() {
            let acc = body.drive(bridge_force * body.drive_participation, self.dt)
                .map_err(|error| BowedRunError::Radiation(error.to_string()))?;
            let velocity = body.volume_velocity();
            let pressure = body.radiate(acc, *density, self.listener_m);
            if !(velocity.is_finite() && pressure.is_finite()) {
                return Err(BowedRunError::Radiation(
                    "compact observer produced non-finite volume velocity or pressure".to_string()));
            }
            (Some(velocity), Some(pressure))
        } else { (None, None) };
        Ok(BowedSample {
            bow_point_velocity_m_s: endpoint_velocity,
            relative_velocity_m_s: self.gesture.v_bow_m_s - endpoint_velocity,
            bridge_force_n: bridge_force, body_volume_velocity_m3_s, radiated_pressure_pa,
            total_modal_energy_j: frame_total_energy(&self.model),
        })
    }

    /// Render a diagnostic window with bounded host callbacks and a short final
    /// block. Cancellation reports the EXACT written prefix and leaves every
    /// suffix slot untouched. The same state resumes with an unrequested gate.
    pub fn render_under_gate(&mut self, gate: &CancelGate, out: &mut [BowedSample], block_len: usize)
        -> Result<BowedRenderOutcome, BowedRunError>
    {
        self.validate_request(block_len)?;
        let count = u64::try_from(out.len()).map_err(|_| invalid("bowed output exceeds the sample clock"))?;
        self.samples.checked_add(count).ok_or_else(|| invalid("bowed sample clock overflow"))?;
        let mut samples = 0;
        for block in out.chunks_mut(block_len) {
            if gate.is_requested() { return Ok(BowedRenderOutcome::Cancelled { samples }); }
            self.block(block)?;
            samples += block.len();
        }
        Ok(BowedRenderOutcome::Completed { samples })
    }
}

/// Samples written by one gated request, not the lifetime clock.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BowedRenderOutcome {
    /// All output samples, including a short final block, were completed.
    Completed { /// Written sample count.
        samples: usize },
    /// Cancelled before the next host callback; state is resumable.
    Cancelled { /// Written sample count.
        samples: usize },
}
