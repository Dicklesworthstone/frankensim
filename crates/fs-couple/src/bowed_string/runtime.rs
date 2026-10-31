//! Persistent bowed-string stepping. The same exact-ZOH/friction loop serves
//! batch runs and interactive performance; stopping a block never resets a mode.

use super::{
    BowGesture, BowedRunConfig, BowedRunError, BowedStringCard, FrictionIsland,
    Termination, frame_total_energy, unit_shape, unit_shape_slope_at_bridge,
};
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAcousticWorkspace, advance_exact_zoh,
};
use crate::stribeck_friction::{
    ContactBranch, ImplicitFrictionError, StribeckLoadLine, solve_implicit_friction,
};
use crate::render::plate::MAX_PLATE_RETARDED_STEPS;
use crate::thin_plate::{CompactBody, HeldStep};
use fs_exec::CancelGate;

/// Sample-clock binding of physical bow gesture tracks to this retained state.
pub mod schedule;

/// Material-derived strings with reciprocal interior body connections.
pub mod reciprocal;

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
    /// Compact-body observer pressure [Pa] at this sample's end, emitted r/c
    /// earlier, when a body is attached. Zero until the first sound arrives.
    pub radiated_pressure_pa: Option<f64>,
    /// Endpoint string modal energy [J]; does not include soundboard energy.
    pub total_modal_energy_j: f64,
    /// Contact receipt of this sample's last substep.
    pub contact: BowContactReceipt,
    /// Substeps of this sample that ended stuck.
    pub stuck_substeps: u32,
}

/// Disclosure of one held-force contact step (one substep).
///
/// `roots` counts the intersections of the discrete load line with the
/// friction graph, the static cone included; `selected` is the zero-based
/// index of the committed one by slip speed, cone first. A lifted bow has no
/// contact and reports zero roots. The viscous falsifier has one root.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BowContactReceipt {
    /// Branch after the step; `Fresh` while lifted.
    pub branch: ContactBranch,
    /// Intersections found, 0..=3.
    pub roots: u8,
    /// Index of the committed intersection.
    pub selected: u8,
}

/// State of the original bowed-string physics across audio callbacks.
///
/// The initial configuration's `steps` is a batch-output length, not a limit on
/// this runtime. `max_block` bounds one callback and hence cancellation latency
/// in samples. Substeps per sample remain explicitly selected by the caller.
/// The string's candidate states reuse a workspace; a plate retains its own
/// allocation behavior. No blanket real-time guarantee is implied.
///
/// Sticking and sliding use the same exact-ZOH port response. One set-valued
/// step per substep ([`StribeckLoadLine`]) chooses endpoint no-slip within the
/// static capacity or the unregularized kinetic law at the endpoint, keeping the
/// contact's current branch where several roots exist.
/// The velocity response is cached from the shared modal transition, with no
/// state subtraction to recover small mobilities. This is not a continuous-time
/// constraint, a continuous-time passivity certificate or a new friction law.
/// The one-way
/// rigid-bridge/body approximation is unchanged. Control changes are physical
/// inputs, not mode resets or frequency assignments.
pub struct BowedStringState {
    card: BowedStringCard,
    island: FrictionIsland,
    kinetic_law: Option<fs_tribo::FrictionLaw>,
    load_line: Option<StribeckLoadLine>,
    branch: ContactBranch,
    receipts: Vec<BowContactReceipt>,
    gesture: BowGesture,
    model: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    shapes_at_bow: Vec<f64>,
    // Exact-ZOH endpoint velocity row [coefficient of q, v, held modal force].
    // Material, basis and substep are immutable; a moving bow changes only phi.
    contact_velocity_rows: Vec<[f64; 3]>,
    bridge_force_factors: Vec<f64>,
    bridge_viscous_factors: Vec<f64>,
    plate: Option<RadiatingBody>,
    listener_m: f64,
    forces: Vec<f64>,
    dt: f64,
    sub_dt: f64,
    subsamples: usize,
    peak_energy_j: f64,
    samples: u64,
    max_block: usize,
    poisoned: bool,
}

fn invalid(what: &'static str) -> BowedRunError {
    BowedRunError::InvalidRequest { what }
}

/// The bridge-driven body under its declared gas. The gas's baffled-piston
/// load acts on the body's motion, and the listener hears it at t - r/c,
/// replayed exactly from the recorded exact-ZOH step (the plate renderer's
/// scheme).
struct RadiatingBody {
    body: CompactBody,
    density: f64,
    /// Whole samples of the r/c delay.
    delay_whole: u64,
    /// Time into step n - delay_whole that the listener hears at sample n [s].
    retarded_elapsed_s: f64,
    /// Ring of the last delay_whole + 1 held steps.
    history: Vec<HeldStep>,
}

fn friction_refusal<E>(error: ImplicitFrictionError<E>, source: impl FnOnce(E) -> BowedRunError) -> BowedRunError {
    match error {
        ImplicitFrictionError::Constitutive(cause) => source(cause),
        ImplicitFrictionError::InvalidInput(what) => BowedRunError::Friction(what),
        ImplicitFrictionError::NonDissipative => BowedRunError::Friction("implicit contact would drive relative slip"),
        ImplicitFrictionError::NoConvergence => BowedRunError::Friction("implicit contact did not resolve the declared sliding law"),
    }
}

fn admit_contact_force(law: &fs_tribo::FrictionLaw, velocity_m_s: f64, normal_n: f64)
    -> Result<(), BowedRunError>
{
    let capacity = law.static_coefficient().map_err(BowedRunError::Interface)? * normal_n;
    if !capacity.is_finite() {
        return Err(BowedRunError::Friction("static friction capacity must be finite"));
    }
    if velocity_m_s != 0.0 {
        let force = law.kinetic_coefficient(velocity_m_s.abs())
            .map_err(BowedRunError::Interface)? * normal_n;
        if !force.is_finite() {
            return Err(BowedRunError::Friction("kinetic friction force must be finite"));
        }
    }
    Ok(())
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
        let kinetic_law = config.island.kinetic_law()?;
        config.island.check_state(config.gesture.v_bow_m_s, config.gesture.normal_force_n)?;
        if let Some(law) = &kinetic_law {
            admit_contact_force(law, config.gesture.v_bow_m_s, config.gesture.normal_force_n)?;
        }
        let plate = match &config.termination {
            Termination::PlateOnePort { body, ambient } => {
                if !(ambient.density.is_finite() && ambient.density > 0.0) {
                    return Err(BowedRunError::InvalidAmbientDensity { density_kg_m3: ambient.density });
                }
                if !(config.listener_m.is_finite() && config.listener_m > 0.0) {
                    return Err(BowedRunError::InvalidListenerDistance { distance_m: config.listener_m });
                }
                if !(ambient.sound_speed.is_finite() && ambient.sound_speed > 0.0) {
                    return Err(invalid("ambient sound speed must be positive and finite"));
                }
                let mut body = body.as_ref().clone();
                // A silent body (zero modal area) has no face to load.
                if body.area_m2 != 0.0 {
                    body.attach_piston_load(ambient, config.card.sample_rate_hz)
                        .map_err(|error| BowedRunError::Radiation(error.to_string()))?;
                }
                // Output n is the pressure at t = (n+1) dt. Its source time is
                // (n+1) dt - r/c = t_{n-m} + (1-f) dt for a delay of m + f samples.
                let rate = f64::from(config.card.sample_rate_hz);
                let delay = config.listener_m * rate / ambient.sound_speed;
                let rows = delay.floor() + 1.0;
                if !(rows.is_finite() && rows <= MAX_PLATE_RETARDED_STEPS as f64) {
                    return Err(BowedRunError::Radiation(
                        "listener distance exceeds the retarded body history budget".to_string()));
                }
                Some(RadiatingBody {
                    body,
                    density: ambient.density,
                    delay_whole: delay.floor() as u64,
                    retarded_elapsed_s: (1.0 - (delay - delay.floor())) / rate,
                    history: vec![HeldStep::REST; rows as usize],
                })
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
        // Probe the linear transition once, rather than copying its oscillator
        // equations or subtracting two nearly equal moving-state predictions.
        // The unit-force response starts at ZERO state; it is not inverse mass.
        let contact_velocity_rows: Vec<[f64; 3]> = model.modes().iter().map(|&mode| {
            let velocity = |q, v, force| advance_exact_zoh(mode, ModalAcousticState {
                displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v,
            }, force, sub_dt).velocity_m_sqrt_kg_per_s;
            [velocity(1.0, 0.0, 0.0), velocity(0.0, 1.0, 0.0), velocity(0.0, 0.0, 1.0)]
        }).collect();
        if contact_velocity_rows.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("bow contact response is not finite"));
        }
        let peak_energy_j = frame_total_energy(&model);
        Ok(Self {
            card: card.clone(), island: config.island.clone(), kinetic_law,
            load_line: None, branch: ContactBranch::Fresh,
            receipts: Vec::with_capacity(subsamples),
            gesture: config.gesture, model, workspace, shapes_at_bow, contact_velocity_rows,
            bridge_force_factors, bridge_viscous_factors, plate,
            listener_m: config.listener_m, forces: vec![0.0; card.mode_count],
            dt, sub_dt, subsamples, peak_energy_j,
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

    /// Replace the string's modal state through the modal model's transactional
    /// restart seam, for example to begin from a known motion. Gesture and clock
    /// are unchanged; the contact history restarts fresh. A refusal leaves the
    /// voice unchanged.
    pub fn restore_string_states(&mut self, states: &[ModalAcousticState]) -> Result<(), BowedRunError> {
        self.check_live()?;
        self.model.restore_states(states).map_err(BowedRunError::Model)?;
        self.branch = ContactBranch::Fresh;
        self.peak_energy_j = self.peak_energy_j.max(frame_total_energy(&self.model));
        Ok(())
    }

    /// Contact receipts of every substep of the last completed sample, in order.
    #[must_use]
    pub fn last_contact_receipts(&self) -> &[BowContactReceipt] { &self.receipts }

    /// Whether this voice has an actual pressure observer.
    #[must_use]
    pub fn has_radiation(&self) -> bool { self.plate.is_some() }

    /// Set physical bow inputs between samples. Signed speed permits reversals;
    /// zero normal force explicitly LIFTS the bow and bypasses friction, including
    /// the viscous falsifier. It does not silence or reset the string or body.
    ///
    /// A speed, station or load change keeps the contact branch; the next substep
    /// evaluates it on the new load line. Lifting or landing the bow clears that
    /// history, so a new touch starts fresh. Refusals leave all state untouched.
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
            if let Some(law) = &self.kinetic_law {
                admit_contact_force(law, velocity_m_s - velocity, normal_force_n)?;
            }
        }
        if normal_force_n == 0.0 || self.gesture.normal_force_n == 0.0 {
            self.branch = ContactBranch::Fresh;
        }
        // Repeating an unchanged station preserves the original reduction bits.
        if station.to_bits() != self.gesture.station_fraction.to_bits() {
            for (k, phi) in self.shapes_at_bow.iter_mut().enumerate() {
                *phi = unit_shape(k, station, self.card.length_m, self.card.linear_density_kg_m);
            }
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

    // Both static and sliding branches condense the SAME conjugate force and
    // velocity port. A material swap changes omega/zeta and phi, hence this
    // response, without any instrument or material-name branch.
    fn contact_response(&self) -> Result<(f64, f64), BowedRunError> {
        let mut free_velocity = 0.0;
        let mut mobility = 0.0;
        for ((state, phi), row) in self.model.states().iter()
            .zip(&self.shapes_at_bow).zip(&self.contact_velocity_rows)
        {
            free_velocity += phi * (row[0] * state.displacement_m_sqrt_kg
                + row[1] * state.velocity_m_sqrt_kg_per_s);
            mobility += phi * phi * row[2];
        }
        let free_relative = self.gesture.v_bow_m_s - free_velocity;
        if !free_relative.is_finite() || !mobility.is_finite() || mobility <= 0.0 {
            return Err(invalid("bow contact needs finite free slip and positive discrete mobility"));
        }
        Ok((free_relative, mobility))
    }

    fn implicit_viscous_traction(&self, viscous_n_s_per_m: f64, free_relative: f64, mobility: f64)
        -> Result<f64, BowedRunError>
    {
        let step = solve_implicit_friction(free_relative, mobility,
            |velocity| Ok::<_, BowedRunError>(viscous_n_s_per_m * velocity))
            .map_err(|error| friction_refusal(error, |source| source))?;
        Ok(step.traction_n)
    }

    // The owner's set-valued step: stick inside the static cone or slide on the
    // unregularized kinetic law, keeping the branch the contact is on. Source
    // applicability is checked at the committed endpoint, zero when sticking.
    fn set_valued_traction(&mut self, free_relative: f64, mobility: f64)
        -> Result<(f64, BowContactReceipt), BowedRunError>
    {
        let normal = self.gesture.normal_force_n;
        if !self.load_line.as_ref().is_some_and(|line| line.matches(normal, mobility)) {
            let law = self.kinetic_law.as_ref().ok_or(BowedRunError::Friction("missing shared dry law"))?;
            let line = StribeckLoadLine::new(law, normal, mobility)
                .map_err(|error| friction_refusal(error, BowedRunError::Interface))?;
            self.load_line = Some(line);
        }
        let line = self.load_line.as_ref().ok_or(BowedRunError::Friction("missing load line"))?;
        let step = line.step(free_relative, self.branch)
            .map_err(|error| friction_refusal(error, BowedRunError::Interface))?;
        self.island.check_state(step.relative_velocity_m_s, normal)?;
        self.branch = step.branch;
        Ok((step.traction_n, BowContactReceipt { branch: step.branch, roots: step.roots, selected: step.selected }))
    }

    fn advance_sample(&mut self) -> Result<BowedSample, BowedRunError> {
        self.receipts.clear();
        for _ in 0..self.subsamples {
            let v_str: f64 = self.model.states().iter().zip(&self.shapes_at_bow)
                .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s).sum();
            let v_rel = self.gesture.v_bow_m_s - v_str;
            // No contact exists while the bow is lifted. In particular the
            // interface's positive-load domain is not queried at zero load.
            let (traction, receipt) = if self.gesture.normal_force_n == 0.0 {
                self.branch = ContactBranch::Fresh;
                (0.0, BowContactReceipt::default())
            } else {
                self.island.check_state(v_rel, self.gesture.normal_force_n)?;
                let (free_relative, mobility) = self.contact_response()?;
                if let FrictionIsland::ViscousOnly { viscous_n_s_per_m } = self.island {
                    let traction = self.implicit_viscous_traction(viscous_n_s_per_m, free_relative, mobility)?;
                    let branch = ContactBranch::SlidingLow { positive: free_relative > 0.0 };
                    (traction, BowContactReceipt { branch, roots: 1, selected: 0 })
                } else {
                    self.set_valued_traction(free_relative, mobility)?
                }
            };
            self.receipts.push(receipt);
            for (q, phi) in self.forces.iter_mut().zip(&self.shapes_at_bow) { *q = traction * phi; }
            let stepped = self.model.step_duration_into(&self.forces, self.sub_dt, &mut self.workspace)
                .map_err(BowedRunError::LimitExceeded)?;
            let energy = stepped.total_modal_energy_j;
            if self.gesture.normal_force_n > 0.0
                && matches!(&self.island, FrictionIsland::InterfaceStribeck { .. } | FrictionIsland::Interface { .. })
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
        let (body_volume_velocity_m3_s, radiated_pressure_pa) = if let Some(plate) = self.plate.as_mut() {
            let n = self.samples;
            let rows = plate.delay_whole + 1;
            // Rows fit MAX_PLATE_RETARDED_STEPS, so the casts are exact.
            let (_, held) = plate.body.drive_held(bridge_force * plate.body.drive_participation, self.dt)
                .map_err(|error| BowedRunError::Radiation(error.to_string()))?;
            plate.history[(n % rows) as usize] = held;
            let velocity = plate.body.volume_velocity();
            let pressure = if n >= plate.delay_whole {
                // The row of step n - m; -m is congruent to 1 modulo m + 1.
                let step = plate.history[((n + 1) % rows) as usize];
                let acc = plate.body.acceleration_within(step, plate.retarded_elapsed_s)
                    .map_err(|error| BowedRunError::Radiation(error.to_string()))?;
                plate.body.radiate(acc, plate.density, self.listener_m)
            } else {
                0.0
            };
            if !(velocity.is_finite() && pressure.is_finite()) {
                return Err(BowedRunError::Radiation(
                    "compact observer produced non-finite volume velocity or pressure".to_string()));
            }
            (Some(velocity), Some(pressure))
        } else { (None, None) };
        let stuck = self.receipts.iter().filter(|r| r.branch == ContactBranch::Stuck).count();
        Ok(BowedSample {
            bow_point_velocity_m_s: endpoint_velocity,
            relative_velocity_m_s: self.gesture.v_bow_m_s - endpoint_velocity,
            bridge_force_n: bridge_force, body_volume_velocity_m3_s, radiated_pressure_pa,
            total_modal_energy_j: frame_total_energy(&self.model),
            contact: self.receipts.last().copied().unwrap_or_default(),
            stuck_substeps: u32::try_from(stuck).map_err(|_| invalid("bow substep count exceeds u32"))?,
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

#[cfg(test)]
mod exact_stiction_tests {
    use super::*;
    use crate::stribeck_friction::StribeckFriction;

    // The unregularized kinetic traction a sliding set-valued step commits.
    fn kinetic_traction(cfg: &BowedRunConfig, slip: f64) -> f64 {
        let law = cfg.island.kinetic_law().unwrap().unwrap();
        slip.signum() * law.kinetic_coefficient(slip.abs()).unwrap() * cfg.gesture.normal_force_n
    }

    // Numerical mechanics fixtures, not claims of measured instrument timbre.
    fn config(length: f64, mu: f64, tension: f64, bending: f64) -> BowedRunConfig {
        BowedRunConfig {
            card: BowedStringCard {
                length_m: length, tension_n: tension, linear_density_kg_m: mu,
                bending_stiffness_n_m2: bending, viscous_bending_n_m2_s: 0.0,
                mode_count: 12, zetas: vec![0.12; 12], sample_rate_hz: 48_000,
            },
            island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.3, 0.1, 0.1).unwrap()),
            gesture: BowGesture::admit(0.01, 100.0, 0.13).unwrap(),
            steps: 0, subsamples: 1, termination: Termination::Rigid, listener_m: 1.0,
        }
    }

    #[test]
    fn exact_stiction_tracks_the_bow_across_string_scales_and_substeps() {
        for (length, mu, tension, bending) in [
            (0.328, 0.003, 50.0, 1.0e-5),
            (0.69, 0.012, 90.0, 1.0e-4),
            (1.04, 0.04, 140.0, 1.0e-3),
        ] {
            for subsamples in [1, 2, 8] {
                let mut cfg = config(length, mu, tension, bending);
                cfg.subsamples = subsamples;
                let mut voice = BowedStringState::new(&cfg, 64).unwrap();
                for _ in 0..64 {
                    let frame = voice.step().unwrap();
                    assert_eq!(voice.branch, ContactBranch::Stuck);
                    assert!(frame.relative_velocity_m_s.abs() < 2.0e-11,
                        "endpoint slip {} with {subsamples} substeps", frame.relative_velocity_m_s);
                }
            }
        }
    }

    #[test]
    fn exact_stiction_reverses_and_release_preserves_free_motion() {
        let cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        let mut voice = BowedStringState::new(&cfg, 1).unwrap();
        voice.step().unwrap();
        // A held contact follows the bow through rest and reversal: the cone
        // still holds, so the branch is kept across the speed changes.
        for velocity in [0.0, -0.01] {
            voice.set_bow(velocity, 100.0, 0.13).unwrap();
            let frame = voice.step().unwrap();
            assert_eq!(voice.branch, ContactBranch::Stuck);
            assert!(frame.relative_velocity_m_s.abs() < 2.0e-11);
        }
        let mut free = voice.model.clone();
        free.step_duration(&vec![0.0; cfg.card.mode_count], voice.sub_dt).unwrap();
        voice.set_bow(0.5, 0.0, 0.13).unwrap();
        voice.step().unwrap();
        assert_eq!(voice.branch, ContactBranch::Fresh);
        assert_eq!(voice.model.states(), free.states());
        assert!(voice.total_modal_energy_j() > 0.0);
    }

    #[test]
    fn exact_stiction_never_overrides_the_static_force_capacity() {
        let mut cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        cfg.gesture.normal_force_n = 1.0e-6;
        let mut voice = BowedStringState::new(&cfg, 1).unwrap();
        let required = voice.model.held_force_for_port_velocity(
            &voice.shapes_at_bow, cfg.gesture.v_bow_m_s, voice.sub_dt,
        ).unwrap();
        assert!(required > 0.8 * cfg.gesture.normal_force_n);
        let frame = voice.step().unwrap();
        assert_ne!(voice.branch, ContactBranch::Stuck);
        assert!(frame.relative_velocity_m_s > 0.009);
        // The actual load is the unregularized kinetic law at the NEW endpoint,
        // here inside the island's ODE ramp scale, which the bow never uses.
        let FrictionIsland::Stribeck(legacy) = &cfg.island else { panic!("legacy fixture") };
        assert!(frame.relative_velocity_m_s < legacy.regularization_m_s);
        let traction = kinetic_traction(&cfg, frame.relative_velocity_m_s);
        for (force, phi) in voice.forces.iter().zip(&voice.shapes_at_bow) {
            assert!((*force - traction * phi).abs() < 1.0e-12);
        }
    }

    #[test]
    fn exact_stiction_force_responds_to_material_mass_and_stiffness() {
        let cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        let mut heavier = cfg.clone();
        // Preserve modal frequencies while tripling mechanical impedance.
        heavier.card.linear_density_kg_m *= 3.0;
        heavier.card.tension_n *= 3.0;
        heavier.card.bending_stiffness_n_m2 *= 3.0;
        let mut a = BowedStringState::new(&cfg, 1).unwrap();
        let mut b = BowedStringState::new(&heavier, 1).unwrap();
        a.step().unwrap();
        b.step().unwrap();
        let force_a = a.forces[0] / a.shapes_at_bow[0];
        let force_b = b.forces[0] / b.shapes_at_bow[0];
        assert!((force_b / force_a - 3.0).abs() < 1.0e-10);
        assert!((b.total_modal_energy_j() / a.total_modal_energy_j() - 3.0).abs() < 1.0e-9);
    }

    #[test]
    fn implicit_sliding_resolves_the_actual_endpoint_across_string_scales() {
        for (length, mu, tension, bending) in [
            (0.328, 0.003, 50.0, 1.0e-5),
            (0.69, 0.012, 90.0, 1.0e-4),
            (1.04, 0.04, 140.0, 1.0e-3),
        ] {
            let mut cfg = config(length, mu, tension, bending);
            cfg.gesture = BowGesture::admit(0.2, 0.02, 0.13).unwrap();
            let mut voice = BowedStringState::new(&cfg, 64).unwrap();
            for sample in 0..64 {
                // Exercise a geometry gesture too; it changes both conjugate
                // projections without resetting the cached material response.
                if sample == 32 { voice.set_bow(0.2, 0.02, 0.19).unwrap(); }
                let frame = voice.step().unwrap();
                assert_ne!(voice.branch, ContactBranch::Stuck);
                let expected = kinetic_traction(&cfg, frame.relative_velocity_m_s);
                for (force, phi) in voice.forces.iter().zip(&voice.shapes_at_bow) {
                    assert!((*force - expected * phi).abs() < 1.0e-10);
                }
            }
        }
    }

    #[test]
    fn hyperbolic_bow_contact_matches_closed_form_and_shared_owner() {
        let law = fs_tribo::FrictionLaw::Hyperbolic {
            static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: 0.2,
        };
        for (length, mu, tension, bending) in [
            (0.328, 0.003, 50.0, 1.0e-5),
            (0.69, 0.012, 90.0, 1.0e-4),
            (1.04, 0.04, 140.0, 1.0e-3),
        ] {
            for speed in [-0.35_f64, 0.35] {
                let mut cfg = config(length, mu, tension, bending);
                cfg.island = FrictionIsland::Constitutive(law.clone());
                cfg.gesture = BowGesture::admit(speed, 0.02, 0.13).unwrap();
                let mut voice = BowedStringState::new(&cfg, 1).unwrap();
                // Measure mobility by a unit-force transition, independently
                // of the cached contact row. On the positive-speed branch,
                // (v - D)(v + v0) + m N K = 0, D = |v_free| - m N mu_k.
                let mut unit = voice.model.clone();
                unit.step_duration(&voice.shapes_at_bow, voice.sub_dt).unwrap();
                let mobility: f64 = unit.states().iter().zip(&voice.shapes_at_bow)
                    .map(|(state, phi)| phi * state.velocity_m_sqrt_kg_per_s).sum();
                let d = speed.abs() - mobility * 0.02 * 0.4;
                let slip = 0.5 * (d - 0.2 + ((d + 0.2).powi(2)
                    - 4.0 * mobility * 0.02 * 0.08).sqrt());
                assert!(speed.abs() > mobility * 0.02 * 0.8);
                let expected_force = speed.signum() * 0.02 * (0.4 + 0.08 / (slip + 0.2));
                let (free, cached_mobility) = voice.contact_response().unwrap();
                let shared = StribeckLoadLine::new(&law, 0.02, cached_mobility).unwrap()
                    .step(free, ContactBranch::Fresh).unwrap();
                let frame = voice.step().unwrap();
                assert!((frame.relative_velocity_m_s - speed.signum() * slip).abs() < 1.0e-12);
                assert_eq!(frame.contact, BowContactReceipt {
                    branch: shared.branch, roots: shared.roots, selected: shared.selected,
                });
                for (force, phi) in voice.forces.iter().zip(&voice.shapes_at_bow) {
                    assert!((*force - expected_force * phi).abs() < 1.0e-12);
                    assert_eq!(*force, shared.traction_n * phi);
                }
            }
        }
    }

    #[test]
    fn shared_law_keeps_legacy_gaussian_output_independent_of_ode_ramp() {
        let mut cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        let mut voices = Vec::new();
        for ramp in [1.0e-6, 10.0] {
            cfg.island = FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.3, 0.1, ramp).unwrap());
            voices.push(BowedStringState::new(&cfg, 1).unwrap());
        }
        cfg.island = FrictionIsland::Constitutive(fs_tribo::FrictionLaw::Stribeck {
            static_mu: 0.8, kinetic_mu: 0.3, characteristic_speed: 0.1, viscous_per_speed: 0.0,
        });
        voices.push(BowedStringState::new(&cfg, 1).unwrap());
        for (speed, load) in [(0.2, 0.02), (0.01, 100.0), (-0.2, 0.02), (0.0, 0.0)] {
            for voice in &mut voices { voice.set_bow(speed, load, 0.13).unwrap(); }
            for _ in 0..32 {
                let expected = voices[0].step().unwrap();
                for voice in &mut voices[1..] { assert_eq!(voice.step().unwrap(), expected); }
            }
        }
    }

    #[test]
    fn hyperbolic_interface_retains_applicability_and_static_release_behavior() {
        use fs_tribo::{ApplicabilityRange as Range, DryFrictionApplicability, DryInterfaceSystemCard,
            FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef};
        let law = FrictionLaw::Hyperbolic {
            static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: 0.2,
        };
        let interface = |min_speed| FrictionIsland::Interface {
            card: Box::new(DryInterfaceSystemCard::new(
                InterfaceSystemRef::new("bow/string", "unworn", "synthetic-test",
                    InputAuthority::SyntheticFixture, InterfaceMedium::Dry).unwrap(),
                law.clone(), DryFrictionApplicability::new(
                    Range::new(280.0, 320.0).unwrap(), Range::new(0.0, 2.0e8).unwrap(),
                    Range::new(min_speed, 1.0).unwrap()).unwrap(),
            ).unwrap()), temperature_kelvin: 300.0, nominal_contact_area_m2: 1.0e-6,
        };
        let mut cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        cfg.gesture = BowGesture::admit(0.35, 0.02, 0.13).unwrap();
        cfg.island = FrictionIsland::Constitutive(law.clone());
        let mut direct = BowedStringState::new(&cfg, 1).unwrap();
        cfg.island = interface(0.0);
        let mut bound = BowedStringState::new(&cfg, 1).unwrap();
        for (speed, load) in [(0.35, 0.02), (0.01, 100.0), (-0.01, 100.0), (0.0, 0.0)] {
            direct.set_bow(speed, load, 0.13).unwrap();
            bound.set_bow(speed, load, 0.13).unwrap();
            for _ in 0..16 {
                let expected = direct.step().unwrap();
                assert_eq!(bound.step().unwrap(), expected);
                if load == 100.0 { assert_eq!(expected.contact.branch, ContactBranch::Stuck); }
                if load == 0.0 { assert_eq!(expected.contact.branch, ContactBranch::Fresh); }
            }
        }
        // The original interface owner still checks the actual solved slip,
        // not merely the input bow speed. A narrow source domain refuses it.
        cfg.island = interface(0.349999);
        let mut narrow = BowedStringState::new(&cfg, 1).unwrap();
        assert!(matches!(narrow.step(), Err(BowedRunError::Interface(
            fs_tribo::TriboError::OutsideApplicability { field: "slip_speed_mps", .. }))));
        cfg.island = interface(0.0);
        if let FrictionIsland::Interface { temperature_kelvin, .. } = &mut cfg.island {
            *temperature_kelvin = 330.0;
        }
        assert!(matches!(BowedStringState::new(&cfg, 1), Err(BowedRunError::Interface(
            fs_tribo::TriboError::OutsideApplicability { field: "temperature_kelvin", .. }))));
    }

    #[test]
    fn implicit_viscous_contact_matches_closed_form_without_explicit_reversal() {
        for (length, mu, tension, bending) in [
            (0.328, 0.003, 50.0, 1.0e-5),
            (0.69, 0.012, 90.0, 1.0e-4),
            (1.04, 0.04, 140.0, 1.0e-3),
        ] {
            let mut cfg = config(length, mu, tension, bending);
            cfg.card.mode_count = 1;
            cfg.card.zetas.truncate(1);
            let resistance = 100_000.0;
            cfg.island = FrictionIsland::ViscousOnly { viscous_n_s_per_m: resistance };
            cfg.gesture = BowGesture::admit(0.2, 1.0, 0.13).unwrap();
            let mut voice = BowedStringState::new(&cfg, 1).unwrap();
            // Independently obtain mobility by committing a real unit-force
            // transition, not by inspecting the cached response rows.
            let mut unit = voice.model.clone();
            unit.step_duration(&voice.shapes_at_bow, voice.sub_dt).unwrap();
            let mobility: f64 = unit.states().iter().zip(&voice.shapes_at_bow)
                .map(|(state, phi)| phi * state.velocity_m_sqrt_kg_per_s).sum();
            assert!(resistance * mobility > 1.0); // old explicit step reverses slip
            let expected = cfg.gesture.v_bow_m_s / (1.0 + resistance * mobility);
            let frame = voice.step().unwrap();
            assert!(frame.relative_velocity_m_s > 0.0);
            assert!((frame.relative_velocity_m_s - expected).abs() < 1.0e-12);
            let force = voice.forces[0] / voice.shapes_at_bow[0];
            assert!((force - resistance * expected).abs() < 2.0e-9);
        }
    }

    #[test]
    fn sliding_callback_partition_and_release_preserve_the_physical_state() {
        let mut cfg = config(0.328, 0.003, 50.0, 1.0e-5);
        cfg.gesture = BowGesture::admit(0.2, 0.02, 0.13).unwrap();
        let mut whole = BowedStringState::new(&cfg, 64).unwrap();
        let mut pieces = BowedStringState::new(&cfg, 64).unwrap();
        let mut a = [BowedSample::default(); 64];
        let mut b = a;
        whole.block(&mut a).unwrap();
        for block in b.chunks_mut(7) { pieces.block(block).unwrap(); }
        assert_eq!(a, b);
        let mut free = pieces.model.clone();
        free.step_duration(&vec![0.0; cfg.card.mode_count], pieces.sub_dt).unwrap();
        pieces.set_bow(-0.3, 0.0, 0.13).unwrap();
        pieces.step().unwrap();
        assert_eq!(pieces.model.states(), free.states());
        assert!(pieces.total_modal_energy_j() > 0.0);
    }
}
