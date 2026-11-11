//! The lip on the shared nonlinear gas bore (bead `frankensim-giuhr`).
//!
//! The same outward-striking lip island as [`super::BrassVoice`] meters the
//! mouth into a bore whose gas obeys the conservative Euler equations of the
//! performance's own NASA-9 mixture ([`NonlinearBore`]). Wave steepening, the
//! Zwikker-Kosten wall memory and the baffled-piston bell all live in that one
//! field; no reflectance line, transfer function or fitted load is used.
//!
//! The lip channel is an isentropic nozzle of area `width * gap` (with the opt-in exponent row, `width * gap * (gap / rest gap)^(q - 1)`;
//! with the opt-in second lip mode, the gap is the rest gap plus both modes' displacements) from the
//! mouth (stagnation gas at the blowing pressure and the room temperature) to
//! the pressure on the bore's first face. Its jet mixes out over that face
//! without pressure recovery, the linear loop's Bernoulli law made
//! compressible. When the cup pressure exceeds the mouth's, gas leaves through
//! the same channel. At every Runge-Kutta stage the face pressure solves mass
//! continuity between the metered channel and the interior's own end wave
//! (this source's shock or fan). The mixed-out face cools by its kinetic head:
//! its static enthalpy plus `u^2 / 2` is the mouth's stagnation enthalpy, using
//! the same source's temperature-dependent heat capacity. A closed lip is a wall.
//!
//! An optional finite mouth plenum stores mass and internal energy. The
//! blowing-pressure control then feeds it through an authored opening, and
//! the evolving plenum pressure drives the lip. Feed and lip exchanges share
//! the bore's SSP-RK2 stages, including its retried trials.
//!
//! The lip steps once per gas step under the cup pressure of the previous
//! step's first stage, and the gap is held through both stages. The radiated
//! observer receives every mechanics sample's mean bell volume flow.
//!
//! Limits: the held feed uses the room's temperature and frozen composition
//! (no exhaled CO2). A finite mouth is rigid, adiabatic and acoustically compact;
//! it is not an anatomical lung model. The existing lip's swept-volume work
//! remains outside this approximation; gas-exchange conservation is not a
//! whole lip/gas passivity claim. Reverse flow discharges the cup's static
//! state; one held bore without valves; the bell is a rigid baffled piston.

use fs_duct::Duct;
use fs_duct::modal::MouthKernel;
use fs_duct::nonlinear::bore::{BoreSpec, NonlinearBore, WallGas, WallMemory, WallReceipt};
use fs_duct::nonlinear::{FaceFlux, FrozenGas, PortStage, Primitive};
use fs_gas::dynamics::Error as GasError;
use fs_phs::linear_exact::{LinearPropagator, LinearSplit, step_split};
use fs_phs::{PortHamiltonian, mass_spring_damper};

use super::{
    BrassBlockDiag, BrassControl, LipIsland, LipMode2, LipSecondMode, TensionWork, prepare_lip, spring_work_j,
    zero_crossing_lock_hz,
};
use crate::turbulence_noise::{
    ConstrictionNoise, ConstrictionNoiseCard, JetNoiseCard, JetNoiseCoverage, LipJetNoise,
};

mod finite_mouth;
pub use finite_mouth::{FiniteMouthReport, FiniteMouthSpec, FiniteMouthState};
use super::bell_shell::{ModalBellWall, ShellMaterial};
use super::bore_wall::{
    BRASS_DENSITY_KG_M3, BRASS_POISSON_RATIO, BRASS_YOUNG_PA, BellShellRecord, BoreWallRecord, ThinBrassWall,
};
use super::performance::BellShellInfo;
use finite_mouth::FiniteMouth;

/// How the nonlinear bore is cut, clocked and walled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearClock {
    /// Largest axial gas cell [m].
    pub maximum_cell_m: f64,
    /// Output sample rate [Hz].
    pub output_rate_hz: u32,
    /// Mechanics (observer) samples per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Gas steps per mechanics sample.
    pub gas_steps_per_mechanics: u32,
    /// Wall memory of every cell.
    pub wall: WallMemory,
}

impl NonlinearClock {
    /// Mechanics (observer) clock [Hz].
    #[must_use]
    pub fn mechanics_rate_hz(&self) -> u32 {
        self.output_rate_hz * self.mechanics_per_output
    }

    /// Fixed gas step [s].
    #[must_use]
    pub fn gas_step_s(&self) -> f64 {
        1.0 / (f64::from(self.mechanics_rate_hz()) * f64::from(self.gas_steps_per_mechanics))
    }
}

/// The face the lip channel and the bore agree on at one stage.
#[derive(Debug, Clone, Copy)]
struct LipFace {
    flux: FaceFlux,
    pressure_pa: f64,
    mass_flow_kg_s: f64,
}

/// The lip island driving a nonlinear gas bore.
pub struct NonlinearLip {
    bore: NonlinearBore,
    gas: FrozenGas,
    room: Primitive,
    lip_sys: PortHamiltonian,
    lip_linear: LinearSplit,
    lip_state: Vec<f64>,
    lip: LipIsland,
    tension: f64,
    tension_work: TensionWork,
    blow_pa: f64,
    mouth: Primitive,
    finite_mouth: Option<FiniteMouth>,
    clock: NonlinearClock,
    face_pa: f64,
    block_index: u64,
    diags: Vec<BrassBlockDiag>,
    /// Opt-in lip jet noise ([`LipJetNoise`]: Stevens' form with Verge's constant, an extrapolation); `None` adds nothing.
    jet_noise: Option<LipJetNoise>,
    /// Opt-in constriction noise ([`ConstrictionNoise`]: Oldham and Ukpoho's printed duct-constriction law in third-octave
    /// bands, an extrapolation to the lip); `None` adds nothing.
    constriction_noise: Option<ConstrictionNoise>,
    /// The bore gas's dynamic viscosity [Pa s] and density [kg/m^3], the engine's own (the wall binding's), read only by the
    /// opt-in noises: the jet noise's Reynolds gate and reference jet speed, the constriction noise's reference jet speed.
    wall_viscosity_pa_s: f64,
    reference_density_kg_m3: f64,
    /// Opt-in second lip mode ([`LipMode2`], the linear loop's own integrator); `None` leaves the lip exactly the one-mode lip.
    lip2: Option<LipSecondMode>,
    /// The opt-in flow-area exponent `q` ([`NonlinearLip::with_lip_area_exponent`]); exactly 1.0 is the plain `width x gap`.
    area_exponent: f64,
    /// The modal bell wall's admission receipt (`bore_wall shell-brass`), if bound.
    bell_shell: Option<BellShellInfo>,
}

impl NonlinearLip {
    /// Fill `duct` with the room gas at rest and attach the lip at its inlet.
    /// `cup_volume_m3` is added to the first gas cell.
    ///
    /// # Errors
    /// A named refusal from the lip island, the gas owner, the mesh, the wall
    /// binding or the clock.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: WallGas,
        lip: LipIsland,
        cup_volume_m3: f64,
        clock: NonlinearClock,
    ) -> Result<Self, String> {
        Self::new_with_mouth(duct, gas, room_temperature_k, room_pressure_pa, wall_gas, lip, cup_volume_m3, clock, None)
    }

    /// As [`Self::new`], with the bore's radiating end optionally the plane-mode equivalent of the complete flanged multimodal
    /// mouth (the `mouth_load modal` row): `mouth` is the causal kernel [`fs_duct::modal::plane_mouth_kernel`] built on this
    /// bore's own duct, with this bore's gas as the bore gas and the exterior it names (the same gas unless the card declares a `bore_gas` row), on this clock's gas
    /// step. `None` is [`Self::new`] bit for bit.
    ///
    /// # Errors
    /// As [`Self::new`], or a kernel that does not fit the bore.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_mouth(
        duct: &Duct,
        gas: FrozenGas,
        room_temperature_k: f64,
        room_pressure_pa: f64,
        wall_gas: WallGas,
        lip: LipIsland,
        cup_volume_m3: f64,
        clock: NonlinearClock,
        mouth: Option<&MouthKernel>,
    ) -> Result<Self, String> {
        if !(1..=16).contains(&clock.mechanics_per_output) || clock.gas_steps_per_mechanics == 0 {
            return Err(
                "nonlinear bore needs 1..=16 mechanics samples per output sample and at least one gas step per mechanics sample"
                    .into(),
            );
        }
        if !(lip.mass_kg > 0.0
            && lip.stiffness_n_m > 0.0
            && lip.damping_n_s_m >= 0.0
            && lip.width_m > 0.0
            && lip.rest_gap_m > 0.0
            && lip.face_area_m2 > 0.0)
        {
            return Err("lip island parameters must be positive".into());
        }
        let (wall_viscosity_pa_s, reference_density_kg_m3) = (wall_gas.dynamic_viscosity_pa_s, wall_gas.density_kg_m3);
        let bore = NonlinearBore::with_mouth_kernel(
            duct,
            gas.clone(),
            room_temperature_k,
            room_pressure_pa,
            Some(wall_gas),
            BoreSpec {
                maximum_cell_m: clock.maximum_cell_m,
                inlet_volume_m3: cup_volume_m3,
                step_s: clock.gas_step_s(),
                cfl: 0.45,
                wall: Some(clock.wall),
            },
            mouth,
        )
        .map_err(|e| format!("nonlinear bore refused: {e}"))?;
        let lip_sys = mass_spring_damper(lip.mass_kg, lip.stiffness_n_m, lip.damping_n_s_m)
            .map_err(|_| "lip island refused admission".to_string())?;
        let lip_linear = prepare_lip(&lip_sys).map_err(|e| e.to_string())?;
        let room = bore.room();
        Ok(Self {
            bore,
            gas,
            room,
            lip_sys,
            lip_linear,
            lip_state: vec![0.0, 0.0],
            lip,
            tension: 1.0,
            tension_work: TensionWork::default(),
            blow_pa: 0.0,
            mouth: room,
            finite_mouth: None,
            clock,
            face_pa: room.pressure,
            block_index: 0,
            diags: Vec::new(),
            jet_noise: None,
            constriction_noise: None,
            wall_viscosity_pa_s,
            reference_density_kg_m3,
            lip2: None,
            area_exponent: 1.0,
            bell_shell: None,
        })
    }

    /// Add the opt-in second lip mode before any gas step: the linear loop's own integrator ([`super::BrassVoice::with_lip_mode2`]
    /// is the reference: a second mass-spring-damper stepped under the SAME transmural force with the sign of its printed
    /// `1/mu`, the aperture being `rest gap + q1 + q2` clipped at zero). A lip without it steps exactly the one-mode lip. Its
    /// published numbers carry no pre-stress law, so it runs only at lip tension 1.0 and [`NonlinearLip::apply`] refuses any
    /// other `SetLipTension` by name.
    ///
    /// # Errors
    /// The lip's tension is not 1.0, or the mode is refused admission.
    pub fn with_lip_mode2(mut self, mode: LipMode2) -> Result<Self, String> {
        if self.tension != 1.0 {
            return Err("the second lip mode (lip_response_mode2) has no printed pre-stress law and runs only at lip tension 1.0".into());
        }
        let sys = mass_spring_damper(mode.mass_kg, mode.stiffness_n_m, mode.damping_n_s_m)
            .map_err(|_| "the second lip mode refused admission".to_string())?;
        let linear = prepare_lip(&sys).map_err(|e| e.to_string())?;
        self.lip2 = Some(LipSecondMode { mode, sys, linear, state: vec![0.0, 0.0] });
        Ok(self)
    }

    /// Bind the opt-in breathing wall (`bore_wall thin-brass`) before any gas step: every cell of the bore gets a radial wall
    /// oscillator of cartridge brass of the declared thickness and breathing-mode quality factor
    /// ([`NonlinearBore::bind_shell`], [`super::bore_wall`]). A lip without it steps exactly the rigid-walled bore, bit for bit.
    ///
    /// # Errors
    /// A thickness or quality factor the wall refuses, or a duct that is not this bore's own.
    pub fn with_bore_wall(mut self, duct: &Duct, record: BoreWallRecord) -> Result<Self, String> {
        let law = ThinBrassWall::new(record.thickness_m)?;
        self.bore
            .bind_shell(duct, &law, record.quality)
            .map_err(|e| format!("bore wall refused: {e}"))?;
        Ok(self)
    }

    /// Bind the opt-in modal bell wall (`bore_wall shell-brass`) before any gas step: the flaring end of the bore is meshed as a
    /// brass shell of revolution, its axisymmetric modes below the cut-off become a [`fs_duct::nonlinear::bore::WallDriver`], and
    /// every bell cell's pressure drives them while they sweep volume from every cell ([`super::bell_shell`]). A lip without it
    /// steps exactly the rigid-walled bore, bit for bit.
    ///
    /// # Errors
    /// A bore without a flaring bell, a shell the owner refuses, or a driver the bore refuses.
    pub fn with_bell_shell(mut self, duct: &Duct, record: BellShellRecord) -> Result<Self, String> {
        let faces = self.bore.solver().mesh().positions().to_vec();
        let wall = ModalBellWall::from_duct(
            duct,
            &faces,
            record,
            ShellMaterial {
                density_kg_m3: BRASS_DENSITY_KG_M3,
                young_pa: BRASS_YOUNG_PA,
                poisson: BRASS_POISSON_RATIO,
            },
            self.clock.gas_step_s(),
        )?;
        let frequencies = wall.frequencies_hz();
        self.bell_shell = Some(BellShellInfo {
            record,
            cells_covered: wall.covered_cells(),
            coupling_modes: frequencies.len(),
            modes_below_cutoff: wall.modes_below_cutoff(),
            first_covered_cell: wall.first_covered_cell().unwrap_or(0),
            lowest_mode_hz: frequencies.first().copied().unwrap_or(0.0),
            highest_mode_hz: frequencies.last().copied().unwrap_or(0.0),
        });
        self.bore
            .bind_wall_driver(Box::new(wall))
            .map_err(|e| format!("bell shell refused: {e}"))?;
        Ok(self)
    }

    /// What binding the modal bell wall admitted, or `None` without one.
    #[must_use]
    pub fn bell_shell_info(&self) -> Option<BellShellInfo> {
        self.bell_shell
    }

    /// Add the opt-in flow-area exponent `q` before any gas step: the aperture the lip channel meters through becomes
    /// `w h (h/H0)^(q-1)` for the lip gap `h`, the rest gap `H0` and the width `w`, as on the linear bore
    /// ([`super::BrassVoice::with_lip_area_exponent`] carries the label and the source: a bracket parameter, not a measured value).
    /// `q = 1` is the plain `w h` and is bit-identical to a lip without the row. Only the flow area changes.
    ///
    /// # Errors
    /// Unless `1 <= q <= 4`.
    pub fn with_lip_area_exponent(mut self, q: f64) -> Result<Self, String> {
        if !(q.is_finite() && (1.0..=4.0).contains(&q)) {
            return Err("the lip area exponent (lip_area_exponent) must be finite and in [1, 4]".into());
        }
        self.area_exponent = q;
        Ok(self)
    }

    // The flow aperture `w h (h/H0)^(q-1)` [m^2] for the lip gap `h`; exactly `w h` when the exponent is 1 or the lip is closed.
    fn flow_width_height(&self, gap: f64) -> f64 {
        let plain = self.lip.width_m * gap;
        if self.area_exponent == 1.0 || gap <= 0.0 {
            plain
        } else {
            plain * (gap / self.lip.rest_gap_m).powf(self.area_exponent - 1.0)
        }
    }

    /// Attach finite upstream gas storage before any gas step.
    ///
    /// The blowing-pressure control becomes the held feed pressure; it never
    /// resets the stored mass, energy or pressure. The actual mouth pressure
    /// drives the lip. A zero feed area gives a sealed, precharged supply.
    /// Every source state uses the bore's exact frozen caloric law.
    ///
    /// # Errors
    /// Invalid storage or source data, or repeated/late attachment.
    ///
    /// The opt-in noise rows stand beside the finite mouth, in either builder order. Their summed pressure is a SERIES offset
    /// between the stored gas and the lip face: both stages of the gas step offer the lip the stored gas at its own temperature
    /// and the pressure raised by the noise (an isothermal modulation, as the original reservoir is rebuilt at the room
    /// temperature), the feed nozzle and the stored inventory never see it, and the source's mechanical work is reported as
    /// [`FiniteMouthReport::series_work_j`] (rejected as heat by the isothermal modulation, so it is not drawn from the storage).
    pub fn with_finite_mouth(mut self, spec: FiniteMouthSpec) -> Result<Self, String> {
        if self.bore.solver().time().value() != 0.0 || self.finite_mouth.is_some() {
            return Err("finite mouth must be attached once before any gas step".into());
        }
        let (supply, primitive) = FiniteMouth::new(&self.gas, spec, self.mouth)
            .map_err(|e| format!("finite mouth refused: {e}"))?;
        self.finite_mouth = Some(supply);
        self.mouth = primitive;
        Ok(self)
    }

    /// Current finite upstream storage; absent for the original held reservoir.
    #[must_use]
    pub fn finite_mouth_state(&self) -> Option<FiniteMouthState> {
        self.finite_mouth.as_ref().map(|supply| supply.state(self.mouth))
    }

    /// Most recent accepted finite-mouth gas exchange, on the bore's clock.
    #[must_use]
    pub fn finite_mouth_report(&self) -> Option<&FiniteMouthReport> {
        self.finite_mouth.as_ref().and_then(FiniteMouth::report)
    }

    /// Attach the opt-in lip jet noise before any step: a series pressure source added to the blowing pressure that drives
    /// the lip (its force) and the mouth reservoir the lip channel meters from, of rms `gain rho U_j^2 / 2` inside the band
    /// `St_h = f h / U_j` in `[st_lo, st_hi]` (a band reaching the gas step's Nyquist frequency keeps only its part below
    /// it), shut while `Re_h = mdot / (w mu)` is below the card's onset or the jet flow is
    /// not forward (see [`LipJetNoise`] for the law, what it rests on and the label). It runs once per gas step on the
    /// previous gas step's metered mass flow and gap. A lip without it is bit-identical to one built before this existed; a
    /// gate that never opens adds exactly zero.
    ///
    /// # Errors
    /// The reason [`JetNoiseCard::validate`] gives.
    pub fn with_jet_noise(mut self, card: JetNoiseCard) -> Result<Self, String> {
        self.jet_noise = Some(LipJetNoise::new(card).map_err(str::to_string)?);
        Ok(self)
    }

    /// The lip jet noise's open-gate gas steps so far by how much of the band each could represent; `None` without the row.
    #[must_use]
    pub fn jet_noise_coverage(&self) -> Option<JetNoiseCoverage> {
        self.jet_noise.as_ref().map(LipJetNoise::coverage)
    }

    /// Attach the opt-in constriction noise before any step: a series pressure source added to the blowing pressure that
    /// drives the lip (its force) and the mouth reservoir the lip channel meters from, beside (not instead of) the jet noise,
    /// of rms `K(St) P_blow` in each third octave below the duct cutoff `100 / r`, `K` the printed trend of Oldham and Ukpoho
    /// and `P_blow` the gesture's blowing pressure (see [`ConstrictionNoise`] for the law, what it rests on and the label).
    /// Its bank runs at the gas step. It runs once per gas step on the previous gas step's metered flow, and is silent while
    /// that flow is not forward. A lip without it is bit-identical to one built before this existed; a card whose cutoff
    /// leaves no driven band adds exactly zero.
    ///
    /// # Errors
    /// The reason [`ConstrictionNoiseCard::validate`] gives.
    pub fn with_constriction_noise(mut self, card: ConstrictionNoiseCard) -> Result<Self, String> {
        self.constriction_noise = Some(ConstrictionNoise::new(card, self.clock.gas_step_s()).map_err(str::to_string)?);
        Ok(self)
    }

    /// Apply a control-rate change between blocks. A held bore has no valves.
    ///
    /// # Errors
    /// An unusable value, a refused lip rebuild, or a valve change.
    pub fn apply(&mut self, control: BrassControl) -> Result<(), String> {
        match control {
            BrassControl::SetBlowingPressure(p) => {
                if !(p.is_finite() && p >= 0.0) {
                    return Err("blowing pressure must be finite non-negative".into());
                }
                let pressure = self.room.pressure + p;
                let density = pressure / (self.gas.gas_constant() * self.room.temperature);
                let state = self
                    .gas
                    .conserved(density, 0.0, pressure)
                    .map_err(|e| format!("mouth reservoir refused: {e}"))?;
                let feed = self
                    .gas
                    .primitive(state, || false)
                    .map_err(|e| format!("mouth reservoir refused: {e}"))?;
                if let Some(supply) = self.finite_mouth.as_mut() {
                    supply.set_feed(feed);
                } else {
                    self.mouth = feed;
                }
                self.blow_pa = p;
            }
            BrassControl::SetLipTension(t) => {
                if !(t.is_finite() && t > 0.0) {
                    return Err("lip tension must be finite positive".into());
                }
                if self.lip2.is_some() && t != 1.0 {
                    return Err("the second lip mode (lip_response_mode2) has no printed pre-stress law and runs only at lip tension 1.0".into());
                }
                // As in the linear loop: stiffness ~ t^2, damping ~ t, same state.
                let lip_sys = mass_spring_damper(
                    self.lip.mass_kg,
                    self.lip.stiffness_n_m * t * t,
                    self.lip.damping_n_s_m * t,
                )
                .map_err(|_| "tensioned lip island refused admission".to_string())?;
                let lip_linear = prepare_lip(&lip_sys).map_err(|e| e.to_string())?;
                let work_j = spring_work_j(
                    &self.lip_sys,
                    &lip_sys,
                    &self.lip_state,
                    self.lip.stiffness_n_m * self.tension * self.tension,
                    self.lip.stiffness_n_m * t * t,
                )?;
                let tension_work = self.tension_work.with_increment([work_j, 0.0])?;
                self.lip_linear = lip_linear;
                self.lip_sys = lip_sys;
                self.tension = t;
                self.tension_work = tension_work;
            }
            BrassControl::SetValve { .. } => {
                return Err("the nonlinear bore is one held bore without valves".into());
            }
        }
        Ok(())
    }

    /// Per-block diagnostics so far.
    #[must_use]
    pub fn diagnostics(&self) -> &[BrassBlockDiag] {
        &self.diags
    }

    /// Signed tension-actuator work [J] by mode since construction, including
    /// accepted controls not yet rendered. The second entry is zero: the second
    /// lip mode runs only at tension 1.0, so no tension change reaches it. It is
    /// not part of gas exchange or a lip/gas passivity claim.
    #[must_use]
    pub fn lip_tension_work_j(&self) -> [f64; 2] {
        self.tension_work.total_j()
    }

    /// The bore's gas field and its receipts.
    #[must_use]
    pub const fn bore(&self) -> &NonlinearBore {
        &self.bore
    }

    /// The wall memory's admission receipt.
    #[must_use]
    pub fn wall_receipt(&self) -> Option<WallReceipt> {
        self.bore.wall_receipt()
    }

    /// The declared clock.
    #[must_use]
    pub const fn clock(&self) -> NonlinearClock {
        self.clock
    }

    /// Render one block: `mouthpiece` receives the cup-face overpressure [Pa]
    /// at each output sample, `flows` the mean bell volume flow [m^3/s] of
    /// every mechanics sample, `mechanics_per_output` per output sample.
    ///
    /// # Errors
    /// A length mismatch, or any refusal of the lip, the lip face or the gas
    /// bore. A refusal leaves the block partly advanced and must end the run.
    pub fn step_block(&mut self, mouthpiece: &mut [f64], flows: &mut [f64]) -> Result<(), String> {
        let ratio = self.clock.mechanics_per_output as usize;
        if mouthpiece.is_empty() || flows.len() != mouthpiece.len() * ratio {
            return Err("nonlinear lip block needs one flow per mechanics sample".into());
        }
        let steps = self.clock.gas_steps_per_mechanics as usize;
        let dt = self.clock.gas_step_s();
        let mut gap_min = f64::INFINITY;
        let mut gap_max = 0.0f64;
        let mut gap_acc = 0.0f64;
        let mut blow_work = 0.0f64;
        let mut bore_work = 0.0f64;
        let mut worst_defect = 0.0f64;
        for (sample, slot) in mouthpiece.iter_mut().enumerate() {
            for frame in 0..ratio {
                let mut flow = 0.0;
                for _ in 0..steps {
                    let (gap, defect, face, outlet, drive, mouth_density) = self.gas_step(dt)?;
                    gap_min = gap_min.min(gap);
                    gap_max = gap_max.max(gap);
                    gap_acc += gap;
                    worst_defect = worst_defect.max(defect);
                    let mouth_volume = face.mass_flow_kg_s.max(0.0) / mouth_density;
                    blow_work += drive * mouth_volume * dt;
                    bore_work += (face.pressure_pa - self.room.pressure) * face.mass_flow_kg_s
                        / self.room.density
                        * dt;
                    flow += outlet;
                }
                flows[sample * ratio + frame] = flow / steps as f64;
            }
            *slot = self.face_pa - self.room.pressure;
        }
        let output_dt = 1.0 / f64::from(self.clock.output_rate_hz);
        let n = (mouthpiece.len() * ratio * steps) as f64;
        let lip_tension_work_j = self.tension_work.finish_block();
        self.diags.push(BrassBlockDiag {
            block: self.block_index,
            lock_hz: zero_crossing_lock_hz(mouthpiece, output_dt),
            gap_min_m: gap_min,
            gap_mean_m: gap_acc / n,
            gap_max_m: gap_max,
            blow_work_j: blow_work,
            bore_work_j: bore_work,
            lip_tension_work_j,
            supply_defect: worst_defect,
        });
        self.block_index += 1;
        Ok(())
    }

    // One gas step: lip under the last cup pressure, then the bore with the
    // lip face at both stages. Returns the gap, the lip supply defect, the
    // first-stage face, the mean bell volume flow and the driving pressure the
    // step used, and the original mouth density for its volume-flow diagnostic.
    fn gas_step(&mut self, dt: f64) -> Result<(f64, f64, LipFace, f64, f64, f64), String> {
        // The opt-in jet noise is a series pressure source in the lip's driving pressure: it adds to the force on the lip and
        // to the stagnation pressure of the mouth reservoir the lip channel meters from, as in the linear loop. It reads the
        // previous gas step's metered flow and gap, and its speed is bounded by the Bernoulli speed of the present pressure
        // drop. The opt-in constriction noise joins the same series drive beside it: it reads the previous gas step's forward
        // flow, its amplitude follows the held blowing pressure, and it needs no gap. Zero is not added and the mouth state is
        // not rebuilt, so a lip without either, or whose gates are shut, keeps `drive == blow_pa` and `mouth` bit for bit.
        let mouth_density = self.mouth.density;
        let mut drive = if self.finite_mouth.is_some() {
            self.mouth.pressure - self.room.pressure
        } else {
            self.blow_pa
        };
        let mut mouth = self.mouth;
        let mut noisy = false;
        // The drive before any noise: the held blow pressure for the original reservoir, the stored gas's gauge pressure for a
        // finite mouth. The jet noise's Bernoulli bound reads the pressure drop of THIS drive.
        let base_gauge = drive;
        // The summed series noise of this gas step [Pa]; exactly zero when no source speaks.
        let mut series_pa = 0.0f64;
        if let Some(source) = self.jet_noise.as_mut() {
            let drop = (base_gauge - (self.face_pa - self.room.pressure)).abs();
            let bound = (2.0 * drop / self.reference_density_kg_m3).sqrt();
            let noise = source.pressure_pa(
                self.lip.width_m,
                bound,
                self.reference_density_kg_m3,
                self.wall_viscosity_pa_s,
                dt,
            );
            if noise != 0.0 {
                drive += noise;
                series_pa += noise;
                noisy = true;
            }
        }
        if let Some(source) = self.constriction_noise.as_mut() {
            let noise = source.pressure_pa(self.blow_pa, self.reference_density_kg_m3);
            if noise != 0.0 {
                drive += noise;
                series_pa += noise;
                noisy = true;
            }
        }
        // The original held reservoir is rebuilt at the noisy stagnation pressure and the room temperature. A finite mouth is
        // NOT rebuilt: its stored gas keeps its own state and the noise is handed to the finite-mouth step as a series offset
        // that both stages of the lip face see and the feed does not (see `FiniteMouthSpec` and `FiniteMouthTrial::flux`).
        if noisy && self.finite_mouth.is_none() {
            let pressure = self.room.pressure + drive;
            let density = pressure / (self.gas.gas_constant() * self.room.temperature);
            let state = self
                .gas
                .conserved(density, 0.0, pressure)
                .map_err(|e| format!("mouth reservoir refused: {e}"))?;
            mouth = self
                .gas
                .primitive(state, || false)
                .map_err(|e| format!("mouth reservoir refused: {e}"))?;
        }
        let force = (drive - (self.face_pa - self.room.pressure)) * self.lip.face_area_m2;
        let record = step_split(
            &self.lip_sys,
            &self.lip_linear,
            LinearPropagator::ExactZoh,
            &self.lip_state,
            &[force],
            dt,
        )
        .map_err(|_| "lip island step refused".to_string())?;
        let mut defect = record.supply_defect().abs();
        // The opt-in second lip mode steps under `sign x force` and its displacement joins the first's; without it `displacement`
        // is `record.x[0]` and every operation below is the one-mode lip's, in the same order.
        let mut displacement = record.x[0];
        let mut second_record = None;
        if let Some(second) = self.lip2.as_ref() {
            let record2 = step_split(
                &second.sys,
                &second.linear,
                LinearPropagator::ExactZoh,
                &second.state,
                &[second.mode.sign * force],
                dt,
            )
            .map_err(|_| "second lip mode step refused".to_string())?;
            defect = defect.max(record2.supply_defect().abs());
            displacement += record2.x[0];
            second_record = Some(record2);
        }
        let gap = (self.lip.rest_gap_m + displacement).max(0.0);
        let opening = self.flow_width_height(gap);
        let inlet = self.bore.inlet_area_m2();
        let gas = &self.gas;
        let mut finite_trial = self.finite_mouth.as_ref().map(|supply| supply.trial(mouth, dt, inlet, series_pa));
        let mut first = None;
        let step = self
            .bore
            .step(|stage: &PortStage| {
                let mut evaluate = |mouth| {
                    let face = lip_face(gas, stage.left.1, mouth, opening, inlet)?;
                    if stage.stage == 0 {
                        first = Some(face);
                    }
                    Ok(face.flux)
                };
                let flux = if let Some(trial) = finite_trial.as_mut() {
                    trial.flux(gas, stage, evaluate)?
                } else {
                    evaluate(mouth)?
                };
                Ok(Some(flux))
            })
            .map_err(|e| format!("nonlinear bore step refused: {e}"))?;
        let face = first.ok_or_else(|| "the bore offered no first stage".to_string())?;
        if let Some(trial) = finite_trial {
            let supply = self.finite_mouth.as_mut()
                .ok_or_else(|| "finite mouth lost its accepted owner".to_string())?;
            self.mouth = supply.commit(trial, &step.report)
                .map_err(|e| format!("finite mouth commit refused: {e}"))?;
        }
        self.lip_state = record.x;
        if let (Some(second), Some(record2)) = (self.lip2.as_mut(), second_record) {
            second.state = record2.x;
        }
        self.face_pa = face.pressure_pa;
        if let Some(source) = self.jet_noise.as_mut() {
            // The metered mass flow as a volume flow at the reference density: `Re_h = Q rho / (w mu) = mdot / (w mu)` needs
            // no density choice, and a negative flow (the cup above the mouth) is silent.
            source.record(face.mass_flow_kg_s / self.reference_density_kg_m3, gap);
        }
        if let Some(source) = self.constriction_noise.as_mut() {
            // The same metered volume flow: the source needs it forward (positive), and a negative flow is silent.
            source.record(face.mass_flow_kg_s / self.reference_density_kg_m3);
        }
        Ok((gap, defect, face, step.outlet_volume_flow_m3_s, drive, mouth_density))
    }
}

// Signed mass flow the lip channel meters at cup-face pressure `p` [kg/s]:
// mouth to cup below the mouth pressure, the cup's static gas back out above it.
fn metered(
    gas: &FrozenGas,
    star: Primitive,
    mouth: Primitive,
    opening: f64,
    p: f64,
) -> Result<f64, GasError> {
    if opening == 0.0 || p == mouth.pressure {
        Ok(0.0)
    } else if p < mouth.pressure {
        let (jet, _) = gas.nozzle_exit(mouth, p, || false)?;
        Ok(jet.density * jet.velocity * opening)
    } else {
        // The star state is the end wave's state AT `p`: the caller solved it for that pressure, so the donor's pressure is `p` by construction. The
        // primitive recovered from the conserved variables carries it only up to rounding. With `p` a few ulp above the mouth pressure (a lip just
        // opening at equilibrium) that recovery can come back below the mouth pressure, and the nozzle then refuses a stationary donor at equality
        // ("stationary nozzle donor or backpressure": NOZ#653 on the Bb4 card, p 1 ulp above the mouth, star 2 ulp below p). The nozzle's own
        // precondition (0 < backpressure <= donor.pressure, donor.velocity == 0) is unchanged.
        let donor = Primitive {
            velocity: 0.0,
            pressure: p,
            ..star
        };
        let (jet, _) = gas.nozzle_exit(donor, mouth.pressure, || false)?;
        Ok(-jet.density * jet.velocity * opening)
    }
}

fn lip_face(
    gas: &FrozenGas,
    interior: Primitive,
    mouth: Primitive,
    opening: f64,
    inlet_area: f64,
) -> Result<LipFace, GasError> {
    // The pressure at which the interior's end wave stops the gas: the face
    // pressure of a closed lip, and one end of the bracket for an open one.
    let mut wall = interior.pressure;
    let mut stopped = false;
    for _ in 0..64 {
        let wave = gas.end_wave(interior, -1.0, wall, || false)?;
        let step = wave.primitive.velocity / wave.velocity_per_pressure;
        wall -= step;
        if !(wall > 0.0 && wall.is_finite()) {
            return Err(GasError::Invalid(
                "lip face closed-end pressure left the gas domain",
            ));
        }
        if step.abs() <= 1e-14 * wall {
            stopped = true;
            break;
        }
    }
    if !stopped {
        return Err(GasError::Invalid(
            "lip face closed-end pressure did not converge",
        ));
    }
    // Continuity G(p) = metered(p) / A - face flux(p) falls strictly with p and
    // changes sign between the closed-end and mouth pressures.
    let residual = |p: f64| -> Result<(f64, Primitive, f64), GasError> {
        let wave = gas.end_wave(interior, -1.0, p, || false)?;
        let mass = metered(gas, wave.primitive, mouth, opening, p)?;
        let face = if wave.primitive.velocity >= 0.0 {
            gas.stagnation_to_static_state(mouth.temperature, p, wave.primitive.velocity, || false)?
        } else {
            wave.primitive
        };
        let value = mass / inlet_area - face.density * face.velocity;
        Ok((value, face, wave.signal_speed))
    };
    let mut pressure = wall;
    let (mut lo, mut hi) = (wall.min(mouth.pressure), wall.max(mouth.pressure));
    if opening > 0.0 && lo < hi {
        // Illinois false position on the bracket [lo, hi], G(lo) >= 0 >= G(hi).
        // At high mouth pressure the auxiliary upper face speed can exhaust
        // the caloric path even though the metered solution is supported.
        // Keep that pressure only as a bound and bisect until a supported
        // nonpositive residual exists; never accept a temperature clamp.
        let cold_trial = |error: &GasError| {
            matches!(
                error,
                GasError::CaloricPath("reservoir acceleration leaves cold source domain")
            )
        };
        let mut g_lo = residual(lo)?.0;
        let mut g_hi = match residual(hi) {
            Ok(value) => Some(value.0),
            Err(error) if wall < mouth.pressure && cold_trial(&error) => None,
            Err(error) => return Err(error),
        };
        let mut side = 0i8;
        let mut converged = false;
        for _ in 0..200 {
            if g_hi.is_none() && hi - lo <= 1e-12 * hi {
                return Err(GasError::CaloricPath(
                    "lip face continuity requires an unsupported cold inflow",
                ));
            }
            let next = match g_hi {
                Some(g_hi) if g_lo > g_hi => (lo * g_hi - hi * g_lo) / (g_hi - g_lo),
                _ => lo.midpoint(hi),
            };
            let next = if next > lo && next < hi {
                next
            } else {
                lo.midpoint(hi)
            };
            let g = match residual(next) {
                Ok(value) => value.0,
                Err(error) if g_hi.is_none() && cold_trial(&error) => {
                    hi = next;
                    side = 0;
                    continue;
                }
                Err(error) => return Err(error),
            };
            pressure = next;
            if g == 0.0 {
                converged = true;
                break;
            }
            if g > 0.0 {
                (lo, g_lo) = (next, g);
                if side == 1 {
                    g_hi = g_hi.map(|value| value * 0.5);
                }
                side = 1;
            } else {
                (hi, g_hi) = (next, Some(g));
                if side == -1 {
                    g_lo *= 0.5;
                }
                side = -1;
            }
            if g_hi.is_some() && hi - lo <= 1e-12 * hi {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(GasError::Invalid("lip face continuity did not converge"));
        }
    }
    let (_, face, signal) = residual(pressure)?;
    let flux = if face.velocity >= 0.0 {
        gas.inflow_flux(face, mouth.temperature, signal)?
    } else {
        gas.face_flux(face, signal)?
    };
    Ok(LipFace {
        flux,
        pressure_pa: pressure,
        mass_flow_kg_s: flux.mass * inlet_area,
    })
}

#[cfg(test)]
#[path = "nonlinear_tests.rs"]
mod tests;
