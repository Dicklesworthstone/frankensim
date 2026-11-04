//! Brass playing loop (music bead `frankensim-music-v8-root-3ez8g.4.3`):
//! the lip island x the multimodal bore, composed — no `Trumpet` type,
//! a composition recipe. The lip is a 1-DOF `fs_phs::mass_spring_damper`
//! island whose parameters come from a reduce-lab valve card (bead
//! 3ez8g.3.3), blowing through `fs_phs::bernoulli_volume_flow` into the
//! D18 mouthpiece-cup junction (bead 3ez8g.4.2) and the per-valve
//! multimodal characteristic lines terminated by the baked bell load
//! (beads 3ez8g.4.1 + zolja).
//!
//! PITCH IS NEVER ASSIGNED. The control surface is blowing pressure,
//! lip tension (a dimensionless pre-stress multiplier on the card's
//! stiffness), and the valve combination — no oscillator, no target
//! note, no frequency input of any kind (the battery greps this file
//! to prove it). The note that sounds is the LOCK the loop finds
//! between the lip resonance and the column's impedance peaks.
//!
//! Sign convention: the lip is OUTWARD-STRIKING — mouth overpressure
//! `p_blow - p_bore` pushes the lip open (`h = h0 + q`), the classic
//! brass embouchure family. The lip's supply defect is watched every
//! sample (the exact held-force quadratic island keeps an honest
//! energy ledger); per-block diagnostics carry the lock estimate, gap
//! statistics, and the power ledger.
//!
//! Honesty ladder: every fixture input is labeled card-backed or
//! authored in [`LipIsland::provenance`]; with authored tissue numbers
//! the whole loop is Estimate (claim-side promotion happens in the
//! brass-gates bead, not here). The lip MESH never appears in this
//! loop — only the reduced card's numbers do.

use fs_material::gas::GasState;
use fs_phs::linear_exact::{
    LinearPropagator, LinearSplit, LinearSplitBudget, split_quadratic, step_split,
};
use fs_phs::{PortHamiltonian, mass_spring_damper};

use crate::mm_line::{
    CupState, MmBankLoad, MmLineBank, MmLineConfig, MmLineError, cup_junction_boundary, cup_port,
};
use crate::turbulence_noise::{ConstrictionNoise, ConstrictionNoiseCard, JetNoiseCard, JetNoiseCoverage, LipJetNoise};
use fs_duct::Duct;

/// File-driven lip-valve performances observed at a baffled microphone.
pub mod performance;

/// The lip on the shared nonlinear gas bore (the `nonlinear-brass` feature).
#[cfg(feature = "nonlinear-brass")]
pub mod nonlinear;

/// The lip island's reduced parameters plus their provenance labels.
#[derive(Debug, Clone)]
pub struct LipIsland {
    /// Effective mass [kg] (reduce-lab card, face-unit normalization).
    pub mass_kg: f64,
    /// Effective stiffness [N/m] at unit tension.
    pub stiffness_n_m: f64,
    /// Effective damping [N s/m] at unit tension.
    pub damping_n_s_m: f64,
    /// Slit width [m] (card, measured from the mesh).
    pub width_m: f64,
    /// Rest gap [m] (card, measured from the mesh).
    pub rest_gap_m: f64,
    /// Pressure-collection face area [m^2] (card orifice face).
    pub face_area_m2: f64,
    /// The provenance sentence: WHICH inputs are card-backed and which
    /// are authored (the receipts-chain requirement).
    pub provenance: String,
}

/// The second, INWARD-striking lip mode (opt-in; see [`BrassVoice::with_lip_mode2`]).
///
/// Lip-reed measurements publish the mechanical response of an embouchure as two resonances: an outward-striking one and
/// an inward-striking one, the played note between them (Newton, Campbell and Gilbert, JASA 123(1) EL14-EL20, 2008, who
/// conclude that "a model incorporating at least two degrees of freedom should be required to simulate computationally
/// realistic oscillations"). Cullen (2000, PhD, University of Edinburgh, Appendix A Table A.2) prints both for twenty
/// embouchures of an artificial lip pair on a trombone-size mouthpiece: `fL1 QL1 1/muL1` and `fL2 QL2 1/muL2`, the
/// second `1/mu` NEGATIVE because that mode moves against the force. The aperture is the sum of the two modes' displacements,
/// each a mass-spring-damper `m = A / |1/mu|`, `k = m (2 pi f)^2`, `c = m (2 pi f) / Q` driven by the same transmural
/// force, the second with the sign of its printed `1/mu`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LipMode2 {
    /// Effective mass magnitude [kg]: `A / |1/mu_2|`.
    pub mass_kg: f64,
    /// Stiffness [N/m] at unit tension: `m (2 pi f_2)^2`.
    pub stiffness_n_m: f64,
    /// Damping [N s/m] at unit tension: `m (2 pi f_2) / Q_2`.
    pub damping_n_s_m: f64,
    /// `-1.0` for an inward-striking mode (a negative printed `1/mu`: the mode moves against the force), `1.0` otherwise.
    pub sign: f64,
}

impl LipMode2 {
    /// The mode from its published quantities and the lip's pressure-collecting face area `A` [m^2].
    ///
    /// # Errors
    /// The reason, when `f_hz`, `q` or `face_area_m2` is not positive finite, or `inverse_mu_m2_per_kg` is zero or not finite.
    pub fn from_response(
        f_hz: f64,
        q: f64,
        inverse_mu_m2_per_kg: f64,
        face_area_m2: f64,
    ) -> Result<LipMode2, &'static str> {
        if !(f_hz.is_finite() && f_hz > 0.0) {
            return Err("the second lip mode's resonance must be positive finite");
        }
        if !(q.is_finite() && q > 0.0) {
            return Err("the second lip mode's quality factor must be positive finite");
        }
        if !(inverse_mu_m2_per_kg.is_finite() && inverse_mu_m2_per_kg != 0.0) {
            return Err("the second lip mode's inverse mu must be finite and nonzero");
        }
        if !(face_area_m2.is_finite() && face_area_m2 > 0.0) {
            return Err("the lip's face area must be positive finite");
        }
        let mass = face_area_m2 / inverse_mu_m2_per_kg.abs();
        let omega = core::f64::consts::TAU * f_hz;
        Ok(LipMode2 {
            mass_kg: mass,
            stiffness_n_m: mass * omega * omega,
            damping_n_s_m: mass * omega / q,
            sign: if inverse_mu_m2_per_kg < 0.0 { -1.0 } else { 1.0 },
        })
    }
}

/// The second lip mode's integrator: its own exact quadratic island and state beside the lip island's.
struct LipSecondMode {
    mode: LipMode2,
    sys: PortHamiltonian,
    linear: LinearSplit,
    state: Vec<f64>,
}

/// Control-rate inputs (applied between blocks; D17).
#[derive(Debug, Clone, Copy)]
pub enum BrassControl {
    /// Steady blowing pressure [Pa].
    SetBlowingPressure(f64),
    /// Dimensionless lip pre-stress multiplier: stiffness scales by
    /// `tension^2` (lip resonance scales by `tension`), damping by
    /// `tension` (constant loss factor).
    SetLipTension(f64),
    /// Valve combination switch with a crossfade length [samples].
    SetValve {
        /// Target combo index in the bank.
        combo: usize,
        /// Crossfade samples (0 = hard swap).
        fade_samples: usize,
    },
}

/// Per-block diagnostics (the bead's rich-logging requirement).
#[derive(Debug, Clone)]
pub struct BrassBlockDiag {
    /// Block index.
    pub block: u64,
    /// Zero-crossing lock estimate over the block [Hz] (0 = silent).
    pub lock_hz: f64,
    /// Lip gap statistics over the block [m].
    pub gap_min_m: f64,
    /// Mean gap [m].
    pub gap_mean_m: f64,
    /// Max gap [m].
    pub gap_max_m: f64,
    /// Pneumatic work done by the mouth this block [J].
    pub blow_work_j: f64,
    /// Wave energy sent into the bore this block [J].
    pub bore_work_j: f64,
    /// Worst lip-island supply defect this block.
    pub supply_defect: f64,
}

/// The composed brass voice.
pub struct BrassVoice {
    bank: MmLineBank,
    lip_sys: PortHamiltonian,
    lip_linear: LinearSplit,
    lip_state: Vec<f64>,
    lip: LipIsland,
    tension: f64,
    blow_pa: f64,
    cup: CupState,
    cup_compliance: f64,
    p_bore_prev: f64,
    rho: f64,
    /// Dynamic viscosity of the ambient gas [Pa s] (only the opt-in jet noise's Reynolds gate reads it).
    mu: f64,
    dt: f64,
    block_index: u64,
    diags: Vec<BrassBlockDiag>,
    /// Opt-in lip jet noise ([`LipJetNoise`]: Stevens' form with Verge's constant, an extrapolation); `None` adds nothing.
    jet_noise: Option<LipJetNoise>,
    /// Opt-in constriction noise ([`ConstrictionNoise`]: Oldham and Ukpoho's printed duct-constriction law in third-octave
    /// bands, an extrapolation to the lip); `None` adds nothing.
    constriction_noise: Option<ConstrictionNoise>,
    /// Opt-in second lip mode ([`LipMode2`]); `None` leaves the loop exactly the one-mode loop.
    lip2: Option<LipSecondMode>,
    /// The opt-in flow-area exponent `q` ([`BrassVoice::with_lip_area_exponent`]); exactly 1.0 is the plain `w h`.
    area_exponent: f64,
}

impl BrassVoice {
    /// Compose the voice: valve-combination ducts + the lip card + a
    /// mouthpiece-cup volume (authored; D18 electrically-short).
    /// A shared mouth load or [`MmBankLoad::PerCombo`] is realized once, then
    /// valve changes carry the outgoing-wave history through the same bank.
    ///
    /// # Errors
    /// [`MmLineError`] from the bank realization or lip admission.
    pub fn new<'a>(
        combos: &[Duct],
        labels: &[&str],
        gas: &GasState,
        load: impl Into<MmBankLoad<'a>>,
        config: &MmLineConfig,
        lip: LipIsland,
        cup_volume_m3: f64,
    ) -> Result<BrassVoice, MmLineError> {
        if !(cup_volume_m3 > 0.0 && cup_volume_m3.is_finite()) {
            return Err(MmLineError::Invalid {
                what: "cup volume must be positive finite",
            });
        }
        if !(lip.mass_kg > 0.0
            && lip.stiffness_n_m > 0.0
            && lip.damping_n_s_m >= 0.0
            && lip.width_m > 0.0
            && lip.rest_gap_m > 0.0
            && lip.face_area_m2 > 0.0)
        {
            return Err(MmLineError::Invalid {
                what: "lip island parameters must be positive",
            });
        }
        let bank = MmLineBank::new(combos, labels, gas, load, config)?;
        let lip_sys = mass_spring_damper(lip.mass_kg, lip.stiffness_n_m, lip.damping_n_s_m)
            .map_err(|_| MmLineError::Invalid {
                what: "lip island refused admission",
            })?;
        let lip_linear = prepare_lip(&lip_sys)?;
        let cup_compliance = cup_volume_m3 / (gas.density * gas.sound_speed * gas.sound_speed);
        Ok(BrassVoice {
            bank,
            lip_sys,
            lip_linear,
            lip_state: vec![0.0, 0.0],
            lip,
            tension: 1.0,
            blow_pa: 0.0,
            cup: CupState::default(),
            cup_compliance,
            p_bore_prev: 0.0,
            rho: gas.density,
            mu: gas.dynamic_viscosity,
            dt: 1.0 / f64::from(config.sample_rate_hz),
            block_index: 0,
            diags: Vec::new(),
            jet_noise: None,
            constriction_noise: None,
            lip2: None,
            area_exponent: 1.0,
        })
    }

    /// Attach the opt-in lip jet noise before any sample is rendered: a series pressure source added to the lip's driving
    /// pressure (the force on the lip and the head of the flow solve alike), of rms `gain rho U_j^2 / 2` inside the band
    /// `St_h = f h / U_j` in `[st_lo, st_hi]` (a band reaching the loop's Nyquist frequency keeps only its part below it),
    /// shut while `Re_h = Q / (w nu)` is below the card's onset or the jet flow is not forward (see [`LipJetNoise`] for the
    /// law, what it rests on and the label). A voice without it is bit-identical to one built before this existed; a gate
    /// that never opens adds exactly zero.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] with the reason [`JetNoiseCard::validate`] gives.
    pub fn with_jet_noise(mut self, card: JetNoiseCard) -> Result<BrassVoice, MmLineError> {
        self.jet_noise = Some(LipJetNoise::new(card).map_err(|what| MmLineError::Invalid { what })?);
        Ok(self)
    }

    /// The lip jet noise's open-gate samples so far by how much of the band each could represent; `None` without the row.
    #[must_use]
    pub fn jet_noise_coverage(&self) -> Option<JetNoiseCoverage> {
        self.jet_noise.as_ref().map(LipJetNoise::coverage)
    }

    /// Add the opt-in second lip mode before any sample is rendered: a second mass-spring-damper (see [`LipMode2`] for what it
    /// is and where its numbers come from) stepped under the SAME transmural force with the sign of its printed `1/mu`,
    /// the aperture being `rest gap + q1 + q2` (clipped at zero). A voice without it is bit-identical to one built before
    /// this existed. Its published numbers carry no pre-stress law, so it runs only at lip tension 1.0 and
    /// [`BrassVoice::apply`] refuses any other `SetLipTension` by name.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] when the voice's tension is not 1.0 or the mode is refused admission.
    pub fn with_lip_mode2(mut self, mode: LipMode2) -> Result<BrassVoice, MmLineError> {
        if self.tension != 1.0 {
            return Err(MmLineError::Invalid {
                what: "the second lip mode (lip_response_mode2) has no printed pre-stress law and runs only at lip tension 1.0",
            });
        }
        let sys = mass_spring_damper(mode.mass_kg, mode.stiffness_n_m, mode.damping_n_s_m).map_err(|_| {
            MmLineError::Invalid {
                what: "the second lip mode refused admission",
            }
        })?;
        let linear = prepare_lip(&sys)?;
        self.lip2 = Some(LipSecondMode {
            mode,
            sys,
            linear,
            state: vec![0.0, 0.0],
        });
        Ok(self)
    }

    /// Add the opt-in flow-area exponent `q` before any sample is rendered: the aperture the Bernoulli flow sees becomes
    /// `w h (h/H0)^(q-1)` for the lip gap `h`, the rest gap `H0` and the width `w` (an area that grows as `h^q`; `q = 1` is the plain
    /// `w h` and is bit-identical to a voice without the row). The normalisation at `H0` is ours: it leaves the aperture at the
    /// rest gap unchanged. Copley and Strong (1996, JASA 99(2) 1219) plot the lip area against its height and describe an
    /// "almost linear relationship between area and height" as a first approximation without printing an exponent, so `q` is a
    /// bracket parameter and not a measured value. Only the flow area changes; the lip's motion is untouched.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] unless `1 <= q <= 4`.
    pub fn with_lip_area_exponent(mut self, q: f64) -> Result<BrassVoice, MmLineError> {
        if !(q.is_finite() && (1.0..=4.0).contains(&q)) {
            return Err(MmLineError::Invalid {
                what: "the lip area exponent (lip_area_exponent) must be finite and in [1, 4]",
            });
        }
        self.area_exponent = q;
        Ok(self)
    }

    /// The flow aperture `w h (h/H0)^(q-1)` [m^2] for the lip gap `h`; exactly `w h` when the exponent is 1 or the lip is closed.
    fn flow_width_height(&self, gap: f64) -> f64 {
        let plain = self.lip.width_m * gap;
        if self.area_exponent == 1.0 || gap <= 0.0 {
            plain
        } else {
            plain * (gap / self.lip.rest_gap_m).powf(self.area_exponent - 1.0)
        }
    }

    /// One lip step under `force` [N]: the first mode always, the second when present (driven by `sign x force`). Returns the
    /// summed displacement `q1 (+ q2)` [m] and the worst supply defect. Without the second mode this is exactly the original
    /// single-mode step (the same operations in the same order).
    fn step_lip(&mut self, force: f64, dt: f64) -> Result<(f64, f64), MmLineError> {
        let record = step_split(
            &self.lip_sys,
            &self.lip_linear,
            LinearPropagator::ExactZoh,
            &self.lip_state,
            &[force],
            dt,
        )
        .map_err(|_| MmLineError::Invalid {
            what: "lip island step refused",
        })?;
        let mut defect = record.supply_defect().abs();
        self.lip_state = record.x;
        let mut displacement = self.lip_state[0];
        if let Some(second) = self.lip2.as_mut() {
            let record2 = step_split(
                &second.sys,
                &second.linear,
                LinearPropagator::ExactZoh,
                &second.state,
                &[second.mode.sign * force],
                dt,
            )
            .map_err(|_| MmLineError::Invalid {
                what: "second lip mode step refused",
            })?;
            defect = defect.max(record2.supply_defect().abs());
            second.state = record2.x;
            displacement += second.state[0];
        }
        Ok((displacement, defect))
    }

    /// Attach the opt-in constriction noise before any sample is rendered: a series pressure source added to the lip's
    /// driving pressure beside (not instead of) the jet noise, of rms `K(St) P_blow` in each third octave below the duct
    /// cutoff, `K` the printed trend of Oldham and Ukpoho and `P_blow` the gesture's blowing pressure (see
    /// [`ConstrictionNoise`] for the law, what it rests on and the label). It is silent while the jet flow is not forward. A
    /// voice without it is bit-identical to one built before this existed.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] with the reason [`ConstrictionNoiseCard::validate`] gives.
    pub fn with_constriction_noise(mut self, card: ConstrictionNoiseCard) -> Result<BrassVoice, MmLineError> {
        self.constriction_noise =
            Some(ConstrictionNoise::new(card, self.dt).map_err(|what| MmLineError::Invalid { what })?);
        Ok(self)
    }

    /// Apply a control-rate delta (BETWEEN blocks; D17 lift for valve
    /// switches lives in the bank).
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] on an unusable value.
    pub fn apply(&mut self, control: BrassControl) -> Result<(), MmLineError> {
        match control {
            BrassControl::SetBlowingPressure(p) => {
                if !(p.is_finite() && p >= 0.0) {
                    return Err(MmLineError::Invalid {
                        what: "blowing pressure must be finite non-negative",
                    });
                }
                self.blow_pa = p;
            }
            BrassControl::SetLipTension(t) => {
                if !(t.is_finite() && t > 0.0) {
                    return Err(MmLineError::Invalid {
                        what: "lip tension must be finite positive",
                    });
                }
                if self.lip2.is_some() && t != 1.0 {
                    return Err(MmLineError::Invalid {
                        what: "the second lip mode (lip_response_mode2) has no printed pre-stress law and runs only at lip tension 1.0",
                    });
                }
                // Stiffness ~ t^2, damping ~ t: rebuild the island with
                // the SAME state (q, p persist through an embouchure
                // move — the lip does not teleport).
                let lip_sys = mass_spring_damper(
                    self.lip.mass_kg,
                    self.lip.stiffness_n_m * t * t,
                    self.lip.damping_n_s_m * t,
                )
                .map_err(|_| MmLineError::Invalid {
                    what: "tensioned lip island refused admission",
                })?;
                let lip_linear = prepare_lip(&lip_sys)?;
                self.lip_sys = lip_sys;
                self.lip_linear = lip_linear;
                self.tension = t;
            }
            BrassControl::SetValve {
                combo,
                fade_samples,
            } => {
                self.bank.switch(combo, fade_samples)?;
            }
        }
        Ok(())
    }

    /// Reset the whole voice to silence (lines, lip, cup) keeping the
    /// realized bank — the probe/battery seam.
    ///
    /// # Errors
    /// As [`MmLineBank::reset`].
    pub fn reset(&mut self) -> Result<(), MmLineError> {
        self.bank.reset()?;
        self.lip_state = vec![0.0, 0.0];
        if let Some(second) = self.lip2.as_mut() {
            second.state = vec![0.0, 0.0];
        }
        self.cup = CupState::default();
        self.p_bore_prev = 0.0;
        if let Some(source) = self.jet_noise.as_mut() {
            source.record(0.0, 0.0);
        }
        if let Some(source) = self.constriction_noise.as_mut() {
            source.record(0.0);
        }
        Ok(())
    }

    /// The bank's lift log (valve switches).
    #[must_use]
    pub fn lift_log(&self) -> &[crate::mm_line::MmLiftRecord] {
        self.bank.lift_log()
    }

    /// Per-block diagnostics so far.
    #[must_use]
    pub fn diagnostics(&self) -> &[BrassBlockDiag] {
        &self.diags
    }

    /// Realization receipts of the bank, one per valve combination.
    #[must_use]
    pub fn receipts(&self) -> &[crate::mm_line::MmComboReceipt] {
        self.bank.receipts()
    }

    /// Plane-mode characteristic impedance of the bore's input plane
    /// [Pa s/m^3]: the port impedance the outgoing wave is referenced to.
    #[must_use]
    pub fn port_impedance(&self) -> f64 {
        self.bank.zc0()
    }

    /// Render one block of mouthpiece pressure [Pa]. The loop each
    /// sample: lip island steps under the transmural force, the
    /// Bernoulli aperture meters the flow, the cup junction turns flow
    /// into waves, the bank returns the bore's reflection.
    ///
    /// # Errors
    /// [`MmLineError`] when any island refuses (non-finite state).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), MmLineError> {
        self.step_block_traced(out, None)
    }

    /// [`Self::step_block`] that also records, per sample, the wave `p+` the
    /// junction sent INTO the bore [Pa]. Together with the port impedance it
    /// is the bore's incident wave, which every downstream transfer of the
    /// same duct (mouth flow, radiation) is referenced to.
    ///
    /// # Errors
    /// As [`Self::step_block`]; also when `outgoing` and `out` differ in length.
    pub fn step_block_outgoing(
        &mut self,
        out: &mut [f64],
        outgoing: &mut [f64],
    ) -> Result<(), MmLineError> {
        if outgoing.len() != out.len() {
            return Err(MmLineError::Invalid {
                what: "outgoing-wave trace must match the block length",
            });
        }
        self.step_block_traced(out, Some(outgoing))
    }

    fn step_block_traced(
        &mut self,
        out: &mut [f64],
        mut outgoing: Option<&mut [f64]>,
    ) -> Result<(), MmLineError> {
        if out.is_empty() {
            return Err(MmLineError::Invalid {
                what: "empty block",
            });
        }
        let dt = self.dt;
        let mut gap_min = f64::INFINITY;
        let mut gap_max = 0.0f64;
        let mut gap_acc = 0.0f64;
        let mut blow_work = 0.0f64;
        let mut bore_work = 0.0f64;
        let mut worst_defect = 0.0f64;
        let zc0 = self.bank.zc0();
        for (index, slot) in out.iter_mut().enumerate() {
            // 1. Lip island under the transmural force (outward-striking:
            //    mouth overpressure opens the lip).
            // The opt-in jet noise is a series pressure source in the lip's driving pressure (force on the lip and head of
            // the flow solve alike, as the voice's glottal sources are). It reads the previous sample's jet flow and lip
            // gap, its speed is bounded by the Bernoulli speed of the present pressure drop, and zero is not added, so a
            // voice without it, or whose gate is shut, keeps `drive == blow_pa` bit for bit.
            let mut drive = self.blow_pa;
            if let Some(source) = self.jet_noise.as_mut() {
                let bound = (2.0 * (self.blow_pa - self.p_bore_prev).abs() / self.rho).sqrt();
                let noise = source.pressure_pa(self.lip.width_m, bound, self.rho, self.mu, dt);
                if noise != 0.0 {
                    drive += noise;
                }
            }
            // The opt-in constriction noise joins the same series drive: it reads the previous sample's jet flow (forward
            // flow only), its amplitude follows the held blowing pressure, and zero is not added, so a voice without it keeps
            // `drive` bit for bit.
            if let Some(source) = self.constriction_noise.as_mut() {
                let noise = source.pressure_pa(self.blow_pa, self.rho);
                if noise != 0.0 {
                    drive += noise;
                }
            }
            let dp_lip = drive - self.p_bore_prev;
            let force = dp_lip * self.lip.face_area_m2;
            let (displacement, defect) = self.step_lip(force, dt)?;
            worst_defect = worst_defect.max(defect);
            let gap = (self.lip.rest_gap_m + displacement).max(0.0);
            gap_min = gap_min.min(gap);
            gap_max = gap_max.max(gap);
            gap_acc += gap;
            // 2+3. IMPLICIT flow/junction solve at the CURRENT sample.
            //    An explicit Bernoulli flow computed from the previous
            //    sample's bore pressure carries a one-sample delay whose
            //    phase flips the aperture's small-signal conductance
            //    ANTI-dissipative above fs/4 — executed: a crook combo
            //    locked onto a 7.3 kHz parasitic scream. Substituting the
            //    (linear) cup-junction relation p = u/a + K into the
            //    Bernoulli law gives a closed-form quadratic in u, so the
            //    junction is solved exactly with no delay in the flow
            //    path (the lip GAP stays explicit: lip mass filters it).
            //    The bore's reflection joins the same solve: the junction,
            //    the flow and the reflection all hold at this sample.
            let boundary = self.bank.boundary_response()?;
            let (a, k_off) = cup_port(boundary, &self.cup, self.cup_compliance, zc0, dt);
            let wh = self.flow_width_height(gap);
            let flow = if wh <= 0.0 {
                0.0
            } else {
                let b = 2.0 * wh * wh / (self.rho * a);
                let head = drive - k_off;
                if head >= 0.0 {
                    let c = 2.0 * wh * wh * head / self.rho;
                    -0.5 * b + (0.25 * b * b + c).sqrt()
                } else {
                    let c = -2.0 * wh * wh * head / self.rho;
                    0.5 * b - (0.25 * b * b + c).sqrt()
                }
            };
            blow_work += drive * flow * dt;
            if let Some(source) = self.jet_noise.as_mut() {
                source.record(flow, gap);
            }
            if let Some(source) = self.constriction_noise.as_mut() {
                source.record(flow);
            }
            let (p_plus, p_minus, p_now) = cup_junction_boundary(
                flow,
                boundary,
                &mut self.cup,
                self.cup_compliance,
                zc0,
                dt,
            );
            let _ = self.bank.push(p_plus)?;
            if let Some(trace) = outgoing.as_deref_mut() {
                trace[index] = p_plus;
            }
            bore_work += (p_plus * p_plus - p_minus * p_minus) / zc0 * dt;
            self.p_bore_prev = p_now;
            *slot = p_now;
        }
        let n = out.len() as f64;
        self.diags.push(BrassBlockDiag {
            block: self.block_index,
            lock_hz: zero_crossing_lock_hz(out, dt),
            gap_min_m: gap_min,
            gap_mean_m: gap_acc / n,
            gap_max_m: gap_max,
            blow_work_j: blow_work,
            bore_work_j: bore_work,
            supply_defect: worst_defect,
        });
        self.block_index += 1;
        Ok(())
    }
}

// Zero-crossing lock estimate against the block mean (the mouthpiece
// pressure rides a DC operating point).
fn zero_crossing_lock_hz(out: &[f64], dt: f64) -> f64 {
    let n = out.len() as f64;
    let mean = out.iter().sum::<f64>() / n;
    let mut crossings_ac = 0usize;
    let mut prev = out[0] - mean;
    for &p in &out[1..] {
        let v = p - mean;
        if prev < 0.0 && v >= 0.0 {
            crossings_ac += 1;
        }
        prev = v;
    }
    crossings_ac as f64 / (n * dt)
}

fn prepare_lip(sys: &PortHamiltonian) -> Result<LinearSplit, MmLineError> {
    split_quadratic(
        sys,
        LinearSplitBudget {
            max_modes: 2,
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .map_err(|_| MmLineError::Invalid {
        what: "exact quadratic lip preparation refused",
    })
}

#[cfg(test)]
pub(crate) mod brass_loop_tests {
    use super::*;
    use crate::mm_line::MmLoad;
    use fs_duct::{LossModel, Segment, Termination};
    use fs_material::gas::GasSpec;

    pub(crate) fn air(temperature_k: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101_325.0).expect("air")
    }

    pub(crate) fn brass_combo(crook_m: f64) -> Duct {
        let mut segments = vec![Segment::Cylinder {
            radius: 0.006,
            length: 0.30,
        }];
        if crook_m > 0.0 {
            segments.push(Segment::Cylinder {
                radius: 0.006,
                length: crook_m,
            });
        }
        segments.push(Segment::Cone {
            inlet_radius: 0.006,
            outlet_radius: 0.012,
            length: 0.25,
        });
        segments.push(Segment::Cone {
            inlet_radius: 0.012,
            outlet_radius: 0.035,
            length: 0.12,
        });
        Duct { segments }
    }

    pub(crate) fn config() -> MmLineConfig {
        MmLineConfig {
            sample_rate_hz: 24_000,
            n_modes: 3,
            extra_slices: 1,
            // Exact Bessel (.76): the line ladder refuses a WideTube/Bessel kink.
            loss: LossModel::Bessel,
        }
    }

    /// AUTHORED lip island for the probe (card-backed lane wired in the
    /// battery proper).
    pub(crate) fn authored_lip() -> LipIsland {
        LipIsland {
            mass_kg: 1.8e-4,
            stiffness_n_m: 350.0,
            damping_n_s_m: 0.05,
            width_m: 0.012,
            rest_gap_m: 5.0e-4,
            face_area_m2: 2.0e-5,
            provenance: "AUTHORED probe values (Estimate)".to_string(),
        }
    }

    /// The probe voice (authored lip with a 1e-4 m^2 face, tension 1.6) with the opt-in lip jet noise, or without it.
    fn jet_voice(card: Option<JetNoiseCard>) -> BrassVoice {
        let gas = air(293.15);
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let voice = BrassVoice::new(&[brass_combo(0.0)], &["open"], &gas, &load, &config(), lip, 1.5e-6).expect("voice");
        let mut voice = match card {
            Some(card) => voice.with_jet_noise(card).expect("a usable card"),
            None => voice,
        };
        voice.apply(BrassControl::SetLipTension(1.6)).expect("tension");
        voice
    }

    /// Three 1200-sample blocks of mouthpiece pressure at a steady `blow` [Pa].
    fn jet_render(voice: &mut BrassVoice, blow: f64) -> Vec<f64> {
        voice.apply(BrassControl::SetBlowingPressure(blow)).expect("blow");
        let (mut out, mut block) = (Vec::new(), vec![0.0f64; 1200]);
        for _ in 0..3 {
            voice.step_block(&mut block).expect("block");
            out.extend_from_slice(&block);
        }
        out
    }

    fn jet_card(onset: f64, seed: u64) -> JetNoiseCard {
        JetNoiseCard { seed, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: onset }
    }

    fn same_bits(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
    }

    /// Gauntlet G0 and G3. A voice without the jet noise, and one whose gate can never open, are bit-identical; with the gate
    /// open (the lip at its rest gap passes `Re_h` 2300 at 3 kPa) the sound changes, is finite, and is seeded.
    #[test]
    fn lip_jet_noise_shut_gate_is_bit_identical_and_an_open_gate_changes_the_sound_seeded() {
        let plain = jet_render(&mut jet_voice(None), 3000.0);
        assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
        let shut = jet_render(&mut jet_voice(Some(jet_card(1.0e12, 7))), 3000.0);
        assert!(same_bits(&plain, &shut), "a gate that never opens adds exactly zero");
        let open = jet_render(&mut jet_voice(Some(jet_card(1700.0, 7))), 3000.0);
        assert!(open.iter().all(|p| p.is_finite()));
        assert!(!same_bits(&plain, &open), "an open gate must reach the pressure");
        assert!(same_bits(&open, &jet_render(&mut jet_voice(Some(jet_card(1700.0, 7))), 3000.0)), "a declared seed repeats bit for bit");
        assert!(!same_bits(&open, &jet_render(&mut jet_voice(Some(jet_card(1700.0, 8))), 3000.0)), "another seed gives another noise");
    }

    /// Gauntlet G3 with a planted negative. At 300 Pa the jet is 22 m/s, and `Re_h = Q / (w nu)` stays under 1700 for any lip
    /// gap this drive can open (it would take more than 1.1 mm): the row with the 1700 gate is bit-identical to no row, while the
    /// same row with NO gate (onset 0) speaks. The second half is what keeps the first from being vacuous. With no mouth
    /// pressure at all there is no forward flow and no source, whatever the gate.
    #[test]
    fn lip_jet_noise_is_gated_on_the_jet_and_has_no_source_without_forward_flow() {
        let plain = jet_render(&mut jet_voice(None), 300.0);
        let gated = jet_render(&mut jet_voice(Some(jet_card(1700.0, 7))), 300.0);
        assert!(same_bits(&plain, &gated), "below the gate's flow threshold the source stays shut");
        let ungated = jet_render(&mut jet_voice(Some(jet_card(0.0, 7))), 300.0);
        assert!(!same_bits(&plain, &ungated), "an ungated source speaks at the same pressure, so the gate is what shut the first");
        let still = jet_render(&mut jet_voice(None), 0.0);
        let none = jet_render(&mut jet_voice(Some(jet_card(0.0, 7))), 0.0);
        assert!(same_bits(&still, &none), "no mouth pressure, no forward flow: no source even with no gate");
    }

    /// The probe voice of [`jet_voice`] with the opt-in constriction noise (cup radius from the card), or without it.
    fn constriction_voice(card: Option<ConstrictionNoiseCard>) -> BrassVoice {
        let gas = air(293.15);
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let voice = BrassVoice::new(&[brass_combo(0.0)], &["open"], &gas, &load, &config(), lip, 1.5e-6).expect("voice");
        let mut voice = match card {
            Some(card) => voice.with_constriction_noise(card).expect("a usable card"),
            None => voice,
        };
        voice.apply(BrassControl::SetLipTension(1.6)).expect("tension");
        voice
    }

    fn constriction_card(seed: u64) -> ConstrictionNoiseCard {
        ConstrictionNoiseCard { seed, duct_radius_m: 8.4e-3 }
    }

    /// Gauntlet G0 and G3. A voice without the constriction noise is the voice built before it existed; with no mouth
    /// pressure there is no forward flow and no source, so the row adds exactly zero; with a blowing pressure the printed law
    /// drives bands (a negative control inside the test: the planted 12 m radius puts the cutoff `100 / r` under the bank's
    /// lowest band, and then the row is bit-identical to none even though the pressure is the same),
    /// the pressure changes, stays finite, repeats bit for bit under a declared seed, and another seed gives another noise.
    #[test]
    fn constriction_noise_is_absent_silent_without_a_jet_and_speaks_seeded_on_the_printed_law() {
        let plain = jet_render(&mut constriction_voice(None), 3000.0);
        assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
        let still = jet_render(&mut constriction_voice(None), 0.0);
        let none = jet_render(&mut constriction_voice(Some(constriction_card(7))), 0.0);
        assert!(same_bits(&still, &none), "no mouth pressure, no forward flow: no source");
        // Planted negative: a 12 m radius puts f_c = 8.3 Hz under the bank's lowest band (10 Hz), so nothing is driven at any
        // pressure. (A smaller "duct" does not do it: the Strouhal length grows with the radius, so a 3 m radius, f_c = 33 Hz,
        // still drives the 10-31.5 Hz bands, and the first run of this test said so.)
        let cut_off = jet_render(
            &mut constriction_voice(Some(ConstrictionNoiseCard { seed: 7, duct_radius_m: 12.0 })),
            3000.0,
        );
        assert!(same_bits(&plain, &cut_off), "no band below the duct cutoff, so the law supplies nothing");
        let open = jet_render(&mut constriction_voice(Some(constriction_card(7))), 3000.0);
        assert!(open.iter().all(|p| p.is_finite()));
        assert!(!same_bits(&plain, &open), "the printed law drives bands at this pressure, so it must reach the pressure");
        assert!(same_bits(&open, &jet_render(&mut constriction_voice(Some(constriction_card(7))), 3000.0)), "a declared seed repeats bit for bit");
        assert!(!same_bits(&open, &jet_render(&mut constriction_voice(Some(constriction_card(8))), 3000.0)), "another seed gives another noise");
    }

    /// G0/G3: distinct referenced mouth loads participate in the live lip/cup
    /// loop. Both hard and faded valve changes retain the actual outgoing waves.
    #[test]
    fn per_combo_strict_loads_render_and_carry_waves_through_valve_changes() {
        use fs_duct::{TabulatedLoad, TabulatedMouth};
        use fs_math::c64::C64;

        let gas = air(293.15);
        let cfg = MmLineConfig {
            sample_rate_hz: 48_000,
            n_modes: 2,
            extra_slices: 1,
            loss: LossModel::Lossless,
        };
        let radius = 0.002;
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let dt = 1.0 / f64::from(cfg.sample_rate_hz);
        let ceiling = crate::driving_point::LineLadder::default().max_taps;
        let resistances = [0.7, 1.3];
        let mut ducts = Vec::new();
        let mut tables = Vec::new();
        for (transit_samples, resistance) in [8.0, 12.0].into_iter().zip(resistances) {
            let length = transit_samples * gas.sound_speed * dt;
            ducts.push(Duct {
                segments: vec![Segment::Cylinder { radius, length }],
            });
            tables.push(
                TabulatedLoad::try_new(
                    vec![
                        core::f64::consts::TAU / (ceiling as f64 * dt),
                        core::f64::consts::PI / dt,
                    ],
                    vec![C64::new(resistance * zc, 0.0); 2],
                    "analytic resistive valve-path control",
                    0.0,
                )
                .unwrap()
                .with_mouth(TabulatedMouth {
                    axial_position_m: length,
                    radius_m: radius,
                    gas,
                    dc_impedance: resistance * zc,
                })
                .unwrap(),
            );
        }
        let loads: Vec<_> = tables.iter().map(MmLoad::Tabulated).collect();
        let mut voice = BrassVoice::new(
            &ducts,
            &["short", "long"],
            &gas,
            MmBankLoad::PerCombo(&loads),
            &cfg,
            authored_lip(),
            1.5e-6,
        )
        .expect("each valve path has its own strict load");
        for (combo, resistance) in resistances.into_iter().enumerate() {
            let receipt = &voice.receipts()[combo];
            assert_eq!(receipt.table_bins, receipt.n_fft / 2);
            assert_eq!(receipt.fallback_bins, 0);
            let dc: f64 = voice.bank.combo_taps(combo).iter().sum();
            let expected_dc = (resistance - 1.0) / (resistance + 1.0);
            assert!((dc - expected_dc).abs() < 2.0e-11);
        }
        assert_eq!(voice.receipts()[0].n_fft, voice.receipts()[1].n_fft);
        voice
            .apply(BrassControl::SetBlowingPressure(1_500.0))
            .unwrap();
        let mut history = Vec::new();
        let mut previous_combo = 0;
        for (combo, fade_samples, samples) in [(0, 0, 256), (1, 16, 64), (0, 0, 64)] {
            let lip_before = voice.lip_state.clone();
            let cup_before = voice.cup;
            voice
                .apply(BrassControl::SetValve {
                    combo,
                    fade_samples,
                })
                .unwrap();
            assert_eq!(voice.lip_state, lip_before);
            assert_eq!(voice.cup.p_prev, cup_before.p_prev);
            assert_eq!(voice.cup.u_prev, cup_before.u_prev);
            for sample in 0..samples {
                let mut pressure = [0.0];
                let mut outgoing = [0.0];
                voice
                    .step_block_outgoing(&mut pressure, &mut outgoing)
                    .unwrap();
                assert!(pressure[0].is_finite() && outgoing[0].is_finite());
                history.push(outgoing[0]);
                // Independent direct convolution of the rendered source trace,
                // including waves emitted before this valve was selected.
                let reflected = |which| -> f64 {
                    voice
                        .bank
                        .combo_taps(which)
                        .iter()
                        .zip(history.iter().rev())
                        .map(|(tap, wave)| tap * wave)
                        .sum()
                };
                let expected = if sample < fade_samples {
                    let weight = (sample + 1) as f64 / fade_samples as f64;
                    (1.0 - weight) * reflected(previous_combo) + weight * reflected(combo)
                } else {
                    reflected(combo)
                };
                let actual = voice.bank.incoming();
                assert!(
                    (actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0),
                    "combo {combo}, sample {sample}: incoming {actual} != history {expected}"
                );
            }
            previous_combo = combo;
        }
        assert!(history.iter().any(|wave| wave.abs() > 1.0));
        assert_eq!(voice.lift_log().len(), 2);
        for lift in voice.lift_log() {
            assert_eq!(lift.carried, voice.receipts()[0].n_fft);
        }
    }

    #[test]
    fn exact_lip_step_and_tension_rebuild_match_the_closed_form() {
        let gas = air(293.15);
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let lip = authored_lip();
        let (mass, stiffness, loss) = (lip.mass_kg, lip.stiffness_n_m, lip.damping_n_s_m);
        let mut voice = BrassVoice::new(
            &[brass_combo(0.0)],
            &["analytic fixture"],
            &gas,
            &load,
            &cfg,
            lip,
            1.5e-6,
        )
        .unwrap();
        voice.lip_state = vec![1e-5, 0.0];
        let initial = voice.lip_state.clone();
        voice.apply(BrassControl::SetLipTension(3.0)).unwrap();
        assert_eq!(
            voice.lip_state, initial,
            "repreparation must not reset momentum or position"
        );
        let alpha = 3.0 * loss / (2.0 * mass);
        let omega2 = 9.0 * stiffness / mass;
        let damped = fs_math::det::sqrt(omega2 - alpha * alpha);
        let decay = fs_math::det::exp(-alpha * voice.dt);
        let sine = fs_math::det::sin(damped * voice.dt);
        let cosine = fs_math::det::cos(damped * voice.dt);
        let expected_q = initial[0] * decay * (cosine + alpha / damped * sine);
        let expected_p = -mass * initial[0] * decay * omega2 / damped * sine;
        voice.step_block(&mut [0.0]).unwrap();
        assert!((voice.lip_state[0] - expected_q).abs() <= 1e-12 * initial[0]);
        assert!(
            (voice.lip_state[1] - expected_p).abs()
                <= 1e-12 * (mass * fs_math::det::sqrt(omega2) * initial[0])
        );
    }

    /// Dominant spectral peak of a block [Hz] (parabolic-refined),
    /// searched above 40 Hz.
    pub(crate) fn dominant_hz(block: &[f64], rate: f64) -> f64 {
        use fs_fft::{C64 as FftC64, Fft};
        let n = block.len().next_power_of_two();
        let mean = block.iter().sum::<f64>() / block.len() as f64;
        let mut buf: Vec<FftC64> = (0..n)
            .map(|k| FftC64::new(block.get(k).map_or(0.0, |p| p - mean), 0.0))
            .collect();
        let mut scratch = vec![FftC64::new(0.0, 0.0); n];
        Fft::new(n).forward(&mut buf, &mut scratch);
        let df = rate / n as f64;
        let k_lo = (40.0 / df).ceil() as usize;
        let mags: Vec<f64> = buf[..n / 2]
            .iter()
            .map(|c| (c.re * c.re + c.im * c.im).sqrt())
            .collect();
        let mut best = k_lo;
        for k in k_lo..mags.len() {
            if mags[k] > mags[best] {
                best = k;
            }
        }
        if best == 0 || best + 1 >= mags.len() || mags[best] <= 0.0 {
            return 0.0;
        }
        let (ya, yb, yc) = (
            mags[best - 1].max(1e-300).ln(),
            mags[best].ln(),
            mags[best + 1].max(1e-300).ln(),
        );
        let denom = ya - 2.0 * yb + yc;
        let shift = if denom.abs() > 1e-12 {
            0.5 * (ya - yc) / denom
        } else {
            0.0
        };
        (best as f64 + shift) * df
    }

    #[test]
    #[ignore = "operating-point probe: bore-vs-lip lock discrimination"]
    fn zz_probe_bore_discrimination() {
        let gas = air(293.15);
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let combos = [brass_combo(0.0), brass_combo(0.08), brass_combo(0.15)];
        let mut voice = BrassVoice::new(
            &combos,
            &["open", "c1", "c2"],
            &gas,
            &load,
            &cfg,
            authored_lip(),
            1.5e-6,
        )
        .expect("voice");
        for combo in 0..3usize {
            for &(tension, blow) in &[
                (0.8f64, 4000.0f64),
                (0.8, 6000.0),
                (1.0, 4000.0),
                (1.3, 6000.0),
                (1.6, 8000.0),
                (2.0, 8000.0),
            ] {
                voice.reset().expect("reset");
                voice
                    .apply(BrassControl::SetValve {
                        combo,
                        fade_samples: 0,
                    })
                    .expect("valve");
                voice
                    .apply(BrassControl::SetLipTension(tension))
                    .expect("t");
                voice
                    .apply(BrassControl::SetBlowingPressure(blow))
                    .expect("p");
                let mut block = vec![0.0f64; 2400];
                for _ in 0..8 {
                    voice.step_block(&mut block).expect("block");
                }
                let d = voice.diagnostics().last().expect("diag").clone();
                let mean = block.iter().sum::<f64>() / block.len() as f64;
                let amp = block
                    .iter()
                    .map(|p| (p - mean).abs())
                    .fold(0.0f64, f64::max);
                let f_dom = dominant_hz(&block, 24_000.0);
                println!(
                    "combo={combo} t={tension:.1} blow={blow:.0} zc={:.1} fft={f_dom:.1} amp={amp:.0}",
                    d.lock_hz
                );
            }
        }
    }

    #[test]
    #[ignore = "operating-point probe: upper-partial reachability"]
    fn zz_probe_upper_partials() {
        let gas = air(293.15);
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let combos = [brass_combo(0.0)];
        for &face in &[5.0e-5f64, 1.0e-4] {
            let mut lip = authored_lip();
            lip.face_area_m2 = face;
            let mut voice =
                BrassVoice::new(&combos, &["open"], &gas, &load, &cfg, lip, 1.5e-6).expect("voice");
            for &tension in &[0.8f64, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2] {
                for &blow in &[3000.0f64, 6000.0] {
                    voice.reset().expect("reset");
                    voice
                        .apply(BrassControl::SetLipTension(tension))
                        .expect("t");
                    voice
                        .apply(BrassControl::SetBlowingPressure(blow))
                        .expect("p");
                    let mut block = vec![0.0f64; 2400];
                    for _ in 0..8 {
                        voice.step_block(&mut block).expect("block");
                    }
                    let mean = block.iter().sum::<f64>() / block.len() as f64;
                    let amp = block
                        .iter()
                        .map(|p| (p - mean).abs())
                        .fold(0.0f64, f64::max);
                    let f_dom = dominant_hz(&block, 24_000.0);
                    let f_lip =
                        (350.0 * tension * tension / 1.8e-4).sqrt() / core::f64::consts::TAU;
                    println!(
                        "face={face:.0e} t={tension:.1} blow={blow:.0} f_lip={f_lip:.0}                          fft={f_dom:.1} amp={amp:.0}"
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "operating-point probe: valve ladder only"]
    fn zz_probe_valve_only() {
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let gas = air(293.15);
        let combos = [brass_combo(0.0), brass_combo(0.08)];
        let mut voice = BrassVoice::new(&combos, &["open", "c1"], &gas, &load, &cfg, lip, 1.5e-6)
            .expect("voice");
        for combo in 0..2usize {
            voice.reset().expect("reset");
            voice
                .apply(BrassControl::SetValve {
                    combo,
                    fade_samples: 0,
                })
                .expect("valve");
            voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
            voice
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut block = vec![0.0f64; 2400];
            for _ in 0..8 {
                voice.step_block(&mut block).expect("block");
            }
            let mean = block.iter().sum::<f64>() / block.len() as f64;
            let amp = block
                .iter()
                .map(|p| (p - mean).abs())
                .fold(0.0f64, f64::max);
            let f = dominant_hz(&block, 24_000.0);
            println!("VALVE combo={combo} lock={f:.1} amp={amp:.0}");
        }
    }

    #[test]
    #[ignore = "operating-point probe: valve ladder, slot walk, temperature"]
    #[allow(clippy::too_many_lines)] // one coherent probe
    fn zz_probe_gates_measurements() {
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        // Valve ladder + slot walk at 293 K.
        let gas = air(293.15);
        let combos = [brass_combo(0.0), brass_combo(0.08), brass_combo(0.15)];
        let mut voice = BrassVoice::new(
            &combos,
            &["open", "c1", "c2"],
            &gas,
            &load,
            &cfg,
            lip.clone(),
            1.5e-6,
        )
        .expect("voice");
        let run = |voice: &mut BrassVoice, combo: usize, tension: f64| -> (f64, f64) {
            voice.reset().expect("reset");
            voice
                .apply(BrassControl::SetValve {
                    combo,
                    fade_samples: 0,
                })
                .expect("valve");
            voice
                .apply(BrassControl::SetLipTension(tension))
                .expect("t");
            voice
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut block = vec![0.0f64; 2400];
            for _ in 0..8 {
                voice.step_block(&mut block).expect("block");
            }
            let mean = block.iter().sum::<f64>() / block.len() as f64;
            let amp = block
                .iter()
                .map(|p| (p - mean).abs())
                .fold(0.0f64, f64::max);
            (dominant_hz(&block, 24_000.0), amp)
        };
        for &blow in &[3000.0f64] {
            for combo in 0..3usize {
                voice.reset().expect("reset");
                voice
                    .apply(BrassControl::SetValve {
                        combo,
                        fade_samples: 0,
                    })
                    .expect("valve");
                voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
                voice
                    .apply(BrassControl::SetBlowingPressure(blow))
                    .expect("p");
                let mut block = vec![0.0f64; 2400];
                for _ in 0..8 {
                    voice.step_block(&mut block).expect("block");
                }
                let mean = block.iter().sum::<f64>() / block.len() as f64;
                let amp = block
                    .iter()
                    .map(|p| (p - mean).abs())
                    .fold(0.0f64, f64::max);
                let f = dominant_hz(&block, 24_000.0);
                println!("VALVE blow={blow:.0} combo={combo} lock={f:.1} amp={amp:.0}");
            }
        }
        for &t in &[1.5f64, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1] {
            let (f, amp) = run(&mut voice, 0, t);
            println!("SLOT t={t:.2} lock={f:.1} amp={amp:.0}");
        }
        // Temperature arm: cold horn.
        let cold = air(277.15);
        let mut voice_cold = BrassVoice::new(
            &[brass_combo(0.0)],
            &["open"],
            &cold,
            &load,
            &cfg,
            lip,
            1.5e-6,
        )
        .expect("voice");
        let (f_cold, amp_cold) = {
            voice_cold
                .apply(BrassControl::SetLipTension(0.8))
                .expect("t");
            voice_cold
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut block = vec![0.0f64; 2400];
            for _ in 0..8 {
                voice_cold.step_block(&mut block).expect("block");
            }
            let mean = block.iter().sum::<f64>() / block.len() as f64;
            let amp = block
                .iter()
                .map(|p| (p - mean).abs())
                .fold(0.0f64, f64::max);
            (dominant_hz(&block, 24_000.0), amp)
        };
        println!("COLD lock={f_cold:.1} amp={amp_cold:.0}");
    }

    #[test]
    #[ignore = "operating-point probe: prints the lock map"]
    fn zz_probe_lock_map() {
        let gas = air(293.15);
        let cfg = config();
        let combos = [brass_combo(0.0)];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        for &tension in &[0.6f64, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0] {
            for &blow in &[1000.0f64, 2000.0, 4000.0] {
                let mut voice = BrassVoice::new(
                    &combos,
                    &["open"],
                    &gas,
                    &load,
                    &cfg,
                    authored_lip(),
                    1.5e-6,
                )
                .expect("voice");
                voice
                    .apply(BrassControl::SetLipTension(tension))
                    .expect("t");
                voice
                    .apply(BrassControl::SetBlowingPressure(blow))
                    .expect("p");
                let mut block = vec![0.0f64; 2400];
                for _ in 0..5 {
                    voice.step_block(&mut block).expect("block");
                }
                let d = voice.diagnostics().last().expect("diag").clone();
                let f_lip = (350.0 * tension * tension / 1.8e-4).sqrt() / core::f64::consts::TAU;
                let amp = block
                    .iter()
                    .map(|p| (p - block.iter().sum::<f64>() / block.len() as f64).abs())
                    .fold(0.0f64, f64::max);
                println!(
                    "t={tension:.1} blow={blow:.0} f_lip={f_lip:.0} lock={:.1} amp={amp:.1} \
                     gap=[{:.2e},{:.2e}] defect={:.1e}",
                    d.lock_hz, d.gap_min_m, d.gap_max_m, d.supply_defect
                );
            }
        }
    }

    /// The probe voice (authored lip with a 1e-4 m^2 face, tension 1) with the opt-in second lip mode, or without it.
    fn second_mode_voice(second: Option<LipMode2>) -> BrassVoice {
        let gas = air(293.15);
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let voice = BrassVoice::new(&[brass_combo(0.0)], &["open"], &gas, &load, &config(), lip, 1.5e-6).expect("voice");
        match second {
            Some(mode) => voice.with_lip_mode2(mode).expect("second mode"),
            None => voice,
        }
    }

    /// Cullen (2000) Table A.2 embouchure 7's inward-striking mode: 244.5 Hz, Q 11.0, 1/mu -0.073 m^2/kg.
    fn embouchure_7_mode2() -> LipMode2 {
        LipMode2::from_response(244.5, 11.0, -0.073, 1.0e-4).expect("the printed numbers")
    }

    /// Displacement of a mass-spring-damper released from rest under a constant force, underdamped closed form.
    fn damped_step_response(force: f64, m: f64, k: f64, c: f64, t: f64) -> f64 {
        let alpha = c / (2.0 * m);
        let wd = (k / m - alpha * alpha).sqrt();
        force / k * (1.0 - (-alpha * t).exp() * ((wd * t).cos() + alpha / wd * (wd * t).sin()))
    }

    #[test]
    fn second_lip_mode_follows_the_published_quantities_and_refuses_unusable_ones() {
        let mode = embouchure_7_mode2();
        let omega = core::f64::consts::TAU * 244.5;
        assert!((mode.mass_kg - 1.0e-4 / 0.073).abs() <= 1e-15 * mode.mass_kg);
        assert!((mode.stiffness_n_m - mode.mass_kg * omega * omega).abs() <= 1e-12 * mode.stiffness_n_m);
        assert!((mode.damping_n_s_m - mode.mass_kg * omega / 11.0).abs() <= 1e-12 * mode.damping_n_s_m);
        assert_eq!(mode.sign, -1.0, "a negative printed 1/mu is an inward-striking mode");
        assert_eq!(LipMode2::from_response(244.5, 11.0, 0.073, 1.0e-4).expect("outward").sign, 1.0);
        for (f, q, inverse_mu, area, needle) in [
            (0.0, 11.0, -0.073, 1.0e-4, "resonance"),
            (f64::NAN, 11.0, -0.073, 1.0e-4, "resonance"),
            (244.5, 0.0, -0.073, 1.0e-4, "quality factor"),
            (244.5, 11.0, 0.0, 1.0e-4, "nonzero"),
            (244.5, 11.0, f64::INFINITY, 1.0e-4, "nonzero"),
            (244.5, 11.0, -0.073, 0.0, "face area"),
        ] {
            let error = LipMode2::from_response(f, q, inverse_mu, area).expect_err(needle);
            assert!(error.contains(needle), "{error}");
        }
    }

    #[test]
    fn second_lip_mode_static_response_is_the_signed_sum_of_two_compliances() {
        let mode = embouchure_7_mode2();
        let mut voice = second_mode_voice(Some(mode));
        let (dt, force) = (1.0 / 24_000.0, 2.5e-2);
        let mut displacement = 0.0f64;
        for _ in 0..48_000 {
            displacement = voice.step_lip(force, dt).expect("step").0;
        }
        let expected = force / voice.lip.stiffness_n_m + mode.sign * force / mode.stiffness_n_m;
        assert!(
            (displacement - expected).abs() <= 1e-9 * expected.abs(),
            "static opening {displacement:.12e} against F (1/k1 - 1/k2) = {expected:.12e}"
        );
        assert!(displacement < force / voice.lip.stiffness_n_m, "the inward-striking mode moves against the force");
    }

    #[test]
    fn second_lip_mode_step_response_matches_the_two_oscillator_closed_form() {
        let mode = embouchure_7_mode2();
        let mut voice = second_mode_voice(Some(mode));
        let (dt, force) = (1.0 / 24_000.0, 2.5e-2);
        let (m1, k1, c1) = (voice.lip.mass_kg, voice.lip.stiffness_n_m, voice.lip.damping_n_s_m);
        let mut at = [0.0f64; 2];
        for n in 1..=120usize {
            let (displacement, _) = voice.step_lip(force, dt).expect("step");
            if n == 7 {
                at[0] = displacement;
            }
            if n == 120 {
                at[1] = displacement;
            }
        }
        for (n, got) in [(7usize, at[0]), (120usize, at[1])] {
            let t = n as f64 * dt;
            let expected = damped_step_response(force, m1, k1, c1, t)
                + mode.sign * damped_step_response(force, mode.mass_kg, mode.stiffness_n_m, mode.damping_n_s_m, t);
            assert!(
                (got - expected).abs() <= 1e-9 * expected.abs().max(1e-12),
                "step {n}: {got:.12e} against the closed form {expected:.12e}"
            );
        }
    }

    #[test]
    fn second_lip_mode_runs_only_at_unit_tension_and_names_the_refusal() {
        let mut voice = second_mode_voice(Some(embouchure_7_mode2()));
        voice.apply(BrassControl::SetLipTension(1.0)).expect("tension 1 is the printed state");
        let refusal = voice.apply(BrassControl::SetLipTension(1.2)).expect_err("another tension has no printed law");
        assert!(format!("{refusal:?}").contains("lip_response_mode2"), "{refusal:?}");
        let mut one_mode = second_mode_voice(None);
        one_mode.apply(BrassControl::SetLipTension(1.2)).expect("the one-mode lip scales as before");
        let mut tensioned = second_mode_voice(None);
        tensioned.apply(BrassControl::SetLipTension(1.3)).expect("tension");
        let refusal = tensioned.with_lip_mode2(embouchure_7_mode2()).err().expect("added after a tension change");
        assert!(format!("{refusal:?}").contains("lip_response_mode2"), "{refusal:?}");
    }

    #[test]
    fn a_voice_without_the_second_mode_steps_exactly_the_original_lip_island() {
        let mut voice = second_mode_voice(None);
        let dt = 1.0 / 24_000.0;
        let reference = step_split(&voice.lip_sys, &voice.lip_linear, LinearPropagator::ExactZoh, &voice.lip_state, &[1.0e-2], dt)
            .expect("reference step");
        let (displacement, defect) = voice.step_lip(1.0e-2, dt).expect("step");
        assert!(voice.lip2.is_none());
        assert_eq!(voice.lip_state, reference.x);
        assert_eq!(displacement.to_bits(), reference.x[0].to_bits());
        assert_eq!(defect.to_bits(), reference.supply_defect().abs().to_bits());
    }

    #[test]
    fn second_lip_mode_changes_the_rendered_sound_and_reset_clears_it() {
        let one = jet_render(&mut second_mode_voice(None), 4000.0);
        let mut voice = second_mode_voice(Some(embouchure_7_mode2()));
        let two = jet_render(&mut voice, 4000.0);
        assert!(two.iter().all(|x| x.is_finite()));
        assert!(!same_bits(&one, &two), "a second mode must reach the mouthpiece pressure");
        assert!(voice.lip2.as_ref().is_some_and(|second| second.state.iter().any(|x| *x != 0.0)));
        voice.reset().expect("reset");
        assert_eq!(voice.lip_state, vec![0.0, 0.0]);
        assert_eq!(voice.lip2.as_ref().expect("second mode").state, vec![0.0, 0.0]);
    }

    #[test]
    fn lip_area_exponent_one_is_bit_identical_and_q_scales_the_aperture_about_the_rest_gap() {
        // q = 1 (the row set to exactly 1) steps exactly the original loop, bit for bit.
        let plain = jet_render(&mut second_mode_voice(None), 4000.0);
        let unit = jet_render(&mut second_mode_voice(None).with_lip_area_exponent(1.0).expect("q = 1"), 4000.0);
        assert!(same_bits(&plain, &unit), "q = 1 must be the plain w h");
        // The aperture: unchanged at the rest gap, multiplied by (h/H0)^(q-1) elsewhere, zero when closed.
        let mut voice = second_mode_voice(None).with_lip_area_exponent(2.0).expect("q = 2");
        let (w, h0) = (voice.lip.width_m, voice.lip.rest_gap_m);
        assert_eq!(voice.flow_width_height(h0).to_bits(), (w * h0).to_bits());
        let wide = 2.0 * h0;
        assert!((voice.flow_width_height(wide) - w * wide * 2.0).abs() <= 1e-15 * w * wide * 2.0);
        let narrow = 0.5 * h0;
        assert!((voice.flow_width_height(narrow) - w * narrow * 0.5).abs() <= 1e-15 * w * narrow);
        assert_eq!(voice.flow_width_height(0.0), 0.0);
        // The exponent reaches the mouthpiece pressure.
        let steeper = jet_render(&mut voice, 4000.0);
        assert!(steeper.iter().all(|x| x.is_finite()));
        assert!(!same_bits(&plain, &steeper), "q = 2 must reach the mouthpiece pressure");
    }

    #[test]
    fn lip_area_exponent_outside_one_to_four_is_refused_by_name() {
        for q in [0.0, 0.99, 4.01, -1.0, f64::NAN, f64::INFINITY] {
            let refusal = second_mode_voice(None).with_lip_area_exponent(q).err().expect("outside [1, 4]");
            assert!(format!("{refusal:?}").contains("lip_area_exponent"), "{q}: {refusal:?}");
        }
        assert!(second_mode_voice(None).with_lip_area_exponent(4.0).is_ok());
    }
}

#[cfg(test)]
mod brass_gate_tests {
    use super::brass_loop_tests::{air, authored_lip, brass_combo, config, dominant_hz};
    use super::*;
    use crate::mm_line::MmLoad;
    use fs_duct::Termination;
    use fs_math::c64::C64;

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    fn cents(a: f64, b: f64) -> f64 {
        1200.0 * (b / a).log2()
    }

    /// Impedance peaks of the CUP-LOADED column `|Z/(1 + i w C Z)|`
    /// from the MM authority (the column the lips actually see).
    fn cup_loaded_peaks_hz(
        duct: &Duct,
        gas: &GasState,
        cfg: &MmLineConfig,
        c_cup: f64,
    ) -> Vec<f64> {
        let mut mags = Vec::new();
        for i in 0..391 {
            let f = 120.0 + (900.0 - 120.0) * f64::from(i) / 390.0;
            let omega = core::f64::consts::TAU * f;
            let z = fs_duct::modal::mm_input_impedance(
                duct,
                gas,
                omega,
                cfg.loss,
                Termination::FlangedOpen,
                cfg.n_modes,
                cfg.extra_slices,
            )
            .expect("authority")
            .plane_impedance;
            let zc = C64::new(z.re, z.im);
            let denom = C64::ONE + C64::new(0.0, -omega * c_cup) * zc;
            mags.push((f, (zc * denom.recip()).abs()));
        }
        let mut peaks = Vec::new();
        for i in 1..mags.len() - 1 {
            if mags[i].1 > mags[i - 1].1 && mags[i].1 > mags[i + 1].1 {
                peaks.push(mags[i].0);
            }
        }
        peaks
    }

    fn run_note(voice: &mut BrassVoice, combo: usize, tension: f64, blow: f64) -> (f64, f64, f64) {
        voice.reset().expect("reset");
        voice
            .apply(BrassControl::SetValve {
                combo,
                fade_samples: 0,
            })
            .expect("valve");
        voice
            .apply(BrassControl::SetLipTension(tension))
            .expect("t");
        voice
            .apply(BrassControl::SetBlowingPressure(blow))
            .expect("p");
        let mut block = vec![0.0f64; 2400];
        for _ in 0..8 {
            voice.step_block(&mut block).expect("block");
        }
        let mean = block.iter().sum::<f64>() / block.len() as f64;
        let amp = block
            .iter()
            .map(|p| (p - mean).abs())
            .fold(0.0f64, f64::max);
        let defect = voice.diagnostics().last().expect("diag").supply_defect;
        (dominant_hz(&block, 24_000.0), amp, defect)
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one shared bank serves five emergent gates
    fn bl_001_005_emergent_gates() {
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        lip.provenance = "card-shaped AUTHORED lip (Estimate): mass/stiffness/gap/width at \
                          reduce-lab scale; mouth-side face area authored above the \
                          orifice-strip area"
            .to_string();
        let gas = air(293.15);
        let combos = [brass_combo(0.0), brass_combo(0.08), brass_combo(0.15)];
        let cup_volume = 1.5e-6;
        let mut voice = BrassVoice::new(
            &combos,
            &["open", "crook-1", "crook-2"],
            &gas,
            &load,
            &cfg,
            lip.clone(),
            cup_volume,
        )
        .expect("voice");
        let c_cup = cup_volume / (gas.density * gas.sound_speed * gas.sound_speed);
        let peaks = cup_loaded_peaks_hz(&combos[0], &gas, &cfg, c_cup);
        assert!(
            peaks.len() >= 2,
            "need at least two column peaks: {peaks:?}"
        );

        // ---- bl-001 LOCK ----
        let (f_lock, amp, defect) = run_note(&mut voice, 0, 0.8, 3000.0);
        let f_lip = (350.0f64 * 0.8 * 0.8 / 1.8e-4).sqrt() / core::f64::consts::TAU;
        let nearest_peak = peaks
            .iter()
            .copied()
            .min_by(|x, y| {
                (x - f_lock)
                    .abs()
                    .partial_cmp(&(y - f_lock).abs())
                    .expect("f")
            })
            .expect("peak");
        let above_peak_cents = cents(nearest_peak, f_lock);
        let lip_distance_cents = cents(f_lip, f_lock).abs();
        let pass1 = amp > 3000.0
            && (0.0..300.0).contains(&above_peak_cents)
            && lip_distance_cents > 150.0
            && defect < 1e-4;
        verdict(
            "bl-001-emergent-lock",
            pass1,
            &format!(
                "lock {f_lock:.1} Hz, amp {amp:.0} Pa; {above_peak_cents:.0} cents above the \
                 cup-loaded peak {nearest_peak:.1} (outward-striking band), {lip_distance_cents:.0} \
                 cents from the lip alone ({f_lip:.0}); worst supply defect {defect:.1e}"
            ),
        );

        // ---- bl-002 HARMONIC SERIES ----
        let (f_slot2, amp2, _) = run_note(&mut voice, 0, 1.7, 3000.0);
        let ratio = f_slot2 / f_lock;
        let peak2 = peaks
            .iter()
            .copied()
            .min_by(|x, y| {
                (x - f_slot2)
                    .abs()
                    .partial_cmp(&(y - f_slot2).abs())
                    .expect("f")
            })
            .expect("peak");
        let above2 = cents(peak2, f_slot2);
        let pass2 = amp2 > 1500.0
            && ratio > 1.9
            && ratio < 2.4
            && (0.0..300.0).contains(&above2)
            && (peak2 / nearest_peak) > 1.8;
        verdict(
            "bl-002-harmonic-series",
            pass2,
            &format!(
                "tension walks the lock {f_lock:.1} -> {f_slot2:.1} Hz (ratio {ratio:.2}) across \
                 distinct column peaks {nearest_peak:.1} -> {peak2:.1}"
            ),
        );

        // ---- bl-003 THE VALVE TEST ----
        let (f_open, a_open, _) = run_note(&mut voice, 0, 0.8, 3000.0);
        let (f_c1, a_c1, _) = run_note(&mut voice, 1, 0.8, 3000.0);
        let (f_c2, a_c2, _) = run_note(&mut voice, 2, 0.8, 3000.0);
        let d1 = cents(f_open, f_c1);
        let d2 = cents(f_open, f_c2);
        let pass3 = a_open > 3000.0
            && a_c1 > 3000.0
            && a_c2 > 3000.0
            && d1 < -60.0
            && d1 > -220.0
            && d2 < d1
            && d2 > -320.0;
        verdict(
            "bl-003-valve-drops-the-series",
            pass3,
            &format!(
                "open {f_open:.1} Hz; +crook1 {f_c1:.1} ({d1:.0} cents); +crook2 {f_c2:.1} \
                 ({d2:.0} cents) — the whole series drops with inserted length, no pitch \
                 was ever assigned"
            ),
        );

        // ---- bl-004 SLOTTING ----
        let (f_a, amp_a, _) = run_note(&mut voice, 0, 1.7, 3000.0);
        let (f_b, amp_b, _) = run_note(&mut voice, 0, 1.9, 3000.0);
        let bend = cents(f_a, f_b);
        let slot_jump = cents(f_lock, f_a);
        let pass4 =
            amp_a > 1500.0 && amp_b > 1500.0 && bend > 30.0 && bend < 150.0 && slot_jump > 900.0;
        verdict(
            "bl-004-slotting",
            pass4,
            &format!(
                "within slot 2 the lip bends the pitch {f_a:.1} -> {f_b:.1} Hz (+{bend:.0} \
                 cents, continuous); between slots it jumps {slot_jump:.0} cents"
            ),
        );

        // ---- bl-005 c(T) ----
        let cold_gas = air(277.15);
        let mut cold_voice = BrassVoice::new(
            &combos[..1],
            &["open"],
            &cold_gas,
            &load,
            &cfg,
            lip,
            cup_volume,
        )
        .expect("cold voice");
        let (f_cold, amp_cold, _) = run_note(&mut cold_voice, 0, 0.8, 3000.0);
        let warp = cents(f_open, f_cold);
        let predicted = 600.0 * (277.15f64 / 293.15).log2();
        let pass5 = amp_cold > 3000.0 && (warp - predicted).abs() < 20.0;
        verdict(
            "bl-005-cold-horn-plays-flat",
            pass5,
            &format!(
                "cold (277 K) lock {f_cold:.1} vs warm {f_open:.1}: {warp:.1} cents \
                 (sound-speed prediction {predicted:.1}); the temperature enters ONLY \
                 through the gas card"
            ),
        );
    }

    #[test]
    fn bl_006_no_click_valve_change() {
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let gas = air(293.15);
        let combos = [brass_combo(0.0), brass_combo(0.08)];
        let run = |hard_cut: Option<bool>| -> (Vec<f64>, usize) {
            let mut voice = BrassVoice::new(
                &combos,
                &["open", "crook-1"],
                &gas,
                &load,
                &cfg,
                lip.clone(),
                1.5e-6,
            )
            .expect("voice");
            voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
            voice
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut signal = Vec::new();
            let mut block = vec![0.0f64; 600];
            let mut carried = 0usize;
            for b in 0..24usize {
                if b == 14 {
                    match hard_cut {
                        None => {}
                        Some(false) => {
                            voice
                                .apply(BrassControl::SetValve {
                                    combo: 1,
                                    fade_samples: 96,
                                })
                                .expect("switch");
                            carried = voice.lift_log().last().expect("lift").carried;
                        }
                        Some(true) => {
                            // The falsifier: a fresh silent voice replaces
                            // the sounding one.
                            voice = BrassVoice::new(
                                &combos,
                                &["open", "crook-1"],
                                &gas,
                                &load,
                                &cfg,
                                lip.clone(),
                                1.5e-6,
                            )
                            .expect("voice");
                            voice
                                .apply(BrassControl::SetValve {
                                    combo: 1,
                                    fade_samples: 0,
                                })
                                .expect("switch");
                            voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
                            voice
                                .apply(BrassControl::SetBlowingPressure(3000.0))
                                .expect("p");
                        }
                    }
                }
                voice.step_block(&mut block).expect("block");
                signal.extend_from_slice(&block);
            }
            (signal, carried)
        };
        let (control, _) = run(None);
        let (lifted, carried) = run(Some(false));
        let (cut, _) = run(Some(true));
        let jump = |p: &[f64]| -> f64 {
            (14 * 600 - 30..14 * 600 + 300)
                .map(|k| (p[k] - p[k - 1]).abs())
                .fold(0.0f64, f64::max)
        };
        let block_rms = |p: &[f64], b: usize| -> f64 {
            let seg = &p[b * 600..(b + 1) * 600];
            let mean = seg.iter().sum::<f64>() / 600.0;
            (seg.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / 600.0).sqrt()
        };
        let pre_rms = block_rms(&lifted, 12);
        // The audible discriminator is the ENVELOPE: the lift keeps the
        // note alive through the change; the hard cut drops out and
        // re-attacks.
        let lifted_min = (14..17)
            .map(|b| block_rms(&lifted, b))
            .fold(f64::INFINITY, f64::min);
        // The re-attack recovers within ~3 blocks; the dropout lives in
        // the FIRST post-switch block.
        let cut_min = block_rms(&cut, 14);
        let control_jump = jump(&control);
        let lifted_jump = jump(&lifted);
        let pass = carried > 0
            && lifted_jump < 3.0 * control_jump
            && lifted_min > 0.75 * pre_rms
            && cut_min < 0.55 * pre_rms;
        verdict(
            "bl-006-no-click-valve",
            pass,
            &format!(
                "carried {carried}; slew control {control_jump:.0} vs lifted {lifted_jump:.0}; \
                 envelope through the change: lifted min RMS {lifted_min:.0} vs pre {pre_rms:.0} \
                 (no dropout), hard-cut falsifier collapses to {cut_min:.0}"
            ),
        );
    }

    #[test]
    fn bl_007_structural_honesty_and_refusals() {
        // PITCH IS NEVER ASSIGNED + no MIDI + no EQ + no mesh: grep the
        // module source (the structurally-enforced claims).
        let full = include_str!("brass_loop.rs");
        // Scan the PRODUCTION module only (the gate's own assertion
        // strings would otherwise match themselves).
        let source = full
            .split("mod brass_loop_tests")
            .next()
            .expect("production half");
        let control_surface = source
            .split("pub enum BrassControl")
            .nth(1)
            .expect("control enum")
            .split("pub struct")
            .next()
            .expect("enum body");
        let no_pitch_inputs = !control_surface.to_lowercase().contains("hz")
            && !control_surface.to_lowercase().contains("freq")
            && !control_surface.to_lowercase().contains("pitch")
            && !control_surface.to_lowercase().contains("note");
        let lower = source.to_lowercase();
        let no_midi = !lower.contains("midi");
        let no_eq = !lower.contains("equalizer");
        let no_mesh = !source.contains("fs_rep_mesh") && !source.contains("Soup");
        // Refusals.
        let cfg = config();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let gas = air(293.15);
        let combos = [brass_combo(0.0)];
        let bad_cup =
            BrassVoice::new(&combos, &["open"], &gas, &load, &cfg, authored_lip(), 0.0).is_err();
        let mut bad_lip = authored_lip();
        bad_lip.mass_kg = -1.0;
        let bad_lip_refused =
            BrassVoice::new(&combos, &["open"], &gas, &load, &cfg, bad_lip, 1.5e-6).is_err();
        let mut voice = BrassVoice::new(
            &combos,
            &["open"],
            &gas,
            &load,
            &cfg,
            authored_lip(),
            1.5e-6,
        )
        .expect("voice");
        let bad_pressure = voice
            .apply(BrassControl::SetBlowingPressure(f64::NAN))
            .is_err();
        let bad_tension = voice.apply(BrassControl::SetLipTension(0.0)).is_err();
        let unknown_valve = voice
            .apply(BrassControl::SetValve {
                combo: 9,
                fade_samples: 0,
            })
            .is_err();
        let empty_block = voice.step_block(&mut []).is_err();
        let pass = no_pitch_inputs
            && no_midi
            && no_eq
            && no_mesh
            && bad_cup
            && bad_lip_refused
            && bad_pressure
            && bad_tension
            && unknown_valve
            && empty_block;
        verdict(
            "bl-007-structural-honesty",
            pass,
            &format!(
                "control surface pitch-free {no_pitch_inputs}; no MIDI {no_midi}; no EQ \
                 {no_eq}; no mesh in the loop {no_mesh}; refusals \
                 {bad_cup}/{bad_lip_refused}/{bad_pressure}/{bad_tension}/{unknown_valve}/{empty_block}"
            ),
        );
    }
}

#[cfg(test)]
mod brass_observer_tests {
    use super::brass_loop_tests::{air, config};
    use fs_bem::helmholtz::Medium;
    use fs_bem::panel3d::SpherePanels;
    use fs_bem::radiation_bake::{bake_radiation_load, lathe_profile};

    /// Scope-4 demonstration: the radiated observer's weight comes from
    /// the SAME Z_L bake's directivity table — what you hear is the
    /// bell's radiated field, evaluated per direction, not the
    /// mouthpiece pressure. At the measured lock frequency the mouth is
    /// compact (ka ~ 0.2), so near-omnidirectional weights ARE the
    /// correct physics; the broadband radiated-observer convolution is
    /// the brass-gates bead's lane.
    #[test]
    fn bl_008_radiated_observer_from_the_bake() {
        let _ = (air(293.15), config());
        let medium = Medium::air();
        let profile = [
            (0.0f64, 0.008f64),
            (0.03, 0.009),
            (0.06, 0.012),
            (0.09, 0.017),
            (0.12, 0.025),
        ];
        let bell = lathe_profile(&profile, 16).expect("lathe");
        let surface = SpherePanels::from_triangles(bell.triangles.clone()).expect("panels");
        // One bake row AT the measured lock frequency (bl-001: 208.9 Hz).
        let omega = core::f64::consts::TAU * 208.9;
        let bake = bake_radiation_load(
            &surface,
            &bell.driven_mask(),
            &[omega],
            medium,
            4,
            "brass-loop/observer/lock-frequency/v1",
        )
        .expect("bake");
        let row = &bake.rows[0];
        let on_axis = row.directivity.evaluate([1.0, 0.0, 0.0]).abs();
        let side = row.directivity.evaluate([0.0, 1.0, 0.0]).abs();
        let back = row.directivity.evaluate([-1.0, 0.0, 0.0]).abs();
        let ratio_side = on_axis / side.max(1e-300);
        let pass = on_axis.is_finite()
            && side.is_finite()
            && on_axis > 0.0
            && row.directivity.captured_fraction > 0.98
            && ratio_side > 0.5
            && ratio_side < 3.0;
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"bl-008-radiated-observer\",\"verdict\":\"{}\",\
             \"detail\":\"observer weights from the bake's SH table at the lock frequency \
             (ka {:.2}): on-axis {on_axis:.3e}, side {side:.3e}, back {back:.3e} \
             (compact-mouth near-omni is the correct physics; captured {:.3})\"}}",
            if pass { "pass" } else { "fail" },
            row.mouth_ka,
            row.directivity.captured_fraction
        );
        assert!(pass, "observer demonstration failed");
    }
}
