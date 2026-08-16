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
//! sample (the Gonzalez-discrete-gradient island keeps an honest
//! energy ledger); per-block diagnostics carry the lock estimate, gap
//! statistics, and the power ledger.
//!
//! Honesty ladder: every fixture input is labeled card-backed or
//! authored in [`LipIsland::provenance`]; with authored tissue numbers
//! the whole loop is Estimate (claim-side promotion happens in the
//! brass-gates bead, not here). The lip MESH never appears in this
//! loop — only the reduced card's numbers do.

use fs_material::gas::GasState;
use fs_phs::{PortHamiltonian, mass_spring_damper, step};

use crate::mm_line::{CupState, MmLineBank, MmLineConfig, MmLineError, MmLoad, cup_junction};
use fs_duct::Duct;

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
    lip_state: Vec<f64>,
    lip: LipIsland,
    tension: f64,
    blow_pa: f64,
    cup: CupState,
    cup_compliance: f64,
    p_bore_prev: f64,
    rho: f64,
    dt: f64,
    block_index: u64,
    diags: Vec<BrassBlockDiag>,
}

impl BrassVoice {
    /// Compose the voice: valve-combination ducts + the lip card + a
    /// mouthpiece-cup volume (authored; D18 electrically-short).
    ///
    /// # Errors
    /// [`MmLineError`] from the bank realization or lip admission.
    pub fn new(
        combos: &[Duct],
        labels: &[&str],
        gas: &GasState,
        load: &MmLoad<'_>,
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
        let cup_compliance = cup_volume_m3 / (gas.density * gas.sound_speed * gas.sound_speed);
        Ok(BrassVoice {
            bank,
            lip_sys,
            lip_state: vec![0.0, 0.0],
            lip,
            tension: 1.0,
            blow_pa: 0.0,
            cup: CupState::default(),
            cup_compliance,
            p_bore_prev: 0.0,
            rho: gas.density,
            dt: 1.0 / f64::from(config.sample_rate_hz),
            block_index: 0,
            diags: Vec::new(),
        })
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
                self.tension = t;
                // Stiffness ~ t^2, damping ~ t: rebuild the island with
                // the SAME state (q, p persist through an embouchure
                // move — the lip does not teleport).
                self.lip_sys = mass_spring_damper(
                    self.lip.mass_kg,
                    self.lip.stiffness_n_m * t * t,
                    self.lip.damping_n_s_m * t,
                )
                .map_err(|_| MmLineError::Invalid {
                    what: "tensioned lip island refused admission",
                })?;
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
        self.cup = CupState::default();
        self.p_bore_prev = 0.0;
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

    /// Render one block of mouthpiece pressure [Pa]. The loop each
    /// sample: lip island steps under the transmural force, the
    /// Bernoulli aperture meters the flow, the cup junction turns flow
    /// into waves, the bank returns the bore's reflection.
    ///
    /// # Errors
    /// [`MmLineError`] when any island refuses (non-finite state).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), MmLineError> {
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
        for slot in out.iter_mut() {
            // 1. Lip island under the transmural force (outward-striking:
            //    mouth overpressure opens the lip).
            let dp_lip = self.blow_pa - self.p_bore_prev;
            let force = dp_lip * self.lip.face_area_m2;
            let record = step(&self.lip_sys, &self.lip_state, &[force], dt)
                .map_err(|_| MmLineError::Invalid {
                    what: "lip island step refused",
                })?;
            worst_defect = worst_defect.max(record.supply_defect().abs());
            self.lip_state = record.x;
            let gap = (self.lip.rest_gap_m + self.lip_state[0]).max(0.0);
            gap_min = gap_min.min(gap);
            gap_max = gap_max.max(gap);
            gap_acc += gap;
            // 2. Bernoulli flow through the aperture.
            let flow =
                fs_phs::bernoulli_volume_flow(self.lip.width_m, gap, dp_lip, self.rho);
            blow_work += self.blow_pa * flow * dt;
            // 3. Cup junction + bore.
            let p_minus = self.bank.incoming();
            let (p_plus, p_now) = cup_junction(
                flow,
                p_minus,
                &mut self.cup,
                self.cup_compliance,
                zc0,
                dt,
            );
            let _ = self.bank.push(p_plus)?;
            bore_work += (p_plus * p_plus - p_minus * p_minus) / zc0 * dt;
            self.p_bore_prev = p_now;
            *slot = p_now;
        }
        let n = out.len() as f64;
        // Zero-crossing lock estimate against the block mean (the
        // mouthpiece pressure rides a DC operating point).
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
        let lock_hz = crossings_ac as f64 / (n * dt);
        self.diags.push(BrassBlockDiag {
            block: self.block_index,
            lock_hz,
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

#[cfg(test)]
mod brass_loop_tests {
    use super::*;
    use crate::mm_line::MmLoad;
    use fs_duct::{LossModel, Segment, Termination};
    use fs_material::gas::GasSpec;

    fn air(temperature_k: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101_325.0).expect("air")
    }

    fn brass_combo(crook_m: f64) -> Duct {
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

    fn config() -> MmLineConfig {
        MmLineConfig {
            sample_rate_hz: 24_000,
            n_modes: 3,
            extra_slices: 1,
            loss: LossModel::WideTube,
        }
    }

    /// AUTHORED lip island for the probe (card-backed lane wired in the
    /// battery proper).
    fn authored_lip() -> LipIsland {
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

    /// Dominant spectral peak of a block [Hz] (parabolic-refined),
    /// searched above 40 Hz.
    fn dominant_hz(block: &[f64], rate: f64) -> f64 {
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
        let shift = if denom.abs() > 1e-12 { 0.5 * (ya - yc) / denom } else { 0.0 };
        (best as f64 - shift) * df
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
                voice.apply(BrassControl::SetLipTension(tension)).expect("t");
                voice.apply(BrassControl::SetBlowingPressure(blow)).expect("p");
                let mut block = vec![0.0f64; 2400];
                for _ in 0..8 {
                    voice.step_block(&mut block).expect("block");
                }
                let d = voice.diagnostics().last().expect("diag").clone();
                let mean = block.iter().sum::<f64>() / block.len() as f64;
                let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
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
            let mut voice = BrassVoice::new(
                &combos,
                &["open"],
                &gas,
                &load,
                &cfg,
                lip,
                1.5e-6,
            )
            .expect("voice");
            for &tension in &[0.8f64, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2] {
                for &blow in &[3000.0f64, 6000.0] {
                    voice.reset().expect("reset");
                    voice.apply(BrassControl::SetLipTension(tension)).expect("t");
                    voice.apply(BrassControl::SetBlowingPressure(blow)).expect("p");
                    let mut block = vec![0.0f64; 2400];
                    for _ in 0..8 {
                        voice.step_block(&mut block).expect("block");
                    }
                    let mean = block.iter().sum::<f64>() / block.len() as f64;
                    let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
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
        let mut voice = BrassVoice::new(
            &combos,
            &["open", "c1"],
            &gas,
            &load,
            &cfg,
            lip,
            1.5e-6,
        )
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
            let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
            let f = dominant_hz(&block, 24_000.0);
            println!("VALVE combo={combo} lock={f:.1} amp={amp:.0}");
        }
    }

    #[test]
    #[ignore = "operating-point probe: valve ladder, slot walk, temperature"]
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
        let mut run = |voice: &mut BrassVoice, combo: usize, tension: f64| -> (f64, f64) {
            voice.reset().expect("reset");
            voice
                .apply(BrassControl::SetValve {
                    combo,
                    fade_samples: 0,
                })
                .expect("valve");
            voice.apply(BrassControl::SetLipTension(tension)).expect("t");
            voice
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut block = vec![0.0f64; 2400];
            for _ in 0..8 {
                voice.step_block(&mut block).expect("block");
            }
            let mean = block.iter().sum::<f64>() / block.len() as f64;
            let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
            (dominant_hz(&block, 24_000.0), amp)
        };
        for &blow in &[1500.0f64, 2000.0, 2500.0, 3000.0] {
            for combo in 0..3usize {
                voice.reset().expect("reset");
                voice
                    .apply(BrassControl::SetValve {
                        combo,
                        fade_samples: 0,
                    })
                    .expect("valve");
                voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
                voice.apply(BrassControl::SetBlowingPressure(blow)).expect("p");
                let mut block = vec![0.0f64; 2400];
                for _ in 0..8 {
                    voice.step_block(&mut block).expect("block");
                }
                let mean = block.iter().sum::<f64>() / block.len() as f64;
                let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
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
            voice_cold.apply(BrassControl::SetLipTension(0.8)).expect("t");
            voice_cold
                .apply(BrassControl::SetBlowingPressure(3000.0))
                .expect("p");
            let mut block = vec![0.0f64; 2400];
            for _ in 0..8 {
                voice_cold.step_block(&mut block).expect("block");
            }
            let mean = block.iter().sum::<f64>() / block.len() as f64;
            let amp = block.iter().map(|p| (p - mean).abs()).fold(0.0f64, f64::max);
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
                voice.apply(BrassControl::SetLipTension(tension)).expect("t");
                voice.apply(BrassControl::SetBlowingPressure(blow)).expect("p");
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
}
