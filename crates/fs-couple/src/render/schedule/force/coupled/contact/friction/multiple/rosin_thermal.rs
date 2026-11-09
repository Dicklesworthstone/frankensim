//! Smith & Woodhouse (2000) thermal rosin on selected ports of the simultaneous
//! owner (`fs_tribo::rosin_thermal`, the paper's Eqs. 6, 10-16 and 21).
//!
//! Friction. Each selected port is one cell of the SAME rosin contact: the
//! points of the string under the hair. Over a step its law is the set-valued
//! Coulomb law `|F_i| <= mu_i N_i`, with `mu_i = pi a^2 k_y(T_i)/N_ref` from the
//! cell's temperature at the START of the step, as the paper's law takes it.
//! `a` is the whole contact's radius at the whole load `N_ref`, and the area
//! grows with load, so `mu_i` does not depend on the load. The coefficient
//! enters the existing joint solve through the port's own Coulomb law.
//!
//! Heat. The heat balance advances only after the joint step is accepted. Each
//! cell's history steps through the law's own `advance` at zero compliance with
//! that port's accepted normal load and slip. There the law's slip is the slip
//! given, so the temperature is the law's Eq. 6 step, not a second balance. A
//! hair point crosses the contact across the string, so its crossing length is
//! the same in every cell: each cell's law has radius `a` at reference load
//! `N_ref/cells`, which gives the whole contact's radius at an even split.
//! Refining the cells then leaves the transit and advection terms unchanged.
//!
//! One stroke per cell: landing and lift. A cell is loaded on an accepted step
//! when its port reports a positive normal load and is not released.
//! - Before its first loaded step the cell is inert, since there is no contact.
//!   Its history starts at that step with the law's own initial state, the
//!   rosin at ambient, as at t = 0 of the paper. Its t = 0 is its first
//!   contact, not the first sample.
//! - Once the load vanishes (the lift), the friction force is mu N = 0 whatever
//!   mu is, so the cell adds nothing. Its history is frozen and reported.
//! - A second landing after a lift refuses by name. What the rosin does off the
//!   string (its cooling) is unsourced, so this law does not cover it.
//! Each cell keeps its own phase; the accepted step indices are in the receipt.
//!
//! Refusals, by name. A cell's second landing, an unloaded cell that still
//! carries traction, a nonfinite or negative normal load, a temperature above
//! the calibration's limit, the law's history capacity, a port that is not
//! plain Coulomb, and a law step other than the sample period. Loaded sliding
//! also refuses negative actual traction-slip work or declared numerical work:
//! a reversal's signed work cannot become heat by taking an absolute slip.
//! A refusal leaves every thermal history unchanged. The mechanical owner has
//! already accepted its frame and remains responsible for joint publication.
//! A stuck cell takes zero slip: its numerical creep is the declared stick work
//! of the receipt, not friction heat.
use super::*;
use fs_tribo::rosin_thermal::{SmithWoodhouse, SmithWoodhouseParams, SmithWoodhouseState};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// Where a cell is on its one stroke. Steps are 0-based accepted-step indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RosinCellPhase {
    /// No loaded step yet: inert, its history not started.
    BeforeLanding,
    /// Loaded from step `landed` on.
    InContact {
        /// The first loaded step: the history's t = 0.
        landed: usize,
    },
    /// Loaded from `landed`, unloaded from `lifted` on: zero force, history frozen.
    Lifted {
        /// The first loaded step.
        landed: usize,
        /// The first unloaded step after it.
        lifted: usize,
    },
}

/// One Smith & Woodhouse history per selected port, all on one declared law.
#[derive(Clone, Debug)]
pub struct RosinThermalPorts {
    law: SmithWoodhouse,
    ports: Vec<usize>,
    states: Vec<SmithWoodhouseState>,
    phases: Vec<RosinCellPhase>,
    // Accepted steps so far, loaded or not.
    accepted: usize,
    // pi a^2 / N_ref of the whole contact [m^2/N].
    area_per_load: f64,
    coefficients: Vec<f64>,
    heat_seconds: f64,
}

impl RosinThermalPorts {
    /// Bind `ports` (distinct, each a plain Coulomb port of `system`) to one law
    /// whose `contact` is the WHOLE contact: `radius_m` at `reference_normal_n`.
    /// Each cell's history uses that radius at `reference_normal_n/cells`. The
    /// law's step must equal the network's sample period. Every history starts
    /// at ambient, and every port takes `pi a^2 k_y(0)/N_ref` until its landing.
    pub fn new(mut params: SmithWoodhouseParams, ports: Vec<usize>, system: &mut MultiFrictionModalSystem)
        -> Result<Self, ModalFrictionError>
    {
        let cells = ports.len();
        if cells == 0 || ports.iter().enumerate().any(|(k, p)| *p >= system.port_count() || ports[..k].contains(p)) {
            return Err(invalid("rosin thermal friction needs distinct existing ports").into());
        }
        if params.step_s.to_bits() != system.sample_period_s().to_bits() {
            return Err(invalid("rosin thermal law step differs from the network sample period").into());
        }
        let whole = params.contact;
        #[allow(clippy::cast_precision_loss)]
        let share = cells as f64;
        params.contact.reference_normal_n = finite(whole.reference_normal_n / share)?;
        let law = SmithWoodhouse::new(params)?;
        let area_per_load = finite(core::f64::consts::PI * whole.radius_m * whole.radius_m / whole.reference_normal_n)?;
        let states = vec![law.initial_state(); cells];
        let mut bound = Self { law, ports, states, phases: vec![RosinCellPhase::BeforeLanding; cells], accepted: 0,
            area_per_load, coefficients: vec![0.0; cells], heat_seconds: 0.0 };
        bound.apply(system)?;
        Ok(bound)
    }

    /// Set every bound port's coefficient from its cell's current temperature.
    /// Call before each joint step. Every port admits its coefficient before
    /// any is set, so a refusal leaves the system and the receipt unchanged.
    pub fn apply(&mut self, system: &mut MultiFrictionModalSystem) -> Result<(), ModalFrictionError> {
        let mut next = Vec::with_capacity(self.ports.len());
        for (&port, state) in self.ports.iter().zip(&self.states) {
            let yield_pa = self.law.yield_stress_pa(state.temperature_rise_k())?;
            let mu = finite(self.area_per_load * yield_pa)?;
            system.admit_port_coulomb_coefficient(port, mu)?;
            next.push(mu);
        }
        for (&port, &mu) in self.ports.iter().zip(&next) {
            system.set_port_coulomb_coefficient(port, mu)?;
        }
        self.coefficients = next;
        Ok(())
    }

    /// Advance every loaded cell's heat on the accepted `frame`, with each
    /// port's own accepted normal load and slip, and move each cell along its
    /// stroke (landing, lift). All cells commit together or none does.
    pub fn accept(&mut self, frame: &MultiFrictionFrame) -> Result<(), ModalFrictionError> {
        let started = std::time::Instant::now();
        let step_s = self.law.params().step_s;
        let step = self.accepted;
        // The law's own force (its whole-contact limit) is never used at zero
        // compliance; the ceiling only admits it.
        let solve = TangentialSolveConfig { max_iterations: 1, maximum_force_n: f64::MAX,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-12 };
        // None: the cell's history is unchanged on this step.
        let mut staged = Vec::with_capacity(self.ports.len());
        let mut phases = Vec::with_capacity(self.ports.len());
        for ((&port, state), &phase) in self.ports.iter().zip(&self.states).zip(&self.phases) {
            let report = frame.ports.get(port).ok_or_else(|| invalid("rosin thermal port is missing from the frame"))?;
            let normal = report.normal_force_n;
            if !normal.is_finite() || normal < 0.0 {
                return Err(invalid("rosin thermal friction needs a finite nonnegative normal load").into());
            }
            if normal == 0.0 || report.regime == ModalFrictionRegime::Released {
                // Inert before the landing; after the lift mu N = 0 and the history is frozen.
                if report.traction_n != 0.0 {
                    return Err(invalid("an unloaded rosin thermal cell carries traction").into());
                }
                staged.push(None);
                phases.push(match phase {
                    RosinCellPhase::InContact { landed } => RosinCellPhase::Lifted { landed, lifted: step },
                    unchanged => unchanged,
                });
                continue;
            }
            let landed = match phase {
                RosinCellPhase::BeforeLanding => step,
                RosinCellPhase::InContact { landed } => landed,
                RosinCellPhase::Lifted { .. } => return Err(invalid(
                    "rosin thermal friction: a cell lands again after its lift; the rosin's cooling off the string is unsourced").into()),
            };
            let slip = match report.regime {
                ModalFrictionRegime::Sliding => {
                    let slip = finite(report.slip_distance_m)?;
                    let work = finite(finite(report.traction_n)? * slip)?;
                    let declared = finite(report.declared_slip_work_j)?;
                    let dissipated = finite(report.dissipation_j)?;
                    // The endpoint branch can slide even when its actual
                    // displacement gives signed, non-dissipative work.
                    // The thermal law uses |slip| and cannot represent that.
                    if work < 0.0 || declared != 0.0 || dissipated < 0.0 {
                        return Err(invalid("rosin thermal friction needs dissipative sliding work").into());
                    }
                    slip
                },
                ModalFrictionRegime::Sticking => 0.0,
                _ => return Err(invalid("a step Coulomb coefficient binds only to a plain Coulomb port").into()),
            };
            let input = TangentialStepInput { normal_n: normal, slip_increment_m: [slip, 0.0],
                dt_s: step_s, contact_temperature_k: None, nominal_pressure_pa: None,
                compliance_m_per_n: 0.0, solve };
            staged.push(Some(self.law.advance(state, &input)?.state));
            phases.push(RosinCellPhase::InContact { landed });
        }
        for (state, next) in self.states.iter_mut().zip(staged) {
            if let Some(next) = next { *state = next; }
        }
        self.phases = phases;
        self.accepted += 1;
        self.heat_seconds += started.elapsed().as_secs_f64();
        Ok(())
    }

    /// The declared per-cell law (its contact carries `N_ref/cells`).
    #[must_use] pub const fn law(&self) -> &SmithWoodhouse { &self.law }
    /// Bound ports, in cell order.
    #[must_use] pub fn ports(&self) -> &[usize] { &self.ports }
    /// Each cell's temperature above ambient over its last heat step [K]: 0
    /// before its landing, frozen after its lift.
    #[must_use] pub fn temperatures_k(&self) -> Vec<f64> { self.states.iter().map(SmithWoodhouseState::temperature_rise_k).collect() }
    /// Coefficients the next step will use, in cell order.
    #[must_use] pub fn coefficients(&self) -> &[f64] { &self.coefficients }
    /// Each cell's place on its stroke, in cell order.
    #[must_use] pub fn phases(&self) -> &[RosinCellPhase] { &self.phases }
    /// Accepted steps so far, loaded or not.
    #[must_use] pub const fn accepted_steps(&self) -> usize { self.accepted }
    /// Heat steps per cell: its loaded accepted steps.
    #[must_use] pub fn heat_steps(&self) -> Vec<usize> { self.states.iter().map(SmithWoodhouseState::steps).collect() }
    /// Wall time spent advancing the heat histories [s]; a receipt, not physics.
    #[must_use] pub const fn heat_seconds(&self) -> f64 { self.heat_seconds }
}

#[cfg(test)]
mod tests;
