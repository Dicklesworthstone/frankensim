//! Several structural modes coupled through ONE slit and pressure/flow port.
//!
//! The gap is H + g.q, pressure forces are -A*dp and swept flow is -A.v.
//! A and g are separate signed projections; neither is a fitted sound gain.
//! Contact and continuous Bernoulli flow remain in their original owners. The
//! scalar solve generalizes the existing independent-blade midpoint condensation.
//! This is a fixed linear modal basis, diagonal viscous loss and a uniform
//! effective slit, not distributed closure, hereditary material or curved flow.

use super::dynamic::{ApertureDrive, ApertureFrame, ApertureState};
use super::paired::{BladeState, PairedApertureLimits};
use super::uniform::{self, MAX_MODES};
use crate::acoustic_realize::AcousticRealizeError;
use crate::unilateral_contact::distributed::ApertureContactStep;
use fs_dcontact::{ContactStorage, Obstacle};
use fs_phs::Storage;

/// One fixed-basis mode, in a physical metre-valued coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ApertureMode {
    /// Positive generalized mass [kg].
    pub mass_kg: f64,
    /// Positive generalized restoring stiffness [N/m].
    pub stiffness_n_m: f64,
    /// Explicit nonnegative diagonal viscous loss [N s/m].
    pub damping_n_s_m: f64,
    /// Signed pressure virtual-work area and swept-volume projection [m²].
    pub pressure_area_m2: f64,
    /// Signed mean-slit displacement per coordinate metre [dimensionless].
    pub gap_projection: f64,
    /// Two independent physical force ports: generalized force is this row
    /// dotted with the supplied newtons. An all-zero column is an absent port.
    pub force_projection: [f64; 2],
}

/// Immutable physical basis, fluid load, clock and existing explicit limits.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalApertureSpec {
    /// Positive effective rest gap [m].
    pub rest_opening_m: f64,
    /// Single shared slit width [m], not multiplied by the mode count.
    pub width_m: f64,
    /// 1..=32 independent modes. Missing modes are NOT inferred or rescaled.
    pub modes: Vec<ApertureMode>,
    /// Actual inlet gas density [kg/m³].
    pub density_kg_m3: f64,
    /// Characteristic pressure/flow load [Pa s/m³].
    pub impedance_pa_s_m3: f64,
    /// Fixed mechanical step [s].
    pub time_step_s: f64,
    /// Total accepted samples, including continuation; at most 2^53.
    pub max_steps: u64,
    /// Bounds apply to every mode; velocity tolerance applies to the total gap.
    pub limits: PairedApertureLimits,
}

/// Completely accepted mechanics and original pressure/flow work account.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalApertureFrame {
    /// Common slit, fluid and complete mechanical/contact energy observation.
    pub aperture: ApertureFrame,
    /// Independent endpoint coordinates in the original basis order.
    pub modes: Vec<BladeState>,
    /// Same-sample pressure difference across the source [Pa].
    pub pressure_drop_pa: f64,
    /// Single physical contact resultant [N]; generalized forces are g_i*F.
    pub contact_force_n: f64,
    /// Independent force-balance residuals, in the original coordinates [N].
    pub momentum_residual_n: Vec<f64>,
}

/// A complete modal state; no acoustic network or filter is duplicated here.
pub struct ModalAperture {
    spec: ModalApertureSpec,
    state: Vec<BladeState>,
    obstacle: Obstacle,
    contact_storage: ContactStorage,
    samples: u64,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn nonlinear(e: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(e.to_string())
}
struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}

impl ModalAperture {
    /// Admit a supplied physical modal basis and ONE effective gap contact.
    /// There is no named-instrument preset or inferred damping. Modal sign and
    /// scaling belong to the supplied projections, including the force ports.
    /// A successful construction does not validate omitted modes or mean-gap
    /// closure. Nonproportional damping and hereditary laws require other images.
    pub fn new(spec: ModalApertureSpec, state: Vec<BladeState>, obstacle: Obstacle)
        -> Result<Self, AcousticRealizeError>
    {
        let l = spec.limits;
        if !(1..=MAX_MODES).contains(&spec.modes.len()) || state.len() != spec.modes.len()
            || ![spec.rest_opening_m, spec.width_m, spec.density_kg_m3,
                spec.impedance_pa_s_m3, spec.time_step_s, l.maximum_abs_displacement_m,
                l.maximum_abs_velocity_m_s, l.maximum_energy_j, l.velocity_tolerance_m_s,
                l.relative_tolerance].iter().all(|x| x.is_finite() && *x > 0.0)
            || l.relative_tolerance >= 1.0 || !(1..=512).contains(&l.max_evaluations)
            || spec.max_steps == 0 || spec.max_steps > (1_u64 << 53)
            || !(spec.time_step_s * spec.max_steps as f64).is_finite()
            || obstacle.provenance().trim().is_empty() {
            return Err(invalid("modal aperture requires bounded modes, matched states and finite physical limits"));
        }
        for mode in &spec.modes {
            if ![mode.mass_kg, mode.stiffness_n_m].iter().all(|x| x.is_finite() && *x > 0.0)
                || !mode.damping_n_s_m.is_finite() || mode.damping_n_s_m < 0.0
                || [mode.pressure_area_m2, mode.gap_projection, mode.force_projection[0],
                    mode.force_projection[1]].iter().any(|x| !x.is_finite()) {
                return Err(invalid("modal aperture has invalid mechanics or signed physical projections"));
            }
        }
        if spec.modes.iter().all(|m| m.gap_projection == 0.0) {
            return Err(invalid("modal aperture basis does not move its declared slit"));
        }
        let gap = spec.rest_opening_m + spec.modes.iter().zip(&state)
            .map(|(m,s)| m.gap_projection*s.displacement_m).sum::<f64>();
        ApertureContactStep::new(&obstacle, gap).map_err(nonlinear)?;
        let contact_storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![obstacle.clone()])
            .map_err(nonlinear)?;
        let model = Self { spec, state, obstacle, contact_storage, samples: 0 };
        model.validate_state(&model.state)?;
        Ok(model)
    }
    /// Original physical inputs and explicit approximation limits.
    #[must_use] pub const fn spec(&self) -> &ModalApertureSpec { &self.spec }
    /// Actual independent physical coordinates, not a single equivalent reed.
    #[must_use] pub fn modal_states(&self) -> &[BladeState] { &self.state }
    /// Accepted mechanical samples only.
    #[must_use] pub const fn accepted_steps(&self) -> u64 { self.samples }
    fn gap(&self, state: &[BladeState]) -> f64 {
        self.spec.rest_opening_m + self.spec.modes.iter().zip(state)
            .map(|(m,s)| m.gap_projection*s.displacement_m).sum::<f64>()
    }
    fn energy(&self, state: &[BladeState]) -> f64 {
        self.contact_storage.hamiltonian(&[self.gap(state), 0.0])
            + self.spec.modes.iter().zip(state).map(|(m,s)|
                0.5*m.mass_kg*s.velocity_m_s*s.velocity_m_s
                + 0.5*m.stiffness_n_m*s.displacement_m*s.displacement_m).sum::<f64>()
    }
    /// Actual mechanical and contact energy [J], not accumulated input work.
    #[must_use] pub fn stored_energy_j(&self) -> f64 { self.energy(&self.state) }
    fn validate_state(&self, state: &[BladeState]) -> Result<(), AcousticRealizeError> {
        let l = self.spec.limits;
        if state.iter().any(|s| !s.displacement_m.is_finite() || !s.velocity_m_s.is_finite()
            || s.displacement_m.abs() > l.maximum_abs_displacement_m
            || s.velocity_m_s.abs() > l.maximum_abs_velocity_m_s) {
            return Err(invalid("modal aperture state exceeds its physical coordinate bounds"));
        }
        let energy = self.energy(state);
        if !self.gap(state).is_finite() || !energy.is_finite() || energy < 0.0 || energy > l.maximum_energy_j {
            return Err(invalid("modal aperture gap or stored energy exceeds its physical domain"));
        }
        Ok(())
    }
    /// Check finite forces and reject an input addressed to an absent force port.
    pub fn validate_forces(&self, forces: [f64; 2]) -> Result<(), AcousticRealizeError> {
        for i in 0..2 {
            if !forces[i].is_finite() || (forces[i] != 0.0
                && self.spec.modes.iter().all(|m| m.force_projection[i] == 0.0)) {
                return Err(invalid("modal force requires finite newtons and an existing physical projection"));
            }
        }
        Ok(())
    }
    /// Extend a finite horizon without changing any physical coordinate.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        if total <= self.spec.max_steps || total > (1_u64 << 53)
            || !(total as f64*self.spec.time_step_s).is_finite() {
            return Err(invalid("modal aperture horizon extension must increase a finite sample budget"));
        }
        self.spec.max_steps = total;
        Ok(())
    }
    /// Advance pressure and all modes through their shared contact and flow solve.
    pub fn step(&mut self, drive: ApertureDrive, forces: [f64; 2])
        -> Result<ModalApertureFrame, AcousticRealizeError>
    {
        let frame = self.preview(drive, forces)?;
        self.accept(&frame);
        Ok(frame)
    }
    pub(super) fn accept(&mut self, frame: &ModalApertureFrame) {
        self.state.copy_from_slice(&frame.modes);
        self.samples = frame.aperture.step;
    }
    pub(super) fn preview(&self, drive: ApertureDrive, forces: [f64; 2])
        -> Result<ModalApertureFrame, AcousticRealizeError>
    {
        self.validate_forces(forces)?;
        if self.samples >= self.spec.max_steps || ![drive.upstream_pressure_pa,
            drive.incoming_pressure_pa, drive.body_flow_m3_s].iter().all(|x| x.is_finite()) {
            return Err(invalid("modal aperture needs finite supplies and an unexhausted horizon"));
        }
        let n = self.state.len();
        let dt = self.spec.time_step_s;
        let h = 0.5*dt;
        let (mut free, mut mobility, mut area, mut shape, mut external) =
            ([0.0; MAX_MODES], [0.0; MAX_MODES], [0.0; MAX_MODES],
             [0.0; MAX_MODES], [0.0; MAX_MODES]);
        for (i, (m,s)) in self.spec.modes.iter().zip(&self.state).enumerate() {
            external[i] = m.force_projection[0]*forces[0] + m.force_projection[1]*forces[1];
            let d = m.mass_kg + h*m.damping_n_s_m + h*h*m.stiffness_n_m;
            mobility[i] = h/d;
            free[i] = (m.mass_kg*s.velocity_m_s - h*m.stiffness_n_m*s.displacement_m + h*external[i])/d;
            area[i] = m.pressure_area_m2; shape[i] = m.gap_projection;
            if !d.is_finite() || d <= 0.0 || !mobility[i].is_finite() || mobility[i] <= 0.0
                || !free[i].is_finite() { return Err(invalid("modal midpoint response is not representable")); }
        }
        let gap = self.gap(&self.state);
        let z = self.spec.impedance_pa_s_m3;
        let trial = uniform::solve(&free[..n], &mobility[..n], &area[..n], &shape[..n],
            gap, dt, self.spec.width_m, self.spec.density_kg_m3, z, drive,
            self.spec.limits, &self.obstacle)?;
        let next: Vec<_> = self.state.iter().zip(&trial.velocity).map(|(s,v)| BladeState {
            displacement_m: s.displacement_m + dt*v, velocity_m_s: 2.0*v - s.velocity_m_s,
        }).collect();
        self.validate_state(&next)?;
        let next_gap = self.gap(&next);
        let vm: f64 = shape[..n].iter().zip(&trial.velocity).map(|(g,v)| g*v).sum();
        let contact = ApertureContactStep::new(&self.obstacle, gap)
            .and_then(|c| c.response(next_gap, vm)).map_err(nonlinear)?;
        let opening = f64::midpoint(gap, next_gap).max(0.0);
        let jet = super::moving::volume_flow(self.spec.width_m, opening, trial.dp, self.spec.density_kg_m3);
        let swept = -area[..n].iter().zip(&trial.velocity).map(|(a,v)| a*v).sum::<f64>();
        let no_drop = drive.upstream_pressure_pa - drive.incoming_pressure_pa;
        let load_pressure = z*(jet + swept + drive.body_flow_m3_s);
        let outgoing = if no_drop.abs().max(trial.dp.abs())
            <= drive.incoming_pressure_pa.abs().max(load_pressure.abs()) {
            no_drop - trial.dp
        } else { drive.incoming_pressure_pa + load_pressure };
        let wave = (outgoing - drive.incoming_pressure_pa)/z;
        let viscous: f64 = self.spec.modes.iter().zip(&trial.velocity)
            .map(|(m,v)| m.damping_n_s_m*v*v).sum();
        let mechanical: f64 = external[..n].iter().zip(&trial.velocity).map(|(f,v)| dt*f*v).sum();
        let energy = self.energy(&next);
        let aperture = ApertureFrame {
            step: self.samples + 1, time_s: (self.samples + 1) as f64*dt,
            state: ApertureState { opening_m: next_gap, opening_velocity_m_s:
                shape[..n].iter().zip(&next).map(|(g,s)| g*s.velocity_m_s).sum() },
            midpoint_opening_m: opening, outgoing_pressure_pa: outgoing,
            bore_pressure_pa: outgoing + drive.incoming_pressure_pa,
            jet_flow_m3_s: jet, swept_flow_m3_s: swept, bore_flow_m3_s: wave,
            flow_residual_m3_s: jet + swept + drive.body_flow_m3_s - wave,
            stored_energy_j: energy, storage_change_j: energy - self.stored_energy_j(),
            dissipated_energy_j: dt*(trial.dp*jet + viscous + contact.dissipated_power),
            pressure_work_j: dt*trial.dp*(wave - drive.body_flow_m3_s), mechanical_work_j: mechanical,
            relaxation_energy_j: 0.0, relaxation_loss_j: 0.0,
        };
        let momentum: Vec<_> = self.spec.modes.iter().enumerate().map(|(i,m)|
            m.mass_kg*(next[i].velocity_m_s - self.state[i].velocity_m_s)/dt
            + m.stiffness_n_m*f64::midpoint(self.state[i].displacement_m, next[i].displacement_m)
            + m.damping_n_s_m*trial.velocity[i] + area[i]*trial.dp
            - shape[i]*contact.force - external[i]).collect();
        let flow_scale = jet.abs() + swept.abs() + drive.body_flow_m3_s.abs() + wave.abs();
        let scale = self.stored_energy_j() + energy + aperture.pressure_work_j.abs()
            + mechanical.abs() + aperture.dissipated_energy_j.abs();
        let allowance = 128.0*self.spec.limits.relative_tolerance.max(f64::EPSILON);
        if ![outgoing, aperture.bore_pressure_pa, jet, swept, wave, scale, flow_scale,
            aperture.balance_residual_j(), aperture.flow_residual_m3_s,
            aperture.dissipated_energy_j, mechanical].iter().all(|x| x.is_finite())
            || momentum.iter().any(|x| !x.is_finite())
            || aperture.balance_residual_j().abs() > allowance*scale.max(f64::MIN_POSITIVE)
            || aperture.dissipated_energy_j < -allowance*scale
            || aperture.flow_residual_m3_s.abs() > allowance*flow_scale.max(f64::MIN_POSITIVE) {
            return Err(invalid("modal aperture refused an unresolved physical work or flow balance"));
        }
        Ok(ModalApertureFrame { aperture, modes: next, pressure_drop_pa: trial.dp,
            contact_force_n: contact.force, momentum_residual_n: momentum })
    }
}
