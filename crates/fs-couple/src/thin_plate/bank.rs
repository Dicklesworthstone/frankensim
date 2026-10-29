//! Reciprocal body/cavity/exterior evolution through the existing prepared owner.
use super::piston::BaffledPistonRadiation;
use super::{
    AcousticRealizeError, CompactBody, GasState, PreparedVkBody, PreparedVkStep, VkBody,
    prepare_render_linear, render_step_error,
};
use fs_exec::CancelGate;
use fs_phs::{PortHamiltonian, PreparedStepRecord};

pub(super) struct PlateCavity {
    sys: PortHamiltonian,
    linear: fs_phs::linear_exact::LinearSplit,
    pub(super) x: Vec<f64>,
}

/// Accepted acoustic observation and complete mechanical/material/air work.
pub type PlateBankStep = PreparedVkStep;

/// Linear compact surfaces and an optional complete nonlinear plate.
///
/// Cavity and continuous piston reactions share the body's exact/Gonzalez solve.
/// All material and acoustic histories remain in the same accepted state. Source
/// fields may be configured before the first step. Independent changes to their
/// state, coefficients or topology after preparation refuse; they cannot silently
/// invalidate the cached physical owner. Observations use their accepted mirrors.
#[derive(Default)]
pub struct PlateBank {
    pub linear: Vec<CompactBody>,
    pub vk: Option<VkBody>,
    pub(super) cavity: Option<PlateCavity>,
    radiation: Vec<FaceRadiation>,
    radiation_definitions: Option<Vec<[u64; 7]>>,
    radiation_clock: Option<(GasState, u32)>,
    radiation_pending: Option<(GasState, u32, Vec<[u64; 7]>)>,
    runtime: Option<BankRuntime>,
}

struct BankRuntime {
    body: PreparedVkBody,
    prototype: PortHamiltonian,
    areas: Vec<f64>,
    definitions: Vec<[u64; 7]>,
    mirrors: Vec<[u64; 2]>,
    roots: Vec<f64>,
    vk_identity: Option<usize>,
    vk_areas: Vec<f64>,
    vk_start: usize,
    vk_extent: usize,
    cavity_start: Option<usize>,
}

struct FaceRadiation {
    source: BaffledPistonRadiation,
    members: Vec<usize>,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn definition(b: &CompactBody) -> [u64; 7] {
    [
        b.area_m2.to_bits(),
        b.mass_kg.to_bits(),
        b.drive_participation.to_bits(),
        b.omega.to_bits(),
        b.zeta.to_bits(),
        b.piston_area_m2.to_bits(),
        b.piston_face
            .as_ref()
            .map_or(0, |face| std::sync::Arc::as_ptr(face).addr() as u64),
    ]
}
fn identity(body: &VkBody) -> usize {
    core::ptr::from_ref(body.sys.storage()).cast::<()>().addr()
}
fn capacity(n: usize) -> Result<(), AcousticRealizeError> {
    if n == 0
        || n.checked_add(2)
            .and_then(|n| n.checked_mul(2))
            .is_none_or(|n| n > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM)
    {
        return Err(invalid(
            "complete plate/air state exceeds the shared dense render capacity",
        ));
    }
    Ok(())
}

impl PlateBank {
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.linear.is_empty() && self.vk.is_none()
    }

    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        self.linear
            .iter()
            .map(CompactBody::volume_velocity)
            .sum::<f64>()
            + self.vk.as_ref().map_or(0.0, VkBody::volume_velocity)
    }

    /// Includes every attached material and acoustic history exactly once.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        if let Some(runtime) = &self.runtime {
            return runtime.body.stored_energy_j();
        }
        self.linear
            .iter()
            .map(|b| 0.5 * b.mass_kg * (b.v * b.v + b.omega * b.omega * b.y * b.y))
            .sum::<f64>()
            + self.vk.as_ref().map_or(0.0, VkBody::stored_energy_j)
            + self
                .cavity
                .as_ref()
                .map_or(0.0, |c| c.sys.hamiltonian(&c.x))
    }

    /// Drive every original force projection and sum compact-monopole pressures
    /// from the jointly accepted interval-mean surface acceleration.
    pub fn drive_and_radiate(
        &mut self,
        force_n: f64,
        dt: f64,
        rho: f64,
        listener_m: f64,
    ) -> Result<f64, AcousticRealizeError> {
        // Independent oscillators do not need a dense coupled owner. Preserve
        // the existing O(number of modes) exact scalar path and its standalone
        // digital-filter semantics when no joint physical load was requested.
        if self.vk.is_none()
            && self.cavity.is_none()
            && self.radiation.is_empty()
            && self.radiation_pending.is_none()
            && self.runtime.is_none()
        {
            if !rho.is_finite() || rho <= 0.0 || !listener_m.is_finite() || listener_m <= 0.0 {
                return Err(invalid(
                    "compact bank observation needs positive finite density and distance",
                ));
            }
            let mut pressure = 0.0;
            for body in &mut self.linear {
                pressure += body.drive_and_radiate(
                    force_n * body.drive_participation,
                    dt,
                    rho,
                    listener_m,
                )?;
            }
            return Ok(pressure);
        }
        self.advance_controlled(force_n, dt, rho, listener_m, || false)
            .map(|r| r.pressure_pa)
    }

    /// A common accepted step and complete work/heat record. Solver refusal,
    /// invalid observation and late cancellation publish no physical history.
    /// Fixed-clock steps reuse the full exponential/Newton workspace; changing
    /// dt prepares a new numerical cache without changing the physical state.
    pub fn advance_controlled(
        &mut self,
        force_n: f64,
        dt: f64,
        rho: f64,
        listener_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlateBankStep, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if !force_n.is_finite()
            || !dt.is_finite()
            || dt <= 0.0
            || !rho.is_finite()
            || rho <= 0.0
            || !listener_m.is_finite()
            || listener_m <= 0.0
        {
            return Err(invalid(
                "plate bank needs finite force and positive duration, density and distance",
            ));
        }
        if let Some((gas, rate, definitions)) = &self.radiation_pending {
            if definitions.len() != self.linear.len()
                || self
                    .linear
                    .iter()
                    .zip(definitions)
                    .any(|(b, d)| definition(b) != *d)
            {
                return Err(invalid(
                    "plate bank changed after its radiation source was declared",
                ));
            }
            let (gas, rate) = (*gas, *rate);
            self.attach_radiation_loads_controlled(&gas, rate, &CancelGate::default())?;
            self.radiation_pending = None;
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
        }
        if let Some((gas, rate)) = self.radiation_clock {
            if rho.to_bits() != gas.density.to_bits() || dt > f64::from(rate).recip() {
                return Err(invalid(
                    "bank radiation needs its source gas and a clock no coarser than the declared period",
                ));
            }
        }
        self.validate_radiation_source()?;
        if self.runtime.is_some() {
            self.validate_runtime()?;
        }
        if self.is_empty() {
            let work = if let Some(cavity) = self.cavity.as_mut() {
                let step = fs_phs::linear_exact::step_render(
                    &cavity.sys,
                    &cavity.linear,
                    &cavity.x,
                    &[0.0],
                    dt,
                )
                .map_err(render_step_error)?;
                if cancelled() {
                    return Err(AcousticRealizeError::Cancelled);
                }
                cavity.x = step.x;
                PreparedStepRecord {
                    delta_h: step.delta_h,
                    dissipated: step.dissipated,
                    supplied: step.supplied,
                    newton_iters: step.newton_iters,
                    solver_residual: step.solver_residual,
                    solver_certificate: step.solver_certificate,
                }
            } else {
                PreparedStepRecord {
                    delta_h: 0.0,
                    dissipated: 0.0,
                    supplied: 0.0,
                    newton_iters: 0,
                    solver_residual: 0.0,
                    solver_certificate: None,
                }
            };
            return Ok(PlateBankStep {
                pressure_pa: 0.0,
                volume_velocity_m3_s: 0.0,
                work,
            });
        }
        if self.runtime.is_none() {
            let runtime = self.prepare_runtime(dt)?;
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
            self.runtime = Some(runtime);
        }
        self.validate_runtime()?;
        let runtime = self.runtime.as_mut().expect("prepared above");
        if runtime.body.time_step_s().to_bits() != dt.to_bits() {
            let sys = runtime.prototype.share_storage();
            let body = VkBody {
                linear: prepare_render_linear(&sys)?,
                sys,
                x: runtime.body.state().to_vec(),
                areas: runtime.areas.clone(),
            }
            .into_prepared(dt)?;
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
            runtime.body = body;
        }
        let record =
            runtime
                .body
                .drive_and_radiate_controlled(force_n, rho, listener_m, &mut cancelled)?;
        // The prepared owner's last cancellation check precedes its commit.
        // Only infallible accepted-state mirroring follows that boundary.
        let x = runtime.body.state();
        for (i, b) in self.linear.iter_mut().enumerate() {
            b.y = x[2 * i] / runtime.roots[i];
            b.v = x[2 * i + 1] / runtime.roots[i];
            runtime.mirrors[i] = [b.y.to_bits(), b.v.to_bits()];
        }
        if let Some(body) = self.vk.as_mut() {
            body.x
                .copy_from_slice(&x[runtime.vk_start..runtime.vk_start + runtime.vk_extent]);
        }
        if let (Some(cavity), Some(start)) = (self.cavity.as_mut(), runtime.cavity_start) {
            cavity.x.copy_from_slice(&x[start..start + 2]);
        }
        Ok(record)
    }

    /// Attach the existing uniform-pressure single-neck owner at rest. Positive
    /// area means OUTWARD surface motion: cavity input is -A v, and pressure
    /// applies +A p through the SAME signed surface map, independently of the
    /// point-force footprint. End corrections and resonance-frozen neck/mouth
    /// resistance remain with the existing cavity owner and are counted once.
    pub fn attach_cavity(
        &mut self,
        cavity: fs_scenario::HelmholtzCavity,
        gas: &GasState,
    ) -> Result<(), AcousticRealizeError> {
        if self.runtime.is_some() || self.cavity.is_some() {
            return Err(invalid(
                "cavity attachment precedes evolution and cannot reset an existing cavity",
            ));
        }
        gas.validate_coherence()
            .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        let resistance = crate::cavity_modal::resonance_neck_resistance(&cavity, gas)?;
        let sys = fs_phs::helmholtz_resonator_flow(
            cavity.volume_m3,
            cavity.neck_radius_m,
            cavity.neck_length_m,
            gas.density,
            gas.sound_speed,
            resistance,
            crate::cavity_modal::neck_ends(&cavity),
        )
        .map_err(render_step_error)?;
        self.cavity = Some(PlateCavity {
            x: vec![0.0; sys.state_dim()],
            linear: prepare_render_linear(&sys)?,
            sys,
        });
        Ok(())
    }

    /// Attach continuous passive piston histories instead of held previous-step
    /// reactions. Source geometry, gas, declared band and sampled fitting errors
    /// remain explicit in BaffledPistonRadiation. Modes harvested from one plate
    /// share one face memory, driven by their SUMMED signed volume flow. Its
    /// cross-mode reaction is retained. Distinct faces have independent loads;
    /// no mutual exterior coupling between separate faces is inferred.
    pub fn attach_radiation_loads(
        &mut self,
        gas: &GasState,
        rate: u32,
    ) -> Result<(), AcousticRealizeError> {
        self.attach_radiation_loads_controlled(gas, rate, &CancelGate::default())
    }

    // Assembly routing can choose a different complete pHS owner. Delay this
    // bank's fit until its own stepping path is actually selected; an unused
    // bank must not impose its dense capacity on another physical owner.
    pub(crate) fn declare_radiation_loads(
        &mut self,
        gas: &GasState,
        rate: u32,
    ) -> Result<(), AcousticRealizeError> {
        gas.validate_coherence()
            .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        if rate == 0 {
            return Err(invalid("radiation source needs a positive sample clock"));
        }
        if self.runtime.is_some()
            || self.radiation_clock.is_some()
            || self.radiation_pending.is_some()
        {
            return Err(invalid(
                "bank radiation source is already declared or evolving",
            ));
        }
        if self.linear.iter().any(|b| b.area_m2 != 0.0) {
            self.radiation_pending =
                Some((*gas, rate, self.linear.iter().map(definition).collect()));
        }
        Ok(())
    }

    /// Uses the original fit window [min(40,.05 Nyquist),max(.4 Nyquist,4 low)]
    /// with257 samples and at most128 poles per face. Complete histories and the
    /// shared total-state capacity are admitted before publication; none is
    /// truncated. Warm digital histories have no inferred energy coordinates and
    /// refuse migration. New acoustic histories explicitly begin quiescent.
    pub fn attach_radiation_loads_controlled(
        &mut self,
        gas: &GasState,
        rate: u32,
        gate: &CancelGate,
    ) -> Result<(), AcousticRealizeError> {
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        gas.validate_coherence()
            .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        if let Some((pending_gas, pending_rate, definitions)) = &self.radiation_pending {
            if pending_gas != gas
                || *pending_rate != rate
                || definitions.len() != self.linear.len()
                || self
                    .linear
                    .iter()
                    .zip(definitions)
                    .any(|(b, d)| definition(b) != *d)
            {
                return Err(invalid(
                    "continuous radiation differs from the bank's declared physical source",
                ));
            }
        }
        if let Some((old_gas, old_rate)) = self.radiation_clock {
            self.validate_radiation_source()?;
            if old_gas == *gas && old_rate == rate {
                return Ok(());
            }
            return Err(invalid(
                "existing bank radiation cannot reset its gas, clock or history",
            ));
        }
        if self.runtime.is_some()
            || rate == 0
            || self
                .linear
                .iter()
                .any(|b| b.rad.is_some() || b.y != 0.0 || b.v != 0.0)
        {
            return Err(invalid(
                "continuous bank radiation needs a positive clock and quiescent bodies without digital histories",
            ));
        }
        let nyquist = core::f64::consts::PI * f64::from(rate);
        let low = 40.0_f64.min(0.05 * nyquist);
        let high = (0.4 * nyquist).max(4.0 * low);
        let mut loads = Vec::with_capacity(self.linear.len());
        let mut extent = 2 * self.linear.len()
            + self.vk.as_ref().map_or(0, |v| v.x.len())
            + self.cavity.as_ref().map_or(0, |c| c.x.len());
        let mut assigned = vec![false; self.linear.len()];
        for (index, body) in self.linear.iter().enumerate() {
            if assigned[index] || body.area_m2 == 0.0 {
                continue;
            }
            let mut members = vec![index];
            assigned[index] = true;
            if let Some(face) = &body.piston_face {
                for j in index + 1..self.linear.len() {
                    let other = &self.linear[j];
                    if other
                        .piston_face
                        .as_ref()
                        .is_some_and(|f| std::sync::Arc::ptr_eq(face, f))
                    {
                        if other.piston_area_m2.to_bits() != body.piston_area_m2.to_bits() {
                            return Err(invalid(
                                "one declared piston face has inconsistent physical area",
                            ));
                        }
                        members.push(j);
                        assigned[j] = true;
                    }
                }
            }
            let radius = (body.piston_area_m2 / core::f64::consts::PI).sqrt();
            let load = BaffledPistonRadiation::new(radius, gas, [low, high], 257, 128, gate)?;
            for &i in &members {
                if self.linear[i].area_m2 != 0.0 {
                    // Passivity of the fitted continuation is not accuracy
                    // authority outside the geometry's admitted source band.
                    load.impedance(self.linear[i].omega)?;
                }
            }
            extent = extent
                .checked_add(2 * load.model().poles.len())
                .ok_or_else(|| invalid("complete piston state capacity overflow"))?;
            capacity(extent)?;
            loads.push(FaceRadiation {
                source: load,
                members,
            });
        }
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if loads.is_empty() {
            return Ok(());
        }
        self.radiation = loads;
        self.radiation_definitions = Some(self.linear.iter().map(definition).collect());
        self.radiation_clock = Some((*gas, rate));
        Ok(())
    }

    fn validate_radiation_source(&self) -> Result<(), AcousticRealizeError> {
        if let Some(definitions) = &self.radiation_definitions {
            if definitions.len() != self.linear.len()
                || self
                    .linear
                    .iter()
                    .zip(definitions)
                    .any(|(b, d)| definition(b) != *d)
            {
                return Err(invalid(
                    "plate bank geometry, modes or face identity changed after radiation attachment",
                ));
            }
        }
        Ok(())
    }

    fn prepare_runtime(&mut self, dt: f64) -> Result<BankRuntime, AcousticRealizeError> {
        let count = self.linear.len();
        if self.linear.iter().any(|b| b.rad.is_some()) {
            return Err(invalid(
                "joint bank cannot reinterpret sampled piston filter histories",
            ));
        }
        if self
            .radiation
            .iter()
            .any(|r| r.members.iter().any(|&i| i >= count))
        {
            return Err(invalid(
                "plate bank topology changed after radiation attachment",
            ));
        }
        for b in &self.linear {
            if [
                b.area_m2,
                b.mass_kg,
                b.drive_participation,
                b.omega,
                b.zeta,
                b.y,
                b.v,
            ]
            .iter()
            .any(|x| !x.is_finite())
                || b.mass_kg <= 0.0
                || b.omega <= 0.0
                || b.zeta < 0.0
            {
                return Err(invalid(
                    "joint bank requires finite physical body coefficients and state",
                ));
            }
        }
        let vk_extent = self.vk.as_ref().map_or(0, |v| v.x.len());
        let total = 2 * count
            + vk_extent
            + self.cavity.as_ref().map_or(0, |c| c.x.len())
            + self
                .radiation
                .iter()
                .map(|r| 2 * r.source.model().poles.len())
                .sum::<usize>();
        capacity(total)?;
        let roots: Vec<_> = self.linear.iter().map(|b| b.mass_kg.sqrt()).collect();
        let mut x = Vec::with_capacity(total);
        let mut areas = Vec::new();
        let mut system = if count == 0 {
            None
        } else {
            let omega: Vec<_> = self.linear.iter().map(|b| b.omega).collect();
            let zeta: Vec<_> = self.linear.iter().map(|b| b.zeta).collect();
            let drive: Vec<_> = self
                .linear
                .iter()
                .zip(&roots)
                .map(|(b, m)| b.drive_participation / m)
                .collect();
            for (b, m) in self.linear.iter().zip(&roots) {
                x.extend_from_slice(&[b.y * m, b.v * m]);
                areas.push(b.area_m2 / m);
            }
            Some(fs_phs::modal_bank_ports(&omega, &zeta, &[&drive]).map_err(render_step_error)?)
        };
        let vk_start = x.len();
        if let Some(body) = self.vk.as_mut() {
            if body.sys.port_dim() != 1
                || body.sys.state_dim() != body.x.len()
                || body.areas.is_empty()
                || 2 * body.areas.len() > body.x.len()
                || body.areas.iter().chain(&body.x).any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "joint nonlinear plate needs its original force port and complete finite state",
                ));
            }
            areas.extend_from_slice(&body.areas);
            x.extend_from_slice(&body.x);
            let copy = body.sys.share_storage();
            system = Some(match system {
                Some(first) => {
                    fs_phs::series_impedance_ports(first, copy).map_err(render_step_error)?
                }
                None => copy,
            });
        }
        let mut sys = system.ok_or_else(|| invalid("joint bank has no mechanical storage"))?;
        let mut cavity_start = None;
        if let Some(cavity) = self.cavity.as_mut() {
            cavity_start = Some(x.len());
            x.extend_from_slice(&cavity.x);
            let mut columns = vec![0.0; sys.state_dim()];
            for (k, &area) in areas.iter().enumerate() {
                columns[2 * k + 1] = -area;
            }
            sys = sys
                .with_parallel_load_columns(cavity.sys.share_storage(), &columns)
                .map_err(render_step_error)?;
        }
        for load in &self.radiation {
            let acoustic = load.source.port_hamiltonian()?;
            let mut columns = vec![0.0; sys.state_dim()];
            for &i in &load.members {
                columns[2 * i + 1] = areas[i];
            }
            x.resize(x.len() + acoustic.state_dim(), 0.0);
            sys = sys
                .with_parallel_load_columns(acoustic, &columns)
                .map_err(render_step_error)?;
        }
        let prototype = sys.share_storage();
        let body = VkBody {
            linear: prepare_render_linear(&sys)?,
            sys,
            x,
            areas: areas.clone(),
        }
        .into_prepared(dt)?;
        Ok(BankRuntime {
            body,
            prototype,
            areas,
            roots,
            vk_start,
            vk_extent,
            cavity_start,
            definitions: self.linear.iter().map(definition).collect(),
            mirrors: self
                .linear
                .iter()
                .map(|b| [b.y.to_bits(), b.v.to_bits()])
                .collect(),
            vk_identity: self.vk.as_ref().map(identity),
            vk_areas: self.vk.as_ref().map_or_else(Vec::new, |v| v.areas.clone()),
        })
    }

    fn validate_runtime(&self) -> Result<(), AcousticRealizeError> {
        let r = self.runtime.as_ref().expect("prepared above");
        if self.linear.len() != r.definitions.len()
            || self.vk.as_ref().map(identity) != r.vk_identity
        {
            return Err(invalid(
                "plate bank physical topology changed after preparation",
            ));
        }
        for (i, b) in self.linear.iter().enumerate() {
            if definition(b) != r.definitions[i]
                || [b.y.to_bits(), b.v.to_bits()] != r.mirrors[i]
                || b.rad.is_some()
            {
                return Err(invalid(
                    "plate bank body coefficients or state changed outside the accepted clock",
                ));
            }
        }
        if let Some(body) = &self.vk {
            if body.areas != r.vk_areas
                || body.x != r.body.state()[r.vk_start..r.vk_start + r.vk_extent]
            {
                return Err(invalid(
                    "plate bank nonlinear source changed outside the accepted clock",
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
