//! Simultaneous prescribed-load friction on ONE retained mechanical network.
//!
//! Every port sees d_i = u_i dt - b_i^T delta_q = d_free_i - sum_j S_ij F_j.
//! The full S includes all body, spring and dashpot reactions. Bounded nonlinear
//! coordinate sweeps reuse the original scalar stick/slip solver, then check
//! ALL laws on the same endpoint before publishing a single mechanical step.
//! No sequential voice advancement, averaged material law or delayed reaction.
//!
//! This permits double-stops and several driven/braking surfaces. Normal loads
//! and surface speeds remain explicit inputs, not a solved bow-hair/normal
//! contact model. Sticking is step-averaged. Nonmonotone laws can have several
//! roots or fail to converge; order-independent/global convergence, continuous
//! no-slip and alias-free or experimentally validated sound are not claimed.

use super::*;
use super::super::network_response;

/// Compliant obstacles in the same joint solve as prescribed-load friction.
pub mod obstacles;

/// Bounded joint setup and nonlinear work, in addition to each port's limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MultiFrictionConfig {
    /// Total ports including the original port, in 1..=64.
    pub max_ports: usize,
    /// Complete nonlinear coordinate sweeps per physical sample, in 1..=128.
    pub max_sweeps: usize,
    /// Setup operations bounded by p*(n*(k+2)+(k+1)^2)+n*p^2, with n modes,
    /// k bilateral connections and p friction ports. Checked before allocation.
    pub max_setup_terms: usize,
}

/// A port's actual work, loss and law residual on the jointly accepted sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionPortFrame {
    /// Applied signed traction [N].
    pub traction_n: f64,
    /// Actual surface-minus-attachment displacement over the step [m].
    pub slip_distance_m: f64,
    /// Work supplied by this prescribed moving surface, F*u*dt [J].
    pub drive_work_j: f64,
    /// Nonnegative sliding loss F*slip, zero on sticking/release [J].
    pub dissipation_j: f64,
    /// Actual signed kinetic force residual [N], zero on other branches.
    pub constitutive_residual_n: f64,
    /// Resolved branch of the original scalar law.
    pub regime: ModalFrictionRegime,
}
impl Default for FrictionPortFrame {
    fn default() -> Self {
        Self { traction_n: 0.0, slip_distance_m: 0.0, drive_work_j: 0.0,
            dissipation_j: 0.0, constitutive_residual_n: 0.0, regime: ModalFrictionRegime::Released }
    }
}

/// One accepted clock and energy balance, not a sum of separately stepped voices.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiFrictionFrame {
    /// Original network's one-based physical clock, including prior motion.
    pub sample: u64,
    /// Original modal or acceleration pressure observation [Pa].
    pub observer_pressure_pa: f64,
    /// Modal and bilateral-spring storage [J].
    pub network_energy_j: f64,
    /// Work of the supplied external modal forces [J].
    pub external_work_j: f64,
    /// Sum of the separately reported surface supplies [J].
    pub drive_work_j: f64,
    /// Original modal and bilateral-dashpot loss [J].
    pub network_dissipation_j: f64,
    /// Sum of actual sliding losses [J].
    pub friction_dissipation_j: f64,
    /// Stored unilateral-contact potential, separate from network storage [J].
    pub obstacle_energy_j: f64,
    /// Nonadhesive loss at the optional compliant obstacles [J].
    pub obstacle_dissipation_j: f64,
    /// Change of ALL storage + all losses - external and surface work [J].
    pub energy_residual_j: f64,
    /// Original network's energy allowance, scaled by absolute surface work.
    pub energy_tolerance_j: f64,
    /// Joint sweeps used; one for the exact original single-port path.
    pub sweeps: usize,
    /// Original port first, then additional ports in their declared order.
    pub ports: Vec<FrictionPortFrame>,
}

struct Point {
    port: ModalFrictionPort,
    config: ModalFrictionConfig,
    column: Vec<f64>,
}

/// Original network and observer plus several jointly solved tangential ports.
/// A refused step retains the accepted state, observations and clock. Scratch
/// is discarded on retry; no warm-start history depends on rejected trials.
pub struct MultiFrictionModalSystem {
    core: FrictionModalSystem,
    points: Vec<Point>,
    config: MultiFrictionConfig,
    compliance: Vec<f64>,
    free_delta: Vec<f64>,
    free_slip: Vec<f64>,
    tractions: Vec<f64>,
    regimes: Vec<ModalFrictionRegime>,
    iterations: Vec<usize>,
    staged_ports: Vec<FrictionPortFrame>,
    frame: MultiFrictionFrame,
    accepted: bool,
    interfaces: Vec<Option<interface::DryFrictionBinding>>,
    obstacles: Option<obstacles::ObstacleSet>,
}
impl FrictionModalSystem {
    /// Add prescribed friction contacts without resetting the original network.
    /// Existing body/cavity vibration, clock, pressure weights and material
    /// bindings stay intact. No additional contact stores mechanical energy.
    /// Empty additions use the original single-port step exactly.
    ///
    /// # Errors
    /// Invalid ports/laws, unrepresentable response, setup limits or cancellation.
    pub fn with_additional_ports(self, additional: Vec<(ModalFrictionPort, ModalFrictionConfig)>,
        config: MultiFrictionConfig, gate: &CancelGate) -> Result<MultiFrictionModalSystem, ModalFrictionError>
    {
        MultiFrictionModalSystem::from_single(self, additional, config, gate)
    }
}
impl MultiFrictionModalSystem {
    fn from_single(core: FrictionModalSystem, additional: Vec<(ModalFrictionPort, ModalFrictionConfig)>,
        config: MultiFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        let p = additional.len().checked_add(1).ok_or_else(|| invalid("friction port count overflow"))?;
        if !(1..=64).contains(&config.max_ports) || p > config.max_ports
            || !(1..=128).contains(&config.max_sweeps) {
            return Err(invalid("simultaneous friction needs bounded nonempty port and sweep counts").into());
        }
        let n = core.mode_count();
        let k = core.network.columns.len();
        let terms = k.checked_add(1).and_then(|v| v.checked_mul(v))
            .and_then(|square| k.checked_add(2).and_then(|v| n.checked_mul(v))
                .and_then(|v| v.checked_add(square)))
            .and_then(|v| v.checked_mul(p))
            .and_then(|v| n.checked_mul(p*p).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("simultaneous friction setup count overflow"))?;
        if terms > config.max_setup_terms {
            return Err(invalid("simultaneous friction exceeds max_setup_terms").into());
        }
        let mut points = Vec::with_capacity(p);
        points.push(Point { port: core.port.clone(), config: core.config, column: core.column.clone() });
        for (port, config) in additional {
            poll(Some(gate))?;
            let column = scalar::admit_port(&core.network, &port, config)?;
            points.push(Point { port, config, column });
        }
        let mut compliance = vec![0.0; p*p];
        for j in 0..p {
            let response = network_response(&core.network, &points[j].column, gate)?;
            for i in 0..p { compliance[i*p+j] = dot(&points[i].column, &response)?; }
            if compliance[j*p+j] <= 0.0 {
                return Err(invalid("each friction port needs positive network compliance").into());
            }
        }
        poll(Some(gate))?;
        Ok(Self { core, points, config, compliance, free_delta: vec![0.0; n],
            free_slip: vec![0.0; p], tractions: vec![0.0; p],
            regimes: vec![ModalFrictionRegime::Released; p], iterations: vec![0; p],
            staged_ports: vec![FrictionPortFrame::default(); p],
            frame: MultiFrictionFrame { ports: vec![FrictionPortFrame::default(); p], ..MultiFrictionFrame::default() },
            accepted: false, interfaces: vec![None; p-1], obstacles: None })
    }
    /// Accepted components in their original basis and order.
    #[must_use] pub fn components(&self) -> &[ModalAcousticTimeModel] { self.core.components() }
    /// External modal-force dimension.
    #[must_use] pub fn mode_count(&self) -> usize { self.core.mode_count() }
    /// Original shared mechanical period [s].
    #[must_use] pub fn sample_period_s(&self) -> f64 { self.core.sample_period_s() }
    /// Preserved original accepted sample clock.
    #[must_use] pub fn samples_rendered(&self) -> u64 { self.core.samples_rendered() }
    /// All prescribed-load friction ports, including released contacts.
    #[must_use] pub fn port_count(&self) -> usize { self.points.len() }
    /// Original geometry and law, with no mutable bypass of admission.
    #[must_use] pub fn port(&self, index: usize) -> Option<&ModalFrictionPort> { self.points.get(index).map(|p| &p.port) }
    /// Actual modal, bilateral and optional obstacle stored energy [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let network = self.core.total_energy_j()?;
        match &self.obstacles {
            Some(obstacles) => finite(network + obstacles.energy(self.components())?),
            None => Ok(network),
        }
    }
    /// Settle the complete network with every driven friction port released.
    /// Independent normal obstacles remain present: their activity and reactions
    /// are solved by the original contact-equilibrium owner, not the dynamic
    /// compliance matrix. Original per-obstacle limits, joint sweep count and
    /// obstacle setup budget apply. Dashpots carry no static load.
    ///
    /// Refuses after any sample or nonzero initial vibration; all sources,
    /// observers and clocks survive. No fallible work follows state publication.
    /// Returns initial storage INCLUDING obstacles, not prior actuator work.
    pub fn initialize_released_equilibrium(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalFrictionError>
    {
        if let Some(obstacles) = &self.obstacles {
            return obstacles.initialize_equilibrium(&mut self.core.network, external,
                self.config.max_sweeps, gate).map_err(Into::into);
        }
        self.core.initialize_released_equilibrium(external, gate)
    }
    /// Last joint report; absent until this wrapper accepts a sample.
    #[must_use] pub fn last_frame(&self) -> Option<&MultiFrictionFrame> { self.accepted.then_some(&self.frame) }
    /// Existing acceleration-to-pressure coefficients, if any, unchanged.
    #[must_use] pub fn acceleration_pressure_weights(&self) -> Option<&[f64]> { self.core.acceleration_pressure_weights() }
    /// Advance every contact and body on the same physical sample.
    pub fn step(&mut self, external: &[f64], drives: &[FrictionDrive]) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        self.step_inner(external, drives, None)
    }
    /// Same transaction, with polls during prediction, sweeps and publication.
    pub fn step_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive], gate: &CancelGate)
        -> Result<&MultiFrictionFrame, ModalFrictionError> { self.step_inner(external, drives, Some(gate)) }

    fn step_inner(&mut self, external: &[f64], drives: &[FrictionDrive], gate: Option<&CancelGate>)
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        poll(gate)?;
        let p = self.points.len();
        if drives.len() != p || drives.iter().any(|d| !d.speed_m_s.is_finite()
            || !d.normal_force_n.is_finite() || d.normal_force_n < 0.0) {
            return Err(invalid("simultaneous friction requires one finite speed and nonnegative load per port").into());
        }
        if p == 1 && self.obstacles.is_none() {
            let f = *self.core.step_inner(external, drives[0], gate)?;
            self.frame.ports[0] = FrictionPortFrame { traction_n: f.traction_n, slip_distance_m: f.slip_distance_m,
                drive_work_j: f.drive_work_j, dissipation_j: f.friction_dissipation_j,
                constitutive_residual_n: f.constitutive_residual_n, regime: f.regime };
            self.frame.sample = f.sample; self.frame.observer_pressure_pa = f.observer_pressure_pa;
            self.frame.network_energy_j = f.network_energy_j; self.frame.external_work_j = f.external_work_j;
            self.frame.drive_work_j = f.drive_work_j; self.frame.network_dissipation_j = f.network_dissipation_j;
            self.frame.friction_dissipation_j = f.friction_dissipation_j;
            self.frame.energy_residual_j = f.energy_residual_j; self.frame.energy_tolerance_j = f.energy_tolerance_j;
            self.frame.sweeps = 1; self.accepted = true;
            return Ok(&self.frame);
        }
        for i in 0..p {
            poll(gate)?;
            interface::check_endpoint(self.interface(i), &self.points[i].column,
                &self.core.network.models, drives[i])?;
        }
        let before = self.total_energy_j()?;
        self.core.network.prepare_prediction(external, gate)?;
        self.free_delta.copy_from_slice(&self.core.network.free_delta);
        if let Some(obstacles) = &mut self.obstacles {
            obstacles.prepare(&self.core.network.old_q, &self.free_delta, gate)?;
        }
        for i in 0..p {
            self.free_slip[i] = finite(finite(drives[i].speed_m_s*self.sample_period_s())?
                - dot(&self.points[i].column, &self.free_delta)?)?;
        }
        self.tractions.fill(0.0); self.iterations.fill(0);
        let mut sweeps = 0;
        for sweep in 1..=self.config.max_sweeps {
            if let Some(obstacles) = &mut self.obstacles {
                obstacles.sweep(&self.tractions, self.core.network.dt, gate)?;
            }
            for i in 0..p {
                poll(gate)?;
                let mut free = self.free_slip[i];
                for j in 0..p { if j != i { free = finite(free - self.compliance[i*p+j]*self.tractions[j])?; } }
                if let Some(obstacles) = &self.obstacles { free = finite(free + obstacles.slip_shift(i)?)?; }
                let point = &self.points[i];
                let (force, regime, iterations) = scalar::solve(&point.port.law, point.config,
                    self.compliance[i*p+i], self.sample_period_s(), free, drives[i].normal_force_n, gate)?;
                self.tractions[i] = force; self.regimes[i] = regime;
                self.iterations[i] += iterations;
            }
            // Joint convergence, not just small force updates or each coordinate's
            // last local root. Temporary physical-limit failures cannot publish.
            let mut failure = self.obstacles.as_ref().and_then(|obstacles|
                obstacles.check_prediction(&self.tractions, self.core.network.dt, gate).err())
                .map(ModalFrictionError::Network);
            for i in 0..p {
                poll(gate)?;
                let mut slip = finite(self.free_slip[i] - dot(&self.compliance[i*p..(i+1)*p], &self.tractions)?)?;
                if let Some(obstacles) = &self.obstacles { slip = finite(slip + obstacles.slip_shift(i)?)?; }
                if let Err(error) = scalar::check(&self.points[i].port.law, self.points[i].config,
                    self.sample_period_s(), drives[i], self.tractions[i], slip, self.regimes[i], self.iterations[i]) {
                    if failure.is_none() { failure = Some(error); }
                }
            }
            match failure {
                None => { sweeps = sweep; break; },
                Some(error) if sweep == self.config.max_sweeps => return Err(error),
                _ => {},
            }
        }
        self.core.forces.copy_from_slice(external);
        for (point, &traction) in self.points.iter().zip(&self.tractions) {
            if traction != 0.0 {
                for (force, &b) in self.core.forces.iter_mut().zip(&point.column) { *force = finite(*force+b*traction)?; }
            }
        }
        if let Some(obstacles) = &self.obstacles { obstacles.add_forces(&mut self.core.forces, gate)?; }
        self.core.network.stage_inner(&self.core.forces, gate)?;
        let (obstacle_energy, obstacle_loss) = match &mut self.obstacles {
            Some(obstacles) => obstacles.stage(&self.core.network.candidates, self.core.network.dt, gate)?,
            None => (0.0, 0.0),
        };
        let (mut loss, mut work, mut absolute_work) = (0.0, 0.0, 0.0);
        for i in 0..p {
            poll(gate)?;
            let distance = finite(drives[i].speed_m_s*self.sample_period_s())?;
            let slip = finite(distance-dot(&self.points[i].column, &self.core.network.free_delta)?)?;
            let (dissipation, residual) = scalar::check(&self.points[i].port.law, self.points[i].config,
                self.sample_period_s(), drives[i], self.tractions[i], slip, self.regimes[i], self.iterations[i])?;
            interface::check_average(self.interface(i), drives[i], slip, self.sample_period_s())?;
            interface::check_endpoint(self.interface(i), &self.points[i].column,
                &self.core.network.candidates, drives[i])?;
            let supply = finite(self.tractions[i]*distance)?;
            loss = finite(loss+dissipation)?; work = finite(work+supply)?;
            absolute_work = finite(absolute_work+supply.abs())?;
            self.staged_ports[i] = FrictionPortFrame { traction_n: self.tractions[i], slip_distance_m: slip,
                drive_work_j: supply, dissipation_j: dissipation, constitutive_residual_n: residual, regime: self.regimes[i] };
        }
        let pressure = self.core.staged_observer_pressure(gate)?;
        let staged = &self.core.network.staged_frame;
        let network_energy = finite(staged.modal_energy_j+staged.connection_energy_j)?;
        let after = if self.obstacles.is_some() { finite(network_energy+obstacle_energy)? } else { network_energy };
        let network_loss = finite(staged.component_dissipation_j+staged.connection_dissipation_j)?;
        let external_work = dot(external, &self.core.network.free_delta)?;
        let residual = if self.obstacles.is_some() {
            finite((after-before)+network_loss+loss+obstacle_loss-external_work-work)?
        } else { finite((after-before)+network_loss+loss-external_work-work)? };
        let scale = before.max(after).max(network_loss).max(loss).max(obstacle_loss).max(external_work.abs()).max(absolute_work);
        if self.obstacles.is_some() {
            limit("total energy including friction obstacles", after, self.core.network.config.maximum_total_energy_j)?;
        }
        let tolerance = finite(self.core.network.config.energy_absolute_tolerance_j
            +self.core.network.config.energy_relative_tolerance*scale)?;
        if loss < 0.0 || residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance }.into());
        }
        let sample = staged.sample;
        poll(gate)?;
        self.core.network.staged_frame.observer_pressure_pa = pressure;
        self.core.network.publish_staged();
        if let Some(obstacles) = &mut self.obstacles { obstacles.publish(); }
        self.frame.sample = sample; self.frame.observer_pressure_pa = pressure;
        self.frame.network_energy_j = network_energy; self.frame.external_work_j = external_work;
        self.frame.drive_work_j = work; self.frame.network_dissipation_j = network_loss;
        self.frame.friction_dissipation_j = loss; self.frame.energy_residual_j = residual;
        self.frame.obstacle_energy_j = obstacle_energy; self.frame.obstacle_dissipation_j = obstacle_loss;
        self.frame.energy_tolerance_j = tolerance; self.frame.sweeps = sweeps;
        self.frame.ports.copy_from_slice(&self.staged_ports); self.accepted = true;
        Ok(&self.frame)
    }
}

impl MultiFrictionModalSystem {
    /// Attach an original interface to one unchanged physical port. Existing
    /// motion, observers, other interfaces and clocks are preserved. A source
    /// can be attached only once and must carry the identical retained law;
    /// previous unbound samples are not retroactively source-validated.
    pub fn with_interface(mut self, index: usize, binding: interface::DryFrictionBinding)
        -> Result<Self, ModalFrictionError>
    {
        let point = self.points.get(index).ok_or_else(|| invalid("interface names an unknown friction port"))?;
        if self.interface(index).is_some() || point.port.law != *binding.card().friction_law() {
            return Err(invalid("interface must match the unchanged physical port law and cannot replace a retained source").into());
        }
        if index == 0 {
            self.core = self.core.with_interface(binding)?;
        } else {
            self.interfaces[index-1] = Some(binding);
        }
        Ok(self)
    }
    /// Source belonging to this physical contact, in original drive order.
    /// The first binding remains in the original scalar owner, including its
    /// exact singleton execution path; it is not copied into a competing slot.
    #[must_use]
    pub fn interface(&self, index: usize) -> Option<&interface::DryFrictionBinding> {
        if index == 0 { self.core.interface() }
        else { self.interfaces.get(index-1).and_then(Option::as_ref) }
    }
}
