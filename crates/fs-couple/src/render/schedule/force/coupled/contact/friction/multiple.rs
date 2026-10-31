//! Simultaneous prescribed- or contact-loaded friction on ONE mechanical network.
//!
//! Every port sees d_i = u_i dt - b_i^T delta_q = d_free_i - sum_j S_ij F_j.
//! The full S includes all body, spring and dashpot reactions. Bounded nonlinear
//! coordinate sweeps reuse the original scalar stick/slip solver, then check
//! ALL laws on the same endpoint before publishing a single mechanical step.
//! No sequential voice advancement, averaged material law or delayed reaction.
//! Bound geometric storage wraps these same trials in the shared DG iteration;
//! normal and tangential laws, material reaction and total work must all agree.
//!
//! Normal loads are prescribed by default. Explicit contact bindings instead
//! use the simultaneously solved compliant reaction, including automatic release
//! on separation. Surface speeds remain prescribed; no bow-hair model is inferred.
//! Sticking is step-averaged. Nonmonotone laws can have several
//! roots or fail to converge; order-independent/global convergence, continuous
//! no-slip and alias-free or experimentally validated sound are not claimed.

use super::*;
use super::super::network_response;
use crate::modal_acoustic_time::ModalAcousticState;

/// Compliant obstacles in the same joint solve as prescribed or contact-loaded friction.
pub mod obstacles;

/// One normal-contact reaction driving the corresponding tangential law.
pub mod normal_load;

/// Reciprocal moving-surface attachments in the original mechanical network.
pub mod surfaces;

mod stretching;

/// A declared compliant drive: a massless series spring-dashpot between a
/// prescribed surface and the contact (Pitteroff & Woodhouse bow hair, M4.4).
pub mod compliant_drive;

/// Quasi-static response of modes omitted from fixed friction attachments.
mod residual;

/// Active-set rescue for a joint solve the Gauss-Seidel sweep cannot converge.
mod active_set;

/// Stateful material shear in the same mechanical and acoustic transaction.
pub mod layer;

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
    /// Normal load used by this law [N]: prescribed or the SAME held contact
    /// reaction from this accepted step. Not the previous step's normal force.
    pub normal_force_n: f64,
    /// Applied signed traction [N].
    pub traction_n: f64,
    /// Actual surface-minus-attachment displacement over the step [m].
    pub slip_distance_m: f64,
    /// Work supplied by this prescribed moving surface, F*u*dt [J].
    pub drive_work_j: f64,
    /// Nonnegative dissipation [J]: dry sliding work or physical layer plastic loss.
    pub dissipation_j: f64,
    /// Physical loss in this port's series-drive dashpot [J], separately from
    /// interface loss. Already included in `friction_dissipation_j`.
    pub series_drive_dissipation_j: f64,
    /// Actual signed kinetic force residual [N], zero on other branches.
    pub constitutive_residual_n: f64,
    /// Dry contact branch, or explicitly resolved elastic/plastic material shear.
    pub regime: ModalFrictionRegime,
    /// Root disclosure for a Stribeck or hyperbolic port; `None` for laws
    /// without a disclosed rule and for a released contact.
    pub roots: Option<FrictionRoots>,
    /// Independent material-layer storage/loss, absent for unchanged dry ports.
    pub shear_layer: Option<layer::ShearLayerFrame>,
}
impl Default for FrictionPortFrame {
    fn default() -> Self {
        Self { normal_force_n: 0.0, traction_n: 0.0, slip_distance_m: 0.0, drive_work_j: 0.0,
            dissipation_j: 0.0, series_drive_dissipation_j: 0.0,
            constitutive_residual_n: 0.0, regime: ModalFrictionRegime::Released,
            roots: None, shear_layer: None }
    }
}

/// One accepted clock and energy balance, not a sum of separately stepped voices.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiFrictionFrame {
    /// Original network's one-based physical clock, including prior motion.
    pub sample: u64,
    /// Original modal or acceleration pressure observation [Pa].
    pub observer_pressure_pa: f64,
    /// Modal, bilateral-spring, material-memory and any bound nonlinear storage [J].
    pub network_energy_j: f64,
    /// Work of the supplied external modal forces [J].
    pub external_work_j: f64,
    /// Sum of the separately reported surface supplies [J].
    pub drive_work_j: f64,
    /// Original modal, bilateral-dashpot and material-relaxation loss [J].
    pub network_dissipation_j: f64,
    /// Sum of physical dry/layer losses and any series-drive dashpot loss [J].
    pub friction_dissipation_j: f64,
    /// Reversible material shear-layer storage, separate from network storage [J].
    pub shear_layer_energy_j: f64,
    /// Stored unilateral-contact potential, separate from network storage [J].
    pub obstacle_energy_j: f64,
    /// Nonadhesive loss at the optional compliant obstacles [J].
    pub obstacle_dissipation_j: f64,
    /// Work supplied by prescribed support translation and attachment motion [J]. This is
    /// separate from bow/surface work and is not counted as contact loss.
    pub obstacle_work_j: f64,
    /// Static storage of modes the network omits, under the obstacles' reactions,
    /// R G R / 2 [J]. Zero unless an obstacle residual compliance is bound.
    pub obstacle_residual_energy_j: f64,
    /// Declared NUMERICAL loss of that quasi-static residual under step-held
    /// reactions, dR G dR / 2 >= 0 [J]; not a physical dissipation.
    pub obstacle_residual_loss_j: f64,
    /// Static storage of omitted modes under the friction tractions, F G F / 2 [J].
    pub friction_residual_energy_j: f64,
    /// NUMERICAL loss of their step-held response, dF G dF / 2 [J]. This is
    /// separate from physical friction and material dissipation.
    pub friction_residual_loss_j: f64,
    /// Change of ALL storage + all losses - external and surface work [J].
    pub energy_residual_j: f64,
    /// Original network's energy allowance, scaled by absolute surface work.
    pub energy_tolerance_j: f64,
    /// Joint sweeps used; one for the exact original single-port path.
    pub sweeps: usize,
    /// Times the active-set rescue resolved a step the Gauss-Seidel sweep could
    /// not converge (frankensim-o7qbp). Zero on every step the sweep converges.
    pub active_set_rescues: usize,
    /// One when a stop's actual endpoint missed its unchanged tolerance and the
    /// step was solved again with headroom; zero on every step that passed first.
    pub obstacle_resolves: usize,
    /// Jointly accepted nonlinear storage and reaction convergence, when bound.
    /// Refused trials leave this report and every port report unchanged.
    pub nonlinear: Option<super::nonlinear::NonlinearFrictionDiagnostics>,
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
    // McIntyre-Schumacher-Woodhouse branch per port: committed, staged by a
    // trial, and the per-port load line and disclosed roots scratch.
    branches: Vec<ContactBranch>,
    staged_branches: Vec<ContactBranch>,
    load_lines: Vec<Option<StribeckLoadLine>>,
    port_roots: Vec<Option<FrictionRoots>>,
    staged_ports: Vec<FrictionPortFrame>,
    frame: MultiFrictionFrame,
    accepted: bool,
    interfaces: Vec<Option<interface::DryFrictionBinding>>,
    // Per-call applicability maps only. The original held map still owns
    // integration/work; these true endpoint maps never become physical history.
    endpoint_controls: Option<Vec<surfaces::geometry::EndpointControls>>,
    obstacles: Option<obstacles::ObstacleSet>,
    contact_loads: Vec<Option<usize>>,
    // Optional massless series spring-dashpot drive per port (bow-hair
    // compliance, M4.4) and its staged stretch for the current trial.
    series_drives: Vec<Option<compliant_drive::SeriesCompliantDrive>>,
    staged_stretch: Vec<f64>,
    friction_residual: Option<residual::FrictionResidual>,
    // Constitutive candidates never replace accepted history until the complete
    // coupled network, all contacts and the energy ledger accept the step.
    shear_layers: Vec<Option<layer::MaterialShearLayer>>,
    // Pivot budget for the active-set rescue (frankensim-o7qbp). Set at
    // construction; a test may lower it to plant a budget-exhaustion negative.
    active_set_budget: usize,
    // Accepted steps solved again after a stop's actual-endpoint miss.
    obstacle_resolves: usize,
}
impl FrictionModalSystem {
    /// Add prescribed friction contacts without resetting the original network.
    /// Existing body/cavity vibration, clock, pressure weights, nonlinear storage
    /// and material bindings stay intact. New ports start with dry laws and no
    /// storage; material shear layers can be explicitly bound before playback.
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
    /// Same original linear network with all dry bows released. Bound normal
    /// contacts, material shear layers and compliant-drive storage refuse;
    /// this read-only AC solve does not invent their operating-point tangent.
    pub fn released_harmonic_response(
        &self, omega_rad_s: f64, force: &[fs_math::c64::C64],
        maximum_work: usize, gate: &CancelGate,
    ) -> Result<super::harmonic::ReleasedHarmonicResponse, ModalFrictionError> {
        self.admit_released_harmonic(gate)?;
        self.core.released_harmonic_response(omega_rad_s, force, maximum_work, gate)
    }

    /// Independent load cases on the same admitted released-contact network,
    /// sharing its original harmonic factorization and preserving RHS order.
    pub fn released_harmonic_responses(
        &self, omega_rad_s: f64, forces: &[&[fs_math::c64::C64]],
        maximum_work: usize, gate: &CancelGate,
    ) -> Result<Vec<super::harmonic::ReleasedHarmonicResponse>, ModalFrictionError> {
        self.admit_released_harmonic(gate)?;
        self.core.released_harmonic_responses(omega_rad_s, forces, maximum_work, gate)
    }

    fn admit_released_harmonic(&self, gate: &CancelGate) -> Result<(), ModalFrictionError> {
        poll(Some(gate))?;
        if self.obstacles.is_some()
            || self.contact_loads.iter().any(Option::is_some)
            || self.series_drives.iter().any(Option::is_some)
            || self.shear_layers.iter().any(Option::is_some)
        {
            return Err(invalid(
                "released harmonic response requires no bound normal contact, shear layer or compliant drive",
            ).into());
        }
        Ok(())
    }

    fn from_single(core: FrictionModalSystem, additional: Vec<(ModalFrictionPort, ModalFrictionConfig)>,
        config: MultiFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if core.adhesive.is_some() {
            return Err(invalid("simultaneous dry ports cannot discard a retained adhesive patch").into());
        }
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
        if !core.network.is_harmonic_only() {
            for j in 0..p {
                let response = network_response(&core.network, &points[j].column, gate)?;
                for i in 0..p { compliance[i*p+j] = dot(&points[i].column, &response)?; }
                if compliance[j*p+j] <= 0.0 {
                    return Err(invalid("each friction port needs positive network compliance").into());
                }
            }
        }
        poll(Some(gate))?;
        let mut branches = vec![ContactBranch::Fresh; p];
        branches[0] = core.branch;
        Ok(Self { core, points, config, compliance, free_delta: vec![0.0; n],
            free_slip: vec![0.0; p], tractions: vec![0.0; p],
            regimes: vec![ModalFrictionRegime::Released; p], iterations: vec![0; p],
            staged_branches: branches.clone(), branches,
            load_lines: vec![None; p], port_roots: vec![None; p],
            staged_ports: vec![FrictionPortFrame::default(); p],
            frame: MultiFrictionFrame { ports: vec![FrictionPortFrame::default(); p], ..MultiFrictionFrame::default() },
            accepted: false, interfaces: vec![None; p-1], endpoint_controls: None,
            obstacles: None, contact_loads: vec![None; p],
            series_drives: vec![None; p], staged_stretch: vec![0.0; p],
            friction_residual: None,
            shear_layers: vec![None; p], active_set_budget: 4 * p + 16, obstacle_resolves: 0 })
    }
    /// Accepted steps whose stop missed at its actual endpoint and was solved
    /// again with headroom, then passed the unchanged tolerance.
    #[must_use] pub fn obstacle_resolves(&self) -> usize { self.obstacle_resolves }
    /// Accepted components in their original basis and order.
    #[must_use] pub fn components(&self) -> &[ModalAcousticTimeModel] { self.core.components() }
    /// Original constitutive parameters, retained through the joint friction solve.
    #[must_use]
    pub fn material_relaxation_branches(&self) -> &[fs_phs::RelaxationBranch] {
        self.core.material_relaxation_branches()
    }
    /// Accepted first-order material history, with no mutable admission bypass.
    #[must_use]
    pub fn material_memory(&self) -> &[f64] { self.core.material_memory() }
    /// Accepted nonlinear constitutive histories, in original binding order.
    pub fn nonlinear_material_memory(&self) -> Vec<&[f64]> {self.core.nonlinear_material_memory()}
    /// Retain the common constitutive state after contact/static preparation.
    /// Contact maps and their response matrices are unchanged by this binding.
    pub fn with_mapped_material_histories(mut self,
        bindings:Vec<super::nonlinear::history::MappedModalHistoryBinding>,
        maximum_states:usize,gate:&CancelGate) -> Result<Self,ModalFrictionError> {
        if self.friction_residual.is_some() && !bindings.is_empty() {
            return Err(invalid("friction residual compliance does not include nonlinear material storage").into());
        }
        self.core=self.core.with_mapped_material_histories(bindings,maximum_states,gate)?;
        Ok(self)
    }
    /// Source admission reads the unchanged original mechanical basis.
    pub(crate) fn original_core(&self) -> &FrictionModalSystem {&self.core}
    /// Work required to observe an original range of accepted material arms.
    pub fn material_relaxation_observation_work(&self, branches: core::ops::Range<usize>)
        -> Result<usize, ModalCouplingError>
    {
        self.core.material_relaxation_observation_work(branches)
    }
    /// Accepted arm storage and instantaneous power after the joint solve.
    pub fn material_relaxation_observation(&self, branches: core::ops::Range<usize>,
        maximum_work: usize, gate: &CancelGate)
        -> Result<MaterialRelaxationObservation, ModalCouplingError>
    {
        self.core.material_relaxation_observation(branches, maximum_work, gate)
    }
    /// Last accepted linear account, including independent material storage/loss.
    #[must_use]
    pub fn coupled_frame(&self) -> Option<&super::super::super::CoupledModalFrame> {
        self.core.coupled_frame()
    }
    /// External modal-force dimension.
    #[must_use] pub fn mode_count(&self) -> usize { self.core.mode_count() }
    /// The shared linear network's admitted stiffness bound [rad/s]; contacts are not in it.
    #[must_use] pub fn linear_frequency_upper_bound_rad_s(&self) -> f64 { self.core.linear_frequency_upper_bound_rad_s() }
    /// Original shared mechanical period [s].
    #[must_use] pub fn sample_period_s(&self) -> f64 { self.core.sample_period_s() }
    /// Preserved original accepted sample clock.
    #[must_use] pub fn samples_rendered(&self) -> u64 { self.core.samples_rendered() }
    /// All prescribed-load friction ports, including released contacts.
    #[must_use] pub fn port_count(&self) -> usize { self.points.len() }
    /// Original geometry and dry construction law, with no mutable admission bypass.
    /// An explicit material_shear_layer replaces that port's active dry mechanism.
    #[must_use] pub fn port(&self, index: usize) -> Option<&ModalFrictionPort> { self.points.get(index).map(|p| &p.port) }
    /// Actual modal, bilateral, material, series-drive and obstacle stored energy [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let network = self.core.total_energy_j()?;
        let network = if self.shear_layers.iter().any(Option::is_some) {
            finite(network + self.shear_layers.iter().flatten().map(layer::MaterialShearLayer::energy_j).sum::<f64>())?
        } else { network };
        let hair: f64 = self.series_drives.iter().flatten()
            .map(|drive| drive.spring_energy_j(drive.stretch_m())).sum();
        let network = match &self.friction_residual {
            Some(residual) => finite(network + residual.energy_j()?)?,
            None => network,
        };
        match &self.obstacles {
            Some(obstacles) => finite(network + obstacles.energy(self.components())? + hair),
            None => finite(network + hair),
        }
    }
    /// Attach a massless series spring-dashpot drive to one or more ports: the
    /// declared compliant drive of M4.4 (bow hair). Each binding is `(port,
    /// stiffness_n_m, damping_kg_s)`; the drive begins relaxed. Bind the
    /// complete set once before playback. Other ports keep the rigid prescribed
    /// drive. The series element enters the same joint load line, free slip and
    /// energy ledger as every other reaction.
    ///
    /// # Errors
    /// A sample already rendered, a rebind, an unknown port, or non-physical
    /// stiffness/damping.
    pub fn with_series_compliant_drives(mut self, drives: Vec<(usize, f64, f64)>, gate: &CancelGate)
        -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if drives.is_empty() { return Ok(self); }
        if self.samples_rendered() != 0 || self.series_drives.iter().any(Option::is_some) {
            return Err(invalid("series compliant drives need one complete binding before playback").into());
        }
        for (port, stiffness_n_m, damping_kg_s) in drives {
            poll(Some(gate))?;
            if port >= self.series_drives.len() || self.series_drives[port].is_some() {
                return Err(invalid("one series compliant drive per actual friction port").into());
            }
            self.series_drives[port] = Some(compliant_drive::SeriesCompliantDrive::new(stiffness_n_m, damping_kg_s)?);
        }
        Ok(self)
    }

    /// Replace one component's modal state before playback, for example to
    /// begin a string from a known Helmholtz motion. The component's model and
    /// its staging candidate are both set through the modal model's own restart
    /// seam, which checks the basis and the displacement/velocity/energy
    /// budgets; `prepare_prediction` re-derives the network's cached positions
    /// from the models on the next step, so no other cache needs touching. The
    /// shared clock is unchanged and any later contact history starts fresh.
    ///
    /// # Errors
    /// A sample already rendered, an unknown component, or a state that the
    /// modal basis or budgets reject.
    pub fn restore_component_states(&mut self, component: usize, states: &[ModalAcousticState])
        -> Result<(), ModalFrictionError>
    {
        if self.samples_rendered() != 0 {
            return Err(invalid("restore component states before any sample is rendered").into());
        }
        let network = &mut self.core.network;
        if component >= network.models.len() {
            return Err(invalid("restore names a component outside the coupled network").into());
        }
        if !network.plane_stress_material_points.is_empty() {
            // A cold restart is still a physical source-domain change. Stage
            // the complete accepted prefix before exposing this replacement.
            network.candidates.clone_from(&network.models);
            network.candidate_material_memory.clone_from(&network.material_memory);
            network.candidates[component].restore_states(states)
                .map_err(|source|ModalCouplingError::Component {component,source})?;
            network.validate_staged_material_points(None)?;
        }
        network.models[component].restore_states(states)
            .map_err(|source| ModalCouplingError::Component { component, source })?;
        network.candidates[component].restore_states(states)
            .map_err(|source| ModalCouplingError::Component { component, source })?;
        Ok(())
    }

    /// Settle the complete network with every driven friction port released.
    /// Independent normal obstacles remain present: their activity and reactions
    /// are solved by the original contact-equilibrium owner, not the dynamic
    /// compliance matrix. Original per-obstacle limits, joint sweep count and
    /// obstacle setup budget apply. Dashpots carry no static load.
    ///
    /// Bound nonlinear storage and normal contacts share one conservative
    /// equilibrium objective and the original nonlinear force/iteration limits.
    /// Stationary contact gradients carry no dynamic dashpot force. Original
    /// normal source domains and force/penetration limits apply to the result.
    /// Refuses after any sample or nonzero initial vibration; all sources,
    /// observers and clocks survive. No fallible work follows state publication.
    /// Returns initial storage INCLUDING obstacles, not prior actuator work.
    pub fn initialize_released_equilibrium(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalFrictionError>
    {
        if self.shear_layers.iter().any(Option::is_some) {
            return Err(invalid("released static preload cannot discard retained material-layer shear").into());
        }
        if let Some(state) = &mut self.core.nonlinear {
            let report = state.initialize_equilibrium(&mut self.core.network, external,
                self.obstacles.as_ref(), gate)?;
            let energy = report.total_energy_j;
            state.preload = Some(report);
            return Ok(energy);
        }
        if let Some(obstacles) = &self.obstacles {
            return obstacles.initialize_equilibrium(&mut self.core.network, external,
                self.config.max_sweeps, gate).map_err(Into::into);
        }
        self.core.initialize_released_equilibrium(external, gate)
    }
    /// Accepted nonlinear initialization on this same mechanical owner. The
    /// initial storage is not a reconstructed pre-window actuator work.
    pub fn nonlinear_preload_diagnostics(&self) -> Option<&super::nonlinear::NonlinearPreloadDiagnostics> {
        self.core.nonlinear_preload_diagnostics()
    }
    /// Last joint report; absent until this wrapper accepts a sample.
    #[must_use] pub fn last_frame(&self) -> Option<&MultiFrictionFrame> { self.accepted.then_some(&self.frame) }
    /// Existing acceleration-to-pressure coefficients, if any, unchanged.
    #[must_use] pub fn acceleration_pressure_weights(&self) -> Option<&[f64]> { self.core.acceleration_pressure_weights() }
    /// Advance every contact and body on the same physical sample.
    pub fn step(&mut self, external: &[f64], drives: &[FrictionDrive]) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        self.step_inner(external, drives, None, None, None)
    }
    /// Same transaction, with polls during prediction, sweeps and publication.
    pub fn step_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive], gate: &CancelGate)
        -> Result<&MultiFrictionFrame, ModalFrictionError> { self.step_inner(external, drives, None, None, Some(gate)) }

    fn step_inner(&mut self, external: &[f64], drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        gate: Option<&CancelGate>)
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        self.step_with_patch_geometry_inner(external, drives, obstacle_offsets, obstacle_geometry, false, gate)
    }

    // Only the complete normal/tangent geometry entry point may transport a
    // contact-loaded patch. Ordinary normal motion must not detach its tangent.
    fn step_with_patch_geometry_inner(&mut self, external: &[f64], drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        contact_patch_motion: bool, gate: Option<&CancelGate>)
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        if (0..self.port_count()).any(|port|
            self.interface(port).is_some_and(|source| source.temperature_source().is_some())) {
            return Err(invalid("temperature-dependent friction requires its finite-heat transaction").into());
        }
        self.prepare_step(drives, obstacle_offsets, obstacle_geometry, contact_patch_motion, gate)?;
        let p = self.points.len();
        if p == 1 && self.obstacles.is_none() && self.series_drives.iter().all(Option::is_none)
            && self.shear_layers.iter().all(Option::is_none) && self.endpoint_controls.is_none()
            && self.friction_residual.is_none() {
            let f = *self.core.step_inner(external, drives[0], gate)?;
            // A later endpoint-bound step uses the joint owner. Preserve the
            // same accepted nonmonotone contact branch across both routes.
            self.branches[0] = self.core.branch;
            self.frame.ports[0] = FrictionPortFrame { normal_force_n: drives[0].normal_force_n, traction_n: f.traction_n, slip_distance_m: f.slip_distance_m,
                drive_work_j: f.drive_work_j, dissipation_j: f.friction_dissipation_j,
                series_drive_dissipation_j: 0.0,
                constitutive_residual_n: f.constitutive_residual_n, regime: f.regime, roots: f.roots, shear_layer: None };
            self.frame.sample = f.sample; self.frame.observer_pressure_pa = f.observer_pressure_pa;
            self.frame.network_energy_j = f.network_energy_j; self.frame.external_work_j = f.external_work_j;
            self.frame.drive_work_j = f.drive_work_j; self.frame.network_dissipation_j = f.network_dissipation_j;
            self.frame.friction_dissipation_j = f.friction_dissipation_j;
            self.frame.energy_residual_j = f.energy_residual_j; self.frame.energy_tolerance_j = f.energy_tolerance_j;
            self.frame.nonlinear = self.core.nonlinear_frame().copied().map(|f| f.diagnostics());
            self.frame.sweeps = 1; self.accepted = true;
            return Ok(&self.frame);
        }
        let accepted = self.stage_constitutive(external, drives, obstacle_offsets, gate)?;
        self.validate_obstacle_sources(gate)?;
        poll(gate)?;
        Ok(self.publish_step(accepted))
    }

    // Both ordinary playback and thermal coefficient trials stage the same
    // complete constitutive owner. Neither caller may bypass its exclusions or
    // replace nonlinear string storage with a linear-only mechanical trial.
    fn stage_constitutive(&mut self, external: &[f64], drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, gate: Option<&CancelGate>)
        -> Result<MultiFrictionFrame, ModalFrictionError> {
        self.core.network.require_time_domain()?;
        if self.core.has_nonlinear_potentials() {
            if self.obstacles.as_ref().is_some_and(obstacles::ObstacleSet::has_residual) {
                return Err(invalid("bound nonlinear storage is not staged with the obstacles' residual compliance").into());
            }
            self.stage_with_nonlinear(external, drives, obstacle_offsets, gate)
        } else { self.stage_multiple(external, drives, obstacle_offsets, gate) }
    }

    // Shared by the ordinary step and the closed layer/enthalpy transaction.
    // Preparation changes only scratch, never accepted physical history.
    fn prepare_step(&mut self, drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        contact_patch_motion: bool, gate: Option<&CancelGate>) -> Result<(), ModalFrictionError> {
        if let Some(state) = &mut self.core.nonlinear { state.discard_staged_history_heat(); }
        poll(gate)?;
        if self.friction_residual.is_some() { self.admit_friction_residual()?; }
        if let Some(offsets) = obstacle_offsets {
            let obstacles = self.obstacles.as_ref()
                .ok_or_else(|| invalid("normal motion requires attached obstacles"))?;
            obstacles.check_motion(offsets, self.sample_period_s())?;
        }
        let p = self.points.len();
        if drives.len() != p || drives.iter().any(|d| !d.speed_m_s.is_finite()
            || !d.normal_force_n.is_finite() || d.normal_force_n < 0.0) {
            return Err(invalid("simultaneous friction requires one finite speed and nonnegative load per port").into());
        }
        self.validate_contact_normal_drives(drives)?;
        if obstacle_geometry.is_some() && self.obstacles.is_none() {
            return Err(invalid("attachment motion requires attached obstacles").into());
        }
        if let Some(obstacles) = &mut self.obstacles {
            // Geometry is fixed for every inner friction/contact and outer DG
            // trial of THIS sample. Accepted columns remain unchanged until
            // the complete mechanical transaction publishes.
            obstacles.prepare_geometry(&self.core.network, &self.points,
                &self.contact_loads, obstacle_geometry, contact_patch_motion, gate)?;
        }
        Ok(())
    }

    // Source admission follows all numerical trials, before any publication.
    // The thermal transaction calls this only at its final accepted endpoint.
    fn validate_obstacle_sources(&self, gate: Option<&CancelGate>) -> Result<(), ModalFrictionError> {
        self.core.network.validate_staged_material_points(gate)?;
        if let Some(obstacles) = &self.obstacles {
            obstacles.validate_sources(&self.core.network.candidates, self.sample_period_s(), gate)?;
        }
        Ok(())
    }

    // No fallible work may follow this common physical publication boundary.
    fn publish_step(&mut self, mut accepted: MultiFrictionFrame) -> &MultiFrictionFrame {
        self.core.network.staged_frame.observer_pressure_pa = accepted.observer_pressure_pa;
        self.core.network.publish_staged();
        if accepted.nonlinear.is_some() {
            if let Some(state)=&mut self.core.nonlinear {state.publish_histories();}
        }
        self.branches.copy_from_slice(&self.staged_branches);
        self.core.branch = self.branches[0];
        for (drive, &stretch) in self.series_drives.iter_mut().zip(&self.staged_stretch) {
            if let Some(drive) = drive { drive.commit(stretch); }
        }
        if let Some(residual) = &mut self.friction_residual { residual.publish(&self.tractions); }
        for layer in self.shear_layers.iter_mut().flatten() { layer.publish(); }
        if let Some(obstacles) = &mut self.obstacles { obstacles.publish(); }
        self.obstacle_resolves += accepted.obstacle_resolves;
        // Reuse the accepted port allocation through repeated constitutive
        // trials. No fallible work follows mechanical publication.
        std::mem::swap(&mut accepted.ports, &mut self.frame.ports);
        accepted.ports.copy_from_slice(&self.staged_ports);
        self.frame = accepted;
        self.accepted = true;
        &self.frame
    }

    // Stribeck and hyperbolic ports take the branch-continuation owner on
    // their own held-force load line, keeping the branch each occupies across
    // the joint sweep; the branch commits only when the step publishes. Every
    // other law, and a released port, keeps the scalar bisection unchanged.
    fn solve_point(&mut self, i: usize, free: f64, normal: f64, gate: Option<&CancelGate>)
        -> Result<(f64, ModalFrictionRegime, usize, Option<FrictionRoots>), ModalFrictionError>
    {
        let law = self.points[i].port.law.clone();
        let config = self.points[i].config;
        let dt = self.sample_period_s();
        let mut compliance = self.friction_compliance(i, i)?;
        if let Some(drive) = self.series_drives[i] {
            compliance = finite(compliance + drive.series_compliance_m_n(dt))?;
        }
        if let Some(layer) = &self.shear_layers[i] {
            poll(gate)?;
            self.staged_branches[i] = ContactBranch::Fresh;
            return layer.solve(free, compliance, normal, dt, config);
        }
        let continued = matches!(law, FrictionLaw::Stribeck { .. } | FrictionLaw::Hyperbolic { .. });
        if !continued || normal == 0.0 {
            self.staged_branches[i] = ContactBranch::Fresh;
            let (force, regime, iterations) = scalar::solve(&law, config, compliance, dt, free, normal, gate)?;
            return Ok((force, regime, iterations, None));
        }
        poll(gate)?;
        let mobility = finite(compliance/dt)?;
        let refusal = |error: ImplicitFrictionError<TriboError>| -> ModalFrictionError { match error {
            ImplicitFrictionError::Constitutive(cause) => cause.into(),
            ImplicitFrictionError::InvalidInput(what) => invalid(what).into(),
            ImplicitFrictionError::NonDissipative => invalid("set-valued friction would drive relative slip").into(),
            ImplicitFrictionError::NoConvergence => invalid("set-valued friction root did not resolve within max_iterations").into(),
        }};
        if !self.load_lines[i].as_ref().is_some_and(|line| line.matches(normal, mobility)) {
            self.load_lines[i] = Some(StribeckLoadLine::new(&law, normal, mobility).map_err(refusal)?);
        }
        let line = self.load_lines[i].as_ref().ok_or_else(|| invalid("missing set-valued load line"))?;
        let (step, iterations) = line.step_within(finite(free/dt)?, self.branches[i], config.max_iterations)
            .map_err(refusal)?;
        self.staged_branches[i] = step.branch;
        let regime = if step.branch == ContactBranch::Stuck {
            ModalFrictionRegime::Sticking
        } else {
            ModalFrictionRegime::Sliding
        };
        Ok((step.traction_n, regime, iterations, Some(FrictionRoots { found: step.roots, selected: step.selected })))
    }

    /// Active-set rescue for a step the Gauss-Seidel sweep cannot converge
    /// (frankensim-o7qbp). Builds the effective load line (bare compliance plus
    /// each port's series/hair compliance on the diagonal; free slip plus each
    /// port's stored-stretch shift), solves the unique consistent stick/slip
    /// partition directly, and on success overwrites the held tractions and
    /// regimes with it. Returns whether a consistent partition was adopted.
    ///
    /// Restricted to the plain friction case: when any port carries a unilateral
    /// obstacle or a material shear layer the rescue declines and the caller
    /// refuses with the sweep's own last failure as the decline's cause.
    /// Override the active-set pivot budget (frankensim-o7qbp). For tests that
    /// plant a budget-exhaustion negative; production uses the construction default.
    pub fn force_active_set_budget(&mut self, budget: usize) { self.active_set_budget = budget; }

    fn active_set_rescue(&mut self, drives: &[FrictionDrive], gate: Option<&CancelGate>,
        cause: ModalFrictionError) -> Result<(), ModalFrictionError>
    {
        let p = self.points.len();
        if self.obstacles.is_some() || self.shear_layers.iter().any(Option::is_some) {
            // Not attempted: the active-set model covers only the plain friction
            // ports. Name it a decline that carries the Gauss-Seidel failure.
            return Err(ModalFrictionError::ActiveSetDeclined { cause: Box::new(cause) });
        }
        let dt = self.sample_period_s();
        let mut a = self.compliance.clone();
        if let Some(residual) = &self.friction_residual {
            for (a, g) in a.iter_mut().zip(residual.compliance()) { *a = finite(*a + g)?; }
        }
        let mut f = vec![0.0; p];
        let mut normal = vec![0.0; p];
        for i in 0..p {
            if let Some(series) = self.series_drives[i] {
                a[i*p+i] = finite(a[i*p+i] + series.series_compliance_m_n(dt))?;
                f[i] = finite(self.free_slip[i] + series.free_slip_shift_m(dt))?;
            } else {
                f[i] = self.free_slip[i];
            }
            normal[i] = self.contact_loaded_drive(i, drives[i])?.normal_force_n;
        }
        let laws: Vec<FrictionLaw> = (0..p).map(|i| self.points[i].port.law.clone()).collect();
        let lawrefs: Vec<&FrictionLaw> = laws.iter().collect();
        let seed = self.regimes.clone();
        match active_set::solve_active_set(&a, &f, &normal, &lawrefs, dt, &seed, self.active_set_budget)? {
            active_set::Resolution::Solved(tractions, regimes) => {
                poll(gate)?;
                self.tractions.copy_from_slice(&tractions);
                self.regimes.copy_from_slice(&regimes);
                for b in &mut self.staged_branches { *b = ContactBranch::Fresh; }
                Ok(())
            }
            active_set::Resolution::Exhausted { pivots, last_set, worst } => {
                Err(ModalFrictionError::ActiveSetExhausted {
                    pivots,
                    budget: self.active_set_budget,
                    last_active_set: last_set,
                    worst_port: worst.port,
                    worst_condition: worst.condition,
                    worst_size: worst.size,
                })
            }
        }
    }

    // Original simultaneous friction/obstacle solve with held modal forces.
    // Only candidates and scratch change. Nonlinear storage can repeat this
    // solve until every bow and the shared constitutive reaction agree.
    //
    // Each stop is checked again at its ACTUAL endpoint after staging, which is
    // rounded along a different path than its prediction. A prediction accepted
    // at the edge of the tolerance can then miss by parts per million. Only then
    // is the same step solved again from the accepted state with headroom and
    // checked with the UNCHANGED tolerance; it refuses if it misses again. A
    // step that passes the first time is unchanged.
    fn stage_multiple(&mut self, external: &[f64], drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, gate: Option<&CancelGate>)
        -> Result<MultiFrictionFrame, ModalFrictionError>
    {
        let first = self.stage_multiple_once(external, drives, obstacle_offsets, gate);
        let Some(obstacles) = self.obstacles.as_mut() else { return first; };
        if !obstacles.take_endpoint_miss() { return first; }
        obstacles.set_headroom(true);
        let second = self.stage_multiple_once(external, drives, obstacle_offsets, gate);
        if let Some(obstacles) = &mut self.obstacles {
            obstacles.set_headroom(false);
            obstacles.take_endpoint_miss();
        }
        let mut frame = second?;
        frame.obstacle_resolves = 1;
        Ok(frame)
    }

    fn stage_multiple_once(&mut self, external: &[f64], drives: &[FrictionDrive],
        obstacle_offsets: Option<&[f64]>, gate: Option<&CancelGate>)
        -> Result<MultiFrictionFrame, ModalFrictionError>
    {
        let p = self.points.len();
        for i in 0..p {
            poll(gate)?;
            if self.contact_loads[i].is_none() {
                self.check_interface_endpoint(i, &self.core.network.models, drives[i], false)?;
            }
        }
        let before = self.total_energy_j()?;
        self.core.network.prepare_prediction(external, gate)?;
        self.free_delta.copy_from_slice(&self.core.network.free_delta);
        if let Some(obstacles) = &mut self.obstacles {
            obstacles.prepare(&self.core.network.old_q, &self.free_delta, obstacle_offsets, gate)?;
        }
        for i in 0..p {
            self.free_slip[i] = finite(finite(drives[i].speed_m_s*self.sample_period_s())?
                - dot(&self.points[i].column, &self.free_delta)?)?;
            if let Some(residual) = &self.friction_residual {
                self.free_slip[i] = finite(self.free_slip[i] + residual.old_displacement(i)?)?;
            }
        }
        self.tractions.fill(0.0); self.iterations.fill(0);
        for r in &mut self.port_roots { *r = None; }
        let mut sweeps = 0;
        let mut rescues = 0usize;
        for sweep in 1..=self.config.max_sweeps {
            if let Some(obstacles) = &mut self.obstacles {
                obstacles.sweep(&self.tractions, self.core.network.dt, gate)?;
            }
            for i in 0..p {
                poll(gate)?;
                let mut free = self.free_slip[i];
                for j in 0..p { if j != i { free = finite(free - self.friction_compliance(i,j)?*self.tractions[j])?; } }
                if let Some(obstacles) = &self.obstacles { free = finite(free + obstacles.slip_shift(i)?)?; }
                if let Some(series) = self.series_drives[i] { free = finite(free + series.free_slip_shift_m(self.sample_period_s()))?; }
                let drive = self.contact_loaded_drive(i, drives[i])?;
                let (force, regime, iterations, roots) = self.solve_point(i, free, drive.normal_force_n, gate)?;
                self.tractions[i] = force; self.regimes[i] = regime;
                self.iterations[i] += iterations; self.port_roots[i] = roots;
            }
            // Joint convergence, not just small force updates or each coordinate's
            // last local root. Temporary physical-limit failures cannot publish.
            let mut failure = self.obstacles.as_ref().and_then(|obstacles|
                obstacles.check_prediction(&self.tractions, self.core.network.dt, gate).err())
                .map(ModalFrictionError::Network);
            for i in 0..p {
                poll(gate)?;
                let mut slip = finite(self.free_slip[i] - dot(&self.compliance[i*p..(i+1)*p], &self.tractions)?)?;
                if let Some(residual) = &self.friction_residual {
                    slip = finite(slip - residual.displacement(i, &self.tractions)?)?;
                }
                if let Some(obstacles) = &self.obstacles { slip = finite(slip + obstacles.slip_shift(i)?)?; }
                if let Some(series) = self.series_drives[i] {
                    let eta = series.staged_stretch_m(self.tractions[i], self.sample_period_s());
                    slip = finite(slip - (eta - series.stretch_m()))?;
                }
                let drive = self.contact_loaded_drive(i, drives[i])?;
                let checked = if let Some(layer) = &self.shear_layers[i] {
                    layer.check(slip, self.tractions[i], drive.normal_force_n, self.sample_period_s(),
                        self.points[i].config, self.iterations[i])
                } else {
                    scalar::check(&self.points[i].port.law, self.points[i].config,
                        self.sample_period_s(), drive, self.tractions[i], slip, self.regimes[i], self.iterations[i]).map(|_| ())
                };
                if let Err(error) = checked {
                    if failure.is_none() { failure = Some(error); }
                }
            }
            match failure {
                None => { sweeps = sweep; break; },
                Some(cause) if sweep == self.config.max_sweeps => {
                    // The Gauss-Seidel sweep could not converge this step
                    // (frankensim-o7qbp). For a well-posed step the consistent
                    // stick/slip partition is unique; solve it directly with the
                    // active-set rescue. If it finds and verifies that partition,
                    // adopt it and accept; otherwise refuse on budget as before.
                    // Adopts the unique partition and continues, or returns a
                    // NAMED active-set error (exhausted with payload, or declined
                    // with the Gauss-Seidel error as its cause).
                    self.active_set_rescue(drives, gate, cause)?;
                    rescues = 1;
                    sweeps = sweep;
                    break;
                },
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
        let (obstacle_energy, obstacle_loss, obstacle_work, absolute_obstacle_work) = match &mut self.obstacles {
            Some(obstacles) => obstacles.stage(&self.core.network.candidates, self.core.network.dt, gate)?,
            None => (0.0, 0.0, 0.0, 0.0),
        };
        let residual_bound = self.obstacles.as_ref().is_some_and(obstacles::ObstacleSet::has_residual);
        let (residual_energy, residual_loss) = match &self.obstacles {
            Some(obstacles) if residual_bound => obstacles.staged_residual_j()?,
            _ => (0.0, 0.0),
        };
        let (mut loss, mut work, mut absolute_work) = (0.0, 0.0, 0.0);
        for i in 0..p {
            poll(gate)?;
            let distance = finite(drives[i].speed_m_s*self.sample_period_s())?;
            let mut slip = finite(distance-dot(&self.points[i].column, &self.core.network.free_delta)?)?;
            if let Some(residual) = &self.friction_residual {
                slip = finite(slip - residual.increment(i, &self.tractions)?)?;
            }
            let mut hair_loss = 0.0;
            if let Some(series) = self.series_drives[i] {
                let eta = series.staged_stretch_m(self.tractions[i], self.sample_period_s());
                slip = finite(slip - (eta - series.stretch_m()))?;
                hair_loss = series.dashpot_loss_j(eta, self.sample_period_s());
                self.staged_stretch[i] = eta;
            }
            let drive = self.contact_loaded_drive(i, drives[i])?;
            let dt = self.sample_period_s();
            let (dissipation, residual, shear_layer) = if let Some(layer) = &mut self.shear_layers[i] {
                let (frame, residual) = layer.stage(slip, self.tractions[i], drive.normal_force_n, dt,
                    self.points[i].config, self.iterations[i])?;
                self.regimes[i] = ModalFrictionRegime::Shearing { plastic_flow: frame.plastic_flow };
                (frame.dissipated_j, residual, Some(frame))
            } else {
                let (loss, residual) = scalar::check(&self.points[i].port.law, self.points[i].config,
                    dt, drive, self.tractions[i], slip, self.regimes[i], self.iterations[i])?;
                (loss, residual, None)
            };
            // A contact-loaded law has no known load until the JOINT solve.
            // Validate all three slip rates with that same held reaction before
            // any candidate, observer, normal offset or clock is published.
            if self.contact_loads[i].is_some() {
                self.check_interface_endpoint(i, &self.core.network.models, drive, false)?;
            }
            interface::check_average(self.interface(i), drive, slip, self.sample_period_s())?;
            self.check_interface_endpoint(i, &self.core.network.candidates, drive, true)?;
            let supply = finite(self.tractions[i]*distance)?;
            loss = finite(loss+dissipation+hair_loss)?; work = finite(work+supply)?;
            absolute_work = finite(absolute_work+supply.abs())?;
            self.staged_ports[i] = FrictionPortFrame { normal_force_n: drive.normal_force_n,
                traction_n: self.tractions[i], slip_distance_m: slip,
                drive_work_j: supply, dissipation_j: dissipation,
                series_drive_dissipation_j: hair_loss, constitutive_residual_n: residual,
                regime: self.regimes[i], roots: self.port_roots[i], shear_layer };
        }
        let mut pressure = self.core.staged_observer_pressure(gate)?;
        if let (Some(obstacles), Some(weights)) = (&self.obstacles, self.core.acceleration_pressure_weights()) {
            pressure = obstacles.endpoint_pressure(pressure, weights, gate)?;
            limit("acceleration observer pressure", pressure.abs(), self.core.network.config.maximum_abs_pressure_pa)?;
        }
        let network_energy = finite(self.core.network.staged_frame.total_energy_j())?;
        let hair_energy: f64 = (0..p).filter_map(|i|
            self.series_drives[i].map(|drive| drive.spring_energy_j(self.staged_stretch[i]))).sum();
        let after = if residual_bound { finite(network_energy+obstacle_energy+residual_energy+hair_energy)? }
            else if self.obstacles.is_some() { finite(network_energy+obstacle_energy+hair_energy)? }
            else { finite(network_energy+hair_energy)? };
        let shear_layer_energy = finite(self.shear_layers.iter().flatten()
            .map(layer::MaterialShearLayer::staged_energy_j).sum())?;
        let after = if self.shear_layers.iter().any(Option::is_some) { finite(after+shear_layer_energy)? }
            else { after };
        let (friction_residual_energy, friction_residual_loss) = match &self.friction_residual {
            Some(residual) => residual.staged_energy_loss_j(&self.tractions)?,
            None => (0.0, 0.0),
        };
        let after = if self.friction_residual.is_some() { finite(after+friction_residual_energy)? }
            else { after };
        let staged = &self.core.network.staged_frame;
        let network_loss = finite(staged.total_dissipation_j())?;
        let external_work = dot(external, &self.core.network.free_delta)?;
        let residual = if residual_bound {
            finite((after-before)+network_loss+loss+obstacle_loss+residual_loss-external_work-work-obstacle_work)?
        } else if self.obstacles.is_some() {
            finite((after-before)+network_loss+loss+obstacle_loss-external_work-work-obstacle_work)?
        } else { finite((after-before)+network_loss+loss-external_work-work)? };
        let residual = if self.friction_residual.is_some() { finite(residual+friction_residual_loss)? }
            else { residual };
        let scale = before.max(after).max(network_loss).max(loss).max(obstacle_loss)
            .max(external_work.abs()).max(absolute_work).max(absolute_obstacle_work);
        if self.obstacles.is_some() || self.shear_layers.iter().any(Option::is_some)
            || self.friction_residual.is_some() {
            limit("total energy including contact storage", after, self.core.network.config.maximum_total_energy_j)?;
        }
        let tolerance = finite(self.core.network.config.energy_absolute_tolerance_j
            +self.core.network.config.energy_relative_tolerance*scale)?;
        if loss < 0.0 || residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance }.into());
        }
        poll(gate)?;
        Ok(MultiFrictionFrame {
            sample: staged.sample, observer_pressure_pa: pressure,
            network_energy_j: network_energy, external_work_j: external_work,
            drive_work_j: work, network_dissipation_j: network_loss,
            friction_dissipation_j: loss, shear_layer_energy_j: shear_layer_energy, energy_residual_j: residual,
            obstacle_energy_j: obstacle_energy, obstacle_dissipation_j: obstacle_loss,
            obstacle_work_j: obstacle_work, obstacle_residual_energy_j: residual_energy,
            obstacle_residual_loss_j: residual_loss, energy_tolerance_j: tolerance, sweeps,
            friction_residual_energy_j: friction_residual_energy,
            friction_residual_loss_j: friction_residual_loss,
            active_set_rescues: rescues, obstacle_resolves: 0,
            nonlinear: None, ports: Vec::new(),
        })
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
        if self.friction_residual.is_some() {
            return Err(invalid("friction residual compliance does not reconstruct source-bound endpoint velocities").into());
        }
        let point = self.points.get(index).ok_or_else(|| invalid("interface names an unknown friction port"))?;
        if self.interface(index).is_some() || self.shear_layers[index].is_some()
            || point.port.law != *binding.card().friction_law() {
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

#[cfg(test)]
mod continuation_tests;
