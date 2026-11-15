//! The prepared linear impact description stepped as a declared graph on the
//! composition kernel (bead rc-1109 .6.4).
//!
//! [`ImpactGraph::new`] takes what [`LinearImpactSystem::new`] takes, the
//! same linear bodies, normal contact, volume springs and limits, and
//! declares it as a graph that fs-system steps, where the module compiles it
//! into its coupled owners. The parts:
//!
//! - A body that no volume touches is fs-nlmodal's `ModalStructure` over the
//!   owner model the module compiles for it: the same exact held-force step
//!   and the same physical damping rates.
//! - The bodies that volumes touch are exact owners, one per set of bodies
//!   the volumes join: fs-system's `ExactPortHamiltonianComponent` over their
//!   coordinates, interleaved `(q, p)` in the original order, with
//!   `H = p.p/2 + q.K q/2` and `K = diag(omega^2) + sum k c c^T`, one term
//!   per volume, the module's reference-area spring `k = bulk A (A / V)` on
//!   `c = areas / A`; each mode's physical rate acts on its momentum.
//! - The contact point is fs-dcontact's `ContactInterface` on the kernel's
//!   load line between the two parts its row touches, end `a` on the part
//!   with the lower coordinates. Each of the two takes the row's entries on
//!   its own coordinates through one port, `contact`, and publishes the end
//!   of its flow, which the interface reads. Its start closure is the row
//!   times the start displacements.
//! - Each part's heat goes to a sink outside the system at `heat_sink_k`.
//!
//! Against the module the coefficients are the same and the rounding
//! differs: the interface carries the closure as its own state, the kernel
//! forms each end's compliance from the part's tangent, and an exact owner
//! prepares its exponential from `H`, apart from the module's coupled block.
//! fs-couple's `tests/impact_graph.rs` and `tests/drum_strike_graph.rs`
//! bound the difference.
//!
//! A frame ([`ImpactGraphFrame`]) reads the kernel's ledger: the stored
//! energy is the parts' stored energy at the build plus the ledger's stored
//! change, a step's dissipation is the heat its sinks received, and its
//! residual is the change of the ledger's first-law residual.
//!
//! Refused by name: a nonlinear body, more than one contact point, a contact
//! row that does not join exactly two parts, a contact with internal loss
//! (the interface has no heat port), and a step the kernel shortened (the
//! clock is fixed). Not declared: viscous dampers
//! ([`LinearImpactSystem::new_with_dampers`]), external generalized forces
//! and the module's whole-system energy ceiling.
//!
//! [`LinearImpactSystem::new`]: super::linear::LinearImpactSystem::new
//! [`LinearImpactSystem::new_with_dampers`]: super::linear::LinearImpactSystem::new_with_dampers

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_nlmodal::component::{ModalAttachment, ModalStructure};
use fs_ports::component::{Component, PortRole, PortSpec};
use fs_ports::{PortKind, StableId};
use fs_system::phs_exact::{
    ExactPortHamiltonianComponent, LinearSplitBudget, PortHamiltonian, QuadraticStorage,
};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, SystemError, heat_sink};

use super::linear::{LinearImpactConfig, VolumeConnection};
use super::{BodyPotential, ImpactBody, ImpactError};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeModel};

/// The port through which a part takes the contact.
const CONTACT_PORT: &str = "contact";

/// The stream key of every step: no part draws random numbers.
const KEY: StreamKey = StreamKey {
    seed: 0,
    kernel_id: 0,
    tile: 0,
    iteration: 0,
};

fn invalid(what: &'static str) -> ImpactError {
    ImpactError::Invalid(what)
}

fn owner(error: impl core::fmt::Display) -> ImpactError {
    ImpactError::Owner(error.to_string())
}

/// For an identifier or sink refusal, whose error has no display form.
fn named(error: impl core::fmt::Debug) -> ImpactError {
    ImpactError::Owner(format!("{error:?}"))
}

fn kernel(error: SystemError) -> ImpactError {
    match error {
        SystemError::Cancelled => ImpactError::Cancelled,
        other => owner(other),
    }
}

/// One accepted step of an [`ImpactGraph`], read from the kernel's ledger.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImpactGraphFrame {
    /// One-based accepted mechanical sample; zero before a step.
    pub sample: u64,
    /// The kernel's time at the sample's end [s].
    pub time_s: f64,
    /// Energy stored in the parts [J]: their stored energy at the build plus
    /// the ledger's stored change.
    pub stored_energy_j: f64,
    /// Heat the sinks received over the step [J].
    pub dissipated_energy_j: f64,
    /// The step's change of the ledger's uncorrected first-law residual [J].
    pub balance_residual_j: f64,
    /// The kernel's iterations for the step.
    pub iterations: u32,
}

/// Which owner a part is.
enum Owner {
    Structure,
    Exact,
}

/// One part of the graph and the original coordinates it holds, increasing.
struct Part {
    node: NodeId,
    owner: Owner,
    coordinates: Vec<usize>,
}

/// The impact description on the composition kernel. See the module
/// documentation.
pub struct ImpactGraph {
    system: System,
    parts: Vec<Part>,
    contact: Option<NodeId>,
    state: Vec<f64>,
    period_s: f64,
    max_steps: u64,
    pool: ArenaPool,
    frame: ImpactGraphFrame,
    stored_at_build_j: f64,
    heat_j: f64,
    residual_j: f64,
}

impl ImpactGraph {
    /// Declare the bodies, contacts and volumes as a graph and seal it, the
    /// parts' heat going to sinks at `heat_sink_k` [K].
    ///
    /// # Errors
    /// [`ImpactError::Invalid`] for the description's own admission (as the
    /// module's) and for what this image refuses by name (the module
    /// documentation); [`ImpactError::Owner`] for any part's or the kernel's
    /// refusal, its message unchanged; [`ImpactError::Cancelled`].
    pub fn new(
        bodies: Vec<ImpactBody>,
        contacts: Vec<Obstacle>,
        volumes: Vec<VolumeConnection>,
        config: LinearImpactConfig,
        heat_sink_k: f64,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let counts: Vec<usize> = bodies.iter().map(|b| b.potential.count()).collect();
        let count = counts
            .iter()
            .try_fold(0usize, |n, &m| n.checked_add(m))
            .ok_or_else(|| invalid("linear impact mode count overflow"))?;
        if count == 0
            || counts.contains(&0)
            || count > config.coupling.max_modes
            || config.sample_rate_hz == 0
            || config.max_steps == 0
            || config.max_steps > (1u64 << 53)
            || volumes.len() > config.coupling.max_connections
            || !heat_sink_k.is_finite()
            || heat_sink_k <= 0.0
        {
            return Err(invalid(
                "the kernel impact image needs bounded bodies and volumes, a positive clock, an exact step budget and a positive sink temperature",
            ));
        }
        let mut omega = Vec::with_capacity(count);
        let mut rate = Vec::with_capacity(count);
        let mut start = Vec::with_capacity(count);
        let body_count = bodies.len();
        let mut first = Vec::with_capacity(body_count);
        for body in bodies {
            let BodyPotential::Linear(w) = &body.potential else {
                return Err(invalid(
                    "nonlinear shell storage cannot be declared on the kernel impact image",
                ));
            };
            if body.initial.len() != w.len()
                || body.damping_per_s.len() != w.len()
                || w.iter().any(|w| !w.is_finite() || *w < 0.0)
                || body
                    .damping_per_s
                    .iter()
                    .any(|d| !d.is_finite() || *d < 0.0)
            {
                return Err(invalid(
                    "linear impact body state/damping dimensions or coefficients are invalid",
                ));
            }
            first.push(omega.len());
            omega.extend_from_slice(w);
            rate.extend_from_slice(&body.damping_per_s);
            start.extend_from_slice(&body.initial);
        }
        let body_of = |coordinate: usize| first.partition_point(|&f| f <= coordinate) - 1;

        // The volume springs, and the sets of bodies they join.
        let mut set: Vec<usize> = (0..body_count).collect();
        let mut exact = vec![false; body_count];
        let mut springs = Vec::with_capacity(volumes.len());
        for volume in volumes {
            let v = volume.spring;
            let area = volume.reference_area_m2;
            if v.areas.len() != count
                || v.areas.iter().any(|a| !a.is_finite())
                || [area, v.bulk_modulus_pa, v.volume_m3]
                    .iter()
                    .any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(invalid(
                    "linear impact volume requires finite positive SI data and every signed modal area",
                ));
            }
            let stiffness = (v.bulk_modulus_pa * area) * (area / v.volume_m3);
            if !stiffness.is_finite() || stiffness <= 0.0 {
                return Err(invalid(
                    "volume reference-area stiffness is not representable",
                ));
            }
            let column: Vec<f64> = v.areas.iter().map(|b| b / area).collect();
            if column
                .iter()
                .zip(&v.areas)
                .any(|(&c, &a)| a != 0.0 && c == 0.0)
            {
                return Err(invalid("volume reference-area participation underflow"));
            }
            let mut touched: Vec<usize> = (0..count)
                .filter(|&c| column[c] != 0.0)
                .map(body_of)
                .collect();
            touched.dedup();
            if let Some(&joined) = touched.iter().map(|&b| &set[b]).min() {
                let old: Vec<usize> = touched.iter().map(|&b| set[b]).collect();
                for label in &mut set {
                    if old.contains(label) {
                        *label = joined;
                    }
                }
                for &b in &touched {
                    exact[b] = true;
                }
            }
            springs.push((stiffness, column));
        }

        // The parts, in the order of their first body.
        let mut parts: Vec<(usize, Vec<usize>, bool)> = Vec::new();
        for (b, &label) in set.iter().enumerate() {
            let range = first[b]..first[b] + counts[b];
            match parts.iter_mut().find(|(l, _, _)| *l == label) {
                Some((_, coordinates, _)) => coordinates.extend(range),
                None => parts.push((label, range.collect(), exact[b])),
            }
        }

        // The contact point: its row and its per-point law, as the module
        // admits them.
        let point_count = contacts.iter().map(Obstacle::n_points).sum::<usize>();
        if point_count > 1 {
            return Err(invalid(
                "the kernel impact image holds one contact point; distributed and several contacts are not declared",
            ));
        }
        let point = match contacts.into_iter().find(|ob| ob.n_points() > 0) {
            None => None,
            Some(ob) => {
                let ob = Obstacle::new(
                    ob.collocation().to_vec(),
                    ob.n_points(),
                    count,
                    ob.gaps().to_vec(),
                    ob.weights().to_vec(),
                    ob.stiffness(),
                    ob.alpha(),
                    ob.provenance().to_string(),
                )
                .and_then(|clean| clean.with_internal_loss(ob.internal_loss()))
                .map_err(owner)?;
                let column = ob.collocation()[..count].to_vec();
                let law = Obstacle::new(
                    vec![-1.0],
                    1,
                    1,
                    vec![ob.gaps()[0]],
                    vec![ob.weights()[0]],
                    ob.stiffness(),
                    ob.alpha(),
                    ob.provenance().to_string(),
                )
                .and_then(|law| law.with_internal_loss(ob.internal_loss()))
                .map_err(owner)?;
                let ends: Vec<usize> = (0..parts.len())
                    .filter(|&p| parts[p].1.iter().any(|&c| column[c] != 0.0))
                    .collect();
                if ends.len() != 2 {
                    return Err(invalid(
                        "a contact row must join exactly two parts of the kernel impact image",
                    ));
                }
                Some((column, law, [ends[0], ends[1]]))
            }
        };

        let rate_hz = config.sample_rate_hz;
        let mut builder = SystemBuilder::new();
        let mut stored_at_build_j = 0.0;
        let mut declared = Vec::with_capacity(parts.len());
        for (index, (_, coordinates, is_exact)) in parts.into_iter().enumerate() {
            if gate.is_requested() {
                return Err(ImpactError::Cancelled);
            }
            // The row's entries on this part, `-` at end `a` and `+` at end
            // `b`: a common flow gives each end minus its part's flow, so the
            // closure rate `flow_a - flow_b` is the row times the velocities.
            let weights = point.as_ref().and_then(|(column, _, ends)| {
                let sign = if ends[0] == index {
                    -1.0
                } else if ends[1] == index {
                    1.0
                } else {
                    return None;
                };
                Some(
                    coordinates
                        .iter()
                        .map(|&c| sign * column[c])
                        .collect::<Vec<f64>>(),
                )
            });
            let name = format!("part{index}");
            let (node, owner_kind) = if is_exact {
                let component = exact_owner(
                    &coordinates,
                    &omega,
                    &rate,
                    &springs,
                    weights.as_deref(),
                    &config,
                )?;
                let state: Vec<f64> = coordinates
                    .iter()
                    .flat_map(|&c| {
                        [
                            start[c].displacement_m_sqrt_kg,
                            start[c].velocity_m_sqrt_kg_per_s,
                        ]
                    })
                    .collect();
                stored_at_build_j += component.stored_energy_j(&state).map_err(owner)?;
                let node = builder
                    .add_with_tangent(&name, component, state)
                    .map_err(owner)?;
                (node, Owner::Exact)
            } else {
                let modes = coordinates
                    .iter()
                    .map(|&c| ModalAcousticMode {
                        angular_frequency_rad_s: omega[c],
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    })
                    .collect();
                let rates: Vec<f64> = coordinates.iter().map(|&c| rate[c]).collect();
                let model = ModalAcousticTimeModel::try_new_with_damping_rates(
                    rate_hz,
                    modes,
                    &rates,
                    config.component,
                )
                .map_err(owner)?;
                let attachments = match weights {
                    None => Vec::new(),
                    Some(weights) => vec![ModalAttachment {
                        id: StableId::new(CONTACT_PORT).map_err(named)?,
                        kind: PortKind::MechanicalForceVelocity,
                        width: NonZeroUsize::MIN,
                        weights,
                    }],
                };
                let structure = ModalStructure::new(model, attachments)
                    .and_then(ModalStructure::publishing_ends)
                    .map_err(owner)?;
                let state: Vec<ModalAcousticState> =
                    coordinates.iter().map(|&c| start[c]).collect();
                stored_at_build_j += structure.stored_energy_j(&state).map_err(owner)?;
                let node = builder
                    .add_with_tangent(&name, structure, state)
                    .map_err(owner)?;
                (node, Owner::Structure)
            };
            let heat = builder.port(node, "heat").map_err(owner)?;
            builder
                .heat_route(
                    heat,
                    &[],
                    Some(heat_sink(&format!("{name}_heat"), heat_sink_k).map_err(named)?),
                    RouteTemperature::Fixed(heat_sink_k),
                )
                .map_err(owner)?;
            declared.push(Part {
                node,
                owner: owner_kind,
                coordinates,
            });
        }

        let contact = match point {
            None => None,
            Some((column, law, ends)) => {
                let limits = config.contact;
                let solve = ContactSolve {
                    max_iterations: limits.max_iterations,
                    maximum_force_n: limits.maximum_force_n,
                    maximum_penetration_m: limits.maximum_penetration_m,
                    force_absolute_tolerance_n: limits.force_absolute_tolerance_n,
                    force_relative_tolerance: limits.force_relative_tolerance,
                };
                let interface = ContactInterface::new(law, solve)
                    .and_then(|interface| interface.reading_ends([true, true]))
                    .map_err(owner)?;
                let closure: f64 = column
                    .iter()
                    .zip(&start)
                    .map(|(c, s)| c * s.displacement_m_sqrt_kg)
                    .sum();
                let state = ContactState::at(closure);
                stored_at_build_j += interface.stored_energy_j(&state).map_err(owner)?;
                let node = builder
                    .add_on_line("contact", interface, state)
                    .map_err(owner)?;
                for (end, id) in [(ends[0], "a"), (ends[1], "b")] {
                    let part = builder
                        .port(declared[end].node, CONTACT_PORT)
                        .map_err(owner)?;
                    let interface_end = builder.port(node, id).map_err(owner)?;
                    builder.common_flow(part, &[interface_end]).map_err(owner)?;
                }
                Some(node)
            }
        };

        let system = builder.seal().map_err(owner)?;
        let state = start
            .iter()
            .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
            .collect();
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        Ok(Self {
            system,
            parts: declared,
            contact,
            state,
            period_s: f64::from(rate_hz).recip(),
            max_steps: config.max_steps,
            pool: ArenaPool::new(ArenaConfig::default()),
            frame: ImpactGraphFrame {
                stored_energy_j: stored_at_build_j,
                ..ImpactGraphFrame::default()
            },
            stored_at_build_j,
            heat_j: 0.0,
            residual_j: 0.0,
        })
    }

    /// Accepted interleaved `q`/`v` values in the original body/mode order,
    /// as [`LinearImpactSystem::state`](super::linear::LinearImpactSystem::state).
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.state
    }

    /// Number of mechanical coordinates.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.state.len() / 2
    }

    /// The fixed mechanical period [s].
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 {
        self.period_s
    }

    /// The last accepted frame, or the stored energy at the build before a
    /// step.
    #[must_use]
    pub const fn frame(&self) -> &ImpactGraphFrame {
        &self.frame
    }

    /// The sealed system, for its ledger and its parts' states.
    #[must_use]
    pub const fn system(&self) -> &System {
        &self.system
    }

    /// The contact's accepted state (closure and the last step's reaction),
    /// if the description has a contact.
    #[must_use]
    pub fn contact(&self) -> Option<&ContactState> {
        self.contact
            .and_then(|node| self.system.state::<ContactState>(node))
    }

    /// Step the system over one period of the fixed clock.
    ///
    /// # Errors
    /// [`ImpactError::Budget`] when the lifetime step budget is spent;
    /// [`ImpactError::Cancelled`]; [`ImpactError::Owner`] for the kernel's
    /// refusal; [`ImpactError::Invalid`] when the kernel shortened the step.
    pub fn step(&mut self, gate: &CancelGate) -> Result<ImpactGraphFrame, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if self.frame.sample >= self.max_steps {
            return Err(ImpactError::Budget);
        }
        let Self {
            system,
            pool,
            period_s,
            ..
        } = self;
        let period = *period_s;
        let report = pool
            .scope(|arena| {
                let cx = Cx::new(gate, arena, KEY, Budget::INFINITE, ExecMode::Deterministic);
                system.step(period, &cx)
            })
            .map_err(kernel)?;
        if report.dt_s.to_bits() != period.to_bits() {
            return Err(invalid(
                "the kernel shortened a step of the fixed mechanics clock",
            ));
        }
        for part in &self.parts {
            match part.owner {
                Owner::Structure => {
                    let modes = self
                        .system
                        .state::<Vec<ModalAcousticState>>(part.node)
                        .ok_or_else(|| invalid("a modal part's state is not modal"))?;
                    for (&c, s) in part.coordinates.iter().zip(modes) {
                        self.state[2 * c] = s.displacement_m_sqrt_kg;
                        self.state[2 * c + 1] = s.velocity_m_sqrt_kg_per_s;
                    }
                }
                Owner::Exact => {
                    let x = self
                        .system
                        .state::<Vec<f64>>(part.node)
                        .ok_or_else(|| invalid("an exact part's state is not a vector"))?;
                    for (&c, pair) in part.coordinates.iter().zip(x.as_chunks::<2>().0) {
                        self.state[2 * c] = pair[0];
                        self.state[2 * c + 1] = pair[1];
                    }
                }
            }
        }
        let ledger = self.system.ledger();
        let heat = ledger.sink_heat_j();
        let residual = ledger.energy_residual_j();
        self.frame = ImpactGraphFrame {
            sample: self.frame.sample + 1,
            time_s: report.time_s,
            stored_energy_j: self.stored_at_build_j + ledger.stored_change_j(),
            dissipated_energy_j: heat - self.heat_j,
            balance_residual_j: residual - self.residual_j,
            iterations: report.iterations,
        };
        self.heat_j = heat;
        self.residual_j = residual;
        Ok(self.frame)
    }
}

/// The exact owner of the bodies a set of volumes joins: `H` over the part's
/// interleaved `(q, p)`, the volumes' springs restricted to its coordinates
/// (a volume touches one part only), each mode's rate on its momentum, and
/// the contact's port if `weights` is given.
fn exact_owner(
    coordinates: &[usize],
    omega: &[f64],
    rate: &[f64],
    springs: &[(f64, Vec<f64>)],
    weights: Option<&[f64]>,
    config: &LinearImpactConfig,
) -> Result<ExactPortHamiltonianComponent, ImpactError> {
    let n = 2 * coordinates.len();
    let m = usize::from(weights.is_some());
    let mut storage = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n * m];
    for (mode, &c) in coordinates.iter().enumerate() {
        let (q, p) = (2 * mode, 2 * mode + 1);
        storage[q * n + q] = omega[c] * omega[c];
        storage[p * n + p] = 1.0;
        j[q * n + p] = 1.0;
        j[p * n + q] = -1.0;
        r[p * n + p] = rate[c];
        if let Some(weights) = weights {
            g[p] = weights[mode];
        }
    }
    for (k, column) in springs {
        for (a, &ca) in coordinates.iter().enumerate() {
            for (b, &cb) in coordinates.iter().enumerate() {
                storage[(2 * a) * n + 2 * b] += k * column[ca] * column[cb];
            }
        }
    }
    let storage = QuadraticStorage::new(storage, n).map_err(owner)?;
    let system = PortHamiltonian::new(n, m, j, r, g, Box::new(storage)).map_err(owner)?;
    let ports = if m == 1 {
        vec![
            PortSpec::power(
                StableId::new(CONTACT_PORT).map_err(named)?,
                PortKind::MechanicalForceVelocity,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )
            .map_err(owner)?,
        ]
    } else {
        Vec::new()
    };
    let budget = LinearSplitBudget {
        max_modes: 2 * config.coupling.max_modes,
        max_omega_dt: config.coupling.nyquist_guard_fraction * core::f64::consts::PI,
    };
    ExactPortHamiltonianComponent::new(system, ports, budget)
        .and_then(ExactPortHamiltonianComponent::publishing_ends)
        .map_err(owner)
}
