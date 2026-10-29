//! Cold exact mechanical propagation with full per-resistance loss integrals.
//!
//! z = [D q, v, memory, u] contains every original mass-normalized coordinate,
//! energy-normalized Maxwell internal strain and held generalized force.
//! D only balances units (sqrt(K_instantaneous_ii), or one
//! on free coordinates); it is invertible and removes no state. The constant
//! augmented generator is the same q'=v, v'=-Kq-Cv+u pencil as ExactBlock.
//! Short-interval Van Loan integrals avoid exp(-A^T dt) overflow for stiff
//! viscous memory: W(2h)=W(h)+F(h)^T W(h) F(h). The integral is linear in
//! the resistance, so the step's total connection loss is ONE integral of the
//! summed dashpots; each physical dashpot keeps its OWN integral for the port
//! ledger, computed when a ledger is first asked for. Maxwell arms use this same
//! exponential and their own positive loss Gramian, without mechanical masses.
//! This is a floating-point
//! matrix exponential, not a nonlinear/continuum accuracy certificate or a
//! real-time timing promise.
//!
//! A link anchored to a constant loss factor (`ConnectionLaw::constant_loss_reference_omega_rad_s`) wrote its dashpot as
//! `c = eta k / omega_ref`. Its resistance acts along the damping column `g(K) b` instead of `b`, with
//! `g(K) = V diag(s_a) V^T` over the coupled undamped modes of the complete stiffness (identity mass) and
//! `s_a = sqrt(omega_ref G(omega_a) / omega_a)`, G being the shell's five-branch unit spectrum
//! (`crate::shear_shell::constant_loss`). Coupled mode a then carries the loss factor `eta share_a G(omega_a)`, where
//! `share_a = k (b.v_a)^2 / omega_a^2` is the link's share of that mode's strain energy: the modal strain-energy rule
//! of the shell's constant loss, applied to one joint. The congruence keeps the resistance positive semidefinite. The
//! spring, its storage and its rest force keep `b`; the loss outside the spectrum's fit band is the fit's own shape.
//! With Maxwell arms, this connection policy uses the equilibrium mechanical
//! stiffness; the material branches retain their own causal frequency dependence.

use super::*;
use fs_phs::linear_exact::{DenseExponentialBudget, UpperBlocks, exponential_with_budget, upper_block_exponential_with_budget};
use crate::port_adoption::{DeferredTransfers, PortTransfer};
use crate::shear_shell::constant_loss::ConstantLossSpectrum;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// Explicit upper bound on complete mechanical modes in the cold dense path.
/// Each mode retains displacement, velocity and its independent held input.
pub const MAX_PREPARED_MODAL_MODES: usize = 1024;
// Each loss integral exponentiates Van Loan's 2d operator, d = 3 coordinates per mode; at the ceiling it must fit the
// budgeted exponential.
const _: () = assert!(6 * MAX_PREPARED_MODAL_MODES <= fs_phs::linear_exact::MAX_PREPARED_EXPONENTIAL_DIM);

/// Additional cold dense resources; original modal/connection/work and
/// physical/energy admissions remain binding. There is deliberately no default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PreparedModalBudget {
    /// Maximum TOTAL original mechanical modes, up to 1024; none are dropped.
    /// With M material memories and N modes, the independent exponential ceiling
    /// still requires 2*(3*N+M) <= 6144 before allocation.
    pub max_modes: usize,
    /// Conservative peak bytes of additional dense preparation/cache storage.
    /// Includes up to 48 scale levels (kept for the per-link integrals), the component and total connection loss matrices
    /// and every per-link loss matrix the port ledger may compute.
    pub max_dense_bytes: usize,
    /// Conservative scalar multiply-add term ceiling, including every shared
    /// exponential's 32 Taylor products and at most 48 squarings.
    pub max_multiply_terms: u64,
}

/// One exact step before any energy accounting: the held input (with the air's held port force, if any), the propagated
/// end state and its finite increments.
struct Propagation {
    x: Vec<f64>,
    u: Vec<f64>,
    z: Vec<f64>,
    endpoint: Vec<f64>,
    y: Vec<f64>,
    supplied: f64,
    /// The air's held force on each source port.
    air_force: Option<Vec<f64>>,
}

/// The loss quadratics z^T W z of ONE state/input vector (component, connection and material Gramians),
/// and the exact bits of that vector. One stage evaluates the identical vector in `advance` and in `losses`; they share this
/// sweep instead of each streaming the dense Gramians again. The values are the very ones a fresh evaluation returns.
struct LossSweep {
    z: Vec<f64>,
    component: f64,
    connection: f64,
    material: f64,
}

/// What the per-link port ledger needs and every step shares: the little of the layout that turns a state into z, and each
/// link's dense Gramian (None without damping). The Gramians are computed when a ledger is first asked for, and the sweep
/// over all of them runs only then.
struct LedgerData {
    count: usize,
    memory_count: usize,
    dt: f64,
    weights: Vec<f64>,
    columns: Vec<Vec<f64>>,
    /// (stiffness [N/m], rest extension [m]) of every link.
    springs: Vec<(f64, f64)>,
    /// Damping [N s/m] of every link.
    dampings: Vec<f64>,
    /// The damping column `g(K) b` of every constant-loss link; None where the dashpot acts along the link's own column.
    damping_columns: Vec<Option<Vec<f64>>>,
    /// What a link's Gramian is computed from; None without a dissipative link.
    integral: Option<LinkIntegral>,
    gramians: OnceLock<Vec<Option<Vec<f64>>>>,
    evaluations: AtomicU64,
}

/// The generator, every scale level (the last is the transition) and the short interval of the network: one link's loss
/// Gramian is `gramian` of these and its own resistance, bit for bit the integral the construction would have computed.
struct LinkIntegral {
    generator: Vec<f64>,
    levels: Vec<Vec<f64>>,
    h: f64,
    budget: DenseExponentialBudget,
}

/// An accepted step's per-link port ledger, kept as its inputs until an observer asks (`DeferredTransfers`): the start state,
/// the held input, the end state and the step's total connection loss.
struct PendingLedger {
    data: Arc<LedgerData>,
    x: Vec<f64>,
    u: Vec<f64>,
    endpoint: Vec<f64>,
    connection_loss: f64,
}

pub(super) struct PreparedNetwork {
    count: usize,
    dt: f64,
    weights: Vec<f64>,
    stiffness: Vec<f64>,
    rest_force: Vec<f64>,
    transition: Vec<f64>,
    component_loss: Option<Vec<f64>>,
    /// Each row c defines the original arm storage (c.q - memory)^2 / 2.
    material_couplings: Vec<Vec<f64>>,
    material_loss: Option<Vec<f64>>,
    /// The Gramian of the summed link resistance: z^T W z of it is the step's total connection loss. None without a
    /// dissipative link. The per-link Gramians are the port ledger's (`ledger`).
    connection_total: Option<Vec<f64>>,
    ledger: Arc<LedgerData>,
    sweep: Mutex<Option<LossSweep>>,
    sweeps: AtomicU64,
    /// Two-way exterior air beside the transition (`air`); None for the mechanics alone.
    air: Option<air::PreparedAir>,
}

const MATRIX_TOL: f64 = 1e-11;
const MAX_SCALES: usize = 48;

impl PreparedModalBudget {
    /// Validate fixed capacity fields without allocating or choosing a backend.
    pub fn validate(self) -> Result<(), ModalCouplingError> {
        if self.max_modes == 0 || self.max_modes > MAX_PREPARED_MODAL_MODES
            || self.max_dense_bytes == 0 || self.max_multiply_terms == 0 {
            return Err(invalid("prepared exact modal byte/work/mode budget"));
        }
        Ok(())
    }

    /// `air` is (ports, poles) of a two-way air bank, admitted with its own entries and setup work.
    fn admit(self, count: usize, losses: usize, links: usize, anchored: usize, air: Option<(usize, usize)>)
        -> Result<DenseExponentialBudget, ModalCouplingError>
    {
        self.admit_with_memory(count, 0, losses, links, anchored, air)
    }

    fn admit_with_memory(self, count: usize, memory: usize, losses: usize, links: usize,
        anchored: usize, air: Option<(usize, usize)>)
        -> Result<DenseExponentialBudget, ModalCouplingError>
    {
        self.validate()?;
        if count == 0 || count > self.max_modes {
            return Err(invalid("complete prepared modal count exceeds max_modes; no modes were dropped"));
        }
        // The hard mode limit makes these small, but every resource product
        // remains checked before the first prepared dense allocation.
        let d = count.checked_mul(3).and_then(|n| n.checked_add(memory))
            .ok_or_else(|| invalid("prepared state extent"))?;
        if d > fs_phs::linear_exact::MAX_PREPARED_EXPONENTIAL_DIM / 2 {
            return Err(invalid("complete prepared mechanical/material state exceeds exponential capacity; no states were dropped"));
        }
        let extent = d.checked_mul(d).ok_or_else(|| invalid("prepared matrix extent"))?;
        let matrices = losses.checked_add(128).ok_or_else(|| invalid("prepared workspace count"))?;
        // A constant-loss link adds one coupled eigensolve of the n x n stiffness (its copy, the eigenvectors and the
        // solver's own two squares) and one damping column.
        let eigensolve_entries = if anchored == 0 { 0 } else { 4 * count * count + anchored * count };
        let (air_entries, air_terms) = match air {
            None => (Some(0), Some(0)),
            Some((ports, poles)) => (air::PreparedAir::dense_entries(count, memory, ports, poles),
                air::PreparedAir::setup_terms(ports, poles)),
        };
        let bytes = extent.checked_mul(matrices)
            .and_then(|n| n.checked_add(eigensolve_entries))
            .and_then(|n| n.checked_add(air_entries?))
            .and_then(|n| n.checked_add(count.checked_mul(memory)?.checked_mul(2)?))
            .and_then(|n| n.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(links.checked_mul(core::mem::size_of::<Option<Vec<f64>>>())?))
            .ok_or_else(|| invalid("prepared dense byte overflow"))?;
        let cube = (d as u64).checked_pow(3).ok_or_else(|| invalid("prepared dense work overflow"))?;
        // A 2d Van Loan exponential per nonzero channel; two d transition
        // exponentials for the independent h versus h/2 check; Gramian
        // doubling, energy-identity validation and all physical assembly.
        let products = (losses as u64).checked_mul(80 * 8 + 2 * MAX_SCALES as u64 + 2)
            .and_then(|n| n.checked_add(2 * 80 + 2 * MAX_SCALES as u64 + 16))
            .ok_or_else(|| invalid("prepared product work overflow"))?;
        // Householder tridiagonalization and QL with vectors are about 9 n^3; each damping column is two n^2 products.
        let eigensolve_terms = if anchored == 0 { Some(0) } else {
            (count as u64).checked_pow(3).and_then(|n| n.checked_mul(9))
                .and_then(|n| n.checked_add((anchored as u64).checked_mul((count as u64).checked_pow(2)?)?.checked_mul(2)?))
        };
        let terms = cube.checked_mul(products)
            .and_then(|n| n.checked_add((links as u64).checked_mul((count as u64).checked_pow(2)?)?.checked_mul(8)?))
            .and_then(|n| n.checked_add((memory as u64).checked_mul((count as u64).checked_pow(2)?)?.checked_mul(16)?))
            .and_then(|n| n.checked_add(eigensolve_terms?))
            .and_then(|n| n.checked_add(air_terms?))
            .ok_or_else(|| invalid("prepared complete work overflow"))?;
        if bytes > self.max_dense_bytes {
            return Err(invalid("prepared exact dense bytes exceed max_dense_bytes"));
        }
        if terms > self.max_multiply_terms {
            return Err(invalid("prepared exact multiply terms exceed max_multiply_terms"));
        }
        Ok(DenseExponentialBudget {
            max_dimension: 2 * d,
            max_dense_bytes: self.max_dense_bytes,
            max_multiply_terms: self.max_multiply_terms,
        })
    }
}

impl PreparedNetwork {
    pub(super) fn new(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        count: usize,
        dt: f64,
        budget: PreparedModalBudget,
        branches: &[fs_phs::RelaxationBranch],
        radiation: Option<&radiation::Binding>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if !dt.is_finite() || dt <= 0.0 || columns.len() != links.len()
            || columns.iter().any(|c| c.len() != count)
            || models.iter().map(|m| m.modes().len()).sum::<usize>() != count {
            return Err(invalid("prepared complete mechanical layout/clock"));
        }
        let air_extent = radiation.map(|binding| (binding.model.ports, binding.model.poles.len()));
        let has_component_loss = models.iter().flat_map(|m| m.modes())
            .any(|m| m.damping_ratio != 0.0 && m.angular_frequency_rad_s != 0.0);
        // One integral for the component loss, one for the summed connection loss, and one per dissipative link that an
        // observed port ledger computes: all are admitted here, before anything is allocated.
        let dissipative = links.iter().filter(|c| c.damping_n_s_m != 0.0).count();
        let anchored = links.iter().filter(|c| c.constant_loss_reference_omega_rad_s.is_some()).count();
        let losses = usize::from(has_component_loss) + usize::from(dissipative > 0) + dissipative
            + usize::from(!branches.is_empty());
        let exp_budget = if branches.is_empty() {
            budget.admit(count, losses, links.len(), anchored, air_extent)?
        } else {
            budget.admit_with_memory(count, branches.len(), losses, links.len(), anchored, air_extent)?
        };
        let state_count = 2 * count + branches.len();
        let d = state_count + count;
        let mut material_couplings = Vec::with_capacity(branches.len());
        for branch in branches {
            poll(Some(gate))?;
            if branch.projection.len() != 2 * count || !branch.stiffness.is_finite()
                || branch.stiffness <= 0.0 || !branch.relaxation_time_s.is_finite()
                || branch.relaxation_time_s <= 0.0 || !(1.0 / branch.relaxation_time_s).is_finite()
                || branch.projection.iter().any(|v| !v.is_finite())
                || branch.projection.iter().skip(1).step_by(2).any(|v| *v != 0.0) {
                return Err(invalid("prepared material requires complete displacement-only passive Maxwell arms"));
            }
            let root = det::sqrt(branch.stiffness);
            material_couplings.push(branch.projection.iter().step_by(2)
                .map(|b| finite(root * b)).collect::<Result<Vec<_>, _>>()?);
        }
        let mut stiffness = vec![0.0; count * count];
        let mut damping = vec![0.0; count * count];
        let mut component_rates = vec![0.0; count];
        let mut rest_force = vec![0.0; count];
        for (i, mode) in models.iter().flat_map(|m| m.modes()).enumerate() {
            stiffness[i * count + i] = finite(mode.angular_frequency_rad_s.powi(2))?;
            component_rates[i] = finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s)?;
            damping[i * count + i] = component_rates[i];
        }
        for (column, link) in columns.iter().zip(links) {
            poll(Some(gate))?;
            // A constant-loss dashpot needs the complete stiffness first; it is added below along its damping column.
            let viscous = link.constant_loss_reference_omega_rad_s.is_none();
            for i in 0..count {
                rest_force[i] = finite(rest_force[i] + link.stiffness_n_m * link.rest_extension_m * column[i])?;
                for j in 0..count {
                    stiffness[i * count + j] = finite(stiffness[i * count + j] + link.stiffness_n_m * column[i] * column[j])?;
                    if viscous {
                        damping[i * count + j] = finite(damping[i * count + j] + link.damping_n_s_m * column[i] * column[j])?;
                    }
                }
            }
        }
        let damping_columns = constant_loss_columns(&stiffness, columns, links, count, gate)?;
        for (own, link) in damping_columns.iter().zip(links) {
            let Some(column) = own else { continue };
            poll(Some(gate))?;
            for i in 0..count {
                for j in 0..count {
                    damping[i * count + j] = finite(damping[i * count + j] + link.damping_n_s_m * column[i] * column[j])?;
                }
            }
        }
        let weights: Vec<_> = (0..count).map(|i| {
            let k = material_couplings.iter().try_fold(stiffness[i * count + i],
                |sum, row| finite(sum + row[i] * row[i]))?;
            Ok(if k > 0.0 { det::sqrt(k) } else { 1.0 })
        }).collect::<Result<_, ModalCouplingError>>()?;
        let mut generator = vec![0.0; d * d];
        for i in 0..count {
            generator[(2 * i) * d + 2 * i + 1] = weights[i];
            generator[(2 * i + 1) * d + state_count + i] = 1.0;
            for j in 0..count {
                generator[(2 * i + 1) * d + 2 * j] = finite(-stiffness[i * count + j] / weights[j])?;
                generator[(2 * i + 1) * d + 2 * j + 1] = -damping[i * count + j];
            }
        }
        let mut material_resistance = (!branches.is_empty()).then(|| vec![0.0; d * d]);
        for (arm, (coupling, branch)) in material_couplings.iter().zip(branches).enumerate() {
            poll(Some(gate))?;
            let memory = 2 * count + arm;
            let rate = 1.0 / branch.relaxation_time_s;
            let scaled: Vec<_> = coupling.iter().zip(&weights)
                .map(|(c, weight)| finite(c / weight)).collect::<Result<_, _>>()?;
            // z' = (c.q-z)/tau; v' receives -c*(c.q-z). Every signed
            // cross-modal entry remains in this original physical generator.
            generator[memory * d + memory] = -rate;
            let resistance = material_resistance.as_mut().unwrap();
            resistance[memory * d + memory] = rate;
            for i in 0..count {
                generator[(2 * i + 1) * d + memory] = coupling[i];
                generator[memory * d + 2 * i] = finite(rate * scaled[i])?;
                let cross = finite(-rate * scaled[i])?;
                resistance[(2 * i) * d + memory] = cross;
                resistance[memory * d + 2 * i] = cross;
                for j in 0..count {
                    let entry = &mut generator[(2 * i + 1) * d + 2 * j];
                    *entry = finite(*entry - coupling[i] * scaled[j])?;
                    let entry = &mut resistance[(2 * i) * d + 2 * j];
                    *entry = finite(*entry + rate * (scaled[i] * scaled[j]))?;
                }
            }
        }
        let (levels, h) = scaled_levels(&generator, d, dt, exp_budget, gate)?;

        let mut component_loss = None;
        if has_component_loss {
            let mut resistance = vec![0.0; d * d];
            for (i, c) in component_rates.iter().enumerate() {
                resistance[(2 * i + 1) * d + 2 * i + 1] = *c;
            }
            component_loss = Some(gramian(&generator, &resistance, &levels, h, d, exp_budget, gate)?);
        }
        // The Gramian is linear in the resistance: ONE integral of the summed dashpots is the total connection loss, however
        // many links there are. Each link's own integral is left to the port ledger.
        let mut connection_resistance: Option<Vec<f64>> = None;
        for ((column, link), damping_column) in columns.iter().zip(links).zip(&damping_columns) {
            poll(Some(gate))?;
            if link.damping_n_s_m == 0.0 {
                continue;
            }
            let own = link_resistance(damping_column.as_deref().unwrap_or(column), link.damping_n_s_m, count, d)?;
            match connection_resistance.as_mut() {
                Some(sum) => for (s, v) in sum.iter_mut().zip(&own) { *s = finite(*s + v)?; },
                None => connection_resistance = Some(own),
            }
        }
        let connection_total = match &connection_resistance {
            Some(resistance) => Some(gramian(&generator, resistance, &levels, h, d, exp_budget, gate)?),
            None => None,
        };
        let material_loss = match &material_resistance {
            Some(resistance) => Some(gramian(&generator, resistance, &levels, h, d, exp_budget, gate)?),
            None => None,
        };
        let transition = levels.last().unwrap().clone();
        let air = radiation.map(|binding| air::PreparedAir::new(binding, count, branches.len(), dt, &weights, &transition, gate))
            .transpose()?;
        let integral = connection_resistance.map(|_| LinkIntegral { generator, levels, h, budget: exp_budget });
        let ledger = Arc::new(LedgerData {
            count, memory_count: branches.len(), dt, weights: weights.clone(), columns: columns.to_vec(),
            springs: links.iter().map(|l| (l.stiffness_n_m, l.rest_extension_m)).collect(),
            dampings: links.iter().map(|l| l.damping_n_s_m).collect(), damping_columns,
            integral, gramians: OnceLock::new(), evaluations: AtomicU64::new(0),
        });
        let result = Self { count, dt, weights, stiffness, rest_force, transition, component_loss, connection_total,
            material_couplings, material_loss, ledger,
            sweep: Mutex::new(None), sweeps: AtomicU64::new(0), air };
        result.validate_energy_identity(gate)?;
        poll(Some(gate))?;
        Ok(result)
    }

    fn state_input(&self, x: &[f64], u: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        state_vector(self.count, self.material_couplings.len(), &self.weights, self.dt, x, u, dt)
    }

    /// Original mechanical storage plus every Maxwell arm's stored strain energy.
    pub(super) fn hamiltonian(&self, x: &[f64]) -> f64 {
        let n = self.count;
        let mut energy = 0.0;
        for i in 0..n {
            energy += 0.5 * x[2 * i + 1].powi(2);
            for j in 0..n {
                energy += 0.5 * x[2 * i] * self.stiffness[i * n + j] * x[2 * j];
            }
        }
        for (arm, row) in self.material_couplings.iter().enumerate() {
            let mismatch = row.iter().enumerate().fold(-x[2 * n + arm],
                |sum, (i, c)| sum + c * x[2 * i]);
            energy += 0.5 * mismatch * mismatch;
        }
        energy
    }

    /// Storage of the complete state: mechanics, material history, then air pairs.
    pub(super) fn total_hamiltonian(&self, x: &[f64]) -> Result<f64, ModalCouplingError> {
        let mechanical = self.hamiltonian(x);
        match &self.air {
            None => Ok(mechanical),
            Some(air) => finite(mechanical + air.energy(x.get(2 * self.count + self.material_couplings.len()..).unwrap_or_default())?),
        }
    }

    /// The acoustic state entries a complete state carries after mechanics and material history.
    pub(super) fn air_state_len(&self) -> usize { self.air.as_ref().map_or(0, air::PreparedAir::state_len) }

    fn propagate(&self, models: &[ModalAcousticTimeModel], memory: &[f64], radiation: &[f64], external: &[f64], dt: f64)
        -> Result<Propagation, ModalCouplingError>
    {
        if external.len() != self.count || memory.len() != self.material_couplings.len()
            || radiation.len() != self.air_state_len() {
            return Err(invalid("prepared external force layout"));
        }
        let mut x = flatten(models);
        x.extend_from_slice(memory);
        let mut u: Vec<_> = external.iter().zip(&self.rest_force)
            .map(|(f, r)| finite(f + r)).collect::<Result<_, _>>()?;
        let mut z = self.state_input(&x, &u, dt)?;
        // The air's held force solves the step's interconnection against this held input (see `air`).
        let air_force = match &self.air {
            None => None,
            Some(air) => {
                let force = air.exchange_force(&air.port_increments(&z, &x, &self.weights)?, &air.free_increments(radiation)?)?;
                for (&mode, f) in air.sources().iter().zip(&force) {
                    u[mode] = finite(u[mode] + f)?;
                }
                z = self.state_input(&x, &u, dt)?;
                Some(force)
            }
        };
        let end = apply(&self.transition, &z)?;
        let mut endpoint = vec![0.0; 2 * self.count + memory.len()];
        let mut y = vec![0.0; self.count];
        let mut supplied = 0.0;
        for i in 0..self.count {
            endpoint[2 * i] = finite(end[2 * i] / self.weights[i])?;
            endpoint[2 * i + 1] = end[2 * i + 1];
            let dq = finite(endpoint[2 * i] - x[2 * i])?;
            y[i] = finite(dq / dt)?;
            supplied = finite(supplied + u[i] * dq)?;
        }
        endpoint[2 * self.count..].copy_from_slice(&end[2 * self.count..2 * self.count + memory.len()]);
        Ok(Propagation { x, u, z, endpoint, y, supplied, air_force })
    }

    /// The end state of the free (contact-free) step, without the loss Gramians and the trajectory energy-residual check.
    /// Every finiteness check of the propagation still runs. With air the end state is the air-loaded one. Returns
    /// (mechanical end state, mechanical start state).
    pub(super) fn predict_free(&self, models: &[ModalAcousticTimeModel], memory: &[f64], radiation: &[f64], external: &[f64], dt: f64)
        -> Result<(Vec<f64>, Vec<f64>), ModalCouplingError>
    {
        let step = self.propagate(models, memory, radiation, external, dt)?;
        Ok((step.endpoint, step.x))
    }

    /// One exact step. With air the record's state and the returned start state carry the air's (q, p) pairs after the
    /// mechanics, the returned held input includes the air's held port force, and the record's storage change, loss and
    /// supplied work are those of mechanics and air together (the interconnection is internal).
    pub(super) fn advance(&self, models: &[ModalAcousticTimeModel], memory: &[f64], radiation: &[f64], external: &[f64], dt: f64)
        -> Result<(StepRecord, Vec<f64>, Vec<f64>), ModalCouplingError>
    {
        let Propagation { x, u, z, endpoint, y, supplied, air_force } = self.propagate(models, memory, radiation, external, dt)?;
        let (component, connection, material) = self.loss_sums(&z)?;
        let dissipated = finite(component + connection + material)?;
        let before = finite(self.hamiltonian(&x))?;
        let after = finite(self.hamiltonian(&endpoint))?;
        let delta_h = finite(after - before)?;
        let scale = before.abs().max(after.abs()).max(supplied.abs()).max(dissipated).max(f64::MIN_POSITIVE);
        let residual = finite(delta_h + dissipated - supplied)?;
        if residual.abs() > MATRIX_TOL * scale {
            return Err(invalid("prepared exact trajectory energy residual"));
        }
        let (Some(air), Some(air_force)) = (&self.air, air_force) else {
            return Ok((StepRecord {
                x: endpoint, y, delta_h, dissipated, supplied,
                newton_iters: 0, solver_residual: residual.abs(),
                solver_certificate: None,
            }, x, u));
        };
        let dq_ports = air.port_increments(&z, &x, &self.weights)?;
        let step = air.advance(radiation, &dq_ports)?;
        // The held force the air's own step implies must be the one the mechanics were held at, and the work each side
        // received must cancel. Both forces are formed from increments, end less start, so each port's defect is held
        // beyond the rounding bound of those differences (see `force_rounding`).
        let exchange = air_force.iter().zip(&dq_ports).try_fold(0.0, |s, (f, dq)| finite(s + f * dq))?;
        let force_scale = air_force.iter().chain(&step.force).fold(0.0_f64, |m, f| m.max(f.abs())).max(f64::MIN_POSITIVE);
        let rounding = air.force_rounding(&z, &x, &self.weights, radiation, &dq_ports)?;
        let force_defect = air_force.iter().zip(&step.force).zip(&rounding)
            .fold(0.0_f64, |m, ((a, b), r)| m.max((a - b).abs() - r));
        let energy_scale = scale.max(exchange.abs()).max(step.work.abs()).max(step.loss);
        if force_defect > 1e-9 * force_scale || (exchange + step.work).abs() > MATRIX_TOL * energy_scale {
            return Err(invalid("prepared air exchange residual"));
        }
        let mut record_x = endpoint;
        record_x.extend_from_slice(&step.state);
        let mut start = x;
        start.extend_from_slice(radiation);
        let supplied = finite(supplied - exchange)?;
        let delta_h = finite(delta_h + step.delta_h)?;
        let dissipated = finite(dissipated + step.loss)?;
        let residual = finite(delta_h + dissipated - supplied)?;
        Ok((StepRecord {
            x: record_x, y, delta_h, dissipated, supplied,
            newton_iters: 0, solver_residual: residual.abs(),
            solver_certificate: None,
        }, start, u))
    }

    /// Every mode's end displacement per unit held force along `column`: the air-loaded response when the network has
    /// air, at rest.
    pub(super) fn response(&self, column: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        let z = self.state_input(&vec![0.0; 2 * self.count + self.material_couplings.len()], column, dt)?;
        let end = apply(&self.transition, &z)?;
        let mut response = (0..self.count).map(|i| finite(end[2 * i] / self.weights[i])).collect::<Result<Vec<_>, _>>()?;
        if let Some(air) = &self.air {
            air.load_response(&mut response)?;
        }
        Ok(response)
    }

    /// Original component and pre-summed connection loss, plus material heat.
    /// The number of dense forms is independent of the number of links/arms;
    /// each link's own form is evaluated only when its ledger is requested.
    fn compute_sweep(&self, z: &[f64]) -> Result<LossSweep, ModalCouplingError> {
        self.sweeps.fetch_add(1, Ordering::Relaxed);
        let (component, connection) = match (&self.component_loss, &self.connection_total) {
            (Some(c), Some(k)) => quadratic_pair(c, k, z)?,
            (Some(c), None) => (quadratic(c, z)?, 0.0),
            (None, Some(k)) => (0.0, quadratic(k, z)?),
            (None, None) => (0.0, 0.0),
        };
        let material = match &self.material_loss {
            Some(w) => quadratic(w, z)?,
            None => 0.0,
        };
        Ok(LossSweep { z: z.to_vec(), component, connection, material })
    }

    /// Run `read` on the sweep for exactly this vector (bitwise), evaluating the Gramians only when it is a new vector.
    fn with_sweep<T>(&self, z: &[f64], read: impl FnOnce(&LossSweep) -> Result<T, ModalCouplingError>)
        -> Result<T, ModalCouplingError>
    {
        let mut slot = self.sweep.lock().map_err(|_| invalid("prepared loss sweep lock poisoned"))?;
        let same = slot.as_ref().is_some_and(|s| s.z.len() == z.len()
            && s.z.iter().zip(z).all(|(a, b)| a.to_bits() == b.to_bits()));
        if !same { *slot = Some(self.compute_sweep(z)?); }
        match slot.as_ref() {
            Some(sweep) => read(sweep),
            None => Err(invalid("prepared loss sweep unavailable")),
        }
    }

    fn loss_sums(&self, z: &[f64]) -> Result<(f64, f64, f64), ModalCouplingError> {
        self.with_sweep(z, |sweep| Ok((sweep.component, sweep.connection, sweep.material)))
    }

    /// The damping column of a constant-loss link (`g(K) b`), or None where the link's dashpot acts along its own column.
    pub(super) fn damping_column(&self, link: usize) -> Option<&[f64]> {
        self.ledger.damping_columns.get(link).and_then(Option::as_deref)
    }

    #[cfg(test)]
    pub(super) fn sweeps(&self) -> u64 { self.sweeps.load(Ordering::Relaxed) }

    #[cfg(test)]
    pub(super) fn ledger_evaluations(&self) -> u64 { self.ledger.evaluations.load(Ordering::Relaxed) }

    /// (component, connection, material, radiated) losses of the step from the complete start state `x` under the held input `u`,
    /// both as `advance` returned them.
    pub(super) fn losses(&self, x: &[f64], u: &[f64], dt: f64, total: f64)
        -> Result<(f64, f64, f64, f64), ModalCouplingError>
    {
        let state_count = 2 * self.count + self.material_couplings.len();
        let mechanical = x.get(..state_count).ok_or_else(|| invalid("prepared complete state layout"))?;
        let z = self.state_input(mechanical, u, dt)?;
        let (component, connection, material) = self.loss_sums(&z)?;
        // The radiated loss of the same step: the transition's own port rows on the same z give the bits `advance` drove
        // the air with.
        let radiated = match &self.air {
            None => 0.0,
            Some(air) => air.loss(&x[state_count..], &air.port_increments(&z, mechanical, &self.weights)?)?,
        };
        let sum = finite(component + connection + material + radiated)?;
        if !total.is_finite() || (sum - total).abs() > MATRIX_TOL * sum.abs().max(total.abs()).max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared independent loss channels disagree with total"));
        }
        Ok((component, connection, material, radiated))
    }

    /// The per-port work of an accepted step, evaluated now: the entry of every caller that is not deferred. Only the
    /// mechanical/material trajectory determines each link's integrated loss.
    pub(super) fn connection_work(
        &self, x: &[f64], u: &[f64], endpoint: &[f64], dt: f64,
        columns: &[Vec<f64>], links: &[ConnectionLaw], connection_loss: f64,
        gate: Option<&CancelGate>,
    ) -> Result<Vec<PortTransfer>, ModalCouplingError> {
        let springs: Vec<_> = links.iter().map(|l| (l.stiffness_n_m, l.rest_extension_m)).collect();
        let n = 2 * self.count + self.material_couplings.len();
        let (Some(x), Some(endpoint)) = (x.get(..n), endpoint.get(..n)) else {
            return Err(invalid("prepared connection observation layout"));
        };
        self.ledger.work(x, u, endpoint, dt, connection_loss, columns, &springs, gate)
    }

    /// The same ledger as a deferred observation: nothing is evaluated until an observer asks for the window's transfers, and
    /// the sweep over every link's Gramian then runs once. The step's mechanical inputs are copied; the Gramians are shared.
    pub(super) fn pending_ledger(&self, x: &[f64], u: &[f64], endpoint: &[f64], connection_loss: f64)
        -> Arc<dyn DeferredTransfers>
    {
        let n = (2 * self.count + self.material_couplings.len()).min(x.len()).min(endpoint.len());
        Arc::new(PendingLedger {
            data: Arc::clone(&self.ledger), x: x[..n].to_vec(), u: u.to_vec(), endpoint: endpoint[..n].to_vec(),
            connection_loss,
        })
    }

    fn validate_energy_identity(&self, gate: &CancelGate) -> Result<(), ModalCouplingError> {
        let n = self.count;
        let state_count = 2 * n + self.material_couplings.len();
        let d = state_count + n;
        let mut storage = vec![0.0; d * d];
        for i in 0..n {
            storage[(2 * i + 1) * d + 2 * i + 1] = 1.0;
            for j in 0..n {
                storage[(2 * i) * d + 2 * j] = finite(self.stiffness[i * n + j] / self.weights[i] / self.weights[j])?;
            }
        }
        for (arm, row) in self.material_couplings.iter().enumerate() {
            poll(Some(gate))?;
            let memory = 2 * n + arm;
            storage[memory * d + memory] = 1.0;
            for i in 0..n {
                let a = finite(row[i] / self.weights[i])?;
                storage[(2 * i) * d + memory] = -a;
                storage[memory * d + 2 * i] = -a;
                for j in 0..n {
                    let entry = &mut storage[(2 * i) * d + 2 * j];
                    *entry = finite(*entry + a * (row[j] / self.weights[j]))?;
                }
            }
        }
        let propagated = sandwich(&self.transition, &storage, d, gate)?;
        let mut defect = vec![0.0; d * d];
        let mut scale = vec![0.0; d * d];
        for i in 0..d * d {
            defect[i] = finite(0.5 * (propagated[i] - storage[i]))?;
            scale[i] = finite(0.5 * (propagated[i].abs() + storage[i].abs()))?;
        }
        for matrix in self.component_loss.iter().chain(&self.connection_total).chain(&self.material_loss) {
            for i in 0..d * d {
                defect[i] = finite(defect[i] + matrix[i])?;
                scale[i] = finite(scale[i] + matrix[i].abs())?;
            }
        }
        // u^T (q(dt)-q(0)) is an EXACT held-force work quadratic. Symmetric
        // input/state cross terms count once in z^T work z.
        for i in 0..n {
            let qi = 2 * i;
            let ui = state_count + i;
            for j in 0..d {
                let initial = if qi == j { 1.0 } else { 0.0 };
                let work = finite(0.5 * (self.transition[qi * d + j] - initial) / self.weights[i])?;
                defect[ui * d + j] = finite(defect[ui * d + j] - work)?;
                defect[j * d + ui] = finite(defect[j * d + ui] - work)?;
                scale[ui * d + j] = finite(scale[ui * d + j] + work.abs())?;
                scale[j * d + ui] = finite(scale[j * d + ui] + work.abs())?;
            }
        }
        if matrix_norm(&defect, d)? > MATRIX_TOL * matrix_norm(&scale, d)?.max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared complete transition/loss matrix energy identity"));
        }
        Ok(())
    }
}

impl LedgerData {
    /// Each link's own loss Gramian (None without damping), computed by the first caller and shared after: one Van Loan
    /// integral per dissipative link, the cost the construction no longer pays. A deferred observation has no gate of its
    /// own and runs the admitted work to completion.
    fn gramians(&self, gate: Option<&CancelGate>) -> Result<&[Option<Vec<f64>>], ModalCouplingError> {
        if let Some(done) = self.gramians.get() {
            return Ok(done);
        }
        let idle = CancelGate::new_clock_free();
        let gate = gate.unwrap_or(&idle);
        let d = 3 * self.count + self.memory_count;
        let mut computed = Vec::with_capacity(self.dampings.len());
        for ((column, &damping), own) in self.columns.iter().zip(&self.dampings).zip(&self.damping_columns) {
            poll(Some(gate))?;
            let Some(integral) = self.integral.as_ref().filter(|_| damping != 0.0) else {
                computed.push(None);
                continue;
            };
            let resistance = link_resistance(own.as_deref().unwrap_or(column), damping, self.count, d)?;
            computed.push(Some(gramian(&integral.generator, &resistance, &integral.levels, integral.h, d,
                integral.budget, gate)?));
        }
        Ok(self.gramians.get_or_init(|| computed))
    }

    /// Work delivered to each port over one accepted step: the spring's stored-energy change, less the port's own viscous loss
    /// z^T W_k z. The per-port losses must add up to the step's total connection loss `connection_loss` (which comes from the
    /// one integral of the summed resistance, a different computation), so this is an independent check of both.
    fn work(
        &self, x: &[f64], u: &[f64], endpoint: &[f64], dt: f64, connection_loss: f64,
        columns: &[Vec<f64>], springs: &[(f64, f64)], gate: Option<&CancelGate>,
    ) -> Result<Vec<PortTransfer>, ModalCouplingError> {
        if springs.len() != self.dampings.len() || columns.len() != springs.len() || endpoint.len() != x.len() {
            return Err(invalid("prepared connection observation layout"));
        }
        self.evaluations.fetch_add(1, Ordering::Relaxed);
        let z = state_vector(self.count, self.memory_count, &self.weights, self.dt, x, u, dt)?;
        let gramians = self.gramians(gate)?;
        let mut transfers = Vec::with_capacity(springs.len());
        let mut total = 0.0;
        for ((column, &(stiffness, rest_extension)), gramian) in columns.iter().zip(springs).zip(gramians) {
            poll(gate)?;
            let coordinate = |state: &[f64]| -> Result<f64, ModalCouplingError> {
                column.iter().enumerate().try_fold(-rest_extension,
                    |s, (i, b)| finite(s + b * state[2 * i]))
            };
            let before = coordinate(x)?;
            let after = coordinate(endpoint)?;
            let loss = match gramian { Some(w) => quadratic(w, &z)?, None => 0.0 };
            total = finite(total + loss)?;
            let work = finite(-0.5 * stiffness * (after - before) * (after + before) - loss)?;
            transfers.push(PortTransfer::work(work, false));
        }
        if (total - connection_loss).abs() > MATRIX_TOL * total.abs().max(connection_loss.abs()).max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared per-port losses do not sum to connection loss"));
        }
        Ok(transfers)
    }
}

impl DeferredTransfers for PendingLedger {
    fn evaluate(&self) -> Result<Vec<PortTransfer>, fs_ports::CoupleError> {
        let data = &self.data;
        data.work(&self.x, &self.u, &self.endpoint, data.dt, self.connection_loss, &data.columns, &data.springs, None)
            .map_err(|error| fs_ports::CoupleError::ModuleObservationFailed {
                owner: "coupled-modal".into(),
                reason: error.to_string(),
            })
    }
}

/// z = [D q, v, memory, u] of the complete state and held input at the network's clock.
fn state_vector(count: usize, memory: usize, weights: &[f64], expected_dt: f64, x: &[f64], u: &[f64], dt: f64)
    -> Result<Vec<f64>, ModalCouplingError>
{
    let state_count = 2 * count + memory;
    if dt.to_bits() != expected_dt.to_bits() || x.len() != state_count || u.len() != count
        || x.iter().chain(u).any(|v| !v.is_finite()) {
        return Err(invalid("prepared exact clock or complete finite state/input changed"));
    }
    let mut z = vec![0.0; state_count + count];
    for i in 0..count {
        z[2 * i] = finite(weights[i] * x[2 * i])?;
        z[2 * i + 1] = x[2 * i + 1];
        z[state_count + i] = u[i];
    }
    z[2 * count..state_count].copy_from_slice(&x[2 * count..state_count]);
    Ok(z)
}

/// The damping column `g(K) b` of every link anchored to a constant loss factor, None for a Kelvin-Voigt link (see the
/// module notes). `stiffness` is the complete n x n stiffness with identity mass. Without an anchored link nothing is
/// solved. K is positive semidefinite: an eigenvalue below zero by more than the solver's roundoff refuses, and one
/// within it is a free coordinate, where `G(omega) / omega` has its finite limit.
fn constant_loss_columns(
    stiffness: &[f64], columns: &[Vec<f64>], links: &[ConnectionLaw], count: usize, gate: &CancelGate,
) -> Result<Vec<Option<Vec<f64>>>, ModalCouplingError> {
    if links.iter().all(|link| link.constant_loss_reference_omega_rad_s.is_none()) {
        return Ok(vec![None; links.len()]);
    }
    poll(Some(gate))?;
    let _admitted = fs_la::eigen::admit_jacobi_eigh(count).map_err(|_| invalid("constant-loss coupled eigensolve size"))?;
    let (values, vectors) = fs_la::eigen::householder_eigh(stiffness, count);
    poll(Some(gate))?;
    let largest = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let floor = 64.0 * count as f64 * f64::EPSILON * largest;
    if values.iter().chain(&vectors).any(|v| !v.is_finite()) || values.iter().any(|&v| v < -floor) {
        return Err(invalid("constant-loss coupled stiffness is not finite positive semidefinite"));
    }
    // G(omega) / omega = sum_j g_j / (rate_j (1 + (omega / rate_j)^2)): finite at a free coordinate.
    let spectrum: Vec<_> = ConstantLossSpectrum::fit().branches().collect();
    let per_omega: Vec<f64> = values.iter().map(|&value| {
        let omega = det::sqrt(value.max(0.0));
        spectrum.iter().map(|&(rate, g)| { let x = omega / rate; g / (rate * (1.0 + x * x)) }).sum::<f64>()
    }).collect();
    let mut result = Vec::with_capacity(links.len());
    for (column, link) in columns.iter().zip(links) {
        poll(Some(gate))?;
        let Some(reference) = link.constant_loss_reference_omega_rad_s else {
            result.push(None);
            continue;
        };
        // d = V diag(s) V^T b: project on the coupled modes, scale each by s_a, and return.
        let mut scaled = vec![0.0; count];
        for (a, (value, ratio)) in scaled.iter_mut().zip(&per_omega).enumerate() {
            let projection = (0..count).fold(0.0, |s, i| s + vectors[i * count + a] * column[i]);
            *value = finite(det::sqrt(reference * ratio) * projection)?;
        }
        let damping_column = (0..count)
            .map(|i| finite(vectors[i * count..(i + 1) * count].iter().zip(&scaled).fold(0.0, |s, (v, x)| s + v * x)))
            .collect::<Result<Vec<_>, _>>()?;
        result.push(Some(damping_column));
    }
    Ok(result)
}

/// One link's own resistance: its dashpot c b b^T in the velocity block.
fn link_resistance(column: &[f64], damping: f64, count: usize, d: usize) -> Result<Vec<f64>, ModalCouplingError> {
    let mut resistance = vec![0.0; d * d];
    for i in 0..count {
        for j in 0..count {
            resistance[(2 * i + 1) * d + 2 * j + 1] = finite(damping * column[i] * column[j])?;
        }
    }
    Ok(resistance)
}

fn matrix_norm(a: &[f64], n: usize) -> Result<f64, ModalCouplingError> {
    let mut norm = 0.0_f64;
    for row in a.chunks_exact(n) {
        norm = norm.max(row.iter().try_fold(0.0, |s, v| finite(s + v.abs()))?);
    }
    Ok(norm)
}

fn matrix_agrees(a: &[f64], b: &[f64], n: usize, what: &'static str) -> Result<(), ModalCouplingError> {
    let mut residual = 0.0_f64;
    for (left, right) in a.chunks_exact(n).zip(b.chunks_exact(n)) {
        residual = residual.max(left.iter().zip(right).try_fold(0.0,
            |s, (x, y)| finite(s + (x - y).abs()))?);
    }
    let scale = matrix_norm(a, n)?.max(matrix_norm(b, n)?).max(f64::MIN_POSITIVE);
    if residual > MATRIX_TOL * scale { return Err(invalid(what)); }
    Ok(())
}

fn product(a: &[f64], b: &[f64], n: usize, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError> {
    let mut result = vec![0.0; n * n];
    // Every entry is the k = 0, 1, ... sum of the textbook i-k-j loop, bit for bit. Four rows share each row of b while it
    // is in cache, and the slice update vectorizes over j.
    for (a_rows, out_rows) in a.chunks(4 * n).zip(result.chunks_mut(4 * n)) {
        poll(Some(gate))?;
        for (k, b_row) in b.chunks_exact(n).enumerate() {
            for (a_row, out_row) in a_rows.chunks_exact(n).zip(out_rows.chunks_exact_mut(n)) {
                let scale = a_row[k];
                for (value, &right) in out_row.iter_mut().zip(b_row) {
                    *value += scale * right;
                }
            }
        }
    }
    if result.iter().any(|v| !v.is_finite()) { return Err(invalid("prepared dense product overflow")); }
    Ok(result)
}

fn sandwich(f: &[f64], w: &[f64], n: usize, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError> {
    let wf = product(w, f, n, gate)?;
    let transpose: Vec<_> = (0..n * n).map(|i| f[(i % n) * n + i / n]).collect();
    product(&transpose, &wf, n, gate)
}

/// The transition of `generator` over `dt` by scaling and squaring: every scale level (the last is the transition) and the
/// short interval h.
fn scaled_levels(generator: &[f64], d: usize, dt: f64, budget: DenseExponentialBudget, gate: &CancelGate)
    -> Result<(Vec<Vec<f64>>, f64), ModalCouplingError>
{
    let norm = matrix_norm(generator, d)?;
    let mut h = dt;
    let mut scales = 0;
    while norm * h > 0.5 {
        if scales == MAX_SCALES {
            return Err(invalid("prepared exact scaling exceeds 48 doublings"));
        }
        h *= 0.5;
        scales += 1;
    }
    let scaled: Vec<_> = generator.iter().map(|v| v * h).collect();
    let mut levels = Vec::with_capacity(scales + 1);
    levels.push(exponential_with_budget(&scaled, d, budget).map_err(ModalCouplingError::Propagation)?);
    for i in 0..scales {
        levels.push(product(&levels[i], &levels[i], d, gate)?);
    }
    // Compare the complete transition against one finer independent series
    // evaluation. Extra squaring roundoff is admitted, not silently trusted.
    let finer: Vec<_> = scaled.iter().map(|v| v * 0.5).collect();
    let mut check = exponential_with_budget(&finer, d, budget).map_err(ModalCouplingError::Propagation)?;
    for _ in 0..=scales { check = product(&check, &check, d, gate)?; }
    matrix_agrees(levels.last().unwrap(), &check, d, "prepared transition refinement residual")?;
    Ok((levels, h))
}

fn gramian(
    generator: &[f64], resistance: &[f64], levels: &[Vec<f64>], h: f64,
    d: usize, budget: DenseExponentialBudget, gate: &CancelGate,
) -> Result<Vec<f64>, ModalCouplingError> {
    poll(Some(gate))?;
    // Van Loan's [[-A^T h, R h], [0, A h]] is block upper-triangular: its exponential is formed by blocks, the same bits as
    // the assembled 2d x 2d exponential at half its products.
    let mut van_loan = UpperBlocks { first: vec![0.0; d * d], cross: vec![0.0; d * d], last: vec![0.0; d * d] };
    for i in 0..d {
        for j in 0..d {
            van_loan.first[i * d + j] = finite(-generator[j * d + i] * h)?;
            van_loan.last[i * d + j] = finite(generator[i * d + j] * h)?;
            van_loan.cross[i * d + j] = finite(resistance[i * d + j] * h)?;
        }
    }
    let UpperBlocks { cross, last: e, .. } = upper_block_exponential_with_budget(&van_loan, d, budget)
        .map_err(ModalCouplingError::Propagation)?;
    matrix_agrees(&e, &levels[0], d, "prepared loss/transition short-interval disagreement")?;
    let transpose: Vec<_> = (0..d * d).map(|i| e[(i % d) * d + i / d]).collect();
    let mut w = product(&transpose, &cross, d, gate)?;
    for level in levels.iter().take(levels.len() - 1) {
        let extra = sandwich(level, &w, d, gate)?;
        for (value, increment) in w.iter_mut().zip(extra) { *value = finite(*value + increment)?; }
    }
    let transpose: Vec<_> = (0..d * d).map(|i| w[(i % d) * d + i / d]).collect();
    matrix_agrees(&w, &transpose, d, "prepared loss Gramian symmetry residual")?;
    // Only antisymmetric floating roundoff is removed. No spectral projection,
    // eigenvalue clamp or omitted channel is used to make a loss look positive.
    for i in 0..d {
        for j in i + 1..d {
            let value = finite(0.5 * w[i * d + j] + 0.5 * w[j * d + i])?;
            w[i * d + j] = value;
            w[j * d + i] = value;
        }
    }
    Ok(w)
}

fn apply(matrix: &[f64], z: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
    // Each entry is its row's left-to-right sum, bit for bit; eight rows advance together so their independent sums overlap.
    // A partial sum that leaves the finite range never returns to it (inf + finite = inf, inf - inf = NaN, NaN stays), so
    // checking each finished row refuses exactly what a check of every partial sum refused.
    let n = z.len();
    let mut result = Vec::with_capacity(matrix.len() / n);
    let mut blocks = matrix.chunks_exact(8 * n);
    for rows in blocks.by_ref() {
        let rows: [&[f64]; 8] = std::array::from_fn(|k| &rows[k * n..(k + 1) * n]);
        let mut sums = [0.0_f64; 8];
        for (j, &value) in z.iter().enumerate() {
            for (sum, row) in sums.iter_mut().zip(&rows) {
                *sum += row[j] * value;
            }
        }
        for sum in sums {
            result.push(finite(sum)?);
        }
    }
    for row in blocks.remainder().chunks_exact(n) {
        result.push(finite(row.iter().zip(z).fold(0.0, |s, (a, b)| s + a * b))?);
    }
    Ok(result)
}

fn quadratic(matrix: &[f64], z: &[f64]) -> Result<f64, ModalCouplingError> {
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for (row, &left) in matrix.chunks_exact(z.len()).zip(z) {
        for (&value, &right) in row.iter().zip(z) {
            let term = left * value * right;
            sum += term;
            magnitude += term.abs();
        }
    }
    checked_loss(sum, magnitude)
}

/// Both loss quadratics of one vector in one pass over the two Gramians. Each keeps its own left-to-right sum and magnitude,
/// so each equals `quadratic` of its matrix bit for bit, refusals included. The four independent sums overlap in the
/// pipeline instead of one quadratic running after the other.
fn quadratic_pair(first: &[f64], second: &[f64], z: &[f64]) -> Result<(f64, f64), ModalCouplingError> {
    let n = z.len();
    let (mut sum_first, mut magnitude_first, mut sum_second, mut magnitude_second) = (0.0, 0.0, 0.0, 0.0);
    for ((row_first, row_second), &left) in first.chunks_exact(n).zip(second.chunks_exact(n)).zip(z) {
        for ((&a, &b), &right) in row_first.iter().zip(row_second).zip(z) {
            let term_first = left * a * right;
            let term_second = left * b * right;
            sum_first += term_first;
            magnitude_first += term_first.abs();
            sum_second += term_second;
            magnitude_second += term_second.abs();
        }
    }
    // The order of the refusals is the order of the two separate evaluations: the first matrix's, then the second's.
    let first = checked_loss(sum_first, magnitude_first)?;
    Ok((first, checked_loss(sum_second, magnitude_second)?))
}

/// The refusals of a finished quadratic. A term or partial sum that leaves the finite range never returns to it
/// (inf + finite = inf, inf - inf = NaN, NaN stays), so a non-finite sum or magnitude refuses exactly what a check of every
/// term and partial sum refused.
fn checked_loss(sum: f64, magnitude: f64) -> Result<f64, ModalCouplingError> {
    let (sum, magnitude) = (finite(sum)?, finite(magnitude)?);
    // A Gramian is PSD in exact arithmetic. Only a negative value below its
    // explicit summation-roundoff allowance is reported as zero physical loss.
    // A larger negative quadratic refuses the trial and preserves accepted state.
    let roundoff = 4096.0 * f64::EPSILON * magnitude.max(f64::MIN_POSITIVE);
    if sum < -roundoff { return Err(invalid("prepared resistance integral is negative beyond roundoff")); }
    Ok(sum.max(0.0))
}

pub(super) mod air;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod constant_loss_tests;
#[cfg(test)]
mod air_tests;
#[cfg(test)]
mod memory_tests;
