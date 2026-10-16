//! Matrix-free version of the existing centroid Burton-Miller discretization.
//!
//! Adjacent-leaf pairs use fs-fmm's exact point kernels. Self integrals and
//! jumps reuse the dense owner's equivalent-disc rule. Static row closures use
//! the SAME FMM discretization as the dynamic operator, preserving N0[1]=0.
//! FMM approximation error is distinct from the recomputed GMRES residual;
//! callers must qualify the explicit order against dense/refinement oracles.

use super::{
    Formulation, HelmholtzError, MAX_DENSE_PANELS, MIN_PANELS_PER_WAVELENGTH, Medium,
    RadiationSolution, alpha_for, characteristic_panel_size, evaluate_radiated_power, kernels,
    self_terms, surface_fingerprint,
};
use crate::{BemError, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_fmm::{
    FmmC64, FmmError, HelmholtzAdjointDouble, HelmholtzDouble, HelmholtzHypersingular,
    HelmholtzSingle,
};
use fs_la::eigen_complex::{LuComplex, lu_complex};
use fs_math::c64::C64;
use fs_solver::{
    krylov::{GmresState, SolveReport},
    op::LinearOp,
};
use std::{
    cell::{Cell, RefCell},
    time::Instant,
};

/// Right preconditioning leaves GMRES's residual in the original equation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preconditioner {
    /// Original equation, without inverse application.
    None,
    /// Exact LU blocks drawn from individual direct-interaction FMM leaves.
    /// Large leaves are split in their stable original point order.
    NearFieldBlockJacobi { max_block_panels: usize },
}

/// Explicit accuracy, restart, FMM and near-field storage budgets.
#[derive(Clone, Debug)]
pub struct IterativeOptions {
    pub tol: f64,
    pub restart: usize,
    pub max_iters: usize,
    /// Explicit interpolation order; never silently increased or reduced.
    pub fmm_order: usize,
    pub leaf_capacity: usize,
    pub precond: Preconditioner,
    /// Maximum stored complex entries across all inverse-preconditioner blocks.
    pub max_near_entries: usize,
}

impl Default for IterativeOptions {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            restart: 60,
            max_iters: 600,
            fmm_order: 4,
            leaf_capacity: 128,
            precond: Preconditioner::NearFieldBlockJacobi {
                max_block_panels: 32,
            },
            max_near_entries: 4_194_304,
        }
    }
}

/// Measured process evidence, separate from numerical convergence authority.
#[derive(Clone, Debug)]
pub struct ResourceReceipt {
    pub host: String,
    /// Preparation, RHS and committed/aborted solve time for this state.
    pub wall_s: f64,
    /// Linux process VmHWM; None where that measurement is unavailable.
    /// This is process-wide high-water RSS, not a per-solve allocation claim.
    pub peak_rss_bytes: Option<u64>,
    pub matvecs: usize,
    pub iterations: usize,
    pub fmm_order: usize,
    pub near_field_bytes: usize,
    pub fmm_stats: [String; 2],
}

/// Refusals retain diagnostics; an unconverged pressure is never a solution.
#[derive(Clone, Debug)]
pub enum IterativeError {
    Geometry(BemError),
    Helmholtz(HelmholtzError),
    Fmm(FmmError),
    Cancelled {
        matvecs: usize,
    },
    Invalid {
        what: &'static str,
    },
    NearFieldBudget {
        required_entries: usize,
        maximum_entries: usize,
    },
    SingularNearBlock {
        block: usize,
    },
    NotConverged {
        pressure: Vec<C64>,
        report: SolveReport,
        receipt: ResourceReceipt,
    },
}
impl core::fmt::Display for IterativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry(e) => e.fmt(f),
            Self::Helmholtz(e) => e.fmt(f),
            Self::Fmm(e) => e.fmt(f),
            Self::Cancelled { matvecs } => {
                write!(f, "FS-BEM-ITERATIVE-CANCELLED: after {matvecs} matvecs")
            }
            Self::Invalid { what } => write!(f, "FS-BEM-ITERATIVE-INVALID: {what}"),
            Self::NearFieldBudget {
                required_entries,
                maximum_entries,
            } => write!(
                f,
                "FS-BEM-ITERATIVE-NEAR-BUDGET: {required_entries} entries exceed {maximum_entries}"
            ),
            Self::SingularNearBlock { block } => {
                write!(f, "FS-BEM-ITERATIVE-SINGULAR-BLOCK: {block}")
            }
            Self::NotConverged { report, .. } => write!(
                f,
                "FS-BEM-ITERATIVE-NOT-CONVERGED: {} iterations, relative residual {}",
                report.iters, report.rel_residual
            ),
        }
    }
}
impl std::error::Error for IterativeError {}
impl From<HelmholtzError> for IterativeError {
    fn from(e: HelmholtzError) -> Self {
        Self::Helmholtz(e)
    }
}
impl From<FmmError> for IterativeError {
    fn from(e: FmmError) -> Self {
        Self::Fmm(e)
    }
}

struct NearBlock {
    indices: Vec<usize>,
    inverse: LuComplex,
}

/// Prepared immutable geometry/operator, reusable across boundary velocities.
/// Only Burton-Miller and its explicit sign-mutation arm are supported. The
/// triangle PlainCbie quadrature is NOT replaced with centroid FMM kernels.
pub struct IterativeRadiation<'a> {
    surface: &'a SpherePanels,
    k: f64,
    medium: Medium,
    alpha: C64,
    options: IterativeOptions,
    ppw: f64,
    double: FmmC64<'static>,
    mixed: FmmC64<'static>,
    diagonal_a: Vec<C64>,
    diagonal_b: Vec<C64>,
    blocks: Vec<NearBlock>,
    near_entries: usize,
    setup_s: f64,
    stamp: u64,
}

/// Clone is a numerical checkpoint at a completed GMRES restart boundary.
/// Wall-time counters are diagnostic and do not participate in replay equality.
#[derive(Clone, Debug)]
pub struct IterativeState {
    gmres: GmresState,
    rhs: Vec<f64>,
    velocity: Vec<C64>,
    stamp: u64,
    matvecs: usize,
    wall_s: f64,
}
impl IterativeState {
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.gmres.iters
    }
    #[must_use]
    pub fn report(&self, tolerance: f64) -> SolveReport {
        SolveReport::from_claim(
            self.gmres.iters,
            self.gmres.residual_claim(),
            tolerance,
            self.gmres.history.clone(),
        )
    }
}

fn poll(gate: &CancelGate) -> Result<(), IterativeError> {
    if gate.is_requested() {
        Err(IterativeError::Cancelled { matvecs: 0 })
    } else {
        Ok(())
    }
}

impl<'a> IterativeRadiation<'a> {
    /// Prepare without allocating any N-by-N matrix. Near blocks are explicitly
    /// bounded, and FMM keeps its original coefficient/translation/P2P budgets.
    pub fn prepare(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        formulation: Formulation,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        let started = Instant::now();
        poll(gate)?;
        surface.validate().map_err(IterativeError::Geometry)?;
        if !k.is_finite()
            || k <= 0.0
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(IterativeError::Invalid {
                what: "positive finite k, density and sound speed required",
            });
        }
        if formulation == Formulation::PlainCbie {
            return Err(IterativeError::Invalid {
                what: "iterative owner requires Burton-Miller; triangle CBIE keeps its dense quadrature",
            });
        }
        if !options.tol.is_finite()
            || options.tol <= 0.0
            || options.tol >= 1.0
            || options.restart == 0
            || options.restart > 256
        {
            return Err(IterativeError::Invalid {
                what: "0 < tol < 1 and 1 <= restart <= 256 required",
            });
        }
        if let Preconditioner::NearFieldBlockJacobi { max_block_panels } = options.precond {
            if !(1..=128).contains(&max_block_panels) {
                return Err(IterativeError::Invalid {
                    what: "near blocks require 1..=128 panels",
                });
            }
        }
        let ppw = 2.0 * core::f64::consts::PI / (k * characteristic_panel_size(surface));
        if ppw < MIN_PANELS_PER_WAVELENGTH {
            return Err(HelmholtzError::TooCoarse {
                panels_per_wavelength: ppw,
            }
            .into());
        }
        let points = || surface.centroids().to_vec();
        let normals = || surface.normals().to_vec();
        let double = FmmC64::new_double(
            HelmholtzDouble { k },
            points(),
            normals(),
            options.fmm_order,
            options.leaf_capacity,
        )?;
        poll(gate)?;
        let mixed = FmmC64::new_hypersingular(
            HelmholtzHypersingular { k },
            points(),
            normals(),
            options.fmm_order,
            options.leaf_capacity,
        )?;
        let weights: Vec<_> = surface.areas().iter().map(|&a| C64::from_re(a)).collect();
        // Both static closures use the same point distribution, traversal and
        // explicit order as the actual operator. No dense row scans occur.
        let d0 = {
            poll(gate)?;
            let op = FmmC64::new_double(
                HelmholtzDouble { k: 0.0 },
                points(),
                normals(),
                options.fmm_order,
                options.leaf_capacity,
            )?;
            op.potentials(&weights)?
        };
        let n0 = {
            poll(gate)?;
            let op = FmmC64::new_hypersingular(
                HelmholtzHypersingular { k: 0.0 },
                points(),
                normals(),
                options.fmm_order,
                options.leaf_capacity,
            )?;
            op.potentials(&weights)?
        };
        poll(gate)?;
        let alpha = alpha_for(formulation, k);
        let mut diagonal_a = Vec::with_capacity(weights.len());
        let mut diagonal_b = Vec::with_capacity(weights.len());
        for i in 0..weights.len() {
            let (s, n_reg) = self_terms(k, surface.areas()[i]);
            let d_self = C64::from_re(-0.5) - d0[i];
            diagonal_a.push(d_self - C64::from_re(0.5) + alpha * (n_reg - n0[i]));
            diagonal_b.push(s + alpha * (d_self + C64::from_re(0.5)));
        }
        let groups: Vec<Vec<usize>> = match options.precond {
            Preconditioner::None => Vec::new(),
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => double
                .near_field_blocks()
                .flat_map(|leaf| leaf.chunks(max_block_panels).map(<[usize]>::to_vec))
                .collect(),
        };
        let near_entries = groups
            .iter()
            .try_fold(0usize, |count, block| {
                block
                    .len()
                    .checked_mul(block.len())
                    .and_then(|entries| count.checked_add(entries))
            })
            .ok_or(IterativeError::Invalid {
                what: "near-field extent overflow",
            })?;
        if near_entries > options.max_near_entries {
            return Err(IterativeError::NearFieldBudget {
                required_entries: near_entries,
                maximum_entries: options.max_near_entries,
            });
        }
        let mut blocks = Vec::with_capacity(groups.len());
        for indices in groups {
            poll(gate)?;
            let m = indices.len();
            let mut matrix = vec![C64::ZERO; m * m];
            for (a, &i) in indices.iter().enumerate() {
                for (b, &j) in indices.iter().enumerate() {
                    matrix[a * m + b] = if i == j {
                        diagonal_a[i]
                    } else {
                        // This is the SAME exact P2P quadrature as dense BM and
                        // FMM's near field, not a second panel integration rule.
                        let (_, d, _, h) = kernels(
                            k,
                            surface.centroids()[i],
                            surface.normals()[i],
                            surface.centroids()[j],
                            surface.normals()[j],
                        );
                        (d + alpha * h).scale(surface.areas()[j])
                    };
                }
            }
            let inverse =
                lu_complex(&matrix, m).map_err(|_| IterativeError::SingularNearBlock {
                    block: blocks.len(),
                })?;
            blocks.push(NearBlock { indices, inverse });
        }
        // Cross-wiring guard for numerical checkpoints, not cryptographic proof.
        let mut stamp = surface_fingerprint(surface);
        let preconditioner_key = match options.precond {
            Preconditioner::None => 0,
            Preconditioner::NearFieldBlockJacobi { max_block_panels } => max_block_panels as u64,
        };
        for word in [
            k.to_bits(),
            alpha.im.to_bits(),
            medium.density.to_bits(),
            medium.sound_speed.to_bits(),
            options.tol.to_bits(),
            options.restart as u64,
            options.max_iters as u64,
            options.fmm_order as u64,
            options.leaf_capacity as u64,
            near_entries as u64,
            preconditioner_key,
        ] {
            stamp = (stamp ^ word).wrapping_mul(0x0000_0100_0000_01b3);
        }
        Ok(Self {
            surface,
            k,
            medium,
            alpha,
            options: options.clone(),
            ppw,
            double,
            mixed,
            diagonal_a,
            diagonal_b,
            blocks,
            near_entries,
            setup_s: started.elapsed().as_secs_f64(),
            stamp,
        })
    }

    /// Compute Bq, with q=i*k*c*rho*v, using the original dense B signs.
    pub fn start(
        &self,
        velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<IterativeState, IterativeError> {
        let started = Instant::now();
        poll(gate)?;
        if velocity.len() != self.surface.areas().len() || velocity.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "one finite complex normal velocity per panel required",
            });
        }
        let q: Vec<_> = velocity
            .iter()
            .map(|&v| v * C64::new(0.0, self.k * self.medium.sound_speed * self.medium.density))
            .collect();
        let weighted: Vec<_> = q
            .iter()
            .zip(self.surface.areas())
            .map(|(&q, &a)| q.scale(a))
            .collect();
        let kernel = HelmholtzSingle { k: self.k };
        let single = FmmC64::new(
            &kernel,
            self.surface.centroids().to_vec(),
            self.options.fmm_order,
            self.options.leaf_capacity,
        )?
        .potentials(&weighted)?;
        poll(gate)?;
        let adjoint = FmmC64::new_adjoint_double(
            HelmholtzAdjointDouble { k: self.k },
            self.surface.centroids().to_vec(),
            self.surface.normals().to_vec(),
            self.options.fmm_order,
            self.options.leaf_capacity,
        )?
        .potentials(&weighted)?;
        poll(gate)?;
        let rhs_complex: Vec<_> = (0..q.len())
            .map(|i| single[i] + self.alpha * adjoint[i] + self.diagonal_b[i] * q[i])
            .collect();
        if rhs_complex.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "nonfinite radiation right-hand side",
            });
        }
        let rhs = stack(&rhs_complex);
        let mut gmres = GmresState::new(&rhs, self.options.restart.min(rhs.len()));
        if rhs.iter().all(|&v| v == 0.0) {
            gmres.run(&ZeroOp(rhs.len()), &rhs, self.options.tol, 1, false);
        }
        Ok(IterativeState {
            gmres,
            rhs,
            velocity: velocity.to_vec(),
            stamp: self.stamp,
            matvecs: 0,
            wall_s: started.elapsed().as_secs_f64(),
        })
    }

    fn inverse(&self, input: &[C64]) -> Vec<C64> {
        if self.blocks.is_empty() {
            return input.to_vec();
        }
        let mut out = vec![C64::ZERO; input.len()];
        for block in &self.blocks {
            let mut values: Vec<_> = block.indices.iter().map(|&i| input[i]).collect();
            block.inverse.solve(&mut values);
            for (&i, value) in block.indices.iter().zip(values) {
                out[i] = value;
            }
        }
        out
    }

    /// Original unpreconditioned complex A application, including self closure.
    pub fn apply(&self, pressure: &[C64], gate: &CancelGate) -> Result<Vec<C64>, IterativeError> {
        poll(gate)?;
        if pressure.len() != self.surface.areas().len() {
            return Err(IterativeError::Invalid {
                what: "pressure/operator dimension mismatch",
            });
        }
        let weighted: Vec<_> = pressure
            .iter()
            .zip(self.surface.areas())
            .map(|(&p, &a)| p.scale(a))
            .collect();
        let d = self.double.potentials(&weighted)?;
        poll(gate)?;
        let h = self.mixed.potentials(&weighted)?;
        poll(gate)?;
        let out: Vec<_> = (0..pressure.len())
            .map(|i| d[i] + self.alpha * h[i] + self.diagonal_a[i] * pressure[i])
            .collect();
        if out.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::Invalid {
                what: "nonfinite Burton-Miller application",
            });
        }
        Ok(out)
    }

    /// Run additional restart cycles; cancellation/refusal discards the current
    /// partial cycle and retains the preceding checkpoint. max_iters counts
    /// actual Arnoldi updates, including a final partial restart cycle.
    pub fn advance(
        &self,
        state: &mut IterativeState,
        cycles: usize,
        gate: &CancelGate,
    ) -> Result<SolveReport, IterativeError> {
        self.advance_controlled(state, cycles, gate, &mut |_| false)
    }

    /// The predicate is polled immediately before each actual matvec, with the
    /// cumulative matvec index. A true result cancels and is latched for this
    /// call. This supports deterministic cancellation without a timing race.
    pub fn advance_controlled(
        &self,
        state: &mut IterativeState,
        cycles: usize,
        gate: &CancelGate,
        cancel: &mut impl FnMut(usize) -> bool,
    ) -> Result<SolveReport, IterativeError> {
        if state.stamp != self.stamp {
            return Err(IterativeError::Invalid {
                what: "checkpoint belongs to a different prepared operator",
            });
        }
        let started = Instant::now();
        let op = RealOp {
            owner: self,
            gate,
            first_error: RefCell::new(None),
            matvecs: Cell::new(0),
            offset: state.matvecs,
            cancel: RefCell::new(cancel),
        };
        for _ in 0..cycles {
            if gate.is_requested() {
                state.matvecs += op.matvecs.get();
                state.wall_s += started.elapsed().as_secs_f64();
                return Err(IterativeError::Cancelled {
                    matvecs: state.matvecs,
                });
            }
            if state.gmres.rel_residual() < self.options.tol
                || state.gmres.iters >= self.options.max_iters
            {
                break;
            }
            let mut trial = state.gmres.clone();
            trial.restart = self
                .options
                .restart
                .min(self.options.max_iters - trial.iters)
                .min(state.rhs.len());
            trial.run(&op, &state.rhs, self.options.tol, 1, false);
            if gate.is_requested() && op.first_error.borrow().is_none() {
                *op.first_error.borrow_mut() = Some(IterativeError::Cancelled { matvecs: 0 });
            }
            if let Some(e) = op.first_error.borrow_mut().take() {
                state.matvecs += op.matvecs.get();
                state.wall_s += started.elapsed().as_secs_f64();
                return Err(match e {
                    IterativeError::Cancelled { .. } => IterativeError::Cancelled {
                        matvecs: state.matvecs,
                    },
                    other => other,
                });
            }
            trial.restart = self.options.restart.min(state.rhs.len());
            state.gmres = trial;
        }
        state.matvecs += op.matvecs.get();
        state.wall_s += started.elapsed().as_secs_f64();
        Ok(state.report(self.options.tol))
    }

    /// Publish only a converged original-equation residual and finite pressure.
    pub fn finish(
        &self,
        state: &IterativeState,
    ) -> Result<(RadiationSolution, SolveReport, ResourceReceipt), IterativeError> {
        if state.stamp != self.stamp {
            return Err(IterativeError::Invalid {
                what: "checkpoint belongs to a different prepared operator",
            });
        }
        let pressure = self.inverse(&unstack(&state.gmres.x));
        let report = state.report(self.options.tol);
        let receipt = self.receipt(state);
        if !report.converged_euclidean() || pressure.iter().any(|v| !finite(*v)) {
            return Err(IterativeError::NotConverged {
                pressure,
                report,
                receipt,
            });
        }
        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(
                &pressure,
                &state.velocity,
                self.surface.areas(),
                self.medium,
            );
        Ok((
            RadiationSolution {
                pressure,
                velocity: state.velocity.clone(),
                k: self.k,
                medium: self.medium,
                surface_fingerprint: surface_fingerprint(self.surface),
                panels_per_wavelength: self.ppw,
                radiated_power,
                radiated_power_roundoff_interval,
                power_diagnostics,
                // Universal condition-number lower bound; no dense probe performed.
                condition_lower_bound: 1.0,
                dense_cap_utilization: self.surface.areas().len() as f64 / MAX_DENSE_PANELS as f64,
            },
            report,
            receipt,
        ))
    }

    #[must_use]
    pub fn receipt(&self, state: &IterativeState) -> ResourceReceipt {
        ResourceReceipt {
            host: host(),
            wall_s: self.setup_s + state.wall_s,
            peak_rss_bytes: peak_rss(),
            matvecs: state.matvecs,
            iterations: state.gmres.iters,
            fmm_order: self.options.fmm_order,
            near_field_bytes: self.near_entries * core::mem::size_of::<C64>()
                + self
                    .blocks
                    .iter()
                    .map(|b| b.indices.len() * core::mem::size_of::<usize>() * 2)
                    .sum::<usize>(),
            fmm_stats: [self.double.stats(), self.mixed.stats()],
        }
    }
}

struct RealOp<'a, 's> {
    owner: &'a IterativeRadiation<'s>,
    gate: &'a CancelGate,
    first_error: RefCell<Option<IterativeError>>,
    matvecs: Cell<usize>,
    offset: usize,
    cancel: RefCell<&'a mut dyn FnMut(usize) -> bool>,
}
impl LinearOp for RealOp<'_, '_> {
    fn n(&self) -> usize {
        2 * self.owner.surface.areas().len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.first_error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        let result = if (self.cancel.borrow_mut())(self.offset + self.matvecs.get()) {
            Err(IterativeError::Cancelled { matvecs: 0 })
        } else {
            poll(self.gate)
        }
        .and_then(|()| {
            self.matvecs.set(self.matvecs.get() + 1);
            self.owner
                .apply(&self.owner.inverse(&unstack(x)), self.gate)
        });
        match result {
            Ok(v) => y.copy_from_slice(&stack(&v)),
            Err(e) => {
                *self.first_error.borrow_mut() = Some(e);
                y.fill(f64::NAN);
            }
        }
    }
    fn apply_transpose(&self, _x: &[f64], y: &mut [f64]) {
        *self.first_error.borrow_mut() = Some(IterativeError::Invalid {
            what: "transpose iterative radiation is not implemented",
        });
        y.fill(f64::NAN);
    }
}
struct ZeroOp(usize);
impl LinearOp for ZeroOp {
    fn n(&self) -> usize {
        self.0
    }
    fn apply(&self, _: &[f64], y: &mut [f64]) {
        y.fill(0.0);
    }
}
fn finite(v: C64) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
fn stack(v: &[C64]) -> Vec<f64> {
    v.iter()
        .map(|z| z.re)
        .chain(v.iter().map(|z| z.im))
        .collect()
}
fn unstack(v: &[f64]) -> Vec<C64> {
    let n = v.len() / 2;
    (0..n).map(|i| C64::new(v[i], v[n + i])).collect()
}
fn host() -> String {
    let name = std::fs::read_to_string("/proc/sys/kernel/hostname").ok();
    let cpu = std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|text| {
            text.lines()
                .find_map(|line| line.strip_prefix("model name").map(str::to_string))
        });
    format!(
        "{} {} {} {}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        name.as_deref().unwrap_or("hostname-unavailable").trim(),
        cpu.as_deref().unwrap_or("cpu-unavailable").trim()
    )
}
fn peak_rss() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| {
            let value = line
                .strip_prefix("VmHWM:")?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()?;
            value.checked_mul(1024)
        })
}

/// One-shot convenience over the same checkpointable operator/state API.
pub fn solve_radiation_iterative(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity: &[C64],
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<(RadiationSolution, SolveReport, ResourceReceipt), IterativeError> {
    let operator = IterativeRadiation::prepare(surface, k, medium, formulation, options, gate)?;
    let mut state = operator.start(velocity, gate)?;
    operator.advance(&mut state, options.max_iters.saturating_add(1), gate)?;
    operator.finish(&state)
}
