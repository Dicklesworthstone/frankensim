//! Additional coupled-block solve modes for the composition kernel:
//!
//! - Derivative-free fixed-point iteration with bounded Aitken dynamic
//!   relaxation or bounded vector IQN-ILS (`Component` mandatory tier only).
//! - Monolithic Newton--Krylov via `fs-solver`'s fallible adapter
//!   ([`fs_solver::solve_fallible_newton_krylov`]) when every component in the
//!   looped block offers [`fs_ports::component::StageResidual`], carrying
//!   `Result` end-to-end without ever encoding refusal as `NaN`.

use core::cell::RefCell;

use fs_exec::Cx;
use fs_solver::{
    FallibleNewtonError, FallibleNonlinearProblem, NewtonKrylovConfig, solve_fallible_newton_krylov,
};

use super::{BlockWork, ROUNDING_UNITS, System, SystemError};
use crate::accelerator::{AitkenConfig, IqnIlsConfig};
use crate::diagnostics::{
    BlockSolveMode, BlockSolverPolicy, BlockStepDiagnostic, FixedPointMethod, MonolithicConfig,
};

impl System {
    /// Choose the solve mode for `block` from the component tiers present and
    /// the configured [`BlockSolverPolicy`].
    #[must_use]
    pub(super) fn select_block_mode(&self, block: usize) -> BlockSolveMode {
        if !self.blocks[block].looped() {
            return BlockSolveMode::OrderedSweep;
        }
        if !self.lines[block].is_empty() {
            return BlockSolveMode::LoadLineNewton;
        }
        let nodes = &self.blocks[block].nodes;
        let all_stage = !nodes.is_empty()
            && nodes
                .iter()
                .all(|&n| self.nodes[n].has_stage_residual() && self.nodes[n].stage_unknowns() > 0);
        let any_tangent = nodes.iter().any(|&n| self.nodes[n].has_tangent());
        match self.solver_policy {
            BlockSolverPolicy::Auto => {
                if all_stage {
                    BlockSolveMode::MonolithicNewtonKrylov
                } else {
                    BlockSolveMode::InterfaceNewton
                }
            }
            BlockSolverPolicy::ByTier { fixed_point, .. } => {
                if all_stage {
                    BlockSolveMode::MonolithicNewtonKrylov
                } else if any_tangent {
                    BlockSolveMode::InterfaceNewton
                } else {
                    match fixed_point {
                        FixedPointMethod::Aitken(_) => BlockSolveMode::FixedPointAitken,
                        FixedPointMethod::IqnIls(_) => BlockSolveMode::FixedPointIqnIls,
                    }
                }
            }
            BlockSolverPolicy::FixedPointAitken(_) => BlockSolveMode::FixedPointAitken,
            BlockSolverPolicy::FixedPointIqnIls(_) => BlockSolveMode::FixedPointIqnIls,
            BlockSolverPolicy::InterfaceNewton => BlockSolveMode::InterfaceNewton,
            BlockSolverPolicy::MonolithicNewtonKrylov(_) => BlockSolveMode::MonolithicNewtonKrylov,
        }
    }

    fn aitken_config(&self) -> AitkenConfig {
        match self.solver_policy {
            BlockSolverPolicy::FixedPointAitken(cfg)
            | BlockSolverPolicy::ByTier {
                fixed_point: FixedPointMethod::Aitken(cfg),
                ..
            } => cfg,
            BlockSolverPolicy::FixedPointIqnIls(cfg)
            | BlockSolverPolicy::ByTier {
                fixed_point: FixedPointMethod::IqnIls(cfg),
                ..
            } => cfg.fallback,
            _ => AitkenConfig::default(),
        }
    }

    fn iqn_config(&self) -> IqnIlsConfig {
        match self.solver_policy {
            BlockSolverPolicy::FixedPointIqnIls(cfg)
            | BlockSolverPolicy::ByTier {
                fixed_point: FixedPointMethod::IqnIls(cfg),
                ..
            } => cfg,
            _ => IqnIlsConfig::default(),
        }
    }

    fn monolithic_config(&self) -> MonolithicConfig {
        match self.solver_policy {
            BlockSolverPolicy::MonolithicNewtonKrylov(cfg)
            | BlockSolverPolicy::ByTier {
                monolithic: cfg, ..
            } => cfg,
            _ => MonolithicConfig::default(),
        }
    }

    /// Solve one looped block using derivative-free fixed-point iteration with
    /// bounded Aitken or vector IQN-ILS acceleration on the torn port variables.
    ///
    /// Mutates only `self.candidate_controller.blocks[block]` and pre-allocated
    /// scratch in `work`; accepted state and accelerator history remain bitwise
    /// untouched until step publication.
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub(super) fn solve_fixed_point_block(
        &mut self,
        block: usize,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        mode: BlockSolveMode,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        let m = unknowns.len();
        for &index in unknowns {
            let start = self.accepted_inputs[index];
            if self.input_is_temperature[index] && !(start.is_finite() && start > 0.0) {
                return Err(SystemError::GuessRequired {
                    port: self.input_label[index].clone(),
                });
            }
            self.inputs[index] = start;
        }
        work.group_scale.resize(self.group_count, 0.0);
        work.allowance.fill(0.0);

        let aitken_cfg = self.aitken_config();
        let iqn_cfg = self.iqn_config();
        let keep_history = match mode {
            BlockSolveMode::FixedPointIqnIls => iqn_cfg.reuse_history_across_steps,
            _ => true,
        };
        let init_omega = match mode {
            BlockSolveMode::FixedPointAitken => aitken_cfg.initial_relaxation,
            _ => iqn_cfg.fallback.initial_relaxation,
        };
        self.candidate_controller.blocks[block].begin_step(keep_history, init_omega);

        let eval_start = self.step_tracker.evaluations;
        let mut iterations = 0_u32;
        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            for &node in nodes {
                let worst_prev = work
                    .residual
                    .iter()
                    .fold(0.0_f64, |acc, &r| acc.max(r.abs()));
                self.charge_eval(
                    self.node_inputs(node).len().max(1) as u64,
                    block,
                    mode,
                    iterations,
                    if iterations == 0 {
                        f64::NAN
                    } else {
                        worst_prev
                    },
                )?;
                if let Some(shorter) = self.propose(node, dt_s, cx)? {
                    return Ok(Err(shorter));
                }
            }

            work.group_scale.fill(0.0);
            let mut l2_sq = 0.0_f64;
            let mut linf = 0.0_f64;
            for (row, &index) in unknowns.iter().enumerate() {
                let driver = &self.drivers[index];
                let target = driver.value(&self.outputs, &self.boundary_values);
                work.fp_current[row] = self.inputs[index];
                work.fp_image[row] = target;
                let res = self.inputs[index] - target;
                work.residual[row] = res;
                l2_sq += res * res;
                linf = linf.max(res.abs());
                let group = self.input_group[index];
                let terms = driver.magnitude(&self.outputs, &self.boundary_values);
                work.group_scale[group] = work.group_scale[group]
                    .max(self.inputs[index].abs())
                    .max(terms);
            }

            let power_scale = self.block_power_scale(nodes);
            let mut worst = (0_usize, 0.0_f64);
            let mut converged = true;
            for (row, &index) in unknowns.iter().enumerate() {
                let group = self.input_group[index];
                let scale = work.group_scale[group];
                let mut allow =
                    self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * scale;
                if self.input_is_energy_rate[index] {
                    allow += ROUNDING_UNITS * f64::EPSILON * power_scale;
                }
                work.allowance[row] = allow;
                let allowed = self.config.relative_tolerance * scale + allow;
                let excess = work.residual[row].abs() - allowed;
                if !(excess <= 0.0) {
                    converged = false;
                    if !(excess <= worst.1) {
                        worst = (row, excess);
                    }
                }
            }

            if converged {
                for (row, &index) in unknowns.iter().enumerate() {
                    self.input_allowance[index] = work.allowance[row];
                }
                self.record_block_diag(BlockStepDiagnostic {
                    block,
                    mode,
                    iterations,
                    evaluations: self.step_tracker.evaluations - eval_start,
                    residual_l2: l2_sq.sqrt(),
                    residual_linf: linf,
                });
                return Ok(Ok(iterations));
            }
            if iterations >= self.config.max_iterations {
                self.record_block_diag(BlockStepDiagnostic {
                    block,
                    mode,
                    iterations,
                    evaluations: self.step_tracker.evaluations - eval_start,
                    residual_l2: l2_sq.sqrt(),
                    residual_linf: linf,
                });
                return Err(self.not_converged(
                    nodes,
                    unknowns,
                    &work.residual,
                    iterations,
                    worst.0,
                ));
            }

            let step_ok = match mode {
                BlockSolveMode::FixedPointAitken => self.candidate_controller.blocks[block]
                    .aitken_step(
                        &work.fp_current[..m],
                        &work.fp_image[..m],
                        aitken_cfg,
                        &mut work.fp_next[..m],
                    )
                    .map(|_| ()),
                BlockSolveMode::FixedPointIqnIls => self.candidate_controller.blocks[block]
                    .iqn_ils_step(
                        &work.fp_current[..m],
                        &work.fp_image[..m],
                        iqn_cfg,
                        &mut work.iqn_scratch,
                        &mut work.fp_next[..m],
                    )
                    .map(|_| ()),
                _ => unreachable!("non-fixed-point mode dispatched to fixed-point solver"),
            };
            if step_ok.is_err() {
                return Err(SystemError::Singular {
                    nodes: self.names_of(nodes),
                });
            }

            // Safeguard temperature inputs against non-positive iterates by damping
            // back toward the admissible current iterate.
            let mut alpha = 1.0_f64;
            for _ in 0..12 {
                let ok = unknowns.iter().enumerate().all(|(row, &index)| {
                    let cand =
                        work.fp_current[row] + alpha * (work.fp_next[row] - work.fp_current[row]);
                    cand.is_finite() && (!self.input_is_temperature[index] || cand > 0.0)
                });
                if ok {
                    break;
                }
                alpha *= 0.5;
            }
            for (row, &index) in unknowns.iter().enumerate() {
                let cand =
                    work.fp_current[row] + alpha * (work.fp_next[row] - work.fp_current[row]);
                if !cand.is_finite() || (self.input_is_temperature[index] && !(cand > 0.0)) {
                    let shorter = 0.5 * dt_s;
                    if shorter > 0.0 && shorter < dt_s {
                        return Ok(Err(shorter));
                    }
                    return Err(SystemError::Singular {
                        nodes: self.names_of(nodes),
                    });
                }
                self.inputs[index] = cand;
            }
            iterations += 1;
        }
    }

    /// Solve one looped block monolithically using `StageResidual` through
    /// `fs-solver`'s fallible Newton--Krylov adapter.
    ///
    /// Component refusals propagate as typed [`SystemError::Component`] values
    /// directly through `Result` and are never encoded as `NaN`.
    #[allow(clippy::too_many_lines, clippy::items_after_statements)]
    pub(super) fn solve_monolithic_block(
        &mut self,
        block: usize,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        let mode = BlockSolveMode::MonolithicNewtonKrylov;
        for &node in nodes {
            if !self.nodes[node].has_stage_residual() || self.nodes[node].stage_unknowns() == 0 {
                return Err(SystemError::Unsupported {
                    node: self.names[node].clone(),
                    needs: "StageResidual with at least one unknown for monolithic Newton-Krylov",
                });
            }
        }
        let m = unknowns.len();
        for &index in unknowns {
            let start = self.accepted_inputs[index];
            if self.input_is_temperature[index] && !(start.is_finite() && start > 0.0) {
                return Err(SystemError::GuessRequired {
                    port: self.input_label[index].clone(),
                });
            }
            self.inputs[index] = start;
        }

        // Partition monolithic state X = [z_0, ..., z_{k-1}, u_0, ..., u_{m-1}].
        let mut z_offsets = Vec::with_capacity(nodes.len() + 1);
        z_offsets.push(0_usize);
        for &node in nodes {
            let last = *z_offsets.last().unwrap();
            z_offsets.push(last + self.nodes[node].stage_unknowns());
        }
        let n_z = *z_offsets.last().unwrap();
        let total_dim = n_z + m;

        let mut x0 = vec![0.0_f64; total_dim];
        for (row, &index) in unknowns.iter().enumerate() {
            x0[n_z + row] = self.inputs[index];
        }
        let mut row_scale = vec![1.0_f64; total_dim];
        for (pos, &node) in nodes.iter().enumerate() {
            let zr = z_offsets[pos]..z_offsets[pos + 1];
            let nz_node = zr.len();
            let inputs = self.node_inputs(node);
            self.nodes[node]
                .stage_initial_guess(
                    self.accepted[node].as_ref(),
                    dt_s,
                    &self.inputs[inputs.clone()],
                    &mut x0[zr.clone()],
                )
                .map_err(|error| self.refused(node, error, cx))?;
            // Estimate natural row scale per internal unknown using a unit `jvp_z` probe
            // so left-diagonal scaling `D^{-1} F(X) = 0` normalizes mixed physical units
            // (e.g. 300 K thermal capacity beside 1 m/s mechanics) cleanly to O(1).
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            let ones = vec![1.0_f64; nz_node];
            let mut dr_probe = vec![0.0_f64; nz_node];
            let mut do_probe = vec![0.0_f64; self.node_outputs(node).len()];
            match self.nodes[node].stage_jvp_z(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs],
                &x0[zr.clone()],
                &ones,
                &mut dr_probe,
                &mut do_probe,
                cx,
            ) {
                Ok(()) => {
                    for k in 0..nz_node {
                        let s = (dr_probe[k].abs().max(1.0)) * x0[zr.start + k].abs().max(1.0);
                        if s.is_finite() && s > 0.0 {
                            row_scale[zr.start + k] = s;
                        }
                    }
                }
                Err(_) if cx.checkpoint().is_err() => return Err(SystemError::Cancelled),
                Err(_) => {}
            }
        }
        for row in 0..m {
            let s = x0[n_z + row].abs().max(1.0);
            if s.is_finite() && s > 0.0 {
                row_scale[n_z + row] = s;
            }
        }

        struct MonolithicScratch {
            inputs: Vec<f64>,
            outputs: Vec<f64>,
            d_outputs: Vec<f64>,
            plus_r: Vec<f64>,
            minus_r: Vec<f64>,
            plus_o: Vec<f64>,
            minus_o: Vec<f64>,
            grad: Vec<(usize, f64)>,
            evals: u64,
            work_units: u64,
            last_linf: f64,
        }

        struct MonolithicProblem<'a> {
            system: &'a System,
            block: usize,
            nodes: &'a [usize],
            unknowns: &'a [usize],
            z_offsets: &'a [usize],
            row_scale: &'a [f64],
            n_z: usize,
            dt_s: f64,
            base_evals: u64,
            base_work: u64,
            scratch: RefCell<MonolithicScratch>,
        }

        impl MonolithicProblem<'_> {
            fn charge(&self, weight: u64) -> Result<(), SystemError> {
                let mut sc = self.scratch.borrow_mut();
                sc.evals = sc.evals.saturating_add(1);
                sc.work_units = sc.work_units.saturating_add(weight.max(1));
                let total_evals = self.base_evals.saturating_add(sc.evals);
                let total_work = self.base_work.saturating_add(sc.work_units);
                self.system.check_budgets_at(
                    total_evals,
                    total_work,
                    self.block,
                    BlockSolveMode::MonolithicNewtonKrylov,
                    0,
                    sc.last_linf,
                )
            }
        }

        impl FallibleNonlinearProblem for MonolithicProblem<'_> {
            type Error = SystemError;

            fn dimension(&self) -> usize {
                self.n_z + self.unknowns.len()
            }

            fn residual(
                &self,
                x: &[f64],
                residual: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), SystemError> {
                let u = &x[self.n_z..];
                for (row, &index) in self.unknowns.iter().enumerate() {
                    if self.system.input_is_temperature[index] && !(u[row] > 0.0) {
                        return Err(SystemError::Singular {
                            nodes: self.system.names_of(self.nodes),
                        });
                    }
                }
                {
                    let mut sc = self.scratch.borrow_mut();
                    for (row, &index) in self.unknowns.iter().enumerate() {
                        sc.inputs[index] = u[row];
                    }
                }
                for (pos, &node) in self.nodes.iter().enumerate() {
                    let zr = self.z_offsets[pos]..self.z_offsets[pos + 1];
                    self.charge((zr.len() + self.system.node_inputs(node).len()) as u64)?;
                    let inp_r = self.system.node_inputs(node);
                    let out_r = self.system.node_outputs(node);
                    let mut sc = self.scratch.borrow_mut();
                    let MonolithicScratch {
                        inputs, outputs, ..
                    } = &mut *sc;
                    self.system.nodes[node]
                        .stage_residual(
                            self.system.accepted[node].as_ref(),
                            self.dt_s,
                            &inputs[inp_r],
                            &x[zr.clone()],
                            &mut residual[zr],
                            &mut outputs[out_r],
                            cx,
                        )
                        .map_err(|error| self.system.refused(node, error, cx))?;
                }
                let mut sc = self.scratch.borrow_mut();
                let mut linf = 0.0_f64;
                for (row, &index) in self.unknowns.iter().enumerate() {
                    let target =
                        self.system.drivers[index].value(&sc.outputs, &self.system.boundary_values);
                    let r = u[row] - target;
                    residual[self.n_z + row] = r;
                    linf = linf.max(r.abs());
                }
                for (r_val, &scale) in residual.iter_mut().zip(self.row_scale) {
                    *r_val /= scale;
                }
                sc.last_linf = linf;
                Ok(())
            }

            fn jacobian_apply(
                &self,
                x: &[f64],
                direction: &[f64],
                output: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), SystemError> {
                let u = &x[self.n_z..];
                let du = &direction[self.n_z..];
                {
                    let mut sc = self.scratch.borrow_mut();
                    for (row, &index) in self.unknowns.iter().enumerate() {
                        sc.inputs[index] = u[row];
                    }
                    sc.d_outputs.fill(0.0);
                }

                for (pos, &node) in self.nodes.iter().enumerate() {
                    let zr = self.z_offsets[pos]..self.z_offsets[pos + 1];
                    let nz_node = zr.len();
                    let inp_r = self.system.node_inputs(node);
                    let out_r = self.system.node_outputs(node);
                    let n_out = out_r.len();

                    // 1. Exact internal derivative (dr/dz * dz, d_outputs/dz * dz) via StageResidual::jvp_z.
                    self.charge((nz_node + inp_r.len()) as u64)?;
                    {
                        let mut sc = self.scratch.borrow_mut();
                        let MonolithicScratch {
                            inputs,
                            outputs,
                            d_outputs,
                            plus_r,
                            ..
                        } = &mut *sc;
                        // Base outputs needed for Driver::gradient on nonlinear drivers.
                        plus_r.resize(nz_node, 0.0);
                        self.system.nodes[node]
                            .stage_residual(
                                self.system.accepted[node].as_ref(),
                                self.dt_s,
                                &inputs[inp_r.clone()],
                                &x[zr.clone()],
                                &mut plus_r[..nz_node],
                                &mut outputs[out_r.clone()],
                                cx,
                            )
                            .map_err(|error| self.system.refused(node, error, cx))?;
                        self.system.nodes[node]
                            .stage_jvp_z(
                                self.system.accepted[node].as_ref(),
                                self.dt_s,
                                &inputs[inp_r.clone()],
                                &x[zr.clone()],
                                &direction[zr.clone()],
                                &mut output[zr.clone()],
                                &mut d_outputs[out_r.clone()],
                                cx,
                            )
                            .map_err(|error| self.system.refused(node, error, cx))?;
                    }

                    // 2. Exact central directional derivative along du_local (when nonzero).
                    let mut du_norm2 = 0.0_f64;
                    let mut u_scale = 1.0_f64;
                    for (row, &index) in self.unknowns.iter().enumerate() {
                        if inp_r.contains(&index) {
                            du_norm2 += du[row] * du[row];
                            u_scale = u_scale.max(u[row].abs());
                        }
                    }
                    if du_norm2 > 0.0 {
                        let du_norm = du_norm2.sqrt();
                        let mut eps = 1.0e-6 * u_scale / du_norm;
                        // Keep temperature perturbations strictly positive.
                        for (row, &index) in self.unknowns.iter().enumerate() {
                            if inp_r.contains(&index) && self.system.input_is_temperature[index] {
                                let step_t = (eps * du[row]).abs();
                                if step_t >= 0.25 * u[row] && du[row] != 0.0 {
                                    eps = eps.min(0.25 * u[row] / du[row].abs());
                                }
                            }
                        }
                        if eps > 0.0 {
                            self.charge((2 * nz_node) as u64)?;
                            let mut sc = self.scratch.borrow_mut();
                            sc.plus_r.resize(nz_node, 0.0);
                            sc.minus_r.resize(nz_node, 0.0);
                            sc.plus_o.resize(n_out, 0.0);
                            sc.minus_o.resize(n_out, 0.0);

                            // +eps
                            for (row, &index) in self.unknowns.iter().enumerate() {
                                if inp_r.contains(&index) {
                                    sc.inputs[index] = u[row] + eps * du[row];
                                }
                            }
                            {
                                let MonolithicScratch {
                                    inputs,
                                    plus_r,
                                    plus_o,
                                    ..
                                } = &mut *sc;
                                self.system.nodes[node]
                                    .stage_residual(
                                        self.system.accepted[node].as_ref(),
                                        self.dt_s,
                                        &inputs[inp_r.clone()],
                                        &x[zr.clone()],
                                        &mut plus_r[..nz_node],
                                        &mut plus_o[..n_out],
                                        cx,
                                    )
                                    .map_err(|error| self.system.refused(node, error, cx))?;
                            }

                            // -eps
                            for (row, &index) in self.unknowns.iter().enumerate() {
                                if inp_r.contains(&index) {
                                    sc.inputs[index] = u[row] - eps * du[row];
                                }
                            }
                            {
                                let MonolithicScratch {
                                    inputs,
                                    minus_r,
                                    minus_o,
                                    ..
                                } = &mut *sc;
                                self.system.nodes[node]
                                    .stage_residual(
                                        self.system.accepted[node].as_ref(),
                                        self.dt_s,
                                        &inputs[inp_r.clone()],
                                        &x[zr.clone()],
                                        &mut minus_r[..nz_node],
                                        &mut minus_o[..n_out],
                                        cx,
                                    )
                                    .map_err(|error| self.system.refused(node, error, cx))?;
                            }

                            // Restore base inputs and accumulate central directional derivative.
                            for (row, &index) in self.unknowns.iter().enumerate() {
                                if inp_r.contains(&index) {
                                    sc.inputs[index] = u[row];
                                }
                            }
                            let inv_2eps = 0.5 / eps;
                            for k in 0..nz_node {
                                output[zr.start + k] += (sc.plus_r[k] - sc.minus_r[k]) * inv_2eps;
                            }
                            for k in 0..n_out {
                                sc.d_outputs[out_r.start + k] +=
                                    (sc.plus_o[k] - sc.minus_o[k]) * inv_2eps;
                            }
                        }
                    }
                }

                // 3. Interface equations: d(u - driver(outputs)) = du - (d driver / d outputs) d_outputs.
                let mut sc = self.scratch.borrow_mut();
                let MonolithicScratch {
                    outputs,
                    d_outputs,
                    grad,
                    ..
                } = &mut *sc;
                for (row, &index) in self.unknowns.iter().enumerate() {
                    self.system.drivers[index].gradient(outputs, grad);
                    let mut d_target = 0.0_f64;
                    for &(src, slope) in grad.iter() {
                        d_target += slope * d_outputs[src];
                    }
                    output[self.n_z + row] = du[row] - d_target;
                }
                for (o_val, &scale) in output.iter_mut().zip(self.row_scale) {
                    *o_val /= scale;
                }
                Ok(())
            }
        }

        let mono_cfg = self.monolithic_config();
        let nk_cfg = NewtonKrylovConfig {
            absolute_tolerance: mono_cfg
                .absolute_tolerance
                .max(self.config.absolute_tolerance),
            relative_tolerance: mono_cfg.relative_tolerance,
            linear_restart: mono_cfg.linear_restart.max(total_dim).min(64),
            max_linear_cycles: mono_cfg.max_linear_cycles,
            forcing_minimum: 1.0e-13,
            forcing_maximum: 1.0e-6,
            ..NewtonKrylovConfig::default()
        };

        let eval_start = self.step_tracker.evaluations;
        let base_evals = self.step_tracker.evaluations;
        let base_work = self.step_tracker.work_units;
        let problem = MonolithicProblem {
            system: self,
            block,
            nodes,
            unknowns,
            z_offsets: &z_offsets,
            row_scale: &row_scale,
            n_z,
            dt_s,
            base_evals,
            base_work,
            scratch: RefCell::new(MonolithicScratch {
                inputs: self.inputs.clone(),
                outputs: self.outputs.clone(),
                d_outputs: vec![0.0; self.outputs.len()],
                plus_r: Vec::new(),
                minus_r: Vec::new(),
                plus_o: Vec::new(),
                minus_o: Vec::new(),
                grad: Vec::new(),
                evals: 0,
                work_units: 0,
                last_linf: f64::NAN,
            }),
        };

        let solved = solve_fallible_newton_krylov(
            &problem,
            x0,
            nk_cfg,
            mono_cfg
                .max_iterations
                .min(self.config.max_iterations as usize),
            cx,
        );
        let (spent_evals, spent_work) = {
            let sc = problem.scratch.borrow();
            (sc.evals, sc.work_units)
        };
        drop(problem);
        self.step_tracker.evaluations = self.step_tracker.evaluations.saturating_add(spent_evals);
        self.step_tracker.work_units = self.step_tracker.work_units.saturating_add(spent_work);

        let sol = match solved {
            Ok(sol) => sol,
            Err(FallibleNewtonError::Cancelled) => return Err(SystemError::Cancelled),
            Err(FallibleNewtonError::Problem(err)) => return Err(err),
            Err(_) => {
                return Err(SystemError::Singular {
                    nodes: self.names_of(nodes),
                });
            }
        };

        let iters_u32 = u32::try_from(sol.report.iterations).unwrap_or(u32::MAX);
        if !sol.report.converged {
            let mut worst_row = 0usize;
            let mut worst_abs = -1.0_f64;
            for row in 0..m {
                let r = sol.residual[n_z + row] * row_scale[n_z + row];
                work.residual[row] = r;
                if r.abs() > worst_abs {
                    worst_abs = r.abs();
                    worst_row = row;
                }
            }
            let l2 = sol
                .residual
                .iter()
                .zip(&row_scale)
                .map(|(&v, &s)| {
                    let u = v * s;
                    u * u
                })
                .sum::<f64>()
                .sqrt();
            let linf = sol
                .residual
                .iter()
                .zip(&row_scale)
                .fold(0.0_f64, |a, (&v, &s)| a.max((v * s).abs()));
            self.record_block_diag(BlockStepDiagnostic {
                block,
                mode,
                iterations: iters_u32,
                evaluations: self.step_tracker.evaluations - eval_start,
                residual_l2: l2,
                residual_linf: linf,
            });
            return Err(self.not_converged(nodes, unknowns, &work.residual, iters_u32, worst_row));
        }

        // Publish converged u into trial inputs, evaluate exact final outputs,
        // and call `stage_accept` on every node to populate candidate & receipts.
        for (row, &index) in unknowns.iter().enumerate() {
            self.inputs[index] = sol.x[n_z + row];
        }
        let mut r_tmp = Vec::new();
        for (pos, &node) in nodes.iter().enumerate() {
            let zr = z_offsets[pos]..z_offsets[pos + 1];
            let inp_r = self.node_inputs(node);
            let out_r = self.node_outputs(node);
            r_tmp.resize(zr.len(), 0.0);
            self.nodes[node]
                .stage_residual(
                    self.accepted[node].as_ref(),
                    dt_s,
                    &self.inputs[inp_r.clone()],
                    &sol.x[zr.clone()],
                    &mut r_tmp,
                    &mut self.outputs[out_r],
                    cx,
                )
                .map_err(|error| self.refused(node, error, cx))?;
            let receipt = self.nodes[node]
                .stage_accept(
                    self.accepted[node].as_ref(),
                    dt_s,
                    &self.inputs[inp_r],
                    &sol.x[zr],
                    self.candidate[node].as_mut(),
                )
                .map_err(|error| self.refused(node, error, cx))?;
            self.receipts[node] = receipt;
        }

        // Compute interface allowances for the energy audit.
        work.group_scale.resize(self.group_count, 0.0);
        work.group_scale.fill(0.0);
        let mut l2_sq = 0.0_f64;
        let mut linf = 0.0_f64;
        for (row, &index) in unknowns.iter().enumerate() {
            let driver = &self.drivers[index];
            let target = driver.value(&self.outputs, &self.boundary_values);
            let res = self.inputs[index] - target;
            work.residual[row] = res;
            l2_sq += res * res;
            linf = linf.max(res.abs());
            let group = self.input_group[index];
            let terms = driver.magnitude(&self.outputs, &self.boundary_values);
            work.group_scale[group] = work.group_scale[group]
                .max(self.inputs[index].abs())
                .max(terms);
        }
        let power_scale = self.block_power_scale(nodes);
        for &index in unknowns {
            let scale = work.group_scale[self.input_group[index]];
            let mut allow = self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * scale;
            if self.input_is_energy_rate[index] {
                allow += ROUNDING_UNITS * f64::EPSILON * power_scale;
            }
            self.input_allowance[index] = allow;
        }
        self.record_block_diag(BlockStepDiagnostic {
            block,
            mode,
            iterations: iters_u32,
            evaluations: self.step_tracker.evaluations - eval_start,
            residual_l2: l2_sq.sqrt(),
            residual_linf: linf,
        });
        Ok(Ok(iters_u32))
    }
}
