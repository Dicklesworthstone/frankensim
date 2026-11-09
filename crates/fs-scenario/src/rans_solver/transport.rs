//! Physical-time, constant-property 1D Launder--Sharma channel transport.
//!
//! The transported dissipation is modified epsilon: the physical k sink is
//! epsilon + D, D=2 nu |d sqrt(k)/dy|^2. Its equation includes the additive
//! E=2 nu nu_t |d²U/dy²|² term. See NASA TM-110446 Part B §3 and the original
//! Launder--Sharma closure (Ceps1=1.44, as on the admitted model card).
//!
//! One first-order segregated step solves momentum, k and modified epsilon
//! with implicit diffusion/destruction and lagged constitutive coefficients.
//! Face eddy work from the accepted momentum trial supplies exactly the same
//! integrated production to k. The optional thermal step transports sensible
//! heat with the prepared endpoint eddy viscosity and can deposit molecular
//! and modeled turbulent dissipation. All fields are accepted together.
//! Neither path supplies porous, buoyancy, transition or
//! three-dimensional flow physics.

use super::{RansChannelGrid, RansFieldState, RansSolver};

mod thermal;
pub use thermal::{
    RansConjugateHeatLoad, RansConjugateReport, RansConjugateWall, RansConjugateWorkspace,
    RansHeatCoupling, RansHeatLoad, RansThermalTransportReport, RansWallLayer,
};
mod adaptive;
pub use adaptive::{
    RansAdaptiveControls, RansAdaptiveError, RansAdaptiveHeatBalance, RansAdaptiveReport,
    RansAdaptiveWorkspace, RansFieldTolerances,
};

/// A refused step leaves every accepted field untouched.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RansTransportError {
    /// Invalid geometry, coefficients, state or physical time interval.
    InvalidInput(&'static str),
    /// The declared workspace node budget is too small.
    NodeBudget {
        /// Nodes required by the requested grid.
        needed: usize,
        /// Caller-declared node ceiling.
        maximum: usize,
    },
    /// Numerical buffers could not be allocated.
    Allocation,
    /// An actual coefficient, pivot, diagnostic or update is unrepresentable.
    Arithmetic {
        /// Equation or diagnostic that could not be represented.
        equation: &'static str,
        /// Affected wall-normal node.
        node: usize,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

impl core::fmt::Display for RansTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for RansTransportError {}

/// Reusable numerical storage: nine f64 fields, no allocation in a step.
pub struct RansTransportWorkspace {
    velocity: Vec<f64>,
    kinetic: Vec<f64>,
    dissipation: Vec<f64>,
    viscosity: Vec<f64>,
    production: Vec<f64>,
    lower: Vec<f64>,
    diagonal: Vec<f64>,
    upper: Vec<f64>,
    rhs: Vec<f64>,
}

impl RansTransportWorkspace {
    /// Allocate bounded storage for the grid's nodes, including both boundaries.
    pub fn try_new(
        grid: &RansChannelGrid,
        maximum_nodes: usize,
    ) -> Result<Self, RansTransportError> {
        let n = grid
            .num_cells
            .checked_add(1)
            .ok_or(RansTransportError::InvalidInput("node count overflow"))?;
        if n > maximum_nodes {
            return Err(RansTransportError::NodeBudget {
                needed: n,
                maximum: maximum_nodes,
            });
        }
        n.checked_mul(9 * core::mem::size_of::<f64>())
            .ok_or(RansTransportError::Allocation)?;
        let field = || {
            let mut values = Vec::new();
            values
                .try_reserve_exact(n)
                .map_err(|_| RansTransportError::Allocation)?;
            values.resize(n, 0.0);
            Ok::<_, RansTransportError>(values)
        };
        Ok(Self {
            velocity: field()?,
            kinetic: field()?,
            dissipation: field()?,
            viscosity: field()?,
            production: field()?,
            lower: field()?,
            diagonal: field()?,
            upper: field()?,
            rhs: field()?,
        })
    }

    /// Numerical payload, excluding Vec headers and allocator bookkeeping.
    #[must_use]
    pub fn numerical_bytes(&self) -> usize {
        self.velocity.len() * 9 * core::mem::size_of::<f64>()
    }

    fn publish_mechanics(&self, state: &mut RansFieldState) {
        state.u.copy_from_slice(&self.velocity);
        state.k.copy_from_slice(&self.kinetic);
        state.eps.copy_from_slice(&self.dissipation);
        state.nu_t.copy_from_slice(&self.viscosity);
    }

    fn solve(
        &mut self,
        equation: &'static str,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RansTransportError> {
        let n = self.rhs.len();
        for i in 1..n {
            checkpoint(i, cancelled)?;
            let factor = self.lower[i] / self.diagonal[i - 1];
            self.diagonal[i] -= factor * self.upper[i - 1];
            self.rhs[i] -= factor * self.rhs[i - 1];
            if !self.diagonal[i].is_finite() || self.diagonal[i] <= 0.0 || !self.rhs[i].is_finite()
            {
                return Err(RansTransportError::Arithmetic { equation, node: i });
            }
        }
        self.rhs[n - 1] /= self.diagonal[n - 1];
        for i in (0..n - 1).rev() {
            checkpoint(i, cancelled)?;
            self.rhs[i] = (self.rhs[i] - self.upper[i] * self.rhs[i + 1]) / self.diagonal[i];
        }
        for (i, value) in self.rhs.iter().enumerate() {
            checkpoint(i, cancelled)?;
            if !value.is_finite() {
                return Err(RansTransportError::Arithmetic { equation, node: i });
            }
        }
        Ok(())
    }
}

/// Discrete mechanical budgets per unit density and unit channel span.
/// Energies have units m³/s²; rates and balance defects have units m³/s³.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansTransportReport {
    /// Common physical interval advanced by all three equations [s].
    pub time_step_s: f64,
    /// Discrete integral of U²/2 before the step.
    pub mean_energy_before: f64,
    /// Discrete integral of U²/2 at the accepted endpoint.
    pub mean_energy_after: f64,
    /// Discrete integral of k before the step.
    pub turbulent_energy_before: f64,
    /// Discrete integral of k at the accepted endpoint.
    pub turbulent_energy_after: f64,
    /// Pressure acceleration times accepted mean volume flux.
    pub pressure_power: f64,
    /// Molecular mean-flow dissipation from shared face gradients.
    pub molecular_dissipation: f64,
    /// Mean-flow eddy loss transferred to the k equation exactly once.
    pub turbulence_production: f64,
    /// Patankar destruction actually used by this first-order k step.
    pub turbulent_destruction: f64,
    /// Diffusive k loss into the prescribed zero-k wall.
    pub turbulent_wall_loss: f64,
    /// Backward-Euler mean-velocity numerical dissipation, always nonnegative.
    pub temporal_dissipation: f64,
    /// Mean-energy rate minus pressure power net of all discrete mean losses.
    pub mean_balance_defect: f64,
    /// Turbulent-energy rate minus production net of destruction and wall loss.
    pub turbulent_balance_defect: f64,
}

impl RansSolver {
    /// Advance mean momentum, k and modified epsilon over one physical interval.
    ///
    /// A stationary wall has U=k=epsilon=0. The centerline has zero normal
    /// flux and its own half control volume; its values are solved, not copied
    /// from the adjacent node. Pressure gradient is held over the interval.
    /// Material coefficients are constant, with no wall-normal mean advection.
    ///
    /// Admission requires either strictly positive k/epsilon at every fluid
    /// node, or the exact globally laminar k=epsilon=0 state. Mixed zero cells
    /// make the low-Re destruction ratios singular and are explicitly refused.
    /// There are no arbitrary turbulence floors or automatically seeded eddies.
    ///
    /// Momentum uses old eddy viscosity. Its accepted face eddy work produces
    /// k; destruction is linearized with old k/epsilon and treated implicitly.
    /// E uses the new velocity curvature and old eddy viscosity. Updated k and
    /// epsilon then determine the returned eddy viscosity. This is first-order
    /// time integration, not a fully converged nonlinear endpoint solve.
    ///
    /// The workspace can be reused. Cancellation is polled at most 64 nodes
    /// apart and once before joint publication. No temperature is advanced.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_launder_sharma(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        pressure_gradient_pa_m: f64,
        time_step_s: f64,
        work: &mut RansTransportWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansTransportReport, RansTransportError> {
        let report = self.prepare_launder_sharma(
            grid,
            state,
            pressure_gradient_pa_m,
            time_step_s,
            false,
            work,
            &mut cancelled,
        )?;
        work.publish_mechanics(state);
        Ok(report)
    }

    // Retain the original numerical path and cancellation order, but defer
    // publication so heat can share the same acceptance boundary. The optional
    // local dissipation account reuses production after the epsilon solve.
    #[allow(clippy::too_many_arguments)]
    fn prepare_launder_sharma(
        &self,
        grid: &RansChannelGrid,
        state: &RansFieldState,
        pressure_gradient_pa_m: f64,
        time_step_s: f64,
        retain_dissipation_heat: bool,
        work: &mut RansTransportWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansTransportReport, RansTransportError> {
        checkpoint(0, &mut cancelled)?;
        let n = grid.num_cells.checked_add(1).filter(|&n| n >= 5).ok_or(
            RansTransportError::InvalidInput("channel needs at least four intervals"),
        )?;
        if work.velocity.len() != n
            || grid.y_nodes.len() != n
            || [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k]
                .iter()
                .any(|v| v.len() != n)
        {
            return Err(RansTransportError::InvalidInput(
                "grid/state/workspace shape mismatch",
            ));
        }
        if [
            self.nu_laminar,
            self.rho,
            self.coeffs.c_eps_1,
            self.coeffs.c_eps_2,
            self.coeffs.sigma_k,
            self.coeffs.sigma_eps,
            time_step_s,
            grid.half_height_m,
        ]
        .iter()
        .any(|&x| !x.is_finite() || x <= 0.0)
            || !self.coeffs.c_mu.is_finite()
            || self.coeffs.c_mu < 0.0
            || !pressure_gradient_pa_m.is_finite()
            || grid.y_nodes[0] != 0.0
            || (grid.y_nodes[n - 1] - grid.half_height_m).abs()
                > 16.0 * f64::EPSILON * grid.half_height_m
            || state.u[0] != 0.0
            || state.k[0] != 0.0
            || state.eps[0] != 0.0
        {
            return Err(RansTransportError::InvalidInput(
                "coefficients, interval or wall trace",
            ));
        }
        let mut laminar = true;
        let mut positive = true;
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            if !grid.y_nodes[i].is_finite()
                || (i > 0 && grid.y_nodes[i] <= grid.y_nodes[i - 1])
                || !state.u[i].is_finite()
                || !state.k[i].is_finite()
                || !state.eps[i].is_finite()
                || !state.nu_t[i].is_finite()
                || state.nu_t[i] < 0.0
                || !state.temp_k[i].is_finite()
                || state.temp_k[i] <= 0.0
            {
                return Err(RansTransportError::InvalidInput(
                    "nonfinite or singular field/geometry",
                ));
            }
            if i > 0 {
                laminar &= state.k[i] == 0.0 && state.eps[i] == 0.0;
                positive &= state.k[i] > 0.0 && state.eps[i] > 0.0;
            }
        }
        if !laminar && !positive {
            return Err(RansTransportError::InvalidInput(
                "mixed zero turbulence is singular",
            ));
        }
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            work.viscosity[i] = eddy_viscosity(self, state.k[i], state.eps[i], i)?;
            work.production[i] = 0.0;
        }
        let acceleration = -pressure_gradient_pa_m / self.rho;
        if !acceleration.is_finite() {
            return Err(RansTransportError::InvalidInput(
                "pressure acceleration overflow",
            ));
        }
        let dt = time_step_s;
        // Backward Euler with shared symmetric face stress and a no-flux center.
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let (lower, upper) = diffusion(grid, &work.viscosity, self.nu_laminar, 1.0, i);
            row(
                work,
                i,
                dt,
                lower,
                upper,
                0.0,
                state.u[i],
                acceleration,
                "momentum",
            )?;
        }
        work.solve("momentum", &mut cancelled)?;
        work.velocity.copy_from_slice(&work.rhs);

        // Each positive face eddy-work contribution is transferred once to k.
        // The wall node has no evolved k inventory, so its face belongs to node1.
        let mut molecular_dissipation = 0.0;
        let mut turbulence_production = 0.0;
        for i in 1..n {
            checkpoint(i, &mut cancelled)?;
            let dy = grid.y_nodes[i] - grid.y_nodes[i - 1];
            let du = work.velocity[i] - work.velocity[i - 1];
            let gradient_work = du * (du / dy);
            let eddy = 0.5 * work.viscosity[i - 1] + 0.5 * work.viscosity[i];
            let production = eddy * gradient_work;
            molecular_dissipation += self.nu_laminar * gradient_work;
            turbulence_production += production;
            if i == 1 {
                work.production[i] += production;
            } else {
                work.production[i - 1] += 0.5 * production;
                work.production[i] += 0.5 * production;
            }
        }
        for i in 1..n {
            checkpoint(i, &mut cancelled)?;
            work.production[i] /= volume(grid, i);
        }

        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let (d, _) = low_re_terms(
                grid,
                &state.k,
                &work.velocity,
                self.nu_laminar,
                &work.viscosity,
                i,
            );
            let sink = if laminar || i == 0 {
                0.0
            } else {
                (state.eps[i] + d) / state.k[i]
            };
            let (lower, upper) = diffusion(
                grid,
                &work.viscosity,
                self.nu_laminar,
                self.coeffs.sigma_k,
                i,
            );
            let production = work.production[i];
            row(work, i, dt, lower, upper, sink, state.k[i], production, "k")?;
        }
        work.solve("k", &mut cancelled)?;
        work.kinetic.copy_from_slice(&work.rhs);

        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let (_, e) = low_re_terms(
                grid,
                &state.k,
                &work.velocity,
                self.nu_laminar,
                &work.viscosity,
                i,
            );
            let (sink, source) = if laminar || i == 0 {
                (0.0, 0.0)
            } else {
                let re_t = (state.k[i] * (state.k[i] / state.eps[i])) / self.nu_laminar;
                let f2 = 1.0 - 0.3 * (-(re_t * re_t)).exp();
                let rate = state.eps[i] / state.k[i];
                (
                    self.coeffs.c_eps_2 * f2 * rate,
                    self.coeffs.c_eps_1 * work.production[i] * rate + e,
                )
            };
            let (lower, upper) = diffusion(
                grid,
                &work.viscosity,
                self.nu_laminar,
                self.coeffs.sigma_eps,
                i,
            );
            row(
                work,
                i,
                dt,
                lower,
                upper,
                sink,
                state.eps[i],
                source,
                "epsilon",
            )?;
        }
        work.solve("epsilon", &mut cancelled)?;
        work.dissipation.copy_from_slice(&work.rhs);

        let mut report = RansTransportReport {
            time_step_s: dt,
            mean_energy_before: 0.0,
            mean_energy_after: 0.0,
            turbulent_energy_before: 0.0,
            turbulent_energy_after: 0.0,
            pressure_power: 0.0,
            molecular_dissipation,
            turbulence_production,
            turbulent_destruction: 0.0,
            turbulent_wall_loss: 0.0,
            temporal_dissipation: 0.0,
            mean_balance_defect: 0.0,
            turbulent_balance_defect: 0.0,
        };
        let mut mean_energy_rate = 0.0;
        let mut turbulent_energy_rate = 0.0;
        if retain_dissipation_heat {
            work.production[0] = 0.0;
        }
        for i in 1..n {
            checkpoint(i, &mut cancelled)?;
            if work.kinetic[i] < 0.0
                || work.dissipation[i] < 0.0
                || (!laminar && (work.kinetic[i] == 0.0 || work.dissipation[i] == 0.0))
            {
                return Err(RansTransportError::Arithmetic {
                    equation: "positive turbulence",
                    node: i,
                });
            }
            let v = volume(grid, i);
            report.mean_energy_before += 0.5 * v * state.u[i] * state.u[i];
            report.mean_energy_after += 0.5 * v * work.velocity[i] * work.velocity[i];
            report.turbulent_energy_before += v * state.k[i];
            report.turbulent_energy_after += v * work.kinetic[i];
            report.pressure_power += v * acceleration * work.velocity[i];
            let du = work.velocity[i] - state.u[i];
            report.temporal_dissipation += 0.5 * v * du * (du / dt);
            mean_energy_rate += 0.5 * v * (work.velocity[i] + state.u[i]) * (du / dt);
            turbulent_energy_rate += v * ((work.kinetic[i] - state.k[i]) / dt);
            let destruction = if !laminar {
                let (d, _) = low_re_terms(
                    grid,
                    &state.k,
                    &work.velocity,
                    self.nu_laminar,
                    &work.viscosity,
                    i,
                );
                v * ((state.eps[i] + d) / state.k[i]) * work.kinetic[i]
            } else {
                0.0
            };
            report.turbulent_destruction += destruction;
            if retain_dissipation_heat {
                // Exact Patankar loss used by the accepted k equation, per
                // unit density and integrated over this mechanical dual cell.
                work.production[i] = destruction;
            }
        }
        report.turbulent_wall_loss = (self.nu_laminar
            + 0.5 * work.viscosity[1] / self.coeffs.sigma_k)
            * (work.kinetic[1] / (grid.y_nodes[1] - grid.y_nodes[0]));
        report.mean_balance_defect = mean_energy_rate
            - (report.pressure_power
                - report.molecular_dissipation
                - report.turbulence_production
                - report.temporal_dissipation);
        report.turbulent_balance_defect = turbulent_energy_rate
            - (report.turbulence_production
                - report.turbulent_destruction
                - report.turbulent_wall_loss);
        if [
            report.mean_energy_before,
            report.mean_energy_after,
            report.turbulent_energy_before,
            report.turbulent_energy_after,
            report.pressure_power,
            report.molecular_dissipation,
            report.turbulence_production,
            report.turbulent_destruction,
            report.turbulent_wall_loss,
            report.temporal_dissipation,
            report.mean_balance_defect,
            report.turbulent_balance_defect,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(RansTransportError::Arithmetic {
                equation: "energy budget",
                node: n - 1,
            });
        }
        if retain_dissipation_heat {
            // Molecular face loss is distributed to the two thermal half
            // cells, including the thermal wall inventory. Eddy work already
            // went to k; neither it, temporal loss nor k flux into the wall is
            // deposited a second time. No additional numerical buffer is used.
            for i in 1..n {
                checkpoint(i, &mut cancelled)?;
                let dy = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let du = work.velocity[i] - work.velocity[i - 1];
                let half = 0.5 * (self.nu_laminar * (du * (du / dy)));
                work.production[i - 1] += half;
                work.production[i] += half;
                if !work.production[i - 1].is_finite() || !work.production[i].is_finite() {
                    return Err(RansTransportError::Arithmetic {
                        equation: "local dissipation heat",
                        node: i,
                    });
                }
            }
        }
        // Finish every fallible operation before publishing the joint endpoint.
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            work.viscosity[i] = eddy_viscosity(self, work.kinetic[i], work.dissipation[i], i)?;
        }
        checkpoint(0, &mut cancelled)?;
        Ok(report)
    }
}

fn checkpoint(i: usize, cancelled: &mut impl FnMut() -> bool) -> Result<(), RansTransportError> {
    if i % 64 == 0 && cancelled() {
        Err(RansTransportError::Cancelled)
    } else {
        Ok(())
    }
}

fn volume(grid: &RansChannelGrid, i: usize) -> f64 {
    let left = grid.y_nodes[i] - grid.y_nodes[i - 1];
    let right = if i == grid.num_cells {
        0.0
    } else {
        grid.y_nodes[i + 1] - grid.y_nodes[i]
    };
    0.5 * (left + right)
}

fn diffusion(grid: &RansChannelGrid, nut: &[f64], nu: f64, sigma: f64, i: usize) -> (f64, f64) {
    if i == 0 {
        return (0.0, 0.0);
    }
    let v = volume(grid, i);
    let lower = ((nu + (0.5 * nut[i - 1] + 0.5 * nut[i]) / sigma)
        / (grid.y_nodes[i] - grid.y_nodes[i - 1]))
        / v;
    let upper = if i == grid.num_cells {
        0.0
    } else {
        ((nu + (0.5 * nut[i] + 0.5 * nut[i + 1]) / sigma) / (grid.y_nodes[i + 1] - grid.y_nodes[i]))
            / v
    };
    (lower, upper)
}

#[allow(clippy::too_many_arguments)]
fn row(
    work: &mut RansTransportWorkspace,
    i: usize,
    dt: f64,
    lower: f64,
    upper: f64,
    sink: f64,
    old: f64,
    source: f64,
    equation: &'static str,
) -> Result<(), RansTransportError> {
    let (lo, diag, hi, rhs) = if i == 0 {
        (0.0, 1.0, 0.0, 0.0)
    } else {
        (
            -dt * lower,
            1.0 + dt * (lower + upper + sink),
            -dt * upper,
            old + dt * source,
        )
    };
    if [lo, diag, hi, rhs, sink].iter().any(|x| !x.is_finite())
        || diag <= 0.0
        || sink < 0.0
        || lower < 0.0
        || upper < 0.0
    {
        return Err(RansTransportError::Arithmetic { equation, node: i });
    }
    work.lower[i] = lo;
    work.diagonal[i] = diag;
    work.upper[i] = hi;
    work.rhs[i] = rhs;
    Ok(())
}

fn eddy_viscosity(
    solver: &RansSolver,
    k: f64,
    eps: f64,
    i: usize,
) -> Result<f64, RansTransportError> {
    if k == 0.0 || solver.coeffs.c_mu == 0.0 {
        return Ok(0.0);
    }
    let ratio = k * (k / eps);
    let re_t = ratio / solver.nu_laminar;
    let inverse = 1.0 / (1.0 + re_t / 50.0);
    let value = solver.coeffs.c_mu * (-3.4 * inverse * inverse).exp() * ratio;
    if !ratio.is_finite() || !value.is_finite() || value < 0.0 {
        Err(RansTransportError::Arithmetic {
            equation: "eddy viscosity",
            node: i,
        })
    } else {
        Ok(value)
    }
}

fn low_re_terms(
    grid: &RansChannelGrid,
    k: &[f64],
    u: &[f64],
    nu: f64,
    nut: &[f64],
    i: usize,
) -> (f64, f64) {
    if i == 0 {
        return (0.0, 0.0);
    }
    let a = grid.y_nodes[i] - grid.y_nodes[i - 1];
    let (gradient, curvature) = if i == grid.num_cells {
        (0.0, 2.0 * ((u[i - 1] - u[i]) / a) / a)
    } else {
        let b = grid.y_nodes[i + 1] - grid.y_nodes[i];
        let left = (k[i].sqrt() - k[i - 1].sqrt()) / a;
        let right = (k[i + 1].sqrt() - k[i].sqrt()) / b;
        (
            (b * left + a * right) / (a + b),
            2.0 * ((u[i + 1] - u[i]) / b - (u[i] - u[i - 1]) / a) / (a + b),
        )
    };
    (
        2.0 * nu * gradient * gradient,
        2.0 * nu * nut[i] * curvature * curvature,
    )
}

#[cfg(test)]
mod tests;
