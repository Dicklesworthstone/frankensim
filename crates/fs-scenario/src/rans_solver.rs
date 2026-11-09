//! Low-Re forced-convection RANS solver implementation (bead `frankensim-extreal-program-f85xj.5.8.2`).
//!
//! Frozen-turbulence one-dimensional mean-flow/scalar screen with:
//! `k` and modified dissipation are supplied fields, not transported unknowns.
//! `advance_launder_sharma` separately transports momentum, k and modified
//! dissipation over a physical interval; `advance_launder_sharma_with_heat`
//! also advances conservative passive heat using the evolving closure.
//! `advance_launder_sharma_adaptive` controls local temporal error over a
//! bounded constant-forcing physical window, with optional passive heat.
//! These bounded entry points do not imply empirically validated general RANS;
//! see the method and crate contracts.
//! - Explicit stretched-grid wall-normal transport; wall resolution is caller supplied
//! - Scalar thermal transport with prescribed wall heat flux
//! - Darcy-Forchheimer porous medium sink terms ($S_p = -\frac{\mu}{\kappa} u - \frac{1}{2} C_F \rho |u| u$)
//! - Boussinesq buoyancy source term ($S_b = \beta g (T - T_{\text{ref}})$)
//! - Fixed-order sweeps and explicit iteration diagnostics
//!
//! Nonzero heating with an insulated centerline has no thermal steady state
//! without axial heat removal. This screen does not add that missing transport.

use crate::rans_card::{BoussinesqOption, LaunderSharmaCoefficients, PorousFinSink, RansModelCard};

mod transport;
pub use transport::{
    RansAdaptiveControls, RansAdaptiveError, RansAdaptiveHeatBalance, RansAdaptiveReport,
    RansAdaptiveWorkspace, RansFieldTolerances, RansHeatLoad, RansThermalTransportReport,
    RansTransportError, RansTransportReport, RansTransportWorkspace,
};

/// Convergence termination reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RansTerminationReason {
    /// Residuals dropped below requested tolerance.
    Converged,
    /// Hit max iteration ceiling before convergence.
    MaxIterationsExceeded,
    /// Cancelled via cancellation scope.
    Cancelled,
    /// Refused due to numerical divergence or invalid initial state.
    Diverged,
}

/// Convergence and conservation history for a RANS solution.
#[derive(Debug, Clone, PartialEq)]
pub struct RansConvergenceReport {
    /// Iteration count performed.
    pub iterations: usize,
    /// Largest accepted velocity update in the last pseudo-time sweep [m/s].
    /// Convergence is tested on the unscaled momentum equation rate.
    pub final_momentum_residual: f64,
    /// Largest accepted temperature update in the last pseudo-time sweep [K].
    /// Convergence is tested on the unscaled temperature equation rate.
    pub final_energy_residual: f64,
    /// Legacy aggregate iteration-change diagnostic, not a conservation audit.
    pub conservation_imbalance: f64,
    /// Termination disposition.
    pub reason: RansTerminationReason,
}

/// 1D Channel discretization grid.
#[derive(Debug, Clone, PartialEq)]
pub struct RansChannelGrid {
    /// Number of wall-normal cells across channel half-height $H$.
    pub num_cells: usize,
    /// Channel half-height $H$ in meters.
    pub half_height_m: f64,
    /// Node coordinates $y_i$ from wall ($y=0$) to centerline ($y=H$).
    pub y_nodes: Vec<f64>,
}

impl RansChannelGrid {
    /// Construct a geometric stretched grid resolving the viscous sublayer.
    ///
    /// # Panics
    /// Panics if `num_cells < 4` or `half_height_m <= 0.0`.
    #[must_use]
    pub fn new_stretched(num_cells: usize, half_height_m: f64, stretch_ratio: f64) -> Self {
        assert!(num_cells >= 4, "num_cells must be >= 4");
        assert!(half_height_m > 0.0, "half_height must be positive");
        let mut y_nodes = Vec::with_capacity(num_cells + 1);
        y_nodes.push(0.0);
        let mut dy = 1.0;
        let mut sum = 0.0;
        let mut dys = Vec::with_capacity(num_cells);
        for _ in 0..num_cells {
            dys.push(dy);
            sum += dy;
            dy *= stretch_ratio;
        }
        let scale = half_height_m / sum;
        let mut current_y = 0.0;
        for d in dys {
            current_y += d * scale;
            y_nodes.push(current_y);
        }
        Self {
            num_cells,
            half_height_m,
            y_nodes,
        }
    }
}

/// Solution state of the low-Re RANS equations.
#[derive(Debug, Clone, PartialEq)]
pub struct RansFieldState {
    /// Streamwise velocity $u(y)$ [m/s].
    pub u: Vec<f64>,
    /// Turbulent kinetic energy $k(y)$ [m²/s²].
    pub k: Vec<f64>,
    /// Frozen modified LS dissipation $\widetilde\varepsilon(y)$ [m²/s³].
    pub eps: Vec<f64>,
    /// Eddy viscosity $\nu_t(y)$ [m²/s].
    pub nu_t: Vec<f64>,
    /// Temperature $T(y)$ [K].
    pub temp_k: Vec<f64>,
}

impl RansFieldState {
    /// Initialize state with laminar profile and ambient temperature.
    #[must_use]
    pub fn new_initial(grid: &RansChannelGrid, u_bulk: f64, t_ambient_k: f64) -> Self {
        let n = grid.num_cells + 1;
        let mut u = vec![0.0; n];
        let mut k = vec![0.0; n];
        let mut eps = vec![0.0; n];
        let nu_t = vec![1e-8; n];
        let temp_k = vec![t_ambient_k; n];

        for i in 1..n {
            let eta = grid.y_nodes[i] / grid.half_height_m;
            u[i] = 1.5 * u_bulk * (2.0 * eta - eta * eta); // Poiseuille shape
            k[i] = 0.01 * u_bulk * u_bulk;
            eps[i] = 0.09 * (k[i] * k[i]) / (1e-5 + 1e-4 * grid.half_height_m);
        }

        Self {
            u,
            k,
            eps,
            nu_t,
            temp_k,
        }
    }
}

/// Frozen-field mean-flow/scalar relaxation screen.
pub struct RansSolver {
    /// Associated admitted model card.
    pub card: RansModelCard,
    /// Closure coefficients.
    pub coeffs: LaunderSharmaCoefficients,
    /// Molecular kinematic viscosity $\nu$ [m²/s].
    pub nu_laminar: f64,
    /// Fluid density $\rho$ [kg/m³].
    pub rho: f64,
    /// Fluid specific heat capacity $c_p$ [J/(kg*K)].
    pub cp: f64,
    /// Fluid thermal conductivity $k$ [W/(m*K)].
    pub k_thermal: f64,
    /// Turbulent Prandtl number $\text{Pr}_t$.
    pub pr_turbulent: f64,
}

impl RansSolver {
    /// Create a new solver bound to an admitted RansModelCard.
    #[must_use]
    pub fn new(card: RansModelCard) -> Self {
        let coeffs = *card.coefficients();
        let pr_turbulent = card.turbulent_prandtl();
        Self {
            card,
            coeffs,
            nu_laminar: 1.5e-5, // Air at room temp m²/s
            rho: 1.2,           // kg/m³
            cp: 1005.0,         // J/(kg*K)
            k_thermal: 0.026,   // W/(m*K)
            pr_turbulent,
        }
    }

    /// Relax the frozen-turbulence mean-flow/scalar screen.
    /// `tol` is an absolute equation-rate threshold, applied separately to
    /// momentum [m/s²] and temperature [K/s]; small pseudo-steps cannot satisfy it.
    /// The in-place sweep and existing boundary reconstruction are retained;
    /// spatial face-flux consistency is not a time-integrated conservation claim.
    ///
    /// # Errors
    /// Returns [`RansTerminationReason::Diverged`] on invalid material/state
    /// inputs or a nonfinite numerical update. Input admission precedes state
    /// mutation; a failed iteration never publishes a nonfinite field value.
    pub fn solve(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        dp_dx: f64,
        wall_heat_flux_w_m2: f64,
        porous: Option<PorousFinSink>,
        buoyancy: Option<BoussinesqOption>,
        max_iters: usize,
        tol: f64,
    ) -> Result<RansConvergenceReport, RansTerminationReason> {
        let n = grid
            .num_cells
            .checked_add(1)
            .filter(|&n| n >= 3)
            .ok_or(RansTerminationReason::Diverged)?;
        // Density is the SAME physical mass storage in pressure acceleration
        // and rho*cp thermal storage. Refuse invalid coefficients before an
        // IEEE max reduction can silently hide their NaN consequences.
        if [
            self.rho,
            self.nu_laminar,
            self.cp,
            self.k_thermal,
            self.pr_turbulent,
            tol,
        ]
        .into_iter()
        .any(|value| !value.is_finite() || value <= 0.0)
            || !self.coeffs.c_mu.is_finite()
            || self.coeffs.c_mu < 0.0
            || !dp_dx.is_finite()
            || !wall_heat_flux_w_m2.is_finite()
            || grid.y_nodes.len() != n
            || grid.y_nodes.iter().any(|value| !value.is_finite())
            || grid.y_nodes.windows(2).any(|pair| pair[0] >= pair[1])
            || [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k]
                .into_iter()
                .any(|values| values.len() != n || values.iter().any(|value| !value.is_finite()))
            || [&state.k, &state.eps, &state.nu_t]
                .into_iter()
                .any(|values| values.iter().any(|value| *value < 0.0))
            || state.temp_k.iter().any(|value| *value <= 0.0)
        {
            return Err(RansTerminationReason::Diverged);
        }
        let pressure_acceleration = -dp_dx / self.rho;
        let thermal_diffusivity = (self.k_thermal / self.rho) / self.cp;
        if !pressure_acceleration.is_finite()
            || !thermal_diffusivity.is_finite()
            || thermal_diffusivity <= 0.0
        {
            return Err(RansTerminationReason::Diverged);
        }
        if let Some(p) = porous.filter(|p| p.enabled) {
            if p.permeability_m2
                .is_some_and(|value| !value.is_finite() || value <= 0.0)
                || p.forchheimer_c_f
                    .is_some_and(|value| !value.is_finite() || value < 0.0)
            {
                return Err(RansTerminationReason::Diverged);
            }
        }
        let mut iter = 0;
        let mut momentum_res = 1.0;
        let mut energy_res = 1.0;
        let mut momentum_rate = f64::INFINITY;
        let mut energy_rate = f64::INFINITY;

        while iter < max_iters && (momentum_rate > tol || energy_rate > tol) {
            iter += 1;

            // 1. Update eddy viscosity with Launder-Sharma damping f_mu.
            // A resolved no-slip wall cannot retain the initializer's eddy seed.
            state.nu_t[0] = 0.0;
            for i in 1..n {
                // LS uses the turbulent Reynolds number k²/(nu * eps),
                // not sqrt(k)*y/nu. NASA TM-110446, Part B §3.
                // Exact laminar/disabled controls need no dimensional floors.
                let k2_over_eps = if state.k[i] == 0.0 || self.coeffs.c_mu == 0.0 {
                    0.0
                } else {
                    if state.eps[i] == 0.0 {
                        return Err(RansTerminationReason::Diverged);
                    }
                    state.k[i] * (state.k[i] / state.eps[i])
                };
                let re_t = k2_over_eps / self.nu_laminar;
                // An overflowing Re_t has the finite high-Re limit f_mu = 1.
                let inverse = 1.0 / (1.0 + re_t / 50.0);
                let f_mu = (-3.4 * inverse * inverse).exp();
                let nu_t = self.coeffs.c_mu * f_mu * k2_over_eps;
                if !k2_over_eps.is_finite() || !nu_t.is_finite() || nu_t < 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                state.nu_t[i] = nu_t;
            }

            // 2. Momentum equation with optional porous sink & buoyancy.
            state.u[0] = 0.0;
            state.u[n - 1] = state.u[n - 2];
            let mut max_u_delta: f64 = 0.0;
            momentum_rate = 0.0;
            for i in 1..n - 1 {
                let dy_m = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let dy_p = grid.y_nodes[i + 1] - grid.y_nodes[i];
                let dy_avg = 0.5 * (dy_m + dy_p);

                // Divergence of the shared face stress, including grad(nu_t).
                // The two physical spacings must remain distinct on stretched grids.
                let lower = ((self.nu_laminar + 0.5 * state.nu_t[i - 1] + 0.5 * state.nu_t[i])
                    / dy_m)
                    / dy_avg;
                let upper = ((self.nu_laminar + 0.5 * state.nu_t[i] + 0.5 * state.nu_t[i + 1])
                    / dy_p)
                    / dy_avg;
                let diffusion =
                    lower * (state.u[i - 1] - state.u[i]) + upper * (state.u[i + 1] - state.u[i]);
                let mut diagonal = lower + upper;
                let mut source = pressure_acceleration;

                // Porous sink: S_p = -(mu/kappa)*u - 0.5*C_F*rho*|u|*u
                if let Some(p) = porous {
                    if p.enabled {
                        if let Some(k_perm) = p.permeability_m2 {
                            let rate = self.nu_laminar / k_perm;
                            source -= rate * state.u[i];
                            diagonal += rate;
                        }
                        if let Some(cf) = p.forchheimer_c_f {
                            let rate = cf * state.u[i].abs();
                            source -= 0.5 * rate * state.u[i];
                            diagonal += rate;
                        }
                    }
                }

                // Buoyancy source: S_b = beta * g * (T - T_ref)
                if let Some(b) = buoyancy {
                    if b.enabled {
                        if let Some(beta) = b.beta_per_k {
                            let g = 9.81;
                            source += beta * g * (state.temp_k[i] - b.reference_temperature_k);
                        }
                    }
                }

                let tendency = diffusion + source;
                if !tendency.is_finite() || !diagonal.is_finite() || diagonal <= 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                // Preserve braking and flow reversal. Bound the explicit relaxation
                // by the diffusion/drag diagonal instead of clipping acceleration.
                let pseudo_dt = 0.0002_f64.min(0.9 / diagonal);
                let relaxed_u = state.u[i] + pseudo_dt * tendency;
                let delta = (relaxed_u - state.u[i]).abs();
                if !relaxed_u.is_finite() || !delta.is_finite() {
                    return Err(RansTerminationReason::Diverged);
                }
                max_u_delta = max_u_delta.max(delta);
                momentum_rate = momentum_rate.max(tendency.abs());
                state.u[i] = relaxed_u;
            }

            // Symmetry at centerline
            state.u[n - 1] = state.u[n - 2];

            // 3. Temperature transport with wall heat flux
            let mut max_t_delta: f64 = 0.0;
            energy_rate = 0.0;
            state.temp_k[n - 1] = state.temp_k[n - 2];
            // Use the SAME face diffusivity as the scalar operator. A molecular
            // ghost gradient combined with turbulent face transport would multiply
            // the imposed heat input by alpha_eff/alpha_molecular.
            let dy_0 = grid.y_nodes[1] - grid.y_nodes[0];
            let wall_diffusivity = thermal_diffusivity
                + (0.5 * state.nu_t[0] + 0.5 * state.nu_t[1]) / self.pr_turbulent;
            let wall_flux = (wall_heat_flux_w_m2 / self.rho) / self.cp;
            if !wall_diffusivity.is_finite() || wall_diffusivity <= 0.0 || !wall_flux.is_finite() {
                return Err(RansTerminationReason::Diverged);
            }
            let wall_temperature = state.temp_k[1] + wall_flux * dy_0 / wall_diffusivity;
            if !wall_temperature.is_finite() || wall_temperature <= 0.0 {
                return Err(RansTerminationReason::Diverged);
            }
            state.temp_k[0] = wall_temperature;

            for i in 1..n - 1 {
                let dy_m = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let dy_p = grid.y_nodes[i + 1] - grid.y_nodes[i];
                let dy_avg = 0.5 * (dy_m + dy_p);

                let lower = ((thermal_diffusivity
                    + (0.5 * state.nu_t[i - 1] + 0.5 * state.nu_t[i]) / self.pr_turbulent)
                    / dy_m)
                    / dy_avg;
                let upper = ((thermal_diffusivity
                    + (0.5 * state.nu_t[i] + 0.5 * state.nu_t[i + 1]) / self.pr_turbulent)
                    / dy_p)
                    / dy_avg;
                let diagonal = lower + upper;
                let tendency = lower * (state.temp_k[i - 1] - state.temp_k[i])
                    + upper * (state.temp_k[i + 1] - state.temp_k[i]);
                if !tendency.is_finite() || !diagonal.is_finite() || diagonal <= 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                let pseudo_dt = 0.01_f64.min(0.9 / diagonal);
                let t_new = state.temp_k[i] + pseudo_dt * tendency;
                if !t_new.is_finite() || t_new <= 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                max_t_delta = max_t_delta.max((t_new - state.temp_k[i]).abs());
                energy_rate = energy_rate.max(tendency.abs());
                state.temp_k[i] = t_new;
            }
            state.temp_k[n - 1] = state.temp_k[n - 2];

            momentum_res = max_u_delta;
            energy_res = max_t_delta;

            if !momentum_res.is_finite() || !energy_res.is_finite() {
                return Err(RansTerminationReason::Diverged);
            }
        }

        let reason = if momentum_rate <= tol && energy_rate <= tol {
            RansTerminationReason::Converged
        } else {
            RansTerminationReason::MaxIterationsExceeded
        };

        Ok(RansConvergenceReport {
            iterations: iter,
            final_momentum_residual: momentum_res,
            final_energy_residual: energy_res,
            conservation_imbalance: (momentum_res + energy_res) * 0.5,
            reason,
        })
    }
}

#[cfg(test)]
mod turbulence_closure_regressions {
    use super::*;
    use crate::rans_card::RansCardDraft;

    fn sample_card() -> RansModelCard {
        RansCardDraft::launder_sharma_channel("electronics-cooling/e10-rans")
            .freeze()
            .unwrap()
    }

    #[test]
    fn solver_uses_the_admitted_turbulence_coefficients_and_prandtl_number() {
        let mut draft = RansCardDraft::launder_sharma_channel("electronics-cooling/custom-rans");
        draft.coefficients.c_mu = 0.095;
        draft.coefficients.c_eps_1 = 1.5;
        draft.coefficients.c_eps_2 = 1.9;
        draft.coefficients.sigma_k = 1.2;
        draft.coefficients.sigma_eps = 1.4;
        draft.turbulent_prandtl = 0.7;
        let solver = RansSolver::new(draft.freeze().unwrap());
        assert_eq!(solver.coeffs.c_mu, 0.095);
        assert_eq!(solver.coeffs.c_eps_1, 1.5);
        assert_eq!(solver.coeffs.c_eps_2, 1.9);
        assert_eq!(solver.coeffs.sigma_k, 1.2);
        assert_eq!(solver.coeffs.sigma_eps, 1.4);
        assert_eq!(solver.pr_turbulent, 0.7);
    }

    #[test]
    fn ls_viscosity_uses_turbulent_reynolds_number_independently_of_wall_distance() {
        let grid = RansChannelGrid::new_stretched(8, 1.0, 1.12);
        let solver = RansSolver::new(sample_card());
        let mut state = RansFieldState::new_initial(&grid, 0.0, 300.0);
        for i in 1..=grid.num_cells {
            state.k[i] = 3e-6;
            state.eps[i] = 2e-8;
        }
        solver
            .solve(&grid, &mut state, 0.0, 0.0, None, None, 1, 1e-12)
            .unwrap();
        // Independent literal equation: Re_t=30, not a wall-distance Reynolds number.
        let expected = 0.09 * (-3.4_f64 / 1.6_f64.powi(2)).exp() * 4.5e-4;
        for &actual in &state.nu_t[1..] {
            assert!((actual / expected - 1.0).abs() < 1e-13);
        }
    }

    #[test]
    fn viscosity_is_invariant_under_consistent_length_and_time_unit_rescaling() {
        let base_grid = RansChannelGrid::new_stretched(8, 1.0, 1.0);
        let solver = RansSolver::new(sample_card());
        let mut base = RansFieldState::new_initial(&base_grid, 0.0, 300.0);
        base.k.fill(3e-6);
        base.eps.fill(2e-8);
        solver
            .solve(&base_grid, &mut base, 0.0, 0.0, None, None, 1, 1e-12)
            .unwrap();
        for (length, time) in [(1e-3_f64, 1e3_f64), (1e3, 1e-3)] {
            let grid = RansChannelGrid::new_stretched(8, length, 1.0);
            let mut scaled_solver = RansSolver::new(sample_card());
            scaled_solver.nu_laminar *= length * length / time;
            let mut scaled = RansFieldState::new_initial(&grid, 0.0, 300.0);
            scaled.k.fill(3e-6 * (length / time).powi(2));
            scaled.eps.fill(2e-8 * length.powi(2) / time.powi(3));
            scaled_solver
                .solve(&grid, &mut scaled, 0.0, 0.0, None, None, 1, 1e-12)
                .unwrap();
            for (&original, &actual) in base.nu_t[1..].iter().zip(&scaled.nu_t[1..]) {
                assert!((actual / (original * length * length / time) - 1.0).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn zero_turbulence_and_disabled_closure_are_exact_without_epsilon_floors() {
        let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
        for disabled in [false, true] {
            let mut solver = RansSolver::new(sample_card());
            let mut state = RansFieldState::new_initial(&grid, 0.0, 300.0);
            state.k.fill(if disabled { 1.0 } else { 0.0 });
            state.eps.fill(0.0);
            if disabled {
                solver.coeffs.c_mu = 0.0;
            }
            solver
                .solve(&grid, &mut state, 0.0, 0.0, None, None, 1, 1e-12)
                .unwrap();
            assert!(state.nu_t[1..].iter().all(|&value| value == 0.0));
        }
        let solver = RansSolver::new(sample_card());
        let mut state = RansFieldState::new_initial(&grid, 0.0, 300.0);
        state.k[1] = 1e-6;
        state.eps[1] = 0.0;
        assert_eq!(
            solver.solve(&grid, &mut state, 0.0, 0.0, None, None, 1, 1e-12),
            Err(RansTerminationReason::Diverged),
        );
    }
}

#[cfg(test)]
mod mean_flow_relaxation_regressions {
    use super::*;
    use crate::rans_card::RansCardDraft;

    fn solver() -> RansSolver {
        RansSolver::new(
            RansCardDraft::launder_sharma_channel("electronics-cooling/e10-rans")
                .freeze()
                .unwrap(),
        )
    }

    fn laminar(grid: &RansChannelGrid, bulk: f64) -> RansFieldState {
        let mut state = RansFieldState::new_initial(grid, bulk, 300.0);
        state.k.fill(0.0);
        state.eps.fill(0.0);
        state.nu_t.fill(0.0);
        state
    }

    #[test]
    fn pressure_reversal_preserves_signed_acceleration() {
        let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
        let solver = solver();
        let mut forward = laminar(&grid, 0.0);
        let mut reverse = forward.clone();
        solver
            .solve(&grid, &mut forward, -2.0, 0.0, None, None, 1, 1e-12)
            .unwrap();
        solver
            .solve(&grid, &mut reverse, 2.0, 0.0, None, None, 1, 1e-12)
            .unwrap();
        assert!((reverse.u[1] + 0.0002 * 2.0 / solver.rho).abs() < 1e-15);
        for (&a, &b) in forward.u.iter().zip(&reverse.u) {
            assert_eq!(a, -b);
        }
    }

    #[test]
    fn stretched_diffusion_preserves_analytic_interior_profiles() {
        let grid = RansChannelGrid::new_stretched(8, 1.0, 1.6);
        let solver = solver();
        for quadratic in [false, true] {
            let mut state = laminar(&grid, 0.0);
            for (u, &y) in state.u.iter_mut().zip(&grid.y_nodes) {
                *u = if quadratic {
                    1.5 * (2.0 * y - y * y)
                } else {
                    y
                };
            }
            let before = state.u.clone();
            let dp = if quadratic {
                -3.0 * solver.nu_laminar * solver.rho
            } else {
                0.0
            };
            solver
                .solve(&grid, &mut state, dp, 0.0, None, None, 1, 1e-12)
                .unwrap();
            // The final interior row touches the retained first-order symmetry
            // reconstruction. All other rows have the exact analytic stencil.
            for i in 1..grid.num_cells - 1 {
                assert!((state.u[i] - before[i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn stiff_passive_drag_brakes_and_small_steps_cannot_fake_convergence() {
        let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
        let solver = solver();
        for speed in [1.0, 1e-6] {
            let mut state = laminar(&grid, 0.0);
            state.u[1..].fill(speed);
            let porous = PorousFinSink {
                enabled: true,
                permeability_m2: Some(1.5e-14),
                forchheimer_c_f: Some(if speed == 1.0 { 2e8 } else { 0.0 }),
            };
            let report = solver
                .solve(&grid, &mut state, 0.0, 0.0, Some(porous), None, 1, 1e-4)
                .unwrap();
            assert!(state.u[1] > 0.0 && state.u[1] < speed);
            assert_eq!(report.reason, RansTerminationReason::MaxIterationsExceeded);
            if speed < 1e-4 {
                assert!(report.final_momentum_residual < 1e-4);
            }
        }
    }

    #[test]
    fn stretched_thermal_diffusion_obeys_the_discrete_maximum_principle() {
        let grid = RansChannelGrid::new_stretched(32, 0.01, 1.15);
        let solver = solver();
        let mut state = laminar(&grid, 0.0);
        state.temp_k[16] = 1000.0;
        let report = solver
            .solve(&grid, &mut state, 0.0, 0.0, None, None, 100, 1e-12)
            .unwrap();
        assert!(report.final_energy_residual.is_finite());
        assert!(state.temp_k.iter().all(|&t| (300.0..=1000.0).contains(&t)));
        assert!(state.temp_k[16] < 1000.0);
        assert_eq!(state.nu_t[0], 0.0);
    }

    #[test]
    fn heated_closed_screen_stays_finite_without_claiming_steady_heat_removal() {
        let grid = RansChannelGrid::new_stretched(16, 0.01, 1.1);
        let solver = solver();
        let mut state = RansFieldState::new_initial(&grid, 1.0, 300.0);
        let report = solver
            .solve(&grid, &mut state, -10.0, 1000.0, None, None, 50, 1e-4)
            .unwrap();
        assert_eq!(report.reason, RansTerminationReason::MaxIterationsExceeded);
        assert!((300.0..500.0).contains(&state.temp_k[0]));
        assert!(state.temp_k.iter().all(|&t| t.is_finite() && t >= 300.0));
    }

    #[test]
    fn imposed_wall_heat_is_not_multiplied_by_turbulent_diffusivity() {
        let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
        let solver = solver();
        let mut state = laminar(&grid, 0.0);
        state.k.fill(1e-3);
        state.eps.fill(1e-3);
        solver
            .solve(&grid, &mut state, 0.0, 40.0, None, None, 1, 1e-12)
            .unwrap();
        assert!(state.nu_t[1] > 0.0);
        // The first interior node initially receives only the imposed face flux:
        // no other temperature gradient or axial transport has been supplied.
        let expected = 300.0 + 0.01 * 40.0 / (solver.rho * solver.cp * 0.25);
        assert!((state.temp_k[1] - expected).abs() < 1e-12);
        let molecular_wall = 300.0 + 40.0 * 0.25 / solver.k_thermal;
        assert!(state.temp_k[0] < molecular_wall);
    }

    #[test]
    fn viscosity_overflow_and_the_finite_high_reynolds_limit_are_distinct() {
        let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
        let mut solver = solver();
        let mut state = laminar(&grid, 0.0);
        state.k[1] = 1e200;
        state.eps[1] = 1.0;
        assert_eq!(
            solver.solve(&grid, &mut state, 0.0, 0.0, None, None, 1, 1e-12),
            Err(RansTerminationReason::Diverged),
        );
        let mut state = laminar(&grid, 0.0);
        solver.nu_laminar = 1e-320;
        state.k.fill(1.0);
        state.eps.fill(1.0);
        solver
            .solve(&grid, &mut state, 0.0, 0.0, None, None, 1, 1e-12)
            .unwrap();
        assert!(state.nu_t[1..].iter().all(|&value| value == 0.09));
        assert_eq!(state.nu_t[0], 0.0);
    }
}
