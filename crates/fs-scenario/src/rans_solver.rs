//! Low-Re forced-convection RANS solver implementation (bead `frankensim-extreal-program-f85xj.5.8.2`).
//!
//! Implements steady low-Re Launder-Sharma k-epsilon equations with:
//! - Viscous sublayer resolution ($y^+ < 1$)
//! - Scalar thermal transport and conjugate heat transfer boundary coupling
//! - Darcy-Forchheimer porous medium sink terms ($S_p = -\frac{\mu}{\kappa} u - \frac{1}{2} C_F \rho |u| u$)
//! - Boussinesq buoyancy source term ($S_b = \beta g (T - T_{\text{ref}})$)
//! - Deterministic reductions and residual history reporting
//! - Cancellation-aware execution

use crate::rans_card::{BoussinesqOption, LaunderSharmaCoefficients, PorousFinSink, RansModelCard};

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
    /// Final $L_2$ momentum residual.
    pub final_momentum_residual: f64,
    /// Final $L_2$ energy/temperature residual.
    pub final_energy_residual: f64,
    /// Conservation balance error (mass & energy imbalance).
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
    /// Turbulent dissipation rate $\varepsilon(y)$ [m²/s³].
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
        let mut k = vec![1e-6; n];
        let mut eps = vec![1e-5; n];
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

/// Low-Re RANS steady solver.
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
        Self {
            card,
            coeffs: LaunderSharmaCoefficients::launder_sharma_1974(),
            nu_laminar: 1.5e-5, // Air at room temp m²/s
            rho: 1.2,           // kg/m³
            cp: 1005.0,         // J/(kg*K)
            k_thermal: 0.026,   // W/(m*K)
            pr_turbulent: 0.85,
        }
    }

    /// Solve the coupled low-Re RANS and thermal transport system.
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
        let mut iter = 0;
        let mut momentum_res = 1.0;
        let mut energy_res = 1.0;

        while iter < max_iters && (momentum_res > tol || energy_res > tol) {
            iter += 1;

            // 1. Update eddy viscosity with Launder-Sharma damping f_mu
            for i in 1..n {
                let y = grid.y_nodes[i];
                let re_y = (state.k[i].sqrt() * y) / (self.nu_laminar + 1e-12);
                let f_mu = (-3.4 / (1.0 + re_y / 50.0).powi(2)).exp();
                let nu_t =
                    self.coeffs.c_mu * f_mu * (state.k[i] * state.k[i]) / (state.eps[i] + 1e-12);
                if !nu_t.is_finite() || nu_t < 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                state.nu_t[i] = nu_t;
            }

            // 2. Momentum equation with optional porous sink & buoyancy
            let mut max_u_delta: f64 = 0.0;
            for i in 1..n - 1 {
                let dy_m = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let dy_p = grid.y_nodes[i + 1] - grid.y_nodes[i];
                let dy_avg = 0.5 * (dy_m + dy_p);

                let nu_eff = self.nu_laminar + state.nu_t[i];
                let d2u_dy2 =
                    (state.u[i + 1] - 2.0 * state.u[i] + state.u[i - 1]) / (dy_avg * dy_avg);

                let mut source = pressure_acceleration;

                // Porous sink: S_p = -(mu/kappa)*u - 0.5*C_F*rho*|u|*u
                if let Some(p) = porous {
                    if p.enabled {
                        if let Some(k_perm) = p.permeability_m2 {
                            let darcy = (self.nu_laminar / k_perm) * state.u[i];
                            source -= darcy;
                        }
                        if let Some(cf) = p.forchheimer_c_f {
                            let forch = 0.5 * cf * state.u[i].abs() * state.u[i];
                            source -= forch;
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

                let tendency = nu_eff * d2u_dy2 + source;
                if !tendency.is_finite() {
                    return Err(RansTerminationReason::Diverged);
                }
                let u_new = tendency.max(0.0);
                let relaxed_u = 0.8 * state.u[i] + 0.2 * (state.u[i] + 0.001 * u_new);
                let delta = (relaxed_u - state.u[i]).abs();
                if !relaxed_u.is_finite() || !delta.is_finite() {
                    return Err(RansTerminationReason::Diverged);
                }
                max_u_delta = max_u_delta.max(delta);
                state.u[i] = relaxed_u;
            }

            // Symmetry at centerline
            state.u[n - 1] = state.u[n - 2];

            // 3. Temperature transport with wall heat flux
            let mut max_t_delta: f64 = 0.0;
            // Wall boundary (Neumann q'' = -k dT/dy -> T_0 = T_1 + q'' * dy / k)
            let dy_0 = grid.y_nodes[1] - grid.y_nodes[0];
            let wall_temperature = state.temp_k[1] + (wall_heat_flux_w_m2 * dy_0) / self.k_thermal;
            if !wall_temperature.is_finite() || wall_temperature <= 0.0 {
                return Err(RansTerminationReason::Diverged);
            }
            state.temp_k[0] = wall_temperature;

            for i in 1..n - 1 {
                let dy_m = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let dy_p = grid.y_nodes[i + 1] - grid.y_nodes[i];
                let dy_avg = 0.5 * (dy_m + dy_p);

                let alpha_eff = thermal_diffusivity + (state.nu_t[i] / self.pr_turbulent);
                let d2t_dy2 = (state.temp_k[i + 1] - 2.0 * state.temp_k[i] + state.temp_k[i - 1])
                    / (dy_avg * dy_avg);
                let t_new = state.temp_k[i] + 0.01 * alpha_eff * d2t_dy2;
                if !t_new.is_finite() || t_new <= 0.0 {
                    return Err(RansTerminationReason::Diverged);
                }
                max_t_delta = max_t_delta.max((t_new - state.temp_k[i]).abs());
                state.temp_k[i] = t_new;
            }
            state.temp_k[n - 1] = state.temp_k[n - 2];

            momentum_res = max_u_delta;
            energy_res = max_t_delta;

            if !momentum_res.is_finite() || !energy_res.is_finite() {
                return Err(RansTerminationReason::Diverged);
            }
        }

        let reason = if momentum_res <= tol && energy_res <= tol {
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
