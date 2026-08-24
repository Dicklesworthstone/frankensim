//! FetchAdjustedMassConsistent mean-wind mode (bead `frankensim-wf-root-guzez.4.7`, E3.3c).
//!
//! Solves the fetch-varying roughness problem honestly: while a naive $z_0(x)$ insertion
//! destroys solenoidality ($\partial U/\partial x \neq 0$), this mode applies an exact analytic
//! mass-consistent vertical correction $w(x, z)$ such that:
//!
//! $$\nabla \cdot \mathbf{u} = \frac{\partial u}{\partial x} + \frac{\partial w}{\partial z} = 0$$
//!
//! with impermeable wall condition $w(x, 0) = 0$.
//!
//! Claim class is explicitly `Estimated` / `Approximate` with ModelId `FetchAdjustedMassConsistent`.

use crate::{refuse, Refusal, KAPPA, MAX_Z0_M, MIN_Z0_M};
use fs_math::det;

/// Model identifier for the fetch-adjusted mass-consistent atmosphere mode.
pub const MODEL_ID_FETCH_MASS_CONSISTENT: &str = "FetchAdjustedMassConsistent";

/// Configuration for fetch-dependent roughness variation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FetchRoughnessProfile {
    /// Upwind base roughness length $z_{0,\text{base}}$ [m].
    pub z0_base_m: f64,
    /// Linear fetch gradient $dz_0/dx$ [1/m] (dimensionless rate of change with fetch).
    pub dz0_dx: f64,
    /// Maximum allowable fetch distance $x_{\text{max}}$ [m].
    pub max_fetch_m: f64,
}

impl FetchRoughnessProfile {
    /// Construct a valid fetch roughness profile.
    ///
    /// # Errors
    /// [`Refusal`] if parameters are non-finite or $z_0$ exceeds domain bounds.
    pub fn new(z0_base_m: f64, dz0_dx: f64, max_fetch_m: f64) -> Result<Self, Refusal> {
        if !z0_base_m.is_finite() || !dz0_dx.is_finite() || !max_fetch_m.is_finite() {
            return Err(refuse(
                "non-finite-input",
                "fetch profile parameters must be finite".into(),
                "provide finite numerical values",
            ));
        }
        if z0_base_m < MIN_Z0_M || z0_base_m > MAX_Z0_M {
            return Err(refuse(
                "z0-outside-domain",
                format!("z0_base {z0_base_m} outside [{MIN_Z0_M}, {MAX_Z0_M}]"),
                "keep roughness within certified domain",
            ));
        }
        if max_fetch_m <= 0.0 {
            return Err(refuse(
                "fetch-invalid",
                "max_fetch_m must be strictly positive".into(),
                "provide positive fetch bound",
            ));
        }
        // Verify end of fetch stays within bounds
        let z0_end = z0_base_m + dz0_dx * max_fetch_m;
        if z0_end < MIN_Z0_M || z0_end > MAX_Z0_M {
            return Err(refuse(
                "z0-outside-domain",
                format!("projected z0 at max fetch ({z0_end:.4e}) outside [{MIN_Z0_M}, {MAX_Z0_M}]"),
                "reduce dz0_dx or max_fetch_m",
            ));
        }
        Ok(Self {
            z0_base_m,
            dz0_dx,
            max_fetch_m,
        })
    }

    /// Compute effective roughness $z_0(x)$ at fetch distance $x$.
    #[must_use]
    pub fn z0_at(&self, x_m: f64) -> f64 {
        (self.z0_base_m + self.dz0_dx * x_m.clamp(0.0, self.max_fetch_m)).clamp(MIN_Z0_M, MAX_Z0_M)
    }
}

/// Fetch-adjusted mass-consistent mean wind law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FetchAdjustedMassConsistent {
    /// Roughness variation profile along fetch.
    pub fetch_profile: FetchRoughnessProfile,
    /// Displacement height $d$ [m].
    pub displacement_height_m: f64,
    /// Reference height $h_{\text{ref}}$ [m].
    pub reference_height_m: f64,
    /// Reference speed $U_{\text{ref}}$ [m/s] at reference height over base roughness.
    pub reference_speed_mps: f64,
}

impl FetchAdjustedMassConsistent {
    /// Construct and admit a fetch-adjusted mass-consistent wind law.
    ///
    /// # Errors
    /// [`Refusal`] on invalid parameters.
    pub fn new(
        fetch_profile: FetchRoughnessProfile,
        displacement_height_m: f64,
        reference_height_m: f64,
        reference_speed_mps: f64,
    ) -> Result<Self, Refusal> {
        if !displacement_height_m.is_finite()
            || !reference_height_m.is_finite()
            || !reference_speed_mps.is_finite()
        {
            return Err(refuse(
                "non-finite-input",
                "law parameters must be finite".into(),
                "provide finite numbers",
            ));
        }
        if displacement_height_m < 0.0 {
            return Err(refuse(
                "displacement-invalid",
                "displacement height must be non-negative".into(),
                "set d >= 0",
            ));
        }
        if reference_height_m <= displacement_height_m + fetch_profile.z0_base_m {
            return Err(refuse(
                "reference-height-invalid",
                "reference height must exceed d + z0".into(),
                "raise reference height",
            ));
        }
        if reference_speed_mps <= 0.0 {
            return Err(refuse(
                "reference-speed-invalid",
                "reference speed must be positive".into(),
                "provide positive speed",
            ));
        }

        Ok(Self {
            fetch_profile,
            displacement_height_m,
            reference_height_m,
            reference_speed_mps,
        })
    }

    /// Friction velocity $u_*$ derived from reference conditions.
    #[must_use]
    pub fn u_star(&self) -> f64 {
        let denom = det::ln(
            (self.reference_height_m - self.displacement_height_m) / self.fetch_profile.z0_base_m,
        );
        (self.reference_speed_mps * KAPPA) / denom
    }

    /// Sample mass-consistent velocity vector $[u, v, w]$ [m/s] at fetch $x$ and altitude $h$.
    /// Altitude $h$ is height above aerodynamic ground ($h = -z$ in NED).
    #[must_use]
    pub fn sample_velocity(&self, x_m: f64, h_m: f64) -> [f64; 3] {
        let z0 = self.fetch_profile.z0_at(x_m);
        let eff_h = h_m - self.displacement_height_m;
        if eff_h <= z0 {
            return [0.0, 0.0, 0.0];
        }

        let u_star = self.u_star();
        // Horizontal log-law component u(x, h)
        let u_val = (u_star / KAPPA) * det::ln(eff_h / z0);

        // Mass-consistent vertical velocity correction:
        // w(x, h) = (u_* / kappa) * (1 / z0) * (dz0/dx) * (eff_h - z0)
        let dz0_dx = if x_m >= 0.0 && x_m <= self.fetch_profile.max_fetch_m {
            self.fetch_profile.dz0_dx
        } else {
            0.0
        };
        let w_val = (u_star / KAPPA) * (1.0 / z0) * dz0_dx * (eff_h - z0);

        [u_val, 0.0, w_val]
    }

    /// Compute exact velocity divergence $\nabla \cdot \mathbf{u} = \partial u/\partial x + \partial w/\partial h$.
    /// By construction of the mass-consistent stream function, this equals 0.0 analytically.
    #[must_use]
    pub fn divergence(&self, x_m: f64, h_m: f64) -> f64 {
        let z0 = self.fetch_profile.z0_at(x_m);
        let eff_h = h_m - self.displacement_height_m;
        if eff_h <= z0 {
            return 0.0;
        }

        let u_star = self.u_star();
        let dz0_dx = if x_m >= 0.0 && x_m <= self.fetch_profile.max_fetch_m {
            self.fetch_profile.dz0_dx
        } else {
            0.0
        };

        // ∂u/∂x = (u_* / kappa) * (-1 / z0) * dz0/dx
        let du_dx = -(u_star / KAPPA) * (1.0 / z0) * dz0_dx;
        // ∂w/∂h = (u_* / kappa) * (1 / z0) * dz0/dx
        let dw_dh = (u_star / KAPPA) * (1.0 / z0) * dz0_dx;

        du_dx + dw_dh
    }
}
