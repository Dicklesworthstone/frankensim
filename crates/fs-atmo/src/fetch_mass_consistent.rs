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
//!
//! [`FetchAdjustedAnelastic`] additionally preserves mass continuity for an explicit
//! height-dependent density profile; the original mode assumes constant density.

use crate::{
    AnelasticFlowSample, FlatSiteLogLaw, FlowSample, KAPPA, MAX_Z0_M, MIN_Z0_M, Refusal, refuse,
};
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
                format!(
                    "projected z0 at max fetch ({z0_end:.4e}) outside [{MIN_Z0_M}, {MAX_Z0_M}]"
                ),
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

/// Authored density `rho(h) = rho_ref * exp(-inverse_scale_height * (h - h_ref))`.
/// The finite height interval is inclusive; zero inverse scale height means exactly
/// constant density. This kinematic profile supplies no gas equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExponentialDensityProfile {
    reference_density_kg_m3: f64,
    reference_height_m: f64,
    inverse_scale_height_per_m: f64,
    height_bounds_m: [f64; 2],
}

impl ExponentialDensityProfile {
    /// Admit a positive density, nonnegative inverse scale height and finite height
    /// interval above ground, containing the reference height.
    ///
    /// # Errors
    /// Refuses invalid inputs or a density that overflows or vanishes in the interval.
    pub fn new(
        reference_density_kg_m3: f64,
        reference_height_m: f64,
        inverse_scale_height_per_m: f64,
        height_bounds_m: [f64; 2],
    ) -> Result<Self, Refusal> {
        if !reference_density_kg_m3.is_finite()
            || reference_density_kg_m3 <= 0.0
            || !reference_height_m.is_finite()
            || !inverse_scale_height_per_m.is_finite()
            || inverse_scale_height_per_m < 0.0
            || !height_bounds_m.iter().all(|h| h.is_finite() && *h >= 0.0)
            || height_bounds_m[1] <= height_bounds_m[0]
            || reference_height_m < height_bounds_m[0]
            || reference_height_m > height_bounds_m[1]
        {
            return Err(refuse(
                "density-profile-invalid",
                "density, inverse scale height or reference/domain heights are invalid".into(),
                "supply positive density, nonnegative inverse scale height and ordered heights",
            ));
        }
        let profile = Self {
            reference_density_kg_m3,
            reference_height_m,
            inverse_scale_height_per_m,
            height_bounds_m,
        };
        // Monotonicity makes the endpoints sufficient to admit the entire interval.
        profile.density_kg_m3(height_bounds_m[0])?;
        profile.density_kg_m3(height_bounds_m[1])?;
        Ok(profile)
    }

    /// Physical density [kg/m³] at an admitted height above ground.
    ///
    /// # Errors
    /// Refuses nonfinite/out-of-domain heights and nonrepresentable density.
    pub fn density_kg_m3(&self, height_m: f64) -> Result<f64, Refusal> {
        if !height_m.is_finite()
            || height_m < self.height_bounds_m[0]
            || height_m > self.height_bounds_m[1]
        {
            return Err(refuse(
                "density-profile-domain",
                "height lies outside the declared density interval".into(),
                "query a finite height inside the source interval",
            ));
        }
        let density = self.reference_density_kg_m3
            * det::exp(-self.inverse_scale_height_per_m * (height_m - self.reference_height_m));
        if !density.is_finite() || density <= 0.0 {
            return Err(refuse(
                "density-profile-range",
                "density overflows or vanishes at the requested height".into(),
                "reduce the source height interval or density contrast",
            ));
        }
        Ok(density)
    }

    /// Compose the original horizontal mean with density-weighted turbulence.
    /// For q = rho_ref * curl(psi), u' = q/rho has
    /// du'_i/dh = (rho_ref/rho) * dv_i/dh + k*u'_i. Consequently
    /// div(rho*u') = 0, while the unchanged horizontal mean separately has
    /// zero mass divergence because rho and its speed depend only on height.
    /// Only the known solenoidal modal owners call this helper; applying it
    /// again to an anelastic field would be physically incorrect.
    pub(crate) fn compose_anelastic(
        &self,
        mean: &FlatSiteLogLaw,
        mut turbulent: FlowSample,
        height_m: f64,
    ) -> Result<AnelasticFlowSample, Refusal> {
        let density_kg_m3 = self.density_kg_m3(height_m)?;
        let k = self.inverse_scale_height_per_m;
        if k != 0.0 {
            let scale = self.reference_density_kg_m3 / density_kg_m3;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(refuse(
                    "anelastic-flow-range",
                    "turbulent density scale is not representable".into(),
                    "reduce the density contrast over the observation interval",
                ));
            }
            for (velocity, gradient) in turbulent.u.iter_mut().zip(&mut turbulent.grad) {
                *velocity *= scale;
                for derivative in &mut *gradient {
                    *derivative *= scale;
                }
                gradient[2] += k * *velocity;
            }
        }
        // Keep the original addition order so k=0 is bitwise legacy flow.
        turbulent.u[0] += mean.speed(height_m);
        turbulent.grad[0][2] += mean.dspeed_dh(height_m);
        let mass_flux_kg_m2_s = turbulent.u.map(|velocity| density_kg_m3 * velocity);
        if turbulent
            .u
            .iter()
            .chain(turbulent.grad.iter().flatten())
            .chain(mass_flux_kg_m2_s.iter())
            .any(|value| !value.is_finite())
        {
            return Err(refuse(
                "anelastic-flow-range",
                "velocity, gradient or mass flux is not representable".into(),
                "reduce the density contrast, wind speed or modal amplitude",
            ));
        }
        Ok(AnelasticFlowSample {
            density_kg_m3,
            flow: turbulent,
            mass_flux_kg_m2_s,
        })
    }
}

/// One consistent observation in the `[fetch, lateral, upward]` frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnelasticFetchSample {
    /// Local physical density [kg/m³].
    pub density_kg_m3: f64,
    /// Velocity [m/s], with zero lateral component.
    pub velocity_mps: [f64; 3],
    /// Advective mass flux `rho * velocity` [kg/(m² s)].
    pub mass_flux_kg_m2_s: [f64; 3],
}

/// Fetch log law with `div(rho * velocity) = 0` for the retained density profile.
/// Above `b = d + z0(x)`, the vertical correction is
/// `w = (u_star/kappa) * (z0'/z0) * integral_b^h(rho(s) ds) / rho(h)`.
/// The original roughness clamps and zero-velocity cutoff are retained. Fetch
/// endpoints and the cutoff are nonsmooth; continuity applies on each smooth piece.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FetchAdjustedAnelastic {
    wind: FetchAdjustedMassConsistent,
    density: ExponentialDensityProfile,
    speed_scale_mps: f64,
}

impl FetchAdjustedAnelastic {
    /// Freeze the wind and density sources, including all roughness cutoff heights.
    ///
    /// # Errors
    /// Refuses invalid wind parameters, nonrepresentable speed scale or a density
    /// domain that excludes any cutoff height along the clamped fetch profile.
    pub fn new(
        wind: FetchAdjustedMassConsistent,
        density: ExponentialDensityProfile,
    ) -> Result<Self, Refusal> {
        let profile = wind.fetch_profile;
        FetchRoughnessProfile::new(profile.z0_base_m, profile.dz0_dx, profile.max_fetch_m)?;
        FetchAdjustedMassConsistent::new(
            profile,
            wind.displacement_height_m,
            wind.reference_height_m,
            wind.reference_speed_mps,
        )?;
        for x_m in [0.0, profile.max_fetch_m] {
            density.density_kg_m3(wind.displacement_height_m + profile.z0_at(x_m))?;
        }
        let speed_scale_mps = wind.u_star() / KAPPA;
        if !speed_scale_mps.is_finite() || speed_scale_mps <= 0.0 {
            return Err(refuse(
                "anelastic-flow-range",
                "reference speed scale is not representable".into(),
                "use a resolved reference height and speed",
            ));
        }
        Ok(Self {
            wind,
            density,
            speed_scale_mps,
        })
    }

    /// Sample physical density, velocity and mass flux without allocation or quadrature.
    /// Outside the fetch interval roughness is clamped and vertical velocity is zero.
    /// At either fetch endpoint the interior derivative convention is retained.
    ///
    /// # Errors
    /// Refuses nonfinite fetch, heights outside the density source, and overflow
    /// in velocity, density-column correction or mass flux.
    pub fn sample(&self, x_m: f64, height_m: f64) -> Result<AnelasticFetchSample, Refusal> {
        if !x_m.is_finite() {
            return Err(refuse(
                "non-finite-input",
                "fetch must be finite".into(),
                "use finite fetch",
            ));
        }
        let density_kg_m3 = self.density.density_kg_m3(height_m)?;
        let profile = self.wind.fetch_profile;
        let z0 = profile.z0_at(x_m);
        let eff_h = height_m - self.wind.displacement_height_m;
        let mut velocity_mps = [0.0; 3];
        if eff_h > z0 {
            velocity_mps[0] = self.speed_scale_mps * det::ln(eff_h / z0);
            let gradient = if x_m >= 0.0 && x_m <= profile.max_fetch_m {
                profile.dz0_dx
            } else {
                0.0
            };
            velocity_mps[2] = self.speed_scale_mps * (1.0 / z0) * gradient * (eff_h - z0);
            if gradient != 0.0 && self.density.inverse_scale_height_per_m > 0.0 {
                let exponent = self.density.inverse_scale_height_per_m * (eff_h - z0);
                // expm1(x)/x resolves weak stratification, including x underflowing to zero.
                let correction = if exponent == 0.0 {
                    1.0
                } else {
                    det::expm1(exponent) / exponent
                };
                velocity_mps[2] *= correction;
            }
        }
        let mass_flux_kg_m2_s = velocity_mps.map(|speed| density_kg_m3 * speed);
        if velocity_mps
            .iter()
            .chain(mass_flux_kg_m2_s.iter())
            .any(|value| !value.is_finite())
        {
            return Err(refuse(
                "anelastic-flow-range",
                "velocity, density-column correction or mass flux overflows".into(),
                "reduce the height interval, density contrast or reference speed",
            ));
        }
        Ok(AnelasticFetchSample {
            density_kg_m3,
            velocity_mps,
            mass_flux_kg_m2_s,
        })
    }
}
