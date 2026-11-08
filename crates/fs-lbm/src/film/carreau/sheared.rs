//! Original-law Couette--Nusselt films with prescribed tangential traction.
//!
//! Steady momentum gives `tau(y) = tau_s + rho*g_t*(h-y)`. The original
//! material slit owner inverts this ONE stress profile and integrates its
//! positive mobility moments. It neither freezes an apparent viscosity nor
//! adds independently solved nonlinear pressure and shear profiles.
//!
//! Mechanical power follows by integrating `tau*du/dy`:
//! `rho*g_t*q + tau_s*u(h) = integral_0^h tau*du/dy dy` for a stationary
//! no-slip wall. Either input port may remove power during counterflow;
//! the original local and integrated viscous dissipation remain nonnegative.
//! The constant-source limit is the native [`crate::film::ShearedFilm`].
//!
//! Thickness, temperature, density, acceleration and surface traction are
//! held. Traction may represent gas drag or a surface-tension gradient, but
//! neither is inferred here. No gas dynamics, stability, contact line, slip,
//! thermal feedback or generalized turbulent threshold is introduced.

use super::{CarreauFilm, CarreauFilmError, CarreauFilmSample, material_call, poll, positive};
use crate::film::{FilmError, FilmPower};
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;
use fs_material::fluid_transport::rheology::slit::{
    CarreauSlitControls, CarreauSlitProfile, LaminarCarreauSlit,
};

/// Actual SI geometry and forcing of a flat film on a stationary no-slip wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauShearedFilmSpec {
    pub density_kg_m3: f64,
    pub thickness_m: f64,
    /// Signed tangential body acceleration [m/s²].
    pub tangential_acceleration_m_s2: f64,
    /// Signed traction applied to the liquid at the outward free surface [Pa].
    pub surface_shear_pa: f64,
}

/// Immutable original-law profile, allowing reversal and zero-net-flux flow.
/// Construction and sampling use the original bounded, allocation-free slit
/// inverse/quadrature and retain its source-rate and temperature support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauShearedFilm {
    spec: CarreauShearedFilmSpec,
    profile: CarreauSlitProfile,
    power: FilmPower,
    work_residual_w_m2: f64,
    surface_mobility_m_s_pa: f64,
    flow_surface_mobility_m2_s_pa: f64,
    flow_gradient_mobility_m3_s_pa: f64,
}

impl CarreauShearedFilm {
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauShearedFilmSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauSlitControls,
    ) -> Result<Self, CarreauFilmError> {
        Self::with_cancel(
            law,
            temperature_k,
            spec,
            maximum_laminar_shear_rate_per_s,
            controls,
            || false,
        )
    }

    /// Resolve the imposed traction directly, without a second wall-speed
    /// root. No film is published on cancellation, unsupported shear or an
    /// unrepresentable response. Temperature stays at the declared value.
    pub fn with_cancel(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauShearedFilmSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauSlitControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CarreauFilmError> {
        poll(&mut cancelled)?;
        if !positive(spec.density_kg_m3)
            || !positive(spec.thickness_m)
            || !spec.tangential_acceleration_m_s2.is_finite()
            || !spec.surface_shear_pa.is_finite()
        {
            return Err(FilmError::InvalidInput("Carreau sheared-film geometry or forcing").into());
        }
        let force = spec.density_kg_m3 * spec.tangential_acceleration_m_s2;
        if !force.is_finite() || (force == 0.0 && spec.tangential_acceleration_m_s2 != 0.0) {
            return Err(FilmError::NonRepresentable.into());
        }
        // Unit width/length are only the native chart for q per width and a
        // pressure gradient; they do not infer a physical vessel span.
        let slit = LaminarCarreauSlit::new(
            law,
            temperature_k,
            1.0,
            spec.thickness_m,
            1.0,
            maximum_laminar_shear_rate_per_s,
        )?;
        let profile = material_call(&mut cancelled, |check| {
            slit.profile_with_upper_traction_check(
                force,
                0.0,
                spec.surface_shear_pa,
                controls,
                check,
            )
        })?;
        let flow = profile.flow();
        let surface_velocity = flow.accepted_wall_velocities_m_s[1];
        let surface_work = spec.surface_shear_pa * surface_velocity;
        let work_residual = (spec.surface_shear_pa - flow.wall_tractions_pa[1]) * surface_velocity
            + flow.work_residual_w;

        // The slit tangents hold BOTH wall velocities fixed. Eliminating the
        // upper velocity at fixed traction adds the positive Schur-complement
        // term A²/R to pressure conductance. A = dQ/dUupper at fixed Ulower.
        // Unit width/length make the resulting values per-width film mobilities.
        let resistance = flow.relative_wall_resistance_n_s_m;
        let a = 0.5 * spec.thickness_m + flow.flow_relative_wall_velocity_derivative_m2;
        let surface_mobility = 1.0 / resistance;
        let flow_surface_mobility = a * surface_mobility;
        let flow_gradient_mobility = flow.pressure_conductance_m3_s_pa + a * flow_surface_mobility;
        if ![
            surface_work,
            work_residual,
            flow.volumetric_flow_m3_s / spec.thickness_m,
        ]
        .into_iter()
        .all(f64::is_finite)
            || ![
                surface_mobility,
                flow_surface_mobility,
                flow_gradient_mobility,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(FilmError::NonRepresentable.into());
        }
        let power = FilmPower {
            body_work_w_m2: flow.pressure_power_w,
            surface_work_w_m2: surface_work,
            viscous_dissipation_w_m2: flow.dissipated_power_w,
        };
        poll(&mut cancelled)?;
        Ok(Self {
            spec,
            profile,
            power,
            work_residual_w_m2: work_residual,
            surface_mobility_m_s_pa: surface_mobility,
            flow_surface_mobility_m2_s_pa: flow_surface_mobility,
            flow_gradient_mobility_m3_s_pa: flow_gradient_mobility,
        })
    }

    pub const fn spec(&self) -> CarreauShearedFilmSpec {
        self.spec
    }
    pub const fn law(&self) -> &ThermalCarreauYasuda {
        self.profile.slit().law()
    }
    pub const fn temperature_k(&self) -> f64 {
        self.profile.slit().temperature_k()
    }
    pub const fn thickness_m(&self) -> f64 {
        self.spec.thickness_m
    }
    pub const fn maximum_laminar_shear_rate_per_s(&self) -> f64 {
        self.profile.slit().maximum_laminar_shear_rate_per_s()
    }
    pub fn flow_per_width_m2_s(&self) -> f64 {
        self.profile.flow().volumetric_flow_m3_s
    }
    pub fn mean_velocity_m_s(&self) -> f64 {
        self.flow_per_width_m2_s() / self.spec.thickness_m
    }
    pub fn surface_velocity_m_s(&self) -> f64 {
        self.profile.flow().accepted_wall_velocities_m_s[1]
    }
    /// Signed fluid-on-wall traction [Pa]; the wall-on-fluid sign is opposite.
    pub fn wall_shear_pa(&self) -> f64 {
        -self.profile.flow().wall_tractions_pa[0]
    }
    pub const fn surface_shear_pa(&self) -> f64 {
        self.spec.surface_shear_pa
    }
    /// Accepted minus prescribed surface traction [Pa], including roundoff.
    pub fn surface_shear_residual_pa(&self) -> f64 {
        self.profile.flow().wall_tractions_pa[1] - self.spec.surface_shear_pa
    }
    pub const fn power(&self) -> FilmPower {
        self.power
    }
    /// Requested input work minus original-profile heat, evaluated without
    /// subtracting potentially large opposing body/surface power terms.
    pub const fn work_residual_w_m2(&self) -> f64 {
        self.work_residual_w_m2
    }
    pub fn quadrature_relative_change(&self) -> f64 {
        self.profile.flow().quadrature_relative_change
    }
    /// d(surface velocity)/d(surface traction), fixed body force [m/(s Pa)].
    pub const fn surface_velocity_traction_derivative_m_s_pa(&self) -> f64 {
        self.surface_mobility_m_s_pa
    }
    /// d(flux per width)/d(surface traction), fixed body force [m²/(s Pa)].
    /// This also equals d(surface velocity)/d(body force per volume).
    pub const fn flow_traction_derivative_m2_s_pa(&self) -> f64 {
        self.flow_surface_mobility_m2_s_pa
    }
    /// d(flux per width)/d(body force per volume), fixed traction [m³/(s Pa)].
    /// Multiply by density for the derivative with respect to acceleration.
    pub const fn flow_force_density_derivative_m3_s_pa(&self) -> f64 {
        self.flow_gradient_mobility_m3_s_pa
    }

    pub fn sample(&self, distance_from_wall_m: f64) -> Result<CarreauFilmSample, CarreauFilmError> {
        self.sample_with_cancel(distance_from_wall_m, || false)
    }

    pub fn sample_with_cancel(
        &self,
        distance_from_wall_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CarreauFilmSample, CarreauFilmError> {
        let sample = material_call(&mut cancelled, |check| {
            self.profile.sample_with_check(distance_from_wall_m, check)
        })?;
        let state = self
            .law()
            .evaluate(self.temperature_k(), sample.velocity_gradient_per_s.abs())?;
        let tangent = state.stress_tangent_pa_s();
        let curvature = -self.spec.density_kg_m3 * self.spec.tangential_acceleration_m_s2 / tangent;
        if !curvature.is_finite() {
            return Err(FilmError::NonRepresentable.into());
        }
        poll(&mut cancelled)?;
        Ok(CarreauFilmSample {
            velocity_m_s: sample.axial_velocity_m_s,
            velocity_gradient_per_s: sample.velocity_gradient_per_s,
            velocity_curvature_per_m_s: curvature,
            shear_stress_pa: sample.shear_stress_pa,
            dynamic_viscosity_pa_s: sample.dynamic_viscosity_pa_s,
            stress_tangent_pa_s: tangent,
            dissipation_w_m3: sample.dissipation_w_m3,
            quadrature_relative_change: sample.quadrature_relative_change,
        })
    }
}

impl CarreauFilm {
    /// Apply signed surface traction at this already resolved film thickness.
    /// The flux changes; the original prescribed-flow thickness is not solved again.
    pub fn with_surface_shear(
        &self,
        surface_shear_pa: f64,
        controls: CarreauSlitControls,
    ) -> Result<CarreauShearedFilm, CarreauFilmError> {
        CarreauShearedFilm::new(
            *self.law(),
            self.temperature_k(),
            CarreauShearedFilmSpec {
                density_kg_m3: self.spec().density_kg_m3,
                thickness_m: self.thickness_m(),
                tangential_acceleration_m_s2: self.spec().tangential_acceleration_m_s2,
                surface_shear_pa,
            },
            self.maximum_laminar_shear_rate_per_s(),
            controls,
        )
    }
}
