//! Fully developed films of the original thermal Carreau-Yasuda law.
//!
//! [`CarreauShearedFilm`] resolves signed surface traction at fixed thickness
//! or prescribed signed flux on an explicitly selected thickness branch;
//! [`CarreauFilm`] below resolves a shear-free film at prescribed positive flux.
//!
//! Momentum fixes tau(y)=rho*g_t*(h-y). A symmetric, stationary-wall slit
//! of gap 2h and pressure gradient rho*g_t has exactly this lower-half flow.
//! The original material slit owner supplies its constitutive inverse, positive
//! quadrature, local viscosity and heat. Only the physical thickness is solved
//! here, using dq/dh=h*gamma_wall. No effective Newtonian viscosity is fitted.

use super::{FilmError, FilmMaterial, NusseltFilm};
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;
use fs_material::fluid_transport::rheology::slit::{
    CarreauSlitControls, CarreauSlitProfile, LaminarCarreauSlit,
};

mod sheared;
pub use sheared::{
    CarreauShearedFilm, CarreauShearedFilmFlowSolution, CarreauShearedFilmFlowSpec,
    CarreauShearedFilmSpec,
};

/// Actual SI forcing and liquid density/capillarity; viscosity belongs to the law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauFilmSpec {
    pub density_kg_m3: f64,
    pub surface_tension_n_m: f64,
    pub flow_per_width_m2_s: f64,
    pub tangential_acceleration_m_s2: f64,
    pub normal_acceleration_m_s2: f64,
}

/// Native constitutive/quadrature work plus a bounded physical thickness solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauFilmControls {
    pub slit: CarreauSlitControls,
    /// Relative accepted flux residual AND thickness Newton correction.
    pub relative_flow_tolerance: f64,
    /// Maximum forward profiles, including the initial bracket, up to 128.
    /// A sheared-film inverse requires at least two for its supplied endpoints;
    /// the shear-free inverse permits one initial upper-bracket profile.
    pub maximum_thickness_iterations: u32,
}
impl Default for CarreauFilmControls {
    fn default() -> Self {
        Self {
            slit: CarreauSlitControls::default(),
            relative_flow_tolerance: 1.0e-9,
            maximum_thickness_iterations: 96,
        }
    }
}

/// No partial film is published on source, work, support or cancellation refusal.
#[derive(Debug)]
pub enum CarreauFilmError {
    Film(FilmError),
    Material(MaterialError),
    OutsideRateSupport,
    /// The supplied positive thickness interval crosses a wall-stress zero,
    /// so its flux is not monotone. Choose one physical thickness branch.
    NonMonotoneThicknessInterval {
        wall_shear_at_bounds_pa: [f64; 2],
    },
    /// The requested signed flux is outside the selected monotone branch.
    FlowOutsideThicknessInterval {
        requested_m2_s: f64,
        flow_at_bounds_m2_s: [f64; 2],
    },
    UnresolvedFlow {
        iterations: u32,
        residual_m2_s: f64,
    },
    Cancelled,
}
impl core::fmt::Display for CarreauFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Carreau film refused: {self:?}")
    }
}
impl std::error::Error for CarreauFilmError {}
impl From<FilmError> for CarreauFilmError {
    fn from(value: FilmError) -> Self {
        Self::Film(value)
    }
}
impl From<MaterialError> for CarreauFilmError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}

/// Original accepted nonlinear profile, law, temperature and source-rate envelope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauFilm {
    spec: CarreauFilmSpec,
    profile: CarreauSlitProfile,
    flow_iterations: u32,
    surface_velocity_m_s: f64,
}

/// One physical point in the same base flow, before nondimensionalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauFilmSample {
    pub velocity_m_s: f64,
    pub velocity_gradient_per_s: f64,
    /// d²u/dy² [1/(m s)], from original momentum and differential viscosity.
    pub velocity_curvature_per_m_s: f64,
    pub shear_stress_pa: f64,
    pub dynamic_viscosity_pa_s: f64,
    /// d(mu*gamma)/d(gamma) [Pa s]; distinct from apparent viscosity.
    pub stress_tangent_pa_s: f64,
    pub dissipation_w_m3: f64,
    pub quadrature_relative_change: f64,
}

impl CarreauFilm {
    pub fn at_flow(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauFilmSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauFilmControls,
    ) -> Result<Self, CarreauFilmError> {
        Self::at_flow_with_cancel(
            law,
            temperature_k,
            spec,
            maximum_laminar_shear_rate_per_s,
            controls,
            || false,
        )
    }

    /// Resolve h at the prescribed positive flux, retaining the actual accepted
    /// flux separately. The law must support zero shear through the declared
    /// rate ceiling. Temperature is fixed and no source support is extended.
    ///
    /// A zero-shear-viscosity Nusselt thickness bounds the shear-thinning root
    /// from above. The original maximum shear stress further bounds that
    /// bracket. Every forward profile and point sample uses the original
    /// coupled slit law; its symmetric zero-stress plane is the free surface.
    pub fn at_flow_with_cancel(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        spec: CarreauFilmSpec,
        maximum_laminar_shear_rate_per_s: f64,
        controls: CarreauFilmControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CarreauFilmError> {
        poll(&mut cancelled)?;
        if !positive(spec.density_kg_m3)
            || !positive(controls.relative_flow_tolerance)
            || controls.relative_flow_tolerance >= 1.0
            || !(1..=128).contains(&controls.maximum_thickness_iterations)
        {
            return Err(
                FilmError::InvalidInput("Carreau film density or thickness controls").into(),
            );
        }
        // Unit width/length are only the original slit chart for a flux per
        // width and a pressure gradient. They do not infer a vessel span.
        let slit = LaminarCarreauSlit::new(
            law,
            temperature_k,
            1.0,
            1.0,
            1.0,
            maximum_laminar_shear_rate_per_s,
        )?;
        let reference = NusseltFilm::new(
            FilmMaterial {
                density_kg_m3: spec.density_kg_m3,
                kinematic_viscosity_m2_s: slit.zero_shear_viscosity_pa_s() / spec.density_kg_m3,
                surface_tension_n_m: spec.surface_tension_n_m,
            },
            spec.flow_per_width_m2_s,
            spec.tangential_acceleration_m_s2,
            spec.normal_acceleration_m_s2,
        )?;
        let force = spec.density_kg_m3 * spec.tangential_acceleration_m_s2;
        let ceiling = law.evaluate(temperature_k, slit.maximum_laminar_shear_rate_per_s())?;
        let maximum_stress = ceiling.dynamic_viscosity_pa_s() * ceiling.shear_rate_per_s();
        let support_height = (maximum_stress / force).next_down();
        if !positive(force) || !positive(support_height) {
            return Err(FilmError::NonRepresentable.into());
        }
        let mut high = reference.thickness_m().min(support_height);
        let mut low = 0.0;
        let mut height = high;
        let mut last_residual = 0.0;
        for iteration in 1..=controls.maximum_thickness_iterations {
            poll(&mut cancelled)?;
            let passage = slit.with_gap(2.0 * height)?;
            let profile = material_call(&mut cancelled, |check| {
                passage.profile_with_check(force, [0.0, 0.0], controls.slit, check)
            })?;
            if profile.stress_offset_pa() != 0.0 {
                return Err(FilmError::InvalidInput(
                    "symmetric slit must retain zero midpoint stress",
                )
                .into());
            }
            let flux = 0.5 * profile.flow().pressure_flow_m3_s;
            let wall = material_call(&mut cancelled, |check| {
                profile.sample_with_check(0.0, check)
            })?;
            let derivative = height * wall.velocity_gradient_per_s;
            if !positive(flux) || !positive(derivative) {
                return Err(FilmError::NonRepresentable.into());
            }
            let residual = flux - spec.flow_per_width_m2_s;
            last_residual = residual;
            let correction = residual / derivative;
            if (residual / spec.flow_per_width_m2_s).abs() <= controls.relative_flow_tolerance
                && (correction / height).abs() <= controls.relative_flow_tolerance
            {
                let surface = material_call(&mut cancelled, |check| {
                    profile.sample_with_check(height, check)
                })?;
                let film = Self {
                    spec,
                    profile,
                    flow_iterations: iteration,
                    surface_velocity_m_s: surface.axial_velocity_m_s,
                };
                if !positive(film.mean_velocity_m_s())
                    || !positive(film.surface_velocity_m_s())
                    || !positive(film.reference_reynolds())
                    || !film.dissipation_w_m2().is_finite()
                {
                    return Err(FilmError::NonRepresentable.into());
                }
                poll(&mut cancelled)?;
                return Ok(film);
            }
            if iteration == 1 && residual < 0.0 {
                return Err(CarreauFilmError::OutsideRateSupport);
            }
            if residual > 0.0 {
                high = height;
            } else {
                low = height;
            }
            let candidate = height - correction;
            height = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                f64::midpoint(low, high)
            };
            if !positive(height) {
                return Err(FilmError::NonRepresentable.into());
            }
        }
        Err(CarreauFilmError::UnresolvedFlow {
            iterations: controls.maximum_thickness_iterations,
            residual_m2_s: last_residual,
        })
    }

    pub const fn spec(&self) -> CarreauFilmSpec {
        self.spec
    }
    pub const fn law(&self) -> &ThermalCarreauYasuda {
        self.profile.slit().law()
    }
    pub const fn temperature_k(&self) -> f64 {
        self.profile.slit().temperature_k()
    }
    pub const fn maximum_laminar_shear_rate_per_s(&self) -> f64 {
        self.profile.slit().maximum_laminar_shear_rate_per_s()
    }
    pub const fn zero_shear_viscosity_pa_s(&self) -> f64 {
        self.profile.slit().zero_shear_viscosity_pa_s()
    }
    pub fn thickness_m(&self) -> f64 {
        0.5 * self.profile.slit().gap_m()
    }
    pub fn flow_per_width_m2_s(&self) -> f64 {
        0.5 * self.profile.flow().pressure_flow_m3_s
    }
    pub fn flow_residual_m2_s(&self) -> f64 {
        self.flow_per_width_m2_s() - self.spec.flow_per_width_m2_s
    }
    pub const fn flow_iterations(&self) -> u32 {
        self.flow_iterations
    }
    pub fn mean_velocity_m_s(&self) -> f64 {
        self.flow_per_width_m2_s() / self.thickness_m()
    }
    pub const fn surface_velocity_m_s(&self) -> f64 {
        self.surface_velocity_m_s
    }
    /// Bookkeeping Reynolds number based on the original zero-shear viscosity;
    /// this is not a universal generalized-Newtonian transition criterion.
    pub fn reference_reynolds(&self) -> f64 {
        self.spec.density_kg_m3 * self.flow_per_width_m2_s() / self.zero_shear_viscosity_pa_s()
    }
    /// Native viscous heat per wall area [W/m²], not deposited into a thermal model.
    pub fn dissipation_w_m2(&self) -> f64 {
        0.5 * self.profile.flow().dissipated_power_w
    }
    pub fn quadrature_relative_change(&self) -> f64 {
        self.profile.flow().quadrature_relative_change
    }

    pub fn sample(&self, distance_from_wall_m: f64) -> Result<CarreauFilmSample, CarreauFilmError> {
        self.sample_with_cancel(distance_from_wall_m, || false)
    }

    /// Resolve the original apparent and differential viscosities at this SAME
    /// local shear. The curvature follows mu_t*u''=-rho*g_t exactly.
    pub fn sample_with_cancel(
        &self,
        distance_from_wall_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CarreauFilmSample, CarreauFilmError> {
        if !distance_from_wall_m.is_finite()
            || distance_from_wall_m < 0.0
            || distance_from_wall_m > self.thickness_m()
        {
            return Err(FilmError::InvalidInput("Carreau film sample outside wall/surface").into());
        }
        let sample = material_call(&mut cancelled, |check| {
            self.profile.sample_with_check(distance_from_wall_m, check)
        })?;
        let state = self
            .law()
            .evaluate(self.temperature_k(), sample.velocity_gradient_per_s.abs())?;
        let tangent = state.stress_tangent_pa_s();
        let curvature = -self.spec.density_kg_m3 * self.spec.tangential_acceleration_m_s2 / tangent;
        if !curvature.is_finite() || curvature >= 0.0 {
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

fn material_call<T>(
    cancelled: &mut impl FnMut() -> bool,
    action: impl FnOnce(&mut dyn FnMut() -> Result<(), MaterialError>) -> Result<T, MaterialError>,
) -> Result<T, CarreauFilmError> {
    let mut observed = false;
    let result = action(&mut || {
        if cancelled() {
            observed = true;
            Err(MaterialError::State {
                what: "Carreau film cancelled".into(),
            })
        } else {
            Ok(())
        }
    });
    if observed {
        Err(CarreauFilmError::Cancelled)
    } else {
        result.map_err(CarreauFilmError::Material)
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CarreauFilmError> {
    if cancelled() {
        Err(CarreauFilmError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
