//! Material-bound cavitation on the shared liquid and film owners.
//!
//! Geometry supplies volume/force ports; one original liquid source supplies
//! density, viscosity and equilibrium vapor pressure. The native film solver
//! owns rupture, transported liquid content and reformation. The binding never
//! substitutes a pressure clamp, a gas EOS or a duplicate thermal integrator.
//!
//! This is a prescribed-configuration film step. A mechanical solver must use
//! the returned pressure forces inside its own coupled solve before committing
//! a candidate. Endpoint forces are not time-integrated mechanical work. The
//! isothermal, negligible-void-mass film approximation does not infer latent
//! heat from changing empty volume or deposit such energy into a thermal body.

#[cfg(feature = "cavitation")]
pub mod bubble;

/// Resonant coherent propagation from the native compressible bubble model.
#[cfg(feature = "cavitation")]
pub mod bubbly;

/// Event-resolved pure-vapor cavities on the original liquid duct state.
#[cfg(feature = "cavitation")]
pub mod duct;

use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_tribo::resistive_film::cavitation::{CavitatingFilm, CavitationLiquid};
pub use fs_tribo::resistive_film::cavitation::{
    CavitatingFilmState, CavitationReport, CavitationSolveConfig, CavitationStep,
};
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmError, FilmLimits, ResistiveFilm};

/// Absolute reservoir boundary and the native film's explicit numerical budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CavitationFilmEnvironment {
    pub ambient_pressure_pa: f64,
    /// Incoming liquid content. A partially filled reservoir must be at the
    /// same cavitation pressure as the film.
    pub ambient_liquid_fraction: f64,
    /// Largest admitted saturated-vapor/liquid density ratio, in (0, 0.1].
    /// Known vapor data must satisfy this negligible-void-mass approximation;
    /// absent optional data leaves its applicability as a caller assumption.
    pub maximum_vapor_density_ratio: f64,
    pub solve: CavitationSolveConfig,
}

/// Failure of the original material query or native film owner.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialCavitationFilmError {
    Material(MaterialError),
    Film(FilmError),
}
impl core::fmt::Display for MaterialCavitationFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "cavitation material: {error}"),
            Self::Film(error) => write!(f, "cavitation film: {error}"),
        }
    }
}
impl std::error::Error for MaterialCavitationFilmError {}
impl From<MaterialError> for MaterialCavitationFilmError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<FilmError> for MaterialCavitationFilmError {
    fn from(error: FilmError) -> Self {
        Self::Film(error)
    }
}

/// An original immutable source and the one liquid inventory solver it drives.
///
/// These material coefficients are frozen for this owner's lifetime. Temperature
/// feedback requires a conservative state transfer before adopting a refreshed
/// density or vapor pressure; rebuilding a model is not itself that transfer.
pub struct MaterialCavitatingFilm<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    film: CavitatingFilm,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialCavitatingFilm<'a, F> {
    /// Bind the source's own rho, mu and p_sat to the existing geometric graph.
    pub fn new(
        source: &'a F,
        cells: Vec<FilmCell>,
        channels: Vec<FilmChannel>,
        ports: usize,
        limits: FilmLimits,
        environment: CavitationFilmEnvironment,
    ) -> Result<Self, MaterialCavitationFilmError> {
        if !environment.maximum_vapor_density_ratio.is_finite()
            || environment.maximum_vapor_density_ratio <= 0.0
            || environment.maximum_vapor_density_ratio > 0.1
        {
            return Err(FilmError(
                "cavitating film maximum vapor/liquid density ratio must be finite and in (0, 0.1]",
            )
            .into());
        }
        let material = source.cavitation_state()?;
        if material.saturated_vapor().is_some_and(|vapor| {
            vapor.density_kg_m3() / material.liquid_density_kg_m3()
                > environment.maximum_vapor_density_ratio
        }) {
            return Err(FilmError(
                "saturated vapor density exceeds the negligible-void-mass film domain",
            )
            .into());
        }
        let geometry = ResistiveFilm::new(
            cells,
            channels,
            ports,
            material.dynamic_viscosity_pa_s(),
            limits,
        )?;
        let film = CavitatingFilm::new(
            geometry,
            CavitationLiquid {
                density_kg_m3: material.liquid_density_kg_m3(),
                cavitation_pressure_pa: material.saturation_pressure_pa(),
                ambient_pressure_pa: environment.ambient_pressure_pa,
                ambient_liquid_fraction: environment.ambient_liquid_fraction,
            },
            environment.solve,
        )?;
        Ok(Self {
            source,
            material,
            film,
        })
    }

    /// Original source, retaining its exact query and property identities.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Coherent frozen material state used by the native pressure equation.
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    /// Native owner, with no second film implementation.
    pub const fn film(&self) -> &CavitatingFilm {
        &self.film
    }

    /// Fully filled initialization at the declared absolute reservoir pressure.
    pub fn initial_state(
        &self,
        q: &[f64],
    ) -> Result<CavitatingFilmState, MaterialCavitationFilmError> {
        Ok(self.film.initial_state(q)?)
    }

    /// Explicit initial liquid content, including declared starvation.
    pub fn initialize(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
    ) -> Result<CavitatingFilmState, MaterialCavitationFilmError> {
        Ok(self.film.initialize(q, liquid_fraction)?)
    }

    /// Propose a conservative native step without mutating accepted state.
    ///
    /// Entrainment is signed full-film Couette volume flow on each original
    /// channel. Upwind liquid content is applied by the native owner exactly
    /// once. A `Cx` caller can supply its cancellation check through `cancelled`.
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        channel_entrainment_m3_s: &[f64],
        dt_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CavitationStep, MaterialCavitationFilmError> {
        Ok(self
            .film
            .advance(state, q_next, channel_entrainment_m3_s, dt_s, cancelled)?)
    }
}
