//! Material-bound cavitation on the shared liquid and film owners.
//!
//! Geometry supplies volume/force ports; one original liquid source supplies
//! density, viscosity and equilibrium vapor pressure. The native film solver
//! owns rupture, transported liquid content and reformation. The binding never
//! substitutes a pressure clamp, a gas EOS or a duplicate thermal integrator.
//!
//! The film binding below is a prescribed-configuration step. A mechanical
//! solver must use the returned pressure forces inside its own coupled solve before committing
//! a candidate. Endpoint forces are not time-integrated mechanical work. The
//! negligible-void-mass film approximation does not infer latent heat from
//! changing empty volume or deposit such energy into a thermal body. Explicit
//! source-temperature continuation refreshes the original material claims and
//! transfers the retained mass through the same native pressure/content solve;
//! the caller's thermal owner supplies the proposed temperature and energy.
//! The separate water-cell owner evolves finite mass and internal energy with
//! actual IF97 phase inventories; the duct and bubble owners retain their own
//! explicitly declared isothermal or frozen-material approximations.

#[cfg(feature = "cavitation")]
pub mod bubble;

/// Interacting radial sources and their individually retarded acoustic paths.
#[cfg(feature = "cavitation")]
pub mod cluster;

/// Resonant coherent propagation from the native compressible bubble model.
#[cfg(feature = "cavitation")]
pub mod bubbly;

/// Event-resolved pure-vapor cavities on the original liquid duct state.
#[cfg(feature = "cavitation")]
pub mod duct;

/// Finite water mass and internal energy with native IF97 flashing/condensation.
/// This thermal inventory uses existing material owners without a bubble solver.
pub mod water_cell;

/// Finite solvent heat and the original two-film dissolved-gas kinetics.
pub mod dissolved_gas;

/// Source-temperature shifts of nonlinear film transport on retained liquid mass.
pub mod power_law;
pub use power_law::{MaterialCavitatingPowerLawContinuation, MaterialCavitatingPowerLawFilm};

/// Pressure-dependent liquid storage with original-source thermal continuation.
pub mod compressible_film;
pub use compressible_film::{
    MaterialCompressibleCavitatingFilm, MaterialCompressibleCavitationContinuation,
};

use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState, ResolvedCavitationState};
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
/// Coefficients stay frozen within each accepted step. A resolved source can
/// continue through `advance_with_temperature_source`, preserving the old mass
/// while solving with the new density, viscosity and vapor pressure. Rebuilding
/// a model or reinitializing filling is not that conservative transition.
pub struct MaterialCavitatingFilm<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    film: CavitatingFilm,
    maximum_vapor_density_ratio: f64,
}

/// Refreshed source binding and its matching conservative native candidate.
/// Accept both together only after the surrounding thermal/mechanical proposal
/// accepts. The source reference retains the complete endpoint material bundle.
pub struct MaterialCavitationContinuation<'a, F: CavitationProperties + ?Sized> {
    pub film: MaterialCavitatingFilm<'a, F>,
    pub step: CavitationStep,
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
        validate_vapor_density(material, environment.maximum_vapor_density_ratio)?;
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
            maximum_vapor_density_ratio: environment.maximum_vapor_density_ratio,
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

impl MaterialCavitatingFilm<'_, ResolvedCavitationState> {
    /// Continue the accepted liquid inventory at a proposed source temperature.
    ///
    /// Resolve `source` through the original `ResolvedCavitationState::at_query_point`
    /// and immutable material card. This binding admits only that same card,
    /// complete source claims, property requirements and temperature axis; every
    /// nonthermal coordinate and quantity convention must remain unchanged.
    /// Dynamic viscosity, density and saturation pressure refresh together;
    /// kinematic-source viscosity uses the NEW density. The original dilute
    /// vapor cap, geometry, reservoir boundary and numerical limits remain fixed.
    ///
    /// The native endpoint solve retains `rho_old*V_old*theta_old`, then uses
    /// the new material for the candidate pressure, filling and face mass. No
    /// accepted source, film or liquid state is mutated. Accept the returned
    /// binding and state together after the original thermal/mechanical proposal
    /// succeeds; a failure of any source, mass, domain or cancellation gate
    /// leaves the complete accepted pair usable for retry.
    ///
    /// This is first-order endpoint constitutive splitting with caller-owned
    /// time refinement. The proposed temperature must come from the original
    /// thermal owner or prescribed boundary data. No heat capacity, latent heat,
    /// vapor mass, dissipated work or integrated mechanical passivity is inferred
    /// from a viscosity/density change or from the endpoint pressure force.
    pub fn advance_with_temperature_source<'b>(
        &self,
        source: &'b ResolvedCavitationState,
        state: &CavitatingFilmState,
        q_next: &[f64],
        channel_entrainment_m3_s: &[f64],
        dt_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        MaterialCavitationContinuation<'b, ResolvedCavitationState>,
        MaterialCavitationFilmError,
    > {
        require_temperature_continuation(self.source, source)?;
        let material = source.state();
        validate_vapor_density(material, self.maximum_vapor_density_ratio)?;
        let next_liquid = CavitationLiquid {
            density_kg_m3: material.liquid_density_kg_m3(),
            cavitation_pressure_pa: material.saturation_pressure_pa(),
            ..self.film.liquid()
        };
        let viscosities = vec![material.dynamic_viscosity_pa_s(); self.film.channel_count()];
        let next = self.film.advance_with_liquid(
            state,
            q_next,
            channel_entrainment_m3_s,
            dt_s,
            next_liquid,
            &viscosities,
            cancelled,
        )?;
        Ok(MaterialCavitationContinuation {
            film: MaterialCavitatingFilm {
                source,
                material,
                film: next.film,
                maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
            },
            step: next.step,
        })
    }
}

fn require_temperature_continuation(
    previous: &ResolvedCavitationState,
    next: &ResolvedCavitationState,
) -> Result<(), FilmError> {
    let old = previous.material();
    let new = next.material();
    let axis = previous.temperature_axis();
    if previous.requirements() != next.requirements()
        || axis != next.temperature_axis()
        || old.card_identity() != new.card_identity()
        || old.axis_quantities() != new.axis_quantities()
        || old.query_point().len() != new.query_point().len()
        || old.query_point().iter().zip(new.query_point()).any(
            |((old_axis, old_value), (new_axis, new_value))| {
                old_axis != new_axis
                    || (old_axis != axis && old_value.to_bits() != new_value.to_bits())
            },
        )
        || old.properties().len() != new.properties().len()
        || old.properties().iter().zip(new.properties()).any(|(a, b)| {
            a.requirement() != b.requirement()
                || a.answer().receipt.selected != b.answer().receipt.selected
        })
    {
        return Err(FilmError(
            "cavitating-film temperature continuation requires the original card, claims and fixed context",
        ));
    }
    Ok(())
}

fn validate_vapor_density(
    material: CavitationState,
    maximum_vapor_density_ratio: f64,
) -> Result<(), FilmError> {
    if material.saturated_vapor().is_some_and(|vapor| {
        vapor.density_kg_m3() / material.liquid_density_kg_m3() > maximum_vapor_density_ratio
    }) {
        return Err(FilmError(
            "saturated vapor density exceeds the negligible-void-mass film domain",
        ));
    }
    Ok(())
}
