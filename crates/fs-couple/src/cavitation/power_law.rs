//! Source-temperature continuation of the original nonlinear cavitating film.
//!
//! The authored reference rheology uses an explicit separable approximation:
//! `K(T) = K_ref * mu_source(T) / mu_source(T_ref)`, at fixed n and core rate.
//! The pinned source supplies density and saturation pressure together with
//! that viscosity ratio. Its Newtonian viscosity is never mistaken for K,
//! and it does not establish the independent rheological shape. The native
//! liquid owner conserves OLD mass at each proposed endpoint material state.

use fs_material::cavitation::{CavitationProperties, CavitationState, ResolvedCavitationState};
use fs_math::det;
use fs_tribo::resistive_film::FilmError;
use fs_tribo::resistive_film::cavitation::{
    CavitatingFilmState, CavitatingPowerLawFilm, CavitationLiquid, CavitationReservoir,
    CavitationStep,
};
use fs_tribo::resistive_film::power_law::{PowerLawControls, RegularizedPowerLaw};

use super::{
    MaterialCavitatingFilm, MaterialCavitationFilmError, require_temperature_continuation,
    validate_vapor_density,
};

/// One retained liquid source and its conservative nonlinear transport owner.
///
/// Reference laws are authored at the original source temperature. The complete
/// reference material values and those laws survive every accepted continuation;
/// K is always recomputed from that reference, never from a previous shifted K.
/// The endpoint source retains its own pinned property and query identities.
///
/// The caller must admit the separable shift over the temperature/shear range.
/// Density remains spatially homogeneous, n and the physical core rate stay
/// fixed, and the thin-gap creeping-flow assumptions still apply. No cavitated
/// tangential traction, heat partition, latent heat, thermal time integration,
/// rheological memory or trapped-liquid compressibility is introduced.
pub struct MaterialCavitatingPowerLawFilm<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    reference_material: CavitationState,
    reference_laws: Vec<RegularizedPowerLaw>,
    film: CavitatingPowerLawFilm,
    maximum_vapor_density_ratio: f64,
}

/// Matching source/constitutive owner and native conservative liquid candidate.
/// The surrounding mechanical and thermal proposal must accept both together.
pub struct MaterialCavitatingPowerLawContinuation<'a, F: CavitationProperties + ?Sized> {
    pub film: MaterialCavitatingPowerLawFilm<'a, F>,
    pub step: CavitationStep,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialCavitatingFilm<'a, F> {
    /// Bind explicit SI rheology at this source's CURRENT reference temperature.
    ///
    /// The source density, vapor limit, numerical budgets, original reservoirs
    /// and accepted liquid state identity remain unchanged. Native conversion
    /// admits the channel laws and refuses incompatible Barus viscosity.
    /// `consistency_pa_s_n` is independent caller data; the source viscosity
    /// supplies only its subsequent dimensionless temperature shift.
    pub fn with_power_law(
        &self,
        laws: Vec<RegularizedPowerLaw>,
        controls: PowerLawControls,
    ) -> Result<MaterialCavitatingPowerLawFilm<'a, F>, MaterialCavitationFilmError> {
        let film = self.film.with_power_law(laws.clone(), controls)?;
        Ok(MaterialCavitatingPowerLawFilm {
            source: self.source,
            material: self.material,
            reference_material: self.material,
            reference_laws: laws,
            film,
            maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
        })
    }
}

impl<'a, F: CavitationProperties + ?Sized> MaterialCavitatingPowerLawFilm<'a, F> {
    /// Current original source and all of its property/query identities.
    pub const fn source(&self) -> &'a F {
        self.source
    }

    /// Coherent endpoint source values used by the current native film.
    pub const fn material(&self) -> CavitationState {
        self.material
    }

    /// Original material values defining T_ref and the source viscosity ratio.
    pub const fn reference_material(&self) -> CavitationState {
        self.reference_material
    }

    /// Original SI K_ref, n and core shear rate, in original channel order.
    pub fn reference_laws(&self) -> &[RegularizedPowerLaw] {
        &self.reference_laws
    }

    /// Native nonlinear transport and conservative liquid inventory owner.
    pub const fn film(&self) -> &CavitatingPowerLawFilm {
        &self.film
    }

    /// Rebind prescribed supplies while retaining reference and current source.
    /// Incoming reservoir liquid has the same homogeneous endpoint density.
    pub fn with_reservoirs(
        &self,
        reservoirs: &[Option<CavitationReservoir>],
    ) -> Result<Self, MaterialCavitationFilmError> {
        Ok(Self {
            source: self.source,
            material: self.material,
            reference_material: self.reference_material,
            reference_laws: self.reference_laws.clone(),
            film: self.film.with_reservoirs(reservoirs)?,
            maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
        })
    }

    pub fn initial_state(
        &self,
        q: &[f64],
    ) -> Result<CavitatingFilmState, MaterialCavitationFilmError> {
        Ok(self.film.initial_state(q)?)
    }

    pub fn initialize(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
    ) -> Result<CavitatingFilmState, MaterialCavitationFilmError> {
        Ok(self.film.initialize(q, liquid_fraction)?)
    }

    /// Native mass/normal-force proposal with the currently accepted source.
    /// Mean entrainment and relative sliding use original passage orientation.
    /// Endpoint normal force is not an integrated mechanical-work receipt.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CavitationStep, MaterialCavitationFilmError> {
        Ok(self.film.advance(
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            cancelled,
        )?)
    }
}

impl MaterialCavitatingPowerLawFilm<'_, ResolvedCavitationState> {
    /// Propose new source temperature, nonlinear viscosity and conserved mass.
    ///
    /// Only the same pinned card/claims and unchanged nonthermal query are
    /// admitted through the original material-continuation check. For a
    /// kinematic-viscosity source, its endpoint DYNAMIC viscosity uses the
    /// endpoint density before forming the reference ratio. The same endpoint
    /// density and saturation pressure enter the native conservative solve.
    ///
    /// This first-order endpoint split requires a temperature proposal and
    /// refinement from the original thermal owner or prescribed boundary.
    /// No density/viscosity change is priced as heat. Old liquid mass, source
    /// and rheology remain usable after any failed proposal or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_temperature_source<'b>(
        &self,
        source: &'b ResolvedCavitationState,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MaterialCavitatingPowerLawContinuation<'b, ResolvedCavitationState>,
        MaterialCavitationFilmError,
    > {
        if cancelled() {
            return Err(FilmError("material power-law cavitation step cancelled").into());
        }
        require_temperature_continuation(self.source, source)?;
        let material = source.state();
        validate_vapor_density(material, self.maximum_vapor_density_ratio)?;
        let reference_viscosity = self.reference_material.dynamic_viscosity_pa_s();
        let viscosity = material.dynamic_viscosity_pa_s();
        let laws = self
            .reference_laws
            .iter()
            .map(|reference| {
                // Preserve the exact native Newtonian source image when its
                // independently authored K_ref is that same dynamic viscosity.
                let direct = if reference.flow_index == 1.0
                    && reference.consistency_pa_s_n == reference_viscosity
                {
                    viscosity
                } else {
                    reference.consistency_pa_s_n * (viscosity / reference_viscosity)
                };
                let consistency = if direct.is_finite() && direct > 0.0 {
                    direct
                } else {
                    det::exp(
                        det::ln(reference.consistency_pa_s_n) + det::ln(viscosity)
                            - det::ln(reference_viscosity),
                    )
                };
                if !consistency.is_finite() || consistency <= 0.0 {
                    return Err(FilmError(
                        "temperature-shifted cavitating consistency is unrepresentable",
                    ));
                }
                Ok(RegularizedPowerLaw {
                    consistency_pa_s_n: consistency,
                    ..*reference
                })
            })
            .collect::<Result<Vec<_>, FilmError>>()?;
        let next = self.film.advance_with_liquid(
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            CavitationLiquid {
                density_kg_m3: material.liquid_density_kg_m3(),
                cavitation_pressure_pa: material.saturation_pressure_pa(),
                ..self.film.liquid()
            },
            laws,
            cancelled,
        )?;
        Ok(MaterialCavitatingPowerLawContinuation {
            film: MaterialCavitatingPowerLawFilm {
                source,
                material,
                reference_material: self.reference_material,
                reference_laws: self.reference_laws.clone(),
                film: next.film,
                maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
            },
            step: next.step,
        })
    }
}
