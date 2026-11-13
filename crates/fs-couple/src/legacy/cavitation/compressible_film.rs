//! Original-source viscosity and thermal density on native compressible films.
//!
//! The source supplies rho(T, p_ref), mu(T, p_ref) and saturation pressure.
//! An independently authored, positive ISOTHERMAL bulk modulus supplies the
//! explicit barotropic closure rho(T,p) = rho(T,p_ref)*exp((p-p_ref)/B).
//! Source sound speed is not silently converted to an isothermal modulus.
//! The pressure coordinate must be named and typed; a saturation density must
//! not be mistaken for density at the film's ambient reference pressure.

use fs_material::cavitation::{CavitationState, ResolvedCavitationState};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_tribo::resistive_film::FilmError;
use fs_tribo::resistive_film::cavitation::{
    CavitatingFilmState, CavitationLiquid, CavitationReservoir, CompressibleCavitatingFilm,
    CompressibleCavitationStep,
};

use super::{
    MaterialCavitatingFilm, MaterialCavitationFilmError, require_temperature_continuation,
    validate_vapor_density,
};

/// One immutable resolved source and the native pressure-dependent inventory.
///
/// The reference pressure remains fixed during source-temperature continuation.
/// The caller can explicitly supply a new independently authored B. The
/// source's local pressure validity does not establish this separate
/// barotropic law over every film pressure. Its
/// admissible range and temporal refinement remain caller-owned, along with
/// the thin-gap creeping-flow and negligible-void-mass approximations.
pub struct MaterialCompressibleCavitatingFilm<'a> {
    source: &'a ResolvedCavitationState,
    material: CavitationState,
    pressure_axis: String,
    film: CompressibleCavitatingFilm,
    maximum_vapor_density_ratio: f64,
}

/// Matching endpoint source/film and its conservative native candidate.
/// Accept both together after the surrounding thermal/mechanical proposal.
pub struct MaterialCompressibleCavitationContinuation<'a> {
    pub film: MaterialCompressibleCavitatingFilm<'a>,
    pub step: CompressibleCavitationStep,
}

impl<'a> MaterialCavitatingFilm<'a, ResolvedCavitationState> {
    /// Bind an independently authored isothermal bulk modulus in pascals.
    ///
    /// `pressure_axis` must identify an explicit static semantic Pressure
    /// coordinate in the original source query, equal to the native ambient
    /// pressure. The source density is then rho(T,p_ref), not rho at saturation
    /// or an unspecified pressure. All original viscosities, pressure supplies,
    /// geometry and numerical limits pass to the native compressible owner.
    ///
    /// Initialize through the returned owner: incompressible state identities
    /// are deliberately incompatible with pressure-dependent stored mass.
    /// Neither source sound speed nor Newtonian viscosity determines B.
    pub fn with_isothermal_bulk_modulus(
        &self,
        pressure_axis: &str,
        bulk_modulus_pa: f64,
    ) -> Result<MaterialCompressibleCavitatingFilm<'a>, MaterialCavitationFilmError> {
        require_reference_pressure(
            self.source,
            pressure_axis,
            self.film.liquid().ambient_pressure_pa,
        )?;
        let film = CompressibleCavitatingFilm::new(self.film.clone(), bulk_modulus_pa)?;
        validate_rupture_vapor_density(self.material, &film, self.maximum_vapor_density_ratio)?;
        Ok(MaterialCompressibleCavitatingFilm {
            source: self.source,
            material: self.material,
            pressure_axis: pressure_axis.into(),
            film,
            maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
        })
    }
}

impl<'a> MaterialCompressibleCavitatingFilm<'a> {
    /// Original source with its exact query, selected claims and card identity.
    pub const fn source(&self) -> &'a ResolvedCavitationState {
        self.source
    }

    /// Coherent source values at the explicitly retained reference pressure.
    pub const fn material(&self) -> CavitationState {
        self.material
    }

    /// Original named pressure coordinate defining rho(T,p_ref).
    pub fn pressure_axis(&self) -> &str {
        &self.pressure_axis
    }

    /// Native barotropic transport and pressure-dependent mass storage owner.
    pub const fn film(&self) -> &CompressibleCavitatingFilm {
        &self.film
    }

    /// Retain the source/reference pressure while prescribing absolute supplies.
    /// Each supply's density follows the same explicit barotropic liquid law.
    pub fn with_reservoirs(
        &self,
        reservoirs: &[Option<CavitationReservoir>],
    ) -> Result<Self, MaterialCavitationFilmError> {
        Ok(Self {
            source: self.source,
            material: self.material,
            pressure_axis: self.pressure_axis.clone(),
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

    /// Explicit initial compression; partial cells must be at saturation.
    pub fn initialize_with_pressure(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
        pressure_pa: &[f64],
    ) -> Result<CavitatingFilmState, MaterialCavitationFilmError> {
        Ok(self
            .film
            .initialize_with_pressure(q, liquid_fraction, pressure_pa)?)
    }

    /// Native mass/normal-force proposal with currently accepted source data.
    /// The native report carries one mass flux for each compressible passage.
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        channel_entrainment_m3_s: &[f64],
        dt_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompressibleCavitationStep, MaterialCavitationFilmError> {
        Ok(self
            .film
            .advance(state, q_next, channel_entrainment_m3_s, dt_s, cancelled)?)
    }

    /// Continue the original COMPRESSED mass at a proposed source temperature.
    ///
    /// The same pinned card, claims, quantity conventions and all nonthermal
    /// coordinates must remain unchanged. Density, dynamic viscosity and
    /// saturation pressure refresh together; kinematic-source conversion uses
    /// the new reference density. The original dilute-vapor cap is rechecked
    /// against the native liquid density at saturation, not its larger density
    /// at the reference pressure.
    /// B stays fixed and supplies retain their absolute pressure and filling.
    ///
    /// Old mass is reconstructed with OLD pressure-dependent density. The
    /// native endpoint solve uses new source density and viscosity without
    /// reinitializing pressure or liquid filling. In a sealed full cell at
    /// constant volume this yields p_new-p_old = B*ln(rho_ref,old/rho_ref,new).
    ///
    /// This first-order endpoint constitutive split is not a thermal clock:
    /// temperature, energy and refinement belong to the original thermal owner.
    /// It supplies no latent heat or compressive-energy budget, and endpoint
    /// normal forces do not establish integrated mechanical work. Source,
    /// domain, mass and cancellation refusals preserve the accepted pair.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_temperature_source<'b>(
        &self,
        source: &'b ResolvedCavitationState,
        state: &CavitatingFilmState,
        q_next: &[f64],
        channel_entrainment_m3_s: &[f64],
        dt_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCompressibleCavitationContinuation<'b>, MaterialCavitationFilmError> {
        self.advance_with_temperature_source_and_bulk_modulus(
            source,
            state,
            q_next,
            channel_entrainment_m3_s,
            dt_s,
            self.film.bulk_modulus_pa(),
            cancelled,
        )
    }

    /// Continue the original source and compressed mass with an explicitly
    /// supplied endpoint isothermal bulk modulus in pascals.
    ///
    /// The pinned source and pressure-reference requirements of
    /// [`Self::advance_with_temperature_source`] still apply. Density,
    /// viscosity and saturation pressure come from that source, while the
    /// caller supplies `bulk_modulus_pa` from an independent constitutive law.
    /// Adiabatic source sound speed is never used to infer this isothermal B.
    ///
    /// Stored OLD mass is evaluated with the OLD B; endpoint storage, reservoir
    /// density and pressure-driven mass transport use the NEW B. A changed B
    /// produces a distinct inventory owner even if all source properties agree.
    /// Sealed pressure can cross saturation, so the native complementarity
    /// solve retains mass by rupturing or refilling the film as necessary.
    /// The dilute-vapor check uses the NEW liquid density at rupture.
    ///
    /// This endpoint constitutive split supplies no latent or compressive
    /// energy budget. Accept the returned source, owner and state together;
    /// rejected sources, constitutive data or cancellation leave them intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_temperature_source_and_bulk_modulus<'b>(
        &self,
        source: &'b ResolvedCavitationState,
        state: &CavitatingFilmState,
        q_next: &[f64],
        channel_entrainment_m3_s: &[f64],
        dt_s: f64,
        bulk_modulus_pa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCompressibleCavitationContinuation<'b>, MaterialCavitationFilmError> {
        if cancelled() {
            return Err(FilmError("material compressible cavitation step cancelled").into());
        }
        require_temperature_continuation(self.source, source)?;
        require_reference_pressure(
            source,
            &self.pressure_axis,
            self.film.liquid().ambient_pressure_pa,
        )?;
        let material = source.state();
        validate_vapor_density(material, self.maximum_vapor_density_ratio)?;
        let viscosities =
            vec![material.dynamic_viscosity_pa_s(); self.film.geometry().channel_count()];
        let next = self.film.advance_with_liquid_and_bulk_modulus(
            state,
            q_next,
            channel_entrainment_m3_s,
            dt_s,
            CavitationLiquid {
                density_kg_m3: material.liquid_density_kg_m3(),
                cavitation_pressure_pa: material.saturation_pressure_pa(),
                ..self.film.liquid()
            },
            &viscosities,
            bulk_modulus_pa,
            cancelled,
        )?;
        validate_rupture_vapor_density(material, &next.film, self.maximum_vapor_density_ratio)?;
        Ok(MaterialCompressibleCavitationContinuation {
            film: MaterialCompressibleCavitatingFilm {
                source,
                material,
                pressure_axis: self.pressure_axis.clone(),
                film: next.film,
                maximum_vapor_density_ratio: self.maximum_vapor_density_ratio,
            },
            step: next.step,
        })
    }
}

fn validate_rupture_vapor_density(
    material: CavitationState,
    film: &CompressibleCavitatingFilm,
    maximum_ratio: f64,
) -> Result<(), FilmError> {
    if material.saturated_vapor().is_some_and(|vapor| {
        vapor.density_kg_m3() / film.density_at_cavitation_kg_m3() > maximum_ratio
    }) {
        return Err(FilmError(
            "saturated vapor density exceeds the compressible film's negligible-void-mass domain at rupture",
        ));
    }
    Ok(())
}

fn require_reference_pressure(
    source: &ResolvedCavitationState,
    axis: &str,
    reference_pressure_pa: f64,
) -> Result<(), FilmError> {
    let expected =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let material = source.material();
    if axis.is_empty()
        || axis == source.temperature_axis()
        || material.axis_quantities().get(axis) != Some(&expected)
        || !material.query_point().iter().any(|(name, pressure)| {
            name == axis
                && pressure.is_finite()
                && *pressure > 0.0
                && *pressure == reference_pressure_pa
        })
    {
        return Err(FilmError(
            "compressible film source needs an explicit static pressure axis at the native reference pressure",
        ));
    }
    Ok(())
}
