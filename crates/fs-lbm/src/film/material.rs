//! Source-resolved liquid properties for the native Newtonian film model.
//!
//! The same material owner already supplies momentum and capillarity to the
//! free-surface LBM. This adapter borrows that complete state, retaining its
//! original viscosity basis, query, property receipts and interface context.
//! It requests no caloric properties and performs no new material selection.

use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

use super::transient::{FilmDrive, FilmStep, FilmStepLimits, TransientFilm, TransientFilmError};
use super::{FilmError, FilmMaterial, NusseltFilm, ShearedFilm};

#[cfg(feature = "thin-film-evolution")]
use super::evolution::adaptive::{
    ThinFilmAdaptiveControls, ThinFilmAdaptiveError, ThinFilmAdaptiveReport, ThinFilmDrive,
    ThinFilmDriveQuery,
};
#[cfg(feature = "thin-film-evolution")]
use super::evolution::{
    ThinFilmControls, ThinFilmError, ThinFilmGeometry, ThinFilmLimits, ThinFilmLine,
    ThinFilmOpenEnds, ThinFilmStep,
};

/// Evolving longitudinal thickness retaining the exact admitted Newtonian
/// liquid source, viscosity basis, query and capillarity declaration.
#[cfg(feature = "thin-film-evolution")]
#[derive(Debug)]
pub struct MaterialThinFilmLine<'a> {
    liquid: MaterialFilmLiquid<'a>,
    flow: ThinFilmLine,
}

#[cfg(feature = "thin-film-evolution")]
impl<'a> MaterialThinFilmLine<'a> {
    /// Original immutable source, including all material/property receipts.
    pub const fn liquid(&self) -> MaterialFilmLiquid<'a> {
        self.liquid
    }

    /// Accepted native height inventory, physical clock and film geometry.
    pub const fn flow(&self) -> &ThinFilmLine {
        &self.flow
    }

    /// Evolve the original thickness inventory under the same source density,
    /// viscosity and surface tension. Native volume, support, work and
    /// cancellation gates decide acceptance. Temperature/composition remain
    /// frozen; no caloric law or rheological approximation is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        target_s: f64,
        acceleration_m_s2: f64,
        surface_shear_pa: &[f64],
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        height_guess_m: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThinFilmStep, ThinFilmError> {
        self.flow.advance(
            target_s,
            acceleration_m_s2,
            surface_shear_pa,
            limits,
            controls,
            height_guess_m,
            cancelled,
        )
    }

    /// Advance both native pressure-driven open connections with this same
    /// original liquid source. Native half-cell loss, reservoir work and signed
    /// volume transport are admitted together; all source receipts are retained.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_open(
        &mut self,
        target_s: f64,
        acceleration_m_s2: f64,
        interior_surface_shear_pa: &[f64],
        ends: ThinFilmOpenEnds,
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        height_guess_m: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThinFilmStep, ThinFilmError> {
        self.flow.advance_open(
            target_s,
            acceleration_m_s2,
            interior_surface_shear_pa,
            ends,
            limits,
            controls,
            height_guess_m,
            cancelled,
        )
    }
    /// Adapt the same source-bound Newtonian film to an authored load schedule.
    /// The original liquid query/receipts remain frozen. Drive temperatures
    /// must be None; a new source state cannot replace this borrowed owner.
    /// Native fine-step accounts, source refusals and whole-call atomicity are
    /// inherited without a separate transport or viscosity implementation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive<E>(
        &mut self,
        target_s: f64,
        event_times_s: &[f64],
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        adaptive: ThinFilmAdaptiveControls,
        resolve: impl FnMut(ThinFilmDriveQuery) -> Result<ThinFilmDrive, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThinFilmAdaptiveReport, ThinFilmAdaptiveError<E>> {
        self.flow.advance_adaptive(
            target_s,
            event_times_s,
            limits,
            controls,
            adaptive,
            resolve,
            cancelled,
        )
    }
}

/// Resolved film inertia retaining the original immutable liquid source.
#[derive(Debug)]
pub struct MaterialTransientFilm<'a> {
    liquid: MaterialFilmLiquid<'a>,
    flow: TransientFilm,
}

impl<'a> MaterialTransientFilm<'a> {
    /// Original source, viscosity basis, query and capillarity declaration.
    pub const fn liquid(&self) -> MaterialFilmLiquid<'a> {
        self.liquid
    }

    /// Accepted native profile, physical clock and original SI coefficients.
    pub const fn flow(&self) -> &TransientFilm {
        &self.flow
    }

    /// Prescribe the original native Navier length without altering the liquid
    /// source receipt. The length is authored boundary data, not inferred from
    /// the material's viscosity, surface tension or contact angle.
    pub fn set_wall_slip_length_m(&mut self, length_m: f64) -> Result<(), TransientFilmError> {
        self.flow.set_wall_slip_length_m(length_m)
    }

    /// Advance the same source-bound velocity history without another material
    /// selection, equilibrium reset, heat deposition or source query.
    pub fn step(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        self.flow.step(duration_s, drive, limits, cancelled)
    }

    /// Export the original interval's local viscous losses for an explicitly
    /// selected thermal owner, retaining this same source-bound momentum.
    /// Cell heat is J/m² of wetted wall area; no heat capacity, source lookup or
    /// thermal inventory is inferred. Flow and output publish together.
    pub fn step_with_heat(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        cell_heat_j_m2: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        self.flow
            .step_with_heat(duration_s, drive, limits, cell_heat_j_m2, cancelled)
    }
}

/// A frozen Newtonian liquid binding for repeated local film evaluations.
///
/// The material owner admits its source conditions. This adapter additionally
/// rejects states requiring evolving shear/pressure rheology, because the
/// Nusselt/Benney model has one constant Newtonian viscosity. Temperature or
/// composition changes require a newly admitted source state; this binding
/// supplies no heat evolution, phase-stability or source-accuracy certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFilmLiquid<'a> {
    material: &'a ResolvedLiquidSurfaceState,
}

impl<'a> MaterialFilmLiquid<'a> {
    /// Bind the exact state used by the existing material free-surface owner.
    /// No properties, source receipts or source coordinates are cloned.
    ///
    /// # Errors
    /// Refuses a state that requires an evolving constitutive owner instead of
    /// silently treating its sampled apparent viscosity as Newtonian.
    pub fn try_new(material: &'a ResolvedLiquidSurfaceState) -> Result<Self, FilmError> {
        if material.rheology_axes().is_some() || material.rheology_pressure_axis().is_some() {
            return Err(FilmError::InvalidInput(
                "Newtonian film cannot freeze an evolving rheology source",
            ));
        }
        Ok(Self { material })
    }

    /// Original source owner, including every retained receipt and query value.
    #[must_use]
    pub const fn material(self) -> &'a ResolvedLiquidSurfaceState {
        self.material
    }

    /// SI coefficients for the existing film and vessel evaluators.
    ///
    /// `None` in the admitted source means capillarity was explicitly neglected;
    /// only that declaration becomes zero. Missing active surface data already
    /// refuses in the material owner. Retain this binding beside numerical
    /// results when their source identity is required.
    #[must_use]
    pub const fn coefficients(self) -> FilmMaterial {
        FilmMaterial {
            density_kg_m3: self.material.density_kg_m3(),
            kinematic_viscosity_m2_s: self.material.kinematic_viscosity_m2_s(),
            surface_tension_n_m: match self.material.surface_tension_n_m() {
                Some(value) => value,
                None => 0.0,
            },
        }
    }

    /// Evaluate the native Nusselt base flow with these SAME source values.
    /// All physical inputs and numerical refusals belong to [`NusseltFilm`].
    /// This fixed-cost call allocates nothing and performs no source queries.
    pub fn film(
        self,
        flow_per_width_m2_s: f64,
        tangential_acceleration_m_s2: f64,
        normal_acceleration_m_s2: f64,
    ) -> Result<NusseltFilm, FilmError> {
        NusseltFilm::new(
            self.coefficients(),
            flow_per_width_m2_s,
            tangential_acceleration_m_s2,
            normal_acceleration_m_s2,
        )
    }

    /// Drive a film of prescribed thickness with signed surface traction,
    /// retaining this same frozen source viscosity and density. Traction must
    /// come from the caller's gas or interface law; no new source is inferred.
    pub fn sheared_film(
        self,
        thickness_m: f64,
        tangential_acceleration_m_s2: f64,
        surface_shear_pa: f64,
    ) -> Result<ShearedFilm, FilmError> {
        ShearedFilm::new(
            self.coefficients(),
            thickness_m,
            tangential_acceleration_m_s2,
            surface_shear_pa,
        )
    }

    /// Prepare a fixed-thickness Newtonian film initially at rest. Subsequent
    /// drive changes evolve actual wall-normal momentum under this source's
    /// viscosity and density. Source temperature/composition remain frozen.
    pub fn transient_film(
        self,
        thickness_m: f64,
        cells: usize,
    ) -> Result<MaterialTransientFilm<'a>, TransientFilmError> {
        Ok(MaterialTransientFilm {
            liquid: self,
            flow: TransientFilm::new(self.coefficients(), thickness_m, cells)?,
        })
    }

    /// Prepare native longitudinal thickness evolution from this exact source.
    /// The Newtonian constructor receives the original SI coefficients directly;
    /// no surrogate Carreau parameters, source selection or heat capacity enters.
    #[cfg(feature = "thin-film-evolution")]
    pub fn thin_film_line(
        self,
        geometry: ThinFilmGeometry,
        heights_m: Vec<f64>,
        time_s: f64,
    ) -> Result<MaterialThinFilmLine<'a>, ThinFilmError> {
        Ok(MaterialThinFilmLine {
            liquid: self,
            flow: ThinFilmLine::newtonian(self.coefficients(), geometry, heights_m, time_s)?,
        })
    }
}
