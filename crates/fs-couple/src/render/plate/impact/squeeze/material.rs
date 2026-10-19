//! Source-retaining gas/liquid transport binding for geometric film resistance.
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use super::{FilmCell, FilmChannel, FilmLimits, ResistiveFilm};
use super::super::{ImpactError, ImpactSystem};

/// One material snapshot and the original Reynolds pressure operator it drives.
/// Geometry and source are independent data. The original source retains its
/// complete card/mixture/query receipts; no nominal fluid or duplicate viscosity
/// is accepted by this binding.
pub struct MaterialSqueezeFilm<'a, F: FluidTransport + ?Sized> {
    source: &'a F,
    transport: FluidTransportState,
    film: ResistiveFilm,
}
impl<'a, F: FluidTransport + ?Sized> MaterialSqueezeFilm<'a, F> {
    /// Resolve the same checked transport view used by flow and heat owners,
    /// then construct the existing geometric pressure graph. Gas snapshots must
    /// retain coherent EOS/caloric/acoustic fields even though this quasistatic
    /// pressure image consumes only their dynamic viscosity.
    ///
    /// The caller declares an incompressible Newtonian, no-slip, thin-gap
    /// regime. The existing pressure/gap bounds apply inside every solve. Gas
    /// pressure excursions must be small relative to ambient; liquid absolute
    /// pressure must stay above cavitation onset. This route does not infer a
    /// vapor-pressure law, rarefied slip, thermal expansion or added fluid mass.
    pub fn new(source: &'a F, cells: Vec<FilmCell>, channels: Vec<FilmChannel>,
        ports: usize, limits: FilmLimits) -> Result<Self, ImpactError> {
        let transport=source.transport_state().map_err(|e|ImpactError::Owner(e.to_string()))?;
        let film=ResistiveFilm::new(cells,channels,ports,transport.dynamic_viscosity_pa_s(),limits)
            .map_err(|e|super::super::invalid(e.0))?;
        Ok(Self {source,transport,film})
    }
    /// Exact immutable source, including its concrete material receipts.
    pub const fn source(&self) -> &'a F { self.source }
    /// Shared snapshot; density, caloric and conductivity data remain available
    /// to the caller's separate thermal or fluid-inertia owner.
    pub const fn transport(&self) -> FluidTransportState { self.transport }
    /// The original pressure law and analytic tangent, without a second solver.
    pub const fn film(&self) -> &ResistiveFilm { &self.film }
    /// Attach to the ordinary joint mechanical equation. This source binding
    /// remains available for inspection while the system owns its immutable
    /// compiled operator. Attach before preparation; the returned insertion
    /// index can select the same island in `SqueezeFilmHeatRoute`.
    pub fn attach(&self, system: ImpactSystem) -> Result<(ImpactSystem, usize), ImpactError> {
        let index=system.squeeze_film_count();
        let system=system.with_squeeze_film(self.film.clone())?;
        Ok((system,index))
    }
}
