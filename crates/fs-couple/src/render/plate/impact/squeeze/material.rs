//! Source-retaining gas/liquid transport binding for geometric film resistance.
use super::super::{ImpactError, ImpactSystem};
use super::{FilmCell, FilmChannel, FilmLimits, ResistiveFilm};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

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
    pub fn new(
        source: &'a F,
        cells: Vec<FilmCell>,
        channels: Vec<FilmChannel>,
        ports: usize,
        limits: FilmLimits,
    ) -> Result<Self, ImpactError> {
        let transport = source
            .transport_state()
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let film = ResistiveFilm::new(
            cells,
            channels,
            ports,
            transport.dynamic_viscosity_pa_s(),
            limits,
        )
        .map_err(|e| super::super::invalid(e.0))?;
        Ok(Self {
            source,
            transport,
            film,
        })
    }
    /// Exact immutable source, including its concrete material receipts.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Shared snapshot; density, caloric and conductivity data remain available
    /// to the caller's separate thermal or fluid-inertia owner.
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// The original pressure law and analytic tangent, without a second solver.
    pub const fn film(&self) -> &ResistiveFilm {
        &self.film
    }
    /// Attach to the ordinary joint mechanical equation. This source binding
    /// remains available for inspection while the system owns its immutable
    /// compiled operator. Attach before preparation; the returned insertion
    /// index can select the same island in `SqueezeFilmHeatRoute`.
    pub fn attach(&self, system: ImpactSystem) -> Result<(ImpactSystem, usize), ImpactError> {
        let index = system.squeeze_film_count();
        let system = system.with_squeeze_film(self.film.clone())?;
        Ok((system, index))
    }
}

/// One retained material evaluation per hydraulic passage of a pressure island.
///
/// Each passage is Newtonian and uniform across its thickness and length; its
/// viscosity may differ from adjacent passages. This represents a resolved
/// piecewise-constant thermal/material field, not a homogenized island viscosity
/// or a model of mixing, species transport, or thermal evolution.
pub struct MaterialChannelSqueezeFilm<'a, F: FluidTransport> {
    sources: &'a [F],
    transport: Vec<FluidTransportState>,
    film: ResistiveFilm,
}

impl<'a, F: FluidTransport> MaterialChannelSqueezeFilm<'a, F> {
    /// Resolve sources in the exact supplied channel order and compile their
    /// dynamic viscosities into the original Reynolds pressure operator.
    ///
    /// A slice of `&dyn FluidTransport` admits heterogeneous source owners while
    /// retaining their complete receipts. Every channel needs exactly one valid
    /// source, including a channel whose current aperture is closed. No source
    /// may be replaced by a nominal value. The geometry and count limits are
    /// admitted before querying sources beyond the first.
    ///
    /// The assumptions and pressure bounds of [`MaterialSqueezeFilm::new`]
    /// apply to every passage. All snapshots are frozen during the mechanical
    /// sample; this binding does not interpolate a viscosity between endpoints.
    pub fn new(
        sources: &'a [F],
        cells: Vec<FilmCell>,
        channels: Vec<FilmChannel>,
        ports: usize,
        limits: FilmLimits,
    ) -> Result<Self, ImpactError> {
        if sources.is_empty() || sources.len() != channels.len() {
            return Err(super::super::invalid(
                "channel film needs exactly one material source per passage",
            ));
        }
        let first = sources[0]
            .transport_state()
            .map_err(|error| ImpactError::Owner(error.to_string()))?;
        let geometry = ResistiveFilm::new(
            cells,
            channels,
            ports,
            first.dynamic_viscosity_pa_s(),
            limits,
        )
        .map_err(|error| super::super::invalid(error.0))?;
        let mut transport = Vec::with_capacity(sources.len());
        transport.push(first);
        for source in &sources[1..] {
            transport.push(
                source
                    .transport_state()
                    .map_err(|error| ImpactError::Owner(error.to_string()))?,
            );
        }
        let viscosities: Vec<_> = transport
            .iter()
            .map(|state| state.dynamic_viscosity_pa_s())
            .collect();
        let film = geometry
            .with_channel_viscosities(&viscosities)
            .map_err(|error| super::super::invalid(error.0))?;
        Ok(Self {
            sources,
            transport,
            film,
        })
    }

    /// Original material owners, in the hydraulic channel order.
    pub const fn sources(&self) -> &'a [F] {
        self.sources
    }

    /// Checked snapshots used to compile each channel, in the same order.
    pub fn transport(&self) -> &[FluidTransportState] {
        &self.transport
    }

    /// The original pressure/force law with its exact local viscosity field.
    pub const fn film(&self) -> &ResistiveFilm {
        &self.film
    }

    /// Attach the heterogeneous pressure island to the ordinary joint mechanical
    /// equation. Pressure, analytic tangent and viscous work all use this same
    /// compiled field; the existing film heat routes can consume its work.
    pub fn attach(&self, system: ImpactSystem) -> Result<(ImpactSystem, usize), ImpactError> {
        let index = system.squeeze_film_count();
        let system = system.with_squeeze_film(self.film.clone())?;
        Ok((system, index))
    }
}
