//! Source-pinned optical feedback from the original accepted cell temperatures.
//!
//! Every cell refreshes the original optical/permeability snapshot at its own
//! thermal temperature. The original coherent film solver then integrates its
//! absorbed work, and the original implicit heat solver advances that work once.
//! Optical coefficients are held during each interval. This is an explicit
//! quasi-steady thermo-optic split; time-step convergence remains caller-owned.
//! It supplies no implicit optical tangent, changing-geometry work, transient
//! electromagnetic storage, thermal emission or phase-change law.

use super::*;
use fs_em::plane_wave::{
    OpticalConductorResponse, OpticalMaterialResponse, SpectralOpticalResponse,
};
use fs_matdb::MaterialCard;
use fs_material::state_point::{
    MaterialStatePointError, SpectralOpticalStatePoint, VisibleConductorStatePoint,
    VisibleOpticalStatePoint,
};

/// Mutually incoherent spectral lines sharing one spatial heat update.
pub mod broadband;
pub use broadband::{
    BroadbandThermoOpticalLimits, BroadbandThermoOpticalWaveHeating,
    BroadbandThermoOpticalWaveStep, ThermoOpticalSpectralLine,
};

mod sealed {
    pub trait Sealed {}
    impl Sealed for fs_material::state_point::VisibleConductorStatePoint {}
    impl Sealed for fs_material::state_point::VisibleOpticalStatePoint {}
    impl Sealed for super::SpectralThermoOpticalSource {}
}

/// Original source families supported by this thermal/Maxwell bridge.
/// Sealed implementations retain the complete source snapshot on refresh;
/// arbitrary callbacks cannot replace source or constitutive admission.
/// The default conductor type preserves existing callers. A
/// `VisibleOpticalStatePoint` additionally supports mixed metal/glass layers;
/// `SpectralThermoOpticalSource` retains explicit finite-band infrared sources.
pub trait ThermoOpticalSource: sealed::Sealed + Clone {
    fn resolved(&self) -> &ResolvedMaterialStatePoint;
    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError>;
    /// Admit the requested frequency as well as the source's constitutive law.
    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError>;
}

impl ThermoOpticalSource for VisibleConductorStatePoint {
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.resolved()
    }
    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        self.at_query_point(card, point)
    }
    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError> {
        let response = OpticalConductorResponse::try_new(self)?;
        response.permittivity_f_m(omega_rad_s)?;
        Ok(DielectricResponse::OpticalConductor(response))
    }
}

impl ThermoOpticalSource for VisibleOpticalStatePoint {
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.resolved()
    }
    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        self.at_query_point(card, point)
    }
    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError> {
        let response = OpticalMaterialResponse::try_new(self)?;
        response.permittivity_f_m(omega_rad_s)?;
        Ok(DielectricResponse::OpticalMaterial(response))
    }
}

/// Original finite-band optics and their retained interpolation admission.
///
/// The source includes same-state magnetic permeability. Its n/k projection or
/// original n/alpha table is consumed directly by the existing Maxwell owner.
/// No source spectrum, nonmagnetic default or additional heat law is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralThermoOpticalSource {
    source: SpectralOpticalStatePoint,
    maximum_knot_span_m: f64,
}

impl SpectralThermoOpticalSource {
    /// Admit original source support, interpolation gaps and permeability.
    /// The maximum gap is a positive finite vacuum length in SI metres;
    /// frequency admission remains the original wave owner's responsibility.
    pub fn try_new(
        source: SpectralOpticalStatePoint,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        SpectralOpticalResponse::try_new((&source).into(), maximum_knot_span_m)?;
        Ok(Self {
            source,
            maximum_knot_span_m,
        })
    }

    /// Original source family, complete material snapshot and selected claims.
    #[must_use]
    pub const fn source(&self) -> &SpectralOpticalStatePoint {
        &self.source
    }

    /// Retained wavelength interpolation limit [vacuum m].
    #[must_use]
    pub const fn maximum_knot_span_m(&self) -> f64 {
        self.maximum_knot_span_m
    }
}

impl ThermoOpticalSource for SpectralThermoOpticalSource {
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.source.resolved()
    }

    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Ok(Self {
            source: self.source.at_query_point(card, point)?,
            maximum_knot_span_m: self.maximum_knot_span_m,
        })
    }

    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError> {
        // The thermal owner calls this after each cell refresh, before applying
        // work and before publishing endpoints. Recheck gap/mu and the actual
        // frequency together so a changed source band cannot silently extend.
        let response =
            SpectralOpticalResponse::try_new((&self.source).into(), self.maximum_knot_span_m)?;
        response.permittivity_f_m(omega_rad_s)?;
        Ok(DielectricResponse::SpectralOptical(response))
    }
}

/// One original optical model and the immutable card that owns its selected
/// claims. The complete snapshot includes the same-state magnetic permeability.
#[derive(Debug)]
pub struct ThermoOpticalLayer<'m, S: ThermoOpticalSource = VisibleConductorStatePoint> {
    pub source: &'m S,
    pub card: &'m MaterialCard,
}
impl<S: ThermoOpticalSource> Copy for ThermoOpticalLayer<'_, S> {}
impl<S: ThermoOpticalSource> Clone for ThermoOpticalLayer<'_, S> {
    fn clone(&self) -> Self {
        *self
    }
}

fn admit_sources<S: ThermoOpticalSource>(
    thermal: &BoundLayeredThermalBody<'_>,
    sources: &[ThermoOpticalLayer<'_, S>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LayeredWaveHeatingError> {
    poll(cancelled)?;
    let body = thermal.prepared_body();
    let layers = body.conduction().layers();
    if body.cell_count() > MAX_LAYERS || sources.len() != layers.len() {
        return Err(LayeredWaveHeatingError::Invalid(
            "thermo-optical heating needs one original source per layer and at most 64 cells",
        ));
    }
    for (source, layer) in sources.iter().zip(layers) {
        poll(cancelled)?;
        let material = source.source.resolved();
        if source.card.content_hash() != material.card_identity()
            || material
                .state_coordinate(layer.temperature_axis())
                .is_none()
        {
            return Err(LayeredWaveHeatingError::Invalid(
                "optical refresh requires its original card and the original thermal temperature coordinate",
            ));
        }
    }
    Ok(())
}

fn refresh_cell_sources<S: ThermoOpticalSource>(
    thermal: &BoundLayeredThermalBody<'_>,
    sources: &[ThermoOpticalLayer<'_, S>],
    state: &DepositionSimulation<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<S>, LayeredWaveHeatingError> {
    let temperatures = thermal.cell_temperatures(state, || checkpoint(cancelled))?;
    let body = thermal.prepared_body();
    let mut refreshed = Vec::with_capacity(temperatures.len());
    for (i, (source, layer)) in sources.iter().zip(body.conduction().layers()).enumerate() {
        for _ in 0..body.layer_cell_counts()[i] {
            poll(cancelled)?;
            let query = source
                .source
                .resolved()
                .query()
                .with(layer.temperature_axis(), temperatures[refreshed.len()])
                .map_err(|_| {
                    LayeredWaveHeatingError::Invalid(
                        "native cell temperature violates the original optical coordinate",
                    )
                })?;
            // Pin every original selected claim and all remaining coordinates,
            // including same-state mu and source-backed wavelength properties.
            refreshed.push(source.source.at_query_point(source.card, &query)?);
        }
    }
    poll(cancelled)?;
    Ok(refreshed)
}

/// Fixed physical geometry, spectral drive and original source selections.
/// Cell temperatures always come from the accepted native thermal state; a
/// caller cannot supply a separate optical temperature or averaged layer state.
pub struct ThermoOpticalWaveHeating<'b, 'a, 'm, S: ThermoOpticalSource = VisibleConductorStatePoint>
{
    thermal: &'b BoundLayeredThermalBody<'a>,
    incident: &'m ResolvedMaterialStatePoint,
    exit: DielectricResponse<'m>,
    omega_rad_s: f64,
    incidence_cosine: f64,
    polarization: Polarization,
    sources: Vec<ThermoOpticalLayer<'m, S>>,
}

/// One jointly admitted thermal endpoint and its optical source states.
/// Both state vectors are in the original thermal CellId/thickness order.
/// Applied states generated this interval's heat; endpoint states describe the
/// returned temperature and are not a second heat credit or an implicit solve.
pub struct ThermoOpticalWaveStep<'t, S: ThermoOpticalSource = VisibleConductorStatePoint> {
    pub step: LayeredWaveStep<'t>,
    pub applied_optical_states: Vec<S>,
    pub endpoint_optical_states: Vec<S>,
}

impl<'b, 'a, 'm, S: ThermoOpticalSource> ThermoOpticalWaveHeating<'b, 'a, 'm, S> {
    /// Bind one original admitted optical source per original thermal layer.
    /// Source cards, shared fixed conditions and temperature-axis conventions
    /// must agree with the existing conductivity/capacity owners. Exterior
    /// media, frequency, incidence and polarization remain held between steps.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: &'m ResolvedMaterialStatePoint,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
        sources: &[ThermoOpticalLayer<'m, S>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        admit_sources(thermal, sources, &mut cancelled)?;
        let mut responses = Vec::with_capacity(sources.len());
        for source in sources {
            poll(&mut cancelled)?;
            responses.push(source.source.wave_response(omega_rad_s)?);
        }
        // Reuse the original geometry/source and propagating-port admission.
        let original = LayeredWaveHeating::try_new(
            thermal,
            incident,
            exit,
            omega_rad_s,
            &responses,
            &mut cancelled,
        )?;
        original.prepare(incidence_cosine, polarization, &mut cancelled)?;
        Ok(Self {
            thermal,
            incident,
            exit,
            omega_rad_s,
            incidence_cosine,
            polarization,
            sources: sources.to_vec(),
        })
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }

    pub fn sources(&self) -> &[ThermoOpticalLayer<'m, S>] {
        &self.sources
    }

    fn refresh_sources(
        &self,
        state: &DepositionSimulation<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<S>, LayeredWaveHeatingError> {
        let refreshed = refresh_cell_sources(self.thermal, &self.sources, state, cancelled)?;
        for source in &refreshed {
            poll(cancelled)?;
            source.wave_response(self.omega_rad_s)?;
        }
        poll(cancelled)?;
        Ok(refreshed)
    }

    /// Refresh the optics at the actual accepted temperatures, then advance
    /// the original heat owner with those held coherent cell powers. Finally
    /// refresh all endpoint optical states BEFORE publishing the joint result.
    /// An optical domain exit, thermal refusal or cancellation leaves initial
    /// unchanged; no partial optical or thermal update is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        drive: LayeredWaveDrive<'_>,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermoOpticalWaveStep<'t, S>, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let applied_optical_states = self.refresh_sources(initial, &mut cancelled)?;
        let step = {
            let responses = applied_optical_states
                .iter()
                .map(|state| state.wave_response(self.omega_rad_s))
                .collect::<Result<Vec<_>, _>>()?;
            let owner = LayeredWaveHeating::try_new_cellwise(
                self.thermal,
                self.incident,
                self.exit,
                self.omega_rad_s,
                &responses,
                &mut cancelled,
            )?;
            let prepared =
                owner.prepare(self.incidence_cosine, self.polarization, &mut cancelled)?;
            prepared.advance_to(
                initial,
                target_s,
                drive,
                conduction,
                controls,
                workspace,
                &mut cancelled,
            )?
        };
        let endpoint_optical_states = self.refresh_sources(&step.thermal, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ThermoOpticalWaveStep {
            step,
            applied_optical_states,
            endpoint_optical_states,
        })
    }
}
