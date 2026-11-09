//! Quasi-steady temperature feedback for coherent, polychromatic film heating.
//!
//! One original homogeneous film shares a uniform, fully solid thermal body.
//! At each explicit feedback interval its complete original material bundle,
//! including permeability, is refreshed at the accepted temperature. The native
//! film, two-face angular-spectrum, mixed-statistics and volume-heat owners then
//! rebuild the applied watt spectrum. Native enthalpy evolution consumes that
//! spectrum once, with its absolute phase clock retained in double-double.
//! Original Debye relaxation, Lorentz resonance and independent Ohmic leakage
//! keep their separate reservoir laws; each refresh retains every selected
//! source role, including Lorentz frequency/damping and same-state permeability.
//!
//! Coefficients are HELD over the requested feedback interval, even when the
//! thermal solver takes smaller steps. The temperature-excursion bound controls
//! this explicit splitting; it is not a coupled residual or an error certificate.
//! Every native thermal sample and the endpoint material law must be admitted
//! before a result is returned. Original inputs are never mutated.
//!
//! This is a quasi-steady constitutive approximation, not a transient Maxwell
//! update. No polarization history, electromagnetic storage change, parameter
//! work or switch-on transient is transferred across a coefficient refresh.
//! Such terms are not deposited as heat. The native Biot gate is necessary for
//! this reduced thermal body but does not establish internal equilibration under
//! a nonuniform standing wave or beam: that separation of time scales remains
//! an explicit modelling assumption. Geometry, density, phase, other films and
//! exterior media are held. A nonlinear trajectory driven by expected heat is
//! not, in general, an ensemble-mean temperature.

use fs_conduction::{
    ConductionError,
    lumped::{
        BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig,
        interference::enthalpy::{
            LumpedEnthalpyPowerBinding, LumpedEnthalpySpectrumControls,
            LumpedEnthalpySpectrumEnergy, LumpedEnthalpySpectrumSample,
            solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable,
        },
    },
};
use fs_em::plane_wave::{
    DielectricResponse, ReflectionError,
    angular_spectrum::{
        AngularSpectrumError, AngularSpectrumExposure, AngularSpectrumExposureEnsemble,
        AngularSpectrumExposureEnsembleLimits, AngularSpectrumExposureHeatChannel,
        AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumFactorStatistics, AngularSpectrumGrid, AngularSpectrumHeatAreaQuadrature,
        AngularSpectrumHeatVolumeQuadrature, AngularSpectrumIlluminationLimits,
        AngularSpectrumIntegratedHeatLimits, AngularSpectrumLimits, PlanarAngularSpectrum,
        PreparedAngularSpectrumIllumination,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_fft::harmonic_window::PreparedHarmonicProductSpectrum;
use fs_matdb::{MaterialCard, PropertyUsageReceipt};
use fs_material::{
    conductor::{ConductorError, ResolvedConductorState},
    field_storage::{
        FieldStorageError,
        dielectric::{ResolvedDebyeCapacitor, lorentz::ResolvedLorentzMaterial},
    },
    phase::{EquilibriumPhaseState, PhaseStateError},
    state_point::{MaterialStatePointError, ResolvedMaterialStatePoint},
};
use fs_math::{c64::C64, dd::Dd};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

/// Peak tangential fields on the two original physical faces, in V/m.
/// Both faces use the original grid axes and the same exp(+i omega t) clock.
/// With an explicit transverse carrier, these are envelope samples with the
/// same exp(-i k_carrier.x) removed on both faces and at every frequency.
/// An explicitly zero right-face array describes left-only illumination;
/// the original propagating two-port admission still applies.
#[derive(Clone, Copy, Debug)]
pub struct ThermalInterferenceLine<'a> {
    pub angular_frequency_rad_s: f64,
    pub left_tangential_e_v_per_m: &'a [[C64; 2]],
    pub right_tangential_e_v_per_m: &'a [[C64; 2]],
}

/// One zero-mean source factor spanning all its frequencies and both faces.
/// Fields already contain its statistical weight. Different factors have
/// the original jointly uncorrelated real/proper second moments.
#[derive(Clone, Copy, Debug)]
pub struct ThermalInterferenceFactor<'a> {
    pub lines: &'a [ThermalInterferenceLine<'a>],
    pub statistics: AngularSpectrumFactorStatistics,
}

/// Authored deterministic mean and explicit statistical factors.
/// Empty banks are allowed; at least one authored line is required to anchor
/// the field preparation. Its amplitudes may all be zero. The clock label is
/// a caller declaration of common phase origin, not an inferred synchronization.
#[derive(Clone, Copy, Debug)]
pub struct ThermalInterferenceDrive<'a> {
    pub source_clock: &'a str,
    pub mean: &'a [ThermalInterferenceLine<'a>],
    pub factors: &'a [ThermalInterferenceFactor<'a>],
}

/// Original owner limits plus bounds on their simultaneous preparation.
/// Native workspace figures are conservatively summed while their owners
/// coexist. Borrowed input arrays, material cards, resolved material-source
/// payloads and allocator overhead are excluded. Native thermal state, output
/// and solve storage remain separately bounded by their receiving owner.
#[derive(Clone, Copy, Debug)]
pub struct ThermalInterferenceLimits {
    pub maximum_lines: usize,
    pub maximum_preparation_workspace_bytes: usize,
    pub angular: AngularSpectrumLimits,
    pub illumination: AngularSpectrumIlluminationLimits,
    pub exposure: AngularSpectrumExposureLimits,
    pub ensemble: AngularSpectrumExposureEnsembleLimits,
    pub layer: AngularSpectrumExposureLayerLimits,
    pub heat: AngularSpectrumIntegratedHeatLimits,
}

#[derive(Debug)]
pub enum ThermalInterferenceError {
    Invalid(&'static str),
    Electromagnetic(AngularSpectrumError),
    Film(ReflectionError),
    Material(MaterialStatePointError),
    Dielectric(FieldStorageError),
    Conductor(ConductorError),
    Phase(PhaseStateError),
    Thermal(ConductionError),
    Cancelled,
}
impl core::fmt::Display for ThermalInterferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quasi-steady interference heating: {self:?}")
    }
}
impl std::error::Error for ThermalInterferenceError {}
impl From<AngularSpectrumError> for ThermalInterferenceError {
    fn from(value: AngularSpectrumError) -> Self {
        Self::Electromagnetic(value)
    }
}
impl From<ReflectionError> for ThermalInterferenceError {
    fn from(value: ReflectionError) -> Self {
        Self::Film(value)
    }
}
impl From<MaterialStatePointError> for ThermalInterferenceError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl From<FieldStorageError> for ThermalInterferenceError {
    fn from(value: FieldStorageError) -> Self {
        Self::Dielectric(value)
    }
}
impl From<ConductorError> for ThermalInterferenceError {
    fn from(value: ConductorError) -> Self {
        Self::Conductor(value)
    }
}
impl From<PhaseStateError> for ThermalInterferenceError {
    fn from(value: PhaseStateError) -> Self {
        Self::Phase(value)
    }
}
impl From<ConductionError> for ThermalInterferenceError {
    fn from(value: ConductionError) -> Self {
        Self::Thermal(value)
    }
}

/// Immutable accepted state on this original body and declared source clock.
#[derive(Clone, Copy)]
pub struct QuasiSteadyFilmState<'a> {
    body: &'a LumpedEnthalpyBody<'a>,
    source_clock: &'a str,
    source_time_s: Dd,
    phase: EquilibriumPhaseState,
}
impl<'a> QuasiSteadyFilmState<'a> {
    pub fn phase_state(self) -> EquilibriumPhaseState {
        self.phase
    }
    pub fn source_time_s_dd(self) -> Dd {
        self.source_time_s
    }
    pub fn source_clock(self) -> &'a str {
        self.source_clock
    }
}

pub(super) enum RefreshedFilmLaw {
    Debye(ResolvedDebyeCapacitor),
    Lorentz(ResolvedLorentzMaterial),
    Ohmic(ResolvedConductorState),
    DebyeWithOhmic(ResolvedDebyeCapacitor, ResolvedConductorState),
    LorentzWithOhmic(ResolvedLorentzMaterial, ResolvedConductorState),
}
impl RefreshedFilmLaw {
    pub(super) fn response(&self) -> DielectricResponse<'_> {
        match self {
            Self::Debye(source) => DielectricResponse::Debye(source),
            Self::Lorentz(source) => DielectricResponse::Lorentz(source),
            Self::Ohmic(source) => DielectricResponse::Ohmic(source),
            Self::DebyeWithOhmic(dielectric, conductor) => DielectricResponse::DebyeWithOhmic {
                dielectric,
                conductor,
            },
            Self::LorentzWithOhmic(dielectric, conductor) => DielectricResponse::LorentzWithOhmic {
                dielectric,
                conductor,
            },
        }
    }
}

/// One changing film, the original immutable stack and its thermal receiver.
/// The full original periodic cell and complete film thickness are integrated.
/// The fixed-density, fully solid phase chart and body mass must describe that
/// same volume; partial masks and independent volume multipliers are refused.
pub struct QuasiSteadyFilmHeating<'a, 'm> {
    original: &'a PlanarDielectricLayers<'m>,
    layer_index: usize,
    card: &'a MaterialCard,
    temperature_axis: &'a str,
    body: &'a LumpedEnthalpyBody<'a>,
    grid: AngularSpectrumGrid,
    carrier_rad_m: [f64; 2],
    drive: ThermalInterferenceDrive<'a>,
    quadrature: AngularSpectrumHeatVolumeQuadrature<'a>,
    limits: ThermalInterferenceLimits,
    density_kg_m3: f64,
    line_count: usize,
}

/// A newly resolved, settled heat load. Its spectrum is already in W.
/// Original owner, exact material snapshot, quadrature and native thermal state
/// remain inspectable after the temporary borrowed Maxwell caches are dropped.
pub struct QuasiSteadyFilmHeatLoad<'p, 'a, 'm> {
    owner: &'p QuasiSteadyFilmHeating<'a, 'm>,
    state: QuasiSteadyFilmState<'a>,
    source: RefreshedFilmLaw,
    spectrum: PreparedHarmonicProductSpectrum,
    volume_m3: Dd,
    preparation_workspace_bytes: usize,
}
impl<'p, 'a, 'm> QuasiSteadyFilmHeatLoad<'p, 'a, 'm> {
    pub fn owner(&self) -> &'p QuasiSteadyFilmHeating<'a, 'm> {
        self.owner
    }
    pub fn state(&self) -> QuasiSteadyFilmState<'a> {
        self.state
    }
    pub fn response(&self) -> DielectricResponse<'_> {
        self.source.response()
    }
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    pub fn represented_volume_m3_dd(&self) -> Dd {
        self.volume_m3
    }
    pub fn preparation_workspace_bytes(&self) -> usize {
        self.preparation_workspace_bytes
    }
}

/// Jointly admitted explicit-feedback interval. Applied heat and the native
/// thermal ledger refer to the initial-temperature law held over this interval.
/// The endpoint response is an admission snapshot for subsequent rebuilding;
/// it supplies no second heat credit or reversible energy transfer.
pub struct QuasiSteadyFilmStep<'p, 'a, 'm> {
    applied_load: QuasiSteadyFilmHeatLoad<'p, 'a, 'm>,
    endpoint: QuasiSteadyFilmState<'a>,
    endpoint_source: RefreshedFilmLaw,
    samples: Vec<LumpedEnthalpySpectrumSample>,
    energy: LumpedEnthalpySpectrumEnergy,
    maximum_biot: f64,
    maximum_temperature_excursion_k: f64,
    attempts: usize,
    spectral_work: usize,
}
impl<'p, 'a, 'm> QuasiSteadyFilmStep<'p, 'a, 'm> {
    pub fn applied_load(&self) -> &QuasiSteadyFilmHeatLoad<'p, 'a, 'm> {
        &self.applied_load
    }
    pub fn endpoint(&self) -> QuasiSteadyFilmState<'a> {
        self.endpoint
    }
    pub fn endpoint_response(&self) -> DielectricResponse<'_> {
        self.endpoint_source.response()
    }
    pub fn samples(&self) -> &[LumpedEnthalpySpectrumSample] {
        &self.samples
    }
    pub fn energy(&self) -> LumpedEnthalpySpectrumEnergy {
        self.energy
    }
    pub fn maximum_biot(&self) -> f64 {
        self.maximum_biot
    }
    pub fn maximum_temperature_excursion_k(&self) -> f64 {
        self.maximum_temperature_excursion_k
    }
    pub fn attempts(&self) -> usize {
        self.attempts
    }
    pub fn spectral_work(&self) -> usize {
        self.spectral_work
    }
}

impl<'a, 'm> QuasiSteadyFilmHeating<'a, 'm> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        original: &'a PlanarDielectricLayers<'m>,
        layer_index: usize,
        card: &'a MaterialCard,
        temperature_axis: &'a str,
        body: &'a LumpedEnthalpyBody<'a>,
        grid: AngularSpectrumGrid,
        drive: ThermalInterferenceDrive<'a>,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'a>,
        limits: ThermalInterferenceLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalInterferenceError> {
        Self::try_new_with_carrier(
            original,
            layer_index,
            card,
            temperature_axis,
            body,
            grid,
            [0.0; 2],
            drive,
            quadrature,
            limits,
            cancelled,
        )
    }

    /// Retain a common physical transverse wavevector [rad/m] through every
    /// frequency solve and temperature rebuild. Source arrays describe its
    /// periodic envelopes, so an arbitrary incidence angle does not require
    /// the spatial grid to sample the carrier's short wavelength.
    ///
    /// This holds wavevector, not angle: each line's original exterior law
    /// and frequency determine its longitudinal wavevector and incidence.
    /// Both faces and every mean/factor line share this carrier, as required
    /// by the native coherent exposure. The full shifted spectrum must still
    /// propagate in both exterior ports when the heat load is prepared.
    /// `ContinuousPeriodicCell` integrates transverse carrier products exactly;
    /// sampled `PeriodicCell` heat still needs its own spatial resolution.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_carrier(
        original: &'a PlanarDielectricLayers<'m>,
        layer_index: usize,
        card: &'a MaterialCard,
        temperature_axis: &'a str,
        body: &'a LumpedEnthalpyBody<'a>,
        grid: AngularSpectrumGrid,
        carrier_rad_m: [f64; 2],
        drive: ThermalInterferenceDrive<'a>,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'a>,
        limits: ThermalInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalInterferenceError> {
        poll(&mut cancelled)?;
        if carrier_rad_m.iter().any(|k| !k.is_finite()) {
            return Err(bad("the transverse carrier must be finite in rad/m"));
        }
        let film = original
            .layer(layer_index)
            .ok_or(bad("target film is absent"))?;
        admit_refresh_source(film.response, card, temperature_axis)?;
        let material = film.response.material();
        if body.phase_curve().material_card_identity() != card.content_hash() {
            return Err(bad(
                "field, thermal curve and original immutable card must agree",
            ));
        }
        for receipt in body
            .phase_curve()
            .source_receipts()
            .iter()
            .chain(body.transport().conductivity_receipts())
        {
            poll(&mut cancelled)?;
            admit_thermal_context(material, temperature_axis, receipt)?;
        }
        let first = body
            .phase_curve()
            .knots()
            .first()
            .ok_or(bad("the original thermal chart is empty"))?;
        let density_kg_m3 = first.bulk_density_kg_m3;
        for knot in body.phase_curve().knots() {
            poll(&mut cancelled)?;
            if knot.liquid_mass_fraction != 0.0 || knot.bulk_density_kg_m3 != density_kg_m3 {
                return Err(bad(
                    "fixed film feedback requires a fully solid, fixed-density thermal chart",
                ));
            }
        }
        if film.thickness_m <= 0.0
            || quadrature.depth_span_m != [0.0, film.thickness_m]
            || !matches!(
                quadrature.transverse,
                AngularSpectrumHeatAreaQuadrature::PeriodicCell
                    | AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell
            )
        {
            return Err(bad(
                "feedback integrates the full film thickness and original periodic cell",
            ));
        }
        let line_count = admit_drive(drive, limits, &mut cancelled)?;
        let result = Self {
            original,
            layer_index,
            card,
            temperature_axis,
            body,
            grid,
            carrier_rad_m,
            drive,
            quadrature,
            limits,
            density_kg_m3,
            line_count,
        };
        // Match the physical geometry before material refresh or field allocation.
        let n = grid
            .width
            .checked_mul(grid.height)
            .ok_or(bad("grid size overflows"))?;
        if n == 0 || grid.spacing_m.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return Err(bad(
                "the original grid must have positive finite cell dimensions",
            ));
        }
        let volume = Dd::from_f64(grid.spacing_m[0])
            * Dd::from_f64(grid.spacing_m[1])
            * Dd::from_f64(n as f64)
            * Dd::from_f64(film.thickness_m);
        result.admit_volume(volume)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    pub fn original_stack(&self) -> &'a PlanarDielectricLayers<'m> {
        self.original
    }
    pub fn thermal_body(&self) -> &'a LumpedEnthalpyBody<'a> {
        self.body
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    /// Authored common wavevector, unchanged across frequencies and refreshes.
    pub const fn carrier_transverse_wavenumber_rad_m(&self) -> [f64; 2] {
        self.carrier_rad_m
    }
    pub fn drive(&self) -> ThermalInterferenceDrive<'a> {
        self.drive
    }
    pub fn quadrature(&self) -> AngularSpectrumHeatVolumeQuadrature<'a> {
        self.quadrature
    }

    /// Admit an initial state through the original thermal chart and the
    /// original complete electromagnetic source at that temperature.
    pub fn initial_state(
        &self,
        source_time_s: Dd,
        specific_enthalpy_j_kg: f64,
    ) -> Result<QuasiSteadyFilmState<'a>, ThermalInterferenceError> {
        let state = QuasiSteadyFilmState {
            body: self.body,
            source_clock: self.drive.source_clock,
            source_time_s: clock(source_time_s)?,
            phase: self
                .body
                .phase_curve()
                .state_at_specific_enthalpy(specific_enthalpy_j_kg)?,
        };
        self.admit_state(state)?;
        self.refresh_source(state.phase.temperature_k())?;
        Ok(state)
    }

    fn admit_state(&self, state: QuasiSteadyFilmState<'_>) -> Result<(), ThermalInterferenceError> {
        if !core::ptr::eq(state.body, self.body)
            || state.source_clock != self.drive.source_clock
            || state.phase.phase_curve_identity() != self.body.phase_curve().identity()
            || !state.phase.is_fully_solid()
            || state.phase.bulk_density_kg_m3() != self.density_kg_m3
        {
            return Err(bad(
                "accepted state must retain the original body, clock, chart, solid phase and density",
            ));
        }
        clock(state.source_time_s)?;
        Ok(())
    }

    fn admit_volume(&self, volume: Dd) -> Result<(), ThermalInterferenceError> {
        let expected = Dd::from_f64(self.body.mass_kg()) / Dd::from_f64(self.density_kg_m3);
        if !finite_dd(volume)
            || !finite_dd(expected)
            || volume.to_f64() <= 0.0
            || expected.to_f64() <= 0.0
            || ((volume - expected) / expected).to_f64().abs() > 64.0 * f64::EPSILON
        {
            return Err(bad(
                "thermal mass divided by original density must match the actual complete film volume",
            ));
        }
        Ok(())
    }

    fn refresh_source(
        &self,
        temperature_k: f64,
    ) -> Result<RefreshedFilmLaw, ThermalInterferenceError> {
        let original = self
            .original
            .layer(self.layer_index)
            .ok_or(bad("target film is absent"))?
            .response;
        RefreshedFilmLaw::at_temperature(original, self.card, self.temperature_axis, temperature_k)
    }

    /// Rebuild all original field owners at this accepted temperature.
    /// Repeated calls neither change source amplitudes nor reset their phase.
    pub fn heat_load<'p>(
        &'p self,
        state: QuasiSteadyFilmState<'a>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<QuasiSteadyFilmHeatLoad<'p, 'a, 'm>, ThermalInterferenceError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        let source = self.refresh_source(state.phase.temperature_k())?;
        let (spectrum, volume_m3, preparation_workspace_bytes) = {
            let factor_count = self.drive.factors.len();
            let metadata = checked_sum(&[
                bytes::<&ThermalInterferenceLine<'_>>(self.line_count)?,
                bytes::<f64>(self.line_count)?,
                bytes::<&PlanarDielectricLayers<'_>>(self.line_count)?,
                bytes::<PlanarAngularSpectrum<'_, '_>>(self.line_count)?,
                bytes::<PreparedAngularSpectrumIllumination<'_, '_, '_>>(self.line_count)?,
                bytes::<AngularSpectrumExposure<'_, '_, '_>>(factor_count + 1)?,
                bytes::<AngularSpectrumFactorStatistics>(factor_count)?,
                bytes::<DielectricLayer<'_>>(self.original.layer_count())?,
            ])?;
            let mut work =
                Workspace::new(metadata, self.limits.maximum_preparation_workspace_bytes)?;
            let mut frequencies = reserve(self.line_count)?;
            frequencies.extend(
                self.drive
                    .mean
                    .iter()
                    .chain(self.drive.factors.iter().flat_map(|factor| factor.lines))
                    .map(|line| line.angular_frequency_rad_s),
            );
            frequencies.sort_unstable_by(f64::total_cmp);
            frequencies.dedup_by(|left, right| *left == *right);
            poll(&mut cancelled)?;
            // Every rebuilt stack has the original layer count and geometry.
            // Its small native cache size is known before constructing it.
            work.charge(
                self.original
                    .workspace_bytes()
                    .checked_mul(frequencies.len())
                    .ok_or(bad("frequency stack cache bytes overflow"))?,
            )?;
            let mut films = reserve(self.original.layer_count())?;
            for index in 0..self.original.layer_count() {
                let original = self
                    .original
                    .layer(index)
                    .ok_or(bad("original film is absent"))?;
                films.push(DielectricLayer {
                    response: if index == self.layer_index {
                        source.response()
                    } else {
                        original.response
                    },
                    thickness_m: original.thickness_m,
                });
            }
            let mut stacks = reserve(frequencies.len())?;
            for &omega in &frequencies {
                poll(&mut cancelled)?;
                stacks.push(PlanarDielectricLayers::try_new_with_incident_response(
                    self.original.incident_response(),
                    self.original.exit_response(),
                    omega,
                    &films,
                )?);
            }
            let mut views = reserve(stacks.len())?;
            views.extend(stacks.iter());
            with_drive_ensemble(
                &views,
                &frequencies,
                self.grid,
                self.carrier_rad_m,
                self.drive,
                self.line_count,
                self.limits,
                &mut work,
                &mut cancelled,
                |ensemble, work, cancelled| -> Result<_, ThermalInterferenceError> {
                    let mut limits = self.limits.layer;
                    limits.maximum_workspace_bytes =
                        work.remaining().min(limits.maximum_workspace_bytes);
                    let layer =
                        ensemble.prepare_layer(self.layer_index, limits, &mut *cancelled)?;
                    work.charge(layer.workspace_bytes())?;
                    let mut limits = self.limits.heat;
                    limits.maximum_workspace_bytes =
                        work.remaining().min(limits.maximum_workspace_bytes);
                    let load = layer.integrated_heat_spectrum(
                        self.quadrature,
                        AngularSpectrumExposureHeatChannel::Total,
                        limits,
                        &mut *cancelled,
                    )?;
                    work.charge(load.preparation_workspace_bytes())?;
                    let volume = load.represented_volume_m3_dd();
                    self.admit_volume(volume)?;
                    Ok((load.into_spectrum(), volume, work.used))
                },
            )??
        };
        poll(&mut cancelled)?;
        Ok(QuasiSteadyFilmHeatLoad {
            owner: self,
            state,
            source,
            spectrum,
            volume_m3,
            preparation_workspace_bytes,
        })
    }

    /// Advance one declared feedback interval using the initial-temperature
    /// spectrum throughout. The native maximum thermal step and ripple bound
    /// may subdivide this interval without refreshing its held material law.
    ///
    /// config.initial_specific_enthalpy_j_kg must equal the accepted state.
    /// Its internal_power_w must be zero, as required by the native spectrum
    /// solver. The accepted DD clock is authoritative; controls.start_s is
    /// replaced by its rounded compatibility view. All thermal samples must
    /// stay within maximum_temperature_change_k of the initial temperature.
    /// Refusal asks the caller to shorten the feedback interval; no clipping,
    /// hidden rescaling, partial commit or implicit iteration occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'p>(
        &'p self,
        initial: QuasiSteadyFilmState<'a>,
        gate: BiotGate,
        config: LumpedEnthalpyMarchConfig,
        mut controls: LumpedEnthalpySpectrumControls,
        maximum_temperature_change_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<QuasiSteadyFilmStep<'p, 'a, 'm>, ThermalInterferenceError> {
        poll(&mut cancelled)?;
        self.admit_state(initial)?;
        if config.initial_specific_enthalpy_j_kg != initial.phase.specific_enthalpy_j_kg()
            || !maximum_temperature_change_k.is_finite()
            || maximum_temperature_change_k < 0.0
        {
            return Err(bad(
                "initial enthalpy must match accepted state and temperature excursion must be finite and nonnegative",
            ));
        }
        let applied_load = self.heat_load(initial, &mut cancelled)?;
        controls.start_s = initial.source_time_s.to_f64();
        let march = solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable(
            LumpedEnthalpyPowerBinding {
                body: self.body,
                spectrum: applied_load.spectrum(),
                input_unit: "W",
                watts_per_spectrum_unit: 1.0,
                source_clock: self.drive.source_clock,
            },
            gate,
            config,
            controls,
            initial.source_time_s,
            &mut cancelled,
        )?;
        let mut endpoint = initial;
        let mut endpoint_source = None;
        let mut maximum_temperature_excursion_k = 0.0_f64;
        // The initial sample is already admitted by heat_load. Every subsequent
        // native sample is checked, not merely the eventual endpoint.
        for sample in march.samples().iter().skip(1) {
            poll(&mut cancelled)?;
            let phase = sample.phase_state();
            let excursion = (phase.temperature_k() - initial.phase.temperature_k()).abs();
            if !excursion.is_finite() || excursion > maximum_temperature_change_k {
                return Err(bad(
                    "accepted thermal samples exceed the explicit feedback temperature-excursion limit",
                ));
            }
            maximum_temperature_excursion_k = maximum_temperature_excursion_k.max(excursion);
            endpoint = QuasiSteadyFilmState {
                body: self.body,
                source_clock: self.drive.source_clock,
                source_time_s: clock(initial.source_time_s + sample.elapsed_s_dd())?,
                phase,
            };
            self.admit_state(endpoint)?;
            endpoint_source = Some(self.refresh_source(phase.temperature_k())?);
        }
        let energy = march.energy();
        let maximum_biot = march.maximum_biot();
        let attempts = march.attempts();
        let spectral_work = march.spectral_work();
        let samples = march.into_samples();
        // A zero-duration result has the same constitutive endpoint. Refresh
        // through the original source owner, without manufacturing field work.
        let endpoint_source = match endpoint_source {
            Some(source) => source,
            None => self.refresh_source(initial.phase.temperature_k())?,
        };
        poll(&mut cancelled)?;
        Ok(QuasiSteadyFilmStep {
            applied_load,
            endpoint,
            endpoint_source,
            samples,
            energy,
            maximum_biot,
            maximum_temperature_excursion_k,
            attempts,
            spectral_work,
        })
    }
}

fn admit_thermal_context(
    material: &ResolvedMaterialStatePoint,
    temperature_axis: &str,
    receipt: &PropertyUsageReceipt,
) -> Result<(), ThermalInterferenceError> {
    for (axis, value) in &receipt.query_point {
        if axis == temperature_axis {
            if receipt.axis_quantities.get(axis) != material.axis_quantities().get(axis) {
                return Err(bad(
                    "thermal and electromagnetic temperature conventions differ",
                ));
            }
        } else if material.state_coordinate(axis) != Some(*value)
            || receipt.axis_quantities.get(axis) != material.axis_quantities().get(axis)
        {
            return Err(bad(
                "thermal and electromagnetic source fixed conditions differ",
            ));
        }
    }
    Ok(())
}

impl RefreshedFilmLaw {
    /// Refresh the complete original causal source, including mu and the
    /// same-state bundle shared by combined dielectric/conductor reservoirs.
    pub(super) fn at_temperature(
        original: DielectricResponse<'_>,
        card: &MaterialCard,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, ThermalInterferenceError> {
        let point = original
            .material()
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|_| bad("accepted temperature violates the original source coordinate"))?;
        Ok(match original {
            DielectricResponse::Debye(source) => {
                RefreshedFilmLaw::Debye(source.at_query_point(card, &point)?)
            }
            DielectricResponse::Lorentz(source) => {
                RefreshedFilmLaw::Lorentz(ResolvedLorentzMaterial::try_from_material(
                    source.material().at_query_point(card, &point)?,
                    source.high_frequency_property(),
                    source.pole_properties(),
                )?)
            }
            DielectricResponse::Ohmic(source) => {
                let material = source.material().at_query_point(card, &point)?;
                RefreshedFilmLaw::Ohmic(source.with_material(material)?)
            }
            DielectricResponse::DebyeWithOhmic {
                dielectric,
                conductor,
            } => {
                if dielectric.material().identity() != conductor.material().identity() {
                    return Err(bad(
                        "Debye and Ohmic laws must retain the same complete original snapshot",
                    ));
                }
                let dielectric = dielectric.at_query_point(card, &point)?;
                let conductor = conductor.with_material(dielectric.material().clone())?;
                RefreshedFilmLaw::DebyeWithOhmic(dielectric, conductor)
            }
            DielectricResponse::LorentzWithOhmic {
                dielectric,
                conductor,
            } => {
                if dielectric.material().identity() != conductor.material().identity() {
                    return Err(bad(
                        "Lorentz and Ohmic laws must retain the same complete original snapshot",
                    ));
                }
                // Resolve the original complete bundle once. The same refreshed
                // snapshot owns oscillator coefficients, permeability and leakage.
                let dielectric = ResolvedLorentzMaterial::try_from_material(
                    dielectric.material().at_query_point(card, &point)?,
                    dielectric.high_frequency_property(),
                    dielectric.pole_properties(),
                )?;
                let conductor = conductor.with_material(dielectric.material().clone())?;
                RefreshedFilmLaw::LorentzWithOhmic(dielectric, conductor)
            }
            _ => {
                return Err(bad(
                    "feedback requires original Debye, Lorentz and/or Ohmic reservoirs",
                ));
            }
        })
    }
}

pub(super) fn admit_refresh_source(
    original: DielectricResponse<'_>,
    card: &MaterialCard,
    temperature_axis: &str,
) -> Result<(), ThermalInterferenceError> {
    if !matches!(
        original,
        DielectricResponse::Debye(_)
            | DielectricResponse::Lorentz(_)
            | DielectricResponse::Ohmic(_)
            | DielectricResponse::DebyeWithOhmic { .. }
            | DielectricResponse::LorentzWithOhmic { .. }
    ) {
        return Err(bad(
            "feedback requires original Debye, Lorentz and/or Ohmic reservoirs",
        ));
    }
    let material = original.material();
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if material.card_identity() != card.content_hash()
        || material.axis_quantities().get(temperature_axis) != Some(&absolute)
        || material.state_coordinate(temperature_axis).is_none()
        || material
            .properties()
            .iter()
            .any(|p| !p.has_only_fixed_external_conditions(Some(temperature_axis), false))
    {
        return Err(bad(
            "feedback needs its original card, an explicit absolute temperature and fixed external conditions",
        ));
    }
    Ok(())
}

pub(super) fn admit_drive(
    drive: ThermalInterferenceDrive<'_>,
    limits: ThermalInterferenceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, ThermalInterferenceError> {
    if drive.source_clock.trim().is_empty()
        || drive.factors.len() > 4096
        || drive.factors.len() > limits.ensemble.maximum_factors
    {
        return Err(bad(
            "an explicit clock and bounded factor count are required",
        ));
    }
    let mut line_count = drive.mean.len();
    for factor in drive.factors {
        poll(&mut *cancelled)?;
        line_count = line_count
            .checked_add(factor.lines.len())
            .ok_or(bad("source line count overflows"))?;
    }
    if line_count == 0 || line_count > 4096 || line_count > limits.maximum_lines {
        return Err(bad(
            "the complete drive needs one to 4096 budgeted authored lines",
        ));
    }
    for line in drive
        .mean
        .iter()
        .chain(drive.factors.iter().flat_map(|f| f.lines))
    {
        poll(&mut *cancelled)?;
        if !line.angular_frequency_rad_s.is_finite() || line.angular_frequency_rad_s <= 0.0 {
            return Err(bad(
                "every source line needs a finite positive native angular frequency",
            ));
        }
    }
    Ok(line_count)
}

/// Reuse the exact original two-face field and mean/real/proper bank preparation.
/// The callback cannot retain these temporary field caches; owned physical
/// results and source snapshots remain with their original coupling owner.
#[allow(clippy::too_many_arguments)]
pub(super) fn with_drive_ensemble<'s, 'm, C: FnMut() -> bool, R>(
    stacks: &[&'s PlanarDielectricLayers<'m>],
    frequencies: &[f64],
    grid: AngularSpectrumGrid,
    carrier_rad_m: [f64; 2],
    drive: ThermalInterferenceDrive<'_>,
    line_count: usize,
    owner_limits: ThermalInterferenceLimits,
    work: &mut Workspace,
    cancelled: &mut C,
    evaluate: impl FnOnce(&AngularSpectrumExposureEnsemble<'_, '_, 's, 'm>, &mut Workspace, &mut C) -> R,
) -> Result<R, ThermalInterferenceError> {
    let factor_count = drive.factors.len();
    let mut lines: Vec<&ThermalInterferenceLine<'_>> = reserve(line_count)?;
    lines.extend(drive.mean.iter());
    for factor in drive.factors {
        poll(&mut *cancelled)?;
        lines.extend(factor.lines.iter());
    }
    let mut beams = reserve(line_count)?;
    for line in &lines {
        poll(&mut *cancelled)?;
        let index = frequencies
            .binary_search_by(|omega| omega.total_cmp(&line.angular_frequency_rad_s))
            .map_err(|_| bad("authored frequency is absent from its native cache"))?;
        let mut limits = owner_limits.angular;
        limits.maximum_workspace_bytes = work.remaining().min(limits.maximum_workspace_bytes);
        let beam = PlanarAngularSpectrum::prepare_with_carrier(
            stacks[index],
            grid,
            carrier_rad_m,
            line.left_tangential_e_v_per_m,
            limits,
            &mut *cancelled,
        )?;
        work.charge(beam.workspace_bytes())?;
        beams.push(beam);
    }
    let mut illumination = reserve(line_count)?;
    for (beam, line) in beams.iter().zip(&lines) {
        poll(&mut *cancelled)?;
        let mut limits = owner_limits.illumination;
        limits.maximum_workspace_bytes = work.remaining().min(limits.maximum_workspace_bytes);
        let prepared =
            beam.prepare_illumination(line.right_tangential_e_v_per_m, limits, &mut *cancelled)?;
        work.charge(prepared.workspace_bytes())?;
        illumination.push(prepared);
    }
    let anchor = beams.first().ok_or(bad("a source line is required"))?;
    let mean_count = drive.mean.len();
    let mut limits = owner_limits.exposure;
    limits.maximum_workspace_bytes = work.remaining().min(limits.maximum_workspace_bytes);
    let mean = AngularSpectrumExposure::try_new_illumination(
        anchor,
        &illumination[..mean_count],
        limits,
        &mut *cancelled,
    )?;
    work.charge(mean.workspace_bytes())?;
    let mut factors = reserve(factor_count)?;
    let mut statistics = reserve(factor_count)?;
    let mut offset = mean_count;
    for factor in drive.factors {
        poll(&mut *cancelled)?;
        let end = offset + factor.lines.len();
        let mut limits = owner_limits.exposure;
        limits.maximum_workspace_bytes = work.remaining().min(limits.maximum_workspace_bytes);
        let prepared = AngularSpectrumExposure::try_new_illumination(
            anchor,
            &illumination[offset..end],
            limits,
            &mut *cancelled,
        )?;
        work.charge(prepared.workspace_bytes())?;
        factors.push(prepared);
        statistics.push(factor.statistics);
        offset = end;
    }
    let mut limits = owner_limits.ensemble;
    limits.maximum_workspace_bytes = work.remaining().min(limits.maximum_workspace_bytes);
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        anchor,
        if mean_count == 0 { None } else { Some(&mean) },
        &factors,
        &statistics,
        limits,
        &mut *cancelled,
    )?;
    work.charge(ensemble.workspace_bytes())?;

    Ok(evaluate(&ensemble, work, cancelled))
}

pub(super) struct Workspace {
    used: usize,
    maximum: usize,
}
impl Workspace {
    pub(super) fn new(used: usize, maximum: usize) -> Result<Self, ThermalInterferenceError> {
        if used > maximum {
            return Err(bad("feedback preparation exceeds its workspace limit"));
        }
        Ok(Self { used, maximum })
    }
    pub(super) fn used(&self) -> usize {
        self.used
    }
    pub(super) fn remaining(&self) -> usize {
        self.maximum - self.used
    }
    pub(super) fn charge(&mut self, bytes: usize) -> Result<(), ThermalInterferenceError> {
        self.used = self
            .used
            .checked_add(bytes)
            .ok_or(bad("feedback workspace size overflows"))?;
        if self.used > self.maximum {
            return Err(bad("feedback preparation exceeds its workspace limit"));
        }
        Ok(())
    }
}
pub(super) fn bytes<T>(count: usize) -> Result<usize, ThermalInterferenceError> {
    count
        .checked_mul(core::mem::size_of::<T>())
        .ok_or(bad("feedback metadata bytes overflow"))
}
pub(super) fn checked_sum(values: &[usize]) -> Result<usize, ThermalInterferenceError> {
    values.iter().try_fold(0usize, |total, &value| {
        total
            .checked_add(value)
            .ok_or(bad("feedback metadata bytes overflow"))
    })
}
pub(super) fn reserve<T>(capacity: usize) -> Result<Vec<T>, ThermalInterferenceError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| bad("feedback allocation failed"))?;
    Ok(result)
}
fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn clock(value: Dd) -> Result<Dd, ThermalInterferenceError> {
    if !finite_dd(value) {
        return Err(bad("absolute source clock must be finite"));
    }
    let normalized = Dd::from_pair(value.hi, value.lo);
    if !finite_dd(normalized) {
        return Err(bad("absolute source clock is unrepresentable"));
    }
    Ok(normalized)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermalInterferenceError> {
    if cancelled() {
        Err(ThermalInterferenceError::Cancelled)
    } else {
        Ok(())
    }
}
fn bad(message: &'static str) -> ThermalInterferenceError {
    ThermalInterferenceError::Invalid(message)
}
