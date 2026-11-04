//! Coherent dielectric absorption on the original spatial thermal timeline.
//!
//! Each homogeneous optical layer is split at its existing thermal cell faces.
//! The inserted equal-medium interfaces are identities. The original Maxwell
//! owner's positive slab loss modes therefore integrate the actual standing
//! wave over EACH control volume, without point quadrature, uniform heat
//! partitioning, or subtraction of nearly equal incoming/outgoing fluxes.
//!
//! Optical coefficients, geometry and two-face source statistics are held over
//! the explicit physical interval. The native heat owner retains nonlinear
//! capacity/conduction, contacts, reservoirs, authored sources and acceptance.
//! The base owner holds the supplied optical coefficients. The thermo_optical
//! wrapper refreshes each original source at its native cell temperature before
//! each interval and admits its endpoint state before publishing the result.
//! Neither owner supplies an implicit thermo-optic solve, coherent thermal
//! emission, broadband emissivity, phase change, or mechanical energy feedback.

use fs_em::plane_wave::{
    DielectricResponse, Polarization, ReflectionError,
    layers::{
        DielectricLayer, IlluminationBudget, LayerIlluminationPower, MAX_LAYERS,
        PlanarDielectricLayers, PreparedLayerIllumination,
    },
};
use fs_matdb::QueryPoint;
use fs_material::{
    state_point::ResolvedMaterialStatePoint,
    thermal::{
        conduction::LayerConductionConfig,
        transport::stack::transient::{
            BoundLayeredThermalBody,
            implicit::{
                LayeredThermalSolveControls, LayeredThermalSolveError, LayeredThermalSolveReport,
                LayeredThermalStep, LayeredThermalWorkspace,
            },
        },
    },
};
use fs_math::{c64::C64, dd::Dd};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

/// Explicit accepted-temperature feedback for multi-frequency Debye/Ohmic film heat.
pub mod interference;
pub mod spatial_interference;

/// Accepted-temperature optical refresh on the same spatial thermal timeline.
pub mod thermo_optical;
pub use thermo_optical::{
    BroadbandThermoOpticalLimits, BroadbandThermoOpticalWaveHeating,
    BroadbandThermoOpticalWaveStep, SpectralThermoOpticalSource, ThermoOpticalLayer,
    ThermoOpticalSource, ThermoOpticalSpectralLine, ThermoOpticalWaveHeating,
    ThermoOpticalWaveStep,
};

/// One monochromatic two-face drive in peak tangential E [V/m]. Both sources
/// share frequency, transverse wavevector, polarization and exp(+i omega t).
#[derive(Clone, Copy, Debug)]
pub enum LayeredWaveDrive<'a> {
    Coherent([C64; 2]),
    /// Explicit spectral-line second moment C=B B*. Each row is one coherent
    /// source factor; different factors are uncorrelated. This is not PSD/Hz.
    Coherence {
        factors: &'a [[C64; 2]],
        budget: IlluminationBudget,
    },
    /// Explicitly independent face sources; unknown correlation is not inferred.
    Independent([C64; 2]),
}

#[derive(Debug)]
pub enum LayeredWaveHeatingError {
    Invalid(&'static str),
    Electromagnetic(ReflectionError),
    Material(fs_material::state_point::MaterialStatePointError),
    Thermal(DepositionError),
    Solve(LayeredThermalSolveError),
    Cancelled,
}
impl core::fmt::Display for LayeredWaveHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-LAYERED-WAVE-HEATING: {self:?}")
    }
}
impl std::error::Error for LayeredWaveHeatingError {}
impl From<ReflectionError> for LayeredWaveHeatingError {
    fn from(value: ReflectionError) -> Self {
        Self::Electromagnetic(value)
    }
}
impl From<fs_material::state_point::MaterialStatePointError> for LayeredWaveHeatingError {
    fn from(value: fs_material::state_point::MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl From<DepositionError> for LayeredWaveHeatingError {
    fn from(value: DepositionError) -> Self {
        match value {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<LayeredThermalSolveError> for LayeredWaveHeatingError {
    fn from(value: LayeredThermalSolveError) -> Self {
        match value {
            LayeredThermalSolveError::Thermal(DepositionError::Cancelled) => Self::Cancelled,
            other => Self::Solve(other),
        }
    }
}

/// Original thermal geometry and source objects with one exact optical film
/// per thermal cell. Preparation is O(cells), capped at the native 64 films
/// BEFORE allocation. Repeated source queries never rebuild a Maxwell solve.
pub struct LayeredWaveHeating<'b, 'a, 'm> {
    thermal: &'b BoundLayeredThermalBody<'a>,
    optical: PlanarDielectricLayers<'m>,
    cell_ids: Vec<CellId>,
}

/// A single prepared source-to-cell loss operator, bound to the original body.
pub struct PreparedLayeredWaveHeating<'p, 'b, 'a, 'm> {
    owner: &'p LayeredWaveHeating<'b, 'a, 'm>,
    illumination: PreparedLayerIllumination<'p, 'm>,
}

/// Physical wave work integrated over one explicit interval. Only absorbed_j
/// enters thermal storage. Outgoing reflected/transmitted energy is not heat;
/// it is never debited again by this consumer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredWaveEnergy {
    pub incoming_j: [f64; 2],
    pub outgoing_j: [f64; 2],
    pub absorbed_j: f64,
    pub residual_j: f64,
}

/// One jointly returned candidate. The native heat account also includes any
/// pre-existing sources/reservoirs; optical is its distinct additional input,
/// not a second heat credit. No accepted history changes on failure.
pub struct LayeredWaveStep<'m> {
    pub thermal: DepositionSimulation<'m>,
    pub heat_report: LayeredThermalSolveReport,
    pub optical: LayeredWaveEnergy,
    pub optical_power: LayerIlluminationPower,
    /// One physical heat contribution [J] per original CellId, in native order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LayeredWaveHeatingError> {
    if cancelled() {
        Err(LayeredWaveHeatingError::Cancelled)
    } else {
        Ok(())
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_context(
    optical: &ResolvedMaterialStatePoint,
    thermal: &QueryPoint,
    temperature_axis: &str,
) -> Result<(), LayeredWaveHeatingError> {
    for (axis, value) in optical.query_point() {
        if let Some(other) = thermal.axes().get(axis) {
            if optical.axis_quantities().get(axis) != thermal.axis_quantities().get(axis)
                || (axis != temperature_axis && value != other)
            {
                return Err(LayeredWaveHeatingError::Invalid(
                    "optical and caloric/conductivity sources disagree on a shared fixed condition",
                ));
            }
        }
    }
    Ok(())
}

impl<'b, 'a, 'm> LayeredWaveHeating<'b, 'a, 'm> {
    /// Bind one held optical law per ORIGINAL material layer, in native order.
    /// Card identities and shared fixed query coordinates/conventions must
    /// agree with conductivity AND capacity. The optical temperature remains
    /// explicit; matching cards do not establish an experimentally matched coupon.
    /// Area, thickness, subdivision and cell identities come only from thermal.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: &'m ResolvedMaterialStatePoint,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        responses: &[DielectricResponse<'m>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        Self::try_new_with_incident_response(
            thermal,
            DielectricResponse::Lossless(incident),
            exit,
            omega_rad_s,
            responses,
            cancelled,
        )
    }

    /// Bind original layer laws with a source-backed dispersive incident
    /// medium. Its original epsilon(omega) and same-state mu determine phase,
    /// transverse wave number and physical port power at this frequency.
    /// The native Maxwell owner requires positive, exactly lossless incidence;
    /// no additional static epsilon, averaged index or exterior heat is inferred.
    pub fn try_new_with_incident_response(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        responses: &[DielectricResponse<'m>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        Self::from_responses(
            thermal,
            incident,
            exit,
            omega_rad_s,
            responses,
            false,
            cancelled,
        )
    }

    /// Bind one held source response per ORIGINAL thermal cell, in thickness
    /// order. This resolves a spatial optical state without replacing the
    /// thermal geometry or averaging different cell temperatures into one film.
    /// Each response must retain its own layer's material card and fixed context.
    pub fn try_new_cellwise(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: &'m ResolvedMaterialStatePoint,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        responses: &[DielectricResponse<'m>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        Self::try_new_cellwise_with_incident_response(
            thermal,
            DielectricResponse::Lossless(incident),
            exit,
            omega_rad_s,
            responses,
            cancelled,
        )
    }

    /// Cellwise counterpart retaining the original incident constitutive law.
    /// Every thermal cell remains its original homogeneous optical sublayer;
    /// source support and lossless propagating-port admission are unchanged.
    pub fn try_new_cellwise_with_incident_response(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        responses: &[DielectricResponse<'m>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        Self::from_responses(
            thermal,
            incident,
            exit,
            omega_rad_s,
            responses,
            true,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_responses(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        responses: &[DielectricResponse<'m>],
        cellwise: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        if !body.conduction().geometry().is_planar() {
            return Err(LayeredWaveHeatingError::Invalid(
                "this planar coupling requires a planar distributed thermal body",
            ));
        }
        let layers = body.conduction().layers();
        let expected = if cellwise {
            body.cell_count()
        } else {
            layers.len()
        };
        if body.cell_count() > MAX_LAYERS || responses.len() != expected {
            return Err(LayeredWaveHeatingError::Invalid(
                "the selected original layer/cell response count and at most 64 thermal cells are required",
            ));
        }
        let mut films = Vec::with_capacity(body.cell_count());
        for (i, layer) in layers.iter().enumerate() {
            let count = body.layer_cell_counts()[i];
            let thickness_m = layer.thickness_m() / count as f64;
            for local in 0..count {
                poll(&mut cancelled)?;
                let response = responses[if cellwise { films.len() } else { i }];
                // A uniform layer needs one source admission. Spatial optical
                // states retain their separate source/condition checks.
                if cellwise || local == 0 {
                    let material = response.material();
                    if material.card_identity() != layer.material_card_identity() {
                        return Err(LayeredWaveHeatingError::Invalid(
                            "optical and thermal layer must retain the same original material card",
                        ));
                    }
                    admit_context(
                        material,
                        &layer.query_endpoints()[0],
                        layer.temperature_axis(),
                    )?;
                    for component in body.storages()[i].components() {
                        if material.card_identity() != component.material_card_identity() {
                            return Err(LayeredWaveHeatingError::Invalid(
                                "optical and thermal capacity material cards differ",
                            ));
                        }
                        admit_context(
                            material,
                            &component.query_endpoints()[0],
                            component.temperature_axis(),
                        )?;
                    }
                }
                films.push(DielectricLayer {
                    response,
                    thickness_m,
                });
            }
        }
        let optical = PlanarDielectricLayers::try_new_with_incident_response(
            incident,
            exit,
            omega_rad_s,
            &films,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            thermal,
            optical,
            cell_ids: thermal.cell_ids(),
        })
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }
    pub const fn optical_layers(&self) -> &PlanarDielectricLayers<'m> {
        &self.optical
    }
    pub fn cell_ids(&self) -> &[CellId] {
        &self.cell_ids
    }
    pub fn area_m2(&self) -> f64 {
        self.thermal.prepared_body().conduction().area_m2()
    }
    /// Two physical propagating exterior ports are required by the existing
    /// illumination owner; lossless/TIR/evanescent port admission is unchanged.
    pub fn prepare(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedLayeredWaveHeating<'_, 'b, 'a, 'm>, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let illumination = self
            .optical
            .prepare_illumination(incidence_cosine, polarization)?;
        poll(&mut cancelled)?;
        Ok(PreparedLayeredWaveHeating {
            owner: self,
            illumination,
        })
    }
}

impl<'p, 'b, 'a, 'm> PreparedLayeredWaveHeating<'p, 'b, 'a, 'm> {
    pub const fn owner(&self) -> &'p LayeredWaveHeating<'b, 'a, 'm> {
        self.owner
    }

    /// Allocation-free physical cell powers [W], staged transactionally.
    /// The returned original Maxwell account is normal power per area [W/m²].
    /// No source covariance is replaced by a diagonal cell-power approximation.
    pub fn cell_power_into(
        &self,
        drive: LayeredWaveDrive<'_>,
        cell_power_w: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerIlluminationPower, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let n = self.owner.cell_ids.len();
        if cell_power_w.len() != n {
            return Err(LayeredWaveHeatingError::Invalid(
                "one power per original cell required",
            ));
        }
        let mut temporary = [0.0; MAX_LAYERS];
        let heat = &mut temporary[..n];
        let power = match drive {
            LayeredWaveDrive::Coherent(input) => self.illumination.evaluate(input, heat)?.power,
            LayeredWaveDrive::Independent(input) => self
                .illumination
                .evaluate_independent_sources(input, heat)?,
            LayeredWaveDrive::Coherence { factors, budget } => {
                let mut interrupted = false;
                let result =
                    self.illumination
                        .evaluate_coherence_modes(factors, budget, heat, || {
                            let stop = cancelled();
                            interrupted |= stop;
                            stop
                        });
                if interrupted {
                    return Err(LayeredWaveHeatingError::Cancelled);
                }
                result?
            }
        };
        for value in heat.iter_mut() {
            *value *= self.owner.area_m2();
            if !value.is_finite() || *value < 0.0 {
                return Err(LayeredWaveHeatingError::Invalid(
                    "physical cell power is unrepresentable",
                ));
            }
        }
        poll(&mut cancelled)?;
        cell_power_w.copy_from_slice(heat);
        Ok(power)
    }

    /// Apply this held monochromatic/statistical illumination until target_s.
    /// Source statistics are stationary over the interval; the harmonic factor
    /// 1/2 is already included by Maxwell and is NEVER applied again here.
    /// Native thermal admission, nonlinear solve, events and energy acceptance
    /// remain authoritative. Every refused trial leaves initial unchanged.
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
    ) -> Result<LayeredWaveStep<'t>, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let step =
            LayeredThermalStep::new(self.owner.thermal, initial, target_s, conduction, || {
                checkpoint(&mut cancelled)
            })?;
        let n = self.owner.cell_ids.len();
        let mut power = [0.0; MAX_LAYERS];
        let power = &mut power[..n];
        let optical_power = self.cell_power_into(drive, power, &mut cancelled)?;
        solve_held_cell_power(
            step,
            self.owner.area_m2(),
            optical_power,
            power,
            controls,
            workspace,
            &mut cancelled,
        )
    }
}

/// Shared monochromatic/broadband handoff. The input buffer contains physical
/// cell watts and is consumed into joules for exactly one native thermal solve.
fn solve_held_cell_power<'t>(
    step: LayeredThermalStep<'_, 't, '_>,
    area_m2: f64,
    optical_power: LayerIlluminationPower,
    cell_power_w: &mut [f64],
    controls: LayeredThermalSolveControls,
    workspace: &mut LayeredThermalWorkspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<LayeredWaveStep<'t>, LayeredWaveHeatingError> {
    let dt = step.duration_s();
    let mut absorbed = Dd::ZERO;
    for value in cell_power_w.iter_mut() {
        *value *= dt;
        if !value.is_finite() || *value < 0.0 {
            return Err(LayeredWaveHeatingError::Invalid(
                "optical cell energy is unrepresentable",
            ));
        }
        absorbed = absorbed + Dd::from_f64(*value);
    }
    let convert = |values: [f64; 2]| values.map(|p| p * area_m2 * dt);
    let incoming_j = convert(optical_power.incoming_w_m2);
    let outgoing_j = convert(optical_power.outgoing_w_m2);
    let absorbed_j = absorbed.to_f64();
    let residual_j = (incoming_j
        .iter()
        .fold(Dd::ZERO, |s, &v| s + Dd::from_f64(v))
        - outgoing_j
            .iter()
            .fold(Dd::ZERO, |s, &v| s + Dd::from_f64(v))
        - absorbed)
        .to_f64();
    if incoming_j.iter().chain(&outgoing_j).any(|v| !v.is_finite())
        || !absorbed_j.is_finite()
        || !residual_j.is_finite()
    {
        return Err(LayeredWaveHeatingError::Invalid(
            "optical interval energy is unrepresentable",
        ));
    }
    let cell_energy_j = step
        .cell_ids()
        .iter()
        .copied()
        .zip(cell_power_w.iter().copied())
        .collect();
    let step = step.with_cell_energy_j(cell_power_w, || checkpoint(cancelled))?;
    let (thermal, _, heat_report) = step.solve(controls, workspace, || checkpoint(cancelled))?;
    poll(cancelled)?;
    Ok(LayeredWaveStep {
        thermal,
        heat_report,
        optical: LayeredWaveEnergy {
            incoming_j,
            outgoing_j,
            absorbed_j,
            residual_j,
        },
        optical_power,
        cell_energy_j,
    })
}
