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
//! This is a frozen-optical-coefficient approximation, not an implicit
//! thermo-optic solve. Refreshing optical sources after temperature/material
//! change is the caller's responsibility. No coherent thermal emission,
//! broadband emissivity, phase change, or mechanical energy feedback is inferred.

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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        let layers = body.conduction().layers();
        if body.cell_count() > MAX_LAYERS || responses.len() != layers.len() {
            return Err(LayeredWaveHeatingError::Invalid(
                "one response per original layer and at most 64 thermal cells required",
            ));
        }
        let mut films = Vec::with_capacity(body.cell_count());
        for (i, (&response, layer)) in responses.iter().zip(layers).enumerate() {
            poll(&mut cancelled)?;
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
            let count = body.layer_cell_counts()[i];
            let thickness_m = layer.thickness_m() / count as f64;
            for _ in 0..count {
                films.push(DielectricLayer {
                    response,
                    thickness_m,
                });
            }
        }
        let optical = PlanarDielectricLayers::try_new(incident, exit, omega_rad_s, &films)?;
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
        let dt = step.duration_s();
        let n = self.owner.cell_ids.len();
        let mut energy = [0.0; MAX_LAYERS];
        let energy = &mut energy[..n];
        let optical_power = self.cell_power_into(drive, energy, &mut cancelled)?;
        let mut absorbed = Dd::ZERO;
        for value in energy.iter_mut() {
            *value *= dt;
            if !value.is_finite() {
                return Err(LayeredWaveHeatingError::Invalid(
                    "optical cell energy is unrepresentable",
                ));
            }
            absorbed = absorbed + Dd::from_f64(*value);
        }
        let convert = |values: [f64; 2]| values.map(|p| p * self.owner.area_m2() * dt);
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
        let step = step.with_cell_energy_j(energy, || checkpoint(&mut cancelled))?;
        let (thermal, _, heat_report) =
            step.solve(controls, workspace, || checkpoint(&mut cancelled))?;
        poll(&mut cancelled)?;
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
            cell_energy_j: self
                .owner
                .cell_ids
                .iter()
                .copied()
                .zip(energy.iter().copied())
                .collect(),
        })
    }
}
