//! Material string + material plate + optional body-air cavity in one network.
//! No authored modal frequencies, effective masses, bridge shapes or output gains.

use super::*;
use crate::cavity_modal::HelmholtzModalCavity;
use crate::thin_plate::{ResolvedPlateSpecimen, certified_radiators};

/// Source-bound receiver and physical bridge/observer geometry.
/// The bridge is a bilateral compliant support on the FULL string span. The
/// plate's force footprint is the existing uniform-plate reducer's footprint.
/// This is not a carved shell, soundpost or bass bar model.
#[derive(Clone, Debug)]
pub struct MaterialPlateBody {
    /// Exact material state, thickness/mass constraint, support and geometry.
    pub specimen: ResolvedPlateSpecimen,
    /// Interior station on the full string, including its afterlength.
    pub bridge_station_fraction: f64,
    /// Transverse connection stiffness [N/m], not inferred from either bulk E.
    pub stiffness_n_m: f64,
    /// Connection damping [N s/m], distinct from the two materials' losses.
    pub damping_n_s_m: f64,
    /// Stress-free string-minus-body displacement [m].
    pub rest_extension_m: f64,
    /// Limits applied to the derived mass-normalized receiver modes.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Observer distance [m] in the existing compact baffled half-space law.
    pub listener_m: f64,
}

/// A source-derived lumped cavity on the back of the radiating plate.
/// The plate's positive displacement is OUTWARD from the cavity. Its volume
/// change and outward neck volume use the same signed area projections.
#[derive(Clone, Debug)]
pub struct MaterialBodyCavity {
    /// Existing Helmholtz law, geometry, selected loss and coherent INTERIOR
    /// gas. A distinct exterior requires an exterior-owned mouth and equal
    /// reference pressure; legacy compact necks retain their one-gas law.
    pub cavity: HelmholtzModalCavity,
    /// Limits for the physical neck inertia. Compression energy belongs to the
    /// network's collective spring and is also checked by the network budget.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Neck-to-listener distance [m], in the compact baffled observation.
    pub listener_m: f64,
}

impl ReciprocalMaterialString {
    /// Legacy air convenience form of [`Self::from_material_plate_with_gas`].
    /// Derive receiver modes/ports/losses from its original material and geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, ambient: AmbientGas,
        body: MaterialPlateBody, law: FrictionLaw, options: ReciprocalStringOptions,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let gas = crate::acoustic_realize::gas_state(ambient).map_err(ReciprocalStringError::Material)?;
        let mut model = Self::from_material_plate_with_gas(cx, specimen, gas, body, law, options, gate)?;
        model.source.ambient = Some(ambient);
        Ok(model)
    }

    /// Derive reciprocal string/plate mechanics and sound in one explicit gas.
    /// `options.connections` must be empty: the physical plate reducer supplies
    /// the bridge footprint, not caller-authored modal shapes. No body frequency
    /// or gain is tuned. Both original material specimens survive scheduling.
    /// The same coherent gas supplies string drag and compact pressure density.
    /// No piston reaction filter, room response or live gas replacement is added.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_with_gas(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, gas: GasState,
        body: MaterialPlateBody, law: FrictionLaw, options: ReciprocalStringOptions,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate(cx, specimen, gas, body, None, law, options, gate)
    }

    /// Body-air feedback in the SAME string, bridge, friction and receiver solve.
    ///
    /// The body mode's actual signed area drives a collective cavity compression
    /// spring. That same spring reacts on the body AND the physical neck mass;
    /// neck resistance dissipates through the original network's dashpot law.
    /// The neck is not a second oscillator with an assigned Helmholtz pitch.
    /// Its resonance, loading of the body, and resulting string motion emerge
    /// from the coupled mass/stiffness/geometry and the supplied gas.
    ///
    /// Output is the compact, zero-retardation sum of actual body and outward
    /// neck accelerations in Pa, including all accepted bridge/cavity reactions.
    /// This is NOT a broadband matched radiation field, spatial body acoustics,
    /// f-hole shape model, carved shell or experimental violin/cello/bass claim.
    /// The chosen lumped law and compact observations retain their validity scope.
    /// Changing data requires a new realization; no vibration history is erased
    /// by an in-place material or medium swap.
    ///
    /// # Errors
    /// Conflicting connection inputs, stale cavity gas, invalid geometry/limits,
    /// unsupported material image, aliasing, source/reduction failure or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_with_cavity(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, gas: GasState,
        body: MaterialPlateBody, cavity: MaterialBodyCavity, law: FrictionLaw,
        options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate(cx, specimen, gas, body, Some(cavity), law, options, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn build_material_plate(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, gas: GasState,
        body: MaterialPlateBody, cavity: Option<MaterialBodyCavity>, law: FrictionLaw,
        mut options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if !options.connections.is_empty() {
            return Err(invalid("material plate construction derives its bridge; supplied connections must be empty"));
        }
        gas.validate_coherence().map_err(|e| ReciprocalStringError::Material(AcousticRealizeError::Ambient(e.to_string())))?;
        let description = body.specimen.plate();
        let count = specimen.string().n_modes.checked_add(description.n_modes)
            .and_then(|n| n.checked_add(usize::from(cavity.is_some())))
            .ok_or_else(|| invalid("material string/plate/cavity mode count overflow"))?;
        let needed_links = cavity.as_ref().map_or(1, |c|
            2 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0));
        if description.geometric_nonlinearity || description.n_modes == 0
            || count > options.network.max_modes || needed_links > options.network.max_connections
            || options.sample_rate_hz == 0 || !interior(body.bridge_station_fraction)
            || !body.listener_m.is_finite() || body.listener_m <= 0.0 {
            return Err(invalid("material string/plate requires linear modes within budget, a positive clock/distance and an interior bridge"));
        }
        if let Some(c) = &cavity {
            c.cavity.validate_exterior_gas(&gas).map_err(ReciprocalStringError::Material)?;
            if !c.listener_m.is_finite() || c.listener_m <= 0.0 {
                return Err(invalid("body cavity requires a positive finite observer distance"));
            }
        }
        let image = reduce_material_body(&body.specimen, &gas, options.sample_rate_hz,
            body.modal_budget, body.listener_m, cavity.as_ref(), gate)?;
        let mut weights = vec![0.0; specimen.string().n_modes];
        weights.extend(image.pressure);
        options.connections.push(InteriorStringConnection {
            station_fraction: body.bridge_station_fraction,
            body: ModalAttachment { component: 0, shapes: image.bridge },
            stiffness_n_m: body.stiffness_n_m, damping_n_s_m: body.damping_n_s_m,
            rest_extension_m: body.rest_extension_m,
        });
        let mut assembled = Self::from_specimen_and_body_network(cx, specimen, gas,
            image.bodies, image.body_connections, law, options, gate)?;
        assembled.system = assembled.system.with_acceleration_pressure(weights)
            .map_err(ReciprocalStringError::Friction)?;
        assembled.source.material_plate = Some(body);
        assembled.source.material_cavity = cavity;
        poll(gate)?;
        Ok(assembled)
    }
}

// One body reduction for one-string and sympathetic multi-string realizations.
// Modal mass, signed bridge/volume shapes, cavity reactions and pressure output
// come from the same physical owners and their explicitly bound gas regions.
pub(super) struct MaterialBodyImage {
    pub(super) bodies: Vec<ModalAcousticTimeModel>,
    pub(super) bridge: Vec<f64>,
    pub(super) pressure: Vec<f64>,
    pub(super) body_connections: Vec<crate::render::schedule::force::coupled::ModalConnection>,
}

pub(super) fn reduce_material_body(
    specimen: &ResolvedPlateSpecimen, gas: &GasState, sample_rate_hz: u32,
    modal_budget: ModalAcousticTimeBudget, listener_m: f64,
    cavity: Option<&MaterialBodyCavity>, gate: &CancelGate,
) -> Result<MaterialBodyImage, ReciprocalStringError> {
    poll(gate)?;
    let description = specimen.plate();
    if description.geometric_nonlinearity || description.n_modes == 0 {
        return Err(invalid("shared material body requires nonempty linear modes"));
    }
    validate_body_medium(gas, sample_rate_hz, listener_m, cavity, gate)?;
    // Fresh reduction only; the original DKT/modal owner is synchronous.
    let radiators = certified_radiators(description).map_err(ReciprocalStringError::Material)?;
    reduce_radiating_body(&radiators, gas, sample_rate_hz, modal_budget, listener_m, cavity, None, gate)
}

// Admission before either uniform or regional geometry enters its eigensolve.
pub(super) fn validate_body_medium(
    gas: &GasState, sample_rate_hz: u32, listener_m: f64,
    cavity: Option<&MaterialBodyCavity>, gate: &CancelGate,
) -> Result<(), ReciprocalStringError> {
    poll(gate)?;
    gas.validate_coherence().map_err(|error| ReciprocalStringError::Material(
        AcousticRealizeError::Ambient(error.to_string())))?;
    if sample_rate_hz == 0 || !listener_m.is_finite() || listener_m <= 0.0 {
        return Err(invalid("shared material body requires a positive clock and observer distance"));
    }
    if let Some(c) = cavity {
        c.cavity.validate_exterior_gas(gas).map_err(ReciprocalStringError::Material)?;
        if !c.listener_m.is_finite() || c.listener_m <= 0.0 {
            return Err(invalid("body cavity requires a positive finite observer distance"));
        }
    }
    Ok(())
}

// A fresh modal image from either geometry owner enters exactly the same
// mass-normalization, signed volume, cavity and acceleration-pressure lowering.
// An explicit cavity row uses the same mass-normalized mechanical basis but may
// address different physical facets. Without it, uniform/chart plates preserve
// the existing coincident front/back midsurface trace.
// No receiver states or radiation-filter history may enter this private path.
#[allow(clippy::too_many_arguments)]
pub(super) fn reduce_radiating_body(
    radiators: &[crate::thin_plate::CompactBody], gas: &GasState, sample_rate_hz: u32,
    modal_budget: ModalAcousticTimeBudget, listener_m: f64,
    cavity: Option<&MaterialBodyCavity>, cavity_volume_shapes: Option<&[f64]>, gate: &CancelGate,
) -> Result<MaterialBodyImage, ReciprocalStringError> {
    validate_body_medium(gas, sample_rate_hz, listener_m, cavity, gate)?;
    if radiators.is_empty() { return Err(invalid("shared body reduction returned no modes")); }
    if cavity_volume_shapes.is_some_and(|row| row.len()!=radiators.len()
        || row.iter().any(|value| !value.is_finite())) {
        return Err(invalid("cavity pressure row must use every finite original body coordinate"));
    }
    let mut modes = Vec::with_capacity(radiators.len());
    let mut bridge = Vec::with_capacity(radiators.len());
    let mut areas = Vec::with_capacity(radiators.len());
    let mut pressure = Vec::with_capacity(radiators.len() + usize::from(cavity.is_some()));
    for radiator in radiators {
        poll(gate)?;
        let root_mass = fs_math::det::sqrt(radiator.mass_kg);
        let participation = radiator.drive_participation / root_mass;
        let area = radiator.area_m2 / root_mass;
        let weight = radiator.radiate_mass_normalized(1.0, gas.density, listener_m);
        if !(root_mass.is_finite() && root_mass > 0.0 && participation.is_finite()
            && area.is_finite() && weight.is_finite()) {
            return Err(invalid("material plate port normalization must remain finite"));
        }
        bridge.push(participation);
        areas.push(area);
        pressure.push(weight);
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: radiator.omega, damping_ratio: radiator.zeta,
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        });
    }
    let receiver = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, modal_budget)
        .map_err(|source| ReciprocalStringError::Material(AcousticRealizeError::Modal(source)))?;
    let mut bodies = vec![receiver];
    let mut body_connections = Vec::new();
    if let Some(c) = cavity {
        bodies.push(c.cavity.neck_component(sample_rate_hz, c.modal_budget)
            .map_err(ReciprocalStringError::Network)?);
        body_connections = c.cavity.connections(0, 1, cavity_volume_shapes.unwrap_or(&areas))
            .map_err(ReciprocalStringError::Network)?;
        // Existing compact half-space observation, with the physical outward
        // neck volume coordinate. No assigned oscillator frequency or gain.
        pressure.push(gas.density * c.cavity.neck_volume_shape_m2_per_sqrt_kg()
            / (core::f64::consts::TAU * c.listener_m));
    }
    poll(gate)?;
    Ok(MaterialBodyImage { bodies, bridge, pressure, body_connections })
}
