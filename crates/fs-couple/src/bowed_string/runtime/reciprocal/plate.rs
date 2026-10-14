//! Material string + material plate, lowered into the existing shared network.
//! No authored modal frequencies, effective mass, bridge modal shape or gain.

use super::*;
use crate::thin_plate::{ResolvedPlateSpecimen, certified_radiators};

/// Source-bound receiver and physical bridge/observer geometry.
/// The bridge is a bilateral compliant support on the FULL string span. The
/// plate's force footprint is the existing uniform-plate reducer's footprint.
/// This is not a carved shell, soundpost, bass bar or body-air cavity model.
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

impl ReciprocalMaterialString {
    /// Derive a two-way radiating string/body performance from material sources.
    ///
    /// `options.connections` must be empty: this constructor derives the one
    /// bridge port from the plate's actual force footprint and modal mass. No
    /// caller-supplied attachment shapes are ignored or substituted. All other
    /// controls and budgets retain the common constructor's meanings.
    ///
    /// Changing the plate's source data changes its mass, stiffness, damping,
    /// retained modes, bridge participation and signed radiating areas together.
    /// The same admitted atmosphere supplies the string's drag and the body's
    /// compact pressure density. Both material specimens survive scheduling.
    ///
    /// The newly reduced plate starts at rest and has no piston reaction filter.
    /// The pressure readout uses actual endpoint acceleration INCLUDING bridge
    /// reactions. It is not a frequency-assigned oscillator bank or one-way body
    /// filter. This does not add fluid back-loading, propagation delay, room
    /// acoustics, nonlinear plates or empirical instrument validation. Changing
    /// material/environment requires a new realization, not live state transfer.
    ///
    /// # Errors
    /// Conflicting supplied connections, nonlinear plates, invalid geometry or
    /// clock, material/reduction refusal, aliasing, budgets and cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, ambient: AmbientGas,
        body: MaterialPlateBody, law: FrictionLaw, mut options: ReciprocalStringOptions,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if !options.connections.is_empty() {
            return Err(invalid("material plate construction derives its bridge; supplied connections must be empty"));
        }
        let description = body.specimen.plate();
        let count = specimen.string().n_modes.checked_add(description.n_modes)
            .ok_or_else(|| invalid("material string/plate mode count overflow"))?;
        if description.geometric_nonlinearity || description.n_modes == 0
            || count > options.network.max_modes || options.network.max_connections == 0
            || options.sample_rate_hz == 0 || !interior(body.bridge_station_fraction)
            || !body.listener_m.is_finite() || body.listener_m <= 0.0 {
            return Err(invalid("material string/plate requires linear modes within budget, a positive clock/distance and an interior bridge"));
        }
        let gas = crate::acoustic_realize::gas_state(ambient).map_err(ReciprocalStringError::Material)?;
        // Only fresh geometry reductions enter here. No existing body state or
        // radiation-filter memory is discarded during conversion.
        let radiators = certified_radiators(description).map_err(ReciprocalStringError::Material)?;
        poll(gate)?; // the existing DKT/modal reduction is synchronous
        let mut modes = Vec::with_capacity(radiators.len());
        let mut bridge = Vec::with_capacity(radiators.len());
        let mut weights = vec![0.0; specimen.string().n_modes];
        for radiator in &radiators {
            let root_mass = fs_math::det::sqrt(radiator.mass_kg);
            let participation = radiator.drive_participation / root_mass;
            // Reuse the pressure owner, including its half-space factor and sign.
            let weight = radiator.radiate_mass_normalized(1.0, gas.density, body.listener_m);
            if !(root_mass.is_finite() && root_mass > 0.0 && participation.is_finite() && weight.is_finite()) {
                return Err(invalid("material plate port normalization must remain finite"));
            }
            bridge.push(participation);
            weights.push(weight);
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: radiator.omega, damping_ratio: radiator.zeta,
                pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
            });
        }
        let receiver = ModalAcousticTimeModel::try_new(options.sample_rate_hz, modes, body.modal_budget)
            .map_err(|source| ReciprocalStringError::Network(ModalCouplingError::Component { component: 1, source }))?;
        options.connections.push(InteriorStringConnection {
            station_fraction: body.bridge_station_fraction,
            body: ModalAttachment { component: 0, shapes: bridge },
            stiffness_n_m: body.stiffness_n_m, damping_n_s_m: body.damping_n_s_m,
            rest_extension_m: body.rest_extension_m,
        });
        let mut assembled = Self::from_specimen(cx, specimen, ambient, vec![receiver], law, options, gate)?;
        assembled.system = assembled.system.with_acceleration_pressure(weights)
            .map_err(ReciprocalStringError::Friction)?;
        assembled.source.material_plate = Some(body);
        poll(gate)?;
        Ok(assembled)
    }
}
