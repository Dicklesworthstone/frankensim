//! Reciprocal uniform-pressure loading of the complete nonlinear plate state.
use super::{AcousticRealizeError, VkBody, prepare_render_linear, render_step_error};
use crate::cavity_modal::HelmholtzModalCavity;
use fs_exec::CancelGate;
use fs_phs::{PortHamiltonian, QuadraticStorage};

/// Independent accepted gas state when connecting the existing cavity owner.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VkCavityInitial {
    /// Internal pressure minus equilibrium pressure [Pa].
    pub pressure_increment_pa: f64,
    /// Neck flow OUTWARD from the cavity [m^3/s].
    pub neck_outflow_m3_s: f64,
}

/// Pure observation of the same accepted mechanical/gas state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VkCavityObservation {
    pub pressure_increment_pa: f64,
    pub neck_outflow_m3_s: f64,
    /// Sum of original surface velocities, positive OUTWARD [m^3/s].
    pub surface_volume_rate_m3_s: f64,
    /// Gas and neck energies are already included in the body's total storage.
    pub compression_energy_j: f64,
    pub neck_kinetic_energy_j: f64,
    /// Selected original neck resistance times squared neck flow [W].
    pub neck_loss_power_w: f64,
}

/// Observation handle for one cavity appended to a complete plate.
///
/// It owns no evolving state. Use it only with the body returned by the same
/// attachment (or that body's prepared state). Further appended loads preserve
/// its original offsets. Signed cavity shapes are independent of exterior
/// radiation areas and of the point-force drive projection.
#[derive(Clone, Debug)]
pub struct VkCavityLoad {
    source: HelmholtzModalCavity,
    offset: usize,
    outward_shapes: Vec<f64>,
    inverse_compliance: f64,
    acoustic_inertia: f64,
    acoustic_resistance: f64,
}

impl VkCavityLoad {
    #[must_use]
    pub fn source(&self) -> &HelmholtzModalCavity {
        &self.source
    }

    /// Read the accepted body or prepared-body state, without advancing it.
    /// This pressure is an internal reaction, not a microphone pressure or an
    /// additional force for the caller to apply again.
    pub fn observe(&self, state: &[f64]) -> Result<VkCavityObservation, AcousticRealizeError> {
        if state.len() < self.offset + 2 || state.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "cavity observation needs its complete finite accepted state",
            ));
        }
        let compression = state[self.offset];
        let momentum = state[self.offset + 1];
        let pressure = self.inverse_compliance * compression;
        let flow = -momentum / self.acoustic_inertia;
        let surface_volume_rate_m3_s = self
            .outward_shapes
            .iter()
            .enumerate()
            .map(|(k, a)| a * state[2 * k + 1])
            .sum();
        let result = VkCavityObservation {
            pressure_increment_pa: pressure,
            neck_outflow_m3_s: flow,
            surface_volume_rate_m3_s,
            compression_energy_j: 0.5 * pressure * compression,
            neck_kinetic_energy_j: -0.5 * flow * momentum,
            neck_loss_power_w: self.acoustic_resistance * flow * flow,
        };
        if [
            result.pressure_increment_pa,
            result.neck_outflow_m3_s,
            result.surface_volume_rate_m3_s,
            result.compression_energy_j,
            result.neck_kinetic_energy_j,
            result.neck_loss_power_w,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid("cavity observation is not representable"));
        }
        Ok(result)
    }
}

impl VkBody {
    /// Connect a material-gas cavity to this complete nonlinear plate.
    ///
    /// `outward_volume_shapes` integrates the signed cavity-side normal shapes
    /// in the SAME mass-normalized basis as the original plate. Positive surface
    /// motion expands the cavity: U_in = -A.v, while pressure applies +A.p to the
    /// plate. These are one work-conjugate pair, so internal work cancels inside
    /// the original pHS step. A sum of areas times the point-force port would
    /// not have this property and is never used.
    ///
    /// All original coordinates, full damping, nonlinear storage and material
    /// histories remain first, unchanged. The two appended coordinates are
    /// compression volume [m^3] and inward acoustic momentum [Pa s]. Their
    /// independent initial pressure/neck flow are explicit; attaching the load
    /// does not reset a moving plate or infer its prior gas history. The returned
    /// handle observes this gas state after ordinary or prepared body steps.
    ///
    /// Compliance, neck inertia, resistance and end corrections come from the
    /// supplied existing owner, counted once. Its uniform small-signal pressure,
    /// fixed geometry/gas and selected loss approximation remain in force. This
    /// does not add cavity standing waves, broadband neck loss or mouth sound.
    /// An ExteriorMouth source still needs its separate exterior neck load.
    ///
    /// # Errors
    /// Refuses malformed signed projections/initial data, total-state capacity,
    /// cancellation, unrepresentable physics or original pHS/render admission.
    pub fn with_cavity_load(
        self,
        source: &HelmholtzModalCavity,
        outward_volume_shapes: &[f64],
        initial: VkCavityInitial,
        maximum_total_states: usize,
        gate: &CancelGate,
    ) -> Result<(Self, VkCavityLoad), AcousticRealizeError> {
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let offset = self.sys.state_dim();
        let states = offset
            .checked_add(2)
            .ok_or_else(|| invalid("plate/cavity state extent overflow"))?;
        let block = states
            .checked_add(self.sys.port_dim())
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("plate/cavity preparation extent overflow"))?;
        if self.sys.port_dim() != 1
            || self.x.len() != offset
            || self.areas.is_empty()
            || self.areas.len().checked_mul(2).is_none_or(|n| n > offset)
            || outward_volume_shapes.len() != self.areas.len()
            || states > maximum_total_states
            || block > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || outward_volume_shapes
                .iter()
                .chain(&self.x)
                .chain(&self.areas)
                .any(|x| !x.is_finite())
            || !initial.pressure_increment_pa.is_finite()
            || !initial.neck_outflow_m3_s.is_finite()
        {
            return Err(invalid(
                "plate/cavity loading requires finite original projections/state within its total-state budget",
            ));
        }
        let area_squared = source.neck_area_m2().powi(2);
        let inverse_compliance = source.compression_stiffness_n_m() / area_squared;
        let acoustic_inertia = source.neck_mass_kg() / area_squared;
        let acoustic_resistance = source.neck_resistance_n_s_m() / area_squared;
        if ![
            inverse_compliance,
            acoustic_inertia,
            acoustic_inertia.recip(),
        ]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
            || !acoustic_resistance.is_finite()
            || acoustic_resistance < 0.0
        {
            return Err(invalid(
                "original cavity coefficients are not representable in pressure/flow coordinates",
            ));
        }
        // Exact flow-driven dual of the original cavity's quadratic storage.
        // This lowering also preserves ExteriorMouth's omitted outer-end mass.
        let load = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, acoustic_resistance],
            vec![1.0, 0.0],
            Box::new(
                QuadraticStorage::new(
                    vec![inverse_compliance, 0.0, 0.0, acoustic_inertia.recip()],
                    2,
                )
                .map_err(render_step_error)?,
            ),
        )
        .map_err(render_step_error)?;
        let mut column = vec![0.0; offset];
        for (k, a) in outward_volume_shapes.iter().enumerate() {
            column[2 * k + 1] = -a;
        }
        let observation = VkCavityLoad {
            source: source.clone(),
            offset,
            outward_shapes: outward_volume_shapes.to_vec(),
            inverse_compliance,
            acoustic_inertia,
            acoustic_resistance,
        };
        let mut state = self.x;
        state.extend([
            initial.pressure_increment_pa / inverse_compliance,
            -acoustic_inertia * initial.neck_outflow_m3_s,
        ]);
        observation.observe(&state)?;
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let sys = self
            .sys
            .with_parallel_load_columns(load, &column)
            .map_err(render_step_error)?;
        let linear = prepare_render_linear(&sys)?;
        if !sys.hamiltonian(&state).is_finite() {
            return Err(invalid("complete plate/cavity storage is not finite"));
        }
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        Ok((
            Self {
                sys,
                linear,
                x: state,
                areas: self.areas,
            },
            observation,
        ))
    }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

#[cfg(test)]
mod tests;
