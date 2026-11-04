//! Participation and physical tuning of complete outgoing-boundary poles.
//!
//! The ordinary physical pole owner and fs-la PreparedPole retain the original
//! closed pencil, its row conditioning, fixed gauge and one bordered LU.

use super::{DecayWavePort, ScatteringModel, ScatteringPole};
use crate::harmonic::response::{
    HarmonicPoleResponseError, HarmonicPoleResponseOptions, HarmonicPortResidue,
    PreparedHarmonicPole,
};
use crate::harmonic::{finite, poll};
use fs_math::c64::C64;

mod derivatives;
pub use derivatives::{
    DecayWaveParameterDerivative, ScatteringPoleDerivative, ScatteringPoleParameterDirection,
};

/// Original storage, material work and outgoing-boundary work in the prepared
/// mode's FIXED pivot gauge, v[pivot]=1. These are instantaneous complex-envelope
/// quantities for exp(s*t), not cycle-averaged powers of a driven response.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringPoleObservation {
    /// Original port order; selected inputs are -y/Z, ordinary inputs are zero.
    pub physical_inputs: Vec<C64>,
    /// Authored wave-port order. Incident waves vanish in a free mode.
    pub waves: Vec<DecayWavePort>,
    pub stored_energy_j: f64,
    /// Original R/T work, including signed incremental work on nonlinear images.
    pub internal_dissipated_power_w: f64,
    pub outgoing_power_w: f64,
    pub total_decay_power_w: f64,
    pub decay_power_residual_w: f64,
    pub physical_power_residual_w: f64,
    pub relative_residual: f64,
}

/// Laurent coefficients for an analytic incident/ordinary excitation at s0.
/// State/effort/output coefficients are in `physical`; actual input residues
/// are -Y_res/Z on selected ports and zero elsewhere. Outgoing residues equal
/// Y_res: analytic incident feedthrough has NO residue. These are not amplitudes
/// or watts and do not include the regular or algebraic transfer-function part.
#[derive(Clone, Debug)]
pub struct ScatteringPortResidue {
    pub physical: HarmonicPortResidue,
    pub physical_inputs: Vec<C64>,
    /// Coefficients in ScatteringModel::wave_ports() order.
    pub outgoing: Vec<C64>,
    /// Independent original homogeneous physical-equation residual.
    pub relative_residual: f64,
}

/// One immutable simple-pole factor, bound to this original loaded model.
/// Arbitrary admitted forcing and parameter batches reuse the same bordered LU.
/// Every descriptor reaction and material memory state remains in the source
/// coordinates; no reduced modal model or separate eigenproblem is assembled.
pub struct PreparedScatteringPole<'a> {
    model: &'a ScatteringModel,
    physical: PreparedHarmonicPole<'a>,
    observation: ScatteringPoleObservation,
    maximum_vectors: usize,
}

fn invalid(what: &'static str) -> HarmonicPoleResponseError {
    HarmonicPoleResponseError::Invalid { what }
}

impl ScatteringModel {
    /// Prepare forcing participation and parameter sensitivities at one finite
    /// loaded pole. The supplied s/state are rechecked against THIS complete
    /// closed model. Cached observations and stale residual fields are ignored.
    ///
    /// The existing physical/numerical pole owner supplies the fixed pivot gauge
    /// and bordered factor; no bare impedance, descriptor mass or eigenbasis is
    /// inverted. The same original outgoing-mode observer independently restores
    /// boundary inputs and material/escaped power in that prepared gauge.
    /// Dark internal poles remain valid and can have zero forced residues.
    ///
    /// The options bound the complete physical payload/work, including this
    /// boundary observer. In addition to the ordinary owner's bound, this
    /// reserves 64n+128m+256k+512+b*(32m+192k+512) bytes and
    /// 64*(b+1)*(n*n+n*m+m+k) work units, for b maximum vectors and k wave ports.
    /// Every full batch fits that declared ceiling; no column is truncated.
    /// Infinite, repeated, defective and unresolved poles retain the existing
    /// numerical refusals. These local checks do not certify a full spectrum,
    /// continuum validity, finite parameter changes or closed-loop passivity.
    pub fn prepare_pole_response(
        &self,
        pole: &ScatteringPole,
        options: HarmonicPoleResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedScatteringPole<'_>, HarmonicPoleResponseError> {
        poll(&mut cancelled)?;
        if !(1..=256).contains(&options.numerical.max_vectors) {
            return Err(invalid("complete bounded loaded pole batch required"));
        }
        let (n, m, k, b) = (
            self.closed.n as u128,
            self.closed.m as u128,
            self.ports.len() as u128,
            options.numerical.max_vectors as u128,
        );
        let bytes = 64 * n + 128 * m + 256 * k + 512 + b * (32 * m + 192 * k + 512);
        let work = 64 * (b + 1) * (n * n + n * m + m + k);
        if bytes > options.max_physical_bytes as u128 || work > options.max_physical_work as u128 {
            return Err(invalid(
                "insufficient complete loaded pole observation resources",
            ));
        }
        let mut numerical = options;
        numerical.max_physical_bytes -= bytes as usize;
        numerical.max_physical_work -= work as u64;
        let physical = self.closed.prepare_pole_response_at(
            pole.laplace_pole_per_s,
            &pole.state,
            numerical,
            &mut cancelled,
        )?;
        let mode = self.observe_outgoing_pole(physical.observed_mode(), &mut cancelled)?;
        let observation = ScatteringPoleObservation {
            physical_inputs: mode.physical_inputs,
            waves: mode.waves,
            stored_energy_j: mode.stored_energy_j,
            internal_dissipated_power_w: mode.internal_dissipated_power_w,
            outgoing_power_w: mode.outgoing_power_w,
            total_decay_power_w: mode.total_decay_power_w,
            decay_power_residual_w: mode.decay_power_residual_w,
            physical_power_residual_w: mode.physical_power_residual_w,
            relative_residual: mode.relative_residual,
        };
        poll(&mut cancelled)?;
        Ok(PreparedScatteringPole {
            model: self,
            physical,
            observation,
            maximum_vectors: options.numerical.max_vectors,
        })
    }
}

impl PreparedScatteringPole<'_> {
    pub fn laplace_pole_per_s(&self) -> C64 {
        self.physical.laplace_pole_per_s()
    }
    pub fn right_state(&self) -> &[C64] {
        self.physical.right_state()
    }
    pub fn mode_effort(&self) -> &[C64] {
        self.physical.mode_effort()
    }
    pub fn mode_output(&self) -> &[C64] {
        self.physical.mode_output()
    }
    /// Gauge is v[pivot]=1, held fixed by every parameter derivative.
    pub fn normalization_pivot(&self) -> usize {
        self.physical.normalization_pivot()
    }
    pub fn mode_observation(&self) -> &ScatteringPoleObservation {
        &self.observation
    }

    /// Laurent residue batches for complex excitations evaluated at this pole:
    /// incident output-wave amplitudes at selected original source slots and
    /// ordinary physical inputs elsewhere. The source owner maps these to
    /// 2a/Z and ordinary inputs before the same conditioned bordered solve.
    ///
    /// Original reaction/output coefficients and nonnormal participation remain.
    /// No nonzero observability threshold is imposed. The result is one pole's
    /// contribution, not the full forced response or a modal truncation. All
    /// direct incident and algebraic feedthrough belongs to the regular part.
    /// Failed or cancelled batches publish nothing and leave the factor reusable.
    pub fn port_residues(
        &self,
        inputs: &[&[C64]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<ScatteringPortResidue>, HarmonicPoleResponseError> {
        poll(&mut cancelled)?;
        if inputs.is_empty()
            || inputs.len() > self.maximum_vectors
            || inputs
                .iter()
                .any(|u| u.len() != self.model.closed.m || u.iter().any(|&z| !finite(z)))
        {
            return Err(invalid(
                "finite complete loaded excitation batch within its ceiling required",
            ));
        }
        let mut effective = Vec::with_capacity(inputs.len());
        for &input in inputs {
            let mut mapped = input.to_vec();
            for &(p, z) in &self.model.ports {
                poll(&mut cancelled)?;
                mapped[p] = input[p].scale(2.0 / z);
            }
            if mapped.iter().any(|&z| !finite(z)) {
                return Err(invalid("loaded pole excitation is not representable"));
            }
            effective.push(mapped);
        }
        let refs: Vec<_> = effective.iter().map(Vec::as_slice).collect();
        let physical = self.physical.port_residues(&refs, &mut cancelled)?;
        let mut result = Vec::with_capacity(physical.len());
        for physical in physical {
            poll(&mut cancelled)?;
            let mut physical_inputs = vec![C64::ZERO; self.model.closed.m];
            let mut outgoing = Vec::with_capacity(self.model.ports.len());
            for &(p, z) in &self.model.ports {
                poll(&mut cancelled)?;
                physical_inputs[p] = physical.output[p].scale(-1.0 / z);
                outgoing.push(physical.output[p]);
            }
            let relative_residual = self.model.original_residual(
                self.laplace_pole_per_s(),
                &physical.state,
                &physical.effort,
                &physical_inputs,
                &mut cancelled,
            )?;
            result.push(ScatteringPortResidue {
                physical,
                physical_inputs,
                outgoing,
                relative_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}
