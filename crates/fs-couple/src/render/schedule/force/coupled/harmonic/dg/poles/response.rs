//! Physical forcing and native parameter response of the complete coupled pencil.
//! The shared bordered-pole owner supplies the true left/right participation;
//! neither a right-mode force projection nor a reduced oscillator is substituted.

use super::*;
pub use fs_la::pole_response::PoleResponseOptions;
use fs_la::pole_response::{PoleResponseError, PreparedPole};
use fs_math::dd::Dd;

mod parameter;
pub use parameter::{DgLoadedModalPoleDerivative, DgPoleParameterDirection};

/// Separate limits for the shared numerical factor and complete physical maps.
/// They apply to preparation and independently to each complete response batch.
#[derive(Clone, Copy, Debug)]
pub struct DgPoleResponseOptions {
    pub numerical: PoleResponseOptions,
    /// Complete complex raw-state observations per batch, in 0..=256.
    pub max_observations: usize,
    /// Conservative physical payload, excluding source caches, caller inputs,
    /// allocator overhead and the separately admitted numerical factor.
    pub max_physical_bytes: usize,
    /// Mapping/audit units: 64*(k+1)*d^2+64*k*(m+1)*d+512*d, with full
    /// dimension d and declared force/observation ceilings k,m; not FLOPs.
    pub max_physical_work: u64,
}

#[derive(Debug)]
pub enum DgPoleResponseError {
    Invalid(&'static str),
    Physical(DgLoadedPoleError),
    Numerical(PoleResponseError),
    Cancelled,
}
impl core::fmt::Display for DgPoleResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete modal/DG pole response refused: {self:?}")
    }
}
impl std::error::Error for DgPoleResponseError {}
impl From<DgLoadedPoleError> for DgPoleResponseError {
    fn from(error: DgLoadedPoleError) -> Self {
        Self::Physical(error)
    }
}
impl From<PoleResponseError> for DgPoleResponseError {
    fn from(error: PoleResponseError) -> Self {
        Self::Numerical(error)
    }
}

/// Laurent coefficients X(s)=state/(s-s0)+regular for one original modal force.
/// All DG, branch, Maxwell, radiation, terminal, gap and algebraic coordinates
/// remain in [`DgPoleCoordinates`] order, with their original physical units.
/// These are not finite harmonic amplitudes, stored energies or mean powers.
/// Under exp(-i*omega*t), omega0=i*s0 and omega-chart residues are i times
/// these Laplace coefficients. A zero residue does not remove the free mode.
#[derive(Clone, Debug)]
pub struct DgLoadedModalForceResidue {
    pub laplace_pole_per_s: C64,
    pub state: Vec<C64>,
    /// s0 times each original mechanical displacement residue.
    pub modal_velocity_m_sqrt_kg_per_s: Vec<C64>,
    pub radiation_velocity: Vec<C64>,
    pub inlet_flow_m3_s: C64,
    pub inlet_pressure_pa: C64,
    /// Bilinear weight^T*state coefficients, in supplied observation order.
    pub observations: Vec<C64>,
    /// state = prepared.right_mode().state * participation in the native
    /// pivot-one gauge; it includes nonnormal left-mode forcing participation.
    pub participation: C64,
    /// Shared numerical owner's original bordered-equation backward residual.
    pub linear_residual: f64,
    /// Independent complete original polynomial-row residual, exactly zero for
    /// a zero residue. This is not a forward shape or transfer-error bound.
    pub original_relative_residual: f64,
}

/// One immutable, reusable factor bound to a complete prepared physical pencil.
/// Only an isolated finite simple pole is admitted. This does not replace the
/// complete transfer function, its regular/feedthrough terms, or repeated poles.
pub struct PreparedDgPoleResponse<'prepared, 'source> {
    source: &'prepared PreparedModalDgPoles<'source>,
    numerical: PreparedPole,
    mode: DgLoadedModalPole,
    options: DgPoleResponseOptions,
}

fn check_cancel(gate: &CancelGate) -> Result<(), DgPoleResponseError> {
    if gate.is_requested() {
        Err(DgPoleResponseError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(d: usize, options: &DgPoleResponseOptions) -> Result<(), DgPoleResponseError> {
    if d > options.numerical.max_dimension
        || !(1..=256).contains(&options.numerical.max_vectors)
        || options.max_observations > 256
    {
        return Err(DgPoleResponseError::Invalid(
            "complete dimension and bounded force/observation counts required",
        ));
    }
    let (d, k, m) = (
        d as u128,
        options.numerical.max_vectors as u128,
        options.max_observations as u128,
    );
    let bytes = 16 * (2 * d * d + k * (5 * d + m) + 12 * d) + 512 * (k + m);
    let work = 64 * (k + 1) * d * d + 64 * k * (m + 1) * d + 512 * d;
    if bytes > options.max_physical_bytes as u128 || work > options.max_physical_work as u128 {
        return Err(DgPoleResponseError::Invalid(
            "complete pole forcing/observation mapping exceeds its physical budget",
        ));
    }
    Ok(())
}

impl<'source> PreparedModalDgPoles<'source> {
    /// Prepare forcing residues and native parameter tangents at a resolved pole.
    /// Cached pole observations/residuals are not trusted: the complete supplied
    /// state is rechecked against THIS pencil, its source bands and original
    /// split coefficient rows. The native pivot-one normalization is retained.
    ///
    /// The shared numerical owner factors one bordered system and obtains true
    /// left/right participation without a separate adjoint solve. Fixed row
    /// equilibration is applied to both the pencil and later mechanical forces;
    /// the physical Laplace variable and all original coordinates stay unchanged.
    /// Source band admission does not certify empirical complex continuation.
    pub fn prepare_pole_response<'prepared>(
        &'prepared self,
        pole: &DgLoadedModalPole,
        options: DgPoleResponseOptions,
        gate: &CancelGate,
    ) -> Result<PreparedDgPoleResponse<'prepared, 'source>, DgPoleResponseError> {
        check_cancel(gate)?;
        admit(self.dimension(), &options)?;
        self.admit_oscillation(pole.laplace_pole_per_s, gate)?;
        if pole.state.len() != self.dimension()
            || pole.state.iter().any(|&value| !complex_finite(value))
        {
            return Err(DgPoleResponseError::Invalid(
                "complete finite original pole state required",
            ));
        }
        let evaluation = self
            .pencil
            .evaluate(pole.laplace_pole_per_s, &mut || gate.is_requested())
            .map_err(DgLoadedPoleError::from)?;
        let numerical = PreparedPole::new(
            pole.laplace_pole_per_s,
            &pole.state,
            &evaluation,
            options.numerical,
            || gate.is_requested(),
        )?;
        let mode = self.observe_pole(
            NonlinearPole {
                value: numerical.value(),
                vector: numerical.right_vector().to_vec(),
                residual: numerical.pole_residual(),
                relative_step: numerical.relative_correction(),
                residual_scale: evaluation.residual_scale,
                iterations: 0,
            },
            options.numerical.pole_residual_tolerance,
            gate,
        )?;
        check_cancel(gate)?;
        Ok(PreparedDgPoleResponse {
            source: self,
            numerical,
            mode,
            options,
        })
    }
}

impl PreparedDgPoleResponse<'_, '_> {
    /// The rechecked physical right mode in the native pivot-one gauge.
    pub fn right_mode(&self) -> &DgLoadedModalPole {
        &self.mode
    }
    pub fn coordinates(&self) -> &DgPoleCoordinates {
        self.source.coordinates()
    }
    pub fn normalization_pivot(&self) -> usize {
        self.numerical.normalization_pivot()
    }

    /// Laurent residues for complete original mechanical-force batches evaluated
    /// at s0, in mass-normalized mode order. A frequency-dependent analytic drive
    /// contributes its value at the pole; its derivative changes regular terms.
    /// No inlet flow, internal state source or terminal drive is inferred.
    ///
    /// Each optional observation is a COMPLETE complex weight vector over the
    /// original mixed-unit state, contracted bilinearly without conjugation or
    /// normalization. Supply source-owned pressure/field columns or analytic
    /// sensor weights at s0; their units and frequency law remain caller-owned.
    /// Mechanical/radiation velocity and inlet pressure/flow are mapped directly.
    ///
    /// One cached border serves every batch. Each RHS costs O(d^2), plus the
    /// original polynomial audit and O(d) per sensor. Exact zero and dark-mode
    /// residues are valid. No state is dropped and no energy/Mach interpretation
    /// is assigned to a residue. Refusal/cancellation publishes no partial batch
    /// and leaves the prepared factor and all physical histories reusable.
    pub fn forcing_residues(
        &self,
        forces: &[&[C64]],
        observations: &[&[C64]],
        gate: &CancelGate,
    ) -> Result<Vec<DgLoadedModalForceResidue>, DgPoleResponseError> {
        check_cancel(gate)?;
        let (d, modes) = (self.source.dimension(), self.source.source.mode_count());
        if forces.is_empty()
            || forces.len() > self.options.numerical.max_vectors
            || observations.len() > self.options.max_observations
            || forces
                .iter()
                .any(|force| force.len() != modes || force.iter().any(|&v| !complex_finite(v)))
            || observations
                .iter()
                .any(|weights| weights.len() != d || weights.iter().any(|&v| !complex_finite(v)))
        {
            return Err(DgPoleResponseError::Invalid(
                "complete finite modal forces and raw-state observations within declared ceilings required",
            ));
        }
        let mut mapped = Vec::with_capacity(forces.len());
        for &force in forces {
            check_cancel(gate)?;
            let mut rhs = vec![C64::ZERO; d];
            rhs[self.source.coordinates.modal_displacement.clone()].copy_from_slice(force);
            self.source
                .pencil
                .scale_rhs(&mut rhs, &mut || gate.is_requested())
                .map_err(DgLoadedPoleError::from)?;
            mapped.push(rhs);
        }
        let refs: Vec<_> = mapped.iter().map(Vec::as_slice).collect();
        let numerical = self
            .numerical
            .forcing_residues(&refs, || gate.is_requested())?;
        let tolerance = self
            .options
            .numerical
            .pole_residual_tolerance
            .min(self.source.source.config.solve_relative_tolerance);
        let mut output = Vec::with_capacity(numerical.len());
        for residue in numerical {
            check_cancel(gate)?;
            // The eigenvector audit requires a nonzero mode; a zero forced
            // Laurent coefficient is valid and satisfies every homogeneous row.
            let relative = if residue.vector.iter().all(|&value| value == C64::ZERO) {
                0.0
            } else {
                self.source
                    .pencil
                    .relative_residual(self.mode.laplace_pole_per_s, &residue.vector, &mut || {
                        gate.is_requested()
                    })
                    .map_err(DgLoadedPoleError::from)?
            };
            if relative > tolerance {
                return Err(DgLoadedPoleError::Residual {
                    relative,
                    tolerance,
                }
                .into());
            }
            let velocity: Vec<_> = self
                .mode
                .modal_velocity_m_sqrt_kg_per_s
                .iter()
                .map(|&value| value * residue.coefficient)
                .collect();
            let radiation: Vec<_> = self
                .mode
                .radiation_velocity
                .iter()
                .map(|&value| value * residue.coefficient)
                .collect();
            let flow = self.mode.inlet_flow_m3_s * residue.coefficient;
            let pressure = self.mode.inlet_pressure_pa * residue.coefficient;
            if velocity
                .iter()
                .chain(&radiation)
                .any(|&value| !complex_finite(value))
                || !complex_finite(flow)
                || !complex_finite(pressure)
            {
                return Err(DgPoleResponseError::Invalid(
                    "physical pole residue observation overflows",
                ));
            }
            let mut observed = Vec::with_capacity(observations.len());
            for &weights in observations {
                observed.push(contract(weights, &residue.vector, gate)?);
            }
            output.push(DgLoadedModalForceResidue {
                laplace_pole_per_s: self.mode.laplace_pole_per_s,
                state: residue.vector,
                modal_velocity_m_sqrt_kg_per_s: velocity,
                radiation_velocity: radiation,
                inlet_flow_m3_s: flow,
                inlet_pressure_pa: pressure,
                observations: observed,
                participation: residue.coefficient,
                linear_residual: residue.linear_residual,
                original_relative_residual: relative,
            });
        }
        check_cancel(gate)?;
        Ok(output)
    }
}

fn contract(weights: &[C64], state: &[C64], gate: &CancelGate) -> Result<C64, DgPoleResponseError> {
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for (i, (&w, &x)) in weights.iter().zip(state).enumerate() {
        if i % 64 == 0 {
            check_cancel(gate)?;
        }
        let (wr, wi, xr, xi) = (
            Dd::from_f64(w.re),
            Dd::from_f64(w.im),
            Dd::from_f64(x.re),
            Dd::from_f64(x.im),
        );
        re = re + wr * xr - wi * xi;
        im = im + wr * xi + wi * xr;
    }
    let result = C64::new(re.to_f64(), im.to_f64());
    if !complex_finite(result) {
        return Err(DgPoleResponseError::Invalid(
            "bilinear pole observation is not representable",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
