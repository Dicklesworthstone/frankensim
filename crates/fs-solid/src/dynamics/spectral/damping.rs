//! Complete-spectrum motion with an admitted classical viscous law.

use super::{
    CompleteSpectralSolid, CompleteSpectralSolidControls, ExactLinearSolidStep,
    LinearSolidImpulseError, LinearSolidImpulseReport, Result, add, bound, check_operator,
    invalid, poll, project, reconstruct, symmetric, zeros,
};
use crate::linear3::TetElasticAssembly;
use fs_exec::Cx;
use fs_modal::time::{LinearForceTransition, ModalAcousticState};
use fs_sparse::Csr;

/// Every original elastic and rigid coordinate, with physical viscous damping.
///
/// The supplied symmetric viscosity C must be diagonal in the complete,
/// mass-orthonormal elastic basis to the explicit operator allowance. Each
/// transformed diagonal must be nonnegative. Off-diagonal row sums must be no
/// larger than `maximum_relative_operator_error` times that row's diagonal;
/// a zero diagonal therefore requires an exactly zero row. Meaningful coupled
/// damping refuses rather than being discarded. Use `ExactLinearSolid` for
/// general coupled viscosity.
///
/// Prepared scalar transitions integrate held or affine forces through
/// underdamped, critical, overdamped and viscously damped rigid motion. Heat
/// is an independent positive velocity integral, including losses below the
/// roundoff of stored-energy subtraction. The original M/K still audit both
/// interval and cumulative energy before publication. A zero C uses the
/// conservative owner's arithmetic unchanged.
///
/// This separate type exposes no conservative cubic-force, moving-support or
/// sustained-contact operation that could accidentally bypass its viscosity.
pub struct DampedSpectralSolid {
    solid: CompleteSpectralSolid,
    transitions: Option<Vec<LinearForceTransition>>,
    damping_rates_per_s: Vec<f64>,
    relative_damping_projection_error: f64,
    cumulative_dissipated_j: f64,
    dissipation_correction_j: f64,
}

impl DampedSpectralSolid {
    /// Admit the original pencil, complete basis and physical viscosity.
    /// All state, force, impulse and viscosity rows use `assembly.free_dofs`.
    ///
    /// The existing controls also bound the damping preparation. Its byte
    /// envelope is `96*n^2 + 2048*n`; scalar-product admission is
    /// `6*n^3 + 64*n^2 + 8192*n + 512*sum(s_i)`, where `s_i` resolves each
    /// dimensionless frequency/drag to at most 1/4 by halving (starting at a
    /// scale of at least one). Warm intervals retain the existing
    /// `16*n^2 + 512*n` bound. No eigensolve or augmented exponential is hidden.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_basis(
        assembly: &TetElasticAssembly,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: CompleteSpectralSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let n = assembly.free_dofs.len();
        let nn = n as u128;
        bound(n as u128, controls.maximum_free_dofs as u128)?;
        bound(
            96 * nn * nn + 2048 * nn,
            controls.maximum_preparation_bytes as u128,
        )?;
        let mut preparation_terms = 6 * nn * nn * nn + 64 * nn * nn + 8192 * nn;
        bound(preparation_terms, controls.maximum_preparation_multiply_terms as u128)?;
        let solid = CompleteSpectralSolid::prepare_with_basis(
            assembly,
            eigenvalues_rad2_s2,
            basis_columns,
            dt_s,
            initial_displacement_m,
            initial_velocity_mps,
            controls,
            cx,
        )?;
        let viscosity = symmetric(
            viscosity_kg_per_s,
            n,
            controls.maximum_relative_operator_error,
            cx,
        )?;
        let mut rates = zeros(n)?;
        let mut relative_error = 0.0_f64;
        let transitions = if viscosity.iter().all(|entry| *entry == 0.0) {
            None
        } else {
            let mut c_phi = zeros(n * n)?;
            let mut modal = zeros(n * n)?;
            for i in 0..n {
                for j in 0..n {
                    poll(cx)?;
                    c_phi[i * n + j] = (0..n).fold(0.0, |sum, k| {
                        viscosity[i * n + k].mul_add(solid.basis[k * n + j], sum)
                    });
                    if !c_phi[i * n + j].is_finite() {
                        return Err(invalid("unrepresentable viscosity projection"));
                    }
                }
            }
            for i in 0..n {
                for j in 0..n {
                    poll(cx)?;
                    modal[i * n + j] = (0..n).fold(0.0, |sum, k| {
                        solid.basis[k * n + i].mul_add(c_phi[k * n + j], sum)
                    });
                    if !modal[i * n + j].is_finite() {
                        return Err(invalid("unrepresentable modal viscosity"));
                    }
                }
                rates[i] = modal[i * n + i];
                if rates[i] < 0.0 {
                    return Err(invalid("nonnegative modal viscosity required"));
                }
            }
            for i in 0..n {
                poll(cx)?;
                let (mut off_diagonal, mut correction) = (0.0, 0.0);
                for j in 0..n {
                    if i != j {
                        // Only the symmetric part contributes physical power.
                        // C was admitted symmetric in original coordinates.
                        let value = 0.5 * modal[i * n + j] + 0.5 * modal[j * n + i];
                        add(&mut off_diagonal, &mut correction, value.abs());
                    }
                }
                check_operator(
                    off_diagonal,
                    controls.maximum_relative_operator_error * rates[i],
                )?;
                if rates[i] != 0.0 {
                    relative_error = relative_error.max(off_diagonal / rates[i]);
                }
            }
            if rates.iter().all(|rate| *rate == 0.0) {
                return Err(invalid("nonzero viscosity projection underflowed"));
            }
            let mut prepared = Vec::new();
            prepared
                .try_reserve_exact(n)
                .map_err(|_| invalid("damped coefficient allocation refused"))?;
            for i in 0..n {
                poll(cx)?;
                let omega = eigenvalues_rad2_s2[i].sqrt();
                let angle = omega * dt_s;
                let drag = rates[i] * dt_s;
                if !angle.is_finite() || !drag.is_finite() {
                    return Err(invalid("unrepresentable damped oscillator phase or drag"));
                }
                let mut scale = angle.max(drag).max(1.0);
                while scale > 0.25 {
                    scale *= 0.5;
                    preparation_terms += 512;
                }
                bound(preparation_terms, controls.maximum_preparation_multiply_terms as u128)?;
                prepared.push(
                    LinearForceTransition::try_new_with_damping_rate(omega, rates[i], dt_s)
                        .map_err(|_| invalid("unrepresentable damped oscillator transition"))?,
                );
            }
            Some(prepared)
        };
        poll(cx)?;
        Ok(Self {
            solid,
            transitions,
            damping_rates_per_s: rates,
            relative_damping_projection_error: relative_error,
            cumulative_dissipated_j: 0.0,
            dissipation_correction_j: 0.0,
        })
    }

    pub fn displacement_m(&self) -> &[f64] {
        self.solid.displacement_m()
    }
    pub fn velocity_mps(&self) -> &[f64] {
        self.solid.velocity_mps()
    }
    pub fn free_dofs(&self) -> &[usize] {
        self.solid.free_dofs()
    }
    pub fn steps(&self) -> u64 {
        self.solid.steps()
    }
    pub fn time_s(&self) -> f64 {
        self.solid.time_s()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.solid.stored_energy_j()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.solid.kinetic_energy_j()
    }
    pub fn cumulative_work_j(&self) -> f64 {
        self.solid.cumulative_work_j()
    }
    pub fn cumulative_dissipated_j(&self) -> f64 {
        self.cumulative_dissipated_j
    }
    /// Actual admitted diagonal rates, including drag on zero-stiffness modes.
    pub fn damping_rates_per_s(&self) -> &[f64] {
        &self.damping_rates_per_s
    }
    /// Maximum off-diagonal absolute row sum divided by that row's diagonal
    /// in the computed symmetric Phi^T C Phi. Zero rows contribute zero.
    ///
    /// For the represented trajectory this bounds omitted cross-term work by
    /// this value times its reported viscous heat, before matrix-projection
    /// roundoff. It is a numerical admission diagnostic, not an interval proof.
    pub fn relative_damping_projection_error(&self) -> f64 {
        self.relative_damping_projection_error
    }

    /// Instantaneous original-coordinate impulse; prior physical heat remains
    /// in the cumulative audit. The event itself has no elapsed time or heat.
    pub fn apply_momentum_impulse(
        &mut self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> core::result::Result<LinearSolidImpulseReport, LinearSolidImpulseError> {
        let candidate = self.solid.trial_momentum_impulse_with_dissipation(
            impulse_kg_m_s,
            maximum_impulse_error_kg_m_s,
            self.cumulative_dissipated_j,
            cx,
        )?;
        Ok(self.solid.accept_momentum_impulse(candidate))
    }

    /// Exact held physical force, retaining every complete-basis mode.
    pub fn advance(&mut self, force_n: &[f64], cx: &Cx<'_>) -> Result<ExactLinearSolidStep> {
        self.advance_affine(force_n, force_n, cx)
    }

    /// Exact linear force interpolation with independent work and viscous heat.
    /// A force jump across intervals creates no impulse. Refusals publish none
    /// of the candidate physical/modal state, clock, work or dissipation.
    pub fn advance_affine(
        &mut self,
        force_start_n: &[f64],
        force_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        poll(cx)?;
        let Some(transitions) = &self.transitions else {
            return self.solid.advance_affine(force_start_n, force_end_n, cx);
        };
        let n = self.free_dofs().len();
        if force_start_n.len() != n
            || force_end_n.len() != n
            || force_start_n.iter().chain(force_end_n).any(|value| !value.is_finite())
        {
            return Err(invalid("complete finite force endpoints required"));
        }
        let (step, time_s) = self.solid.next_interval()?;
        let force_start = project(&self.solid.basis, force_start_n, cx)?;
        let force_end = project(&self.solid.basis, force_end_n, cx)?;
        let (mut q, mut v) = (zeros(n)?, zeros(n)?);
        let (mut work, mut work_correction) = (0.0, 0.0);
        let (mut heat, mut heat_correction) = (0.0, 0.0);
        for i in 0..n {
            poll(cx)?;
            let frame = transitions[i]
                .advance_with_energy(
                    ModalAcousticState {
                        displacement_m_sqrt_kg: self.solid.modal_displacement[i],
                        velocity_m_sqrt_kg_per_s: self.solid.modal_velocity[i],
                    },
                    force_start[i],
                    force_end[i],
                )
                .map_err(|_| invalid("unrepresentable damped motion or work"))?;
            q[i] = frame.state.displacement_m_sqrt_kg;
            v[i] = frame.state.velocity_m_sqrt_kg_per_s;
            add(&mut work, &mut work_correction, frame.input_work_j);
            add(&mut heat, &mut heat_correction, frame.viscous_dissipation_j);
        }
        let displacement = reconstruct(&self.solid.basis, &q, cx)?;
        let velocity = reconstruct(&self.solid.basis, &v, cx)?;
        let (mut cumulative_heat, mut correction) =
            (self.cumulative_dissipated_j, self.dissipation_correction_j);
        add(&mut cumulative_heat, &mut correction, heat);
        if ![work, heat, cumulative_heat, correction].into_iter().all(f64::is_finite) {
            return Err(invalid("unrepresentable accumulated damped work or heat"));
        }
        let result = self.solid.accept_force_path_with_dissipation(
            q, v, displacement, velocity, work, heat, cumulative_heat, step, time_s, cx,
        )?;
        self.cumulative_dissipated_j = cumulative_heat;
        self.dissipation_correction_j = correction;
        Ok(result)
    }
}
