//! Continuous harmonic loads on every admitted physical solid coordinate.

use super::{
    CompleteSpectralSolidControls, DampedSpectralSolid, ExactLinearSolidStep,
    LinearSolidImpulseError, LinearSolidImpulseReport, Result, TetElasticAssembly, add, bound,
    invalid, poll, project, zeros,
};
use fs_exec::Cx;
use fs_modal::time::{ModalAcousticState, PreparedHarmonicForce};
use fs_sparse::Csr;

/// Complete solid motion under one continuous physical forcing frequency.
///
/// Every free coordinate has its own cosine and sine force amplitudes in N.
/// Their projection uses the original complete basis, including rigid modes.
/// The shared scalar owner retains the initial transient and resonant growth;
/// it does not replace forcing by samples or a steady-state approximation.
/// Viscosity obeys the same classical-law admission as `DampedSpectralSolid`.
/// Independent force work and positive viscous heat enter the same original
/// M/K and cumulative energy gates before state or source phase can change.
///
/// Frequencies are fixed at preparation. This is one simultaneous tone, with
/// arbitrary spatial phase; separate-tone heats cannot be added to represent
/// a multi-frequency load. There is no acoustic or sparse-scale claim.
pub struct HarmonicDampedSpectralSolid {
    solid: DampedSpectralSolid,
    transitions: Vec<PreparedHarmonicForce>,
    forcing_angular_frequency_rad_s: f64,
}

impl HarmonicDampedSpectralSolid {
    /// Prepare exact continuous force motion without a large augmented matrix.
    ///
    /// Supply every eigenpair and the original viscosity in the same order as
    /// `DampedSpectralSolid::prepare_with_basis`. The forcing frequency is
    /// finite and nonnegative; zero denotes DC and reuses held-force arithmetic.
    /// No Nyquist cutoff or artificial damping changes a prescribed tone.
    /// The caller still resolves any intermediate observations it requires.
    ///
    /// For nonzero forcing frequency, preparation admits `96*n^2 + 4096*n`
    /// bytes and `6*n^3 + 64*n^2 + 32768*n + 8192*sum(s_i)` product terms.
    /// Each `s_i` halves `max(omega_i*h, c_i*h, Omega*h, 1)` to at most 1/8.
    /// These bounds include the underlying damped preparation. Warm intervals
    /// keep its `16*n^2 + 512*n` bound, including continuous-source rotation.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_basis(
        assembly: &TetElasticAssembly,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        forcing_angular_frequency_rad_s: f64,
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: CompleteSpectralSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        if !forcing_angular_frequency_rad_s.is_finite() || forcing_angular_frequency_rad_s < 0.0 {
            return Err(invalid("finite nonnegative forcing frequency required"));
        }
        let n = assembly.free_dofs.len();
        let nn = n as u128;
        let phase = forcing_angular_frequency_rad_s * dt_s;
        let mut preparation_terms = 6 * nn * nn * nn + 64 * nn * nn + 32768 * nn;
        if forcing_angular_frequency_rad_s != 0.0 {
            if !phase.is_finite() {
                return Err(invalid("representable harmonic forcing phase required"));
            }
            bound(n as u128, controls.maximum_free_dofs as u128)?;
            bound(
                96 * nn * nn + 4096 * nn,
                controls.maximum_preparation_bytes as u128,
            )?;
            bound(
                preparation_terms,
                controls.maximum_preparation_multiply_terms as u128,
            )?;
        }
        let solid = DampedSpectralSolid::prepare_with_basis(
            assembly,
            viscosity_kg_per_s,
            eigenvalues_rad2_s2,
            basis_columns,
            dt_s,
            initial_displacement_m,
            initial_velocity_mps,
            controls,
            cx,
        )?;
        let mut transitions = Vec::new();
        if forcing_angular_frequency_rad_s != 0.0 {
            transitions
                .try_reserve_exact(n)
                .map_err(|_| invalid("harmonic coefficient allocation refused"))?;
            for (i, &lambda) in eigenvalues_rad2_s2.iter().enumerate() {
                poll(cx)?;
                let omega = lambda.sqrt();
                let rate = solid.damping_rates_per_s()[i];
                let mut scale = (omega * dt_s).max(rate * dt_s).max(phase).max(1.0);
                if !scale.is_finite() {
                    return Err(invalid("representable harmonic oscillator scale required"));
                }
                while scale > 0.125 {
                    scale *= 0.5;
                    preparation_terms += 8192;
                }
                bound(
                    preparation_terms,
                    controls.maximum_preparation_multiply_terms as u128,
                )?;
                transitions.push(
                    PreparedHarmonicForce::try_new_with_damping_rate(
                        omega,
                        rate,
                        forcing_angular_frequency_rad_s,
                        dt_s,
                    )
                    .map_err(|_| {
                        invalid("unrepresentable harmonic motion or power coefficients")
                    })?,
                );
            }
        }
        poll(cx)?;
        Ok(Self {
            solid,
            transitions,
            forcing_angular_frequency_rad_s,
        })
    }

    pub fn forcing_angular_frequency_rad_s(&self) -> f64 {
        self.forcing_angular_frequency_rad_s
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
        self.solid.cumulative_dissipated_j()
    }
    pub fn damping_rates_per_s(&self) -> &[f64] {
        self.solid.damping_rates_per_s()
    }
    pub fn relative_damping_projection_error(&self) -> f64 {
        self.solid.relative_damping_projection_error()
    }

    /// Original-coordinate instantaneous impulse, retaining the earlier heat
    /// and unchanged source phase. It consumes no harmonic interval.
    pub fn apply_momentum_impulse(
        &mut self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> core::result::Result<LinearSolidImpulseReport, LinearSolidImpulseError> {
        self.solid
            .apply_momentum_impulse(impulse_kg_m_s, maximum_impulse_error_kg_m_s, cx)
    }

    /// Pure next-interval source rotation: `(Fc*cos + Fs*sin, Fs*cos - Fc*sin)`.
    /// Use only after accepting an interval, or use `advance_continuous`.
    /// Repeated rotations accumulate floating-point phase/amplitude roundoff;
    /// callers may instead form coefficients from their own absolute clock.
    /// This method adds no normalization, damping, or accepted-state change.
    pub fn advance_force_phase(&self, cosine_n: f64, sine_n: f64) -> Result<(f64, f64)> {
        if !cosine_n.is_finite() || !sine_n.is_finite() {
            return Err(invalid("finite harmonic force quadratures required"));
        }
        if self.forcing_angular_frequency_rad_s == 0.0 {
            return Ok((cosine_n, sine_n));
        }
        self.transitions[0]
            .motion()
            .advance_force_phase(cosine_n, sine_n)
            .map_err(|_| invalid("unrepresentable rotated physical force"))
    }

    /// Advance `F(tau)=Fc*cos(Omega*tau)+Fs*sin(Omega*tau)` on this interval.
    /// Both physical force arrays use `free_dofs()` order and an interval-local
    /// time origin. Changed coefficients describe a finite source jump, not an
    /// impulse. For DC only the cosine coefficient contributes physical force.
    ///
    /// Work and heat include the initial transient and both force quadratures.
    /// Refusal preserves every accepted state, clock and energy account.
    pub fn advance(
        &mut self,
        cosine_force_n: &[f64],
        sine_force_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.admit_forces(cosine_force_n, sine_force_n, cx)?;
        if self.forcing_angular_frequency_rad_s == 0.0 {
            return self.solid.advance(cosine_force_n, cx);
        }
        let (step, time_s) = self.solid.solid.next_interval()?;
        let cosine = project(&self.solid.solid.basis, cosine_force_n, cx)?;
        let sine = project(&self.solid.solid.basis, sine_force_n, cx)?;
        let n = self.free_dofs().len();
        let (mut q, mut v) = (zeros(n)?, zeros(n)?);
        let (mut work, mut work_correction) = (0.0, 0.0);
        let (mut heat, mut heat_correction) = (0.0, 0.0);
        for (i, transition) in self.transitions.iter().enumerate() {
            poll(cx)?;
            let frame = transition
                .advance_with_energy(
                    ModalAcousticState {
                        displacement_m_sqrt_kg: self.solid.solid.modal_displacement[i],
                        velocity_m_sqrt_kg_per_s: self.solid.solid.modal_velocity[i],
                    },
                    cosine[i],
                    sine[i],
                )
                .map_err(|_| invalid("unrepresentable harmonic motion, work or heat"))?;
            q[i] = frame.state.displacement_m_sqrt_kg;
            v[i] = frame.state.velocity_m_sqrt_kg_per_s;
            add(&mut work, &mut work_correction, frame.input_work_j);
            add(&mut heat, &mut heat_correction, frame.viscous_dissipation_j);
        }
        self.solid
            .accept_modal_path(q, v, work, heat, step, time_s, cx)
    }

    /// Continue the same physical tone and rotate its caller-owned quadratures
    /// only after accepting motion. Cancellation, phase overflow, work failure
    /// or any other refusal leaves BOTH arrays and the complete solid unchanged.
    /// Checkpoint these arrays together with the physical state for replay.
    pub fn advance_continuous(
        &mut self,
        cosine_force_n: &mut [f64],
        sine_force_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.admit_forces(cosine_force_n, sine_force_n, cx)?;
        let n = self.free_dofs().len();
        let (mut cosine_next, mut sine_next) = (zeros(n)?, zeros(n)?);
        for i in 0..n {
            poll(cx)?;
            (cosine_next[i], sine_next[i]) =
                self.advance_force_phase(cosine_force_n[i], sine_force_n[i])?;
        }
        let result = self.advance(cosine_force_n, sine_force_n, cx)?;
        cosine_force_n.copy_from_slice(&cosine_next);
        sine_force_n.copy_from_slice(&sine_next);
        Ok(result)
    }

    fn admit_forces(&self, cosine: &[f64], sine: &[f64], cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        let n = self.free_dofs().len();
        if cosine.len() != n
            || sine.len() != n
            || cosine.iter().chain(sine).any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite physical force quadratures required",
            ));
        }
        Ok(())
    }
}
