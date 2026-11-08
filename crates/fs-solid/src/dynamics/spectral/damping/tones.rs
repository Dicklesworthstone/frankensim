//! Simultaneous continuous tones and affine envelopes on the complete solid.

use super::{
    CompleteSpectralSolidControls, DampedSpectralSolid, ExactLinearSolidStep,
    LinearSolidImpulseError, LinearSolidImpulseReport, Result, SpectralParameterUpdate,
    TetElasticAssembly, add, bound, invalid, poll, zeros,
};
use fs_exec::Cx;
use fs_modal::time::{
    MAX_HARMONIC_TONES, ModalAcousticState, ModalAcousticTimeError, ModalToneBudget,
    PreparedModalTones,
};
use fs_sparse::Csr;

/// Every original mode under simultaneous physical force tones.
///
/// Forces are DOF-major: entry `i*K+j` is tone j on `solid().free_dofs()[i]`,
/// in N. Exactly repeated frequencies coalesce in the shared forcing owner;
/// close frequencies remain distinct and retain beats. Each mode has ONE
/// total state and complete cross-frequency work/heat, including its initial
/// transient. Summing separate-tone dissipations is never used.
///
/// The original complete M/K and classical viscosity admission are unchanged.
/// No acoustic model, sample Nyquist cutoff or artificial damping is introduced.
/// Exact endpoint integration does not resolve continuous peaks or remove
/// aliasing from observations sampled too slowly.
pub struct MultiToneDampedSpectralSolid {
    solid: DampedSpectralSolid,
    tones: PreparedModalTones,
    force_count: usize,
}

impl MultiToneDampedSpectralSolid {
    /// Bind all original modes and a fixed list of physical source frequencies.
    ///
    /// With `affine_envelopes=false`, quadratures are constant on each interval.
    /// With `true`, their start/end values interpolate linearly in the SAME
    /// interval-start carrier basis. DC is an exact constant/affine force.
    /// Both paths include positive drag on an exactly rigid coordinate.
    ///
    /// The overall solid byte cap reserves `96*n^2+4096*n+32*n*K` for
    /// physical matrices and projection scratch, then clamps the native tone
    /// storage budget to the remaining bytes. The caller's smaller tone-storage
    /// limit still applies. Work caps are separate: `controls` covers original
    /// preparation and `16*n^2+512*n+P*n^2*K` warm projection terms (P=2 or 4);
    /// `tone_budget` additionally bounds ALL modes' native preparation/step
    /// work scores. Neither work budget is reset per mode.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_basis(
        assembly: &TetElasticAssembly,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        forcing_angular_frequencies_rad_s: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        affine_envelopes: bool,
        controls: CompleteSpectralSolidControls,
        mut tone_budget: ModalToneBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let n = assembly.free_dofs.len();
        let k = forcing_angular_frequencies_rad_s.len();
        if k == 0
            || k > MAX_HARMONIC_TONES
            || k > tone_budget.maximum_tones
            || forcing_angular_frequencies_rad_s
                .iter()
                .any(|w| !w.is_finite() || *w < 0.0)
        {
            return Err(invalid(
                "bounded nonempty finite nonnegative tone frequencies required",
            ));
        }
        let force_count = n
            .checked_mul(k)
            .ok_or_else(|| invalid("tone force shape overflow"))?;
        let (nn, kk) = (n as u128, k as u128);
        bound(nn, controls.maximum_free_dofs as u128)?;
        let physical_bytes = 96 * nn * nn + 4096 * nn + 32 * nn * kk;
        bound(physical_bytes, controls.maximum_preparation_bytes as u128)?;
        tone_budget.maximum_storage_bytes = tone_budget
            .maximum_storage_bytes
            .min((controls.maximum_preparation_bytes as u128 - physical_bytes) as usize);
        let columns = if affine_envelopes { 4 } else { 2 };
        bound(
            16 * nn * nn + 512 * nn + columns * nn * nn * kk,
            controls.maximum_step_multiply_terms as u128,
        )?;
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
        Self::prepare_source(
            solid,
            forcing_angular_frequencies_rad_s,
            affine_envelopes,
            tone_budget,
            force_count,
            cx,
        )
    }

    fn prepare_source(
        solid: DampedSpectralSolid,
        forcing_angular_frequencies_rad_s: &[f64],
        affine_envelopes: bool,
        tone_budget: ModalToneBudget,
        force_count: usize,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        let dt_s = solid.solid.dt_s;
        let mut omega = zeros(solid.free_dofs().len())?;
        for (target, oscillator) in omega.iter_mut().zip(&solid.solid.oscillators) {
            poll(cx)?;
            *target = oscillator.lambda.sqrt();
        }
        let prepared = PreparedModalTones::try_new_with_damping_rates_controlled(
            &omega,
            solid.damping_rates_per_s(),
            forcing_angular_frequencies_rad_s,
            dt_s,
            tone_budget,
            affine_envelopes,
            || cx.checkpoint().is_err(),
        );
        let tones = tone_result(prepared, cx)?;
        poll(cx)?;
        Ok(Self {
            solid,
            tones,
            force_count,
        })
    }

    /// Rebuild all physical modes and their combined source-power kernels.
    ///
    /// This has the fixed-mass and signed-work meaning of
    /// `DampedSpectralSolid::replace_parameters`. Original tone order,
    /// frequencies, duration and constant/affine envelope mode remain fixed.
    /// Caller-owned quadratures are not touched or rotated by this event.
    /// Both the complete physical candidate and its exact multi-tone caches
    /// must pass before either is published.
    ///
    /// The global solid byte cap reserves `128*n^2+8192*n+32*n*K`, PLUS the
    /// old bank's admitted payload/scratch envelope, before clamping the new
    /// `tone_budget.maximum_storage_bytes` to what remains. Thus old and new
    /// banks are covered simultaneously. Native preparation/step work remains
    /// an additional aggregate `tone_budget`; no allowance resets per mode.
    #[allow(clippy::too_many_arguments)]
    pub fn replace_parameters(
        &mut self,
        stiffness_n_per_m: &Csr,
        viscosity_kg_per_s: &Csr,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        maximum_parameter_work_j: f64,
        mut tone_budget: ModalToneBudget,
        cx: &Cx<'_>,
    ) -> Result<SpectralParameterUpdate> {
        poll(cx)?;
        let (n, k) = (
            self.solid.free_dofs().len() as u128,
            self.tones.tone_count(),
        );
        if k > tone_budget.maximum_tones {
            return Err(invalid("replacement tone budget excludes retained sources"));
        }
        let controls = self.solid.solid.controls;
        let reserved = 128 * n * n
            + 8192 * n
            + 32 * n * k as u128
            + self.tones.preparation_storage_bound_bytes();
        bound(reserved, controls.maximum_preparation_bytes as u128)?;
        tone_budget.maximum_storage_bytes = tone_budget
            .maximum_storage_bytes
            .min((controls.maximum_preparation_bytes as u128 - reserved) as usize);
        let (solid, report) = self.solid.parameter_candidate(
            stiffness_n_per_m,
            viscosity_kg_per_s,
            eigenvalues_rad2_s2,
            basis_columns,
            maximum_parameter_work_j,
            cx,
        )?;
        let candidate = Self::prepare_source(
            solid,
            self.tones.angular_frequencies_rad_s(),
            self.tones.is_affine(),
            tone_budget,
            self.force_count,
            cx,
        )?;
        poll(cx)?;
        *self = candidate;
        Ok(report)
    }

    /// Borrow all physical state, free-coordinate labels and work/heat accounts.
    /// The immutable view cannot advance motion without its admitted source law.
    pub fn solid(&self) -> &DampedSpectralSolid {
        &self.solid
    }

    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.tones.angular_frequencies_rad_s()
    }
    pub fn is_affine(&self) -> bool {
        self.tones.is_affine()
    }

    /// Pure phase rotation for one original caller-ordered tone.
    /// Repeated rotations accumulate roundoff; absolute-clock quadratures are
    /// also valid inputs. Neither amplitudes nor phase are normalized.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine_n: f64,
        sine_n: f64,
    ) -> Result<(f64, f64)> {
        self.tones
            .advance_force_phase(tone, cosine_n, sine_n)
            .map_err(|_| invalid("invalid tone index or unrepresentable force phase"))
    }

    /// Instantaneous physical impulse; prior heat and all source phases remain.
    pub fn apply_momentum_impulse(
        &mut self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> core::result::Result<LinearSolidImpulseReport, LinearSolidImpulseError> {
        self.solid
            .apply_momentum_impulse(impulse_kg_m_s, maximum_impulse_error_kg_m_s, cx)
    }

    /// One constant-quadrature interval, measured from this interval's start.
    pub fn advance(
        &mut self,
        cosine_n: &[f64],
        sine_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.advance_source(cosine_n, sine_n, None, cx)
    }

    /// Exact linear quadrature envelopes; endpoints share the start's carrier
    /// basis. A finite change between intervals is not a momentum impulse.
    pub fn advance_affine(
        &mut self,
        cosine_start_n: &[f64],
        sine_start_n: &[f64],
        cosine_end_n: &[f64],
        sine_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.advance_source(
            cosine_start_n,
            sine_start_n,
            Some((cosine_end_n, sine_end_n)),
            cx,
        )
    }

    /// Continue constant tones, publishing their rotated phases only on success.
    pub fn advance_continuous(
        &mut self,
        cosine_n: &mut [f64],
        sine_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.advance_continuous_source(cosine_n, sine_n, None, cx)
    }

    /// After an accepted affine interval, replace start arrays by the rotated
    /// END quadratures. Every refusal preserves both arrays and physical state.
    pub fn advance_affine_continuous(
        &mut self,
        cosine_start_n: &mut [f64],
        sine_start_n: &mut [f64],
        cosine_end_n: &[f64],
        sine_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.advance_continuous_source(
            cosine_start_n,
            sine_start_n,
            Some((cosine_end_n, sine_end_n)),
            cx,
        )
    }

    fn advance_continuous_source(
        &mut self,
        cosine: &mut [f64],
        sine: &mut [f64],
        ends: Option<(&[f64], &[f64])>,
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.admit_source(cosine, sine, ends, cx)?;
        let (c, s) = ends.unwrap_or((&*cosine, &*sine));
        let (mut next_c, mut next_s) = (zeros(self.force_count)?, zeros(self.force_count)?);
        let k = self.tones.tone_count();
        for i in 0..self.force_count {
            if i % k == 0 {
                poll(cx)?;
            }
            (next_c[i], next_s[i]) = self.advance_force_phase(i % k, c[i], s[i])?;
        }
        let step = self.advance_source(cosine, sine, ends, cx)?;
        cosine.copy_from_slice(&next_c);
        sine.copy_from_slice(&next_s);
        Ok(step)
    }

    fn advance_source(
        &mut self,
        cosine: &[f64],
        sine: &[f64],
        ends: Option<(&[f64], &[f64])>,
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        self.admit_source(cosine, sine, ends, cx)?;
        let (step, time_s) = self.solid.solid.next_interval()?;
        let (c, s) = (
            self.project_forces(cosine, cx)?,
            self.project_forces(sine, cx)?,
        );
        let projected_ends = if let Some((ce, se)) = ends {
            Some((self.project_forces(ce, cx)?, self.project_forces(se, cx)?))
        } else {
            None
        };
        let n = self.solid.free_dofs().len();
        let k = self.tones.tone_count();
        let (mut q, mut v) = (zeros(n)?, zeros(n)?);
        let (mut work, mut work_correction, mut heat, mut heat_correction) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            poll(cx)?;
            let state = ModalAcousticState {
                displacement_m_sqrt_kg: self.solid.solid.modal_displacement[i],
                velocity_m_sqrt_kg_per_s: self.solid.solid.modal_velocity[i],
            };
            let range = i * k..(i + 1) * k;
            let result = if let Some((ce, se)) = &projected_ends {
                self.tones.advance_mode_affine_with_energy(
                    i,
                    state,
                    &c[range.clone()],
                    &s[range.clone()],
                    &ce[range.clone()],
                    &se[range],
                )
            } else {
                self.tones
                    .advance_mode_with_energy(i, state, &c[range.clone()], &s[range])
            };
            let frame = tone_result(result, cx)?;
            q[i] = frame.state.displacement_m_sqrt_kg;
            v[i] = frame.state.velocity_m_sqrt_kg_per_s;
            add(&mut work, &mut work_correction, frame.input_work_j);
            add(&mut heat, &mut heat_correction, frame.viscous_dissipation_j);
        }
        self.solid
            .accept_modal_path(q, v, work, heat, step, time_s, cx)
    }

    fn admit_source(
        &self,
        cosine: &[f64],
        sine: &[f64],
        ends: Option<(&[f64], &[f64])>,
        cx: &Cx<'_>,
    ) -> Result<()> {
        poll(cx)?;
        if self.tones.is_affine() != ends.is_some() {
            return Err(invalid("source envelopes must match the prepared tone law"));
        }
        for field in [cosine, sine]
            .into_iter()
            .chain(ends.into_iter().flat_map(|(c, s)| [c, s]))
        {
            if field.len() != self.force_count || field.iter().any(|value| !value.is_finite()) {
                return Err(invalid(
                    "complete finite DOF-major tone force arrays required",
                ));
            }
        }
        Ok(())
    }

    fn project_forces(&self, physical: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
        let n = self.solid.free_dofs().len();
        let k = self.tones.tone_count();
        let mut modal = zeros(self.force_count)?;
        for mode in 0..n {
            poll(cx)?;
            for row in 0..n {
                let phi = self.solid.solid.basis[row * n + mode];
                for tone in 0..k {
                    modal[mode * k + tone] =
                        phi.mul_add(physical[row * k + tone], modal[mode * k + tone]);
                }
            }
        }
        if modal.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable physical tone force projection"));
        }
        Ok(modal)
    }
}

fn tone_result<T>(
    value: core::result::Result<T, ModalAcousticTimeError>,
    cx: &Cx<'_>,
) -> Result<T> {
    value.map_err(|_| {
        poll(cx)
            .err()
            .unwrap_or_else(|| invalid("simultaneous-tone forcing refused"))
    })
}
