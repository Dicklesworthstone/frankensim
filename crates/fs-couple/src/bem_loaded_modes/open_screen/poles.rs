//! Loaded open-shell modes in the original displacement/pressure-jump equations.
//!
//! The boundary's analytic cosine-plus-Gram continuation matches its original
//! real-frequency discretization. A contour retains the full pressure block and
//! every original Maxwell state; neither H nor a material denominator is
//! eliminated into a meromorphic spectral quotient. Fixed balancing,
//! original-row audits and both nonlinear eigenvalue algorithms are shared with
//! the existing closed-body owner. These are poles of the selected finite
//! operator, not a continuum completeness or time-domain causality certificate.
use super::super::{PhysicalPencil, pencil::physical_eigenstate};
use super::*;
use fs_bem::reciprocal_radiation::open_screen::ScreenComplexLimits;
use fs_la::nonlinear_eig::beyn::{BeynOptions, BeynReport, Ellipse, beyn_with_cancel};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, NonlinearPole, PoleEvaluation, PoleOptions,
    newton_pole_with_cancel,
};

/// Explicit complex-domain, callback, projection and output limits. The native
/// Newton/Beyn options separately bound their numerical work and scratch.
#[derive(Clone, Copy, Debug)]
pub struct OpenScreenResonanceOptions {
    /// Original H/H' assembly, both half-planes, with a centered growth bound.
    pub boundary: ScreenComplexLimits,
    /// Includes EVERY source coordinate, interior pressure-jump unknown and
    /// retained Maxwell state: d = modes + pressure_unknowns + memory_count.
    pub maximum_dimension: usize,
    /// Conservative callback/seed/audit payload:
    /// 224*d^2 + 512*d + 8*m*p + 16*vertices + 24*m bytes.
    /// Source geometry, prepared BEM cache/workspace, native eigensolver scratch
    /// and returned-result payload have separate limits; excludes allocator/RSS.
    pub maximum_dense_bytes: usize,
    /// The exact cold source pairing uses this same explicit allowance. Each
    /// subsequent coupled callback requires 16*d^2 + 8*m*p units, in addition
    /// to the independently bounded boundary assembly. No mode is truncated.
    pub maximum_projection_work: u64,
    /// Complete boundary evaluations in each query, including initialization
    /// and final physical acceptance. Cold reference balancing costs one more
    /// separately bounded evaluation when this immutable pencil is prepared.
    pub maximum_operator_evaluations: usize,
    /// Logical output ceiling, checked BEFORE native iterations: one result for
    /// Newton, at most `probe_columns` for Beyn. Includes full nodal jump arrays
    /// and every material state. No partial spectrum is returned.
    pub maximum_result_bytes: usize,
    /// Positive real [rad/s] reference for fixed two-sided balancing. Changing
    /// the spectral point never changes the balance or source coordinates.
    pub reference_angular_frequency_rad_s: f64,
    /// Original mechanical, boundary and memory componentwise backward error,
    /// in (0,1).
    pub physical_residual_tolerance: f64,
}

/// Source-bound operator on [modal displacement, interior pressure jump, memory].
/// The source, full K/C, original Maxwell arms, surface quadrature and exact P1
/// work pairing are fixed. Equal-time material histories are never merged.
pub struct OpenScreenResonancePencil<'p, 's> {
    prepared: &'p PreparedOpenScreenModes<'s>,
    /// Mode-major exact weak rows G_mode[pressure_unknown].
    pairing: Vec<Vec<f64>>,
    row_scales: Vec<f64>,
    column_scales: Vec<f64>,
    settings: OpenScreenResonanceOptions,
}

/// An accepted ORIGINAL-coordinate eigenstate, with arbitrary amplitude fixed
/// by maximum structural displacement (or remaining field for zero motion).
#[derive(Clone, Debug)]
pub struct OpenScreenModePole {
    /// exp(-i*omega*t): a decaying mode has Im(omega)<0.
    pub angular_frequency_rad_s: C64,
    /// Complete unchanged source order; no diagonal or bright-mode reduction.
    pub modal_displacement: Vec<C64>,
    /// Complete original acoustic vertex order, including exact free-edge zeros.
    pub pressure_jump_pa: Vec<C64>,
    /// Original arm order, e_j = c_j.q - z_j [sqrt(J)]. The corresponding
    /// real-state spring energy is e_j^2/2; complex pole amplitudes are arbitrary.
    pub elastic_memory_sqrt_j: Vec<C64>,
    pub mechanical_residual: f64,
    pub boundary_residual: f64,
    pub memory_residual: f64,
    /// Final native Newton scaled correction and iteration count.
    pub relative_step: f64,
    pub iterations: usize,
}
impl OpenScreenModePole {
    pub fn decay_rate_per_s(&self) -> f64 {
        -self.angular_frequency_rad_s.im
    }
    /// No finite Q is invented for an undamped or nonoscillatory pole.
    pub fn quality_factor(&self) -> Option<f64> {
        let w = self.angular_frequency_rad_s;
        (w.im < 0.0 && w.re != 0.0).then(|| w.re.abs() / (-2.0 * w.im))
    }
}

/// Native contour/refinement evidence plus original physical modes. Numerical
/// rank/count agreement is not an interval completeness or continuum certificate.
#[derive(Clone, Debug)]
pub struct OpenScreenModeSpectrum {
    pub poles: Vec<OpenScreenModePole>,
    pub unrefined_values: Vec<C64>,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub argument_count: usize,
    pub rank: usize,
    pub nodes: usize,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
    pub operator_evaluations: usize,
}

impl<'s> PreparedOpenScreenModes<'s> {
    /// Prepare loaded-mode analysis of this SAME source and real-frequency
    /// boundary discretization. The open H is analytic at k=0; unlike closed
    /// Burton--Miller, no 1/k singularity or automatic exclusion is introduced.
    /// A multiple rigid zero pole can still refuse the native simple-pole solve.
    ///
    /// This is a bounded holomorphic finite-operator chart in both half-planes.
    /// Its outgoing angular Gram approximation requires independent mesh,
    /// quadrature and angular refinement. It is not the separate damped
    /// point-kernel `solve_laplace` image, an all-RHP positive-real certificate,
    /// a convolution-quadrature realization or a new transient integrator.
    pub fn resonance_pencil<'p>(
        &'p self,
        settings: OpenScreenResonanceOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OpenScreenResonancePencil<'p, 's>, BemModesError> {
        poll(&mut cancelled)?;
        let (m, p) = (self.source.mode_count(), self.screen.unknowns());
        let d = m
            .checked_add(p)
            .and_then(|d| d.checked_add(self.source.memory_count()))
            .ok_or_else(|| bad("open resonance dimension overflow"))?;
        let vertices = self.source.surface.vertices().len();
        let bytes = 224 * (d as u128).pow(2)
            + 512 * d as u128
            + 8 * m as u128 * p as u128
            + 16 * vertices as u128
            + 24 * m as u128;
        let work = 16 * (d as u128).pow(2) + 8 * m as u128 * p as u128;
        if d > MAX_NONLINEAR_POLE_DIMENSION
            || d > settings.maximum_dimension
            || bytes > settings.maximum_dense_bytes as u128
            || work > settings.maximum_projection_work as u128
            || settings.maximum_operator_evaluations == 0
            || !settings.reference_angular_frequency_rad_s.is_finite()
            || settings.reference_angular_frequency_rad_s <= 0.0
            || !settings.physical_residual_tolerance.is_finite()
            || settings.physical_residual_tolerance <= 0.0
            || settings.physical_residual_tolerance >= 1.0
        {
            return Err(bad(
                "open resonance requires complete dimensions, valid reference/accuracy and sufficient work/storage",
            ));
        }
        self.screen
            .admit_complex_domain(
                settings.reference_angular_frequency_rad_s / self.source.medium.sound_speed,
                0.0,
                settings.boundary,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        let rows: Vec<_> = self.source.weights.iter().map(Vec::as_slice).collect();
        let pairing = self
            .screen
            .velocity_port_pairing(
                &rows,
                m * p,
                settings.maximum_projection_work,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        let mut out = OpenScreenResonancePencil {
            prepared: self,
            pairing,
            row_scales: Vec::new(),
            column_scales: Vec::new(),
            settings,
        };
        let raw = out.assemble(
            C64::from_re(settings.reference_angular_frequency_rad_s),
            &mut 0,
            &mut cancelled,
        )?;
        (out.row_scales, out.column_scales) = raw.fixed_scales(d, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(out)
    }
}

impl OpenScreenResonancePencil<'_, '_> {
    pub fn source(&self) -> &OpenScreenLoadedModes {
        self.prepared.source
    }
    pub fn mode_count(&self) -> usize {
        self.source().mode_count()
    }
    pub fn pressure_unknowns(&self) -> usize {
        self.prepared.screen.unknowns()
    }
    pub fn memory_count(&self) -> usize {
        self.source().memory_count()
    }
    pub fn dimension(&self) -> usize {
        self.mode_count() + self.pressure_unknowns() + self.memory_count()
    }
    /// Original mode-major exact P1 work rows G, including signed cross terms.
    pub fn velocity_pressure_pairing(&self) -> &[Vec<f64>] {
        &self.pairing
    }
    /// Exact vertex-to-unknown map; a free-edge pressure jump is zero.
    pub fn pressure_dof_map(&self) -> &[Option<usize>] {
        self.prepared.screen.pressure_dof_map()
    }

    fn assemble(
        &self,
        omega: C64,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PhysicalPencil, BemModesError> {
        poll(cancelled)?;
        if !finite(omega) {
            return Err(bad("open resonance angular frequency must be finite"));
        }
        if *evaluations >= self.settings.maximum_operator_evaluations {
            return Err(BemModesError::EvaluationBudget {
                maximum: self.settings.maximum_operator_evaluations,
            });
        }
        *evaluations += 1;
        let source = self.source();
        let boundary = self
            .prepared
            .screen
            .complex_operator(
                divide(omega, source.medium.sound_speed),
                self.settings.boundary,
                &mut *cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        let (m, p, d) = (
            self.mode_count(),
            self.pressure_unknowns(),
            self.dimension(),
        );
        if boundary.dimension != p
            || boundary.matrix.len() != p * p
            || boundary.derivative_wavenumber.len() != p * p
        {
            return Err(bad(
                "open analytic boundary cardinality disagrees with the retained source",
            ));
        }
        let mut matrix = vec![C64::ZERO; d * d];
        let mut derivative = vec![C64::ZERO; d * d];
        let mut envelope = vec![0.0; d * d];
        let s = C64::new(omega.im, -omega.re);
        let w2 = omega * omega;
        let aw = omega.abs();
        for i in 0..m {
            poll(cancelled)?;
            for j in 0..m {
                let (k, c) = (source.stiffness[i * m + j], source.viscosity[i * m + j]);
                matrix[i * d + j] = C64::from_re(k) + s.scale(c);
                derivative[i * d + j] = C64::new(0.0, -c);
                envelope[i * d + j] = k.abs() + aw * c.abs();
                if i == j {
                    matrix[i * d + j] = matrix[i * d + j] - w2;
                    derivative[i * d + j] = derivative[i * d + j] - omega.scale(2.0);
                    envelope[i * d + j] += aw * aw;
                }
            }
            for j in 0..p {
                let g = self.pairing[i][j];
                matrix[i * d + m + j] = C64::from_re(g);
                envelope[i * d + m + j] = g.abs();
                matrix[(m + j) * d + i] = w2.scale(-source.medium.density * g);
                derivative[(m + j) * d + i] = omega.scale(-2.0 * source.medium.density * g);
                envelope[(m + j) * d + i] = source.medium.density * aw * aw * g.abs();
            }
        }
        for i in 0..p {
            poll(cancelled)?;
            for j in 0..p {
                let at = (m + i) * d + m + j;
                matrix[at] = boundary.matrix[i * p + j];
                derivative[at] = divide(
                    boundary.derivative_wavenumber[i * p + j],
                    source.medium.sound_speed,
                );
                envelope[at] = boundary.matrix[i * p + j].abs();
            }
        }
        source.memory.add_to_pencil(
            omega,
            m,
            p,
            &mut matrix,
            &mut derivative,
            &mut envelope,
            cancelled,
        )?;
        if matrix.iter().chain(&derivative).any(|&z| !finite(z))
            || envelope.iter().any(|v| !v.is_finite())
        {
            return Err(bad(
                "open coupled operator or analytic derivative overflows",
            ));
        }
        let scale = envelope
            .chunks(d)
            .map(|row| row.iter().sum::<f64>())
            .fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(bad("open coupled coefficient envelope is invalid"));
        }
        Ok(PhysicalPencil {
            value: PoleEvaluation {
                matrix,
                derivative,
                residual_scale: scale,
            },
            envelope,
        })
    }

    /// Unscaled SI-coordinate matrix and analytic d/domega, in the complete
    /// [displacement, interior jump, elastic memory] order. Fixed balancing is
    /// not applied. The full pencil is analytic at bare material rates, which
    /// may still be genuine eigenvalues for mechanically cancelling histories.
    pub fn original_pencil(
        &self,
        omega: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoleEvaluation, BemModesError> {
        Ok(self.assemble(omega, &mut 0, &mut cancelled)?.value)
    }

    fn scaled(
        &self,
        omega: C64,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, BemModesError> {
        self.assemble(omega, evaluations, cancelled)?
            .into_balanced(&self.row_scales, &self.column_scales)
    }

    fn admit_result_storage(&self, modes: usize) -> Result<(), BemModesError> {
        let one = 16
            * (self.mode_count() as u128
                + self.source().surface.vertices().len() as u128
                + self.memory_count() as u128)
            + 256;
        let bytes = modes as u128 * (one + 16) + 512;
        if bytes > self.settings.maximum_result_bytes as u128 {
            return Err(bad(
                "complete open resonance output exceeds its full-nodal result allowance",
            ));
        }
        Ok(())
    }

    fn observe(
        &self,
        pole: NonlinearPole,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<OpenScreenModePole, BemModesError> {
        let m = self.mode_count();
        let vector = physical_eigenstate(&pole, &self.column_scales, m)?;
        let raw = self.assemble(pole.value, evaluations, cancelled)?;
        let residuals = raw.physical_residuals(
            m,
            self.pressure_unknowns(),
            &vector,
            None,
            self.settings.physical_residual_tolerance,
            cancelled,
        )?;
        if pole.value.im > 0.0 {
            return Err(BemModesError::GrowingMode {
                angular_frequency_rad_s: pole.value,
            });
        }
        let memory = vector[m + self.pressure_unknowns()..].to_vec();
        let interior = &vector[m..m + self.pressure_unknowns()];
        let mut jump = Vec::with_capacity(self.pressure_dof_map().len());
        for (vertex, &dof) in self.pressure_dof_map().iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cancelled)?;
            }
            jump.push(dof.map_or(C64::ZERO, |i| interior[i]));
        }
        poll(cancelled)?;
        Ok(OpenScreenModePole {
            angular_frequency_rad_s: pole.value,
            modal_displacement: vector[..m].to_vec(),
            pressure_jump_pa: jump,
            elastic_memory_sqrt_j: memory,
            mechanical_residual: residuals[0],
            boundary_residual: residuals[1],
            memory_residual: residuals[2],
            relative_step: pole.relative_step,
            iterations: pole.iterations,
        })
    }

    /// Local simple-pole refinement by the existing bordered Newton owner.
    /// Only the initializer eliminates H and the scalar material equations to
    /// seed the full state from the supplied motion. A singular isolated H or
    /// exact bare material rate refuses this initializer. Contour discovery
    /// retains all states and needs neither inverse; it also admits modes with
    /// zero motion and nonzero mechanically cancelling material histories.
    pub fn resonance_pole(
        &self,
        omega0: C64,
        modal_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OpenScreenModePole, BemModesError> {
        poll(&mut cancelled)?;
        let (m, p, d) = (
            self.mode_count(),
            self.pressure_unknowns(),
            self.dimension(),
        );
        if modal_seed.len() != m
            || modal_seed.iter().any(|&v| !finite(v))
            || modal_seed.iter().all(|v| v.abs() == 0.0)
        {
            return Err(bad(
                "open resonance seed requires every finite source coordinate and nonzero motion",
            ));
        }
        self.admit_result_storage(1)?;
        let mut evaluations = 0;
        let raw = self.assemble(omega0, &mut evaluations, &mut cancelled)?;
        let mut boundary = vec![C64::ZERO; p * p];
        let mut pressure = vec![C64::ZERO; p];
        for i in 0..p {
            poll(&mut cancelled)?;
            for j in 0..p {
                boundary[i * p + j] = raw.value.matrix[(m + i) * d + m + j];
            }
            for j in 0..m {
                pressure[i] = pressure[i] - raw.value.matrix[(m + i) * d + j] * modal_seed[j];
            }
        }
        let lu = lu_complex_with_cancel(&boundary, p, &mut cancelled)
            .map_err(BemModesError::SeedSolve)?;
        lu.solve_with_cancel(&mut pressure, &mut cancelled)
            .map_err(BemModesError::SeedSolve)?;
        let memory = self
            .source()
            .memory
            .seed(omega0, modal_seed, &mut cancelled)?;
        let mut seed = modal_seed.to_vec();
        seed.extend(pressure);
        seed.extend(memory);
        for (value, &scale) in seed.iter_mut().zip(&self.column_scales) {
            *value = divide(*value, scale);
        }
        drop(raw);
        drop(lu);
        drop(boundary);
        let mut origin = None;
        let result = newton_pole_with_cancel(
            d,
            omega0,
            &seed,
            options,
            |w, cancel| {
                self.scaled(w, &mut evaluations, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: w }
                })
            },
            &mut cancelled,
        );
        if let Some(error) = origin {
            return Err(error);
        }
        self.observe(
            result.map_err(BemModesError::Pole)?,
            &mut evaluations,
            &mut cancelled,
        )
    }

    /// Discover roots of the complete holomorphic pressure/structure/material pencil.
    /// The ENTIRE ellipse and its interior must fit the original complex band
    /// and centered growth chart, using conservative geometric bounds before
    /// any contour solve. No pressure DOF, zero-opening mode or material
    /// history is discarded. Crossing a bare Maxwell rate does not subtract a
    /// meromorphic pole from the argument count.
    ///
    /// Native nested moments, argument count, conditioning and final original
    /// physical-row checks remain separate evidence. Physical identification
    /// and mesh/source/triangle/angular convergence are still required; a small
    /// residual is not a proof that every discovered representation root is a
    /// physical continuum resonance. Growing accepted roots are refused.
    pub fn resonance_spectrum(
        &self,
        region: Ellipse,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OpenScreenModeSpectrum, BemModesError> {
        poll(&mut cancelled)?;
        if !finite(region.center)
            || !region.real_radius.is_finite()
            || region.real_radius <= 0.0
            || !region.imaginary_radius.is_finite()
            || region.imaginary_radius <= 0.0
        {
            return Err(bad(
                "open resonance contour requires finite positive ellipse geometry",
            ));
        }
        self.admit_result_storage(options.probe_columns)?;
        let c = self.source().medium.sound_speed;
        self.prepared
            .screen
            .admit_complex_domain(
                (region.center.abs() + region.real_radius.max(region.imaginary_radius)) / c,
                (region.center.im.abs() + region.imaginary_radius) / c,
                self.settings.boundary,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        let mut evaluations = 0;
        let mut origin = None;
        let result = beyn_with_cancel(
            self.dimension(),
            region,
            options,
            |w, cancel| {
                self.scaled(w, &mut evaluations, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: w }
                })
            },
            &mut cancelled,
        );
        if let Some(error) = origin {
            return Err(error);
        }
        let BeynReport {
            poles,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
        } = result.map_err(BemModesError::Contour)?;
        self.admit_result_storage(poles.len())?;
        let mut observed = Vec::with_capacity(poles.len());
        for pole in poles {
            observed.push(self.observe(pole, &mut evaluations, &mut cancelled)?);
        }
        poll(&mut cancelled)?;
        Ok(OpenScreenModeSpectrum {
            poles: observed,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
            operator_evaluations: evaluations,
        })
    }
}

#[cfg(test)]
mod tests;
