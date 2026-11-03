//! Free modes of the original DG generator, including every material memory.
//!
//! The shared nonlinear-eigenvalue owner solves `(s I - L) x = 0`; no impedance
//! fit, frequency-independent loss substitution, or state reduction is used.
//! Physical fields and independent loss squares are observed in the original
//! coordinates. Their Hermitian energy identity describes a complex modal
//! envelope, not a finite-cycle average of a decaying real trajectory.

use super::{DgBoundary, DgContinuous, DgDuct, DgMedium, DuctError, PortId, finite, invalid, poll};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynReport, beyn_with_cancel};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, NonlinearPole, PoleEvaluation,
    newton_pole_with_cancel,
};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

pub use fs_la::nonlinear_eig::PoleOptions;
pub use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse as ResonanceRegion};

/// Original pressure participation of one homogeneous free mode.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPolePort {
    /// Original declared inflow port; its forcing is zero in this free mode.
    pub id: PortId,
    /// Pressure from the original power-conjugate observation column [Pa].
    pub pressure_pa: C64,
}

/// One finite free mode of the complete frozen DG spatial/material model.
#[derive(Clone, Debug, PartialEq)]
pub struct DgResonancePole {
    /// Laplace pole in exp(s t) [1/s]. Growth is retained, never clamped.
    pub laplace_pole_per_s: C64,
    /// omega = i*s under exp(-i omega t); decay has negative imaginary part.
    pub angular_frequency_rad_s: C64,
    /// Absolute oscillation frequency |Im(s)| [rad/s]; zero for real relaxation.
    pub oscillation_frequency_rad_s: f64,
    /// Signed amplitude decay rate -Re(s) [1/s]; negative denotes growth.
    pub amplitude_decay_rate_per_s: f64,
    /// Full original energy-coordinate order, unit numerical Euclidean norm.
    pub energy_coordinates: Vec<C64>,
    /// Physical pressure at original DG nodes for this modal amplitude [Pa].
    pub pressure_pa: Vec<C64>,
    /// Physical axial volume flow at original DG nodes [m3/s].
    pub volume_flow_m3_s: Vec<C64>,
    /// Original port observations; zero participation does not remove a mode.
    pub ports: Vec<DgPolePort>,
    /// Hermitian modal-envelope storage x* x / 2 [J], not a cycle average.
    pub stored_energy_j: f64,
    /// Included storage in every original wall/molecular auxiliary state [J].
    pub auxiliary_stored_energy_j: f64,
    /// Sum of original positive wall-loss squares for Re(x) and Im(x) [W].
    pub wall_dissipation_w: f64,
    /// Original molecular exchange loss, separate from wall dissipation [W].
    pub molecular_dissipation_w: f64,
    /// Original matched-boundary absorption for the complex modal envelope [W].
    pub matched_absorption_w: f64,
    /// Independently evaluated 2 Re(s) H + D_wall + D_molecular + D_boundary.
    pub decay_power_residual_w: f64,
    /// Infinity-norm backward residual of the original `(s I - L)x = 0`.
    pub relative_residual: f64,
    /// Completed evaluations of the shared local pole refinement.
    pub iterations: usize,
}

/// Numerically resolved modes inside one declared Laplace-plane ellipse.
/// Nested quadrature and count/rank agreement are convergence indicators,
/// not a rigorous enclosure or a proof that every physical mode was resolved.
#[derive(Clone, Debug, PartialEq)]
pub struct DgResonanceSpectrum {
    /// Complete physical observations, including zero-port-participation modes.
    pub poles: Vec<DgResonancePole>,
    /// Original spectral seeds before Newton refinement, for quadrature ladders.
    pub unrefined_laplace_poles_per_s: Vec<C64>,
    /// Numerically resolved argument-principle count, equal to the moment rank.
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    /// Smallest sampled numerical conditioning margin on the contour.
    pub minimum_reciprocal_condition: f64,
}

/// A failed refinement never returns a partial mode or changes source history.
#[derive(Clone, Debug, PartialEq)]
pub enum DgPoleError {
    /// Original source band, physical observation, capacity or cancellation.
    Source(DuctError),
    /// Shared local eigensolver admission, convergence or cancellation.
    Solve(NonlinearEigError),
    /// Shared contour admission, quadrature, count, conditioning or refinement.
    Contour(BeynError),
    /// A numerically near-real contour candidate did not refine locally as a
    /// real eigenpair of the original equations; its imaginary part is not cut.
    RealModeUnresolved {
        candidate: C64,
        real_refinement: C64,
    },
    /// Independent original-operator residual exceeds the requested tolerance.
    Residual {
        /// Observed infinity-norm backward residual.
        relative: f64,
        /// Original caller's required residual tolerance.
        tolerance: f64,
    },
}
impl From<DuctError> for DgPoleError {
    fn from(value: DuctError) -> Self {
        Self::Source(value)
    }
}
impl From<NonlinearEigError> for DgPoleError {
    fn from(value: NonlinearEigError) -> Self {
        Self::Solve(value)
    }
}
impl From<BeynError> for DgPoleError {
    fn from(value: BeynError) -> Self {
        Self::Contour(value)
    }
}
impl core::fmt::Display for DgPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-DUCT-POLE: {self:?}")
    }
}
impl std::error::Error for DgPoleError {}

/// Borrowed complete DG operator with fixed equilibration and cached norm.
/// Preparation scans the source once without copying its dense generator.
/// Repeated local/contour queries preserve all original coordinates and losses.
pub struct PreparedDgPoles<'a, M: DgMedium = GasState> {
    view: DgContinuous<'a, M>,
    scale: f64,
    norm_l: f64,
}

impl<M: DgMedium> DgDuct<M> {
    /// Prepare the complete frozen pole pencil once, without allocating another
    /// dense generator. `reference_rate_per_s` is a positive numerical scale,
    /// not a material coefficient. The scan is bounded by `maximum_work` using
    /// six work units per original matrix entry and the shared dense dimension
    /// cap. Each subsequent query has its own explicit solve/observation budget.
    pub fn prepare_resonance_poles(
        &self,
        reference_rate_per_s: f64,
        maximum_work: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDgPoles<'_, M>, DgPoleError> {
        poll(&mut cancelled)?;
        let view = self.continuous_operator();
        let n = view.dimension();
        if n == 0 || n > MAX_NONLINEAR_POLE_DIMENSION {
            return Err(invalid("complete DG pole dimension exceeds the shared capacity").into());
        }
        if !reference_rate_per_s.is_finite() || reference_rate_per_s <= 0.0 {
            return Err(invalid("DG pole reference rate must be positive and finite").into());
        }
        if 6 * (n as u128) * (n as u128) > maximum_work as u128 {
            return Err(invalid("DG pole preparation exceeds its work capacity").into());
        }
        let mut scale = reference_rate_per_s;
        for row in view.generator_per_s().chunks_exact(n) {
            poll(&mut cancelled)?;
            for &value in row {
                if !value.is_finite() {
                    return Err(invalid("DG pole generator must be finite").into());
                }
                scale = scale.max(value.abs());
            }
        }
        let mut square_l = 0.0;
        for row in view.generator_per_s().chunks_exact(n) {
            poll(&mut cancelled)?;
            square_l += row.iter().map(|a| (a / scale) * (a / scale)).sum::<f64>();
        }
        poll(&mut cancelled)?;
        Ok(PreparedDgPoles {
            view,
            scale,
            norm_l: det::sqrt(square_l),
        })
    }

    /// Refine one free mode of the original frozen spatial generator, with all
    /// declared inflows zero and the original boundary conditions unchanged.
    /// A nearby forced harmonic state is a useful seed for observable modes.
    /// No response-bandwidth Q is converted into a decay rate.
    ///
    /// Uses the shared bounded/cancellable Newton pole solver. Constant scalar
    /// equilibration leaves the original state coordinates and analytic pencil
    /// derivative intact. Options include conservative observation payload/work
    /// in addition to the shared solver; the borrowed source is not reallocated.
    ///
    /// Oscillatory seeds and results must lie within any original wall-law
    /// frequency band. Pure real material-relaxation modes remain available,
    /// with zero oscillation frequency and no assigned Q. A wall law's real-axis
    /// frequency admission does not certify its continuation to complex poles.
    /// These are numerical modes of the retained spatial/material realization,
    /// not a continuum eigenvalue enclosure, complete spectrum, transfer-pole
    /// observability claim, or nonlinear/transient resonance prediction.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgResonancePole, DgPoleError> {
        self.prepare_resonance_poles(options.spectral_scale, options.max_work, &mut cancelled)?
            .resonance_pole(laplace_seed_per_s, state_seed, options, cancelled)
    }

    /// Discover complete free modes in a declared Laplace-plane region without
    /// individual pole or state seeds. See `PreparedDgPoles::resonance_spectrum`
    /// for numerical count, physical-band and resource semantics.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgResonanceSpectrum, DgPoleError> {
        self.prepare_resonance_poles(
            options.refinement.spectral_scale,
            options.max_work,
            &mut cancelled,
        )?
        .resonance_spectrum(region, options, cancelled)
    }
}

impl<M: DgMedium> PreparedDgPoles<'_, M> {
    /// Number of retained pressure, flow and original material-memory states.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.view.dimension()
    }

    fn evaluate(
        &self,
        s: C64,
        spectral_scale: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        let n = self.dimension();
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        for i in 0..n {
            if cancelled() {
                return Err(NonlinearEigError::Cancelled);
            }
            for j in 0..n {
                matrix[i * n + j] =
                    C64::from_re(-self.view.generator_per_s()[i * n + j] / self.scale);
            }
            matrix[i * n + i] = matrix[i * n + i] + C64::new(s.re / self.scale, s.im / self.scale);
            derivative[i * n + i] = C64::from_re(1.0 / self.scale);
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale: self.norm_l
                + s.abs().max(spectral_scale) / self.scale * det::sqrt(n as f64),
        })
    }

    /// Refine a seeded free mode using the prepared original pencil and the
    /// same independent physical observations as `DgDuct::resonance_pole`.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgResonancePole, DgPoleError> {
        poll(&mut cancelled)?;
        let view = self.view;
        if !finite(laplace_seed_per_s) {
            return Err(invalid("DG pole seed must be finite").into());
        }
        if laplace_seed_per_s.im != 0.0 {
            view.admit_frequency(laplace_seed_per_s.im.abs())?;
        }
        let n = view.dimension();
        let ports = view.ports().len();
        if n > options.max_dimension {
            return Err(NonlinearEigError::BadInput {
                what: "complete DG pole dimension exceeds the declared/shared capacity",
            }
            .into());
        }
        let bytes = (n as u128) * 128 + (ports as u128) * 64;
        // Include the source scan, callback assembly and physical observations.
        // The scan charge also covers the one-shot DgDuct convenience method.
        let work = (18 + 6 * options.max_iterations as u128) * (n as u128) * (n as u128)
            + 8 * (n as u128) * (ports as u128);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(
                invalid("complete DG pole observations exceed payload/work capacity").into(),
            );
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let pole = newton_pole_with_cancel(
            n,
            laplace_seed_per_s,
            state_seed,
            &numerical,
            |s, cancel| self.evaluate(s, options.spectral_scale, cancel),
            &mut cancelled,
        )?;
        self.observe(pole, options.residual_tolerance, &mut cancelled)
    }

    /// Discover finite free modes with shared Beyn contour integration, then
    /// independently check/observe every original DG eigenvector. The complete
    /// linear pencil is entire: no wall-memory denominators are eliminated.
    /// The fixed preparation scale preserves analyticity at all contour nodes.
    ///
    /// Every returned oscillatory pole must pass the original wall frequency
    /// band. Pure real relaxation modes remain available; contour evaluations
    /// are analytic numerical probes, not physical-frequency predictions.
    /// A candidate within the pole correction tolerance of the real axis is
    /// refined again from its real state on the original real pencil. Only an
    /// independently converged local real eigenpair receives zero frequency.
    /// No invalid/out-of-band root is silently removed. An unresolved count,
    /// rank, cluster, contour, budget, cancellation or physical mode refuses the
    /// whole query, without modifying source state. Repeated roots can refuse.
    ///
    /// Physical output payload and callback/observation work are reserved in
    /// addition to the shared solver's resources, before allocating matrices.
    /// Numerical quadrature agreement is not a complete-spectrum certificate,
    /// continuum error enclosure, or evidence that an internal mode is audible.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgResonanceSpectrum, DgPoleError> {
        poll(&mut cancelled)?;
        let n = self.dimension();
        if n > options.max_dimension
            || n > options.refinement.max_dimension
            || options.probe_columns == 0
            || options.probe_columns > n
        {
            return Err(invalid("complete DG spectrum dimension or probe capacity refused").into());
        }
        let (n, ports, modes) = (
            n as u128,
            self.view.ports().len() as u128,
            options.probe_columns as u128,
        );
        let refinements = modes * options.refinement.max_iterations as u128;
        // Reserve a second local refinement for every potential real mode.
        // Complex contour roundoff must not fabricate an oscillation frequency.
        let evaluations = options.nodes as u128 + 2 * refinements;
        let bytes = modes * (128 * n + 64 * ports) + 64 * modes;
        let work = (6 + 6 * evaluations + 12 * modes) * n * n
            + 8 * modes * n * ports
            + refinements * (n + 1) * (n + 1) * (n + 17);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(
                invalid("complete DG spectrum observations exceed payload/work capacity").into(),
            );
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
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
        } = beyn_with_cancel(
            self.dimension(),
            region,
            &numerical,
            |s, cancel| self.evaluate(s, options.refinement.spectral_scale, cancel),
            &mut cancelled,
        )?;
        let mut physical = Vec::with_capacity(poles.len());
        for mut pole in poles {
            let correction_limit = options.refinement.step_tolerance
                * options.refinement.spectral_scale.max(pole.value.abs());
            if pole.value.im != 0.0 && pole.value.im.abs() <= correction_limit {
                let real_seed: Vec<_> = pole.vector.iter().map(|x| C64::from_re(x.re)).collect();
                let mut real = newton_pole_with_cancel(
                    self.dimension(),
                    C64::from_re(pole.value.re),
                    &real_seed,
                    &options.refinement,
                    |s, cancel| self.evaluate(s, options.refinement.spectral_scale, cancel),
                    &mut cancelled,
                )?;
                if real.value.im != 0.0 || (real.value - pole.value).abs() > correction_limit {
                    return Err(DgPoleError::RealModeUnresolved {
                        candidate: pole.value,
                        real_refinement: real.value,
                    });
                }
                if !region.contains(real.value) {
                    return Err(BeynError::PoleOutside { value: real.value }.into());
                }
                real.iterations += pole.iterations;
                pole = real;
            }
            physical.push(self.observe(
                pole,
                options.refinement.residual_tolerance,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(DgResonanceSpectrum {
            poles: physical,
            unrefined_laplace_poles_per_s: unrefined_values,
            argument_count,
            resolved_rank: rank,
            quadrature_nodes: nodes,
            argument_integral,
            coarse_argument_integral,
            moment_relative_change,
            minimum_reciprocal_condition,
        })
    }

    fn observe(
        &self,
        pole: NonlinearPole,
        tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgResonancePole, DgPoleError> {
        let view = self.view;
        let duct = view.duct;
        let n = self.dimension();
        let ports = view.ports().len();
        let generator = view.generator_per_s();
        let scale = self.scale;
        let s = pole.value;
        if s.im != 0.0 {
            view.admit_frequency(s.im.abs())?;
        }
        let x = pole.vector;
        let norm_x = x.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        let (mut norm_a, mut norm_r) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut residual = C64::new(s.re / scale, s.im / scale) * x[i];
            let mut row = s.abs() / scale;
            for j in 0..n {
                let a = generator[i * n + j] / scale;
                residual = residual - x[j].scale(a);
                row += a.abs();
            }
            norm_a = norm_a.max(row);
            norm_r = norm_r.max(residual.abs());
        }
        let relative = norm_r / (norm_a * norm_x);
        if !relative.is_finite() || relative > tolerance {
            return Err(DgPoleError::Residual {
                relative,
                tolerance,
            });
        }
        let mut native = duct.zero_state();
        for (to, from) in native.tilde.iter_mut().zip(&x) {
            *to = from.re;
        }
        let (real_p, real_q) = duct.from_energy_coords(&native.tilde);
        let real_wall = duct.wall_dissipated_power_w(&native)?;
        let real_molecular = duct.molecular_dissipated_power_w(&native)?;
        poll(&mut cancelled)?;
        for (to, from) in native.tilde.iter_mut().zip(&x) {
            *to = from.im;
        }
        let (imag_p, imag_q) = duct.from_energy_coords(&native.tilde);
        let wall = real_wall + duct.wall_dissipated_power_w(&native)?;
        let molecular = real_molecular + duct.molecular_dissipated_power_w(&native)?;
        let pressure: Vec<_> = real_p
            .into_iter()
            .zip(imag_p)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let flow: Vec<_> = real_q
            .into_iter()
            .zip(imag_q)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let mut boundary = 0.0;
        for (bc, node, gas, area) in [
            (duct.inlet, 0, &duct.gas, duct.inlet_area_m2()),
            (
                duct.outlet,
                duct.field_dofs - 1,
                &duct.outlet_gas,
                duct.outlet_area_m2(),
            ),
        ] {
            if bc == DgBoundary::Matched {
                let z = gas.density_kg_m3() * gas.sound_speed_m_s() / area;
                boundary += 0.5 * (pressure[node].norm_sq() / z + z * flow[node].norm_sq());
            }
        }
        let mut port_outputs = Vec::with_capacity(ports);
        for (id, column) in view.ports() {
            poll(&mut cancelled)?;
            let pressure = column
                .iter()
                .zip(&x)
                .fold(C64::ZERO, |a, (&b, &x)| a + x.scale(b));
            port_outputs.push(DgPolePort {
                id: *id,
                pressure_pa: pressure,
            });
        }
        let energy = 0.5 * x.iter().map(|x| x.norm_sq()).sum::<f64>();
        let auxiliary = 0.5
            * x[2 * duct.field_dofs..]
                .iter()
                .map(|x| x.norm_sq())
                .sum::<f64>();
        let balance = 2.0 * s.re * energy + wall + molecular + boundary;
        if pressure.iter().chain(&flow).any(|&x| !finite(x))
            || port_outputs.iter().any(|p| !finite(p.pressure_pa))
            || [energy, auxiliary, wall, molecular, boundary, balance]
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(invalid("DG pole physical field, storage or loss overflow").into());
        }
        poll(&mut cancelled)?;
        Ok(DgResonancePole {
            laplace_pole_per_s: s,
            angular_frequency_rad_s: C64::new(-s.im, s.re),
            oscillation_frequency_rad_s: s.im.abs(),
            amplitude_decay_rate_per_s: -s.re,
            energy_coordinates: x,
            pressure_pa: pressure,
            volume_flow_m3_s: flow,
            ports: port_outputs,
            stored_energy_j: energy,
            auxiliary_stored_energy_j: auxiliary,
            wall_dissipation_w: wall,
            molecular_dissipation_w: molecular,
            matched_absorption_w: boundary,
            decay_power_residual_w: balance,
            relative_residual: relative,
            iterations: pole.iterations,
        })
    }
}
