//! Free modes of the original DG generator, including every material memory.
//!
//! The shared nonlinear-eigenvalue owner solves `(s I - L) x = 0`; no impedance
//! fit, frequency-independent loss substitution, or state reduction is used.
//! Physical fields and independent loss squares are observed in the original
//! coordinates. Their Hermitian energy identity describes a complex modal
//! envelope, not a finite-cycle average of a decaying real trajectory.

use super::{DgBoundary, DgDuct, DgMedium, DuctError, PortId, finite, invalid, poll};
use fs_la::nonlinear_eig::{NonlinearEigError, PoleEvaluation, newton_pole_with_cancel};
use fs_math::{c64::C64, det};

pub use fs_la::nonlinear_eig::PoleOptions;

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

/// A failed refinement never returns a partial mode or changes source history.
#[derive(Clone, Debug, PartialEq)]
pub enum DgPoleError {
    /// Original source band, physical observation, capacity or cancellation.
    Source(DuctError),
    /// Shared local eigensolver admission, convergence or cancellation.
    Solve(NonlinearEigError),
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
impl core::fmt::Display for DgPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-DUCT-POLE: {self:?}")
    }
}
impl std::error::Error for DgPoleError {}

impl<M: DgMedium> DgDuct<M> {
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
        poll(&mut cancelled)?;
        let view = self.continuous_operator();
        if !finite(laplace_seed_per_s) {
            return Err(invalid("DG pole seed must be finite").into());
        }
        if laplace_seed_per_s.im != 0.0 {
            view.admit_frequency(laplace_seed_per_s.im.abs())?;
        }
        let n = view.dimension();
        let ports = view.ports().len();
        if n > options.max_dimension || n > fs_la::nonlinear_eig::MAX_NONLINEAR_POLE_DIMENSION {
            return Err(NonlinearEigError::BadInput {
                what: "complete DG pole dimension exceeds the declared/shared capacity",
            }
            .into());
        }
        let bytes = (n as u128) * 128 + (ports as u128) * 64;
        let work = 12 * (n as u128) * (n as u128) + 8 * (n as u128) * (ports as u128);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(
                invalid("complete DG pole observations exceed payload/work capacity").into(),
            );
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let generator = view.generator_per_s();
        let scale = generator
            .iter()
            .map(|a| a.abs())
            .fold(options.spectral_scale, f64::max);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("DG pole operator scale must be positive and finite").into());
        }
        let mut square_l = 0.0;
        for row in generator.chunks_exact(n) {
            poll(&mut cancelled)?;
            square_l += row.iter().map(|a| (a / scale) * (a / scale)).sum::<f64>();
        }
        let norm_l = det::sqrt(square_l);
        let pole = newton_pole_with_cancel(
            n,
            laplace_seed_per_s,
            state_seed,
            &numerical,
            |s, cancel| {
                let mut matrix = vec![C64::ZERO; n * n];
                let mut derivative = vec![C64::ZERO; n * n];
                for i in 0..n {
                    if cancel() {
                        return Err(NonlinearEigError::Cancelled);
                    }
                    for j in 0..n {
                        matrix[i * n + j] = C64::from_re(-generator[i * n + j] / scale);
                    }
                    matrix[i * n + i] = matrix[i * n + i] + C64::new(s.re / scale, s.im / scale);
                    derivative[i * n + i] = C64::from_re(1.0 / scale);
                }
                Ok(PoleEvaluation {
                    matrix,
                    derivative,
                    residual_scale: norm_l
                        + s.abs().max(options.spectral_scale) / scale * det::sqrt(n as f64),
                })
            },
            &mut cancelled,
        )?;
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
        if !relative.is_finite() || relative > options.residual_tolerance {
            return Err(DgPoleError::Residual {
                relative,
                tolerance: options.residual_tolerance,
            });
        }
        let mut native = self.zero_state();
        for (to, from) in native.tilde.iter_mut().zip(&x) {
            *to = from.re;
        }
        let (real_p, real_q) = self.from_energy_coords(&native.tilde);
        let real_wall = self.wall_dissipated_power_w(&native)?;
        let real_molecular = self.molecular_dissipated_power_w(&native)?;
        poll(&mut cancelled)?;
        for (to, from) in native.tilde.iter_mut().zip(&x) {
            *to = from.im;
        }
        let (imag_p, imag_q) = self.from_energy_coords(&native.tilde);
        let wall = real_wall + self.wall_dissipated_power_w(&native)?;
        let molecular = real_molecular + self.molecular_dissipated_power_w(&native)?;
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
            (self.inlet, 0, &self.gas, self.inlet_area_m2()),
            (
                self.outlet,
                self.field_dofs - 1,
                &self.outlet_gas,
                self.outlet_area_m2(),
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
            * x[2 * self.field_dofs..]
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
