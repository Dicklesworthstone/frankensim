//! Derivatives of complete forced, autonomous and conservative HB waveforms.
//!
//! At fixed retained order, implicit differentiation gives `R_X X_p = -R_p`.
//! The original real AFT Jacobian and fs-la LU are the numerical owners. One
//! factor serves many physical parameter directions without nonlinear solves.
//! Every DC/real/imaginary coordinate participates, including mode coupling,
//! damping and explicit material memory. No diagonal-mode approximation,
//! regularization or stability inference is introduced.

use crate::{
    AftWorkspace, ForcedParameterPoint, HbAnchor, HbBudget, HbOrbit, HbScales, OrbitError,
    OrbitProblem, Packing, TAU, hb_residual,
    phase::PhaseCondition,
    scaling::{AdmittedScales, length},
    tangent,
};
use fs_la::factor::{Lu, lu};
use fs_math::{c64::C64, dd::Dd, det};

mod autonomous;
pub use autonomous::{
    AutonomousMeanObjective, AutonomousObjectiveGradient, AutonomousOrbitSensitivity,
    AutonomousSensitivityPoint, PreparedAutonomousAdjoint, PreparedAutonomousSensitivity,
};
mod adjoint;
pub use adjoint::{ForcedObjectiveGradient, PreparedForcedAdjoint};
mod objective;
pub use objective::{ForcedMeanObjective, ForcedMeanObjectiveControls};
mod backbone;
use backbone::{AmplitudeClosure, SensitivityFactor};
pub use backbone::{
    BackboneSensitivityAnchor, BackboneSensitivityDirection, PreparedBackboneSensitivity,
};

/// Explicit resources and independently authored numerical acceptance limits.
#[derive(Clone, Copy, Debug)]
pub struct ForcedSensitivityControls {
    /// Coarse AFT quadrature. An independent `2*samples+1` grid is also used.
    /// Must contain at least `4*N+4` samples for the retained order N.
    pub samples_per_period: usize,
    /// Full real unknown dimension, including autonomous frequency when present,
    /// also hard-capped at 1024.
    pub maximum_unknowns: usize,
    /// Conservative peak numerical payload for preparation and one direction.
    /// Excludes caller/source storage, allocator metadata and callback internals.
    pub maximum_payload_bytes: usize,
    /// Conservative arithmetic-work allowance for preparation and one direction;
    /// each reused direction has this same allowance. Callback work is excluded.
    pub maximum_work: u64,
    /// Maximum normalized Parseval norm of the original retained base balance
    /// on EITHER grid, including phase/amplitude rows when present.
    /// Previous convergence fields are never trusted.
    pub maximum_base_residual: f64,
    /// Frobenius norm of the difference of the two row/column-scaled exact
    /// coefficient Jacobians, including frequency/phase/amplitude borders.
    /// Agreement is measured, not an aliasing proof.
    pub maximum_tangent_sampling_defect: f64,
    /// Maximum fs-la Hager/Higham 1-norm condition ESTIMATE of the scaled
    /// Jacobian. This is not a certified upper bound on forward sensitivity.
    pub maximum_condition_estimate: f64,
    /// Maximum componentwise backward error of square differentiated systems.
    /// For the rectangular conservative system, maximum absolute row defect
    /// divided by maximum absolute row action; independent physical residuals
    /// below still bound every row, including weak/zero harmonic equations.
    pub maximum_relative_residual: f64,
    /// Maximum Parseval norm of the differentiated equations on BOTH grids,
    /// in equation-scale units per declared parameter-scale increment,
    /// combined with normalized phase/amplitude defects when present. Forced
    /// objective adjoints use this same absolute allowance for their normalized
    /// dual equations, direction-column sampling change and objective-gradient
    /// sampling change; objective units are supplied by `prepare_adjoint`.
    pub maximum_derivative_residual: f64,
}

/// Physical parameter chart for one directional derivative.
#[derive(Clone, Copy, Debug)]
pub struct ForcedSensitivityDirection {
    /// Current physical parameter value, forwarded unchanged to the callback.
    pub parameter: f64,
    /// Positive physical parameter increment used only for conditioning and
    /// absolute derivative-residual acceptance. Output is still per parameter unit.
    pub parameter_scale: f64,
    /// Derivative of the authored forcing clock in radians/time per parameter.
    /// The COMPLETE residual callback must include this changing phase/time chart.
    pub angular_frequency_derivative: f64,
}

/// Complete derivative of the retained real periodic waveform.
#[derive(Clone, Debug)]
pub struct ForcedOrbitSensitivity {
    /// Base angular frequency, not a differentiated state coordinate.
    pub omega: f64,
    /// Authored clock derivative; no frequency law is inferred from the source.
    pub angular_frequency_derivative: f64,
    /// `coeffs[state][harmonic] = dX_n/dp` in original state/parameter units.
    /// DC is real; positive frequencies are Fourier HALF-amplitudes as in HbOrbit.
    /// These are derivatives at fixed phase, not at fixed absolute time when
    /// the frequency changes. They are not derivatives of peak magnitudes.
    pub coeffs: Vec<Vec<C64>>,
    /// Maximum componentwise backward error of the fine-grid linear solve.
    pub relative_residual: f64,
    /// Maximum normalized differentiated-equation norm over the two grids.
    pub derivative_residual: f64,
}

/// No partially accepted derivative is returned on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ForcedSensitivityError {
    /// Invalid domain, shape, scale or unrepresentable arithmetic.
    Invalid(&'static str),
    /// An exact complete original island state tangent is required.
    TangentUnavailable,
    /// Original AFT tangent failure, propagated without finite-difference fallback.
    Orbit(OrbitError),
    /// Explicit numerical payload or operation allowance was exceeded.
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    /// Original retained base equations fail independent reevaluation.
    BaseUnresolved { residual: f64, limit: f64 },
    /// The two independently sampled exact Jacobians disagree.
    SamplingUnresolved { defect: f64, limit: f64 },
    /// Original retained Jacobian is singular; no regularization is attempted.
    Singular,
    /// Scaled condition estimate exceeds the authored limit or is nonfinite.
    IllConditioned { estimate: f64, limit: f64 },
    /// A complete residual-direction declaration refused or was incomplete.
    DirectionRefused { samples: usize },
    /// An autonomous residual-frequency column refused or was incomplete.
    FrequencyDirectionRefused { samples: usize },
    /// A sampled original island-direction declaration refused or was incomplete.
    IslandDirectionRefused { sample: usize, samples: usize },
    /// A scalar observation or its full state derivative refused or was incomplete.
    ObjectiveRefused { sample: usize, samples: usize },
    /// Cycle means or their complete coefficient gradients disagree across grids.
    ObjectiveSamplingUnresolved {
        value_defect: f64,
        gradient_defect: f64,
        value_limit: f64,
        gradient_limit: f64,
    },
    /// Retained backward error or independent differentiated equations failed.
    DerivativeUnresolved {
        relative_residual: f64,
        derivative_residual: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for ForcedSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-ORBIT-FORCED-SENSITIVITY: {self:?}")
    }
}
impl std::error::Error for ForcedSensitivityError {}
impl From<OrbitError> for ForcedSensitivityError {
    fn from(value: OrbitError) -> Self {
        Self::Orbit(value)
    }
}
type Result<T> = core::result::Result<T, ForcedSensitivityError>;

/// Cached exact original AFT linearization and a single conditioned factor.
///
/// This is a snapshot: subsequent directions must describe the SAME original
/// source, coefficients, period and phase origin used in preparation. Keeping
/// the source unchanged is the caller's obligation. No base-orbit refitting or
/// parameter mutation occurs here. Resource accounting excludes responses retained
/// by a caller across separate queries. Preparation has bounded noninterruptible
/// source/tangent and LU stages; cancellation is polled between those stages.
pub struct PreparedForcedSensitivity {
    pack: Packing,
    coefficients: Vec<f64>,
    omega: f64,
    coarse: AftWorkspace,
    fine: AftWorkspace,
    coarse_jacobian: Vec<f64>,
    fine_jacobian: Vec<f64>,
    normalized: Vec<f64>,
    factor: SensitivityFactor,
    scales: AdmittedScales,
    phase: Option<PhaseCondition>,
    amplitude: Option<AmplitudeClosure>,
    controls: ForcedSensitivityControls,
    /// Maximum recomputed normalized original retained base balance.
    pub base_residual: f64,
    /// Measured difference of the two conditioned coefficient Jacobians.
    pub tangent_sampling_defect: f64,
    /// Scaled 1-norm estimate. It is diagnostic, not a conditioning certificate.
    pub condition_estimate: f64,
}

impl PreparedForcedSensitivity {
    /// Prepare all coupled coefficient derivatives of an explicitly forced orbit.
    /// Uses the shared exact AFT tangent, the shared physical scaling chart and
    /// fs-la partial-pivot LU. Autonomous phase-neutral cycles are refused.
    ///
    /// This differentiates the FIXED-ORDER sampled HB equations. Independently
    /// refine the base and derivative harmonic orders for physical conclusions;
    /// neither two-grid agreement nor this condition estimate bounds the infinite
    /// harmonic tail, proves branch uniqueness, or establishes Floquet stability.
    pub fn new<P: OrbitProblem>(
        problem: &P,
        orbit: &HbOrbit,
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            problem, orbit, controls, scales, None, None, None, cancelled,
        )
    }

    // Forced, isolated-cycle and conservative-family sensitivities share this
    // original AFT owner. An independent frequency column and phase/amplitude
    // closures change the factor shape, never remove original balance rows.
    #[allow(clippy::too_many_arguments)] // complete source and optional physical chart closures
    fn prepare<P: OrbitProblem>(
        problem: &P,
        orbit: &HbOrbit,
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        phase_reference: Option<&[f64]>,
        amplitude_anchor: Option<BackboneSensitivityAnchor>,
        mut frequency_column: Option<
            &mut dyn FnMut(AutonomousSensitivityPoint<'_>, &mut [f64]) -> bool,
        >,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let (pack, coefficients, admitted) = admit(
            problem,
            orbit,
            controls,
            scales,
            phase_reference.is_some(),
            amplitude_anchor.is_some(),
        )?;
        let amplitude = amplitude_anchor
            .map(|anchor| AmplitudeClosure::new(anchor, &pack, &coefficients, &admitted))
            .transpose()?;
        let phase = phase_reference
            .map(|reference| {
                PhaseCondition::new(
                    problem,
                    reference,
                    &HbBudget {
                        harmonics: pack.n_harm,
                        max_newton: 1,
                        tolerance: 0.5,
                        deflate_equilibrium: false,
                    },
                    scales,
                )
            })
            .transpose()?;
        let phase_defect = phase.as_ref().map_or(0.0, |p| p.value(&coefficients));
        let closure_defect = amplitude.as_ref().map_or(phase_defect.abs(), |a| {
            det::hypot(phase_defect, a.base_defect)
        });
        if !closure_defect.is_finite() || closure_defect > controls.maximum_base_residual {
            return Err(ForcedSensitivityError::BaseUnresolved {
                residual: closure_defect,
                limit: controls.maximum_base_residual,
            });
        }
        let callback = problem
            .island_state_jacobian()
            .ok_or(ForcedSensitivityError::TangentUnavailable)?;
        let mut tangent = tangent::Workspace::new(&pack)?;
        let mut build = |samples| -> Result<_> {
            poll(&mut cancelled)?;
            let mut aft = AftWorkspace::with_samples(&pack, samples, 0)?;
            let mut residual = vec![f64::NAN; pack.coeff_len()];
            // Validate full frequency-port matrices before the original residual
            // indexes them. This is source admission, never a replacement port.
            for n in 0..=pack.n_harm {
                let port = problem.port(C64::new(0.0, n as f64 * orbit.omega));
                if port.len() != pack.d * pack.d
                    || port.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
                {
                    return Err(invalid("complete finite original port matrices required"));
                }
            }
            hb_residual(
                &CheckedSource(problem),
                &pack,
                &coefficients,
                orbit.omega,
                &mut aft,
                &mut residual,
            );
            if residual.iter().any(|v| !v.is_finite()) {
                return Err(invalid("original base equations refused"));
            }
            let base = det::hypot(admitted.row_norm(&residual), closure_defect);
            if !base.is_finite() || base > controls.maximum_base_residual {
                return Err(ForcedSensitivityError::BaseUnresolved {
                    residual: base,
                    limit: controls.maximum_base_residual,
                });
            }
            poll(&mut cancelled)?;
            tangent.prepare(problem, callback, &pack, orbit.omega, &aft)?;
            let nc = pack.coeff_len();
            let n = admitted.columns.len();
            let rows = n + usize::from(amplitude.is_some());
            let mut original = vec![0.0; rows * n];
            let mut normalized = vec![0.0; rows * n];
            for row in 0..nc {
                poll(&mut cancelled)?;
                for col in 0..nc {
                    let entry = tangent.coefficient_entry(row, col, &pack);
                    original[row * n + col] = entry;
                    normalized[row * n + col] =
                        scaled(entry, admitted.columns[col], admitted.rows[row])?;
                }
            }
            if let Some(phase) = &phase {
                let mut column = vec![f64::NAN; nc];
                let point = AutonomousSensitivityPoint {
                    angular_frequency: orbit.omega,
                    harmonics: pack.n_harm,
                    aft_samples: samples,
                    coefficients: &coefficients[..nc],
                };
                poll(&mut cancelled)?;
                if !frequency_column
                    .as_mut()
                    .is_some_and(|f| f(point, &mut column))
                    || column.iter().any(|v| !v.is_finite())
                {
                    return Err(ForcedSensitivityError::FrequencyDirectionRefused { samples });
                }
                for row in 0..nc {
                    original[row * n + nc] = column[row];
                    normalized[row * n + nc] =
                        scaled(column[row], admitted.columns[nc], admitted.rows[row])?;
                    // The phase row is dimensionless. Keep it in its normalized
                    // chart instead of forming possibly overflowing 1/state units.
                    normalized[nc * n + row] = phase.scaled_derivative(row, admitted.columns[row]);
                }
            }
            if let Some(amplitude) = &amplitude {
                amplitude.write_row(&mut normalized[n * n..]);
            }
            Ok((aft, original, normalized, base))
        };
        let (coarse, coarse_jacobian, coarse_normalized, coarse_base) =
            build(controls.samples_per_period)?;
        let (fine, fine_jacobian, normalized, fine_base) =
            build(2 * controls.samples_per_period + 1)?;
        let tangent_sampling_defect = length(
            normalized
                .iter()
                .zip(&coarse_normalized)
                .map(|(a, b)| a - b),
        );
        if !tangent_sampling_defect.is_finite()
            || tangent_sampling_defect > controls.maximum_tangent_sampling_defect
        {
            return Err(ForcedSensitivityError::SamplingUnresolved {
                defect: tangent_sampling_defect,
                limit: controls.maximum_tangent_sampling_defect,
            });
        }
        poll(&mut cancelled)?;
        let (factor, condition_estimate) = SensitivityFactor::new(
            &normalized,
            admitted.columns.len(),
            controls.maximum_condition_estimate,
            amplitude.is_some(),
        )?;
        if !condition_estimate.is_finite()
            || condition_estimate > controls.maximum_condition_estimate
        {
            return Err(ForcedSensitivityError::IllConditioned {
                estimate: condition_estimate,
                limit: controls.maximum_condition_estimate,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            pack,
            coefficients,
            omega: orbit.omega,
            coarse,
            fine,
            coarse_jacobian,
            fine_jacobian,
            normalized,
            factor,
            scales: admitted,
            phase,
            amplitude,
            controls,
            base_residual: coarse_base.max(fine_base),
            tangent_sampling_defect,
            condition_estimate,
        })
    }

    /// Number of complete original state coordinates.
    pub fn dimension(&self) -> usize {
        self.pack.d
    }
    /// Retained positive harmonic order.
    pub fn harmonics(&self) -> usize {
        self.pack.n_harm
    }
    /// Original forcing clock in radians per unit time.
    pub fn angular_frequency(&self) -> f64 {
        self.omega
    }

    /// Solve an exact COMPLETE residual-parameter column using the cached factor.
    ///
    /// The callback follows [`crate::ForcedParameterJacobian`]'s convention:
    /// hold every Fourier coefficient fixed and differentiate all original source,
    /// frequency-port and clock dependencies, in physical residual/parameter units.
    /// It is called once on each declared AFT grid. An existing family's exact
    /// callback can be passed as `|point,out| callback(&family,point,out)`.
    /// No parameter probes or additional nonlinear orbit solves are performed.
    pub fn solve(
        &self,
        direction: ForcedSensitivityDirection,
        residual_direction: impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        self.solve_with_amplitude(direction, 0.0, residual_direction, cancelled)
    }

    fn solve_with_amplitude(
        &self,
        direction: ForcedSensitivityDirection,
        amplitude_derivative: f64,
        mut residual_direction: impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        admit_direction(direction)?;
        let mut columns = Vec::with_capacity(2);
        for samples in [self.coarse.samples, self.fine.samples] {
            poll(&mut cancelled)?;
            let mut column = vec![f64::NAN; self.pack.coeff_len()];
            let point = ForcedParameterPoint {
                parameter: direction.parameter,
                angular_frequency: self.omega,
                harmonics: self.pack.n_harm,
                aft_samples: samples,
                coefficients: &self.coefficients[..self.pack.coeff_len()],
            };
            if !residual_direction(point, &mut column) || column.iter().any(|v| !v.is_finite()) {
                return Err(ForcedSensitivityError::DirectionRefused { samples });
            }
            columns.push(column);
        }
        self.solve_columns(
            direction,
            amplitude_derivative,
            &columns[0],
            &columns[1],
            &mut cancelled,
        )
    }

    /// Differentiate ONLY the original time-local island at a FIXED clock and
    /// fixed linear frequency port. `island_direction(t,x,out)` writes the full
    /// partial derivative at the unchanged original physical state. It must
    /// include every source/material-memory row and must not include `J*x_p`.
    ///
    /// Cached original waveforms and the SAME centered AFT projection form `R_p`;
    /// no baseline source or Hessian callbacks are repeated. Nonzero authored clock
    /// derivatives refuse here; use [`Self::solve`] for that complete derivative.
    pub fn solve_island_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        self.solve_island_with_amplitude(direction, 0.0, island_direction, cancelled)
    }

    fn solve_island_with_amplitude(
        &mut self,
        direction: ForcedSensitivityDirection,
        amplitude_derivative: f64,
        mut island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        admit_direction(direction)?;
        if direction.angular_frequency_derivative != 0.0 {
            return Err(invalid(
                "island-only direction requires a fixed physical clock",
            ));
        }
        let coarse = project_direction(
            &self.pack,
            self.omega,
            &mut self.coarse,
            &mut island_direction,
            &mut cancelled,
        )?;
        let fine = project_direction(
            &self.pack,
            self.omega,
            &mut self.fine,
            &mut island_direction,
            &mut cancelled,
        )?;
        self.solve_columns(
            direction,
            amplitude_derivative,
            &coarse,
            &fine,
            &mut cancelled,
        )
    }

    fn solve_columns(
        &self,
        direction: ForcedSensitivityDirection,
        amplitude_derivative: f64,
        coarse: &[f64],
        fine: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        poll(cancelled)?;
        let nc = self.pack.coeff_len();
        let n = self.scales.columns.len();
        let rows = n + usize::from(self.amplitude.is_some());
        let mut rhs = vec![0.0; rows];
        for row in 0..nc {
            rhs[row] = scaled(-fine[row], direction.parameter_scale, self.scales.rows[row])?;
        }
        if let Some(amplitude) = &self.amplitude {
            rhs[n] = amplitude.rhs(amplitude_derivative, direction.parameter_scale)?;
        }
        let mut solution = self.factor.solve(&rhs)?;
        let mut correction = vec![0.0; rows];
        let mut relative_residual = f64::INFINITY;
        // Reuse the same original factor for bounded residual correction.
        for iteration in 0..=3 {
            poll(cancelled)?;
            if solution.iter().any(|v| !v.is_finite()) {
                return Err(invalid("finite differentiated coefficients required"));
            }
            relative_residual = linear_residual(&self.normalized, &solution, &rhs, &mut correction);
            if relative_residual <= self.controls.maximum_relative_residual || iteration == 3 {
                break;
            }
            let change = self.factor.solve(&correction)?;
            for (value, change) in solution.iter_mut().zip(&change) {
                *value += change;
            }
        }
        let physical: Vec<_> = solution
            .iter()
            .zip(&self.scales.columns)
            .map(|(value, column)| scaled(*value, *column, direction.parameter_scale))
            .collect::<Result<_>>()?;
        let mut derivative_residual: f64 = 0.0;
        // Independently contract the ORIGINAL physical matrices and the rounded
        // published derivative, rather than just reassessing the normalized solve.
        for (matrix, column) in [(&self.coarse_jacobian, coarse), (&self.fine_jacobian, fine)] {
            let mut norm = 0.0;
            for row in 0..nc {
                poll(cancelled)?;
                let mut sum = Dd::from_f64(column[row]);
                for col in 0..n {
                    sum = sum + Dd::from_f64(matrix[row * n + col]) * Dd::from_f64(physical[col]);
                }
                norm = det::hypot(
                    norm,
                    scaled(
                        sum.to_f64(),
                        direction.parameter_scale,
                        self.scales.rows[row],
                    )?,
                );
            }
            derivative_residual = derivative_residual.max(norm);
        }
        if let Some(phase) = &self.phase {
            let mut defect = Dd::from_f64(0.0);
            // Check the published physical coefficients in the SAME fixed phase
            // chart. Normalizing first avoids unrepresentable reciprocal units.
            for col in 0..nc {
                let value = scaled(
                    physical[col],
                    direction.parameter_scale,
                    self.scales.columns[col],
                )?;
                defect = defect
                    + Dd::from_f64(phase.scaled_derivative(col, self.scales.columns[col]))
                        * Dd::from_f64(value);
            }
            derivative_residual = det::hypot(derivative_residual, defect.to_f64());
        }
        if let Some(amplitude) = &self.amplitude {
            derivative_residual = det::hypot(
                derivative_residual,
                amplitude.derivative_defect(
                    &physical,
                    &self.scales.columns,
                    amplitude_derivative,
                    direction.parameter_scale,
                )?,
            );
        }
        if !relative_residual.is_finite()
            || !derivative_residual.is_finite()
            || relative_residual > self.controls.maximum_relative_residual
            || derivative_residual > self.controls.maximum_derivative_residual
        {
            return Err(ForcedSensitivityError::DerivativeUnresolved {
                relative_residual,
                derivative_residual,
            });
        }
        poll(cancelled)?;
        let (coeffs, frequency) = self.pack.unpack(&physical);
        Ok(ForcedOrbitSensitivity {
            omega: self.omega,
            angular_frequency_derivative: if self.pack.omega_slot.is_some() {
                frequency
            } else {
                direction.angular_frequency_derivative
            },
            coeffs,
            relative_residual,
            derivative_residual,
        })
    }
}

fn invalid(what: &'static str) -> ForcedSensitivityError {
    ForcedSensitivityError::Invalid(what)
}

// Keep the original residual owner while detecting incomplete source writes
// at EVERY sample; a stale row must not masquerade as a differentiated source.
struct CheckedSource<'a, P>(&'a P);
impl<P: OrbitProblem> OrbitProblem for CheckedSource<'_, P> {
    fn dim(&self) -> usize {
        self.0.dim()
    }
    fn autonomous(&self) -> bool {
        self.0.autonomous()
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        self.0.island(time, state, out);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        self.0.port(s)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(ForcedSensitivityError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn scaled(value: f64, numerator: f64, denominator: f64) -> Result<f64> {
    let result = (value / denominator) * numerator;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(invalid("physical sensitivity scaling is not representable"))
    } else {
        Ok(result)
    }
}
fn admit_direction(direction: ForcedSensitivityDirection) -> Result<()> {
    if !direction.parameter.is_finite()
        || !positive(direction.parameter_scale)
        || !direction.angular_frequency_derivative.is_finite()
    {
        Err(invalid(
            "finite parameter/clock direction and positive parameter scale required",
        ))
    } else {
        Ok(())
    }
}
fn admit<P: OrbitProblem>(
    problem: &P,
    orbit: &HbOrbit,
    controls: ForcedSensitivityControls,
    scales: &HbScales,
    autonomous: bool,
    backbone: bool,
) -> Result<(Packing, Vec<f64>, AdmittedScales)> {
    let d = problem.dim();
    let width = orbit.coeffs.first().map_or(0, Vec::len);
    if (backbone && !autonomous)
        || problem.autonomous() != autonomous
        || d == 0
        || orbit.coeffs.len() != d
        || width < 2
        || width > 257
        || orbit.coeffs.iter().any(|row| {
            row.len() != width
                || row[0].im != 0.0
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        })
        || !positive(orbit.omega)
        || !positive(TAU / orbit.omega)
        || !positive((width - 1) as f64 * orbit.omega)
        || controls.samples_per_period < 4 * (width - 1) + 4
        || controls
            .samples_per_period
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .is_none()
        || !positive(controls.maximum_base_residual)
        || !positive(controls.maximum_tangent_sampling_defect)
        || !positive(controls.maximum_condition_estimate)
        || controls.maximum_condition_estimate < 1.0
        || !positive(controls.maximum_relative_residual)
        || controls.maximum_relative_residual >= 1.0
        || !positive(controls.maximum_derivative_residual)
    {
        return Err(invalid(
            "complete finite orbit, matching clock kind and valid sensitivity controls required",
        ));
    }
    let nc = d
        .checked_mul(2 * width - 1)
        .ok_or_else(|| invalid("coefficient count overflows"))?;
    let n = nc
        .checked_add(usize::from(autonomous))
        .ok_or_else(|| invalid("sensitivity dimension overflows"))?;
    let bounds = [(
        "real sensitivity dimension",
        n as u128,
        controls.maximum_unknowns.min(1024) as u128,
    )];
    for (quantity, requested, cap) in bounds {
        if requested > cap {
            return Err(ForcedSensitivityError::Budget {
                quantity,
                requested,
                cap,
            });
        }
    }
    let (n, d, m, h) = (
        n as u128,
        d as u128,
        controls.samples_per_period as u128,
        (width - 1) as u128,
    );
    // Includes both AFT grids/bases, tangent caches, all retained raw/scaled
    // matrices, LU copies/trailing scratch, condition estimation and one result.
    // GEMM's current fixed panels occupy (128+512)*256 f64 values.
    let transform_payload = crate::transform::payload(m as usize, h as usize, d as usize)? as u128
        + crate::transform::payload((2 * m + 1) as usize, h as usize, d as usize)? as u128;
    let rows = n + u128::from(backbone);
    let payload = transform_payload
        + 8 * (8 * rows * n
            + 16 * n
            + 4 * n * u128::from(backbone)
            + 16 * d * d * (2 * h + 1)
            + 8 * (3 * m + 1) * (d + 1)
            + 163_840);
    let work = 16 * rows * n * n + 128 * rows * n + 128 * (3 * m + 1) * (d * d + d) * (2 * h + 1);
    for (quantity, requested, cap) in [
        (
            "numerical payload bytes",
            payload,
            controls.maximum_payload_bytes as u128,
        ),
        ("numerical work", work, controls.maximum_work as u128),
    ] {
        if requested > cap {
            return Err(ForcedSensitivityError::Budget {
                quantity,
                requested,
                cap,
            });
        }
    }
    let pack = Packing {
        d: d as usize,
        n_harm: h as usize,
        omega_slot: autonomous.then_some(nc),
    };
    let anchor = if autonomous {
        HbAnchor::Autonomous {
            omega_guess: orbit.omega,
        }
    } else {
        HbAnchor::Forced { omega: orbit.omega }
    };
    let coefficients = orbit.packed_seed(&anchor);
    let mut admitted = AdmittedScales::new(
        pack.d,
        anchor,
        &coefficients,
        &HbBudget {
            harmonics: pack.n_harm,
            max_newton: 1,
            tolerance: 0.5,
            deflate_equilibrium: false,
        },
        scales,
    )?;
    if autonomous {
        admitted.rows.push(1.0);
    }
    Ok((pack, coefficients, admitted))
}

fn project_direction(
    pack: &Packing,
    omega: f64,
    aft: &mut AftWorkspace,
    callback: &mut impl FnMut(f64, &[f64], &mut [f64]) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    for sample in 0..aft.samples {
        poll(cancelled)?;
        aft.island.fill(f64::NAN);
        let time = (TAU / omega) * sample as f64 / aft.samples as f64;
        if !callback(
            time,
            &aft.xs[sample * pack.d..(sample + 1) * pack.d],
            &mut aft.island,
        ) || aft.island.iter().any(|v| !v.is_finite())
        {
            return Err(ForcedSensitivityError::IslandDirectionRefused {
                sample,
                samples: aft.samples,
            });
        }
        aft.gs[sample * pack.d..(sample + 1) * pack.d].copy_from_slice(&aft.island);
    }
    aft.center_island(pack.d);
    aft.project_island(pack);
    let mut column = vec![0.0; pack.coeff_len()];
    for n in 0..=pack.n_harm {
        poll(cancelled)?;
        for state in 0..pack.d {
            let coefficient = aft.island_coefficient(pack, n, state);
            column[pack.re_index(n, state)] = -coefficient.re;
            if n > 0 {
                column[pack.im_index(n, state)] = -coefficient.im;
            }
        }
    }
    if column.iter().any(|v| !v.is_finite()) {
        return Err(invalid("finite projected source direction required"));
    }
    Ok(column)
}

fn linear_residual(matrix: &[f64], solution: &[f64], rhs: &[f64], correction: &mut [f64]) -> f64 {
    let n = solution.len();
    let mut maximum: f64 = 0.0;
    let mut maximum_defect: f64 = 0.0;
    let mut maximum_action: f64 = 0.0;
    for row in 0..rhs.len() {
        let mut sum = Dd::from_f64(-rhs[row]);
        let mut scale = Dd::from_f64(rhs[row].abs());
        for col in 0..n {
            let term = Dd::from_f64(matrix[row * n + col]) * Dd::from_f64(solution[col]);
            sum = sum + term;
            scale = scale + Dd::from_f64(term.to_f64().abs());
        }
        let remainder = sum.to_f64();
        let scale = scale.to_f64();
        if !remainder.is_finite() || !scale.is_finite() {
            return f64::INFINITY;
        }
        correction[row] = -remainder;
        maximum_defect = maximum_defect.max(remainder.abs());
        maximum_action = maximum_action.max(scale);
        let relative = if scale == 0.0 && remainder == 0.0 {
            0.0
        } else {
            remainder.abs() / scale
        };
        if !relative.is_finite() {
            return f64::INFINITY;
        }
        maximum = maximum.max(relative);
    }
    if rhs.len() == n {
        maximum
    } else if maximum_action == 0.0 {
        0.0
    } else {
        // A conservative rectangular system has a redundant balance row.
        // Analytically zero harmonic actions can carry only quadrature noise;
        // dividing that row's defect by that same noise is not a useful solve
        // diagnostic. Its absolute physical equation defect remains separately
        // enforced on BOTH grids, including the rounded published derivative.
        maximum_defect / maximum_action
    }
}
