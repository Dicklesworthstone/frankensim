//! Multiplicative, isochoric J2 plasticity with quadratic Hencky elasticity.
//!
//! `F = Fe Fp`, `Cp = Fp^T Fp`, and the retained history is `Cp^-1, alpha`.
//! The yield surface uses **Kirchhoff** stress, not Cauchy stress:
//! `q = sqrt(3/2) ||dev(tau)|| = sigma_y + H alpha`. All energies and
//! dissipation are per reference volume. The exact scalar exponential return
//! is `dgamma = max(q_trial-sigma_y-H alpha_old,0)/(3 mu+H)`.
//!
//! This specializes the logarithmic-elastic multiplicative return-map family
//! of Simo, CMAME 99 (1992), 61–112, doi:10.1016/0045-7825(92)90123-2.
//! The tangent differentiates the discrete update through the Frechet
//! derivative of the matrix logarithm, including coincident eigenvalues.
//! It is neither a rotated infinitesimal-strain update nor a finite-difference
//! tangent. Large rotation and plastic stretch are supported within explicit
//! elastic, plastic and accumulated-history bounds. This law is isothermal,
//! rate independent and isotropically hardening; it does not assert cyclic
//! kinematic hardening, damage, thermal expansion or a mixed-pressure element.

use core::fmt;
use fs_math::det;

/// Explicit-time Perzyna kinetics on this same finite return and tangent owner.
pub mod rate;
/// Strict, explicit model-card binding for this exact Kirchhoff-yield law.
pub mod source;
/// Affine temperature dependence and multiplicative isotropic thermal stretch.
pub mod thermal;

const IDENTITY6: [f64; 6] = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
const ROUND_OFF: f64 = 2048.0 * f64::EPSILON;

/// Independently specified, fixed-temperature coefficients in Pa.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2Parameters {
    /// Elastic shear modulus, positive and finite.
    pub shear_modulus_pa: f64,
    /// Elastic bulk modulus, positive and finite.
    pub bulk_modulus_pa: f64,
    /// Initial yield radius in Kirchhoff equivalent stress.
    pub yield_stress_pa: f64,
    /// Energetic linear isotropic modulus: stored hardening energy is H alpha²/2.
    pub hardening_modulus_pa: f64,
}

/// Explicit physical applicability bounds; there are no material defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2Limits {
    /// Strict lower bound on det(F); must lie strictly between zero and one.
    pub minimum_jacobian: f64,
    /// Maximum Frobenius norm of the **returned** elastic logarithmic strain.
    /// A plastic predictor may exceed this bound before a valid return.
    pub maximum_elastic_log_strain: f64,
    /// Maximum Frobenius norm of `0.5 log(Cp^-1)`.
    pub maximum_plastic_log_strain: f64,
    /// Bounds accumulated flow even when reversal restores a small net metric.
    pub maximum_accumulated_plastic_strain: f64,
}

/// Immutable law; coefficients and domain cannot change underneath a history.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2 {
    parameters: FiniteJ2Parameters,
    limits: FiniteJ2Limits,
}

/// Complete native material-point history, with no hidden deformation reset.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2State {
    /// Symmetric `Cp^-1`, tensor-Voigt order `[xx, yy, zz, xy, yz, zx]`.
    /// It must be positive definite and isochoric, up to numerical roundoff.
    pub inverse_plastic_cauchy_green: [f64; 6],
    /// Accumulated equivalent plastic strain; nonnegative and nondecreasing.
    pub alpha: f64,
}

/// An immutable candidate: the caller commits `state` only after step acceptance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2Response {
    /// Complete candidate history, committed only when the surrounding step succeeds.
    pub state: FiniteJ2State,
    /// Work conjugate to F, row-major tensor components, per reference area.
    pub first_piola_pa: [f64; 9],
    /// Kirchhoff stress `tau = J sigma`, in row-major tensor components.
    pub kirchhoff_pa: [f64; 9],
    /// Exact dP/dF at frozen old history, indexed `[P_component][F_component]`.
    /// At the yield kink the numerical zero-overstress branch is elastic.
    pub tangent_pa: [[f64; 9]; 9],
    /// Reduced incremental potential: stored energy + sigma_y dgamma for the
    /// rate-independent law. Explicit kinetics also adds its integrated convex
    /// rate potential, which is distinct from physical viscous dissipation.
    /// Its derivative is P; it is not the stored energy of the candidate.
    pub incremental_potential_j_m3: f64,
    /// Actual elastic plus energetic-hardening storage at the returned endpoint.
    pub stored_energy_j_m3: f64,
    /// Physical dissipation only: sigma_y dgamma plus endpoint viscous work
    /// when kinetics is explicitly selected, without numerical work defect.
    pub plastic_dissipation_j_m3: f64,
    /// Nonnegative increment of accumulated equivalent plastic strain.
    pub plastic_increment: f64,
    /// Returned `sqrt(3/2) ||dev(tau)||`, independent of volumetric stress.
    pub equivalent_kirchhoff_pa: f64,
    /// Current-to-reference volume ratio `det(F)`.
    pub jacobian: f64,
}

/// Allows a strict source-resolved owner to expose the same immutable native law.
pub trait FiniteJ2Material {
    /// Borrow the admitted native coefficients and applicability bounds.
    fn finite_j2(&self) -> &FiniteJ2;
}

impl FiniteJ2Material for FiniteJ2 {
    fn finite_j2(&self) -> &FiniteJ2 {
        self
    }
}

impl<T: FiniteJ2Material + ?Sized> FiniteJ2Material for &T {
    fn finite_j2(&self) -> &FiniteJ2 {
        (*self).finite_j2()
    }
}

/// Rejections leave the borrowed accepted state untouched.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteJ2Error {
    /// Invalid constitutive coefficients, domain bounds or caller tolerance.
    Parameters(&'static str),
    /// Invalid or unreachable retained material history or stretch metric.
    State(&'static str),
    /// A declared strain or Jacobian bound was exceeded.
    Domain {
        /// Quantity and, for the Jacobian, lower-bound direction.
        quantity: &'static str,
        /// Rejected value.
        value: f64,
        /// Declared applicability bound.
        limit: f64,
    },
    /// An inspected accepted endpoint lies outside the yield surface.
    YieldResidual {
        /// Equivalent Kirchhoff stress minus the current yield radius.
        residual_pa: f64,
        /// Explicit caller tolerance for the accepted-state check.
        tolerance_pa: f64,
    },
    /// The checked update cannot be represented with finite numerical outputs.
    Numerical(&'static str),
}

impl fmt::Display for FiniteJ2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameters(what) => write!(f, "finite J2 parameters: {what}"),
            Self::State(what) => write!(f, "finite J2 history: {what}"),
            Self::Domain {
                quantity,
                value,
                limit,
            } => {
                write!(f, "finite J2 domain: {quantity}={value}, bound={limit}")
            }
            Self::YieldResidual {
                residual_pa,
                tolerance_pa,
            } => write!(
                f,
                "finite J2 accepted-state yield residual {residual_pa} Pa exceeds {tolerance_pa} Pa"
            ),
            Self::Numerical(what) => write!(f, "finite J2 numerical refusal: {what}"),
        }
    }
}

impl std::error::Error for FiniteJ2Error {}

impl FiniteJ2 {
    /// Check explicit coefficients and bounds; no empirical coefficients inferred.
    pub fn new(
        parameters: FiniteJ2Parameters,
        limits: FiniteJ2Limits,
    ) -> Result<Self, FiniteJ2Error> {
        let p = &parameters;
        if [p.shear_modulus_pa, p.bulk_modulus_pa, p.yield_stress_pa]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !p.hardening_modulus_pa.is_finite()
            || p.hardening_modulus_pa < 0.0
            || !(3.0 * p.shear_modulus_pa + p.hardening_modulus_pa).is_finite()
        {
            return Err(FiniteJ2Error::Parameters(
                "finite positive mu, kappa and yield stress, nonnegative H, and representable 3mu+H required",
            ));
        }
        if !limits.minimum_jacobian.is_finite()
            || limits.minimum_jacobian <= 0.0
            || limits.minimum_jacobian >= 1.0
            || [
                limits.maximum_elastic_log_strain,
                limits.maximum_plastic_log_strain,
                limits.maximum_accumulated_plastic_strain,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(FiniteJ2Error::Parameters(
                "0 < minimum_jacobian < 1 and finite positive strain limits required",
            ));
        }
        if !(p.yield_stress_pa + p.hardening_modulus_pa * limits.maximum_accumulated_plastic_strain)
            .is_finite()
        {
            return Err(FiniteJ2Error::Parameters(
                "yield radius must remain representable throughout the accumulated-strain domain",
            ));
        }
        Ok(Self { parameters, limits })
    }

    /// Original immutable constitutive coefficients.
    #[must_use]
    pub fn parameters(&self) -> &FiniteJ2Parameters {
        &self.parameters
    }

    /// Original immutable applicability bounds.
    #[must_use]
    pub fn limits(&self) -> &FiniteJ2Limits {
        &self.limits
    }

    /// Virgin history: identity plastic metric and no accumulated plastic flow.
    #[must_use]
    pub fn initial_state(&self) -> FiniteJ2State {
        FiniteJ2State {
            inverse_plastic_cauchy_green: IDENTITY6,
            alpha: 0.0,
        }
    }

    /// Validate SPD, plastic incompressibility, applicability and reachability.
    /// The history-distance bound follows from the triangle inequality of the
    /// congruence-invariant SPD metric: each exponential flow increment has
    /// distance `sqrt(3/2) dgamma` in the half-log convention.
    pub fn validate_history(&self, state: &FiniteJ2State) -> Result<(), FiniteJ2Error> {
        if !state.alpha.is_finite() || state.alpha < 0.0 {
            return Err(FiniteJ2Error::State("alpha must be finite and nonnegative"));
        }
        upper_bound(
            "accumulated plastic strain",
            state.alpha,
            self.limits.maximum_accumulated_plastic_strain,
        )?;
        // The exact identity needs neither a spectral decomposition nor a
        // reachability test, including a closed plastic path with alpha > 0.
        if state.inverse_plastic_cauchy_green == IDENTITY6 {
            return Ok(());
        }
        if state.alpha == 0.0 {
            return Err(FiniteJ2Error::State(
                "zero accumulated flow requires the virgin plastic metric",
            ));
        }
        let spectrum = Spectrum::spd(&expand(&state.inverse_plastic_cauchy_green))?;
        let strain = spectrum.log_values.map(|x| 0.5 * x);
        let magnitude = norm(&strain);
        if spectrum.log_values.iter().sum::<f64>().abs() > ROUND_OFF * (1.0 + magnitude) {
            return Err(FiniteJ2Error::State(
                "plastic metric must have determinant one",
            ));
        }
        upper_bound(
            "plastic logarithmic strain",
            magnitude,
            self.limits.maximum_plastic_log_strain,
        )?;
        if magnitude > det::sqrt(1.5) * state.alpha + ROUND_OFF * (1.0 + magnitude + state.alpha) {
            return Err(FiniteJ2Error::State(
                "plastic metric is unreachable from the retained accumulated flow",
            ));
        }
        Ok(())
    }

    /// Inspect the actual stored endpoint energy, without a second plastic return.
    pub fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        let endpoint = self.elastic_endpoint(f, state)?;
        self.check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
        self.energy(endpoint.log_j, &endpoint.deviator, state.alpha)
    }

    /// Check that an accepted deformation/history is inside the yield surface.
    /// This never silently returns or repairs an inconsistent checkpoint.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "yield tolerance must be finite and nonnegative",
            ));
        }
        let endpoint = self.elastic_endpoint(f, state)?;
        self.check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
        self.energy(endpoint.log_j, &endpoint.deviator, state.alpha)?;
        let q = det::sqrt(6.0) * self.parameters.shear_modulus_pa * norm(&endpoint.deviator);
        let residual = q - self.yield_radius(state.alpha);
        if !residual.is_finite() {
            return Err(FiniteJ2Error::Numerical(
                "accepted yield residual is not finite",
            ));
        }
        if residual > yield_tolerance_pa {
            return Err(FiniteJ2Error::YieldResidual {
                residual_pa: residual,
                tolerance_pa: yield_tolerance_pa,
            });
        }
        Ok(())
    }

    /// Frozen-history exponential return and analytic, full 9-by-9 tangent.
    /// Caller-owned input history is unchanged on success, rejection or retry.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        let returned = self.return_map(f, old)?;
        self.response_from_return(f, old, returned)
    }

    fn response_from_return(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
        returned: ReturnMap,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        let tangent = self.return_tangent(f, old, &returned)?;
        Ok(FiniteJ2Response {
            state: returned.state,
            first_piola_pa: returned.first_piola,
            kirchhoff_pa: returned.tau,
            tangent_pa: tangent,
            incremental_potential_j_m3: returned.potential,
            stored_energy_j_m3: returned.stored,
            plastic_dissipation_j_m3: returned.dissipation,
            plastic_increment: returned.dgamma,
            equivalent_kirchhoff_pa: returned.scale * returned.q_trial,
            jacobian: returned.trial.jacobian,
        })
    }

    /// One checked Piola stress and reduced potential evaluation, without
    /// assembling dP/dF. Newton line searches use the same return map and
    /// retained/returned-history domain checks as [`Self::response`].
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        let returned = self.return_map(f, old)?;
        Ok((returned.first_piola, returned.potential))
    }

    fn return_map(&self, f: &[f64; 9], old: &FiniteJ2State) -> Result<ReturnMap, FiniteJ2Error> {
        self.return_map_with_rate(f, old, None)
    }

    fn return_map_with_rate(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
        kinetics: Option<&rate::FiniteOverstressJ2>,
    ) -> Result<ReturnMap, FiniteJ2Error> {
        let trial = self.elastic_endpoint(f, old)?;
        let p = &self.parameters;
        let trial_stress = trial.deviator.map(|x| 2.0 * p.shear_modulus_pa * x);
        let q_trial = det::sqrt(1.5) * norm(&trial_stress);
        if !q_trial.is_finite() {
            return Err(FiniteJ2Error::Numerical(
                "trial equivalent stress is not finite",
            ));
        }
        let radius_old = self.yield_radius(old.alpha);
        let excess = q_trial - radius_old;
        // Reconstructing an accepted metric and taking its logarithm can leave
        // a stress-scale roundoff residual at yield. Do not turn an exact hold
        // or rigid rotation into a fictitious retained plastic increment.
        // This numerical zero band is not a user-supplied overstress law.
        let yield_roundoff =
            64.0 * f64::EPSILON * (2.0 * p.shear_modulus_pa).max(q_trial).max(radius_old);
        let plastic = excess > yield_roundoff;
        let denominator = 3.0 * p.shear_modulus_pa + p.hardening_modulus_pa;
        let flow = if plastic {
            kinetics
                .map(|law| law.solve(excess, denominator))
                .transpose()?
        } else {
            None
        };
        let dgamma = if let Some(flow) = &flow {
            flow.increment
        } else if plastic {
            excess / denominator
        } else {
            0.0
        };
        let alpha = old.alpha + dgamma;
        if plastic && (dgamma == 0.0 || alpha == old.alpha) {
            return Err(FiniteJ2Error::Numerical(
                "plastic increment is below retained-history resolution",
            ));
        }
        upper_bound(
            "accumulated plastic strain",
            alpha,
            self.limits.maximum_accumulated_plastic_strain,
        )?;
        // This expression avoids subtracting nearly equal values for a large
        // plastic step or the perfectly plastic limit H=0.
        let scale = if let Some(flow) = &flow {
            (self.yield_radius(alpha) + flow.overstress_pa) / q_trial
        } else if plastic {
            self.yield_radius(alpha) / q_trial
        } else {
            1.0
        };
        let returned_deviator = trial.deviator.map(|x| scale * x);
        self.check_elastic_domain(trial.log_j, &returned_deviator)?;
        let mut tau = trial_stress.map(|x| scale * x);
        for i in [0, 4, 8] {
            tau[i] += p.bulk_modulus_pa * trial.log_j;
        }
        let inverse = inverse(f, trial.jacobian)?;
        let inverse_transpose = transpose(&inverse);
        let first_piola = multiply(&tau, &inverse_transpose);
        let state = if plastic {
            let mean = trial.spectrum.log_values.iter().sum::<f64>() / 6.0;
            let elastic_roots = trial
                .spectrum
                .log_values
                .map(|x| det::exp(2.0 * (trial.log_j / 3.0 + scale * (0.5 * x - mean))));
            if elastic_roots.iter().any(|x| !x.is_finite() || *x <= 0.0) {
                return Err(FiniteJ2Error::Numerical(
                    "returned elastic metric is not representable",
                ));
            }
            let be = trial.spectrum.reconstruct(&elastic_roots);
            let cp_inv = multiply(&multiply(&inverse, &be), &inverse_transpose);
            let candidate = FiniteJ2State {
                inverse_plastic_cauchy_green: pack(&cp_inv),
                alpha,
            };
            self.validate_history(&candidate)?;
            candidate
        } else {
            // Elastic and rigid increments preserve the accepted metric bit for
            // bit; unnecessary spectral reconstruction would drift the history.
            old.clone()
        };
        let stored = self.energy(trial.log_j, &returned_deviator, alpha)?;
        let yield_dissipation = p.yield_stress_pa * dgamma;
        let (dissipation, potential) = if let Some(flow) = &flow {
            let viscous_work = flow.overstress_pa * dgamma;
            (
                yield_dissipation + viscous_work,
                stored + yield_dissipation + flow.potential_fraction * viscous_work,
            )
        } else {
            (yield_dissipation, stored + yield_dissipation)
        };
        // Positive stored/rate contributions avoid cancellation after large
        // plastic steps; physical viscous work remains separate from its integral.
        if tau.iter().chain(first_piola.iter()).any(|x| !x.is_finite())
            || !potential.is_finite()
            || !dissipation.is_finite()
            || (plastic && dissipation == 0.0)
        {
            return Err(FiniteJ2Error::Numerical(
                "stress, potential or physical dissipation is not representable",
            ));
        }
        Ok(ReturnMap {
            trial,
            state,
            first_piola,
            tau,
            inverse_transpose,
            scale,
            q_trial,
            dgamma,
            stored,
            dissipation,
            potential,
            flow,
        })
    }

    fn return_tangent(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
        returned: &ReturnMap,
    ) -> Result<[[f64; 9]; 9], FiniteJ2Error> {
        let p = &self.parameters;
        let trial = &returned.trial;
        let inverse_transpose = &returned.inverse_transpose;
        let first_piola = &returned.first_piola;
        let scale = returned.scale;
        let q_trial = returned.q_trial;
        let plastic = returned.dgamma > 0.0;
        let radius_old = self.yield_radius(old.alpha);
        let denominator = 3.0 * p.shear_modulus_pa + p.hardening_modulus_pa;
        let trial_stress = trial.deviator.map(|x| 2.0 * p.shear_modulus_pa * x);
        let cp_inv = expand(&old.inverse_plastic_cauchy_green);
        let f_cp = multiply(f, &cp_inv);
        let log_derivative = trial.spectrum.divided_differences();
        let mut tangent = [[0.0; 9]; 9];
        for column in 0..9 {
            let row = column / 3;
            let col = column % 3;
            let db: [f64; 9] = core::array::from_fn(|ij| {
                let (i, j) = (ij / 3, ij % 3);
                (if i == row { f_cp[3 * j + col] } else { 0.0 })
                    + (if j == row { f_cp[3 * i + col] } else { 0.0 })
            });
            let dlog = trial.spectrum.frechet(&db, &log_derivative);
            let dtrial = deviator(&dlog).map(|x| p.shear_modulus_pa * x);
            let mut dtau = dtrial.map(|x| scale * x);
            if plastic {
                // n = s_trial/q avoids q² overflow and a norm derivative at q=0.
                let direction = trial_stress.map(|x| x / q_trial);
                let dq = 1.5 * dot(&direction, &dtrial);
                let coefficient = if let Some(flow) = &returned.flow {
                    // Differentiate the implicit scalar flow, not its secant
                    // approximation. The stable positive numerator avoids
                    // cancellation in 1/(A+f') - dgamma/q near perfect flow.
                    -(3.0 * p.shear_modulus_pa / (denominator + flow.tangent_modulus_pa))
                        * ((radius_old + flow.transverse_overstress_pa) / q_trial)
                        * dq
                } else {
                    -(3.0 * p.shear_modulus_pa / denominator) * (radius_old / q_trial) * dq
                };
                for i in 0..9 {
                    dtau[i] += coefficient * direction[i];
                }
            }
            for i in [0, 4, 8] {
                dtau[i] += p.bulk_modulus_pa * inverse_transpose[column];
            }
            let dp = multiply(&dtau, inverse_transpose);
            for ij in 0..9 {
                // dP = dtau F^-T - P dF^T F^-T.
                tangent[ij][column] =
                    dp[ij] - first_piola[3 * (ij / 3) + col] * inverse_transpose[3 * row + ij % 3];
            }
        }
        if tangent.iter().flatten().any(|x| !x.is_finite()) {
            return Err(FiniteJ2Error::Numerical(
                "consistent tangent is not representable",
            ));
        }
        Ok(tangent)
    }

    fn yield_radius(&self, alpha: f64) -> f64 {
        self.parameters.yield_stress_pa + self.parameters.hardening_modulus_pa * alpha
    }

    fn check_elastic_domain(&self, log_j: f64, dev: &[f64; 9]) -> Result<(), FiniteJ2Error> {
        let magnitude = norm(&[norm(dev), log_j / det::sqrt(3.0)]);
        upper_bound(
            "elastic logarithmic strain",
            magnitude,
            self.limits.maximum_elastic_log_strain,
        )
    }

    fn energy(&self, log_j: f64, dev: &[f64; 9], alpha: f64) -> Result<f64, FiniteJ2Error> {
        let p = &self.parameters;
        let energy = 0.5 * p.bulk_modulus_pa * log_j * log_j
            + p.shear_modulus_pa * dot(dev, dev)
            + 0.5 * p.hardening_modulus_pa * alpha * alpha;
        if !energy.is_finite() {
            return Err(FiniteJ2Error::Numerical("stored energy is not finite"));
        }
        Ok(energy)
    }

    fn elastic_endpoint(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
    ) -> Result<ElasticEndpoint, FiniteJ2Error> {
        self.validate_history(state)?;
        if f.iter().any(|x| !x.is_finite()) {
            return Err(FiniteJ2Error::Numerical("F must be finite"));
        }
        let jacobian = determinant(f);
        if !jacobian.is_finite() || jacobian <= self.limits.minimum_jacobian {
            return Err(FiniteJ2Error::Domain {
                quantity: "det(F), strict lower bound",
                value: jacobian,
                limit: self.limits.minimum_jacobian,
            });
        }
        let cp_inv = expand(&state.inverse_plastic_cauchy_green);
        let be = multiply(&multiply(f, &cp_inv), &transpose(f));
        let spectrum = Spectrum::spd(&be)?;
        let logarithm = spectrum.reconstruct(&spectrum.log_values);
        let deviator = deviator(&logarithm).map(|x| 0.5 * x);
        Ok(ElasticEndpoint {
            jacobian,
            log_j: det::ln(jacobian),
            spectrum,
            deviator,
        })
    }
}

struct ElasticEndpoint {
    jacobian: f64,
    log_j: f64,
    spectrum: Spectrum,
    deviator: [f64; 9],
}

/// Checked constitutive candidate, with no tangent storage or evaluation.
/// The spectral trial data are reused only when the caller requests dP/dF.
struct ReturnMap {
    trial: ElasticEndpoint,
    state: FiniteJ2State,
    first_piola: [f64; 9],
    tau: [f64; 9],
    inverse_transpose: [f64; 9],
    scale: f64,
    q_trial: f64,
    dgamma: f64,
    stored: f64,
    dissipation: f64,
    potential: f64,
    flow: Option<rate::FlowReturn>,
}

/// Fixed 3-by-3 data; only the existing native eigensolver owns a small workspace.
struct Spectrum {
    values: [f64; 3],
    log_values: [f64; 3],
    vectors: [f64; 9],
}

impl Spectrum {
    fn spd(matrix: &[f64; 9]) -> Result<Self, FiniteJ2Error> {
        if matrix.iter().any(|x| !x.is_finite()) {
            return Err(FiniteJ2Error::Numerical(
                "elastic or plastic metric must be finite",
            ));
        }
        let (values, vectors) = fs_la::eigen::jacobi_eigh(matrix, 3);
        if values.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || values[0] <= 128.0 * f64::EPSILON * values[2]
        {
            return Err(FiniteJ2Error::State(
                "elastic or plastic metric must be numerically resolved positive definite",
            ));
        }
        Ok(Self {
            values: [values[0], values[1], values[2]],
            log_values: [det::ln(values[0]), det::ln(values[1]), det::ln(values[2])],
            vectors: core::array::from_fn(|i| vectors[i]),
        })
    }

    fn reconstruct(&self, diagonal: &[f64; 3]) -> [f64; 9] {
        core::array::from_fn(|ij| {
            (0..3)
                .map(|k| {
                    self.vectors[3 * (ij / 3) + k] * diagonal[k] * self.vectors[3 * (ij % 3) + k]
                })
                .sum()
        })
    }

    fn divided_differences(&self) -> [f64; 9] {
        core::array::from_fn(|ij| log_divided_difference(self.values[ij / 3], self.values[ij % 3]))
    }

    fn frechet(&self, direction: &[f64; 9], divided_differences: &[f64; 9]) -> [f64; 9] {
        let vt = transpose(&self.vectors);
        let rotated = multiply(&multiply(&vt, direction), &self.vectors);
        let weighted: [f64; 9] = core::array::from_fn(|i| rotated[i] * divided_differences[i]);
        multiply(&multiply(&self.vectors, &weighted), &vt)
    }
}

// Exact coalescence limit and a stable analytic atanh series close to it.
// With |r|<1/8, the first omitted term after r^20/21 is <6e-23.
fn log_divided_difference(a: f64, b: f64) -> f64 {
    if a == b {
        return 1.0 / a;
    }
    let sum = 0.5 * a + 0.5 * b;
    let r = (0.5 * a - 0.5 * b) / sum;
    if r.abs() < 0.125 {
        let squared = r * r;
        let mut power = 1.0;
        let mut series = 1.0;
        for n in 1..=10 {
            power *= squared;
            series += power / f64::from(2 * n + 1);
        }
        series / sum
    } else {
        (det::ln(a) - det::ln(b)) / (a - b)
    }
}

fn upper_bound(quantity: &'static str, value: f64, limit: f64) -> Result<(), FiniteJ2Error> {
    if !value.is_finite() || value > limit {
        return Err(FiniteJ2Error::Domain {
            quantity,
            value,
            limit,
        });
    }
    Ok(())
}

fn norm<const N: usize>(v: &[f64; N]) -> f64 {
    let scale = v.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    if !scale.is_finite() {
        return f64::INFINITY;
    }
    scale * det::sqrt(v.iter().map(|x| (x / scale) * (x / scale)).sum())
}

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn deviator(a: &[f64; 9]) -> [f64; 9] {
    let mean = (a[0] + a[4] + a[8]) / 3.0;
    core::array::from_fn(|i| a[i] - if i % 4 == 0 { mean } else { 0.0 })
}

fn expand(a: &[f64; 6]) -> [f64; 9] {
    [a[0], a[3], a[5], a[3], a[1], a[4], a[5], a[4], a[2]]
}

fn pack(a: &[f64; 9]) -> [f64; 6] {
    [
        a[0],
        a[4],
        a[8],
        0.5 * (a[1] + a[3]),
        0.5 * (a[5] + a[7]),
        0.5 * (a[6] + a[2]),
    ]
}

fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| a[3 * (i % 3) + i / 3])
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn determinant(a: &[f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}

fn inverse(a: &[f64; 9], jacobian: f64) -> Result<[f64; 9], FiniteJ2Error> {
    let cof = |r: usize, c: usize| {
        let (r0, r1) = ((r + 1) % 3, (r + 2) % 3);
        let (c0, c1) = ((c + 1) % 3, (c + 2) % 3);
        (a[3 * r0 + c0] * a[3 * r1 + c1] - a[3 * r0 + c1] * a[3 * r1 + c0]) / jacobian
    };
    let result = core::array::from_fn(|ij| cof(ij % 3, ij / 3));
    if result.iter().any(|x| !x.is_finite()) {
        return Err(FiniteJ2Error::Numerical("inverse F is not representable"));
    }
    Ok(result)
}
