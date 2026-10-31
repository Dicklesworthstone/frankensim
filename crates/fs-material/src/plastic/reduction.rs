//! Stress-constrained reductions of the native three-dimensional material law.
//!
//! A free-sided fiber has uniaxial STRESS, not uniaxial strain. Its transverse
//! contraction and incompressible plastic flow must survive the reduction.
//! The analytical axial restriction below is exact for isotropic J2 with linear
//! isotropic/Prager hardening. The checked plane-stress adapter solves the
//! transverse stress constraints against the original frozen history and
//! condenses the original algorithmic tangent after convergence.

use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use super::voce::VoceJ2;
use super::{J2Plasticity, J2State};
use crate::fiber::Uniaxial;
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;

/// Exact axial stress reduction of the native small-strain, rate-independent J2
/// law. The independent isotropic and Prager moduli remain distinct on reversal.
/// This is not a finite-strain rod, damage law or calibrated cyclic material.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialJ2 {
    backbone: J2Plasticity,
}

/// Committed axial plastic history. Private fields keep this history on the
/// uniaxial-flow subspace; general six-component histories use the 3-D owner.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialJ2State {
    plastic_strain: f64,
    equivalent_plastic_strain: f64,
}
impl UniaxialJ2State {
    /// Signed axial plastic strain; it can decrease during reversed plastic flow.
    pub const fn plastic_strain(&self) -> f64 {
        self.plastic_strain
    }
    /// Nondecreasing accumulated equivalent plastic strain.
    pub const fn equivalent_plastic_strain(&self) -> f64 {
        self.equivalent_plastic_strain
    }
    /// The original tensor-Voigt history, with trace-free lateral plastic flow.
    pub fn tensor_history(&self) -> J2State {
        let p = self.plastic_strain;
        J2State {
            plastic_strain: [p, -0.5 * p, -0.5 * p, 0.0, 0.0, 0.0],
            alpha: self.equivalent_plastic_strain,
        }
    }
}

/// One axial return with the full compatible strain and retained native history.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialJ2Response {
    /// Total tensor strain, including elastic contraction and lateral plasticity.
    pub strain: Voigt,
    /// Axial Cauchy stress in the declared small-strain model [Pa].
    pub stress_pa: f64,
    /// Exact derivative with respect to axial total strain at frozen old history.
    pub tangent_pa: f64,
    /// Candidate committed history; evaluating does not mutate accepted history.
    pub state: UniaxialJ2State,
}

impl UniaxialJ2 {
    /// Admit the same finite elastic/yield/hardening coefficients as the native
    /// energetic J2 owner. No additional material parameter is introduced.
    pub fn new(backbone: J2Plasticity) -> Result<Self, MaterialError> {
        backbone
            .validate_overstress_history(&[0.0; 6], &backbone.initial_state())
            .map_err(|e| MaterialError::Parameters {
                what: e.to_string(),
            })?;
        let denominator = backbone.elastic.youngs + backbone.combined_hardening();
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "uniaxial J2 return denominator must be finite and positive".into(),
            });
        }
        Ok(Self { backbone })
    }

    /// Original native coefficients; source-bound consumers retain their source
    /// bundle alongside this immutable law.
    pub const fn backbone(&self) -> &J2Plasticity {
        &self.backbone
    }

    /// Evaluate the exact restriction of associative J2 to free-sided axial
    /// loading. Here q=|sigma-H_k eps_p|, delta_alpha=|delta_eps_p| and the
    /// returned axial tangent is E(H_iso+H_k)/(E+H_iso+H_k).
    ///
    /// Like the existing `Uniaxial` trait this raw kernel is infallible; callers
    /// must admit finite inputs and the native small-strain domain. The checked
    /// energy entry point below applies the original full-tensor admission.
    pub fn response(&self, strain: f64, state: &UniaxialJ2State) -> UniaxialJ2Response {
        let e = self.backbone.elastic.youngs;
        let h = self.backbone.combined_hardening();
        let trial = e * (strain - state.plastic_strain)
            - self.backbone.kinematic_hardening * state.plastic_strain;
        let excess = trial.abs()
            - (self.backbone.yield_stress
                + self.backbone.hardening * state.equivalent_plastic_strain);
        let mut next = state.clone();
        let tangent_pa = if excess > 0.0 {
            let delta = excess / (e + h);
            next.plastic_strain += trial.signum() * delta;
            next.equivalent_plastic_strain += delta;
            e * (h / (e + h))
        } else {
            e
        };
        let elastic_strain = strain - next.plastic_strain;
        let transverse =
            -self.backbone.elastic.poisson * elastic_strain - 0.5 * next.plastic_strain;
        UniaxialJ2Response {
            strain: [strain, transverse, transverse, 0.0, 0.0, 0.0],
            stress_pa: e * elastic_strain,
            tangent_pa,
            state: next,
        }
    }

    /// Audit axial work, recoverable hardening and plastic heat through the
    /// original 3-D energetic owner. Its tensor result carries all lateral
    /// history; no independent heat fraction or energy formula is introduced.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialJ2State,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        let old_transverse = -self.backbone.elastic.poisson
            * (old_strain - old_state.plastic_strain)
            - 0.5 * old_state.plastic_strain;
        let old_tensor = [old_strain, old_transverse, old_transverse, 0.0, 0.0, 0.0];
        let next = self.response(strain, old_state);
        self.backbone.advance_with_energy(
            &old_tensor,
            &old_state.tensor_history(),
            &next.strain,
            tolerances,
        )
    }
}

impl Uniaxial for UniaxialJ2 {
    type State = UniaxialJ2State;
    fn initial_state(&self) -> Self::State {
        UniaxialJ2State {
            plastic_strain: 0.0,
            equivalent_plastic_strain: 0.0,
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state).stress_pa
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state).tangent_pa
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state).state
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.backbone.card()
    }
}

/// Exact free-sided axial reduction of the native Voce saturation law. The
/// three-dimensional owner supplies the nonlinear plastic increment; the
/// axial reduction changes its elastic correction from `3 mu` to `E`.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialVoce {
    backbone: VoceJ2,
}

impl UniaxialVoce {
    /// Admit the same coefficients as the native Voce and axial J2 owners.
    pub fn new(backbone: VoceJ2) -> Result<Self, MaterialError> {
        UniaxialJ2::new(backbone.backbone().clone())?;
        let stiffness = backbone.backbone().elastic.youngs
            + backbone.backbone().kinematic_hardening
            + backbone.isotropic_hardening_derivative_pa(0.0);
        if !stiffness.is_finite() || stiffness <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "axial Voce return stiffness must be finite and positive".into(),
            });
        }
        Ok(Self { backbone })
    }

    /// Original native law, with independent linear and saturating coefficients.
    pub const fn backbone(&self) -> &VoceJ2 {
        &self.backbone
    }

    /// Pure candidate at fixed accepted history. The returned tensor includes
    /// both Poisson contraction and incompressible lateral plastic flow.
    /// This raw kernel follows `Uniaxial` admission; checked structural or
    /// energetic entry points reject nonfinite results and domain exits.
    pub fn response(&self, strain: f64, state: &UniaxialJ2State) -> UniaxialJ2Response {
        let linear = self.backbone.backbone();
        if self.backbone.saturation_pa() == 0.0 {
            return UniaxialJ2 {
                backbone: linear.clone(),
            }
            .response(strain, state);
        }
        let e = linear.elastic.youngs;
        let trial =
            e * (strain - state.plastic_strain) - linear.kinematic_hardening * state.plastic_strain;
        let delta =
            self.backbone
                .plastic_increment(trial.abs(), state.equivalent_plastic_strain, e);
        let mut next = state.clone();
        next.plastic_strain += trial.signum() * delta;
        next.equivalent_plastic_strain += delta;
        let tangent_pa = if delta == 0.0 {
            e
        } else {
            let hardening = linear.kinematic_hardening
                + self
                    .backbone
                    .isotropic_hardening_derivative_pa(next.equivalent_plastic_strain);
            e * (hardening / (e + hardening))
        };
        let elastic_strain = strain - next.plastic_strain;
        let transverse = -linear.elastic.poisson * elastic_strain - 0.5 * next.plastic_strain;
        UniaxialJ2Response {
            strain: [strain, transverse, transverse, 0.0, 0.0, 0.0],
            stress_pa: e * elastic_strain,
            tangent_pa,
            state: next,
        }
    }

    /// Preserve the native energetic interpretation of saturation. Recoverable
    /// hardening and discrete integration loss are accounted separately from heat.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialJ2State,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        let transverse = -self.backbone.backbone().elastic.poisson
            * (old_strain - old_state.plastic_strain)
            - 0.5 * old_state.plastic_strain;
        let old_tensor = [old_strain, transverse, transverse, 0.0, 0.0, 0.0];
        let next = self.response(strain, old_state);
        self.backbone.advance_with_energy(
            &old_tensor,
            &old_state.tensor_history(),
            &next.strain,
            tolerances,
        )
    }
}

impl Uniaxial for UniaxialVoce {
    type State = UniaxialJ2State;
    fn initial_state(&self) -> Self::State {
        UniaxialJ2State {
            plastic_strain: 0.0,
            equivalent_plastic_strain: 0.0,
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state).stress_pa
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state).tangent_pa
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state).state
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.backbone.card()
    }
}

/// Numerical controls for the bounded plane-stress local solve. These are
/// numerical acceptance limits, not a material calibration or accuracy claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReductionConfig {
    /// Maximum local Newton corrections, from 1 through 64.
    pub max_iterations: usize,
    /// Maximum absolute constrained stress component [Pa].
    pub stress_tolerance_pa: f64,
    /// Maximum full tensor Frobenius strain norm, in (0, 0.05]. The caller must
    /// choose a bound no larger than the underlying material's admitted domain.
    pub strain_limit: f64,
}

/// Full accepted strain and original three-dimensional internal variables.
#[derive(Clone, Debug, PartialEq)]
pub struct ReducedState<S> {
    /// Full tensor strain, retaining transverse contraction and shear.
    pub strain: Voigt,
    /// Original constitutive history. No transverse plastic components are lost.
    pub history: S,
}

/// In-plane response in tensor-shear order [xx, yy, xy].
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneStressResponse<S> {
    /// Candidate accepted full strain and original internal variables.
    pub state: ReducedState<S>,
    /// In-plane stress [Pa] in order [xx, yy, xy].
    pub stress_pa: [f64; 3],
    /// Derivative with respect to [epsilon_xx, epsilon_yy, epsilon_xy].
    /// Shear is TENSOR shear, so work is sx*ex + sy*ey + 2*txy*exy.
    pub tangent_pa: [[f64; 3]; 3],
    /// All endpoint stresses, including the actual constrained residuals [Pa].
    pub full_stress_pa: Voigt,
    /// Number of local Newton corrections; material time advances only once.
    pub iterations: usize,
}

/// Explicit refusal; the input strain/history is unchanged on every failure.
#[derive(Clone, Debug, PartialEq)]
pub enum ReductionError {
    /// Invalid controls, nonfinite arithmetic or an exceeded strain domain.
    Invalid(&'static str),
    /// The shared dense solver refused the transverse tangent or workspace.
    Linear(fs_la::LuWorkspaceError),
    /// A finite local iterate did not meet the caller's stress allowance.
    NonConvergence {
        /// Number of corrections attempted.
        iterations: usize,
        /// Largest absolute constrained stress [Pa].
        transverse_residual_pa: f64,
    },
}
impl core::fmt::Display for ReductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material plane-stress reduction: {self:?}")
    }
}
impl std::error::Error for ReductionError {}
impl From<fs_la::LuWorkspaceError> for ReductionError {
    fn from(error: fs_la::LuWorkspaceError) -> Self {
        Self::Linear(error)
    }
}

/// Reusable three-by-three factorization storage. Prepare once for repeated
/// quadrature-point evaluations; no allocation occurs in local Newton loops.
#[derive(Debug)]
pub struct PlaneStressWorkspace {
    transverse: fs_la::LuWorkspace,
}
impl PlaneStressWorkspace {
    /// Allocate the fixed local workspace through the existing dense-solve owner.
    pub fn new() -> Result<Self, ReductionError> {
        Ok(Self {
            transverse: fs_la::LuWorkspace::new(3)?,
        })
    }
}

/// Checked plane-stress reduction of any compatible native small-strain law.
///
/// Solves sigma_zz=sigma_yz=sigma_zx=0 with frozen OLD material history, then
/// returns C_aa-C_ab C_bb^-1 C_ba. Stress, tangent and commit all use the SAME
/// converged strain and incoming history, including for Perzyna kinetics. This
/// wrapper neither alters the constitutive law nor adds elapsed physical time.
/// The generic caller remains responsible for the underlying law/state domain.
/// Reference: Simo and Taylor, plane-stress elastoplastic return mapping,
/// https://escholarship.org/uc/item/9919000n .
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneStressReduction<L> {
    law: L,
    config: ReductionConfig,
}

const IN_PLANE: [usize; 3] = [0, 1, 3];
const TRANSVERSE: [usize; 3] = [2, 4, 5];

impl<L: SmallStrainLaw> PlaneStressReduction<L> {
    /// Select explicit numerical controls without inventing material parameters.
    pub fn new(law: L, config: ReductionConfig) -> Result<Self, ReductionError> {
        if config.max_iterations == 0
            || config.max_iterations > 64
            || !config.stress_tolerance_pa.is_finite()
            || config.stress_tolerance_pa <= 0.0
            || !config.strain_limit.is_finite()
            || config.strain_limit <= 0.0
            || config.strain_limit > 0.05
        {
            return Err(ReductionError::Invalid(
                "invalid iteration, stress or small-strain allowance",
            ));
        }
        Ok(Self { law, config })
    }
    /// Original immutable law, including any source bundle or physical duration.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Virgin material history at zero total strain. A law with eigenstrain may
    /// need an initial constrained evaluation to establish a stress-free shape.
    pub fn initial_state(&self) -> ReducedState<L::State> {
        ReducedState {
            strain: [0.0; 6],
            history: self.law.initial_state(),
        }
    }

    /// One checked local solve. Repeated element/section assembly can reuse
    /// `evaluate_with_workspace` to avoid repeated workspace allocation.
    pub fn evaluate(
        &self,
        strain: [f64; 3],
        old: &ReducedState<L::State>,
    ) -> Result<PlaneStressResponse<L::State>, ReductionError> {
        self.evaluate_with_workspace(strain, old, &mut PlaneStressWorkspace::new()?)
    }

    /// Evaluate a candidate from the same accepted history throughout Newton,
    /// line search and condensation. No input state is mutated on success/error.
    pub fn evaluate_with_workspace(
        &self,
        in_plane: [f64; 3],
        old: &ReducedState<L::State>,
        workspace: &mut PlaneStressWorkspace,
    ) -> Result<PlaneStressResponse<L::State>, ReductionError> {
        if !self.admitted_strain(&old.strain) || in_plane.iter().any(|v| !v.is_finite()) {
            return Err(ReductionError::Invalid(
                "nonfinite or out-of-domain incoming strain",
            ));
        }
        let mut strain = old.strain;
        for i in 0..3 {
            strain[IN_PLANE[i]] = in_plane[i];
        }
        // The old transverse guess may be too large after a large in-plane
        // reversal even though the new constrained solution is in-domain.
        if !self.admitted_strain(&strain) {
            for i in TRANSVERSE {
                strain[i] = 0.0;
            }
        }
        if !self.admitted_strain(&strain) {
            return Err(ReductionError::Invalid(
                "in-plane strain exceeds the declared small-strain domain",
            ));
        }
        for iteration in 0..=self.config.max_iterations {
            let stress = self.law.stress(&strain, &old.history);
            let residual = transverse_residual(&stress)?;
            let tangent = self.law.tangent(&strain, &old.history);
            if tangent.iter().flatten().any(|v| !v.is_finite()) {
                return Err(ReductionError::Invalid(
                    "nonfinite native algorithmic tangent",
                ));
            }
            let block: [f64; 9] =
                core::array::from_fn(|k| tangent[TRANSVERSE[k / 3]][TRANSVERSE[k % 3]]);
            workspace.transverse.factor(&block)?;
            if residual <= self.config.stress_tolerance_pa {
                let mut reduced = [[0.0; 3]; 3];
                for j in 0..3 {
                    let rhs =
                        core::array::from_fn::<_, 3, _>(|i| tangent[TRANSVERSE[i]][IN_PLANE[j]]);
                    let mut eliminated = [0.0; 3];
                    workspace
                        .transverse
                        .solve_factored_into(&rhs, &mut eliminated)?;
                    for i in 0..3 {
                        reduced[i][j] = tangent[IN_PLANE[i]][IN_PLANE[j]];
                        for k in 0..3 {
                            reduced[i][j] -= tangent[IN_PLANE[i]][TRANSVERSE[k]] * eliminated[k];
                        }
                    }
                }
                if reduced.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(ReductionError::Invalid("unrepresentable condensed tangent"));
                }
                return Ok(PlaneStressResponse {
                    state: ReducedState {
                        strain,
                        history: self.law.update_state(&strain, &old.history),
                    },
                    stress_pa: IN_PLANE.map(|i| stress[i]),
                    tangent_pa: reduced,
                    full_stress_pa: stress,
                    iterations: iteration,
                });
            }
            if iteration == self.config.max_iterations {
                return Err(ReductionError::NonConvergence {
                    iterations: iteration,
                    transverse_residual_pa: residual,
                });
            }
            let rhs = TRANSVERSE.map(|i| -stress[i]);
            let mut correction = [0.0; 3];
            workspace
                .transverse
                .solve_factored_into(&rhs, &mut correction)?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..20 {
                let mut candidate = strain;
                for i in 0..3 {
                    candidate[TRANSVERSE[i]] += fraction * correction[i];
                }
                if self.admitted_strain(&candidate) {
                    let trial = self.law.stress(&candidate, &old.history);
                    if let Ok(candidate_residual) = transverse_residual(&trial) {
                        if candidate_residual <= self.config.stress_tolerance_pa
                            || candidate_residual < (1.0 - 1e-4 * fraction) * residual
                        {
                            accepted = Some(candidate);
                            break;
                        }
                    }
                }
                fraction *= 0.5;
            }
            match accepted {
                Some(next) => strain = next,
                None => {
                    return Err(ReductionError::NonConvergence {
                        iterations: iteration + 1,
                        transverse_residual_pa: residual,
                    });
                }
            }
        }
        unreachable!("bounded Newton exits on its final iteration")
    }

    fn admitted_strain(&self, strain: &Voigt) -> bool {
        strain.iter().all(|v| v.is_finite())
            && crate::tensor::contract(strain, strain)
                <= self.config.strain_limit * self.config.strain_limit
    }
}

fn transverse_residual(stress: &Voigt) -> Result<f64, ReductionError> {
    if stress.iter().any(|v| !v.is_finite()) {
        return Err(ReductionError::Invalid("nonfinite native stress"));
    }
    Ok(TRANSVERSE
        .into_iter()
        .map(|i| stress[i].abs())
        .fold(0.0, f64::max))
}
