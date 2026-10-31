//! Independent nonlinear kinematic hardening with retained tensor history.
//!
//! Each Armstrong-Frederick branch follows
//! `d a_i = (2/3) C_i d eps_p - gamma_i a_i d alpha`; summing branches gives
//! the isothermal Chaboche law. See the
//! [NEML Chaboche description](https://neml.readthedocs.io/en/stable/hardening/non/chaboche.html)
//! for the model family (NEML uses the opposite hardening-force sign).
//! The existing Voce owner supplies isotropic hardening and optional linear
//! Prager hardening. Coefficients are fixed throughout each material history.
//!
//! Backward Euler eliminates every branch as
//! `a_i,new = (a_i,old + (2/3) C_i delta_eps_p)/(1+gamma_i delta_alpha)`.
//! The remaining consistency problem is scalar, but the flow direction changes
//! with that scalar during nonproportional loading. The tangent differentiates
//! both dependencies analytically and is generally NONSYMMETRIC. A radial
//! tangent with frozen old backstress would be incorrect.
//!
//! Branch storage is `3 a_i:a_i/(4 C_i)`. Dynamic recovery dissipates
//! `3 gamma_i delta_alpha a_i,new:a_i,new/(2 C_i)` physically; the separate
//! backward-Euler remainder is `3 delta_a_i:delta_a_i/(4 C_i)`. No part of
//! recoverable hardening or numerical integration loss is relabeled as heat.

use super::J2State;
use super::energy::{EnergyError, EnergyTolerances, IncrementEnergy};
use super::voce::VoceJ2;
use crate::tensor::{contract, deviator, trace};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};
use fs_math::det;

/// Exact free-sided axial restriction using the same scalar consistency solve.
pub mod axial;
/// Fixed-condition source binding with independently selected branch coefficients.
pub mod source;

/// Maximum simultaneously retained backstresses; storage is fixed on the stack.
pub const MAX_BACKSTRESS_BRANCHES: usize = 8;

/// Independently selected branch stiffness and dynamic recovery rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackstressBranch {
    modulus_pa: f64,
    recovery_per_strain: f64,
}

impl BackstressBranch {
    /// C >= 0 [Pa], gamma >= 0 [per equivalent plastic strain]. A zero C
    /// branch is absent; gamma=0 gives an additional linear Prager branch.
    pub fn new(modulus_pa: f64, recovery_per_strain: f64) -> Result<Self, MaterialError> {
        if !modulus_pa.is_finite()
            || modulus_pa < 0.0
            || !recovery_per_strain.is_finite()
            || recovery_per_strain < 0.0
        {
            return Err(MaterialError::Parameters {
                what: "backstress C and dynamic recovery gamma must be finite and nonnegative"
                    .into(),
            });
        }
        Ok(Self {
            modulus_pa,
            recovery_per_strain,
        })
    }

    /// Kinematic modulus C [Pa], under the tensor convention `(2/3) C d eps_p`.
    pub const fn modulus_pa(&self) -> f64 {
        self.modulus_pa
    }
    /// Dynamic recovery per unit accumulated equivalent plastic strain.
    pub const fn recovery_per_strain(&self) -> f64 {
        self.recovery_per_strain
    }
}

/// Original plastic strain and equivalent strain plus independent backstresses.
/// Their deviatoric components must survive unloading, reversal and reduction.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheState<const N: usize = 1> {
    /// Existing J2 plastic strain and monotonically accumulated equivalent strain.
    pub plastic: J2State,
    /// N physical stress-space translations in tensor-Voigt convention [Pa].
    pub backstress_pa: [Voigt; N],
}

/// Isothermal rate-independent Chaboche hardening; N=1 is Armstrong-Frederick.
/// N=0 or all-zero C delegates exactly to the existing Voce law.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheJ2<const N: usize = 1> {
    backbone: VoceJ2,
    branches: [BackstressBranch; N],
}

/// One checked candidate at frozen incoming history, including its exact tangent.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheResponse<const N: usize = 1> {
    /// Candidate total tensor strain.
    pub strain: Voigt,
    /// Endpoint physical stress [Pa].
    pub stress_pa: Voigt,
    /// Derivative with respect to tensor strain at fixed incoming history [Pa].
    pub tangent_pa: Tangent6,
    /// Candidate committed history; incoming history is never mutated.
    pub state: ChabocheState<N>,
    /// Number of iterations of this nonlinear-backstress consistency solve.
    pub iterations: usize,
}

/// Checked endpoint work, distinguishing physical recovery from numerical loss.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheEnergyIncrement<const N: usize = 1> {
    /// Candidate total tensor strain.
    pub strain: Voigt,
    /// Endpoint physical stress [Pa].
    pub stress_pa: Voigt,
    /// Candidate retained history.
    pub state: ChabocheState<N>,
    /// Work and storage; physical dissipation ALREADY includes recovery below.
    pub energy: IncrementEnergy,
    /// Dynamic-recovery contribution to physical dissipation [J/m3], not extra heat.
    pub recovery_dissipation_j_m3: f64,
}

/// Explicit parameter/history, local convergence or energy refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ChabocheError {
    /// Invalid state, numerical domain or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The bounded scalar return did not meet its roundoff-scaled residual gate.
    LocalSolve {
        /// Completed consistency iterations.
        iterations: usize,
        /// Signed consistency residual [Pa].
        residual_pa: f64,
    },
    /// Original energetic history or independently selected balance gate failed.
    Energy(EnergyError),
}
impl core::fmt::Display for ChabocheError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Chaboche J2: {self:?}")
    }
}
impl std::error::Error for ChabocheError {}
impl From<EnergyError> for ChabocheError {
    fn from(value: EnergyError) -> Self {
        Self::Energy(value)
    }
}

// Shared with exact axial restrictions: change only 3 mu to E. Direction has
// unit Frobenius norm; return_modulus_pa is -d(consistency)/d(delta_alpha).
#[derive(Clone, Copy, Debug)]
pub(super) struct ReturnSolution {
    pub(super) delta_alpha: f64,
    pub(super) direction: Voigt,
    pub(super) return_modulus_pa: f64,
    pub(super) vector_norm_pa: f64,
    pub(super) direction_drift_pa: Voigt,
    pub(super) iterations: usize,
}

impl<const N: usize> ChabocheJ2<N> {
    /// Admit up to eight explicit branches over the original Voce/Prager law.
    pub fn new(backbone: VoceJ2, branches: [BackstressBranch; N]) -> Result<Self, MaterialError> {
        if N > MAX_BACKSTRESS_BRANCHES {
            return Err(MaterialError::Parameters {
                what: "at most eight Chaboche backstresses are supported".into(),
            });
        }
        let (_, mu) = backbone.backbone().elastic.lame();
        let mut modulus = 3.0 * mu
            + backbone.backbone().kinematic_hardening
            + backbone.isotropic_hardening_derivative_pa(0.0);
        for branch in &branches {
            modulus += branch.modulus_pa;
        }
        if !modulus.is_finite() || modulus <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "combined Chaboche return modulus is not representable".into(),
            });
        }
        Ok(Self { backbone, branches })
    }

    /// Original Voce isotropic and optional linear Prager owner.
    pub const fn backbone(&self) -> &VoceJ2 {
        &self.backbone
    }
    /// Independent branch parameters; their order is also the state order.
    pub const fn branches(&self) -> &[BackstressBranch; N] {
        &self.branches
    }

    fn absent(&self) -> bool {
        self.branches.iter().all(|b| b.modulus_pa == 0.0)
    }

    /// Admit finite native history and physically reachable backstress bounds.
    /// No instantaneous-yield restriction is imposed on the incoming state.
    pub fn validate_history(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
    ) -> Result<(), ChabocheError> {
        self.backbone
            .backbone()
            .validate_overstress_history(strain, &state.plastic)?;
        self.validate_branches(state)
    }

    fn validate_branches(&self, state: &ChabocheState<N>) -> Result<(), ChabocheError> {
        if !state.plastic.alpha.is_finite() || state.plastic.alpha < 0.0 {
            return Err(ChabocheError::Invalid("invalid accumulated plastic strain"));
        }
        for (branch, a) in self.branches.iter().zip(&state.backstress_pa) {
            if a.iter().any(|x| !x.is_finite()) {
                return Err(ChabocheError::Invalid("nonfinite independent backstress"));
            }
            let norm = det::sqrt(contract(a, a));
            let roundoff = 128.0 * f64::EPSILON;
            let q = det::sqrt(1.5) * norm;
            if !q.is_finite()
                || trace(a).abs() > roundoff * norm.max(f64::MIN_POSITIVE)
                || (branch.modulus_pa == 0.0 && a.iter().any(|x| *x != 0.0))
                || q > branch.modulus_pa * state.plastic.alpha * (1.0 + roundoff)
                || q * branch.recovery_per_strain > branch.modulus_pa * (1.0 + roundoff)
            {
                return Err(ChabocheError::Invalid(
                    "backstress violates its deviatoric, accumulated-strain or saturation bound",
                ));
            }
        }
        Ok(())
    }

    /// Actual effective-stress yield residual including every retained branch.
    pub fn yield_function_with_state(&self, stress: &Voigt, state: &ChabocheState<N>) -> f64 {
        let mut effective = *stress;
        for a in &state.backstress_pa {
            for i in 0..6 {
                effective[i] -= a[i];
            }
        }
        self.backbone
            .yield_function_with_state(&effective, &state.plastic)
    }

    // v(delta) = s_trial - a_Prager,old - sum a_i,old/(1+gamma_i delta).
    // Consistency: sqrt(3/2)|v| - (M+Hk)delta - sum C_i delta/(1+gamma_i delta)
    //              - sigma_y - H alpha_new - R_Voce(alpha_new) = 0.
    fn consistency(
        &self,
        trial: &Voigt,
        state: &ChabocheState<N>,
        modulus: f64,
        delta: f64,
    ) -> Result<(f64, ReturnSolution), ChabocheError> {
        let mut vector = *trial;
        let mut drift = [0.0; 6];
        let mut correction = (modulus + self.backbone.backbone().kinematic_hardening) * delta;
        let mut denominator = modulus
            + self.backbone.backbone().kinematic_hardening
            + self
                .backbone
                .isotropic_hardening_derivative_pa(state.plastic.alpha + delta);
        for (branch, a) in self.branches.iter().zip(&state.backstress_pa) {
            let inverse = 1.0 / (1.0 + branch.recovery_per_strain * delta);
            correction += branch.modulus_pa * (delta * inverse);
            for i in 0..6 {
                vector[i] -= a[i] * inverse;
                drift[i] += (branch.recovery_per_strain * a[i]) * inverse * inverse;
            }
        }
        let norm = det::sqrt(contract(&vector, &vector));
        let direction = if norm == 0.0 {
            [0.0; 6]
        } else {
            core::array::from_fn(|i| vector[i] / norm)
        };
        for (branch, a) in self.branches.iter().zip(&state.backstress_pa) {
            let inverse = 1.0 / (1.0 + branch.recovery_per_strain * delta);
            // Sum each branch's net derivative before adding the elastic
            // modulus, preserving that modulus near backstress saturation.
            let recovery = det::sqrt(1.5) * branch.recovery_per_strain * contract(&direction, a);
            denominator += (branch.modulus_pa - recovery) * inverse * inverse;
        }
        let radius = self.backbone.backbone().yield_stress
            + self
                .backbone
                .isotropic_hardening_pa(state.plastic.alpha + delta);
        let residual = det::sqrt(1.5) * norm - correction - radius;
        if !residual.is_finite()
            || !denominator.is_finite()
            || denominator <= 0.0
            || drift.iter().any(|x| !x.is_finite())
        {
            return Err(ChabocheError::Invalid(
                "nonrepresentable local backstress consistency equation",
            ));
        }
        Ok((
            residual,
            ReturnSolution {
                delta_alpha: delta,
                direction,
                return_modulus_pa: denominator,
                vector_norm_pa: norm,
                direction_drift_pa: drift,
                iterations: 0,
            },
        ))
    }

    pub(super) fn solve_increment(
        &self,
        trial: &Voigt,
        state: &ChabocheState<N>,
        elastic_correction_pa: f64,
    ) -> Result<ReturnSolution, ChabocheError> {
        self.validate_branches(state)?;
        if !elastic_correction_pa.is_finite()
            || elastic_correction_pa <= 0.0
            || trial.iter().any(|x| !x.is_finite())
        {
            return Err(ChabocheError::Invalid("invalid local elastic predictor"));
        }
        let (initial, elastic) = self.consistency(trial, state, elastic_correction_pa, 0.0)?;
        if initial <= 0.0 {
            return Ok(elastic);
        }
        let scale = (det::sqrt(1.5) * elastic.vector_norm_pa)
            .max(self.backbone.backbone().yield_stress)
            .max(f64::MIN_POSITIVE);
        let tolerance = 32.0 * f64::EPSILON * scale;
        let minimum_modulus = elastic_correction_pa + self.backbone.backbone().combined_hardening();
        let mut lower = 0.0;
        let mut upper = initial / minimum_modulus;
        if !upper.is_finite() || upper <= 0.0 {
            return Err(ChabocheError::Invalid(
                "nonrepresentable local plastic bracket",
            ));
        }
        // Reachable branch bounds imply F' <= -(M+Hk+H), hence this is
        // an upper bracket. Expand only for finite-precision endpoint noise.
        for _ in 0..8 {
            if self
                .consistency(trial, state, elastic_correction_pa, upper)?
                .0
                <= 0.0
            {
                break;
            }
            upper *= 2.0;
        }
        let (upper_residual, _) = self.consistency(trial, state, elastic_correction_pa, upper)?;
        if upper_residual > 0.0 {
            return Err(ChabocheError::LocalSolve {
                iterations: 0,
                residual_pa: upper_residual,
            });
        }
        let mut delta = (initial / elastic.return_modulus_pa).min(upper);
        let mut residual = initial;
        for iteration in 1..=96 {
            let (value, mut returned) =
                self.consistency(trial, state, elastic_correction_pa, delta)?;
            residual = value;
            returned.iterations = iteration;
            if value.abs() <= tolerance && delta > 0.0 {
                return Ok(returned);
            }
            if value > 0.0 {
                lower = delta;
            } else {
                upper = delta;
            }
            let newton = delta + value / returned.return_modulus_pa;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == delta || next == lower || next == upper {
                break;
            }
            delta = next;
        }
        Err(ChabocheError::LocalSolve {
            iterations: 96,
            residual_pa: residual,
        })
    }

    pub(super) fn returned_state(
        &self,
        state: &ChabocheState<N>,
        solution: &ReturnSolution,
    ) -> ChabocheState<N> {
        if solution.delta_alpha == 0.0 {
            return state.clone();
        }
        let k = det::sqrt(1.5) * solution.delta_alpha;
        let increment: Voigt = core::array::from_fn(|i| k * solution.direction[i]);
        let mut next = state.clone();
        for i in 0..6 {
            next.plastic.plastic_strain[i] += increment[i];
        }
        next.plastic.alpha += solution.delta_alpha;
        for (branch, a) in self.branches.iter().zip(&mut next.backstress_pa) {
            let divisor = 1.0 + branch.recovery_per_strain * solution.delta_alpha;
            for i in 0..6 {
                a[i] = (a[i] + (2.0 / 3.0) * branch.modulus_pa * increment[i]) / divisor;
            }
        }
        next
    }

    fn consistent_tangent(&self, solution: &ReturnSolution) -> Tangent6 {
        let elastic = &self.backbone.backbone().elastic;
        let mut tangent = elastic.stiffness();
        if solution.delta_alpha == 0.0 {
            return tangent;
        }
        let (_, mu) = elastic.lame();
        let n = solution.direction;
        for j in 0..6 {
            let mut unit = [0.0; 6];
            unit[j] = 1.0;
            let projected = deviator(&unit);
            let ds: Voigt = core::array::from_fn(|i| 2.0 * mu * projected[i]);
            let delta_derivative = det::sqrt(1.5) * contract(&n, &ds) / solution.return_modulus_pa;
            let dv: Voigt =
                core::array::from_fn(|i| ds[i] + solution.direction_drift_pa[i] * delta_derivative);
            let projection = contract(&n, &dv);
            for i in 0..6 {
                let direction_derivative = (dv[i] - n[i] * projection) / solution.vector_norm_pa;
                let plastic_derivative = det::sqrt(1.5)
                    * (n[i] * delta_derivative + solution.delta_alpha * direction_derivative);
                tangent[i][j] -= 2.0 * mu * plastic_derivative;
            }
        }
        tangent
    }

    /// Checked candidate response with frozen incoming history and no allocation.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        self.response_impl(strain, state, true)
    }

    // Stress/history-only trait calls do not need to assemble 36 tangent
    // entries. Public response always returns the full consistent tangent.
    fn response_impl(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
        with_tangent: bool,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        self.backbone.backbone().admit_energy_strain(strain)?;
        self.backbone
            .backbone()
            .validate_plastic_history(&state.plastic)?;
        self.validate_branches(state)?;
        if self.absent() {
            let next = ChabocheState {
                plastic: self.backbone.update_state(strain, &state.plastic),
                backstress_pa: state.backstress_pa,
            };
            let stress_pa = self.backbone.stress(strain, &state.plastic);
            let tangent_pa = if with_tangent {
                self.backbone.tangent(strain, &state.plastic)
            } else {
                [[0.0; 6]; 6]
            };
            self.validate_history(strain, &next)?;
            if stress_pa
                .iter()
                .chain(tangent_pa.iter().flatten())
                .any(|x| !x.is_finite())
            {
                return Err(ChabocheError::Invalid("nonfinite original Voce response"));
            }
            return Ok(ChabocheResponse {
                strain: *strain,
                stress_pa,
                tangent_pa,
                state: next,
                iterations: 0,
            });
        }
        let linear = self.backbone.backbone();
        let (_, mu) = linear.elastic.lame();
        let elastic_strain = subtract(strain, &state.plastic.plastic_strain);
        let trial_dev = deviator(&elastic_strain);
        let prager = linear.backstress_pa(&state.plastic);
        let trial = core::array::from_fn(|i| 2.0 * mu * trial_dev[i] - prager[i]);
        let solution = self.solve_increment(&trial, state, 3.0 * mu)?;
        let next = self.returned_state(state, &solution);
        self.validate_history(strain, &next)?;
        let stress_pa = linear
            .elastic
            .stress(&subtract(strain, &next.plastic.plastic_strain), &());
        let tangent_pa = if with_tangent {
            self.consistent_tangent(&solution)
        } else {
            [[0.0; 6]; 6]
        };
        if stress_pa
            .iter()
            .chain(tangent_pa.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(ChabocheError::Invalid(
                "nonfinite returned stress or tangent",
            ));
        }
        Ok(ChabocheResponse {
            strain: *strain,
            stress_pa,
            tangent_pa,
            state: next,
            iterations: solution.iterations,
        })
    }

    /// Audit native elastic/Voce/Prager storage and independent recovery energy.
    /// Coefficients and their energy interpretation must be unchanged over this step.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ChabocheEnergyIncrement<N>, ChabocheError> {
        self.validate_history(old_strain, old_state)?;
        if self.absent() {
            let next = self.backbone.advance_with_energy(
                old_strain,
                &old_state.plastic,
                strain,
                tolerances,
            )?;
            return Ok(ChabocheEnergyIncrement {
                strain: next.strain,
                stress_pa: next.stress_pa,
                state: ChabocheState {
                    plastic: next.state,
                    backstress_pa: old_state.backstress_pa,
                },
                energy: next.energy,
                recovery_dissipation_j_m3: 0.0,
            });
        }
        if !tolerances.energy_j_m3.is_finite()
            || tolerances.energy_j_m3 <= 0.0
            || !tolerances.yield_pa.is_finite()
            || tolerances.yield_pa <= 0.0
        {
            return Err(ChabocheError::Invalid(
                "energy and yield allowances must be finite and positive",
            ));
        }
        let next = self.response_impl(strain, old_state, false)?;
        let delta = next.state.plastic.alpha - old_state.plastic.alpha;
        if delta < 0.0 || delta == 0.0 && next.state != *old_state {
            return Err(ChabocheError::Invalid(
                "nonrepresentable accumulated history increment",
            ));
        }
        let residual = self.yield_function_with_state(&next.stress_pa, &next.state);
        if !residual.is_finite()
            || residual > tolerances.yield_pa
            || delta > 0.0 && residual.abs() > tolerances.yield_pa
        {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            }
            .into());
        }
        let linear = self.backbone.backbone();
        let old_e = subtract(old_strain, &old_state.plastic.plastic_strain);
        let new_e = subtract(strain, &next.state.plastic.plastic_strain);
        let delta_p = subtract(
            &next.state.plastic.plastic_strain,
            &old_state.plastic.plastic_strain,
        );
        let mut initial = linear.energy_of(
            &old_e,
            old_state.plastic.alpha,
            &old_state.plastic.plastic_strain,
        )?;
        let mut final_state = linear.energy_of(
            &new_e,
            next.state.plastic.alpha,
            &next.state.plastic.plastic_strain,
        )?;
        let numerical = linear.energy_of(&subtract(&new_e, &old_e), delta, &delta_p)?;
        initial.hardening_j_m3 += self.backbone.saturation_storage(old_state.plastic.alpha);
        final_state.hardening_j_m3 += self.backbone.saturation_storage(next.state.plastic.alpha);
        let mut numerical_dissipation = numerical.total_j_m3()
            + self
                .backbone
                .saturation_numerical_loss(old_state.plastic.alpha, delta);
        let mut recovery_dissipation = 0.0;
        for ((branch, a0), a1) in self
            .branches
            .iter()
            .zip(&old_state.backstress_pa)
            .zip(&next.state.backstress_pa)
        {
            if branch.modulus_pa == 0.0 {
                continue;
            }
            let initial_branch = 0.75 * (contract(a0, a0) / branch.modulus_pa);
            let final_branch = 0.75 * (contract(a1, a1) / branch.modulus_pa);
            let da = subtract(a1, a0);
            initial.kinematic_hardening_j_m3 += initial_branch;
            final_state.kinematic_hardening_j_m3 += final_branch;
            numerical_dissipation += 0.75 * (contract(&da, &da) / branch.modulus_pa);
            recovery_dissipation += 2.0 * branch.recovery_per_strain * delta * final_branch;
        }
        let mut energy = IncrementEnergy {
            initial,
            final_state,
            mechanical_work_j_m3: contract(&next.stress_pa, &subtract(strain, old_strain)),
            physical_dissipation_j_m3: linear.yield_stress * delta + recovery_dissipation,
            numerical_dissipation_j_m3: numerical_dissipation,
            residual_j_m3: 0.0,
        };
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - energy.stored_change_j_m3()
            - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3;
        if [
            initial.total_j_m3(),
            final_state.total_j_m3(),
            recovery_dissipation,
            energy.physical_dissipation_j_m3,
            energy.numerical_dissipation_j_m3,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || !energy.residual_j_m3.is_finite()
        {
            return Err(ChabocheError::Invalid(
                "nonrepresentable stored or dissipated energy",
            ));
        }
        if energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            }
            .into());
        }
        Ok(ChabocheEnergyIncrement {
            strain: *strain,
            stress_pa: next.stress_pa,
            state: next.state,
            energy,
            recovery_dissipation_j_m3: recovery_dissipation,
        })
    }
}

impl<const N: usize> SmallStrainLaw for ChabocheJ2<N> {
    type State = ChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        ChabocheState {
            plastic: self.backbone.initial_state(),
            backstress_pa: [[0.0; 6]; N],
        }
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.response_impl(strain, state, false)
            .map_or([f64::NAN; 6], |r| r.stress_pa)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.response_impl(strain, state, false).map_or(
            ChabocheState {
                plastic: J2State {
                    plastic_strain: [f64::NAN; 6],
                    alpha: f64::NAN,
                },
                backstress_pa: [[f64::NAN; 6]; N],
            },
            |r| r.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        if self.absent() {
            return self.backbone.admissibility();
        }
        MaterialAdmissibility {
            tangent_symmetric: self
                .branches
                .iter()
                .all(|b| b.modulus_pa == 0.0 || b.recovery_per_strain == 0.0),
            failure_envelope: "J2 with Voce radius and independent Armstrong-Frederick/Chaboche backstresses",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        if self.absent() {
            return self.backbone.card();
        }
        ModelCard::new(
            "material.j2-chaboche",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain, fixed elastic and hardening coefficients".into(),
                "associative J2 plastic flow with nonassociative dynamic backstress recovery"
                    .into(),
                "backward-Euler branch evolution and analytic generally nonsymmetric tangent"
                    .into(),
                "energetic Voce/Prager and branch storage; physical dynamic-recovery dissipation"
                    .into(),
            ],
            self.backbone.card().validity,
            vec![
                "no static recovery, creep, temperature evolution or finite-strain plasticity"
                    .into(),
                "no experimental ratcheting, damage or fatigue-life calibration".into(),
            ],
            0.05,
        )
    }
}

fn subtract(a: &Voigt, b: &Voigt) -> Voigt {
    core::array::from_fn(|i| a[i] - b[i])
}
