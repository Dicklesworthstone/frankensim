//! Joint Chaboche/Voce hardening and associative Perzyna kinetics.
//!
//! The native consistency F(d) includes the d-dependent nonproportional flow
//! direction and every implicit branch. This owner solves F(d)-f_v(d)=0,
//! f_v=s0[eta d/(dt s0)]^(1/n), then adds f_v/(n d) to the native return
//! modulus. Consequently the native tensor tangent retains direction drift;
//! freezing the backstress direction or sequentially correcting viscosity
//! would produce a different law. The inverse kinetic function is shared with
//! `VoceOverstressJ2`, including its exact n=1 and no-branch reductions.
//!
//! All native storage and numerical-loss terms remain unchanged. The native
//! energy owner adds only f_v d to physical yield/recovery dissipation.
//! Fixed-condition small-strain coefficients and one actual physical interval
//! per law; no sub-yield creep, static recovery or empirical calibration.

use super::{
    ChabocheEnergyIncrement, ChabocheError, ChabocheJ2, ChabocheResponse, ChabocheState,
    ReturnSolution,
};
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::rate::VoceOverstressJ2;
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};
use fs_math::det;

/// Native N-branch Chaboche history with joint, explicitly timed kinetics.
/// State and response use fixed arrays, with no point-response heap allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheOverstressJ2<const N: usize = 1> {
    backbone: ChabocheJ2<N>,
    kinetics: VoceOverstressJ2,
}

impl<const N: usize> ChabocheOverstressJ2<N> {
    /// alpha_dot=max(f,0)/eta with positive viscosity [Pa s] and actual dt [s].
    pub fn new_linear(
        backbone: ChabocheJ2<N>,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let kinetics =
            VoceOverstressJ2::new_linear(backbone.backbone().clone(), viscosity_pa_s, duration_s)?;
        Self::with_kinetics(backbone, kinetics)
    }

    /// alpha_dot=(s0/eta)max(f/s0,0)^n with n>=1 and positive s0 [Pa].
    pub fn new_power_law(
        backbone: ChabocheJ2<N>,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let kinetics = VoceOverstressJ2::new_power_law(
            backbone.backbone().clone(),
            viscosity_pa_s,
            stress_scale_pa,
            exponent,
            duration_s,
        )?;
        Self::with_kinetics(backbone, kinetics)
    }

    fn with_kinetics(
        backbone: ChabocheJ2<N>,
        kinetics: VoceOverstressJ2,
    ) -> Result<Self, MaterialError> {
        let linear = backbone.backbone().backbone();
        let (_, mu) = linear.elastic.lame();
        let mut slope = 3.0 * mu
            + linear.kinematic_hardening
            + backbone.backbone().isotropic_hardening_derivative_pa(0.0);
        for branch in backbone.branches() {
            slope += branch.modulus_pa();
        }
        if !(slope + kinetics.linear_overstress_modulus_pa()).is_finite() {
            return Err(MaterialError::Parameters {
                what: "combined Chaboche and overstress modulus is not representable".into(),
            });
        }
        Ok(Self { backbone, kinetics })
    }

    pub const fn backbone(&self) -> &ChabocheJ2<N> {
        &self.backbone
    }
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.kinetics.viscosity_pa_s()
    }
    pub const fn duration_s(&self) -> f64 {
        self.kinetics.duration_s()
    }
    pub const fn exponent(&self) -> f64 {
        self.kinetics.exponent()
    }
    pub const fn stress_scale_pa(&self) -> Option<f64> {
        self.kinetics.stress_scale_pa()
    }

    /// Native checked stress/history and generally nonsymmetric full tangent.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        self.backbone
            .response_with_rate(strain, state, true, Some(&self.kinetics))
    }

    /// One candidate, with physical yield, viscous and branch-recovery heat.
    /// Recovery is already included in the returned physical-dissipation total.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ChabocheEnergyIncrement<N>, ChabocheError> {
        self.backbone.advance_with_rate_energy(
            old_strain,
            old_state,
            strain,
            tolerances,
            Some(&self.kinetics),
        )
    }
}

impl<const N: usize> ChabocheJ2<N> {
    pub(super) fn solve_rate_increment(
        &self,
        trial: &Voigt,
        state: &ChabocheState<N>,
        elastic_correction_pa: f64,
        kinetics: &VoceOverstressJ2,
    ) -> Result<ReturnSolution, ChabocheError> {
        // Linear viscosity is exactly an additional consistency modulus, not
        // stored hardening. Reuse the original scalar Newton arithmetic.
        if kinetics.exponent() == 1.0 {
            return self.solve_increment(
                trial,
                state,
                elastic_correction_pa + kinetics.linear_overstress_modulus_pa(),
            );
        }
        self.validate_branches(state)?;
        if !elastic_correction_pa.is_finite()
            || elastic_correction_pa <= 0.0
            || trial.iter().any(|x| !x.is_finite())
        {
            return Err(ChabocheError::Invalid(
                "invalid rate-dependent elastic predictor",
            ));
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
        let upper_delta = initial / minimum_modulus;
        if !upper_delta.is_finite() || upper_delta <= 0.0 {
            return Err(ChabocheError::Invalid(
                "nonrepresentable rate-dependent plastic bracket",
            ));
        }
        let evaluate = |log_delta: f64| -> Result<(f64, ReturnSolution, f64, f64), ChabocheError> {
            let delta = det::exp(log_delta);
            let (native, solution) =
                self.consistency(trial, state, elastic_correction_pa, delta)?;
            let viscous = kinetics.overstress_from_log_increment(log_delta);
            let residual = native - viscous;
            // Derivative with respect to log(delta), avoiding the singular
            // power-law slope at zero during bracket evaluation.
            let derivative = solution.return_modulus_pa * delta + viscous / kinetics.exponent();
            if !residual.is_finite() || !derivative.is_finite() || derivative < 0.0 {
                return Err(ChabocheError::Invalid(
                    "nonrepresentable joint Chaboche rate equation",
                ));
            }
            Ok((residual, solution, derivative, viscous))
        };
        let mut lower = det::ln(f64::from_bits(1));
        let mut upper = det::ln(upper_delta);
        if lower > upper || evaluate(lower)?.0 < 0.0 {
            return Err(ChabocheError::Invalid(
                "joint Chaboche rate increment is below representable range",
            ));
        }
        // Reachable branch bounds give the same upper bracket as the native
        // solve; viscosity can only lower the residual further.
        for _ in 0..8 {
            if evaluate(upper)?.0 <= 0.0 {
                break;
            }
            upper += det::ln(2.0);
        }
        if evaluate(upper)?.0 > 0.0 {
            return Err(ChabocheError::LocalSolve {
                iterations: 0,
                residual_pa: evaluate(upper)?.0,
            });
        }
        let mut position = upper;
        let mut residual = initial;
        for iteration in 1..=96 {
            let (value, mut solution, derivative, viscous) = evaluate(position)?;
            residual = value;
            if value.abs() <= tolerance && solution.delta_alpha > 0.0 {
                solution.return_modulus_pa +=
                    (viscous / solution.delta_alpha) / kinetics.exponent();
                if !solution.return_modulus_pa.is_finite() {
                    return Err(ChabocheError::Invalid(
                        "joint Chaboche rate tangent is not representable",
                    ));
                }
                solution.iterations = iteration;
                return Ok(solution);
            }
            if value > 0.0 {
                lower = position;
            } else {
                upper = position;
            }
            let newton = position + value / derivative;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == position || next == lower || next == upper {
                break;
            }
            position = next;
        }
        Err(ChabocheError::LocalSolve {
            iterations: 96,
            residual_pa: residual,
        })
    }
}

impl<const N: usize> SmallStrainLaw for ChabocheOverstressJ2<N> {
    type State = ChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.backbone
            .response_with_rate(strain, state, false, Some(&self.kinetics))
            .map_or([f64::NAN; 6], |r| r.stress_pa)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.backbone
            .response_with_rate(strain, state, false, Some(&self.kinetics))
            .map_or(
                ChabocheState {
                    plastic: crate::plastic::J2State {
                        plastic_strain: [f64::NAN; 6],
                        alpha: f64::NAN,
                    },
                    backstress_pa: [[f64::NAN; 6]; N],
                },
                |r| r.state,
            )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        if self.backbone.absent() {
            return self.kinetics.admissibility();
        }
        MaterialAdmissibility {
            failure_envelope: "small-strain joint Chaboche/Voce J2 overstress, n >= 1; physical recovery",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        if self.backbone.absent() {
            return self.kinetics.card();
        }
        ModelCard::new("material.j2-chaboche-overstress", "0.1.0", Ambition::Solid,
            vec!["small strain, fixed elastic/hardening coefficients".into(),
                "joint Perzyna kinetics and implicit nonlinear kinematic recovery".into(),
                "full generally nonsymmetric tangent including direction drift".into(),
                "energetic Voce/Prager and branch storage; separate physical/numerical dissipation".into(),
                "backward Euler at explicit physical dt".into()], self.backbone.card().validity,
            vec!["no sub-yield creep, static recovery, finite strain, damage or empirical calibration".into()], 0.05)
    }
}
