//! Explicit-time linear static recovery of independent Chaboche backstresses.
//!
//! In addition to the original dynamic recovery, each branch obeys
//! `dX_i/dt|static = -r_i X_i`, with explicitly selected `r_i >= 0` in s^-1.
//! This is the linear-in-backstress case of the time-recovery family described
//! in <https://neml.readthedocs.io/en/stable/hardening/non/chaboche.html>;
//! r here is the COMPLETE coefficient of -X, not NEML's differently normalized A.
//!
//! Backward Euler gives
//! `X1 = (X0 + (2/3) C deps_p)/(1 + gamma dalpha + r dt)`.
//! With a=1+r dt, the ORIGINAL implicit Chaboche solver evaluates the exactly
//! equivalent coefficients C/a, gamma/a and incoming X0/a. This is simultaneous
//! recovery and plastic flow, not a recovery step followed by a separate return.
//! It preserves the original nonproportional flow direction and exact tangent.
//!
//! These effective coefficients are only an algebraic solve preparation.
//! Recoverable storage remains 3 X:X/(4 C) with ORIGINAL C. The original energy
//! owner accounts physical recovery as 3 (gamma dalpha+r dt) X1:X1/(2 C),
//! separately from the backward-Euler remainder 3 (X1-X0):(X1-X0)/(4 C).
//! No heap allocation or material query occurs in a response. Coefficients are
//! fixed over this interval. Nonlinear static recovery, independent isotropic
//! annealing, recrystallization and finite deformation are outside this law.

use super::axial::{UniaxialChaboche, UniaxialChabocheResponse, UniaxialChabocheState};
use super::rate::ChabocheOverstressJ2;
use super::{
    BackstressBranch, ChabocheEnergyIncrement, ChabocheError, ChabocheJ2, ChabocheResponse,
    ChabocheState,
};
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::rate::VoceOverstressJ2;
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};

/// SI inverse-time dimension of the COMPLETE coefficient in dX/dt = -r X.
pub const STATIC_RECOVERY_RATE_DIMS: fs_qty::Dims = fs_qty::Dims([0, 0, -1, 0, 0, 0]);

/// One immutable physical interval, retaining every original branch and its
/// independent time recovery rate. Repeated proposals borrow the SAME history.
#[derive(Clone, Debug, PartialEq)]
pub struct ChabocheStaticRecoveryJ2<const N: usize = 1> {
    backbone: ChabocheJ2<N>,
    axial: UniaxialChaboche<N>,
    kinetics: Option<VoceOverstressJ2>,
    rates_per_s: [f64; N],
    recovery_increments: [f64; N],
    attenuation: [f64; N],
    duration_s: f64,
}

impl<const N: usize> ChabocheStaticRecoveryJ2<N> {
    /// Select independent linear time recovery and rate-independent yielding.
    /// Explicit zeros retain inactive recovery; missing data imply no rate.
    pub fn new(
        backbone: ChabocheJ2<N>,
        rates_per_s: [f64; N],
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        Self::prepare(backbone, None, rates_per_s, duration_s)
    }

    /// Retain the selected native linear/power-law Perzyna kinetics and its
    /// ACTUAL duration while adding independently supplied time-recovery rates.
    pub fn with_overstress(
        rate: ChabocheOverstressJ2<N>,
        rates_per_s: [f64; N],
    ) -> Result<Self, MaterialError> {
        Self::prepare(
            rate.backbone().clone(),
            Some(rate.kinetics().clone()),
            rates_per_s,
            rate.duration_s(),
        )
    }

    fn prepare(
        backbone: ChabocheJ2<N>,
        kinetics: Option<VoceOverstressJ2>,
        rates_per_s: [f64; N],
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || rates_per_s.iter().any(|r| !r.is_finite() || *r < 0.0)
        {
            return Err(MaterialError::Parameters {
                what: "static recovery needs finite nonnegative rates [s^-1] and a positive physical interval [s]".into(),
            });
        }
        let mut recovery_increments = [0.0; N];
        let mut attenuation = [1.0; N];
        let mut branches = *backbone.branches();
        for i in 0..N {
            let increment = rates_per_s[i] * duration_s;
            let divisor = 1.0 + increment;
            if !increment.is_finite()
                || !divisor.is_finite()
                || rates_per_s[i] > 0.0 && increment == 0.0
            {
                return Err(MaterialError::Parameters {
                    what: "static-recovery rate times physical duration must be representable"
                        .into(),
                });
            }
            recovery_increments[i] = increment;
            attenuation[i] = 1.0 / divisor;
            let original = backbone.branches()[i];
            let modulus = original.modulus_pa() / divisor;
            let gamma = original.recovery_per_strain() / divisor;
            if original.modulus_pa() > 0.0 && modulus == 0.0
                || original.recovery_per_strain() > 0.0 && gamma == 0.0
            {
                return Err(MaterialError::Parameters {
                    what: "static recovery would erase an unresolved active branch coefficient"
                        .into(),
                });
            }
            branches[i] = BackstressBranch::new(modulus, gamma)?;
        }
        let effective = ChabocheJ2::new(backbone.backbone().clone(), branches)?;
        let axial = UniaxialChaboche::new(effective)?;
        Ok(Self {
            backbone,
            axial,
            kinetics,
            rates_per_s,
            recovery_increments,
            attenuation,
            duration_s,
        })
    }

    /// Original constitutive moduli and free energy, never effective storage.
    pub const fn backbone(&self) -> &ChabocheJ2<N> {
        &self.backbone
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn rates_per_s(&self) -> &[f64; N] {
        &self.rates_per_s
    }

    fn prepared_history(&self, state: &ChabocheState<N>) -> ChabocheState<N> {
        let mut prepared = state.clone();
        for i in 0..N {
            if self.recovery_increments[i] != 0.0 {
                for component in &mut prepared.backstress_pa[i] {
                    *component *= self.attenuation[i];
                }
            }
        }
        prepared
    }

    fn response_impl(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
        with_tangent: bool,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        self.backbone.validate_history(strain, state)?;
        let prepared = self.prepared_history(state);
        let candidate = self.axial.backbone().response_with_rate(
            strain,
            &prepared,
            with_tangent,
            self.kinetics.as_ref(),
        )?;
        self.backbone.validate_history(strain, &candidate.state)?;
        Ok(candidate)
    }

    /// Exact simultaneous return, including dwell-driven relaxation/re-yielding
    /// and the original full generally nonsymmetric consistent tangent.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        self.response_impl(strain, state, true)
    }

    /// Exact free-sided axial restriction of the SAME prepared scalar return.
    /// Native branch order, lateral contraction and small-strain bounds persist.
    pub fn uniaxial_response(
        &self,
        strain: f64,
        state: &UniaxialChabocheState<N>,
    ) -> Result<UniaxialChabocheResponse<N>, ChabocheError> {
        self.backbone
            .validate_history(&[0.0; 6], state.tensor_history())?;
        let prepared = UniaxialChabocheState::from_transported_history(
            self.prepared_history(state.tensor_history()),
        );
        let candidate =
            self.axial
                .response_with_kinetics(strain, &prepared, self.kinetics.as_ref())?;
        self.backbone
            .validate_history(&candidate.strain, candidate.state.tensor_history())?;
        Ok(candidate)
    }

    /// Audit original branch free energy, physical recovery/yield/viscous heat
    /// and numerical work separately. A candidate does not advance its input.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ChabocheEnergyIncrement<N>, ChabocheError> {
        let candidate = self.response_impl(strain, old_state, false)?;
        self.backbone.audit_increment(
            old_strain,
            old_state,
            strain,
            candidate,
            tolerances,
            self.kinetics.as_ref(),
            &self.recovery_increments,
        )
    }
}

impl<const N: usize> SmallStrainLaw for ChabocheStaticRecoveryJ2<N> {
    type State = ChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.response_impl(strain, state, false)
            .map_or([f64::NAN; 6], |candidate| candidate.stress_pa)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |candidate| candidate.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.response_impl(strain, state, false).map_or_else(
            |_| ChabocheState {
                plastic: crate::plastic::J2State {
                    plastic_strain: [f64::NAN; 6],
                    alpha: f64::NAN,
                },
                backstress_pa: [[f64::NAN; 6]; N],
            },
            |candidate| candidate.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "small-strain joint Chaboche flow and explicit linear static recovery",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-chaboche-static-recovery",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain; original Voce, Prager and independent branch storage".into(),
                "linear static backstress recovery at explicit physical dt".into(),
                "simultaneous native implicit flow, direction drift and exact tangent".into(),
                "physical recovery and optional Perzyna heat exclude numerical loss".into(),
            ],
            self.backbone.card().validity,
            vec![
                "no nonlinear static recovery, isotropic annealing or recrystallization".into(),
                "no finite deformation or empirical cyclic/dwell calibration".into(),
            ],
            0.05,
        )
    }
}
