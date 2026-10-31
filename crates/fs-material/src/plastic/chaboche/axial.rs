//! Exact free-sided axial restriction of nonlinear kinematic hardening.
//!
//! The tensor owner eliminates the independent backstresses and solves the
//! same scalar consistency equation with E replacing 3 mu. Lateral elastic
//! contraction and plastic incompressibility remain in the returned strain.
//! No second return map or independently fitted one-dimensional law is used.

use super::{ChabocheEnergyIncrement, ChabocheError, ChabocheJ2, ChabocheState};
use crate::fiber::Uniaxial;
use crate::plastic::energy::EnergyTolerances;
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;

/// Exact axial-stress restriction; one branch is Armstrong-Frederick.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialChaboche<const N: usize = 1> {
    backbone: ChabocheJ2<N>,
}

/// Independent axial backstresses retained on the uniaxial-flow subspace.
/// Private history prevents importing an arbitrary multiaxial tensor state.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialChabocheState<const N: usize = 1> {
    history: ChabocheState<N>,
}

impl<const N: usize> UniaxialChabocheState<N> {
    /// Signed plastic strain in the fiber direction.
    pub const fn plastic_strain(&self) -> f64 {
        self.history.plastic.plastic_strain[0]
    }
    /// Monotonically accumulated equivalent plastic strain.
    pub const fn equivalent_plastic_strain(&self) -> f64 {
        self.history.plastic.alpha
    }
    /// Original tensor history, including all independent backstresses.
    pub const fn tensor_history(&self) -> &ChabocheState<N> {
        &self.history
    }
    /// Physical axial translations of the yield interval [Pa].
    pub fn backstress_pa(&self) -> [f64; N] {
        core::array::from_fn(|i| 1.5 * self.history.backstress_pa[i][0])
    }
}

/// Pure checked axial candidate and its exact consistent tangent.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialChabocheResponse<const N: usize = 1> {
    /// Full compatible total strain, including free-sided lateral contraction.
    pub strain: Voigt,
    /// Axial physical stress [Pa].
    pub stress_pa: f64,
    /// Derivative of axial stress at frozen incoming history [Pa].
    pub tangent_pa: f64,
    /// Candidate accepted history; the input remains unchanged.
    pub state: UniaxialChabocheState<N>,
}

impl<const N: usize> UniaxialChaboche<N> {
    /// Admit the original law and its representable axial return modulus.
    pub fn new(backbone: ChabocheJ2<N>) -> Result<Self, MaterialError> {
        let linear = backbone.backbone().backbone();
        let mut modulus = linear.elastic.youngs
            + linear.kinematic_hardening
            + backbone.backbone().isotropic_hardening_derivative_pa(0.0);
        for branch in backbone.branches() {
            modulus += branch.modulus_pa();
        }
        if !modulus.is_finite() || modulus <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "axial Chaboche return modulus must be finite and positive".into(),
            });
        }
        Ok(Self { backbone })
    }

    /// Original tensor owner, with independent isotropic and backstress data.
    pub const fn backbone(&self) -> &ChabocheJ2<N> {
        &self.backbone
    }

    fn compatible_strain(&self, strain: f64, plastic: f64) -> Voigt {
        let poisson = self.backbone.backbone().backbone().elastic.poisson;
        let transverse = -poisson * (strain - plastic) - 0.5 * plastic;
        [strain, transverse, transverse, 0.0, 0.0, 0.0]
    }

    /// Solve the native backward-Euler update with exactly zero lateral stress.
    /// Numerical refusal or a full-tensor small-strain domain exit is explicit.
    pub fn response(
        &self,
        strain: f64,
        state: &UniaxialChabocheState<N>,
    ) -> Result<UniaxialChabocheResponse<N>, ChabocheError> {
        let linear = self.backbone.backbone().backbone();
        let e = linear.elastic.youngs;
        let trial = e * (strain - state.plastic_strain())
            - linear.kinematic_hardening * state.plastic_strain();
        // q of this tensor is |trial|; branch shifts are subtracted by the
        // shared owner, including their dependence on the plastic increment.
        let effective = [2.0 * trial / 3.0, -trial / 3.0, -trial / 3.0, 0.0, 0.0, 0.0];
        let returned = self
            .backbone
            .solve_increment(&effective, &state.history, e)?;
        let next = self.backbone.returned_state(&state.history, &returned);
        let tensor = self.compatible_strain(strain, next.plastic.plastic_strain[0]);
        self.backbone.validate_history(&tensor, &next)?;
        let stress_pa = e * (strain - next.plastic.plastic_strain[0]);
        let tangent_pa = if returned.delta_alpha == 0.0 {
            e
        } else {
            e * ((returned.return_modulus_pa - e) / returned.return_modulus_pa)
        };
        if !stress_pa.is_finite() || !tangent_pa.is_finite() {
            return Err(ChabocheError::Invalid(
                "nonrepresentable axial stress or tangent",
            ));
        }
        Ok(UniaxialChabocheResponse {
            strain: tensor,
            stress_pa,
            tangent_pa,
            state: UniaxialChabocheState { history: next },
        })
    }

    /// Audit axial work through the original full-tensor energy owner. Recovery
    /// dissipation is already included in the physical heat-eligible total.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialChabocheState<N>,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<ChabocheEnergyIncrement<N>, ChabocheError> {
        let old_tensor = self.compatible_strain(old_strain, old_state.plastic_strain());
        let next = self.response(strain, old_state)?;
        self.backbone
            .advance_with_energy(&old_tensor, &old_state.history, &next.strain, tolerances)
    }
}

impl<const N: usize> Uniaxial for UniaxialChaboche<N> {
    type State = UniaxialChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        UniaxialChabocheState {
            history: self.backbone.initial_state(),
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |r| r.stress_pa)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |r| r.tangent_pa)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state).map_or_else(
            |_| UniaxialChabocheState {
                history: ChabocheState {
                    plastic: crate::plastic::J2State {
                        plastic_strain: [f64::NAN; 6],
                        alpha: f64::NAN,
                    },
                    backstress_pa: [[f64::NAN; 6]; N],
                },
            },
            |r| r.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        // A scalar axial tangent is symmetric even when the full tensor law is not.
        let mut admission = self.backbone.admissibility();
        admission.tangent_symmetric = true;
        admission
    }
    fn card(&self) -> ModelCard {
        self.backbone.card()
    }
}
