//! Exact free-sided axial restriction of tabulated Perzyna plasticity.
//!
//! The original scalar consistency owner uses `E` instead of `3 mu`; the
//! measured resistance, rate kinetics and Prager translation remain unchanged.
//! Poisson contraction and incompressible lateral plastic flow reconstruct the
//! full strain. The native tensor owner supplies stored energy, the reduced
//! potential and the checked work/heat account, without another energy model.
//!
//! One law represents one actual physical interval. Rebinding its duration
//! preserves the original calibration and optional source evidence; accepted
//! axial history is supplied separately and is never reset or advanced by it.

use super::TabulatedOverstressJ2;
use super::source::ResolvedTabulatedOverstressJ2;
use crate::fiber::Uniaxial;
use crate::plastic::energy::{EnergyError, EnergyIncrement, EnergyTolerances, StoredEnergy};
use crate::plastic::reduction::{UniaxialJ2, UniaxialJ2Response, UniaxialJ2State};
use crate::plastic::tabulated::{TabulatedError, TabulatedJ2};
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;

#[derive(Clone, Debug, PartialEq)]
enum Backbone {
    Authored(TabulatedOverstressJ2),
    Sourced(ResolvedTabulatedOverstressJ2),
}

/// Native tabulated rate law on the original axial J2 history subspace.
/// Source-bound laws retain the original curve and scalar-property receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialTabulatedOverstress {
    backbone: Backbone,
}

impl UniaxialTabulatedOverstress {
    /// Restrict an admitted tensor law without fitting axial coefficients.
    pub fn new(backbone: TabulatedOverstressJ2) -> Result<Self, MaterialError> {
        UniaxialJ2::new(backbone.backbone().backbone().clone())?;
        Ok(Self {
            backbone: Backbone::Authored(backbone),
        })
    }

    /// Select linear viscosity at the actual physical interval in seconds.
    pub fn new_linear(
        backbone: TabulatedJ2,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        Self::new(TabulatedOverstressJ2::new_linear(
            backbone,
            viscosity_pa_s,
            duration_s,
        )?)
    }

    /// Select power-law viscosity with the original positive stress scale and
    /// exponent `n >= 1`, at the actual physical interval in seconds.
    pub fn new_power_law(
        backbone: TabulatedJ2,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        Self::new(TabulatedOverstressJ2::new_power_law(
            backbone,
            viscosity_pa_s,
            stress_scale_pa,
            exponent,
            duration_s,
        )?)
    }

    /// Preserve the original source selection, context, curve and kinetics.
    /// Its immutable source is already shared; no duplicate receipt is created.
    pub fn from_source(source: ResolvedTabulatedOverstressJ2) -> Result<Self, MaterialError> {
        UniaxialJ2::new(source.law().backbone().backbone().clone())?;
        Ok(Self {
            backbone: Backbone::Sourced(source),
        })
    }

    /// Original tensor law for this physical interval.
    pub fn backbone(&self) -> &TabulatedOverstressJ2 {
        match &self.backbone {
            Backbone::Authored(law) => law,
            Backbone::Sourced(source) => source.law(),
        }
    }

    /// Complete original source evidence, when the law was source-bound.
    pub fn source(&self) -> Option<&ResolvedTabulatedOverstressJ2> {
        match &self.backbone {
            Backbone::Authored(_) => None,
            Backbone::Sourced(source) => Some(source),
        }
    }

    /// Original scalar-property bundle; the measured curve stays in `source`.
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source().map(ResolvedTabulatedOverstressJ2::material)
    }

    /// Actual backward-Euler interval, never an equilibrium-iteration count.
    pub fn duration_s(&self) -> f64 {
        self.backbone().duration_s()
    }

    /// Prepare another physical interval using the same immutable calibration.
    /// The original axial history must be passed to the new law's response;
    /// preparing a duration performs no constitutive advance or source query.
    pub fn for_duration(&self, duration_s: f64) -> Result<Self, MaterialError> {
        match &self.backbone {
            Backbone::Sourced(source) => {
                Self::from_source(source.for_duration(duration_s).map_err(|error| {
                    MaterialError::Parameters {
                        what: error.to_string(),
                    }
                })?)
            }
            Backbone::Authored(law) => match law.stress_scale_pa() {
                Some(stress_scale_pa) => Self::new_power_law(
                    law.backbone().clone(),
                    law.viscosity_pa_s(),
                    stress_scale_pa,
                    law.exponent(),
                    duration_s,
                ),
                None => Self::new_linear(law.backbone().clone(), law.viscosity_pa_s(), duration_s),
            },
        }
    }

    /// Reconstruct an accepted free-sided strain without advancing history.
    /// Positive viscous overstress is admissible between physical steps.
    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<Voigt, EnergyError> {
        let poisson = self.backbone().backbone().backbone().elastic.poisson;
        let plastic = state.plastic_strain;
        let transverse = -poisson * (axial_strain - plastic) - 0.5 * plastic;
        let strain = [axial_strain, transverse, transverse, 0.0, 0.0, 0.0];
        self.backbone()
            .validate_history(&strain, &state.tensor_history())?;
        Ok(strain)
    }

    /// Check the original native history, full tensor strain and table support
    /// without imposing an instantaneous rate-independent yield surface.
    pub fn validate_history(
        &self,
        strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<(), EnergyError> {
        self.tensor_strain(strain, state).map(|_| ())
    }

    /// Pure candidate with zero lateral stress, using the native bounded rate
    /// consistency solve. The axial tangent is `E H*/(E+H*)`, where
    /// `H* = H_k + Y'(alpha_new) + f_v/(n delta_alpha)` on plastic loading.
    /// At a measured hardening kink the native table chooses its right slope.
    pub fn response(
        &self,
        strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<UniaxialJ2Response, TabulatedError> {
        let backbone = self.backbone();
        let linear = backbone.backbone().backbone();
        let old_tensor = state.tensor_history();
        linear
            .validate_plastic_history(&old_tensor)
            .map_err(|_| TabulatedError::State("invalid axial plastic history"))?;
        if !strain.is_finite() {
            return Err(TabulatedError::State("axial strain must be finite"));
        }
        let e = linear.elastic.youngs;
        let trial =
            e * (strain - state.plastic_strain) - linear.kinematic_hardening * state.plastic_strain;
        let delta = backbone.plastic_increment(trial.abs(), state.equivalent_plastic_strain, e)?;
        let mut next = state.clone();
        let tangent_pa = if delta == 0.0 {
            e
        } else {
            next.plastic_strain += trial.signum() * delta;
            next.equivalent_plastic_strain += delta;
            if next.plastic_strain == state.plastic_strain
                || next.equivalent_plastic_strain <= state.equivalent_plastic_strain
            {
                return Err(TabulatedError::Numerical(
                    "axial rate increment is below retained-history resolution",
                ));
            }
            let hardening = linear.kinematic_hardening
                + backbone
                    .backbone()
                    .hardening()
                    .tangent_pa(next.equivalent_plastic_strain)?;
            let viscous = (backbone.endpoint_overstress_pa(delta)? / delta) / backbone.exponent();
            let stiffness = hardening + viscous;
            let denominator = e + stiffness;
            if !stiffness.is_finite() || !denominator.is_finite() {
                return Err(TabulatedError::Numerical(
                    "axial rate return modulus is not representable",
                ));
            }
            e * (stiffness / denominator)
        };
        let tensor = self
            .tensor_strain(strain, &next)
            .map_err(|_| TabulatedError::State("axial rate return leaves the native J2 domain"))?;
        // The same tensor is passed to the native work/potential owner. Admit
        // its incoming elastic predictor as well as its returned history.
        backbone
            .validate_history(&tensor, &old_tensor)
            .map_err(|_| {
                TabulatedError::State("axial rate predictor leaves the native J2 domain")
            })?;
        let stress_pa = e * (strain - next.plastic_strain);
        if !stress_pa.is_finite() || !tangent_pa.is_finite() {
            return Err(TabulatedError::Numerical(
                "axial rate stress or tangent is not representable",
            ));
        }
        Ok(UniaxialJ2Response {
            strain: tensor,
            stress_pa,
            tangent_pa,
            state: next,
        })
    }

    /// Original elastic, exact table and independent Prager storage. Viscosity
    /// stores no recoverable energy, even when a retained state is overstressed.
    pub fn stored_energy(
        &self,
        strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<StoredEnergy, EnergyError> {
        self.backbone()
            .stored_energy(&self.tensor_strain(strain, state)?, &state.tensor_history())
    }

    /// Native discrete potential evaluated on the exact free-sided candidate.
    /// Its axial derivative is the physical axial stress; kinetic potential
    /// and physical viscous dissipation retain the tensor owner's distinction.
    pub fn stress_and_potential(
        &self,
        strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<(f64, f64), TabulatedError> {
        let next = self.response(strain, state)?;
        let (stress, potential) = self
            .backbone()
            .stress_and_potential(&next.strain, &state.tensor_history())?;
        Ok((stress[0], potential))
    }

    /// Original full-tensor work audit. Yield and viscosity contribute physical
    /// heat; exact table/Prager storage and backward-Euler loss remain separate.
    /// The returned native state is checked against the axial candidate under
    /// the supplied stress allowance, including every free-side component.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialJ2State,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        let old_tensor = self.tensor_strain(old_strain, old_state)?;
        let next = self.response(strain, old_state).map_err(|_| {
            EnergyError::Invalid("axial tabulated rate return is outside its native domain")
        })?;
        let increment = self.backbone().advance_with_energy(
            &old_tensor,
            &old_state.tensor_history(),
            &next.strain,
            tolerances,
        )?;
        let mut residual = (increment.stress_pa[0] - next.stress_pa).abs();
        for stress in &increment.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = self.backbone().backbone().backbone().elastic.youngs;
        let expected = next.state.tensor_history();
        for (actual, expected) in increment
            .state
            .plastic_strain
            .iter()
            .zip(expected.plastic_strain)
        {
            residual = residual.max(e * (actual - expected).abs());
        }
        residual = residual.max(e * (increment.state.alpha - expected.alpha).abs());
        if !residual.is_finite() || residual > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        Ok(increment)
    }
}

impl Uniaxial for UniaxialTabulatedOverstress {
    type State = UniaxialJ2State;

    fn initial_state(&self) -> Self::State {
        UniaxialJ2State {
            plastic_strain: 0.0,
            equivalent_plastic_strain: 0.0,
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
        self.response(strain, state).map_or(
            UniaxialJ2State {
                plastic_strain: f64::NAN,
                equivalent_plastic_strain: f64::NAN,
            },
            |r| r.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone().admissibility()
    }
    fn card(&self) -> ModelCard {
        self.backbone().card()
    }
}
