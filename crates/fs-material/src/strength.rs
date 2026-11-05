//! Explicit isotropic maximum-principal-stress (Rankine tensile) assessment.
//!
//! The tensile demand is max(0, largest principal Cauchy stress), including
//! hydrostatic tension. Compressive principal stresses never consume tensile
//! capacity. The caller selects the actual named strength claim and the
//! isotropic criterion; a modulus, hardness, yield observation, or typical
//! ultimate strength is not silently substituted for a design allowable.
//!
//! Use with small-strain stresses, or an explicitly compatible Cauchy-strength
//! observation. Engineering ultimate tensile strength is based on original
//! area and is not a true-stress failure law for a necking finite-strain body.
//! This check supplies no anisotropic failure surface, plastic evolution,
//! toughness, crack growth, fatigue, fracture energy, or post-failure transfer.
//! A crossed limit reports a model limit; it never deletes physical state.
//!
//! The source resolver retains condition and claim receipts. Strength is held
//! at those external conditions until explicitly refreshed. Strain/rate/frequency
//! dependence with no physical owner refuses. Reference: COMSOL Structural
//! Mechanics theory, “Rankine Cutoff”:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.033.html>.

mod finite;
pub use finite::ResolvedTensileHyperelastic;

use crate::Voigt;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};

/// An explicitly selected isotropic tensile limit with its original source state.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicTensileLimit {
    resolved: ResolvedMaterialStatePoint,
    limit_property: String,
    limit_pa: f64,
}

impl IsotropicTensileLimit {
    /// Resolve the caller's named strength in Pa at one complete physical point.
    /// Source ambiguity, missing data, nonpositive strength and invalid domains
    /// are errors; no numeric fallback or implicit safety factor is supplied.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        limit_property: &str,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let requirement = ScalarPropertyRequirement::try_new(
            limit_property,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let resolved = resolve_material_state_point(card, point, &[requirement], selection)?;
        Self::try_from_resolved(resolved, limit_property)
    }

    /// Bind one named member of an existing state bundle, preserving every
    /// source receipt. Calling this constructor explicitly selects an isotropic
    /// tensile criterion; it does not infer isotropy from a scalar datum.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        limit_property: &str,
    ) -> Result<Self, MaterialStatePointError> {
        let limit = resolved.property(limit_property)
            .ok_or_else(|| invalid("tensile assessment requires its named strength claim"))?;
        if limit.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS) {
            return Err(invalid("tensile strength must have the declared pressure quantity in Pa"));
        }
        if !limit.is_static_at_fixed_external_conditions() {
            return Err(invalid("tensile strength requires static fixed external conditions; unowned strain, rate or frequency refuses"));
        }
        let limit_pa = limit.value_si();
        if !limit_pa.is_finite() || limit_pa <= 0.0 {
            return Err(invalid("tensile strength must be positive and finite"));
        }
        Ok(Self {
            resolved,
            limit_property: limit_property.into(),
            limit_pa,
        })
    }

    /// Re-evaluate the same original claims on the same material card at new
    /// external conditions. A failed refresh leaves this immutable state intact.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(self.resolved.at_query_point(card, point)?, &self.limit_property)
    }

    /// Complete original material bundle and condition/source identities.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Exact strength property selected by the caller.
    #[must_use]
    pub fn limit_property(&self) -> &str {
        &self.limit_property
    }

    /// Strength at the retained physical state [Pa], without a safety factor.
    #[must_use]
    pub const fn limit_pa(&self) -> f64 {
        self.limit_pa
    }

    /// Original selected scalar claim and its usage receipt.
    #[must_use]
    pub fn source(&self) -> &ResolvedScalarProperty {
        self.resolved.property(&self.limit_property).expect("admitted tensile source")
    }

    /// Compare actual tensile-positive stress [Pa] in tensor-shear Voigt order
    /// [xx, yy, zz, xy, yz, zx]. This bounded calculation allocates nothing.
    /// Nonfinite stress, unrepresentable eigenvalues or utilization refuse.
    pub fn assess(&self, stress_pa: &Voigt) -> Result<TensileAssessment, MaterialStatePointError> {
        let principal_stresses_pa = principal_stresses(stress_pa)?;
        let tensile_stress_pa = principal_stresses_pa[2].max(0.0);
        let utilization = tensile_stress_pa / self.limit_pa;
        if !utilization.is_finite() || (tensile_stress_pa > 0.0 && utilization == 0.0) {
            return Err(invalid("tensile utilization is unrepresentable"));
        }
        Ok(TensileAssessment {
            principal_stresses_pa,
            tensile_stress_pa,
            utilization,
        })
    }
}

/// Local stress/strength observation under the explicitly selected criterion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensileAssessment {
    principal_stresses_pa: [f64; 3],
    tensile_stress_pa: f64,
    utilization: f64,
}

impl TensileAssessment {
    /// All three principal stresses [Pa], ascending and tension positive.
    #[must_use]
    pub const fn principal_stresses_pa(&self) -> [f64; 3] {
        self.principal_stresses_pa
    }

    /// Positive part of the largest principal stress [Pa].
    #[must_use]
    pub const fn tensile_stress_pa(&self) -> f64 {
        self.tensile_stress_pa
    }

    /// Tensile demand divided by the caller's sourced limit.
    #[must_use]
    pub const fn utilization(&self) -> f64 {
        self.utilization
    }

    /// Equality or exceedance under floating-point arithmetic, without tolerance.
    /// A crossed limit does not establish fracture initiation or its evolution.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.utilization >= 1.0
    }
}

/// Principal values of a finite tensor-shear Voigt stress, in ascending order.
/// Reuses the allocation-free, scaled fixed-size fs-la Jacobi kernel.
pub fn principal_stresses(stress_pa: &Voigt) -> Result<[f64; 3], MaterialStatePointError> {
    let [xx, yy, zz, xy, yz, zx] = *stress_pa;
    fs_la::eigen::symmetric_eigenvalues3(&[[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
        .map_err(|_| invalid("principal stress requires a finite tensor and representable converged eigenvalues"))
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
