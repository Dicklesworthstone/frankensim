//! Directional tensile limits in the material axes of an orthotropic solid.
//!
//! These are the three tensile branches of the maximum-stress criterion:
//! `max_i(max(0, sigma_ii) / X_ti)`. The axes belong to the material, not to
//! the principal stresses of the current load. Rotating the load relative to
//! wood grain or reinforcement can therefore change the governing limit.
//!
//! Compression, shear failure and interaction between stress components are
//! deliberately not assessed. In particular, zero tensile utilization under
//! material-axis pure shear is not a statement that the specimen is safe.
//! No damage evolution, fracture energy or automatic stiffness loss is implied.
//! Reference: Abaqus, “Generalized Damage Initiation Criteria for Anisotropic
//! Materials”, maximum stress criterion:
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-damageinitgeneral.htm>.

use super::invalid;
use crate::Voigt;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};

/// Three explicitly selected tensile limits in the original material axes.
///
/// Each axis requires its own named source on the same material card and at
/// the same complete physical query. A missing transverse strength is never
/// inferred from a longitudinal value or from the elastic moduli. The caller
/// declares how axes 1, 2 and 3 correspond to the source's physical directions.
#[derive(Clone, Debug, PartialEq)]
pub struct OrthotropicTensileLimit {
    resolved: ResolvedMaterialStatePoint,
    limit_properties: [String; 3],
    limits_pa: [f64; 3],
}

impl OrthotropicTensileLimit {
    /// Resolve three named tensile-strength claims [Pa], in material-axis order.
    /// No source selection, safety factor or engineering-to-true-stress
    /// conversion is inferred. Use small-strain stresses or explicitly
    /// compatible Cauchy strengths, just as for the isotropic tensile owner.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        limit_properties: [&str; 3],
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        let requirement = |name| {
            ScalarPropertyRequirement::try_new(
                name,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
        };
        let requirements = [
            requirement(limit_properties[0])?,
            requirement(limit_properties[1])?,
            requirement(limit_properties[2])?,
        ];
        let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
        Self::try_from_resolved(resolved, limit_properties)
    }

    /// Bind the three tensile claims from an existing immutable state bundle.
    /// The original complete query and every source receipt remain intact.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        limit_properties: [&str; 3],
    ) -> Result<Self, MaterialStatePointError> {
        if limit_properties[0] == limit_properties[1]
            || limit_properties[1] == limit_properties[2]
            || limit_properties[2] == limit_properties[0]
        {
            return Err(invalid(
                "orthotropic tensile axes require three distinct named claims",
            ));
        }
        let mut limits_pa = [0.0; 3];
        for (axis, name) in limit_properties.iter().enumerate() {
            let source = resolved.property(name).ok_or_else(|| {
                invalid("every material axis requires its named tensile-strength claim")
            })?;
            if source.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS) {
                return Err(invalid(
                    "directional tensile strength must have the declared pressure quantity in Pa",
                ));
            }
            if !source.is_static_at_fixed_external_conditions() {
                return Err(invalid(
                    "directional tensile strength requires fixed external conditions; unowned strain, rate or frequency refuses",
                ));
            }
            let value = source.value_si();
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid(
                    "directional tensile strength must be positive and finite",
                ));
            }
            limits_pa[axis] = value;
        }
        Ok(Self {
            resolved,
            limit_properties: limit_properties.map(str::to_owned),
            limits_pa,
        })
    }

    /// Refresh the exact original claims at new supported external conditions.
    /// Unsupported queries or a different material card refuse atomically.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(
            self.resolved.at_query_point(card, point)?,
            self.limit_properties(),
        )
    }

    /// Complete retained material state and source selections.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Named tensile properties in the original material-axis order.
    #[must_use]
    pub fn limit_properties(&self) -> [&str; 3] {
        self.limit_properties.each_ref().map(String::as_str)
    }

    /// Tensile limits [Pa] along material axes 1, 2 and 3.
    #[must_use]
    pub const fn limits_pa(&self) -> [f64; 3] {
        self.limits_pa
    }

    /// Original source receipts in material-axis order.
    #[must_use]
    pub fn sources(&self) -> [&ResolvedScalarProperty; 3] {
        self.limit_properties.each_ref().map(|name| {
            self.resolved
                .property(name)
                .expect("admitted directional tensile source")
        })
    }

    /// Assess tension-positive tensor-shear stress in the MATERIAL frame,
    /// ordered `[11, 22, 33, 12, 23, 31]`, without allocation or eigensolves.
    ///
    /// A spatial consumer must transform its stress into these same axes.
    /// All components must be finite, even though only normal tension is
    /// assessed. A nonzero demand whose utilization cannot be represented
    /// refuses instead of silently reporting zero.
    pub fn assess_material_stress(
        &self,
        stress_pa: &Voigt,
    ) -> Result<OrthotropicTensileAssessment, MaterialStatePointError> {
        if stress_pa.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "directional tensile assessment requires a finite material-frame stress",
            ));
        }
        let mut tensile_utilizations = [0.0; 3];
        let mut governing_axis = 0;
        for axis in 0..3 {
            let demand = stress_pa[axis].max(0.0);
            let utilization = demand / self.limits_pa[axis];
            if !utilization.is_finite() || demand > 0.0 && utilization == 0.0 {
                return Err(invalid(
                    "directional tensile utilization is unrepresentable",
                ));
            }
            tensile_utilizations[axis] = utilization;
            if utilization > tensile_utilizations[governing_axis] {
                governing_axis = axis;
            }
        }
        Ok(OrthotropicTensileAssessment {
            material_stress_pa: *stress_pa,
            tensile_utilizations,
            governing_axis,
        })
    }
}

/// Tensile branches of the maximum-stress criterion at one material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicTensileAssessment {
    material_stress_pa: Voigt,
    tensile_utilizations: [f64; 3],
    governing_axis: usize,
}

impl OrthotropicTensileAssessment {
    /// Original material-frame tensor-shear stress [Pa], including unassessed
    /// compression and shear, in `[11, 22, 33, 12, 23, 31]` order.
    #[must_use]
    pub const fn material_stress_pa(&self) -> Voigt {
        self.material_stress_pa
    }

    /// Signed normal stresses [Pa] in the original material axes.
    #[must_use]
    pub fn material_normal_stresses_pa(&self) -> [f64; 3] {
        [
            self.material_stress_pa[0],
            self.material_stress_pa[1],
            self.material_stress_pa[2],
        ]
    }

    /// Positive normal stress divided by its own directional limit.
    #[must_use]
    pub const fn tensile_utilizations(&self) -> [f64; 3] {
        self.tensile_utilizations
    }

    /// Zero-based material axis with the largest utilization; the first wins
    /// exact ties, including the all-compressive/zero-tension case.
    #[must_use]
    pub const fn governing_axis(&self) -> usize {
        self.governing_axis
    }

    /// Largest of the three tensile utilizations. This is not a combined
    /// tensile/compressive/shear failure index.
    #[must_use]
    pub fn utilization(&self) -> f64 {
        self.tensile_utilizations[self.governing_axis]
    }

    /// Whether a directional tensile limit is reached, including equality.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.utilization() >= 1.0
    }
}
