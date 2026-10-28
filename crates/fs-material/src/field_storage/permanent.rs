//! Reversible isotropic recoil: B = mu H + Br n at a fixed material state.
//!
//! Absolute recoil permeability, remanence and ALL operating limits are source
//! properties. Coercivity is not substituted for a reversible knee or used to
//! invent permeability. Direction belongs to specimen geometry, not chemistry.
//! Temperature refresh constructs a new pinned snapshot; it does not account
//! for parameter work or irreversible demagnetization during a live trajectory.

use super::{FieldStorageError, MAGNETIC_PERMEABILITY_DIMS, positive};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

/// Source-declared absolute reversible recoil slope [H/m]. A generic magnetic
/// permeability or Br/BHmax-derived model coefficient is not substituted.
pub const RECOIL_PERMEABILITY_PROPERTY: &str = "recoil_magnetic_permeability";
/// Remanent flux-density magnitude [T], never a relative permeability.
pub const REMANENT_FLUX_DENSITY_PROPERTY: &str = "remanent_flux_density";
/// Tesla in the canonical dimensional basis.
pub const FLUX_DENSITY_DIMS: Dims = Dims([0, 1, -2, 0, -1, 0]);
/// Signed lower parallel H limit of this reversible recoil model [A/m].
pub const RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY: &str = "recoil_minimum_parallel_field_strength";
/// Signed upper parallel H limit of this reversible recoil model [A/m].
pub const RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY: &str = "recoil_maximum_parallel_field_strength";
/// Maximum transverse H magnitude admitted by this recoil model [A/m].
pub const RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY: &str =
    "recoil_maximum_transverse_field_strength";
/// Ampere per metre in the canonical dimensional basis.
pub const FIELD_STRENGTH_DIMS: Dims = Dims([-1, 0, 0, 0, 1, 0]);

/// Immutable source-selected recoil law and its independently supplied limits.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedRecoilMaterial {
    material: ResolvedMaterialStatePoint,
    permeability_h_m: f64,
    remanence_t: f64,
    parallel_range_a_m: [f64; 2],
    maximum_transverse_a_m: f64,
}

/// Physical H components and signed distances to the original source limits.
/// A nonnegative margin is source-model admissibility, not a certified bound
/// on the numerical field or a guarantee against irreversible material change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecoilOperatingPoint {
    /// H projected along the prescribed remanence direction [A/m].
    pub parallel_field_a_m: f64,
    /// Magnitude of the perpendicular H component [A/m].
    pub transverse_field_a_m: f64,
    /// Distance above the source's lower parallel limit [A/m].
    pub lower_parallel_margin_a_m: f64,
    /// Distance below the source's upper parallel limit [A/m].
    pub upper_parallel_margin_a_m: f64,
    /// Distance below the source's transverse-magnitude limit [A/m].
    pub transverse_margin_a_m: f64,
}

impl ResolvedRecoilMaterial {
    /// Select the complete recoil law at one original physical query point.
    /// A remanence-only catalog entry is insufficient: all five quantities,
    /// their SI units and their fixed-condition domains must be supported.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FieldStorageError> {
        let requirements = [
            (
                RECOIL_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                REMANENT_FLUX_DENSITY_PROPERTY,
                FLUX_DENSITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::Finite,
            ),
            (
                RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::Finite,
            ),
            (
                RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::NonNegative,
            ),
        ]
        .into_iter()
        .map(|(name, dims, admissibility)| {
            ScalarPropertyRequirement::try_new(name, dims, admissibility)
        })
        .collect::<Result<Vec<_>, _>>()?;
        Self::try_from_material(resolve_material_state_point(
            card,
            point,
            &requirements,
            selection,
        )?)
    }

    /// Adopt a complete material bundle containing explicitly declared recoil
    /// data. Its slope is never inferred from generic permeability, Br/Hc,
    /// maximum energy product or a material name.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, FieldStorageError> {
        let scalar = |name: &str, dims: Dims| -> Result<f64, FieldStorageError> {
            let p = material.property(name).ok_or(FieldStorageError::Invalid(
                "incomplete or ambiguous permanent-magnet recoil data",
            ))?;
            if p.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !p.is_static_at_fixed_external_conditions()
                || !p.value_si().is_finite()
            {
                return Err(FieldStorageError::Invalid(
                    "recoil data must be finite absolute SI scalars at supported fixed external conditions",
                ));
            }
            Ok(p.value_si())
        };
        let permeability_h_m = positive(
            scalar(RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)?,
            "absolute recoil permeability must be positive and finite",
        )?;
        let remanence_t = positive(
            scalar(REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS)?,
            "remanent flux density must be positive and finite",
        )?;
        let parallel_range_a_m = [
            scalar(RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY, FIELD_STRENGTH_DIMS)?,
            scalar(RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY, FIELD_STRENGTH_DIMS)?,
        ];
        let maximum_transverse_a_m = scalar(
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
        )?;
        if parallel_range_a_m[0] >= parallel_range_a_m[1]
            || parallel_range_a_m[0] > 0.0
            || parallel_range_a_m[1] < 0.0
            || maximum_transverse_a_m < 0.0
        {
            return Err(FieldStorageError::Invalid(
                "recoil envelope needs ordered parallel limits containing H=0 and a nonnegative transverse limit",
            ));
        }
        Ok(Self {
            material,
            permeability_h_m,
            remanence_t,
            parallel_range_a_m,
            maximum_transverse_a_m,
        })
    }
    /// Original card identity, source choices, complete conditions and receipts.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Absolute reversible recoil slope [H/m].
    pub const fn permeability_h_m(&self) -> f64 {
        self.permeability_h_m
    }
    /// Original remanent-flux magnitude [T].
    pub const fn remanence_t(&self) -> f64 {
        self.remanence_t
    }
    /// Signed parallel H limits [A/m].
    pub const fn parallel_range_a_m(&self) -> [f64; 2] {
        self.parallel_range_a_m
    }
    /// Maximum perpendicular H magnitude [A/m].
    pub const fn maximum_transverse_a_m(&self) -> f64 {
        self.maximum_transverse_a_m
    }

    /// Re-evaluate all properties together, retaining every original source
    /// claim and coordinate meaning. This returns a new fixed-state law only.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_material(self.material.at_query_point(card, point)?)
    }

    /// Admit a continuing electrical/material state without changing the
    /// original reversible recoil law or any of its operating limits.
    /// The complete original bundle, source claims and coordinate meanings
    /// remain pinned. Conductivity may change with the owned temperature, but
    /// a temperature-dependent recoil property requires its own magnetic
    /// parameter-work and entropy law, even if two endpoint values coincide.
    pub fn admit_fixed_material_refresh(
        &self,
        refreshed: &ResolvedMaterialStatePoint,
    ) -> Result<(), FieldStorageError> {
        super::admit_source_temperature_refresh(&self.material, refreshed)?;
        for name in [
            RECOIL_PERMEABILITY_PROPERTY,
            REMANENT_FLUX_DENSITY_PROPERTY,
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
        ] {
            let original = self
                .material
                .property(name)
                .expect("prepared recoil source");
            let next = refreshed.property(name).ok_or(FieldStorageError::Invalid(
                "refreshed material must retain every original recoil property",
            ))?;
            if original.value_si().to_bits() != next.value_si().to_bits()
                || self
                    .material
                    .query_point()
                    .iter()
                    .zip(refreshed.query_point())
                    .any(|((axis, before), (_, after))| {
                        before.to_bits() != after.to_bits()
                            && original.source_axes().iter().any(|(name, _)| name == axis)
                    })
            {
                return Err(FieldStorageError::Invalid(
                    "temperature-dependent recoil requires magnetic parameter work and entropy",
                ));
            }
        }
        Ok(())
    }

    /// Check the computed field in specimen-aligned coordinates. A point outside
    /// any source limit refuses; no magnetization state is clamped or changed.
    pub fn admit_operating_point(
        &self,
        parallel_field_a_m: f64,
        transverse_field_a_m: f64,
    ) -> Result<RecoilOperatingPoint, FieldStorageError> {
        if !parallel_field_a_m.is_finite()
            || !transverse_field_a_m.is_finite()
            || transverse_field_a_m < 0.0
            || parallel_field_a_m < self.parallel_range_a_m[0]
            || parallel_field_a_m > self.parallel_range_a_m[1]
            || transverse_field_a_m > self.maximum_transverse_a_m
        {
            return Err(FieldStorageError::Invalid(
                "computed magnetic field leaves the source-backed reversible recoil envelope",
            ));
        }
        let margins = [
            parallel_field_a_m - self.parallel_range_a_m[0],
            self.parallel_range_a_m[1] - parallel_field_a_m,
            self.maximum_transverse_a_m - transverse_field_a_m,
        ];
        if !margins.iter().all(|x| x.is_finite()) {
            return Err(FieldStorageError::Invalid(
                "recoil operating margin is unrepresentable",
            ));
        }
        Ok(RecoilOperatingPoint {
            parallel_field_a_m,
            transverse_field_a_m,
            lower_parallel_margin_a_m: margins[0],
            upper_parallel_margin_a_m: margins[1],
            transverse_margin_a_m: margins[2],
        })
    }
}
