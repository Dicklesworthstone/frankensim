//! Explicit physical volumetric viscosity for the native D3Q19 liquid.

use super::{FreeSurface3, MaterialFreeSurface3, MaterialFreeSurfaceScaling, MaterialLbmError};
use crate::d3q19::CollisionModel3;
use crate::d3q19::freesurface3::{FreeSurfaceError3, transport::MassLocation3};
use fs_exec::CancelGate;
use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
use fs_material::state_point::{ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility};
use fs_qty::QuantitySpec;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub(super) struct BulkViscosity<'a> {
    viscosity_pa_s: f64,
    relaxation_rate: f64,
    higher_order_rate: f64,
    // Borrow the original complete source and exact selected property; neither
    // coefficient evaluation nor thermal feedback makes a new claim selection.
    source: Option<(&'a ResolvedMaterialStatePoint, &'a ResolvedScalarProperty)>,
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    /// Author a uniform physical bulk viscosity [Pa s], independent of the
    /// liquid source's shear viscosity. Every 3-D constructor retains this
    /// coefficient through local shear/temperature updates and uses the same
    /// coefficient for `zeta (div u)^2` when molecular heating is requested.
    ///
    /// This explicitly supplied scalar is not inferred from numerical damping,
    /// nor a source law varying with temperature, pressure or composition.
    /// It controls volumetric viscous stress, not reversible pressure work or
    /// compressible thermodynamics. Higher moments use the authored numerical
    /// rate; their damping is not deposited as heat. The native strict `(0,2)`
    /// window requires positive representable zeta; zero refuses without a
    /// viscosity floor. The 2-D and SGS owners refuse this collision family.
    pub fn with_bulk_viscosity(
        mut self,
        bulk_viscosity_pa_s: f64,
        higher_order_rate: f64,
    ) -> Result<Self, MaterialLbmError> {
        let tau = crate::material::physical_bulk_tau(
            self.material.density_kg_m3(),
            self.velocity_scale,
            self.cell_size_m,
            bulk_viscosity_pa_s,
        ).map_err(MaterialLbmError::Invalid)?;
        let relaxation_rate = 1.0 / tau;
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0 / self.flow.tau,
            bulk_rate: relaxation_rate,
            higher_order_rate,
        }.validate().map_err(|_| MaterialLbmError::Invalid(
            "physical bulk collision rates are unrepresentable",
        ))?;
        self.bulk_viscosity = Some(BulkViscosity {
            viscosity_pa_s: bulk_viscosity_pa_s,
            relaxation_rate,
            higher_order_rate,
            source: None,
        });
        Ok(self)
    }

    /// Uniform authored or reference source bulk viscosity [Pa s]. A selected
    /// [`Self::bulk_viscosity_property`] evolves at each material query; this
    /// reference value is not its current local coefficient. None retains
    /// `zeta=2mu/3` from the current molecular shear viscosity.
    #[must_use]
    pub const fn bulk_viscosity_pa_s(&self) -> Option<f64> {
        match self.bulk_viscosity {
            Some(bulk) => Some(bulk.viscosity_pa_s),
            None => None,
        }
    }

    /// Bind one already selected positive Pa s property from this SAME source
    /// bundle. The original name, quantity, claim and complete query remain
    /// authoritative. Thermal/rheological owners refresh zeta with mu at every
    /// actual temperature, pressure and deviatoric-shear query; unsupported
    /// evolving coordinates refuse in frozen consumers. No value, source table
    /// or bulk/shear ratio is inferred when a claim is absent.
    pub fn with_bulk_viscosity_source(
        self,
        property: &str,
        higher_order_rate: f64,
    ) -> Result<Self, MaterialLbmError> {
        let material = self.material.material();
        let selected = material.property(property).ok_or(MaterialLbmError::Invalid(
            "bulk viscosity requires its selected property in the original material bundle",
        ))?;
        let requirement = selected.requirement();
        if requirement.quantity() != QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS)
            || requirement.admissibility() != ScalarAdmissibility::StrictlyPositive
            || requirement.hardness_test().is_some()
            || requirement.elastic_component().is_some()
            || requirement.strain_component().is_some()
            || requirement.stress_component().is_some()
        {
            return Err(MaterialLbmError::Invalid(
                "bulk viscosity requires an exact positive scalar Pa s source requirement",
            ));
        }
        let mut result = self.with_bulk_viscosity(selected.value_si(), higher_order_rate)?;
        result.bulk_viscosity.as_mut().expect("bound bulk").source = Some((material, selected));
        Ok(result)
    }

    /// Exact retained bulk-viscosity property name, when source evolution is selected.
    #[must_use]
    pub fn bulk_viscosity_property(&self) -> Option<&'a str> {
        self.bulk_viscosity.and_then(|bulk| bulk.source.map(|(_, source)| source.requirement().name()))
    }

    pub(in crate::material::free_surface) fn bulk_viscosity_at(
        &self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<Option<f64>, MaterialLbmError> {
        let Some(bulk) = self.bulk_viscosity else { return Ok(None) };
        let Some((original, source)) = bulk.source else { return Ok(Some(bulk.viscosity_pa_s)) };
        let value = material.property(source.requirement().name()).filter(|value| {
            material.card_identity() == original.card_identity()
                && value.requirement() == source.requirement()
                && value.answer().receipt.selected == source.answer().receipt.selected
        }).ok_or(MaterialLbmError::Invalid(
            "local bulk viscosity must retain the original material requirement and selected claim",
        ))?.value_si();
        if !value.is_finite() || value <= 0.0 {
            return Err(MaterialLbmError::Invalid("source bulk viscosity must be finite and positive"));
        }
        Ok(Some(value))
    }

    pub(in crate::material::free_surface) fn bulk_relaxation_rate_at(
        &self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<Option<f64>, MaterialLbmError> {
        let Some(bulk) = self.bulk_viscosity else { return Ok(None) };
        if bulk.source.is_none() { return Ok(Some(bulk.relaxation_rate)) }
        let value = self.bulk_viscosity_at(material)?.expect("bound bulk source");
        let tau = crate::material::physical_bulk_tau(
            self.material.density_kg_m3(), self.velocity_scale, self.cell_size_m, value,
        ).map_err(MaterialLbmError::Invalid)?;
        let rate = 1.0 / tau;
        if !rate.is_finite() || rate <= 0.0 || rate >= 2.0 {
            return Err(MaterialLbmError::Invalid("source bulk relaxation is unrepresentable"));
        }
        Ok(Some(rate))
    }

    pub(super) fn admit_frozen_bulk_source(&self) -> Result<(), MaterialLbmError> {
        if let Some((_, source)) = self.bulk_viscosity.and_then(|bulk| bulk.source) {
            fs_material::liquid::surface::require_frozen_coefficient_support(source)
                .map_err(|_| MaterialLbmError::Invalid(
                    "source bulk viscosity requires its evolving temperature/shear/pressure owner",
                ))?;
        }
        Ok(())
    }

    pub(super) fn source_bulk_relaxations<'m>(
        &self,
        materials: impl IntoIterator<Item = (MassLocation3, &'m ResolvedMaterialStatePoint)>,
        gate: &CancelGate,
    ) -> Result<Option<BTreeMap<(u64, usize), f64>>, MaterialLbmError> {
        if self.bulk_viscosity_property().is_none() { return Ok(None) }
        let mut field = BTreeMap::new();
        for (location, material) in materials {
            if gate.is_requested() { return Err(FreeSurfaceError3::Cancelled.into()) }
            if let MassLocation3::Cell(key, lane) = location {
                let rate = self.bulk_relaxation_rate_at(material)?.expect("bound bulk source");
                field.insert((key, lane), rate);
            }
        }
        if gate.is_requested() { return Err(FreeSurfaceError3::Cancelled.into()) }
        Ok(Some(field))
    }

    pub(super) fn apply_source_bulk_field(
        &self,
        lattice: &mut FreeSurface3,
        field: BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        let bulk = self.bulk_viscosity.ok_or(MaterialLbmError::Invalid("missing source bulk binding"))?;
        if bulk.source.is_none() {
            return Err(MaterialLbmError::Invalid("local bulk field requires its retained material source"));
        }
        lattice.select_bulk_relaxations(bulk.relaxation_rate, bulk.higher_order_rate, field, gate)?;
        Ok(())
    }

    pub(super) fn apply_collision(&self, lattice: &mut FreeSurface3) -> Result<(), MaterialLbmError> {
        if let Some(bulk) = self.bulk_viscosity {
            lattice.set_central_moment_bulk(bulk.relaxation_rate, bulk.higher_order_rate)
                .map_err(|_| MaterialLbmError::Invalid("physical bulk collision rates are unrepresentable"))?;
        }
        Ok(())
    }
}

impl MaterialFreeSurface3<'_> {
    /// Bind independent physical bulk viscosity without advancing the native
    /// trajectory. Retains shear viscosity, source receipts, populations, mass
    /// and force frame. The new coefficient changes how retained trace stress
    /// is interpreted; use the scaling builder before initialization when an
    /// authored initial velocity gradient must be preserved.
    /// Invalid coefficients leave both physical and native bindings unchanged.
    pub fn set_bulk_viscosity(
        &mut self,
        bulk_viscosity_pa_s: f64,
        higher_order_rate: f64,
    ) -> Result<(), MaterialLbmError> {
        let scaling = self.scaling.with_bulk_viscosity(bulk_viscosity_pa_s, higher_order_rate)?;
        scaling.apply_collision(&mut self.lattice)?;
        self.scaling = scaling;
        Ok(())
    }

    /// Bind a retained bulk claim to a frozen trajectory. Claims requiring
    /// evolving shear or pressure refuse; use the corresponding material owner.
    /// Existing populations are preserved, with the same trace-rebind semantics
    /// as [`Self::set_bulk_viscosity`].
    pub fn set_bulk_viscosity_source(
        &mut self,
        property: &str,
        higher_order_rate: f64,
    ) -> Result<(), MaterialLbmError> {
        let scaling = self.scaling.with_bulk_viscosity_source(property, higher_order_rate)?;
        scaling.admit_frozen_bulk_source()?;
        scaling.apply_collision(&mut self.lattice)?;
        self.scaling = scaling;
        Ok(())
    }
}
