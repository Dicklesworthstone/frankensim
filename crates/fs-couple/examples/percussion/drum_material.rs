//! A frozen isotropic head material from the shared source resolver.
use super::{Error, Head};
use fs_couple::{material_source::{MaterialRecord, plate_section}, thin_plate::PlateMaterialModel};
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::state_point::{ResolvedMaterialStatePoint, DENSITY_PROPERTY,
    POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY};
use std::sync::Arc;

/// Original immutable card, query and resolved coefficients used by one head.
/// Geometry and installed tension remain independent explicit specimen inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct HeadMaterial {
    record: MaterialRecord,
    pack: Arc<NormalizedMaterialCardPack>,
    state: ResolvedMaterialStatePoint,
}

impl HeadMaterial {
    pub(super) fn resolve(record: MaterialRecord, pack: Arc<NormalizedMaterialCardPack>,
        thickness_m: f64, tension_n_m: f64) -> Result<(Head, Self), Error>
    {
        let card = record.card(&pack)?;
        let point = record.query_point(card)?;
        let (_, state) = plate_section(card, &point, PlateMaterialModel::Isotropic, thickness_m, 0.0)?;
        // A single frozen elastic pencil cannot stand in for a dispersive or
        // strain/rate-dependent law. Temperature/moisture are fixed conditions.
        if state.properties().iter().any(|p| !p.is_static_at_fixed_external_conditions()) {
            return Err("head_material requires static isotropic density/moduli at fixed external conditions; frequency, stress, strain and rate dependence need their own law".into());
        }
        let source = Self { record, pack, state };
        let head = Head { thickness_m, tension_n_m, damping_ratio: 0.0,
            young_pa: source.value(YOUNG_MODULUS_PROPERTY)?,
            poisson: source.value(POISSON_RATIO_PROPERTY)?,
            density_kg_m3: source.value(DENSITY_PROPERTY)? };
        Ok((head, source))
    }

    fn value(&self, name: &str) -> Result<f64, Error> {
        self.state.property(name).map(|p| p.value_si())
            .ok_or_else(|| format!("head material lost resolved property {name}").into())
    }

    /// Retained source declaration, including explicit physical query axes.
    pub fn record(&self) -> &MaterialRecord { &self.record }
    /// Exact input pack, including original claims and model associations.
    pub fn pack(&self) -> &NormalizedMaterialCardPack { &self.pack }
    /// The one state used for both head inertia and elastic stiffness.
    pub fn state(&self) -> &ResolvedMaterialStatePoint { &self.state }
    pub fn receipt(&self) -> String { self.record.provenance(&self.pack) }

    pub(super) fn validate(&self, head: &Head) -> Result<(), Error> {
        for (name, value) in [(YOUNG_MODULUS_PROPERTY, head.young_pa),
            (POISSON_RATIO_PROPERTY, head.poisson), (DENSITY_PROPERTY, head.density_kg_m3)] {
            if value.to_bits() != self.value(name)?.to_bits() {
                return Err(format!("head_material {name} differs from its original source; resolve the replacement card instead of overwriting a coefficient").into());
            }
        }
        if head.damping_ratio != 0.0 {
            return Err("head_material cannot acquire an authored modal damping ratio".into());
        }
        Ok(())
    }
}
