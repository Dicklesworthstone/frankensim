//! Principal material conductivity on the original finite thermal links.
//!
//! The authored unit vector is the retained temperature-gradient direction in
//! that cell's material axes. Its Galerkin coefficient is n^T K n. This is not
//! the effective conductivity of an arbitrarily oriented, laterally insulated
//! anisotropic rod, nor a replacement for the spatial tensor conduction owner.
//! Scalar source queries, thermal state and series resistance retain their
//! existing owners; this module only binds source roles and link geometry.

use super::*;
use fs_blake3::ContentHash;

/// Principal source roles and physical gradient orientation on one link side.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalConductionSide {
    /// Three distinct scalar static conductivity keys, in material-axis order.
    pub principal_properties: [PropertyKey; 3],
    /// Explicit original claim choices for those three properties.
    pub selection: MaterialPropertySelection,
    /// Unit gradient direction in the declared orthonormal material frame.
    /// Reversal leaves the conductance unchanged. No normalization is inferred.
    pub direction_in_material_axes: [f64; 3],
    /// Nonzero identity of the physical principal axes, shared by this cell's links.
    pub material_frame: ContentHash,
}

/// A source-derived replacement for one original link's two scalar k values.
/// Distances, endpoint order and area remain with the original thermal network.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalConductionLink {
    pub link: u64,
    /// Original thermal-link a, b order, independent of material/cell sorting.
    pub sides: [DirectionalConductionSide; 2],
}

pub(super) struct DirectionalPlan {
    links: BTreeMap<u64, (DirectionalConductionLink, [CellId; 2])>,
    cells: BTreeMap<CellId, DirectionalConductionSide>,
}

impl DirectionalPlan {
    pub(super) fn new(
        network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
        inputs: &[DirectionalConductionLink], cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        if inputs.is_empty() || inputs.len() > network.links.len() {
            return Err(invalid("directional conduction needs a bounded nonempty set of original links"));
        }
        let mut material_links = std::collections::BTreeSet::new();
        for link in network.links {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            material_links.insert(link.link);
        }
        let mut links = BTreeMap::new();
        let mut cells = BTreeMap::<CellId, DirectionalConductionSide>::new();
        for input in inputs {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            let link = thermal.links().get(&input.link)
                .ok_or(invalid("directional conduction link is inactive"))?;
            if !material_links.contains(&input.link)
                || links.insert(input.link, (input.clone(), [link.a, link.b])).is_some() {
                return Err(invalid("directional conduction needs unique existing material link bindings"));
            }
            for (id, side) in [link.a, link.b].into_iter().zip(&input.sides) {
                let norm = side.direction_in_material_axes.iter().map(|x| x * x).sum::<f64>();
                if side.material_frame == ContentHash([0; 32]) || !norm.is_finite()
                    || (norm - 1.0).abs() > 256.0 * f64::EPSILON {
                    return Err(invalid("directional conduction requires an identified material frame and a finite unit gradient direction"));
                }
                for (i, key) in side.principal_properties.iter().enumerate() {
                    let quantity = key.quantity();
                    if (quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
                        && quantity != QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity, ValueForm::Static)))
                        || side.principal_properties[..i].iter().any(|other| other.name() == key.name()) {
                        return Err(invalid("three distinct scalar static thermal conductivity roles are required"));
                    }
                }
                if let MaterialPropertySelection::PinnedByProperty(pins) = &side.selection {
                    if pins.len() != 3 || side.principal_properties.iter().any(|key|
                        pins.iter().filter(|(name, _)| name == key.name()).count() != 1) {
                        return Err(invalid("principal conductivity selection must pin exactly its three source roles"));
                    }
                }
                if let Some(original) = cells.get(&id) {
                    if original.principal_properties != side.principal_properties
                        || original.selection != side.selection || original.material_frame != side.material_frame {
                        return Err(invalid("all links of one cell must retain the same principal conductivity roles, claims and frame"));
                    }
                } else {
                    cells.insert(id, side.clone());
                }
            }
        }
        Ok(Self { links, cells })
    }

    pub(super) fn admit_topology(
        &self, thermal: &DepositionSimulation<'_>,
        check: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialNetworkError> {
        for (&id, (_, cells)) in &self.links {
            check()?;
            if thermal.links().get(&id).is_none_or(|link| [link.a, link.b] != *cells) {
                return Err(MaterialNetworkError::Invalid("directional conductivity cannot migrate to different material endpoints"));
            }
        }
        Ok(())
    }

    pub(super) fn requests(
        &self, cell: CellId, requested: &mut Vec<(ScalarPropertyRequirement, MaterialPropertySelection)>,
    ) -> Result<bool, MaterialNetworkError> {
        let Some(side) = self.cells.get(&cell) else { return Ok(false); };
        for key in &side.principal_properties {
            if requested.iter().any(|(r, _)| r.name() == key.name()) {
                return Err(MaterialNetworkError::Invalid("principal conductivity roles must be distinct from other transport roles"));
            }
            let requirement = ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                .map_err(|e| MaterialNetworkError::Material(ThermalMaterialError::Material(e)))?;
            // The canonical resolver requires each single-property selection to
            // contain exactly that role; final complete bundles retain all pins.
            let selection = match &side.selection {
                MaterialPropertySelection::PinnedByProperty(pins) =>
                    MaterialPropertySelection::PinnedByProperty(pins.iter()
                        .filter(|(name, _)| name == key.name()).cloned().collect()),
                other => other.clone(),
            };
            requested.push((requirement, selection));
        }
        Ok(true)
    }

    pub(super) fn project(
        &self, link: u64, side_index: usize, material: &ResolvedMaterialStatePoint,
    ) -> Result<Option<f64>, MaterialNetworkError> {
        let Some((binding, _)) = self.links.get(&link) else { return Ok(None); };
        let side = &binding.sides[side_index];
        let mut projected = 0.0;
        for (key, direction) in side.principal_properties.iter().zip(side.direction_in_material_axes) {
            let conductivity = material.property(key.name())
                .ok_or(MaterialNetworkError::Invalid("directional conductivity lost its original principal source"))?.value_si();
            if !conductivity.is_finite() || conductivity < 0.0 {
                return Err(MaterialNetworkError::Invalid("principal conductivity must remain finite nonnegative"));
            }
            projected = (direction * direction).mul_add(conductivity, projected);
        }
        if !projected.is_finite() || projected < 0.0 {
            return Err(MaterialNetworkError::Invalid("unrepresentable projected thermal conductivity"));
        }
        Ok(Some(projected))
    }

    pub(super) fn bindings(&self) -> impl Iterator<Item = &DirectionalConductionLink> {
        self.links.values().map(|(binding, _)| binding)
    }

    pub(super) fn retained_bindings(
        &self, pins: &BTreeMap<CellId, MaterialPropertySelection>,
        check: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<DirectionalConductionLink>, MaterialNetworkError> {
        self.links.values().map(|(binding, cells)| {
            check()?;
            let mut retained = binding.clone();
            for (cell, side) in cells.iter().zip(&mut retained.sides) {
                let Some(MaterialPropertySelection::PinnedByProperty(original)) = pins.get(cell) else {
                    return Err(MaterialNetworkError::Invalid("retained principal conductivity claims are absent"));
                };
                side.selection = MaterialPropertySelection::PinnedByProperty(original.iter()
                    .filter(|(name, _)| side.principal_properties.iter().any(|key| key.name() == name))
                    .cloned().collect());
            }
            Ok(retained)
        }).collect()
    }
}
