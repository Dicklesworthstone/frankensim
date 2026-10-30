//! Source principal conductivities expressed in the existing Cartesian heat
//! stencil. A general rotation needs cross-gradient fluxes that stencil does
//! not own, so only an explicitly declared signed permutation is admitted.

use fs_lbm::conjugate::ChtError;
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::state_point::{
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use super::super::{J, array, fields, number, quote, values};
use super::invalid;

pub(in crate::electrothermal_command) struct Conductivity {
    requirements: Vec<ScalarPropertyRequirement>,
    /// Principal axes are rows, expressed in grid coordinates. None is the
    /// existing isotropic declaration, never inferred from a material name.
    axes: Option<[[f64; 3]; 3]>,
    principal_to_grid: [usize; 3],
}

impl Conductivity {
    pub(super) fn load(
        pack: &NormalizedMaterialCardPack,
        input: Option<&J>,
    ) -> Result<Self, ChtError> {
        let Some(input) = input else {
            return Ok(Self {
                requirements: vec![requirement(pack, THERMAL_CONDUCTIVITY_PROPERTY)?],
                axes: None,
                principal_to_grid: [0, 1, 2],
            });
        };
        let parse = || -> Result<(Vec<&str>, [[f64; 3]; 3], [usize; 3]), String> {
            fields(input, &["principal_properties", "principal_axes"])?;
            let properties = array(input, "principal_properties", 3)?;
            let rows = array(input, "principal_axes", 3)?;
            if properties.len() != 3 || rows.len() != 3 {
                return Err(
                    "declare three principal property names and three principal axes".into(),
                );
            }
            let names = properties
                .iter()
                .map(|p| {
                    p.as_str().filter(|s| !s.trim().is_empty()).ok_or_else(|| {
                        "principal property names must be nonblank strings".to_string()
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            if names[0] == names[1] || names[0] == names[2] || names[1] == names[2] {
                return Err("principal properties must name three distinct source roles".into());
            }
            let mut axes = [[0.0; 3]; 3];
            let mut mapping = [0; 3];
            let mut used = [false; 3];
            for (principal, row) in rows.iter().enumerate() {
                let entries = row
                    .as_array()
                    .filter(|r| r.len() == 3)
                    .ok_or("each principal axis needs three grid components")?;
                let mut aligned = None;
                for (grid, entry) in entries.iter().enumerate() {
                    let value = number(entry)?;
                    axes[principal][grid] = value;
                    if value == 1.0 || value == -1.0 {
                        if aligned.replace(grid).is_some() {
                            return Err("principal axes must be signed Cartesian unit axes".into());
                        }
                    } else if value != 0.0 {
                        return Err("off-grid principal axes need a full tensor flux operator; rotation is not projected onto the diagonal".into());
                    }
                }
                let grid = aligned.ok_or("principal axes must be signed Cartesian unit axes")?;
                if used[grid] {
                    return Err(
                        "principal axes must form a signed permutation of the grid axes".into(),
                    );
                }
                used[grid] = true;
                mapping[principal] = grid;
            }
            Ok((names, axes, mapping))
        };
        let (names, axes, principal_to_grid) =
            parse().map_err(|e| invalid("material.thermal_conductivity", e))?;
        let requirements = names
            .into_iter()
            .map(|name| requirement(pack, name))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            requirements,
            axes: Some(axes),
            principal_to_grid,
        })
    }

    pub(super) fn requirements(&self) -> &[ScalarPropertyRequirement] {
        &self.requirements
    }

    pub(super) fn isotropic_requirement(&self) -> Result<&ScalarPropertyRequirement, ChtError> {
        if self.axes.is_some() {
            return Err(invalid(
                "fluid.thermal_conductivity",
                "the native fluid reference requires isotropic conductivity; directional solid materials remain supported",
            ));
        }
        Ok(&self.requirements[0])
    }

    pub(super) fn grid_values(&self, state: &ResolvedMaterialStatePoint) -> [f64; 3] {
        let value = |i: usize| {
            state
                .property(self.requirements[i].name())
                .expect("required source conductivity was resolved")
                .value_si()
        };
        if self.axes.is_none() {
            return [value(0); 3];
        }
        let mut grid = [0.0; 3];
        for principal in 0..3 {
            grid[self.principal_to_grid[principal]] = value(principal);
        }
        grid
    }

    pub(in crate::electrothermal_command) fn render(&self) -> String {
        match self.axes {
            None => format!(
                "{{\"kind\":\"isotropic\",\"property\":{}}}",
                quote(self.requirements[0].name())
            ),
            Some(axes) => format!(
                "{{\"kind\":\"grid-aligned-orthotropic\",\"principal_properties\":[{}],\"principal_axes\":[{}]}}",
                self.requirements
                    .iter()
                    .map(|r| quote(r.name()))
                    .collect::<Vec<_>>()
                    .join(","),
                axes.into_iter().map(values).collect::<Vec<_>>().join(","),
            ),
        }
    }
}

fn requirement(
    pack: &NormalizedMaterialCardPack,
    name: &str,
) -> Result<ScalarPropertyRequirement, ChtError> {
    let thermal_kind = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let candidates = pack.card().claims_for(name);
    let (_, first) = candidates.first().ok_or_else(|| {
        invalid(
            "material.thermal_conductivity",
            format!("card has no conductivity claim for {name}"),
        )
    })?;
    let key = &first.key;
    if (key.quantity() != thermal_kind
        && key.quantity() != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS))
        || candidates.iter().any(|(_, claim)| claim.key != *key)
    {
        return Err(invalid(
            "material.thermal_conductivity",
            "unsupported or conflicting conductivity quantity schemas",
        ));
    }
    ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
        .map_err(|e| invalid("material.thermal_conductivity", e))
}
