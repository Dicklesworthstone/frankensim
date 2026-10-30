//! Source-backed thermal symmetry, lowered into the original tensor owner.

use fs_blake3::ContentHash;
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, PropertyKey, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

use super::{SpatialMaterialSource, SpatialSolidEnthalpyError as Error, poll};

/// Thermal symmetry and source roles for one original native material entry.
/// All laws use that entry's card, temperature interval and complete fixed
/// context. This frame is independent of elastic and expansion symmetry.
#[derive(Clone, Debug)]
pub enum SpatialConductivitySource {
    /// Use `SpatialMaterialSource::conductivity` and its original selection.
    Isotropic,
    /// Three independently sourced positive principal conductivities. The
    /// scalar conductivity property is not queried for this variant.
    OrientedOrthotropic {
        principal_properties: [PropertyKey; 3],
        principal_selections: [ClaimSelection; 3],
        /// Proper rotation whose COLUMNS are the principal axes in mesh space.
        principal_to_mesh: [[f64; 3]; 3],
        /// Nonzero identity of the authored physical principal-axis assignment.
        principal_frame: ContentHash,
        /// Exact destination frame supplied to the spatial heat/solid owner.
        mesh_frame: ContentHash,
    },
}

pub(super) fn prepare(
    source: &mut SpatialMaterialSource<'_>,
    reference: &QueryPoint,
    queries: &[QueryPoint; 2],
    mut binding: SpatialConductivitySource,
    mesh_frame: ContentHash,
    cx: &Cx<'_>,
) -> Result<(ConductivityModel, SpatialConductivitySource), Error> {
    let model = match &mut binding {
        SpatialConductivitySource::Isotropic => {
            let (table, selection) = bind_table(
                source,
                &source.conductivity,
                source.conductivity_selection,
                reference,
                queries,
                cx,
            )?;
            source.conductivity_selection = selection;
            ConductivityModel::isotropic(table)
        }
        SpatialConductivitySource::OrientedOrthotropic {
            principal_properties,
            principal_selections,
            principal_to_mesh,
            principal_frame,
            mesh_frame: destination,
        } => {
            if *principal_frame == ContentHash([0; 32]) || *destination != mesh_frame {
                return Err(Error::Invalid(
                    "nonzero conductivity principal frame and exact mesh frame required",
                ));
            }
            validate_rotation(principal_to_mesh)?;
            let mut tables = Vec::with_capacity(3);
            for axis in 0..3 {
                let (table, selection) = bind_table(
                    source,
                    &principal_properties[axis],
                    principal_selections[axis],
                    reference,
                    queries,
                    cx,
                )?;
                principal_selections[axis] = selection;
                tables.push(table);
            }
            // The native tensor owner takes ROWS; the public physical frame
            // follows the same column convention as stiffness and expansion.
            let axes =
                core::array::from_fn(|i| core::array::from_fn(|j| principal_to_mesh[j][i]));
            let tables = tables.try_into().map_err(|_| {
                Error::Invalid("three principal conductivity tables required")
            })?;
            ConductivityModel::orthotropic(axes, tables)?
        }
    };
    Ok((model, binding))
}

fn bind_table(
    source: &SpatialMaterialSource<'_>,
    property: &PropertyKey,
    selection: ClaimSelection,
    reference: &QueryPoint,
    queries: &[QueryPoint; 2],
    cx: &Cx<'_>,
) -> Result<(ConductivityTable, ClaimSelection), Error> {
    poll(cx)?;
    let table = match selection {
        ClaimSelection::Policy(policy) => ConductivityTable::from_claims_at_query_points(
            source.point.material.claims(),
            property,
            "T",
            queries,
            policy,
        )?,
        ClaimSelection::Pinned(pin) => ConductivityTable::from_claims_pinned_at_query_points(
            source.point.material.claims(),
            property,
            "T",
            queries,
            pin,
        )?,
    };
    // Native tables retain all source knots and prove continuous T support.
    // Admit the selected claim's OTHER dependency axes at the original point;
    // a temperature slice cannot silently freeze an evolving coordinate.
    let pin = table
        .receipts()
        .first()
        .ok_or(Error::Invalid("source conductivity receipts required"))?
        .selected;
    let requirements = [ScalarPropertyRequirement::try_with_key(
        property,
        ScalarAdmissibility::StrictlyPositive,
    )?];
    let resolved = resolve_material_state_point(
        source.point.material,
        reference,
        &requirements,
        MaterialPropertySelection::PinnedByProperty(vec![(property.name().to_owned(), pin)]),
    )?;
    for property in resolved.properties() {
        source.point.admit_static(property)?;
    }
    Ok((table, ClaimSelection::Pinned(pin)))
}

fn validate_rotation(r: &[[f64; 3]; 3]) -> Result<(), Error> {
    if r.iter().flatten().any(|x| !x.is_finite()) {
        return Err(Error::Invalid(
            "finite conductivity principal rotation required",
        ));
    }
    for i in 0..3 {
        for j in 0..3 {
            let dot: f64 = (0..3).map(|k| r[k][i] * r[k][j]).sum();
            if (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1e-12 {
                return Err(Error::Invalid(
                    "orthonormal conductivity principal rotation required",
                ));
            }
        }
    }
    let determinant = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    if (determinant - 1.0).abs() > 1e-12 {
        return Err(Error::Invalid(
            "right-handed conductivity principal rotation required",
        ));
    }
    Ok(())
}
