//! Original shell-skin, finite-cell face and contact-footprint identity.

use super::*;

pub(super) fn validate_bindings<'a, S: ThermalSectionSourceFamily<'a>, E>(
    shell: &ThermoplasticShell<'a, S>,
    thermal: &DepositionSimulation<'_>,
    contacts: &PlasticShellContacts<'_>,
    ordered: &[&ShellContactHeatBinding<'_>],
    tolerance_m: f64,
    cx: &Cx<'_>,
) -> Result<Vec<usize>, ShellContactHeatError<E>> {
    let invalid = ShellContactHeatError::Invalid;
    let witnesses = contacts.witnesses();
    let mut by_id = BTreeMap::new();
    for (i, witness) in witnesses.iter().enumerate() {
        if by_id.insert(witness.id, i).is_some() {
            return Err(invalid("original shell contact IDs must be unique"));
        }
    }
    let shell_cells: BTreeSet<_> = shell
        .accepted
        .iter()
        .flatten()
        .map(|cell| cell.id)
        .collect();
    let mut used = BTreeSet::new();
    let mut faces = BTreeSet::new();
    let mut indices = Vec::with_capacity(ordered.len());
    for binding in ordered {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        let &i = by_id
            .get(&binding.contact_id)
            .ok_or(invalid("heat route names an absent original shell contact"))?;
        let witness = &witnesses[i];
        if !used.insert(i)
            || binding.shell_side > 1
            || !positive(binding.temperature_guess_k)
            || !same_measure(witness.area_m2, binding.route.contact_area_m2())
        {
            return Err(invalid(
                "one finite-temperature heat route per original physical patch required",
            ));
        }
        let facet = shell.facets.get(witness.triangle).ok_or(invalid(
            "thermal contact refers to an absent original facet",
        ))?;
        let geometry = shell
            .mesh
            .facet(witness.triangle)
            .map_err(|_| invalid("invalid original shell contact facet"))?;
        let normal = geometry.frame[2];
        let low = facet
            .layers
            .first()
            .expect("admitted nonempty shell layers")
            .bottom_z_m;
        let high = facet
            .layers
            .last()
            .expect("admitted nonempty shell layers")
            .top_z_m;
        let at_low = near(witness.offset_m, normal.map(|n| low * n), tolerance_m);
        let at_high = near(witness.offset_m, normal.map(|n| high * n), tolerance_m);
        if at_low == at_high || !faces.insert((witness.triangle, at_high)) {
            return Err(invalid(
                "one unambiguous original heat footprint per exposed facet skin required",
            ));
        }
        let layer_index = if at_high { facet.layers.len() - 1 } else { 0 };
        let layer = &facet.layers[layer_index];
        let original = &shell.accepted[witness.triangle][layer_index];
        let ids = binding.route.cells();
        let other_id = ids[1 - binding.shell_side];
        let other = thermal
            .cells()
            .get(&other_id)
            .ok_or(invalid("inactive original obstacle thermal cell"))?;
        if ids[binding.shell_side] != original.id
            || shell_cells.contains(&other_id)
            || other.body_id == original.body_id
            || witness.area_m2 > facet.area_m2 * (1.0 + 512.0 * f64::EPSILON)
        {
            return Err(invalid(
                "heat route must join the actual outer shell layer to a distinct obstacle",
            ));
        }
        // The mechanical normal is the force direction ON the shell. On the
        // upper skin this points opposite the original outward layer normal.
        let inward = normal.map(|n| if at_high { -n } else { n });
        if !near(witness.axes_world[0], inward, 512.0 * f64::EPSILON) {
            return Err(invalid(
                "contact normal must point into its original shell skin",
            ));
        }
        let triangle = shell.mesh.tris[witness.triangle];
        let centroid: [f64; 3] = core::array::from_fn(|d| {
            triangle
                .iter()
                .map(|&n| shell.mesh.nodes[n][d])
                .sum::<f64>()
                / 3.0
        });
        let lower = core::array::from_fn(|d| centroid[d] + layer.bottom_z_m * normal[d]);
        let upper = core::array::from_fn(|d| centroid[d] + layer.top_z_m * normal[d]);
        let cell_geometry = &original.geometry;
        if !(near(cell_geometry.start_m, lower, tolerance_m)
            && near(cell_geometry.end_m, upper, tolerance_m)
            || near(cell_geometry.start_m, upper, tolerance_m)
                && near(cell_geometry.end_m, lower, tolerance_m))
        {
            return Err(invalid(
                "shell layer thermal cell must retain its original normal prism geometry",
            ));
        }
        let skin = if at_high { upper } else { lower };
        let point: [f64; 3] = core::array::from_fn(|d| {
            triangle
                .iter()
                .zip(witness.barycentric)
                .map(|(&n, w)| w * shell.mesh.nodes[n][d])
                .sum::<f64>()
                + witness.offset_m[d]
        });
        if !near(point, skin, tolerance_m)
            || !near(witness.reference_obstacle_point_m, skin, tolerance_m)
            || !(near(other.geometry.start_m, skin, tolerance_m)
                || near(other.geometry.end_m, skin, tolerance_m))
        {
            return Err(invalid(
                "mechanical witness and both original thermal faces must coincide",
            ));
        }
        indices.push(i);
    }
    if used.len() != witnesses.len() {
        return Err(invalid(
            "all original shell contact patches need their own heat route",
        ));
    }
    Ok(indices)
}

fn near(a: [f64; 3], b: [f64; 3], tolerance: f64) -> bool {
    a.iter().chain(&b).all(|v| v.is_finite())
        && a.iter()
            .zip(b)
            .map(|(&x, y)| (x - y).abs())
            .fold(0.0, f64::max)
            <= tolerance
}
