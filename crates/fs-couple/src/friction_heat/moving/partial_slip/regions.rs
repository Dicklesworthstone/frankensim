//! Exact endpoint intersections of transported memory, both original caloric
//! grids and the supplied pressure field. Monotone cursors avoid a Cartesian
//! product of the two surfaces; a native strain is copied, never interpolated.

use super::*;

fn shifted(span: [f64; 2], offset: f64) -> Result<[f64; 2], DepositionError> {
    Ok([finite(span[0] + offset)?, finite(span[1] + offset)?])
}

pub(super) fn partition<E>(
    contact: &MovingPartialSlip,
    profile: Option<&MovingPressureProfile>,
    transported: &[MaterialSlipParcel],
    end: MovingContactPosition,
    patches: &[SweptContactPatch],
    maximum_regions: usize,
    checkpoint: &impl Fn() -> Result<(), DepositionError>,
) -> Result<Vec<MovingSlipRegion>, MovingPartialSlipError<E>> {
    let route = &contact.route;
    let length = route.footprint_m[1] - route.footprint_m[0];
    let area = route.footprint_area_m2();
    let footprint_shift = finite(end.footprint_offset_m - end.body_offsets_m[0])?;
    let b_shift = finite(end.body_offsets_m[1] - end.body_offsets_m[0])?;
    let mut pressure = Vec::new();
    if let Some(profile) = profile {
        for cell in profile.cells() {
            checkpoint()?;
            pressure.push((
                shifted(cell.span_m, footprint_shift)?,
                cell.normal_pressure_pa,
                profile.normalized_patch(cell.normal_pressure_pa, area)?,
            ));
        }
    } else {
        pressure.push((
            shifted(route.footprint_m, footprint_shift)?,
            finite(contact.reference.normal_patch().normal_load_n() / area)?,
            contact.reference.normal_patch().clone(),
        ));
    }
    let mut interfaces = BTreeMap::new();
    for (index, patch) in patches.iter().enumerate() {
        checkpoint()?;
        interfaces.insert(patch.cells, index);
    }
    let mut regions = Vec::new();
    let (mut a_index, mut b_index, mut p_index) = (0, 0, 0);
    for parcel in transported {
        let mut left = parcel.span_m[0];
        while left < parcel.span_m[1] {
            checkpoint()?;
            while a_index < route.surfaces[0].len()
                && route.surfaces[0][a_index].face.span_m[1] <= left
            {
                checkpoint()?;
                a_index += 1;
            }
            while b_index < route.surfaces[1].len()
                && finite(route.surfaces[1][b_index].face.span_m[1] + b_shift)? <= left
            {
                checkpoint()?;
                b_index += 1;
            }
            while p_index < pressure.len() && pressure[p_index].0[1] <= left {
                checkpoint()?;
                p_index += 1;
            }
            let missing = || {
                MovingPartialSlipError::Invalid(
                    "both original surface grids and pressure must cover each material region",
                )
            };
            let a = route.surfaces[0].get(a_index).ok_or_else(missing)?;
            let b = route.surfaces[1].get(b_index).ok_or_else(missing)?;
            let p = pressure.get(p_index).ok_or_else(missing)?;
            let b_span = shifted(b.face.span_m, b_shift)?;
            if a.face.span_m[0] > left || b_span[0] > left || p.0[0] > left {
                return Err(missing());
            }
            let right = parcel.span_m[1]
                .min(a.face.span_m[1])
                .min(b_span[1])
                .min(p.0[1]);
            if right <= left {
                return Err(missing());
            }
            if regions.len() >= maximum_regions {
                return Err(MovingPartialSlipError::RegionBudget);
            }
            let pair = [a.face.cell, b.face.cell];
            let endpoint_interface_index =
                *interfaces
                    .get(&pair)
                    .ok_or(MovingPartialSlipError::Invalid(
                        "endpoint material pair must retain its original swept thermal interface",
                    ))?;
            let area_fraction = finite((right - left) / length)?;
            regions.push(MovingSlipRegion {
                material_span_m: [left, right],
                surface_a_cell: pair[0],
                surface_b_cell: pair[1],
                endpoint_interface_index,
                area_fraction,
                normal_pressure_pa: p.1,
                normal_load_n: finite(area_fraction * p.2.normal_load_n())?,
                state: parcel.state.clone(),
                normalized_patch: p.2.clone(),
            });
            left = right;
        }
    }
    Ok(regions)
}
