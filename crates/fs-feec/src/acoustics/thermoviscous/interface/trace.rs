//! Shared physical P2 interface: geometry, work-dual loads and active traces.
use super::*;
use std::collections::BTreeMap;

type PointKey = [u64; 3];
type FaceKey = [PointKey; 3];
type DofKey = (u8, PointKey, PointKey);

fn point_key(point: [f64; 3]) -> PointKey {
    point.map(|x| {
        if x == 0.0 {
            0.0_f64.to_bits()
        } else {
            x.to_bits()
        }
    })
}
fn face_key(model: &ThermoviscousFe<'_>, face: [u32; 3]) -> FaceKey {
    let mut key = face.map(|v| point_key(model.fe.positions()[v as usize]));
    key.sort_unstable();
    key
}
struct Surface {
    faces: BTreeMap<FaceKey, ([f64; 3], f64)>,
    dofs: BTreeMap<DofKey, usize>,
    selected: Vec<usize>,
}

fn surface(
    model: &ThermoviscousFe<'_>,
    selected_faces: &[[u32; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Surface, ThermoviscousError> {
    if selected_faces.is_empty() || selected_faces.len() > model.space.complex.faces.len() {
        return Err(invalid("nonempty original boundary interface required"));
    }
    let mut faces = BTreeMap::new();
    let mut dofs = BTreeMap::new();
    let mut insert = |key, global| -> Result<(), ThermoviscousError> {
        if dofs
            .insert(key, global)
            .is_some_and(|prior| prior != global)
        {
            return Err(invalid(
                "coincident distinct original interface coordinates refused",
            ));
        }
        Ok(())
    };
    for &face in selected_faces {
        poll(cancelled)?;
        let (_, normal, area) = model.fe.boundary_geometry(face)?;
        if model.wall_faces.contains(&canonical(face))
            || faces
                .insert(face_key(model, face), (normal, area))
                .is_some()
        {
            return Err(invalid("unique original free interface faces required"));
        }
        for &v in &face {
            let point = point_key(model.fe.positions()[v as usize]);
            insert((0, point, point), v as usize)?;
        }
        for a in 0..3 {
            for b in a + 1..3 {
                let mut edge = [face[a], face[b]];
                edge.sort_unstable();
                let id = model
                    .space
                    .complex
                    .edges
                    .binary_search(&edge)
                    .map_err(|_| invalid("original P2 interface edge required"))?;
                let mut points = edge.map(|v| point_key(model.fe.positions()[v as usize]));
                points.sort_unstable();
                insert((1, points[0], points[1]), model.space.edge_off + id)?;
            }
        }
    }
    let mut selected = Vec::new();
    for (mode, traces) in model.physical_modes.iter().enumerate() {
        let overlaps = traces
            .iter()
            .filter(|trace| faces.contains_key(&trace.vertices))
            .count();
        if overlaps > 0 {
            if overlaps != traces.len() {
                return Err(invalid(
                    "one pressure mode cannot straddle internal and external faces",
                ));
            }
            selected.push(mode);
        }
    }
    if selected.is_empty() {
        return Err(invalid(
            "original interface pressure-face ownership required",
        ));
    }
    Ok(Surface {
        faces,
        dofs,
        selected,
    })
}

/// Original weak loads, before any harmonic solve or transient mass inverse.
pub(in crate::acoustics::thermoviscous) struct TraceBinding {
    pub layout: ConformingThermoviscousInterface,
    pub selected: [Vec<usize>; 2],
    /// Row-major physical weak loads. Left positive, right negative.
    pub loads: [Vec<f64>; 2],
    pub heat_flux_areas_m2: Vec<f64>,
}

pub(in crate::acoustics::thermoviscous) fn prepare(
    models: [&ThermoviscousFe<'_>; 2],
    faces: [&[[u32; 3]]; 2],
    maximum_trace_dofs: usize,
    maximum_dense_bytes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TraceBinding, ThermoviscousError> {
    poll(cancelled)?;
    let surfaces = [
        surface(models[0], faces[0], cancelled)?,
        surface(models[1], faces[1], cancelled)?,
    ];
    if surfaces[0].faces.keys().ne(surfaces[1].faces.keys())
        || surfaces[0].dofs.keys().ne(surfaces[1].dofs.keys())
    {
        return Err(invalid(
            "complete original P2 interface triangulations must match",
        ));
    }
    for (key, &(normal, area)) in &surfaces[0].faces {
        poll(cancelled)?;
        let (other, other_area) = surfaces[1].faces[key];
        let dot: f64 = normal.iter().zip(other).map(|(a, b)| a * b).sum();
        if (dot + 1.0).abs() > 1e-12 || (area - other_area).abs() > 1e-12 * area {
            return Err(invalid(
                "opposite interface normals and equal original areas required",
            ));
        }
    }
    let trace_dofs: [Vec<_>; 2] =
        std::array::from_fn(|side| surfaces[side].dofs.values().copied().collect());
    let m = trace_dofs[0].len();
    for side in 0..2 {
        if let Some(wall) = &models[side].moving_wall {
            for &dof in &trace_dofs[side] {
                if wall
                    .velocity_shapes
                    .iter()
                    .any(|shape| shape[dof].iter().any(|x| *x != 0.0))
                {
                    return Err(invalid(
                        "moving structural and fluid interface traces require one kinematic owner",
                    ));
                }
            }
        }
    }
    let traction_indices: Vec<_> = (0..m)
        .filter(|&i| (0..2).any(|side| models[side].velocity_map[trace_dofs[side][i]].is_some()))
        .collect();
    let heat_flux_indices: Vec<_> = (0..m)
        .filter(|&i| (0..2).any(|side| models[side].temperature_map[trace_dofs[side][i]].is_some()))
        .collect();
    let mechanical = 3 * traction_indices.len();
    let count = mechanical + heat_flux_indices.len();
    if count == 0 {
        return Err(invalid(
            "nonzero complete original conforming trace rank required",
        ));
    }
    if count > maximum_trace_dofs {
        return Err(ThermoviscousError::Budget {
            what: "conforming interface coordinates",
            requested: count,
            limit: maximum_trace_dofs,
        });
    }
    let n = models[0].n + models[1].n;
    let bytes = m
        .checked_mul(m)
        .and_then(|v| v.checked_mul(16))
        .and_then(|v| {
            n.checked_mul(count)
                .and_then(|k| k.checked_mul(8))
                .and_then(|k| v.checked_add(k))
        })
        .ok_or(ThermoviscousError::NonFinite)?;
    if bytes > maximum_dense_bytes {
        return Err(ThermoviscousError::Budget {
            what: "conforming physical trace bytes",
            requested: bytes,
            limit: maximum_dense_bytes,
        });
    }
    // Compare complete original surface masses before sharing the left one.
    // No lumping, fitted area, reduced trace or imposed pressure continuity.
    let mut mass = [vec![0.0; m * m], vec![0.0; m * m]];
    for side in 0..2 {
        let source = models[side].fe.boundary_mass(faces[side])?;
        let indices: BTreeMap<_, _> = trace_dofs[side]
            .iter()
            .enumerate()
            .map(|(i, &g)| (g, i))
            .collect();
        for (i, &global) in trace_dofs[side].iter().enumerate() {
            poll(cancelled)?;
            let (columns, values) = source.row(global);
            for (&column, &value) in columns.iter().zip(values) {
                if let Some(&j) = indices.get(&column) {
                    mass[side][i * m + j] = value;
                }
            }
        }
    }
    for (a, b) in mass[0].iter().zip(&mass[1]) {
        if !a.is_finite()
            || !b.is_finite()
            || (a - b).abs() > 512.0 * f64::EPSILON * (a.abs() + b.abs())
        {
            return Err(invalid(
                "original P2 physical surface mass must agree without truncation",
            ));
        }
    }
    let mut loads: [Vec<f64>; 2] = std::array::from_fn(|side| vec![0.0; models[side].n * count]);
    for side in 0..2 {
        let model = models[side];
        let sign = if side == 0 { 1.0 } else { -1.0 };
        for (i, &global) in trace_dofs[side].iter().enumerate() {
            poll(cancelled)?;
            if let Some(v) = model.velocity_map[global] {
                for c in 0..3 {
                    let row = c * model.velocity_count + v;
                    for (j, &trace) in traction_indices.iter().enumerate() {
                        loads[side][row * count + c * traction_indices.len() + j] =
                            sign * mass[0][i * m + trace];
                    }
                }
            }
            if let Some(t) = model.temperature_map[global] {
                let row = model.temperature_start + t;
                for (j, &trace) in heat_flux_indices.iter().enumerate() {
                    // q [W/m2] is dual to theta=T'/T0, with no extra T0.
                    loads[side][row * count + mechanical + j] = sign * mass[0][i * m + trace];
                }
            }
        }
    }
    let (areas, _) = models[0].fe.boundary_load(faces[0])?;
    let heat_flux_areas_m2 = heat_flux_indices
        .iter()
        .map(|&i| areas[trace_dofs[0][i]])
        .collect();
    Ok(TraceBinding {
        layout: ConformingThermoviscousInterface {
            trace_dofs,
            traction_indices,
            heat_flux_indices,
        },
        selected: surfaces.map(|surface| surface.selected),
        loads,
        heat_flux_areas_m2,
    })
}
