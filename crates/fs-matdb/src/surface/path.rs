//! Continuous support on an affine path through the original source cells.

use super::{contains, coordinate_admitted, expanded_determinant, opposite, vertices};
use crate::MatDbError;
use fs_ivl::{Sign, orient2d};
use std::collections::BTreeMap;

/// Return every interior source-edge crossing as an affine path parameter.
/// Connectivity uses original shared vertex/edge identities, not rounded cut
/// coordinates. An admitted triangle intersects a line in one closed interval;
/// a connected chain of these intervals covering both endpoints covers the
/// entire segment. This never fills a hole using a convex hull or tolerance.
pub(crate) fn segment_knots(
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
    start: [f64; 2],
    end: [f64; 2],
) -> Result<Vec<f64>, MatDbError> {
    let refuse = |reason| MatDbError::UnsupportedEvaluation { reason };
    if start == end {
        return Ok(Vec::new());
    }
    let sides: Vec<_> = samples
        .iter()
        .map(|(point, _)| orient2d(start, end, *point))
        .collect();
    let active: Vec<_> = triangles
        .iter()
        .map(|cell| {
            cell.iter().any(|&v| sides[v as usize] != Sign::Positive)
                && cell.iter().any(|&v| sides[v as usize] != Sign::Negative)
        })
        .collect();
    let mut visited = vec![false; triangles.len()];
    let mut pending = Vec::new();
    for (index, &cell) in triangles.iter().enumerate() {
        if active[index] && contains(vertices(samples, cell), start) {
            visited[index] = true;
            pending.push(index);
        }
    }
    while let Some(index) = pending.pop() {
        for (other, cell) in triangles.iter().enumerate() {
            if visited[other] || !active[other] {
                continue;
            }
            let mut shared = [0_u32; 2];
            let mut count = 0;
            for &vertex in &triangles[index] {
                if cell.contains(&vertex) {
                    // Duplicate source cells are already inadmissible.
                    shared[count] = vertex;
                    count += 1;
                }
            }
            let touches = match count {
                1 => sides[shared[0] as usize] == Sign::Zero,
                2 => {
                    let a = sides[shared[0] as usize];
                    let b = sides[shared[1] as usize];
                    a == Sign::Zero || b == Sign::Zero || opposite(a, b)
                }
                _ => false,
            };
            if touches {
                visited[other] = true;
                pending.push(other);
            }
        }
    }
    if !triangles
        .iter()
        .enumerate()
        .any(|(index, &cell)| visited[index] && contains(vertices(samples, cell), end))
    {
        return Err(refuse(
            "declared surface cells do not continuously cover the requested affine path",
        ));
    }

    let varying = usize::from((end[1] - start[1]).abs() > (end[0] - start[0]).abs());
    let mut cuts = BTreeMap::new();
    for (index, cell) in triangles.iter().enumerate() {
        if !visited[index] {
            continue;
        }
        for edge in 0..3 {
            let i = cell[edge].min(cell[(edge + 1) % 3]);
            let j = cell[edge].max(cell[(edge + 1) % 3]);
            let a = samples[i as usize].0;
            let b = samples[j as usize].0;
            for (vertex, point) in [(i, a), (j, b)] {
                if sides[vertex as usize] == Sign::Zero
                    && start[varying].min(end[varying]) < point[varying]
                    && point[varying] < start[varying].max(end[varying])
                {
                    let t = (point[varying] - start[varying]) / (end[varying] - start[varying]);
                    if !(0.0 < t && t < 1.0) {
                        return Err(refuse(
                            "source vertex collapses onto a distinct path endpoint",
                        ));
                    }
                    cuts.insert((vertex, vertex), t);
                }
            }
            if cuts.contains_key(&(i, j))
                || !opposite(sides[i as usize], sides[j as usize])
                || !opposite(orient2d(a, b, start), orient2d(a, b, end))
            {
                // An endpoint cut is already represented by parameter 0 or 1.
                continue;
            }
            let low = expanded_determinant(a, b, start);
            let high = expanded_determinant(a, b, end);
            let t = low / (low - high);
            if !(0.0 < t && t < 1.0) {
                return Err(refuse("source edge crossing is numerically unresolved"));
            }
            let point = [0, 1].map(|axis| {
                if t <= 0.5 {
                    (end[axis] - start[axis]).mul_add(t, start[axis])
                } else {
                    (start[axis] - end[axis]).mul_add(1.0 - t, end[axis])
                }
            });
            if !point.into_iter().all(coordinate_admitted) || point == start || point == end {
                return Err(refuse(
                    "source edge crossing collapses in physical coordinates",
                ));
            }
            cuts.insert((i, j), t);
        }
    }
    let mut ordered: Vec<_> = cuts.into_values().collect();
    ordered.sort_by(f64::total_cmp);
    if ordered.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(refuse(
            "distinct source path crossings collapse to one parameter",
        ));
    }
    Ok(ordered)
}
