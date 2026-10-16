//! Exact-sign topology; f64 only chooses metric refinement points.
use super::*;
use fs_ivl::incircle;

pub(super) fn edges(tris: &[[usize; 3]]) -> BTreeMap<(usize, usize), Vec<(usize, usize)>> {
    let mut out = BTreeMap::new();
    for (t, tri) in tris.iter().enumerate() {
        for i in 0..3 {
            out.entry(edge(tri[i], tri[(i + 1) % 3]))
                .or_insert_with(Vec::new)
                .push((t, tri[(i + 2) % 3]));
        }
    }
    out
}
fn ccw(p: &[[f64; 2]], mut t: [usize; 3]) -> Result<[usize; 3], PlanarError> {
    match orient2d(p[t[0]], p[t[1]], p[t[2]]) {
        Sign::Zero => Err(PlanarError::Precision("zero-area triangulation cell")),
        Sign::Negative => {
            t.swap(1, 2);
            Ok(t)
        }
        Sign::Positive => Ok(t),
    }
}
fn contains(p: &[[f64; 2]], t: [usize; 3], q: [f64; 2]) -> bool {
    (0..3).all(|i| orient2d(p[t[i]], p[t[(i + 1) % 3]], q) != Sign::Negative)
}
fn ear_clip(p: &[[f64; 2]], mut polygon: Vec<usize>) -> Result<Vec<[usize; 3]>, PlanarError> {
    if polygon.len() < 3 {
        return Ok(Vec::new());
    }
    if area(&polygon.iter().map(|&i| p[i]).collect::<Vec<_>>()) < 0.0 {
        polygon.reverse();
    }
    let mut out = Vec::new();
    while polygon.len() > 3 {
        let n = polygon.len();
        let ear = (0..n)
            .find(|&i| {
                let t = [polygon[(i + n - 1) % n], polygon[i], polygon[(i + 1) % n]];
                orient2d(p[t[0]], p[t[1]], p[t[2]]) == Sign::Positive
                    && !polygon
                        .iter()
                        .any(|&v| !t.contains(&v) && contains(p, t, p[v]))
            })
            .ok_or(PlanarError::Topology(
                "constraint cavity has no exact convex ear",
            ))?;
        out.push([
            polygon[(ear + n - 1) % n],
            polygon[ear],
            polygon[(ear + 1) % n],
        ]);
        polygon.remove(ear);
    }
    out.push(ccw(p, [polygon[0], polygon[1], polygon[2]])?);
    Ok(out)
}

pub(super) fn initial(
    input: &[[f64; 2]],
    segments: &[Segment],
    rings: &[Vec<[f64; 2]>],
    cx: &Cx<'_>,
) -> Result<Vec<[usize; 3]>, PlanarError> {
    let mut p = input.to_vec();
    let mut low = input[0];
    let mut high = low;
    for point in input {
        for k in 0..2 {
            low[k] = low[k].min(point[k]);
            high[k] = high[k].max(point[k]);
        }
    }
    let centre = [
        f64::midpoint(low[0], high[0]),
        f64::midpoint(low[1], high[1]),
    ];
    let size = (high[0] - low[0]).max(high[1] - low[1]);
    if !size.is_finite() || size <= 0.0 {
        return Err(PlanarError::InvalidInput("zero geometric scale"));
    }
    let n = input.len();
    p.extend([
        [centre[0] - 64.0 * size, centre[1] - 32.0 * size],
        [centre[0] + 64.0 * size, centre[1] - 32.0 * size],
        [centre[0], centre[1] + 64.0 * size],
    ]);
    let mut tris = vec![[n, n + 1, n + 2]];
    for id in 0..n {
        checkpoint(cx)?;
        let cavity: BTreeSet<_> = tris
            .iter()
            .enumerate()
            .filter(|(_, t)| incircle(p[t[0]], p[t[1]], p[t[2]], p[id]) == Sign::Positive)
            .map(|(i, _)| i)
            .collect();
        if cavity.is_empty() {
            return Err(PlanarError::Precision(
                "Delaunay insertion has an empty cavity",
            ));
        }
        let mut boundary = BTreeMap::<(usize, usize), (usize, (usize, usize))>::new();
        for &i in &cavity {
            for k in 0..3 {
                let a = tris[i][k];
                let b = tris[i][(k + 1) % 3];
                let e = boundary.entry(edge(a, b)).or_insert((0, (a, b)));
                e.0 += 1;
            }
        }
        tris = tris
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !cavity.contains(i))
            .map(|(_, t)| t)
            .collect();
        for (_, (count, (a, b))) in boundary {
            if count == 1 && orient2d(p[a], p[b], p[id]) != Sign::Zero {
                tris.push(ccw(&p, [a, b, id])?);
            }
        }
    }
    tris.retain(|t| t.iter().all(|&i| i < n));
    // Recover a segment by deleting its crossed triangle strip and ear-clipping
    // the two resulting polygons. No flip heuristic can strand a nonconvex
    // strip; already recovered, noncrossing constraints remain on its boundary.
    for s in segments {
        checkpoint(cx)?;
        let adj = edges(&tris);
        if adj.contains_key(&edge(s.a, s.b)) {
            continue;
        }
        let mut cavity = BTreeSet::new();
        for (&(a, b), incident) in &adj {
            if opposite(
                orient2d(p[s.a], p[s.b], p[a]),
                orient2d(p[s.a], p[s.b], p[b]),
            ) && opposite(orient2d(p[a], p[b], p[s.a]), orient2d(p[a], p[b], p[s.b]))
            {
                cavity.extend(incident.iter().map(|&(t, _)| t));
            }
        }
        if cavity.is_empty() {
            return Err(PlanarError::Topology("constraint strip is empty"));
        }
        let selected: Vec<_> = cavity.iter().map(|&i| tris[i]).collect();
        let strip_edges = edges(&selected);
        let mut graph = BTreeMap::<usize, Vec<usize>>::new();
        for ((a, b), inc) in strip_edges {
            if inc.len() == 1 {
                graph.entry(a).or_default().push(b);
                graph.entry(b).or_default().push(a);
            }
        }
        if graph.values().any(|v| v.len() != 2)
            || !graph.contains_key(&s.a)
            || !graph.contains_key(&s.b)
        {
            return Err(PlanarError::Topology(
                "constraint strip boundary is not a simple cycle",
            ));
        }
        let mut cycle = vec![s.a];
        let mut previous = s.a;
        let mut current = graph[&s.a][0];
        while current != s.a {
            cycle.push(current);
            if cycle.len() > graph.len() {
                return Err(PlanarError::Topology("constraint strip boundary cycles"));
            }
            let next = *graph[&current].iter().find(|&&i| i != previous).unwrap();
            previous = current;
            current = next;
        }
        let at = cycle
            .iter()
            .position(|&i| i == s.b)
            .ok_or(PlanarError::Topology("constraint endpoint missing"))?;
        let first = cycle[..=at].to_vec();
        let second = std::iter::once(s.a)
            .chain(cycle[at..].iter().rev().copied())
            .collect();
        let mut replacement = ear_clip(&p, first)?;
        replacement.extend(ear_clip(&p, second)?);
        tris = tris
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !cavity.contains(i))
            .map(|(_, t)| t)
            .collect();
        tris.extend(replacement);
    }
    let protected: BTreeSet<_> = segments.iter().map(|s| edge(s.a, s.b)).collect();
    legalize(&p, &mut tris, &protected, cx)?;
    // Seed each exterior/hole component by its winding classification, then
    // flood ONLY through unconstrained edges. This retains concave boundaries.
    let adj = edges(&tris);
    let mut removed = vec![false; tris.len()];
    let mut stack = Vec::new();
    for (i, t) in tris.iter().enumerate() {
        let q = [
            (p[t[0]][0] + p[t[1]][0] + p[t[2]][0]) / 3.0,
            (p[t[0]][1] + p[t[1]][1] + p[t[2]][1]) / 3.0,
        ];
        if !inside(q, &rings[0]) || rings.iter().skip(1).any(|r| inside(q, r)) {
            removed[i] = true;
            stack.push(i);
        }
    }
    while let Some(i) = stack.pop() {
        for k in 0..3 {
            let e = edge(tris[i][k], tris[i][(k + 1) % 3]);
            if !protected.contains(&e) {
                for &(j, _) in &adj[&e] {
                    if !removed[j] {
                        removed[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
    }
    tris = tris
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed[*i])
        .map(|(_, t)| t)
        .collect();
    if tris.is_empty() {
        return Err(PlanarError::Topology("empty retained domain"));
    }
    Ok(tris)
}

pub(super) fn insert(
    p: &[[f64; 2]],
    tris: &mut Vec<[usize; 3]>,
    id: usize,
    _protected: &BTreeSet<(usize, usize)>,
    cx: &Cx<'_>,
) -> Result<(), PlanarError> {
    checkpoint(cx)?;
    let t = tris
        .iter()
        .position(|&t| contains(p, t, p[id]))
        .ok_or(PlanarError::Precision(
            "refinement point is outside material",
        ))?;
    let original = tris[t];
    let on =
        (0..3).find(|&i| orient2d(p[original[i]], p[original[(i + 1) % 3]], p[id]) == Sign::Zero);
    if let Some(i) = on {
        let (a, b) = (original[i], original[(i + 1) % 3]);
        let incident = edges(tris)[&edge(a, b)].clone();
        let indices: BTreeSet<_> = incident.iter().map(|&(t, _)| t).collect();
        let mut replacement = Vec::new();
        for (_, c) in incident {
            replacement.push(ccw(p, [a, id, c])?);
            replacement.push(ccw(p, [id, b, c])?);
        }
        *tris = tris
            .iter()
            .enumerate()
            .filter(|(i, _)| !indices.contains(i))
            .map(|(_, t)| *t)
            .collect();
        tris.extend(replacement);
    } else {
        tris[t] = ccw(p, [original[0], original[1], id])?;
        tris.push(ccw(p, [original[1], original[2], id])?);
        tris.push(ccw(p, [original[2], original[0], id])?);
    }
    Ok(())
}

pub(super) fn legalize(
    p: &[[f64; 2]],
    tris: &mut [[usize; 3]],
    protected: &BTreeSet<(usize, usize)>,
    cx: &Cx<'_>,
) -> Result<(), PlanarError> {
    let limit = tris
        .len()
        .saturating_mul(tris.len())
        .saturating_mul(8)
        .max(1);
    for _ in 0..limit {
        checkpoint(cx)?;
        let mut flipped = false;
        for ((a, b), incident) in edges(tris) {
            if incident.len() != 2 || protected.contains(&edge(a, b)) {
                continue;
            }
            let (i, c) = incident[0];
            let (j, d) = incident[1];
            if !opposite(orient2d(p[c], p[d], p[a]), orient2d(p[c], p[d], p[b])) {
                continue;
            }
            let t = ccw(p, [a, b, c])?;
            if incircle(p[t[0]], p[t[1]], p[t[2]], p[d]) == Sign::Positive {
                tris[i] = ccw(p, [c, d, a])?;
                tris[j] = ccw(p, [d, c, b])?;
                flipped = true;
                break;
            }
        }
        // Exact cocircular Zero is never flipped, preventing tie cycles.
        if !flipped {
            return Ok(());
        }
    }
    Err(PlanarError::Topology(
        "Delaunay legalization exceeded its finite flip bound",
    ))
}
pub(super) fn circumcentre(p: [[f64; 2]; 3]) -> Result<[f64; 2], PlanarError> {
    let a = sub(p[1], p[0]);
    let b = sub(p[2], p[0]);
    let d = 2.0 * (a[0] * b[1] - a[1] * b[0]);
    if !d.is_finite() || d == 0.0 {
        return Err(PlanarError::Precision("unrepresentable circumcentre"));
    }
    let aa = dot(a, a);
    let bb = dot(b, b);
    Ok([
        p[0][0] + (aa * b[1] - bb * a[1]) / d,
        p[0][1] + (a[0] * bb - b[0] * aa) / d,
    ])
}
