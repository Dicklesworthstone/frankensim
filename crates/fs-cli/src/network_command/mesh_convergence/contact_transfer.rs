//! Transfer a contact's physical support, not a cached integration stencil.
//! Nonmatching traces refine independently; their common refinement is rebuilt
//! by ThermalInterfaces on the admitted child mesh. Do not weld the sides,
//! force one-to-one pairs, rescale R'', or enlarge geometry tolerances/budgets.
use super::*;

pub(super) fn refine(
    cx: &Cx<'_>, contact: &mut J, split: &split::Split, old_vertices: usize,
) -> Result<()> {
    match (contact.get("face_pairs"), contact.get("nonmatching")) {
        (Some(_), None) => refine_matching(cx, contact, split, old_vertices),
        (None, Some(_)) => refine_nonmatching(cx, contact, split, old_vertices),
        _ => Err(bad("refined contact requires exactly one matching or nonmatching declaration")),
    }
}

fn refine_matching(
    cx: &Cx<'_>, contact: &mut J, split: &split::Split, old_vertices: usize,
) -> Result<()> {
    // Keep the original matching path and its deterministic pairing order.
    let mut children = Vec::new();
    for pair in array(get(contact, "face_pairs")?, "contact face pairs", 200_000)? {
        poll(cx)?;
        let a = indices::<3>(get(pair, "side_a")?, "side_a", old_vertices)?;
        let b = indices::<3>(get(pair, "side_b")?, "side_b", old_vertices)?;
        for (a, b) in split.contact_children(a, b)? {
            children.push(J::Object(vec![
                ("side_a".into(), jindices(&a)), ("side_b".into(), jindices(&b)),
            ]));
        }
    }
    replace(contact, "face_pairs", J::Array(children))
}

fn refine_nonmatching(
    cx: &Cx<'_>, contact: &mut J, split: &split::Split, old_vertices: usize,
) -> Result<()> {
    let declaration = get(contact, "nonmatching")?;
    let a = side_children(cx, declaration, "side_a_faces", split, old_vertices)?;
    let b = side_children(cx, declaration, "side_b_faces", split, old_vertices)?;
    // Match the actual binder's budget: cross-side tests PLUS same-side
    // disjointness tests. The global confirmation consumes this budget too.
    let faces = a.len().checked_add(b.len())
        .ok_or_else(|| exhausted("refined contact face count overflow"))?;
    let tests = faces.checked_mul(faces.saturating_sub(1)).map(|n| n / 2)
        .ok_or_else(|| exhausted("refined contact pair-test count overflow"))?;
    let cap = count(get(declaration, "max_pair_tests")?, "contact max_pair_tests", 1_000_000)?;
    if tests > cap {
        return Err(exhausted(format!(
            "refined nonmatching contact needs {tests} pair tests for {} + {} faces, exceeding its unchanged max_pair_tests={cap}; no convergence result published",
            a.len(), b.len(),
        )));
    }
    poll(cx)?;
    let declaration = member_mut(contact, "nonmatching")?;
    replace(declaration, "side_a_faces", J::Array(a))?;
    replace(declaration, "side_b_faces", J::Array(b))?;
    // All remaining fields (including R'', plane/coverage tolerances, overlap
    // output limit and source provenance) remain verbatim in the cloned input.
    Ok(())
}

fn side_children(
    cx: &Cx<'_>, declaration: &J, key: &str, split: &split::Split, old_vertices: usize,
) -> Result<Vec<J>> {
    let mut children = BTreeSet::new();
    for face in array(get(declaration, key)?, key, 200_000)? {
        poll(cx)?;
        let face = indices::<3>(face, key, old_vertices)?;
        for mut child in split.face_children(face)? {
            child.sort_unstable();
            if !children.insert(child) {
                return Err(producer("refinement assigned a contact child face twice"));
            }
        }
    }
    if children.is_empty() {
        return Err(producer("refinement lost an entire nonmatching contact side"));
    }
    Ok(children.iter().map(jindices).collect())
}
