//! Registered skins: exact one-to-one XY-neighborhood correspondence, not a
//! nearest-point weld. Ambiguity refuses even when greedy matching could hide it.
use super::{Spec, ObjDocument, MAX_NODES, chart};
use std::collections::{BTreeMap, BTreeSet};

const MAX_CANDIDATE_CHECKS: usize = MAX_NODES * 64;
fn cell(p: [f64; 3], tolerance: f64) -> Result<(i64, i64), String> {
    let q = [(p[0]/tolerance).floor(), (p[1]/tolerance).floor()];
    // Stay well inside the exact-integer region and leave room for neighbors.
    if q.iter().any(|v| !v.is_finite() || v.abs() > (1u64 << 50) as f64) {
        return Err("pairing grid unresolved at this frame origin and tolerance".into());
    }
    Ok((q[0] as i64, q[1] as i64))
}
pub(super) fn projected(doc: &ObjDocument, spec: &Spec, upper: &BTreeSet<usize>,
    lower: &BTreeSet<usize>, tolerance: f64) -> Result<BTreeMap<usize, usize>, String> {
    if upper.len() != lower.len() || !upper.is_disjoint(lower) {
        return Err("projected pairing requires equally sized disjoint skin vertex sets".into());
    }
    let mut grid = BTreeMap::<(i64, i64), Vec<(usize, [f64; 3])>>::new();
    for &id in lower {
        let p = chart(doc, id, spec)?;
        grid.entry(cell(p, tolerance)?).or_default().push((id, p));
    }
    let mut pairs = BTreeMap::new();
    let mut used = BTreeSet::new();
    let mut checks = 0usize;
    for &id in upper {
        let p = chart(doc, id, spec)?;
        let (x, y) = cell(p, tolerance)?;
        let mut found = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(candidates) = grid.get(&(x+dx, y+dy)) {
                    for &(other, q) in candidates {
                        checks += 1;
                        if checks > MAX_CANDIDATE_CHECKS {
                            return Err("projected pairing candidate budget exceeded; refine the selection or supply explicit pairs".into());
                        }
                        let distance2 = (p[0]-q[0]).powi(2)+(p[1]-q[1]).powi(2);
                        if distance2 <= tolerance*tolerance {
                            if found.replace(other).is_some() {
                                return Err(format!("ambiguous projected correspondence at upper vertex {id}; no nearest-neighbor guessing"));
                            }
                        }
                    }
                }
            }
        }
        let other = found.ok_or_else(|| format!("no lower vertex within the declared projected tolerance of upper vertex {id}"))?;
        if !used.insert(other) { return Err("projected correspondence is not one-to-one".into()); }
        pairs.insert(id, other);
    }
    // Equal cardinality and uniqueness imply full coverage, without node welds.
    Ok(pairs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{import, tests::{fixture, SPEC, PAIRS}};
    use std::fmt::Write;

    const AUTO: &str = "frankensim-board-skins-v1\nlower,lower\nthickness,geometry\nthickness-range,0.001,0.03\npairing,projected,1e-7\n";
    #[test]
    fn registered_skins_match_explicit_pairs_without_requiring_source_vertex_order() {
        let obj = fixture(0.015, 1.);
        let original = import(&obj, SPEC, PAIRS).unwrap();
        assert_eq!(original.fsb, import(&obj, SPEC, AUTO).unwrap().fsb);
        let rows: Vec<_> = obj.lines().collect();
        // Reorder only lower vertices, retaining their actual coordinates/faces.
        let order = [7usize, 5, 9, 6, 8];
        let mut remap: BTreeMap<_, _> = (1..=5).map(|i| (i,i)).collect();
        let mut changed = rows[..5].join("\n")+"\n";
        for (i, &old) in order.iter().enumerate() {
            writeln!(changed,"{}",rows[old]).unwrap(); remap.insert(old+1,i+6);
        }
        for row in &rows[10..] {
            if let Some(face) = row.strip_prefix("f ") {
                let t: Vec<_> = face.split_whitespace().map(|s| remap[&s.parse::<usize>().unwrap()]).collect();
                writeln!(changed,"f {} {} {}",t[0],t[1],t[2]).unwrap();
            } else { writeln!(changed,"{row}").unwrap(); }
        }
        assert_eq!(original.fsb, import(&changed, SPEC, AUTO).unwrap().fsb);
    }
    #[test]
    fn rotated_millimetre_assets_retain_the_same_crown_mass_and_sections() {
        let obj = fixture(0.015, 1.);
        let mut rotated = String::new();
        for row in obj.lines() {
            if let Some(p) = row.strip_prefix("v ") {
                let p: Vec<f64> = p.split_whitespace().map(|s|s.parse().unwrap()).collect();
                writeln!(rotated,"v {:.17e} {:.17e} {:.17e}",10.+1000.*p[2],20.+1000.*p[0],30.+1000.*p[1]).unwrap();
            } else { writeln!(rotated,"{row}").unwrap(); }
        }
        let spec = SPEC.replace("units,1\n","units,0.001\n")
            .replace("frame,0,0,0,1,0,0,0,1,0","frame,10,20,30,0,1,0,0,0,1");
        let a = import(&obj,SPEC,AUTO).unwrap(); let b = import(&rotated,&spec,AUTO).unwrap();
        let parse = |s: &str| s.lines().filter(|r|r.starts_with("node,") || r.starts_with("triangle,"))
            .map(|r|r.split(',').skip(1).map(|v|v.parse::<f64>().unwrap()).collect::<Vec<_>>()).collect::<Vec<_>>();
        let pa = parse(&a.fsb); let pb = parse(&b.fsb);
        assert_eq!(pa.len(),pb.len());
        for (a,b) in pa.iter().zip(&pb) { for (a,b) in a.iter().zip(b) {
            assert!((a-b).abs()<1e-12*a.abs().max(1.));
        }}
    }
    #[test]
    fn ambiguity_or_missing_correspondence_does_not_become_a_greedy_weld() {
        let mut doc = fs_io::obj::read_obj_document(&fixture(0.,1.)).unwrap();
        let spec = Spec::read(SPEC).unwrap();
        let upper = (1..=5).collect(); let lower = (6..=10).collect();
        doc.soup.positions[6].x = 1e-9;
        let e = projected(&doc,&spec,&upper,&lower,1e-7).unwrap_err();
        assert!(e.contains("ambiguous"),"{e}");
        doc.soup.positions[6].x = 1.0001;
        assert!(projected(&doc,&spec,&upper,&lower,1e-7).unwrap_err().contains("no lower vertex"));
        for spec in [format!("{AUTO}pair,1,6\n"),AUTO.replace("1e-7","1e-2"),
            AUTO.replace("1e-7","NaN"),format!("{AUTO}pairing,projected,1e-7\n")] {
            assert!(import(&fixture(0.,1.),SPEC,&spec).is_err());
        }
    }
}
