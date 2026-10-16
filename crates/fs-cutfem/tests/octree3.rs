use std::ops::ControlFlow;
use fs_cutfem::octree3::{Octree3,OctreeError3};
fn go()->ControlFlow<()> {ControlFlow::Continue(())}
#[test]
fn local_refinement_preserves_leaves_and_exact_face_cover() {
    let tree=Octree3::uniform(1,4,1000).unwrap();
    let mark=*tree.leaves().iter().next().unwrap();
    let next=tree.refined(&[mark,mark],go).unwrap();
    assert_eq!(tree.leaves().len(),8);assert_eq!(next.leaves().len(),15);
    for c in tree.leaves() {if *c!=mark {assert!(next.leaves().contains(c));}}
    let faces=next.faces(go).unwrap();
    let cells:Vec<_>=next.leaves().iter().copied().collect();let mut expected=0;
    for (i,&a) in cells.iter().enumerate() {for &b in &cells[i+1..] {
        let (al,sa)=next.lattice_box(a);let (bl,sb)=next.lattice_box(b);
        let overlap:[i64;3]=std::array::from_fn(|d|i64::from((al[d]+sa).min(bl[d]+sb))-i64::from(al[d].max(bl[d])));
        if overlap.iter().all(|n|*n>=0)&&overlap.iter().filter(|&&n|n==0).count()==1 {expected+=1;}
    }}
    assert_eq!(faces.len(),expected);
    assert_eq!(faces.iter().filter(|f|f.lower.level()!=f.upper.level()).count(),12);
}
#[test]
fn recursive_refinement_closes_edge_and_vertex_balance() {
    let mut tree=Octree3::uniform(1,5,4000).unwrap();
    for _ in 0..3 {
        let mark=*tree.leaves().iter().max_by_key(|c|(c.level(),std::cmp::Reverse(c.index()))).unwrap();
        tree=tree.refined(&[mark],go).unwrap();
    }
    for &a in tree.leaves() {for &b in tree.leaves() {
        let (al,sa)=tree.lattice_box(a);let (bl,sb)=tree.lattice_box(b);
        if (0..3).all(|d|al[d]<=bl[d]+sb&&bl[d]<=al[d]+sa) {
            assert!(a.level().abs_diff(b.level())<=1);
        }
    }}
}
#[test]
fn hanging_rows_reproduce_all_trilinear_monomials() {
    let tree=Octree3::uniform(1,4,1000).unwrap();
    let tree=tree.refined(&[*tree.leaves().iter().next().unwrap()],go).unwrap();
    let rows=tree.constraints(go).unwrap();assert!(!rows.is_empty());
    for (&node,row) in &rows {
        assert_eq!(row.iter().map(|(_,w)|w).sum::<f64>(),1.0);
        for mask in 0..8 {
            let eval=|p:[u32;3]|(0..3).map(|a|if mask&(1<<a)==0 {1.0}else{f64::from(p[a])/f64::from(tree.extent())}).product::<f64>();
            let actual:f64=row.iter().map(|&(p,w)|w*eval(p)).sum();
            assert!((actual-eval(node)).abs()<1e-14);
        }
    }
}
#[test]
fn refusals_preserve_original_and_mark_order_is_irrelevant() {
    let tree=Octree3::uniform(1,3,15).unwrap();let marks:Vec<_>=tree.leaves().iter().copied().take(2).collect();
    assert_eq!(tree.refined(&marks,go).unwrap_err(),OctreeError3::LeafBudget);
    assert_eq!(tree.refined(&marks,||ControlFlow::Break(())).unwrap_err(),OctreeError3::Cancelled);
    assert_eq!(tree.leaves().len(),8);
    let a=Octree3::uniform(1,3,1000).unwrap();let mut reversed=marks.clone();reversed.reverse();
    assert_eq!(a.refined(&marks,go).unwrap().leaves(),a.refined(&reversed,go).unwrap().leaves());
    let terminal=Octree3::uniform(0,0,8).unwrap();let mark=*terminal.leaves().iter().next().unwrap();
    assert_eq!(terminal.refined(&[mark],go).unwrap_err(),OctreeError3::LevelBudget);
}
