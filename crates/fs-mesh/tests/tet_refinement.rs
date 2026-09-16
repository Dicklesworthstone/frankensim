use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_mesh::{TetRefinement,TetRefinementError,TetRefinementLimits};
use std::collections::BTreeMap;

fn with_cx<R>(f:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|f(&Cx::new(&gate,arena,
        StreamKey{seed:41,kernel_id:718,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn limits()->TetRefinementLimits {TetRefinementLimits{max_vertices:1000,max_tetrahedra:1000}}
fn mesh()->(Vec<[f64;3]>,Vec<[u32;4]>) {
    (vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]],
     vec![[0,1,2,3],[0,2,1,4]])
}
fn volume(p:&[[f64;3]],t:[u32;4])->f64 {
    let a: [f64;3]=std::array::from_fn(|i|p[t[1] as usize][i]-p[t[0] as usize][i]);
    let b: [f64;3]=std::array::from_fn(|i|p[t[2] as usize][i]-p[t[0] as usize][i]);
    let c: [f64;3]=std::array::from_fn(|i|p[t[3] as usize][i]-p[t[0] as usize][i]);
    (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])).abs()/6.
}
fn integral(p:&[[f64;3]],ts:&[[u32;4]],values:&[f64])->f64 {
    ts.iter().map(|&t|volume(p,t)*t.iter().map(|&v|values[v as usize]).sum::<f64>()/4.).sum()
}
#[test]
fn refinement_shares_face_midpoints_and_preserves_parent_cell_order() {
    with_cx(|cx| {
        let(p,t)=mesh();let fine=TetRefinement::build(cx,&p,&t,limits()).unwrap();
        assert_eq!(fine.positions().len(),14);assert_eq!(fine.tetrahedra().len(),16);
        assert_eq!(&fine.positions()[..p.len()],p.as_slice());
        for (parent,children) in t.iter().zip(fine.tetrahedra().chunks_exact(8)) {
            let sum:f64=children.iter().map(|&child|volume(fine.positions(),child)).sum();
            assert!((sum-volume(&p,*parent)).abs()<1e-14);
        }
        let mut counts=BTreeMap::new();
        for cell in fine.tetrahedra() {
            for omit in 0..4 {
                let mut face:Vec<u32>=cell.iter().enumerate().filter_map(|(i,&v)|(i!=omit).then_some(v)).collect();
                face.sort_unstable();*counts.entry(face).or_insert(0)+=1;
            }
        }
        for child in fine.face_children([0,1,2]).unwrap() {
            let mut face=child.to_vec();face.sort_unstable();assert_eq!(counts[&face],2);
        }
    });
}
#[test]
fn prolongation_preserves_piecewise_linear_source_integrals_without_renormalizing() {
    with_cx(|cx| {
        let(mut p,mut t)=mesh();let mut field=vec![1.,5.,0.,3.,-4.];
        let original=integral(&p,&t,&field);
        for _ in 0..3 {
            let fine=TetRefinement::build(cx,&p,&t,TetRefinementLimits{max_vertices:10000,max_tetrahedra:10000}).unwrap();
            field=fine.prolongate(cx,&field).unwrap();p=fine.positions().to_vec();t=fine.tetrahedra().to_vec();
            assert!((integral(&p,&t,&field)-original).abs()<1e-12);
        }
    });
}
#[test]
fn affine_fields_interpolate_at_the_generated_coordinates() {
    with_cx(|cx| {
        let(p,t)=mesh();let f=|p:[f64;3]|1.+2.*p[0]-3.*p[1]+4.*p[2];
        let fine=TetRefinement::build(cx,&p,&t,limits()).unwrap();
        let values=fine.prolongate(cx,&p.iter().copied().map(f).collect::<Vec<_>>()).unwrap();
        for (&point,&value) in fine.positions().iter().zip(&values){assert_eq!(value,f(point));}
    });
}
#[test]
fn coincident_contact_sides_stay_separate_and_transfer_their_own_fields() {
    with_cx(|cx| {
        let p=vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],
                   [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.]];
        let t=vec![[0,1,2,3],[4,6,5,7]];
        let fine=TetRefinement::build(cx,&p,&t,limits()).unwrap();
        assert_eq!(fine.positions().len(),20);
        let values=fine.prolongate(cx,&[0.,0.,0.,0.,10.,10.,10.,10.]).unwrap();
        for (a,b) in fine.face_children([0,1,2]).unwrap().iter().zip(fine.face_children([4,5,6]).unwrap()) {
            for (&va,vb) in a.iter().zip(b) {
                assert_ne!(va,vb);assert_eq!(fine.positions()[va as usize],fine.positions()[vb as usize]);
                assert_eq!(values[va as usize],0.);assert_eq!(values[vb as usize],10.);
            }
        }
    });
}
#[test]
fn output_limits_and_malformed_fields_refuse() {
    with_cx(|cx| {
        let(p,t)=mesh();
        assert!(matches!(TetRefinement::build(cx,&p,&t,TetRefinementLimits{max_vertices:13,..limits()}),Err(TetRefinementError::OutputLimit)));
        assert!(matches!(TetRefinement::build(cx,&p,&t,TetRefinementLimits{max_tetrahedra:15,..limits()}),Err(TetRefinementError::OutputLimit)));
        assert!(TetRefinement::build(cx,&p,&[[0,0,2,3]],limits()).is_err());
        let fine=TetRefinement::build(cx,&p,&t,limits()).unwrap();
        assert!(fine.prolongate(cx,&[0.]).is_err());assert!(fine.prolongate(cx,&[f64::NAN;5]).is_err());
        assert!(fine.face_children([1,3,4]).is_err());
    });
}
#[test]
fn cancellation_precedes_refinement() {
    let gate=CancelGate::new_clock_free();gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey{seed:41,kernel_id:718,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let(p,t)=mesh();assert!(matches!(TetRefinement::build(&cx,&p,&t,limits()),Err(TetRefinementError::Cancelled)));
    });
}
