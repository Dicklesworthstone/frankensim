use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_mesh::{MarkedTetRefinement as Split,TetRefinementError as Error,TetRefinementLimits};
use std::collections::BTreeMap;
fn with_cx<R>(run:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|run(&Cx::new(&gate,arena,
        StreamKey{seed:41,kernel_id:991,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn limits()->TetRefinementLimits { TetRefinementLimits{max_vertices:1000,max_tetrahedra:4000} }
fn mesh()->(Vec<[f64;3]>,Vec<[u32;4]>) {
    (vec![[0.,0.,0.],[2.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.],
        [10.,0.,0.],[11.,0.,0.],[10.,1.,0.],[10.,0.,1.]],
        vec![[0,1,2,3],[0,2,1,4],[5,6,7,8]])
}
fn volume(p:&[[f64;3]],t:[u32;4])->f64 {
    let a:[f64;3]=std::array::from_fn(|i|p[t[1] as usize][i]-p[t[0] as usize][i]);
    let b:[f64;3]=std::array::from_fn(|i|p[t[2] as usize][i]-p[t[0] as usize][i]);
    let c:[f64;3]=std::array::from_fn(|i|p[t[3] as usize][i]-p[t[0] as usize][i]);
    (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])).abs()/6.
}
fn integral(p:&[[f64;3]],t:&[[u32;4]],f:&[f64])->f64 {
    t.iter().map(|&t|volume(p,t)*t.iter().map(|&v|f[v as usize]).sum::<f64>()/4.).sum()
}
#[test]
fn full_edge_star_closes_the_shared_face_without_refining_a_remote_cell() {
    with_cx(|cx| {
        let(p,t)=mesh();let s=Split::build(cx,&p,&t,&[0],&[],limits()).unwrap();
        assert_eq!(s.tetrahedra().len(),5); assert_eq!(s.positions().len(),10);
        assert_eq!(s.tetrahedra()[2],t[2]);
        let mut faces=BTreeMap::new();
        for cell in s.tetrahedra() { for omitted in 0..4 {
            let mut face:Vec<_>=cell.iter().enumerate().filter_map(|(i,&v)|(i!=omitted).then_some(v)).collect();
            face.sort_unstable();*faces.entry(face).or_insert(0)+=1;
        } }
        for mut face in s.face_children([0,1,2]).unwrap() {
            face.sort_unstable(); assert_eq!(faces[&face.to_vec()],2);
        }
        for edge in s.midpoint_parents() {
            assert!(!s.tetrahedra().iter().any(|t|t.contains(&edge[0])&&t.contains(&edge[1])),"hanging original edge");
        }
        assert_eq!(s.parent_elements().iter().filter(|&&i|i==0).count(),2);
        assert_eq!(s.parent_elements().iter().filter(|&&i|i==1).count(),2);
        assert_eq!(s.parent_elements().iter().filter(|&&i|i==2).count(),1);
    });
}
#[test]
fn arbitrary_signed_p1_sources_and_parent_volumes_survive_several_local_rounds() {
    with_cx(|cx| {
        let(mut p,mut t)=mesh();let mut f=vec![1.,-3.,2.,5.,-7.,0.,2.,1.,-1.];
        let original=integral(&p,&t,&f);
        for _ in 0..5 {
            let s=Split::build(cx,&p,&t,&[0,1],&[],limits()).unwrap();
            let next=s.prolongate(cx,&f).unwrap();
            assert!((integral(s.positions(),s.tetrahedra(),&next)-original).abs()<1e-12);
            let mut volumes=vec![0.;t.len()];
            for (&tet,&parent) in s.tetrahedra().iter().zip(s.parent_elements()) {volumes[parent]+=volume(s.positions(),tet);}
            for (i,&tet) in t.iter().enumerate() {assert!((volumes[i]-volume(&p,tet)).abs()<1e-12);}
            f=next;p=s.positions().to_vec();t=s.tetrahedra().to_vec();
        }
    });
}
#[test]
fn contact_constraints_preserve_separate_traces_with_scrambled_vertex_numbers() {
    with_cx(|cx| {
        let p=vec![[0.,0.,0.],[2.,0.,0.],[0.,1.,0.],[0.,0.,1.],
            [0.,1.,0.],[0.,0.,0.],[2.,0.,0.],[0.,0.,-1.]];
        let t=vec![[0,1,2,3],[5,6,4,7]];
        let links=[([0,1],[5,6]),([0,2],[5,4]),([1,2],[6,4])];
        let s=Split::build(cx,&p,&t,&[0],&links,limits()).unwrap();
        assert_eq!(s.positions().len(),10); assert_eq!(s.tetrahedra().len(),4);
        assert_eq!(s.positions()[8],s.positions()[9]);
        let a=s.face_children([0,1,2]).unwrap();let b=s.face_children([5,6,4]).unwrap();
        assert_eq!(a.len(),2);assert_eq!(b.len(),2);
        for face in a {assert!(b.iter().any(|other|face.iter().all(|&v|
            other.iter().any(|&w|s.positions()[v as usize]==s.positions()[w as usize]))));}
        let values=s.prolongate(cx,&[300.,300.,300.,300.,350.,350.,350.,350.]).unwrap();
        assert!([values[8],values[9]].contains(&300.));assert!([values[8],values[9]].contains(&350.));
    });
}
#[test]
fn mark_order_and_duplicates_do_not_change_the_output() {
    with_cx(|cx| {
        let(p,t)=mesh();let a=Split::build(cx,&p,&t,&[0,1],&[],limits()).unwrap();
        let b=Split::build(cx,&p,&t,&[1,0,1],&[],limits()).unwrap();
        assert_eq!(a.positions(),b.positions());assert_eq!(a.tetrahedra(),b.tetrahedra());
        assert_eq!(a.parent_elements(),b.parent_elements());
    });
}
#[test]
fn invalid_constraints_counts_and_cancellation_refuse_the_complete_operation() {
    with_cx(|cx| {
        let(p,t)=mesh();
        assert_eq!(Split::build(cx,&p,&t,&[],&[],limits()).unwrap_err(),Error::InvalidMesh);
        assert_eq!(Split::build(cx,&p,&t,&[99],&[],limits()).unwrap_err(),Error::InvalidMesh);
        assert_eq!(Split::build(cx,&p,&t,&[0],&[([0,1],[5,6])],limits()).unwrap_err(),Error::InvalidMesh);
        for limited in [TetRefinementLimits{max_vertices:9,..limits()},TetRefinementLimits{max_tetrahedra:4,..limits()}] {
            assert_eq!(Split::build(cx,&p,&t,&[0],&[],limited).unwrap_err(),Error::OutputLimit);
        }
    });
    let gate=CancelGate::new_clock_free();gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey{seed:1,kernel_id:991,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let(p,t)=mesh();assert_eq!(Split::build(&cx,&p,&t,&[0],&[],limits()).unwrap_err(),Error::Cancelled);
    });
}
