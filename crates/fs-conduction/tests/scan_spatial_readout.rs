//! Spatial readouts must sample the actual P1 field and transpose the same row.
mod support;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::spatial::ScanReadout;
use fs_conduction::{ConductionMesh, fixtures::box_grid};
use support::{with_cx, with_cancelled_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2, 1, 1], [2.0, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-13*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}"); }

#[test]
fn g1_interior_probe_reproduces_an_affine_field_without_node_snapping() {
    let mesh = mesh();
    with_cx(|cx| {
        let readout = ScanReadout::in_tetrahedron(cx, &mesh, 0, [1.0,2.0,3.0,4.0]).unwrap();
        let affine = |p: [f64;3]| 20.0+p[0]+2.0*p[1]+3.0*p[2];
        let field: Vec<_> = mesh.positions().iter().map(|&p| affine(p)).collect();
        let mut position = [0.0;3];
        for (&vertex, weight) in mesh.complex().tets[0].iter().zip([0.1,0.2,0.3,0.4]) {
            for axis in 0..3 { position[axis] += weight*mesh.positions()[vertex as usize][axis]; }
        }
        near(readout.evaluate(cx, &field).unwrap(), affine(position));
        assert_eq!(readout.weights().len(), 4);
        assert!(std::ptr::eq(readout.mesh(), &mesh));
        near(ScanReadout::vertex(cx, &mesh, 0).unwrap().evaluate(cx, &field).unwrap(), field[0]);
    });
}

#[test]
fn g1_region_mean_uses_element_measure_and_not_vertex_count() {
    let mesh = mesh();
    with_cx(|cx| {
        let elements = [0,1,2];
        let readout = ScanReadout::volume_mean(cx, &mesh, &elements, 3).unwrap();
        let field: Vec<_> = (0..mesh.vertex_count()).map(|i| (i*i) as f64).collect();
        let total: f64 = elements.iter().map(|&e| mesh.element_volume(e)).sum();
        let expected: f64 = elements.iter().map(|&e| mesh.element_volume(e)*
            mesh.complex().tets[e].iter().map(|&v| field[v as usize]/4.0).sum::<f64>()).sum::<f64>()/total;
        near(readout.evaluate(cx, &field).unwrap(), expected);
        let node_average = readout.weights().iter().map(|&(v,_)| field[v]).sum::<f64>()/readout.weights().len() as f64;
        assert!((expected-node_average).abs()>1e-3);
        let masses = ScanReadout::weighted_mean(cx, &mesh, &[(0,2.0),(mesh.vertex_count()-1,1.0)], 2).unwrap();
        let mut fractions = vec![0.0; mesh.vertex_count()];
        fractions[mesh.vertex_count()-1] = 1.0;
        near(masses.evaluate(cx, &fractions).unwrap(), 1.0/3.0);
    });
}

#[test]
fn g3_spatial_row_transpose_and_weight_rescaling_are_consistent() {
    let mesh = mesh();
    with_cx(|cx| {
        let a = ScanReadout::weighted_mean(cx,&mesh,&[(0,1.0),(1,0.0),(2,3.0)],3).unwrap();
        let b = ScanReadout::weighted_mean(cx,&mesh,&[(0,1e300),(1,0.0),(2,3e300)],3).unwrap();
        assert_eq!(a.weights(), b.weights());
        assert_eq!(a.weights().len(), 2);
        let field: Vec<_> = (0..mesh.vertex_count()).map(|i| i as f64-3.0).collect();
        for seed in [-4.0,0.0,2.5] {
            let gradient = a.pullback(cx, seed).unwrap();
            near(seed*a.evaluate(cx,&field).unwrap(),gradient.iter().zip(&field).map(|(a,b)|a*b).sum());
            assert_eq!(gradient[1],0.0);
        }
    });
}

#[test]
fn g0_spatial_readouts_refuse_bad_support_and_honor_cancellation() {
    let mesh = mesh();
    with_cx(|cx| {
        for weights in [vec![],vec![(0,0.0)],vec![(0,-1.0)],vec![(0,f64::NAN)],
            vec![(0,1.0),(0,2.0)],vec![(mesh.vertex_count(),1.0)],vec![(0,f64::from_bits(1)),(1,1e300)]]
        { assert!(ScanReadout::weighted_mean(cx,&mesh,&weights,10).is_err()); }
        assert!(ScanReadout::weighted_mean(cx,&mesh,&[(0,1.0),(1,1.0)],1).is_err());
        assert!(ScanReadout::in_tetrahedron(cx,&mesh,mesh.element_count(),[1.0;4]).is_err());
        assert!(ScanReadout::in_tetrahedron(cx,&mesh,0,[-1.0,1.0,1.0,1.0]).is_err());
        assert!(ScanReadout::volume_mean(cx,&mesh,&[0,0],2).is_err());
        assert!(ScanReadout::volume_mean(cx,&mesh,&[0,1],1).is_err());
        let readout = ScanReadout::vertex(cx,&mesh,0).unwrap();
        assert!(readout.evaluate(cx,&[]).is_err());
        assert!(readout.pullback(cx,f64::NAN).is_err());
    });
    with_cancelled_cx(|cx| assert!(ScanReadout::vertex(cx,&mesh,0).is_err()));
}
