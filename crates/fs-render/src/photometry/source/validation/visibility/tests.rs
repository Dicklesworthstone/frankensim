use super::*;
use super::super::{EnclosedIlluminanceJob, PhotometricScene, direct_bound, prepare, Prepared};
use super::super::tests::{contains, probe, source, with_cx};
use crate::instances::GeometryInstance;
use crate::photometry::{CalculationGrid, IlluminanceProbe, TransportOptions};
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::spectral::SpectralReflectance;
use crate::tracer::Material;
use fs_geom::{Point3, Vec3};
use std::f64::consts::PI;

fn point(value: [f64; 3]) -> Vector { value.map(Interval::point) }
fn screen(z: f64) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(vec![[-2.0,-1.0,z],[2.0,-1.0,z],[0.0,2.0,z]],vec![[0,1,2]])),
        material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.3).unwrap() },
        emission: None,
    }
}
fn index(primitives: &[Primitive], cx: &Cx<'_>) -> VisibilityIndex {
    VisibilityIndex::prepare(primitives,&mut DomainHasher::new("test.visibility"),cx).unwrap().unwrap()
}

#[test]
fn independent_triangle_cases_keep_winding_edges_and_endpoints_honest() {
    let triangle = [[-2.0,-1.0,2.0],[2.0,-1.0,2.0],[0.0,2.0,2.0]].map(point);
    let origin = point([0.0,0.0,4.0]);
    let direction = point([0.0,0.0,-4.0]);
    let end = Interval::point(0.999);
    assert_eq!(triangle_visibility(origin,direction,end,triangle),Visibility::Blocked);
    assert_eq!(triangle_visibility(origin,direction,end,[triangle[0],triangle[2],triangle[1]]),Visibility::Blocked);
    assert_eq!(triangle_visibility(origin,point([3.0,0.0,-4.0]),end,triangle),Visibility::Clear);
    // The ray goes exactly through vertex (2,-1,2): do not guess either side.
    assert_eq!(triangle_visibility(origin,point([4.0,-2.0,-4.0]),end,triangle),Visibility::Uncertain);
    assert_eq!(triangle_visibility(origin,direction,Interval::point(0.5),triangle),Visibility::Uncertain);
    assert_eq!(triangle_visibility(origin,direction,Interval::point(0.49),triangle),Visibility::Clear);
    assert_eq!(triangle_visibility(origin,direction,Interval::point(0.51),triangle),Visibility::Blocked);
    assert_eq!(triangle_visibility(point([0.0,0.0,2.0]),point([1.0,0.0,0.0]),end,triangle),Visibility::Uncertain);
}

#[test]
fn native_shadow_and_receiver_floor_enclose_independent_point_lux() {
    with_cx(|cx,_| {
        let geometry = [screen(2.0),screen(0.0),screen(-1.0),screen(5.0)];
        let sources = [source(64.0*PI,4.0,AngularDistribution::isotropic())];
        for sigma in [0.0,0.25] {
            let scene = PhotometricScene::try_new(&geometry,&sources,
                HomogeneousMedium::try_new(sigma,0.0,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
            let prepared = prepare(&scene,cx).unwrap();
            let Prepared::Available { visibility, .. } = &prepared else { panic!("native meshes must be admitted"); };
            for (x,reference) in [(0.0,0.0),(3.0,0.512_f64*(-5.0_f64*sigma).exp())] {
                let p = probe(x,0.0);
                let report = scene.illuminance_with_enclosure(cx,&p).unwrap();
                contains(report,reference);
                // Check the independent interval BEFORE nominal-value hulling.
                let (raw,uncertain) = direct_bound(&scene,cx,&p,visibility).unwrap().unwrap();
                assert!(raw.contains(reference),"{raw:?}, reference={reference}");
                assert_eq!(uncertain,0);
                assert_eq!(report.uncertain_visibility_sources,Some(0));
                assert_eq!(report.estimate.lux.to_bits(),scene.illuminance(cx,&p).unwrap().lux.to_bits());
                assert!(report.numerical_error_bound_lux.unwrap() < 1e-10);
            }
        }
        // A probe directly on a floor is not its own blocker: the endpoint
        // exclusion is honored, not replaced by a normal offset or a full ray.
        let floors = [screen(0.0),screen(-1.0),screen(5.0)];
        let scene = PhotometricScene::try_new(&floors,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        contains(scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap(),1.0);
    });
}

#[test]
fn uncertain_vertex_widens_lux_and_csv_without_losing_a_later_blocker() {
    with_cx(|cx,_| {
        let geometry = [screen(2.0)];
        let sources = [source(64.0*PI,4.0,AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let p = probe(4.0,-2.0);
        let report = scene.illuminance_with_enclosure(cx,&p).unwrap();
        // r=6, cos(theta)=2/3, I=16: unblocked E = 8/27 lux.
        assert_eq!(report.uncertain_visibility_sources,Some(1));
        assert_eq!(report.estimate.numerical.lo,0.0);
        assert!(report.estimate.numerical.hi >= 8.0/27.0);
        assert!(report.estimate.total_error_bound_lux.is_none());
        let mut csv = Vec::new();
        report.write_csv_row(0,&mut csv).unwrap();
        assert!(String::from_utf8(csv).unwrap().ends_with(",2,1\n"));
        // Add a triangle centered on a subsequent point of that same segment.
        let mut blocker = screen(1.0);
        if let Shape::Mesh(mesh) = &mut blocker.shape {
            // Rebuild the native BVH after changing the coordinates.
            let vertices = mesh.vertices.iter().map(|v| [v[0]+3.0,v[1]-1.5,v[2]]).collect();
            *mesh = TriMesh::new(vertices,mesh.triangles.clone());
        }
        let both = [screen(2.0),blocker];
        let scene = PhotometricScene::try_new(&both,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx,&p).unwrap();
        contains(report,0.0);
        assert_eq!(report.estimate.numerical.hi,0.0);
        assert_eq!(report.uncertain_visibility_sources,Some(0));
    });
}

#[test]
fn placed_mesh_geometry_and_not_caller_labels_bind_the_identity() {
    with_cx(|cx,_| {
        let supplied = fs_blake3::hash_domain("test.geometry",b"same-caller-label");
        let make = |height: f64| {
            let primitive = screen(height);
            let Shape::Mesh(mesh) = primitive.shape else { unreachable!() };
            // Exact binary unit quaternion cyclically permutes x,y,z.
            Primitive {
                shape: Shape::Instance(GeometryInstance::try_new(1,supplied,SharedGeometry::mesh(mesh),
                    RigidTransform::try_new([0.5;4],[0.0;3]).unwrap()).unwrap()),
                material: primitive.material,
                emission: None,
            }
        };
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0,0.0,0.0,1.0],[4.0,0.0,0.0]).unwrap(),
            64.0*PI,1.0,AngularDistribution::isotropic(),SourceSpectrum::equal_energy()).unwrap()];
        let p = IlluminanceProbe::try_new(Point3::new(0.0,0.0,0.0),Vec3::new(1.0,0.0,0.0)).unwrap();
        let a = [make(2.0)];
        let b = [make(5.0)];
        let sa = PhotometricScene::try_new(&a,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let sb = PhotometricScene::try_new(&b,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let ea = sa.illuminance_with_enclosure(cx,&p).unwrap();
        let eb = sb.illuminance_with_enclosure(cx,&p).unwrap();
        contains(ea,0.0);
        contains(eb,1.0);
        assert_ne!(ea.numerical_model_identity,eb.numerical_model_identity);
        assert_eq!(ea.uncertain_visibility_sources,Some(0));
        assert_eq!(eb.uncertain_visibility_sources,Some(0));
    });
}

#[test]
fn stackless_hierarchy_matches_all_native_triangles_without_dropping_leaves() {
    with_cx(|cx,_| {
        let mut geometry: Vec<_> = (0..257).map(|i| screen(f64::from(i)-128.0)).collect();
        // Exercise unsorted primitive order and a hierarchy with several levels.
        geometry.rotate_left(93);
        let prepared = index(&geometry,cx);
        assert_eq!(prepared.nodes[0].escape,prepared.nodes.len());
        assert_eq!(prepared.nodes.iter().map(|n| n.count).sum::<usize>(),257);
        for x in [-10.0,-3.0,0.0,3.0,10.0] {
            for y in [-2.0,0.0,2.0] {
                let origin = point([x,y,4.5]);
                let direction = point([0.0,0.0,-4.0]);
                let end = Interval::point(0.99);
                let mut brute = Visibility::Clear;
                for primitive in &geometry {
                    let triangle = world_triangle(primitive,0).unwrap();
                    // Same narrow phase, no hierarchy: tests pruning, while
                    // the other tests and rational oracle test the predicate.
                    if !Bounds::triangle(triangle).unwrap().may_intersect(origin,direction,end) { continue; }
                    match triangle_visibility(origin,direction,end,triangle) {
                        Visibility::Blocked => { brute = Visibility::Blocked; break; }
                        Visibility::Uncertain => brute = Visibility::Uncertain,
                        Visibility::Clear => {}
                    }
                }
                assert_eq!(prepared.visibility(&geometry,origin,direction,end,cx).unwrap(),brute);
            }
        }
        // Interval direction includes zero in x: this box is still provably
        // unreachable, whereas an x slab touching the swept range is retained.
        let bounds = Bounds { lo:[2.0,-1.0,-1.0],hi:[3.0,1.0,1.0] };
        assert!(!bounds.may_intersect(point([0.0;3]),[Interval::new(-0.1,0.1),Interval::point(0.0),Interval::point(1.0)],Interval::point(1.0)));
    });
}

#[test]
fn bounded_admission_and_cancelled_mesh_jobs_preserve_the_accepted_prefix() {
    with_cx(|cx,gate| {
        let mut too_large = screen(2.0);
        if let Shape::Mesh(mesh) = &mut too_large.shape {
            // Preparation must reject the metadata size before consulting the
            // old BVH. This fixture is never passed to nominal transport.
            mesh.triangles = vec![[0,1,2]; MAX_VISIBILITY_TRIANGLES+1];
        }
        assert!(matches!(VisibilityIndex::prepare(&[too_large],&mut DomainHasher::new("test.budget"),cx).unwrap(),
            Err(EnclosureUnavailable::GeometryBudgetExceeded)));
        let geometry = [screen(2.0),screen(0.0)];
        let sources = [source(64.0*PI,4.0,AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let grid = CalculationGrid::try_new(Point3::new(-3.0,-2.0,0.0),Vec3::new(0.5,0.0,0.0),Vec3::new(0.0,0.5,0.0),13,9).unwrap();
        let mut full = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
        full.advance(cx,usize::MAX).unwrap();
        let mut split = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
        split.advance(cx,7).unwrap();
        gate.request();
        assert!(matches!(split.advance(cx,9),Err(PhotometryError::Cancelled)));
        assert_eq!(split.results(),&full.results()[..7]);
        assert!(matches!(EnclosedIlluminanceJob::for_grid(&scene,grid,cx),Err(PhotometryError::Cancelled)));
        with_cx(|resume,_| {
            while !split.is_complete() { split.advance(resume,11).unwrap(); }
        });
        assert_eq!(split.results(),full.results());
    });
}
