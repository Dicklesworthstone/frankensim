//! Direct numerical-model regression references; no total physical-accuracy claim.
use super::*;
use super::super::{AngularDistribution, SourceSpectrum};
use crate::instances::RigidTransform;
use crate::photometry::TransportOptions;
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::f64::consts::{FRAC_PI_2, PI};

pub(super) fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 173, kernel_id: 0x5048_454e, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        operation(&cx, &gate)
    })
}

pub(super) fn probe(x: f64, y: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

pub(super) fn source(flux: f64, height: f64, law: AngularDistribution) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 0.0, height]).unwrap(),
        flux, 1.0, law, SourceSpectrum::equal_energy()).unwrap()
}

pub(super) fn contains(report: EnclosedIlluminance, reference: f64) {
    assert_eq!(report.unavailable, None);
    assert_eq!(report.estimate.numerical.kind, NumericalKind::Enclosure);
    assert!(report.estimate.numerical.lo <= reference, "{report:?}, reference={reference}");
    assert!(report.estimate.numerical.hi >= reference, "{report:?}, reference={reference}");
    assert!((report.estimate.lux-reference).abs() <= report.numerical_error_bound_lux.unwrap());
    assert!(report.estimate.total_error_bound_lux.is_none());
    assert!(report.numerical_model_identity.is_some());
}

#[test]
fn cie171_diffuse_point_case_has_an_enclosure_not_a_total_accuracy_claim() {
    // Same publicly described CIE 171 section 5.2 analytical subcase as
    // tests/photometric_path/cie171_analytical.rs. No reference table copied.
    // I=1000*cos(theta), h=3: E=9000/(9+x*x+y*y)^2.
    with_cx(|cx, _| {
        let sources = [source(1000.0*PI, 3.0,
            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
            TransportOptions::default(), cx).unwrap();
        for [x,y] in [[0.0,0.0],[0.5,0.0],[1.0,0.0],[1.5,0.0],[0.5,0.5],
            [1.0,0.5],[1.5,0.5],[1.0,1.0],[1.5,1.0],[1.5,1.5]] {
            let p = probe(x,y);
            let report = scene.illuminance_with_enclosure(cx, &p).unwrap();
            let squared = 9.0+x*x+y*y;
            contains(report, 9000.0/(squared*squared));
            assert_eq!(report.estimate.lux.to_bits(), scene.illuminance(cx, &p).unwrap().lux.to_bits());
            assert!(report.numerical_error_bound_lux.unwrap() < 1e-9);
        }
    });
}

#[test]
fn gray_haze_encloses_absorption_plus_out_scattering_and_positive_underflow() {
    with_cx(|cx, _| {
        let sources = [source(1000.0*PI, 5.0, AngularDistribution::isotropic())];
        for (a,s) in [(0.2,0.0),(0.0,0.2),(0.07,0.13)] {
            let medium = HomogeneousMedium::try_new(a,s,PhaseFunction::HenyeyGreenstein {g:0.7}).unwrap();
            let scene = PhotometricScene::try_new(&[], &sources, medium,
                TransportOptions::default(), cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx, &probe(0.0,0.0)).unwrap();
            contains(report, 3.678_794_411_714_423_3); // independent 10/e
            assert!(report.numerical_error_bound_lux.unwrap() < 1e-10);
        }
        let scene = PhotometricScene::try_new(&[], &sources,
            HomogeneousMedium::try_new(200.0,0.0,PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(), cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx, &probe(0.0,0.0)).unwrap();
        assert_eq!(report.unavailable, None);
        assert!(report.estimate.numerical.hi > 0.0);
        assert!(report.estimate.numerical.hi < f64::MIN_POSITIVE);
        assert!(report.numerical_error_bound_lux.unwrap() > 0.0);
        assert!(report.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn axisymmetric_beam_tables_and_uncertain_cutoffs_keep_their_support() {
    with_cx(|cx, _| {
        let sources = [source(1000.0, 3.0, AngularDistribution::beam_field(0.2,0.4,0.6).unwrap())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
            TransportOptions::default(), cx).unwrap();
        for x in [0.0,0.1,0.3,0.5,0.7,1.0] {
            let p = probe(x,0.0);
            let report = scene.illuminance_with_enclosure(cx, &p).unwrap();
            contains(report, scene.illuminance(cx, &p).unwrap().lux);
            assert!(report.numerical_error_bound_lux.unwrap() < 1e-8);
        }
    });
    let peak = Interval::point(20.0);
    assert_eq!(cutoff(peak, Interval::new(0.5_f64.next_down(),0.5_f64.next_up()),0.5),
        Interval::new(0.0,20.0));
    assert_eq!(cutoff(peak,Interval::point(0.5),0.5),peak);
    assert_eq!(cutoff(peak,Interval::point(0.5_f64.next_down()),0.5),Interval::point(0.0));
}

#[test]
fn supported_mesh_shadows_do_not_promote_unsupported_angular_laws() {
    use crate::charts::TriMesh;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, Primitive, Shape};
    with_cx(|cx, _| {
        let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
        let geometry = [Primitive {
            shape: Shape::Mesh(TriMesh::new(vec![[-1.0,-1.0,1.0],[1.0,-1.0,1.0],[0.0,1.0,1.0]],vec![[0,1,2]])),
            material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() },
            emission: None,
        }];
        let scene = PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
        contains(report, 0.0);
        assert_eq!(report.estimate.lux,0.0);
        assert_eq!(report.estimate.numerical.lo,0.0);
        assert_eq!(report.estimate.numerical.hi,0.0);
        assert_eq!(report.uncertain_visibility_sources,Some(0));
        let sources = [source(1000.0,3.0,AngularDistribution::cosine_power(0.5,FRAC_PI_2).unwrap())];
        let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
        assert!(report.estimate.lux > 0.0);
        assert_eq!(report.unavailable,Some(EnclosureUnavailable::SourceLaw {source:0}));
        assert_eq!(report.estimate.numerical.kind,NumericalKind::NoClaim);
        assert!(report.uncertain_visibility_sources.is_none());
        assert!(report.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn progressive_grid_cancel_resume_preserves_values_bounds_and_identities() {
    with_cx(|cx,gate| {
        let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let grid = CalculationGrid::try_new(Point3::new(-1.0,-1.0,0.0),Vec3::new(0.5,0.0,0.0),Vec3::new(0.0,0.5,0.0),5,5).unwrap();
        let mut full = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
        full.advance(cx,usize::MAX).unwrap();
        let mut partial = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
        partial.advance(cx,3).unwrap();
        assert_eq!(partial.results(),&full.results()[..3]);
        gate.request();
        assert!(matches!(partial.advance(cx,5),Err(PhotometryError::Cancelled)));
        assert_eq!(partial.results(),&full.results()[..3]);
        with_cx(|resume,_| {
            while !partial.is_complete() { partial.advance(resume,7).unwrap(); }
        });
        assert_eq!(partial.results(),full.results());
    });
}

#[test]
fn source_and_receiver_changes_change_the_numerical_model_identity() {
    with_cx(|cx,_| {
        let a = [source(1000.0,3.0,AngularDistribution::isotropic())];
        let b = [source(2000.0,3.0,AngularDistribution::isotropic())];
        let sa = PhotometricScene::try_new(&[], &a, HomogeneousMedium::vacuum(),TransportOptions::default(), cx).unwrap();
        let sb = PhotometricScene::try_new(&[], &b, HomogeneousMedium::vacuum(),TransportOptions::default(), cx).unwrap();
        let ea = sa.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
        let eb = sb.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
        let ec = sa.illuminance_with_enclosure(cx,&probe(1.0,0.0)).unwrap();
        assert_ne!(ea.numerical_model_identity,eb.numerical_model_identity);
        assert_ne!(ea.numerical_model_identity,ec.numerical_model_identity);
        assert_eq!(ea,sa.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap());
    });
}

#[test]
fn csv_keeps_numerical_and_unknown_physical_errors_separate() {
    with_cx(|cx,_| {
        let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
        let mut bytes = Vec::new();
        EnclosedIlluminance::write_csv_header(&mut bytes).unwrap();
        report.write_csv_row(0,&mut bytes).unwrap();
        let csv = String::from_utf8(bytes).unwrap();
        let rows: Vec<_> = csv.lines().collect();
        let columns: Vec<_> = rows[1].split(',').collect();
        assert_eq!(columns.len(),13);
        assert_eq!(rows[0].split(',').count(),13);
        assert_eq!(columns[5],"enclosure");
        assert_eq!(columns[6].len(),64);
        assert_eq!(&columns[7..11],&["unknown";4]);
        assert_eq!(columns[11],"3");
        assert_eq!(columns[12],"0");
        assert_eq!(columns[1].parse::<f64>().unwrap().to_bits(),report.estimate.lux.to_bits());
        struct Refuse;
        impl Write for Refuse {
            fn write(&mut self,_: &[u8])->io::Result<usize> { Err(io::Error::other("test sink")) }
            fn flush(&mut self)->io::Result<()> { Ok(()) }
        }
        assert!(report.write_csv_row(0,&mut Refuse).is_err());
    });
}
