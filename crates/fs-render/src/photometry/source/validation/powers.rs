//! Real cosine-power references through original fixture, mesh and medium owners.
//! No measured-fixture, constructor-normalization or total-accuracy claim.

use super::*;
use super::tests::{contains, probe, source, with_cx};
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
use crate::photometry::source::{AngularDistribution, SourceSpectrum, SPECTRAL_BINS, photopic_weight_nm, wavelength_nm};
use crate::photometry::transport::PhotometricMedium;
use crate::photometry::TransportOptions;
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn established_small_integer_interval_bits_are_unchanged() {
    for exponent in 0..=64 {
        for mu in [Interval::point(0.0), Interval::point(0.5), Interval::new(0.1,0.9), Interval::point(1.0)] {
            let mut original = Interval::point(1.0);
            for _ in 0..exponent { original = nonnegative(original * nonnegative(mu)); }
            let actual = cosine_power_range(mu, f64::from(exponent)).unwrap();
            assert_eq!(actual.lo().to_bits(), original.lo().to_bits());
            assert_eq!(actual.hi().to_bits(), original.hi().to_bits());
        }
    }
}

#[test]
fn fractional_fixture_lux_matches_independent_half_integer_formulae() {
    with_cx(|cx, _| {
        for integer in [0, 1, 65] {
            let exponent = f64::from(integer) + 0.5;
            // Hemisphere flux=2*pi*I0/(p+1), I0=1000 cd. At (3,0,0)
            // beneath a source at (0,0,4), r=5, mu=4/5, E=32*(4/5)^p.
            let sources = [source(2000.0*PI/(exponent+1.0), 4.0,
                AngularDistribution::cosine_power(exponent, FRAC_PI_2).unwrap())];
            let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
                TransportOptions::default(), cx).unwrap();
            let p = probe(3.0,0.0);
            let report = scene.illuminance_with_enclosure(cx,&p).unwrap();
            let mut reference = 32.0 * 0.8_f64.sqrt();
            for _ in 0..integer { reference *= 0.8; }
            contains(report,reference);
            assert!(report.numerical_error_bound_lux.unwrap() < reference * 1e-9);
            assert_eq!(report.estimate.lux.to_bits(),scene.illuminance(cx,&p).unwrap().lux.to_bits());
        }
    });
}

#[test]
fn narrow_and_above_u32_exponents_have_bounded_work_and_preserve_admission() {
    with_cx(|cx, _| {
        for exponent in [73.25, 1_000_000.5, 4_294_967_296.5, 1e12] {
            let law = AngularDistribution::cosine_power(exponent,FRAC_PI_2).unwrap();
            let sources = [source(2000.0*PI/(exponent+1.0), 3.0,law)];
            let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),
                TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            contains(report,1000.0/9.0);
            assert!(report.numerical_error_bound_lux.unwrap().is_finite());
            assert_eq!(report.uncertain_visibility_sources,Some(0));
        }
        // The source owner, not the enclosure, decides whether the original
        // beam width can be represented. No rejection is bypassed here.
        assert!(AngularDistribution::cosine_power(f64::MAX,FRAC_PI_2).is_err());
    });
}

#[test]
fn fractional_cutoff_ambiguity_contains_both_branches_and_underflow_is_not_black() {
    let mu = Interval::new(0.5_f64.next_down(),0.5_f64.next_up());
    let value = cosine_power_range(mu,0.5).unwrap();
    let branch = cutoff(value,mu,0.5);
    assert_eq!(branch.lo(),0.0);
    assert!(branch.contains(0.5_f64.sqrt()));
    assert_eq!(cutoff(value,Interval::point(0.5_f64.next_down()),0.5),Interval::point(0.0));
    let tail = cosine_power_range(Interval::point(0.5),2000.5).unwrap();
    assert_eq!(tail.lo(),0.0);
    assert!(tail.hi()>0.0 && tail.hi()<f64::MIN_POSITIVE);
    assert_eq!(cosine_power_range(Interval::point(0.0),0.5),Some(Interval::point(0.0)));
    assert_eq!(cosine_power_range(Interval::point(0.0),0.0),Some(Interval::point(1.0)));
}

#[test]
fn exponents_beyond_old_limit_bind_provenance_and_csv_without_total_error() {
    with_cx(|cx, _| {
        let mut identities = Vec::new();
        for exponent in [64.25,64.5] {
            let sources = [source(1000.0,3.0,
                AngularDistribution::cosine_power(exponent,FRAC_PI_2).unwrap())];
            let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),
                TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.25,0.0)).unwrap();
            assert_eq!(report.unavailable,None);
            identities.push(report.numerical_model_identity);
            let mut bytes = Vec::new();report.write_csv_row(0,&mut bytes).unwrap();
            let csv=String::from_utf8(bytes).unwrap();
            let fields: Vec<_>=csv.trim_end().split(',').collect();
            assert_eq!(fields[5],"enclosure");
            assert_eq!(&fields[7..11],&["unknown";4]);
            assert_eq!(fields[11],"5");
            assert!(report.estimate.total_error_bound_lux.is_none());
        }
        assert_ne!(identities[0],identities[1]);
    });
}

#[test]
fn fractional_beams_share_colored_grid_shadows_and_progressive_cancel_resume() {
    use crate::charts::TriMesh;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material,Primitive,Shape};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid,AbsorbingGridBudget,MassAbsorptionSpectrum};
    use crate::volumes::transport::scattering::{MassScatteringSpectrum,ScatteringGrid};
    with_cx(|cx,gate| {
        let density=[1.0,2.0,4.0];
        let absorption=AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,3],&density,[-0.5,-0.5,1.0],[1.0;3]),
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)]).unwrap(),None,
            AbsorbingGridBudget{max_grid_cells:3,max_steps_per_segment:128}).unwrap();
        let grid=ScatteringGrid::try_new(cx,absorption,
            MassScatteringSpectrum::rayleigh(400.0,0.5).unwrap(),PhaseFunction::Rayleigh).unwrap();
        let exponent=73.5;
        let sources=[source(2000.0*PI/(exponent+1.0),5.0,
            AngularDistribution::cosine_power(exponent,FRAC_PI_2).unwrap())];
        let material=Material::SpectralDiffuse{reflectance:SpectralReflectance::constant(0.0).unwrap()};
        let geometry=[
            Primitive{shape:Shape::Mesh(TriMesh::new(
                vec![[-0.1,-0.1,2.5],[0.1,-0.1,2.5],[0.0,0.1,2.5]],vec![[0,1,2]])),material,emission:None},
            Primitive{shape:Shape::Mesh(TriMesh::new(
                vec![[-2.0,-2.0,0.0],[2.0,-2.0,0.0],[0.0,2.0,0.0]],vec![[0,1,2]])),material,emission:None},
        ];
        let scene=PhotometricScene::try_with_medium(&geometry,&sources,PhotometricMedium::Scattering(&grid),
            TransportOptions::default(),cx).unwrap();
        let probes=[probe(0.0,0.0),probe(0.4,0.0),probe(0.5,0.0)];
        let mut full=EnclosedIlluminanceJob::new(&scene,&probes,cx).unwrap();
        full.advance(cx,usize::MAX).unwrap();
        let mut split=EnclosedIlluminanceJob::new(&scene,&probes,cx).unwrap();
        split.advance(cx,1).unwrap();
        gate.request();assert!(matches!(split.advance(cx,1),Err(PhotometryError::Cancelled)));
        assert_eq!(split.results(),&full.results()[..1]);
        with_cx(|resume,_| {split.advance(resume,usize::MAX).unwrap();});
        assert_eq!(split.results(),full.results());
        assert_eq!(full.results()[0].estimate.numerical.lo,0.0);
        assert_eq!(full.results()[0].estimate.numerical.hi,0.0);
        for (index,x) in [(1,0.4_f64),(2,0.5)] {
            let radius=(25.0+x*x).sqrt();let cosine=5.0/radius;
            // Independent finite layer sum and half-integer angular law.
            // No nominal transport, native pow, or enclosure helper is used.
            let mut beam=cosine.sqrt();for _ in 0..73 {beam*=cosine;}
            let mut spectral=0.0;
            for bin in 0..SPECTRAL_BINS {
                let lambda=wavelength_nm(bin);let ratio=400.0/lambda;
                let extinction=(lambda-380.0)/400.0+0.5*ratio*ratio*ratio*ratio;
                spectral+=sources[0].spectrum().bins_w_per_lumen_nm()[bin]
                    *photopic_weight_nm(bin)*(-7.0*radius/5.0*extinction).exp();
            }
            let reference=1000.0*beam*cosine/(radius*radius)*683.0*spectral;
            let report=full.results()[index];contains(report,reference);
            assert_eq!(report.uncertain_visibility_sources,Some(0));
            assert!(report.numerical_error_bound_lux.unwrap()<reference*1e-8);
        }
    });
}

#[test]
fn noninteger_source_and_receiver_rotations_reuse_native_fixture_coordinates() {
    use crate::instances::RigidTransform;
    use fs_geom::{Point3,Vec3};
    with_cx(|cx,_| {
        let q=0.5_f64.sqrt();
        let sources=[PhotometricSource::try_new(
            RigidTransform::try_new([0.0,q,0.0,q],[0.0,0.0,0.0]).unwrap(),
            2000.0*PI/1.5,1.0,AngularDistribution::cosine_power(0.5,FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy()).unwrap()];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),
            TransportOptions::default(),cx).unwrap();
        let receiver=IlluminanceProbe::try_new(Point3::new(4.0,3.0,0.0),Vec3::new(-1.0,0.0,0.0)).unwrap();
        let report=scene.illuminance_with_enclosure(cx,&receiver).unwrap();
        contains(report,32.0*0.8_f64.sqrt());
        assert!(report.numerical_error_bound_lux.unwrap()<1e-9);
        assert_eq!(report.estimate.lux.to_bits(),scene.illuminance(cx,&receiver).unwrap().lux.to_bits());
    });
}
