use super::*;
use super::super::tests::{contains, probe, source, with_cx};
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
use crate::photometry::source::{AngularDistribution, SPECTRAL_BINS, photopic_weight_nm, wavelength_nm};
use crate::photometry::source::validation::EnclosedIlluminanceJob;
use crate::photometry::{PhotometricScene, TransportOptions};
use crate::photometry::transport::PhotometricMedium;
use std::f64::consts::{FRAC_PI_2, PI};

fn asymmetric() -> AngularDistribution {
    AngularDistribution::tabulated(
        vec![0.0, FRAC_PI_2, PI],
        vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
        vec![1.0; 5].into_iter()
            .chain([1.0, 3.0, 1.0, 2.0, 1.0])
            .chain([1.0; 5]).collect(),
    ).unwrap()
}

#[test]
fn seam_arcs_remain_local_and_poles_retain_all_azimuths() {
    let arcs = azimuth_arcs(Interval::point(1.0), Interval::new(-1e-12, 1e-12));
    assert_eq!(arcs[0].unwrap().lo(), 0.0);
    assert!(arcs[0].unwrap().hi() < 2e-12);
    assert!(arcs[1].unwrap().lo() > TAU - 2e-12);
    assert_eq!(arcs[1].unwrap().hi(), TAU);
    assert_eq!(
        azimuth_arcs(Interval::new(-1e-12, 1e-12), Interval::point(0.0)),
        [Some(Interval::new(0.0, TAU)), None],
    );
    with_cx(|cx, _| {
        let table = AngularTable {
            cosine: vec![1.0, 0.0, -1.0],
            phi: vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
            values: [vec![1.0; 5], vec![2.0, 8.0, 16.0, 4.0, 2.0], vec![1.0; 5]].concat(),
        };
        let local = [Interval::point(1.0), Interval::new(-1e-12, 1e-12), Interval::point(0.0)];
        let value = table_enclosure(&table, 16.0, local, Interval::point(0.0), cx).unwrap();
        assert!(value.contains(2.0));
        assert!(value.width() < 2e-10, "seam included unrelated lobes: {value:?}");
        let pole = table_enclosure(&table, 16.0, [Interval::point(0.0); 3], Interval::point(1.0), cx).unwrap();
        assert_eq!(pole, Interval::point(1.0));
    });
}

#[test]
fn table_cell_work_limit_returns_a_conservative_envelope_and_polls_cancellation() {
    with_cx(|cx, gate| {
        let columns = MAX_ANGULAR_CELLS_PER_ENCLOSURE + 2;
        let table = AngularTable {
            cosine: vec![1.0, -1.0],
            phi: (0..columns).map(|i| TAU * i as f64 / (columns - 1) as f64).collect(),
            values: vec![7.0; 2 * columns],
        };
        let result = table_enclosure(&table, 7.0, [Interval::point(0.0); 3], Interval::new(-1.0, 1.0), cx).unwrap();
        assert_eq!(result, Interval::new(0.0, 7.0));
        gate.request();
        assert!(matches!(table_enclosure(&table, 7.0, [Interval::point(0.0); 3], Interval::point(1.0), cx), Err(PhotometryError::Cancelled)));
    });
}

#[test]
fn bilinear_ranges_preserve_exact_black_and_large_constant_cells() {
    for value in [0.0, 0.25, f64::MAX] {
        for t in [0.0, 0.125, 0.5, 0.875, 1.0] {
            assert_eq!(convex(Interval::point(value), Interval::point(value), t), Interval::point(value));
        }
    }
    // Independent saddle f(t,u) = 10*(t+u-2*t*u), whose extrema on
    // [1/4,3/4]^2 are 3.75 and 6.25, not its midpoint value 5.
    let mut range: Option<Interval> = None;
    for u in [0.25, 0.75] {
        for t in [0.25, 0.75] {
            let v = convex(convex(Interval::point(0.0), Interval::point(10.0), u),
                convex(Interval::point(10.0), Interval::point(0.0), u), t);
            range = Some(range.map_or(v, |old| old.hull(v)));
        }
    }
    assert!(range.unwrap().contains(3.75) && range.unwrap().contains(6.25));
}

#[test]
fn asymmetric_fixture_lux_matches_independent_solid_angle_normalization() {
    with_cx(|cx, _| {
        let sources = [source(5500.0 * PI, 3.0, asymmetric())];
        for extinction in [0.0, 0.125] {
            let medium = HomogeneousMedium::try_new(extinction, 0.0, PhaseFunction::Isotropic).unwrap();
            let scene = PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx).unwrap();
            for (x, y, equator) in [(3.0, 0.0, 1.0), (-3.0, 0.0, 1.0), (0.0, -3.0, 3.0), (0.0, 3.0, 2.0)] {
                // Independent full-sphere integral: each pole has azimuth
                // integral 2*pi; equator integral is 3.5*pi. Linear polar
                // integration gives total 5.5*pi regardless of middle mu.
                let r = 18.0_f64.sqrt();
                let profile = 1.0 + (1.0 - 3.0 / r) * (equator - 1.0);
                let reference = 1000.0 * profile * 3.0 / (18.0 * r) * (-extinction * r).exp();
                let p = probe(x, y);
                let report = scene.illuminance_with_enclosure(cx, &p).unwrap();
                contains(report, reference);
                assert!(report.numerical_error_bound_lux.unwrap() < 1e-8);
                assert_eq!(report.estimate.lux.to_bits(), scene.illuminance(cx, &p).unwrap().lux.to_bits());
            }
        }
    });
}

#[test]
fn asymmetric_jobs_replay_and_bind_azimuth_data() {
    with_cx(|cx, gate| {
        let sources = [source(5500.0 * PI, 3.0, asymmetric())];
        let probes = [probe(0.0, -3.0), probe(0.0, 3.0), probe(3.0, 1e-12), probe(3.0, -1e-12)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let mut full = EnclosedIlluminanceJob::new(&scene, &probes, cx).unwrap();
        let mut split = EnclosedIlluminanceJob::new(&scene, &probes, cx).unwrap();
        full.advance(cx, usize::MAX).unwrap();
        split.advance(cx, 1).unwrap();
        gate.request();
        assert!(split.advance(cx, usize::MAX).is_err());
        assert_eq!(split.results(), &full.results()[..1]);
        with_cx(|resume, _| {
            split.advance(resume, 2).unwrap(); split.advance(resume, 1).unwrap();
            assert_eq!(split.results(), full.results());
            let mut changed = asymmetric();
            if let crate::photometry::source::AngularLaw::Table(table) = &mut changed.law {
                // Preserve poles, seam, peak and integrated flux; only swap
                // two interior lobes. These bytes must participate in identity.
                std::sync::Arc::get_mut(table).unwrap().values.swap(6, 8);
            }
            let other_sources = [source(5500.0 * PI, 3.0, changed)];
            let other = PhotometricScene::try_new(&[], &other_sources, HomogeneousMedium::vacuum(), TransportOptions::default(), resume).unwrap();
            let result = other.illuminance_with_enclosure(resume, &probes[0]).unwrap();
            assert_ne!(result.numerical_model_identity, full.results()[0].numerical_model_identity);
            assert!(result.estimate.lux < full.results()[0].estimate.lux);
        });
    });
}

#[test]
fn dense_axisymmetric_tables_are_admitted_without_a_sixty_four_row_limit() {
    with_cx(|cx, _| {
        let law = AngularDistribution::tabulated(
            (0..=128).map(|i| PI * f64::from(i) / 128.0).collect(),
            vec![0.0, TAU], vec![1.0; 258],
        ).unwrap();
        let sources = [source(1000.0 * PI, 5.0, law)];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let report = scene.illuminance_with_enclosure(cx, &probe(0.0, 0.0)).unwrap();
        contains(report, 10.0);
        assert!(report.numerical_error_bound_lux.unwrap() < 1e-9);
    });
}

#[test]
fn asymmetric_sources_share_native_colored_haze_blockers_and_receiving_floor() {
    use crate::charts::TriMesh;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, Primitive, Shape};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 4.0];
        let grid = AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1, 1, 3], &density, [-4.0, -4.0, 0.0], [8.0, 8.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(), None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 128 }).unwrap();
        let sources = [source(5500.0 * PI, 3.0, asymmetric())];
        let material = Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() };
        let geometry = [Primitive { material, emission: None, shape: Shape::Mesh(TriMesh::new(
            vec![[-1.0,-1.0,1.5], [1.0,-1.0,1.5], [0.0,1.0,1.5]], vec![[0,1,2]])) },
            Primitive { material, emission: None, shape: Shape::Mesh(TriMesh::new(
            vec![[-4.0,-4.0,0.0], [4.0,-4.0,0.0], [4.0,4.0,0.0], [-4.0,4.0,0.0]], vec![[0,1,2],[0,2,3]])) }];
        let scene = PhotometricScene::try_with_medium(&geometry, &sources, PhotometricMedium::Absorbing(&grid), TransportOptions::default(), cx).unwrap();
        contains(scene.illuminance_with_enclosure(cx, &probe(0.0, 0.0)).unwrap(), 0.0);
        let report = scene.illuminance_with_enclosure(cx, &probe(0.0, -3.0)).unwrap();
        let r = 18.0_f64.sqrt();
        let incident = 1000.0 * (3.0 - 6.0/r) * 3.0 / (18.0*r);
        let mut transfer = 0.0;
        for bin in 0..SPECTRAL_BINS {
            transfer += 683.0 * sources[0].spectrum().bins_w_per_lumen_nm()[bin] * photopic_weight_nm(bin)
                * (-7.0 * 2.0_f64.sqrt() * (wavelength_nm(bin)-380.0)/400.0).exp();
        }
        contains(report, incident * transfer);
        assert_eq!(report.uncertain_visibility_sources, Some(0));
        assert!(report.numerical_error_bound_lux.unwrap() < 1e-7);
    });
}

#[test]
fn fixture_roll_rotates_the_asymmetric_beam_not_its_calibration() {
    use crate::instances::RigidTransform;
    use crate::photometry::source::{PhotometricSource, SourceSpectrum};
    with_cx(|cx, _| {
        let q = 0.5_f64.sqrt();
        let original = [source(5500.0 * PI, 3.0, asymmetric())];
        let rolled = [PhotometricSource::try_new(
            RigidTransform::try_new([q, q, 0.0, 0.0], [0.0, 0.0, 3.0]).unwrap(),
            5500.0 * PI, 1.0, asymmetric(), SourceSpectrum::equal_energy()).unwrap()];
        let a = PhotometricScene::try_new(&[], &original, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let b = PhotometricScene::try_new(&[], &rolled, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let before = a.illuminance_with_enclosure(cx, &probe(0.0, -3.0)).unwrap();
        let after = b.illuminance_with_enclosure(cx, &probe(3.0, 0.0)).unwrap();
        let r = 18.0_f64.sqrt();
        let reference = 1000.0 * (3.0 - 6.0/r) * 3.0 / (18.0*r);
        contains(before, reference); contains(after, reference);
        assert!(after.estimate.lux > b.illuminance(cx, &probe(0.0, -3.0)).unwrap().lux);
        assert_eq!(original[0].flux_lm(), rolled[0].flux_lm());
        assert_ne!(before.numerical_model_identity, after.numerical_model_identity);
    });
}

#[test]
fn azimuth_knot_changes_are_bound_even_for_a_flat_intensity_profile() {
    with_cx(|cx, _| {
        let mut previous = None;
        for middle in [FRAC_PI_2, PI] {
            let law = AngularDistribution::tabulated(
                vec![0.0, FRAC_PI_2, PI], vec![0.0, middle, TAU], vec![1.0; 9],
            ).unwrap();
            let sources = [source(1000.0 * PI, 5.0, law)];
            let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx, &probe(0.0, 0.0)).unwrap();
            contains(report, 10.0);
            if let Some(identity) = previous { assert_ne!(report.numerical_model_identity, Some(identity)); }
            previous = report.numerical_model_identity;
        }
    });
}
