use super::*;
use super::super::tests::{probe, source, with_cx};
use crate::photometry::medium::HomogeneousMedium;
use crate::photometry::source::{AngularDistribution, SourceSpectrum, photopic_weight_nm, wavelength_nm};
use crate::photometry::transport::PhotometricMedium;
use crate::photometry::TransportOptions;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_geom::{Point3, Vec3};
use std::f64::consts::PI;

// These are authored assumptions, NOT measured fixture calibration records.
fn bound(lo: f64, hi: f64) -> Option<FluxScaleBound> {
    Some(FluxScaleBound::try_new(lo, hi, ContentHash([1; 32])).unwrap())
}

fn run(scene: &PhotometricScene<'_>, probes: &[IlluminanceProbe], bounds: &[Option<FluxScaleBound>], cx: &Cx<'_>) -> Vec<FluxBoundedIlluminance> {
    let mut job = FluxEnclosureJob::new(scene, probes, bounds, cx).unwrap();
    job.advance(cx, usize::MAX).unwrap();
    job.results().to_vec()
}

#[test]
fn declarations_reject_invalid_ranges_and_unnamed_evidence() {
    for (lo, hi) in [(-0.1, 1.0), (2.0, 1.0), (0.0, f64::INFINITY), (f64::NAN, 1.0), (0.0, f64::NAN)] {
        assert_eq!(FluxScaleBound::try_new(lo, hi, ContentHash([1; 32])), Err(FluxBoundError::InvalidRange));
    }
    assert_eq!(FluxScaleBound::try_new(1.0, 1.0, ContentHash([0; 32])), Err(FluxBoundError::MissingEvidence));
    let zero = FluxScaleBound::try_new(-0.0, -0.0, ContentHash([2; 32])).unwrap();
    assert_eq!(zero.scale().lo().to_bits(), 0);
    assert_eq!(zero.evidence_identity(), ContentHash([2; 32]));
}

#[test]
fn independent_two_source_extrema_preserve_nominal_and_total_error_scope() {
    with_cx(|cx, _| {
        // Independent inverse-square references: 100 lx and 50 lx.
        let sources = [source(400.0 * PI, 1.0, AngularDistribution::isotropic()),
            source(800.0 * PI, 2.0, AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let probes = [probe(0.0, 0.0)];
        let report = run(&scene, &probes, &[bound(0.9, 1.1), bound(0.8, 1.05)], cx)[0];
        let range = report.flux_conditioned_lux.unwrap();
        assert!(range.contains(130.0) && range.contains(162.5));
        assert!(range.width() < 32.5 + 1e-9);
        let delta = report.flux_perturbation_bound_lux.unwrap();
        assert!(delta >= 20.0 && delta < 20.0 + 1e-9);
        assert_eq!(report.nominal, scene.illuminance_with_enclosure(cx, &probes[0]).unwrap());
        assert_eq!(report, scene.illuminance_with_flux_bounds(cx, &probes[0], &[bound(0.9, 1.1), bound(0.8, 1.05)]).unwrap());
        assert_eq!(report.unbounded_flux_sources, Some(0));
        assert!(report.numerical_and_flux_error_bound_lux().unwrap() < 20.0 + 1e-9);
        assert!(report.nominal.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn shared_calibration_error_is_not_divided_by_sqrt_source_count() {
    with_cx(|cx, _| {
        let sources = [source(400.0 * PI, 1.0, AngularDistribution::isotropic()),
            source(800.0 * PI, 2.0, AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        // Same evidence ID: correlated common gain is explicitly permitted.
        let report = run(&scene, &[probe(0.0, 0.0)], &[bound(0.9, 1.1); 2], cx)[0];
        let range = report.flux_conditioned_lux.unwrap();
        assert!(range.contains(135.0) && range.contains(165.0));
        assert!(report.flux_perturbation_bound_lux.unwrap() >= 15.0);
        assert!(report.flux_perturbation_bound_lux.unwrap() < 15.0 + 1e-9);
    });
}

#[test]
fn bias_missing_bounds_and_explicit_zero_are_different_results() {
    with_cx(|cx, _| {
        let sources = [source(400.0 * PI, 1.0, AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let probes = [probe(0.0, 0.0)];
        let biased = run(&scene, &probes, &[bound(1.25, 1.5)], cx)[0];
        assert!(biased.flux_conditioned_lux.unwrap().contains(125.0));
        assert!(biased.flux_conditioned_lux.unwrap().contains(150.0));
        assert!(!biased.flux_conditioned_lux.unwrap().contains(biased.nominal.estimate.lux));
        let missing = run(&scene, &probes, &[None], cx)[0];
        assert_eq!(missing.unbounded_flux_sources, Some(1));
        assert!(missing.flux_conditioned_lux.is_none());
        assert!(missing.flux_perturbation_bound_lux.is_none());
        assert!(missing.numerical_and_flux_error_bound_lux().is_none());
        let zero = run(&scene, &probes, &[bound(0.0, 0.0)], cx)[0];
        assert_eq!(zero.flux_conditioned_lux, Some(Interval::point(0.0)));
        assert!(zero.flux_perturbation_bound_lux.unwrap() >= 100.0);
        let fixed = run(&scene, &probes, &[bound(1.0, 1.0)], cx)[0];
        assert_eq!(fixed.flux_perturbation_bound_lux, Some(0.0));
        assert_eq!(fixed.nominal, biased.nominal);
        assert_eq!(fixed.nominal, missing.nominal);
        assert_eq!(fixed.nominal, zero.nominal);
        assert!(zero.nominal.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn unknown_provably_blocked_or_exactly_off_sources_do_not_invent_light() {
    use crate::charts::TriMesh;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, Primitive, Shape};
    with_cx(|cx, _| {
        let sources = [source(400.0 * PI, 3.0, AngularDistribution::isotropic()),
            source(0.0, 4.0, AngularDistribution::isotropic())];
        let geometry = [Primitive {
            shape: Shape::Mesh(TriMesh::new(vec![[-1.0,-1.0,1.0],[1.0,-1.0,1.0],[0.0,1.0,1.0]], vec![[0,1,2]])),
            material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() }, emission: None,
        }];
        let scene = PhotometricScene::try_new(&geometry, &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let report = run(&scene, &[probe(0.0, 0.0)], &[None; 2], cx)[0];
        assert_eq!(report.unbounded_flux_sources, Some(0));
        assert_eq!(report.flux_conditioned_lux, Some(Interval::point(0.0)));
        assert_eq!(report.flux_perturbation_bound_lux, Some(0.0));
        assert_eq!(report.numerical_and_flux_error_bound_lux(), Some(0.0));
        assert!(report.nominal.estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn spectral_haze_uses_actual_source_bins_and_budget_failure_discards_partial_accounts() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 4.0];
        let grid = AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,3], &density, [-0.5,-0.5,1.0], [1.0;3]),
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)]).unwrap(), None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 128 }).unwrap();
        for bin in [0, 79] {
            let mut bins = [0.0;80]; bins[bin] = 1.0;
            let mut fixture = source(1000.0 * PI, 5.0, AngularDistribution::isotropic());
            // Construct with the original placement/law and explicit single-bin spectrum.
            fixture = crate::photometry::source::PhotometricSource::try_new(
                fixture.transform(), 1000.0 * PI, 1.0, AngularDistribution::isotropic(),
                SourceSpectrum::try_from_bins(bins).unwrap()).unwrap();
            let sources = [fixture];
            let scene = PhotometricScene::try_with_medium(&[], &sources, PhotometricMedium::Absorbing(&grid), TransportOptions::default(), cx).unwrap();
            let probes = [probe(0.0,0.0)]; let bounds = [bound(0.8,1.2)];
            let report = run(&scene, &probes, &bounds, cx)[0];
            let mass = 683.0 * sources[0].spectrum().bins_w_per_lumen_nm()[bin] * photopic_weight_nm(bin);
            let reference = 10.0 * mass * (-7.0 * (wavelength_nm(bin)-380.0)/400.0).exp();
            assert!(report.flux_conditioned_lux.unwrap().contains(0.8 * reference));
            assert!(report.flux_conditioned_lux.unwrap().contains(1.2 * reference));
            assert!(report.flux_perturbation_bound_lux.unwrap() >= 0.2 * reference);
            assert_eq!(report.nominal, scene.illuminance_with_enclosure(cx, &probes[0]).unwrap());
            let mut exhausted = FluxEnclosureJob::new(&scene, &probes, &bounds, cx).unwrap();
            exhausted.set_grid_enclosure_budget(GridEnclosureBudget { max_nodes_per_segment: 1 }).unwrap();
            exhausted.advance(cx, 1).unwrap();
            let failed = exhausted.results()[0];
            assert_eq!(failed.nominal.unavailable, Some(super::super::EnclosureUnavailable::MediumBudgetExceeded));
            assert!(failed.flux_conditioned_lux.is_none());
            assert!(failed.flux_perturbation_bound_lux.is_none());
            assert_eq!(failed.unbounded_flux_sources, None);
            assert!(failed.combined_identity.is_none());
            assert_eq!(failed.nominal.estimate.lux.to_bits(), report.nominal.estimate.lux.to_bits());
            assert!(exhausted.set_grid_enclosure_budget(GridEnclosureBudget::default()).is_err());
        }
    });
}

#[test]
fn grid_cancel_resume_replays_and_declarations_bind_evidence_and_order() {
    with_cx(|cx, gate| {
        let sources = [source(400.0 * PI, 3.0, AngularDistribution::isotropic()),
            source(800.0 * PI, 4.0, AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        let grid = CalculationGrid::try_new(Point3::new(-0.5,-0.5,0.0), Vec3::new(0.25,0.0,0.0), Vec3::new(0.0,0.25,0.0), 5, 5).unwrap();
        let bounds = [bound(0.8,1.1),bound(0.9,1.2)];
        let mut full = FluxEnclosureJob::for_grid(&scene, grid, &bounds, cx).unwrap();
        full.advance(cx,usize::MAX).unwrap();
        let mut split = FluxEnclosureJob::for_grid(&scene, grid, &bounds, cx).unwrap();
        split.advance(cx,3).unwrap();
        gate.request();
        assert!(matches!(split.advance(cx,4),Err(PhotometryError::Cancelled)));
        assert_eq!(split.results(), &full.results()[..3]);
        with_cx(|resume,_| {
            split.advance(resume,usize::MAX).unwrap();
            assert_eq!(split.results(),full.results());
            let probes = [grid.probe(0).unwrap()];
            let swapped = run(&scene,&probes,&[bounds[1],bounds[0]],resume)[0];
            let evidence = [Some(FluxScaleBound::try_new(0.8,1.1,ContentHash([2;32])).unwrap()),bounds[1]];
            let changed = run(&scene,&probes,&evidence,resume)[0];
            assert_ne!(changed.combined_identity,full.results()[0].combined_identity);
            assert_ne!(swapped.combined_identity,full.results()[0].combined_identity);
            assert_eq!(changed.flux_conditioned_lux,full.results()[0].flux_conditioned_lux);
            assert_eq!(changed.nominal,full.results()[0].nominal);
            assert!(FluxEnclosureJob::new(&scene,&probes,&[bounds[0]],resume).is_err());
        });
    });
}

#[test]
fn accumulator_contains_extreme_ranges_without_fabricating_finite_bounds() {
    let mut underflow = FluxAccount::new();
    underflow.add(Interval::point(f64::from_bits(1)),bound(0.5,0.5));
    assert!(underflow.lux.hi() > 0.0);
    let mut overflow = FluxAccount::new();
    overflow.add(Interval::point(f64::MAX),bound(2.0,2.0));
    assert!(!finite(overflow.lux));
    let mut zero = FluxAccount::new();
    zero.add(Interval::point(f64::MAX),bound(0.0,0.0));
    assert_eq!(zero.lux,Interval::point(0.0));
    let mut fixed = FluxAccount::new();
    fixed.add(Interval::new(0.0,2.0),bound(1.0,1.0));
    assert_eq!(fixed.perturbation,0.0);
    assert!(fixed.lux.contains(0.0) && fixed.lux.contains(2.0));
}

#[test]
fn csv_preserves_unknowns_round_trip_values_and_writer_refusal() {
    with_cx(|cx,_| {
        let sources = [source(400.0 * PI, 3.0, AngularDistribution::isotropic())];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx).unwrap();
        for bounds in [[bound(0.9,1.1)],[None]] {
            let report = run(&scene,&[probe(0.0,0.0)],&bounds,cx)[0];
            let mut bytes = Vec::new();
            FluxBoundedIlluminance::write_csv_header(&mut bytes).unwrap();
            report.write_csv_row(0,&mut bytes).unwrap();
            let csv = String::from_utf8(bytes).unwrap();
            let rows: Vec<_> = csv.lines().collect();
            let fields: Vec<_> = rows[1].split(',').collect();
            assert_eq!(fields.len(),17);
            assert_eq!(rows[0].split(',').count(),17);
            assert_eq!(&fields[13..16],&["unknown";3]);
            assert_eq!(fields[1].parse::<f64>().unwrap().to_bits(),report.nominal.estimate.lux.to_bits());
            if bounds[0].is_none() {
                assert_eq!(&fields[5..9],&["unknown";4]);
                assert_eq!(fields[10],"missing-flux-bounds");
            } else { assert_eq!(fields[10],"conditional-model-bound"); }
            struct Refuse;
            impl Write for Refuse {
                fn write(&mut self,_: &[u8]) -> io::Result<usize> { Err(io::Error::other("test writer")) }
                fn flush(&mut self) -> io::Result<()> { Ok(()) }
            }
            assert!(report.write_csv_row(0,&mut Refuse).is_err());
        }
    });
}
