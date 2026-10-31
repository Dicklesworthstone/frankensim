use super::*;
use crate::instances::RigidTransform;
use crate::photometry::source::{AngularDistribution, SourceSpectrum};
use crate::photometry::source::validation::{EnclosedIlluminanceJob, EnclosureUnavailable};
use crate::photometry::{IlluminanceProbe, PhotometricScene, TransportOptions};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x5350_4543, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}

fn fixture(spectrum: SourceSpectrum) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0,0.0,0.0,1.0], [0.0,0.0,5.0]).unwrap(),
        1000.0 * std::f64::consts::PI, 1.0, AngularDistribution::isotropic(), spectrum).unwrap()
}

fn probe(x: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x,0.0,0.0), Vec3::new(0.0,0.0,1.0)).unwrap()
}

fn medium<'a>(cx: &Cx<'_>, density: &'a [f64;3]) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(cx, VolumeGrid::new([1,1,3], density, [-0.5,-0.5,1.0], [1.0;3]),
        MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)]).unwrap(), None,
        AbsorbingGridBudget {max_grid_cells:3,max_steps_per_segment:128}).unwrap()
}

#[test]
fn layered_colored_attenuation_matches_an_independent_discrete_spectral_integral() {
    with_cx(|cx,_| {
        let density = [1.0,2.0,4.0];
        let grid = medium(cx,&density);
        let mut values = Vec::new();
        let mut identities = Vec::new();
        for bin in [0, SPECTRAL_BINS-1] {
            let mut bins = [0.0;SPECTRAL_BINS]; bins[bin] = 1.0;
            let sources = [fixture(SourceSpectrum::try_from_bins(bins).unwrap())];
            let scene = PhotometricScene::try_with_medium(&[], &sources, PhotometricMedium::Absorbing(&grid),
                TransportOptions::default(), cx).unwrap();
            let p = probe(0.0);
            let report = scene.illuminance_with_enclosure(cx,&p).unwrap();
            // Independent integral: column mass 1+2+4=7; affine k=(lambda-380)/400.
            // Retained source/observer weights are inputs, not an assumed exact
            // normalization or a claim about continuous CIE metrology.
            let mass = 683.0 * sources[0].spectrum().bins_w_per_lumen_nm()[bin] * photopic_weight_nm(bin);
            let reference = 10.0 * mass * (-7.0 * (wavelength_nm(bin)-380.0)/400.0).exp();
            assert_eq!(report.unavailable,None);
            assert_eq!(report.estimate.numerical.kind,NumericalKind::Enclosure);
            assert!(report.estimate.numerical.lo <= reference && report.estimate.numerical.hi >= reference);
            assert!(report.numerical_error_bound_lux.unwrap() < 1e-9);
            assert_eq!(report.estimate.lux.to_bits(),scene.illuminance(cx,&p).unwrap().lux.to_bits());
            assert!(report.estimate.total_error_bound_lux.is_none());
            values.push(report.estimate.lux); identities.push(report.numerical_model_identity);
            let mut csv = Vec::new(); report.write_csv_row(0,&mut csv).unwrap();
            let csv = String::from_utf8(csv).unwrap();
            assert!(csv.contains(",unknown,unknown,unknown,unknown,3,0"));
        }
        assert!(values[0] > 900.0 * values[1], "colored attenuation was incorrectly collapsed to gray");
        assert_ne!(identities[0],identities[1]);
    });
}

#[test]
fn native_scattering_grid_phase_does_not_change_direct_extinction() {
    with_cx(|cx,_| {
        let density = [1.0,2.0,4.0];
        let sources = [fixture(SourceSpectrum::equal_energy())];
        let mut previous = None;
        for phase in [PhaseFunction::Isotropic, PhaseFunction::HenyeyGreenstein {g:0.7}] {
            let grid = ScatteringGrid::try_new(cx, medium(cx,&density),
                MassScatteringSpectrum::rayleigh(550.0,0.5).unwrap(),phase).unwrap();
            let scene = PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),
                TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.0)).unwrap();
            let mut reference = 0.0;
            for bin in 0..SPECTRAL_BINS {
                let wavelength = wavelength_nm(bin);
                let ratio = 550.0/wavelength;
                let extinction = (wavelength-380.0)/400.0 + 0.5*ratio*ratio*ratio*ratio;
                reference += sources[0].spectrum().bins_w_per_lumen_nm()[bin]
                    * photopic_weight_nm(bin) * (-7.0*extinction).exp();
            }
            reference *= 10.0*683.0;
            assert_eq!(report.unavailable,None);
            assert!(report.estimate.numerical.lo <= reference && report.estimate.numerical.hi >= reference);
            assert!(report.estimate.total_error_bound_lux.is_none());
            if let Some(previous) = previous { assert_eq!(report.estimate,previous); }
            previous = Some(report.estimate);
        }
    });
}

#[test]
fn optional_medium_budget_failure_keeps_nominal_lux_and_cannot_change_a_published_job() {
    with_cx(|cx,_| {
        let density = [1.0,2.0,4.0]; let grid = medium(cx,&density);
        let sources = [fixture(SourceSpectrum::equal_energy())]; let probes = [probe(0.0)];
        let scene = PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),cx).unwrap();
        let mut job = EnclosedIlluminanceJob::new(&scene,&probes,cx).unwrap();
        job.set_grid_enclosure_budget(GridEnclosureBudget {max_nodes_per_segment:1}).unwrap();
        job.advance(cx,1).unwrap();
        let report = job.results()[0];
        assert_eq!(report.unavailable,Some(EnclosureUnavailable::MediumBudgetExceeded));
        assert_eq!(report.estimate.numerical.kind,NumericalKind::NoClaim);
        assert!(report.numerical_error_bound_lux.is_none());
        assert_eq!(report.estimate.lux.to_bits(),scene.illuminance(cx,&probes[0]).unwrap().lux.to_bits());
        assert!(job.set_grid_enclosure_budget(GridEnclosureBudget::default()).is_err());
        assert_eq!(job.results(),&[report]);
    });
}

#[test]
fn spectral_jobs_replay_after_cancel_and_bind_density_changes() {
    with_cx(|cx,gate| {
        let density = [1.0,2.0,4.0]; let grid = medium(cx,&density);
        let sources = [fixture(SourceSpectrum::equal_energy())]; let probes = [probe(0.0),probe(0.1),probe(0.2)];
        let scene = PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),cx).unwrap();
        let mut full = EnclosedIlluminanceJob::new(&scene,&probes,cx).unwrap(); full.advance(cx,usize::MAX).unwrap();
        let mut split = EnclosedIlluminanceJob::new(&scene,&probes,cx).unwrap(); split.advance(cx,1).unwrap();
        gate.request(); assert!(split.advance(cx,1).is_err());
        assert_eq!(split.results(),&full.results()[..1]);
        with_cx(|resume,_| {
            split.advance(resume,usize::MAX).unwrap(); assert_eq!(split.results(),full.results());
            let changed_density = [1.0,3.0,4.0]; let changed_grid = medium(resume,&changed_density);
            let changed = PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Absorbing(&changed_grid),
                TransportOptions::default(),resume).unwrap();
            let report = changed.illuminance_with_enclosure(resume,&probes[0]).unwrap();
            assert_ne!(report.numerical_model_identity,full.results()[0].numerical_model_identity);
            assert!(report.estimate.lux < full.results()[0].estimate.lux);
        });
    });
}

#[test]
fn indexed_spectral_haze_and_native_shadow_geometry_share_the_same_result() {
    use crate::charts::TriMesh;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, Primitive, Shape};
    use crate::volumes::transport::{CellPalette, CellPaletteBudget};
    use crate::volumes::transport::scattering::ScatteringOptics;
    use std::sync::Arc;
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 4.0];
        let budget = CellPaletteBudget { max_models: 2, max_cells: 3 };
        let absorption = CellPalette::try_new(cx, vec![
            MassAbsorptionSpectrum::constant(0.25).unwrap(),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(),
        ], Arc::from([0, 1, 0]), budget).unwrap();
        let grid = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([1, 1, 3], &density, [-0.5, -0.5, 1.0], [1.0; 3]),
            absorption, None, AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 128 }).unwrap();
        let scattering = CellPalette::try_new(cx, vec![
            ScatteringOptics { scattering: MassScatteringSpectrum::rayleigh(400.0, 0.5).unwrap(), phase: PhaseFunction::Rayleigh },
            ScatteringOptics { scattering: MassScatteringSpectrum::constant(1.0).unwrap(), phase: PhaseFunction::Isotropic },
        ], Arc::from([1, 0, 1]), budget).unwrap();
        let grid = ScatteringGrid::try_new_indexed(cx, grid, scattering).unwrap();
        let sources = [fixture(SourceSpectrum::equal_energy())];
        let material = Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() };
        let geometry = [
            Primitive { shape: Shape::Mesh(TriMesh::new(
                vec![[-0.1, -0.1, 2.5], [0.1, -0.1, 2.5], [0.0, 0.1, 2.5]], vec![[0, 1, 2]])), material, emission: None },
            Primitive { shape: Shape::Mesh(TriMesh::new(
                vec![[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [0.0, 2.0, 0.0]], vec![[0, 1, 2]])), material, emission: None },
        ];
        let scene = PhotometricScene::try_with_medium(&geometry, &sources, PhotometricMedium::Scattering(&grid),
            TransportOptions::default(), cx).unwrap();
        let probes = [probe(0.0), probe(0.4)];
        let mut job = EnclosedIlluminanceJob::new(&scene, &probes, cx).unwrap();
        job.advance(cx, 1).unwrap();
        assert_eq!(job.results()[0].estimate.numerical.lo, 0.0);
        assert_eq!(job.results()[0].estimate.numerical.hi, 0.0);
        job.advance(cx, 1).unwrap();
        let report = job.results()[1];
        let r = (25.0_f64 + 0.4 * 0.4).sqrt();
        let mut expected = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let lambda = wavelength_nm(bin);
            let ratio = 400.0 / lambda;
            // Independent layer sum. The first/third cells use (a,s)=(1/4,1);
            // only the middle cell has affine absorption and molecular scattering.
            let depth = (6.25 + 2.0 * ((lambda - 380.0) / 400.0 + 0.5 * ratio * ratio * ratio * ratio)) * r / 5.0;
            expected += sources[0].spectrum().bins_w_per_lumen_nm()[bin] * photopic_weight_nm(bin) * (-depth).exp();
        }
        expected *= 683.0 * 1250.0 / (r * r * r);
        assert_eq!(report.unavailable, None);
        assert_eq!(report.uncertain_visibility_sources, Some(0));
        assert!(report.estimate.numerical.lo <= expected && report.estimate.numerical.hi >= expected);
        assert!(report.numerical_error_bound_lux.unwrap() < 1e-9);
        assert_eq!(report.estimate.lux.to_bits(), scene.illuminance(cx, &probes[1]).unwrap().lux.to_bits());
        assert!(report.estimate.total_error_bound_lux.is_none());
    });
}
