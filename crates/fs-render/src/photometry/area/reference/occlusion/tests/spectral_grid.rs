//! Independent spectral layer integrals through the original density/material maps.
use super::*;
use crate::photometry::source::wavelength_nm;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, CellPalette,
    CellPaletteBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction,
    ScatteringGrid, ScatteringOptics};
use std::sync::Arc;

fn medium<'a>(cx: &Cx<'_>, density: &'a [f64]) -> ScatteringGrid<'a> {
    let limits = CellPaletteBudget { max_models: 2, max_cells: 4 };
    let absorption = CellPalette::try_new(cx, vec![
        MassAbsorptionSpectrum::constant(0.125).unwrap(),
        MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(),
    ], Arc::from([0, 1, 1, 0]), limits).unwrap();
    let absorption = AbsorbingGrid::try_new_indexed(cx,
        VolumeGrid::new([2, 1, 2], density, [-2.0, -2.0, 0.5], [2.0, 4.0, 0.5]),
        absorption, None,
        AbsorbingGridBudget { max_grid_cells: 4, max_steps_per_segment: 32 }).unwrap();
    let scattering = CellPalette::try_new(cx, vec![
        ScatteringOptics { scattering: MassScatteringSpectrum::rayleigh(400.0, 0.125).unwrap(),
            phase: PhaseFunction::Rayleigh },
        ScatteringOptics { scattering: MassScatteringSpectrum::constant(0.25).unwrap(),
            phase: PhaseFunction::Isotropic },
    ], Arc::from([1, 0, 1, 0]), limits).unwrap();
    ScatteringGrid::try_new_indexed(cx, absorption, scattering).unwrap()
}

fn spectral_source(blue_fraction: f64) -> PhotometricAreaEmitter {
    let mut bins = [0.0; SPECTRAL_BINS];
    // Desired relative PHOTOPIC weights at 402.5 and 602.5 nm. These are
    // different from equal radiant powers. SourceSpectrum remains authoritative.
    bins[4] = blue_fraction / observer::xyz_at_bin(4)[1];
    bins[44] = (1.0-blue_fraction) / observer::xyz_at_bin(44)[1];
    PhotometricAreaEmitter::try_new(0, 1000.0, 1.0, SourceSpectrum::try_from_bins(bins).unwrap(),
        EmissionSidedness::Front).unwrap()
}

fn target(source: &PhotometricAreaEmitter, shadow: bool) -> f64 {
    // Independently integrated at 70 digits with 24/32-point surface rules:
    // integral 4/r^4 * exp(-r*k_side(lambda)) dx dy, L=1000/(4*pi).
    // The left/right optical depths are found from the two physical layers:
    // k_left=((a0+s1)+2*(a1+s1))/4;
    // k_right=(4*(a1+s0)+8*(a0+s0))/4. The strip removes x in [-.5,.5].
    // No projected-triangle formula, native traversal, or midpoint transmittance
    // supplies this oracle. Retained spectral normalization is applied below.
    let values = if shadow { [10.819_576_768_452_719, 6.415_476_278_584_347] }
        else { [24.759_454_271_862_572, 14.911_501_983_981_905] };
    [4, 44].into_iter().zip(values).map(|(bin, value)| value
        * source.spectrum().bins_w_per_lumen_nm()[bin] * observer::xyz_at_bin(bin)[1]
        * BIN_WIDTH_NM * LUMINOUS_EFFICACY_LM_W).sum()
}

#[test]
fn clear_panel_through_independent_optical_maps_encloses_two_spectral_references() {
    with_cx(|cx, _| {
        let density = [1.0, 4.0, 2.0, 8.0];
        let grid = medium(cx, &density);
        let primitives = [panel(-1.0, 1.0, -1.0, 1.0, 2.0), panel(-4.0, 4.0, -4.0, 4.0, 0.0)];
        let mut values = Vec::new();
        for fraction in [1.0, 0.0] {
            let sources = [spectral_source(fraction)];
            let reference = OccludedAreaReference::new_with_medium(&primitives, &sources,
                PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
                GridEnclosureBudget::default(), cx).unwrap();
            let budget = AreaOcclusionBudget { absolute_error_lux: 2.0, max_patch_evaluations: 4096,
                ..AreaOcclusionBudget::default() };
            let coarse = reference.illuminance(cx, &probe(0.0, 0.0),
                AreaOcclusionBudget { max_patch_evaluations: 2, ..budget }).unwrap();
            let fine = reference.illuminance(cx, &probe(0.0, 0.0), budget).unwrap();
            check(coarse, target(&sources[0], false)); check(fine, target(&sources[0], false));
            assert_eq!(coarse.termination, AreaOcclusionTermination::PatchBudget);
            assert_eq!(fine.termination, AreaOcclusionTermination::ToleranceMet);
            assert!(fine.numerical.lo >= coarse.numerical.lo && fine.numerical.hi <= coarse.numerical.hi);
            assert!(fine.varying_transmittance_patches > 0);
            assert_eq!(fine.absorption_per_m, None); assert_eq!(fine.scattering_per_m, None);
            assert_eq!(fine.physical_grid_identity, Some(grid.content_identity()));
            assert_eq!(fine.numerical_model_identity, coarse.numerical_model_identity);
            values.push(fine);
        }
        assert!(values[0].numerical.lo > values[1].numerical.hi);
        assert_eq!(wavelength_nm(4), 402.5); assert_eq!(wavelength_nm(44), 602.5);
    });
}

#[test]
fn partial_shadow_and_spectral_extinction_share_one_error_budget() {
    with_cx(|cx, _| {
        let density = [1.0, 4.0, 2.0, 8.0]; let grid = medium(cx, &density);
        let primitives = strip_scene(); let sources = [spectral_source(0.75)];
        let reference = OccludedAreaReference::new_with_medium(&primitives, &sources,
            PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
            GridEnclosureBudget::default(), cx).unwrap();
        let budget = AreaOcclusionBudget { absolute_error_lux: 2.0, max_patch_evaluations: 4096,
            ..AreaOcclusionBudget::default() };
        let coarse = reference.illuminance(cx, &probe(0.0, 0.0),
            AreaOcclusionBudget { max_patch_evaluations: 2, ..budget }).unwrap();
        let fine = reference.illuminance(cx, &probe(0.0, 0.0), budget).unwrap();
        check(coarse, target(&sources[0], true)); check(fine, target(&sources[0], true));
        assert_eq!(coarse.termination, AreaOcclusionTermination::PatchBudget);
        assert_eq!(fine.termination, AreaOcclusionTermination::ToleranceMet);
        assert!(fine.varying_transmittance_patches > 0 && fine.uncertain_visibility_patches > 0);
        assert!(fine.evaluated_patches > 2 && fine.evaluated_patches <= 4096);
    });
}

#[test]
fn empty_grid_keeps_vacuum_arithmetic_and_homogeneous_replacement_releases_grid_state() {
    with_cx(|cx, _| {
        let density = [0.0; 4]; let grid = medium(cx, &density);
        let primitives = strip_scene(); let sources = [spectral_source(0.75)];
        let vacuum = OccludedAreaReference::new(&primitives, &sources, AreaReferenceBudget::default(), cx).unwrap();
        let physical = OccludedAreaReference::new_with_medium(&primitives, &sources,
            PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
            GridEnclosureBudget::default(), cx).unwrap();
        let budget = AreaOcclusionBudget { max_patch_evaluations: 128, ..AreaOcclusionBudget::default() };
        let original = vacuum.illuminance(cx, &probe(0.0, 0.0), budget).unwrap();
        let transparent = physical.illuminance(cx, &probe(0.0, 0.0), budget).unwrap();
        assert_eq!(original.numerical, transparent.numerical);
        assert_eq!(original.lux.to_bits(), transparent.lux.to_bits());
        assert_eq!(original.evaluated_patches, transparent.evaluated_patches);
        assert_eq!(transparent.varying_transmittance_patches, 0);
        assert_ne!(original.numerical_model_identity, transparent.numerical_model_identity);
        let restored = physical.with_homogeneous_medium(&HomogeneousMedium::vacuum());
        assert_eq!(original, restored.illuminance(cx, &probe(0.0, 0.0), budget).unwrap());
    });
}

#[test]
fn progressive_grid_results_keep_snapshot_identity_csv_scope_and_cancel_resume() {
    with_cx(|cx, gate| {
        let density = [1.0, 4.0, 2.0, 8.0]; let grid = medium(cx, &density);
        let primitives = strip_scene(); let sources = [spectral_source(0.75)];
        let reference = OccludedAreaReference::new_with_medium(&primitives, &sources,
            PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
            GridEnclosureBudget::default(), cx).unwrap();
        let probes = [probe(0.0, 0.0), probe(0.25, 0.0), probe(-0.25, 0.0)];
        let budget = AreaOcclusionBudget { max_patch_evaluations: 128, ..AreaOcclusionBudget::default() };
        let mut full = OccludedAreaJob::new(&reference, &probes, budget, cx).unwrap();
        full.advance(cx, usize::MAX).unwrap();
        let mut split = OccludedAreaJob::new(&reference, &probes, budget, cx).unwrap();
        split.advance(cx, 1).unwrap();
        gate.request(); assert!(matches!(split.advance(cx, 1), Err(AreaOcclusionError::Cancelled)));
        assert_eq!(split.results(), &full.results()[..1]);
        with_cx(|resume, _| { split.advance(resume, usize::MAX).unwrap(); });
        assert_eq!(split.results(), full.results());
        let report = split.results()[0];
        let mut bytes = Vec::new(); OccludedAreaIlluminance::write_csv_header(&mut bytes).unwrap();
        report.write_csv_row(0, &mut bytes).unwrap();
        let csv = String::from_utf8(bytes).unwrap(); let lines: Vec<_> = csv.lines().collect();
        let fields: Vec<_> = lines[1].split(',').collect();
        assert_eq!(lines[0].split(',').count(), 18); assert_eq!(fields.len(), 18);
        assert_eq!(fields[1], "opaque-mesh-spectral-grid-direct-area");
        assert_eq!(fields[12], "not-applicable"); assert_eq!(fields[13], "not-applicable");
        assert_eq!(fields[15], "unknown"); assert_eq!(fields[16], "3");
        assert_eq!(fields[17], grid.content_identity().to_string());
        assert_eq!(fields[2].parse::<f64>().unwrap().to_bits(), report.lux.to_bits());
    });
}

#[test]
fn physical_grid_refusal_never_publishes_a_partial_receiver() {
    with_cx(|cx, _| {
        let density = [1.0, 4.0, 2.0, 8.0]; let grid = medium(cx, &density);
        let primitives = [panel(-1.0, 1.0, -1.0, 1.0, 2.0)]; let sources = [spectral_source(0.75)];
        let identity = grid.content_identity();
        assert!(matches!(OccludedAreaReference::new_with_medium(&primitives, &sources,
            PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
            GridEnclosureBudget { max_nodes_per_segment: 0 }, cx), Err(AreaOcclusionError::InvalidOptions)));
        let reference = OccludedAreaReference::new_with_medium(&primitives, &sources,
            PhotometricMedium::Scattering(&grid), AreaReferenceBudget::default(),
            GridEnclosureBudget { max_nodes_per_segment: 1 }, cx).unwrap();
        let probes = [probe(0.0, 0.0)];
        let mut job = OccludedAreaJob::new(&reference, &probes, AreaOcclusionBudget::default(), cx).unwrap();
        assert!(matches!(job.advance(cx, 1), Err(AreaOcclusionError::Medium(VolumeTransportError::Budget))));
        assert!(job.results().is_empty()); assert!(!job.is_complete());
        assert_eq!(grid.content_identity(), identity);
    });
}

#[test]
fn independently_prepared_emitters_keep_their_own_spectral_weights() {
    with_cx(|cx, _| {
        let density = [1.0; 4]; let grid = medium(cx, &density);
        let blue = spectral_source(1.0); let red = spectral_source(0.0);
        let sources = [
            PhotometricAreaEmitter::try_new(4, blue.flux_lm(), 1.0, blue.spectrum().clone(), EmissionSidedness::Front).unwrap(),
            PhotometricAreaEmitter::try_new(2, red.flux_lm(), 1.0, red.spectrum().clone(), EmissionSidedness::Front).unwrap(),
        ];
        let prepared = PreparedAreaGrid::new(PhotometricMedium::Scattering(&grid), &sources,
            GridEnclosureBudget::default(), cx).unwrap();
        let a = [point(-0.5), point(0.0), point(2.0)];
        let b = [point(0.0); 3];
        let blue_t = prepared.photopic_transfer(cx, 4, a, b).unwrap();
        let red_t = prepared.photopic_transfer(cx, 2, a, b).unwrap();
        let length = 4.25_f64.sqrt();
        // Both segments stay in the LEFT cells: rho=1 in both layers.
        for (lambda, actual) in [(402.5_f64, blue_t), (602.5_f64, red_t)] {
            let tau = length / 4.0 * (0.125 + (lambda-380.0)/400.0 + 0.5);
            assert!(actual.contains((-tau).exp()));
            assert!(actual.width() < 1e-8);
        }
        assert!(blue_t.lo() > red_t.hi());
        assert!(matches!(prepared.photopic_transfer(cx, 3, a, b),
            Err(AreaOcclusionError::Reference(AreaReferenceError::InvalidGeometry))));
        let changed = [1.0, 4.0, 2.0, 8.0]; let other = medium(cx, &changed);
        let next = PreparedAreaGrid::new(PhotometricMedium::Scattering(&other), &sources,
            GridEnclosureBudget::default(), cx).unwrap();
        assert_ne!(prepared.identity(), next.identity());
    });
}
