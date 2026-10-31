use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x4752_4944, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate));
}

fn grid<'a>(cx: &Cx<'_>, dims: [usize; 3], data: &'a [f64], spectrum: MassAbsorptionSpectrum) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(cx, VolumeGrid::new(dims, data, [0.0; 3], [1.0; 3]),
        spectrum, None, AbsorbingGridBudget { max_grid_cells: data.len(), max_steps_per_segment: 256 }).unwrap()
}

fn tight(value: Interval, reference: f64) {
    assert!(value.lo() <= reference && value.hi() >= reference, "{value:?} excludes {reference}");
    assert!(value.width() <= reference.abs().max(1.0) * 2e-11, "{value:?} is unexpectedly loose");
}

#[test]
fn heterogeneous_layers_have_independent_column_mass_and_spectral_attenuation() {
    with_cx(|cx, _| {
        let data = [1.0, 2.0, 4.0];
        let medium = grid(cx, [1, 1, 3], &data,
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (580.0, 0.5), (780.0, 1.0)]).unwrap());
        let a = Point3::new(0.25, 0.25, -1.0);
        let b = Point3::new(0.25, 0.25, 4.0);
        let column = medium.column_density_enclosure(cx, a, b, GridEnclosureBudget::default()).unwrap();
        tight(column.kg_per_m2, 7.0); // three one-metre layers: 1 + 2 + 4
        assert_eq!(column.possible_cells, 3);
        assert_eq!(column.uncertain_parallel_cells, 0);
        let t = medium.transmittance_enclosure(cx, a, b, &[380.0, 580.0, 780.0], GridEnclosureBudget::default()).unwrap();
        assert_eq!(t[0], Interval::point(1.0));
        tight(t[1], 0.030_197_383_422_318_5); // independent exp(-7/2)
        tight(t[2], 0.000_911_881_965_554_516_2); // exp(-7)
        tight(medium.column_density_enclosure(cx, b, a, GridEnclosureBudget::default()).unwrap().kg_per_m2, 7.0);
    });
}

#[test]
fn oblique_segment_lengths_use_metric_distance_and_original_cell_order() {
    with_cx(|cx, _| {
        // Along (1,0,1), the two occupied cells are indices 0 and 3.
        let data = [1.0, 90.0, 80.0, 3.0];
        let medium = grid(cx, [2, 1, 2], &data, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let a = Point3::new(0.0, 0.25, 0.0);
        let b = Point3::new(2.0, 0.25, 2.0);
        let column = medium.column_density_enclosure(cx, a, b, GridEnclosureBudget::default()).unwrap();
        tight(column.kg_per_m2, 5.656_854_249_492_381); // (1+3)*sqrt(2)
        // Endpoint splitting is additive even when the split lies on cell planes.
        let mid = Point3::new(1.0, 0.25, 1.0);
        let left = medium.column_density_enclosure(cx, a, mid, GridEnclosureBudget::default()).unwrap();
        let right = medium.column_density_enclosure(cx, mid, b, GridEnclosureBudget::default()).unwrap();
        assert!((left.kg_per_m2 + right.kg_per_m2).intersect(column.kg_per_m2).is_some());
    });
}

#[test]
fn parallel_internal_faces_do_not_double_count_a_certain_lower_bound() {
    with_cx(|cx, _| {
        let data = [2.0, 5.0];
        let medium = grid(cx, [2, 1, 1], &data, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let a = Point3::new(1.0, 0.0, 0.25);
        let b = Point3::new(1.0, 1.0, 0.25);
        let column = medium.column_density_enclosure(cx, a, b, GridEnclosureBudget::default()).unwrap();
        // Exact half-open ownership selects x=1, density 5. Rounded cell planes
        // need not decide that ownership, but cannot exclude its true integral.
        assert!(column.kg_per_m2.lo() <= 5.0 && column.kg_per_m2.hi() >= 5.0);
        assert_eq!(column.uncertain_parallel_cells, 2);
        assert_eq!(column.kg_per_m2.lo(), 0.0);
    });
}

#[test]
fn vacuum_misses_and_zero_length_are_exact_while_finite_underflow_is_not_zero() {
    with_cx(|cx, _| {
        let data = [1.0];
        let medium = grid(cx, [1; 3], &data, MassAbsorptionSpectrum::constant(1000.0).unwrap());
        let a = Point3::new(0.25, 0.25, 0.0);
        let b = Point3::new(0.25, 0.25, 1.0);
        let t = medium.transmittance_enclosure(cx, a, b, &[550.0], GridEnclosureBudget::default()).unwrap()[0];
        assert!(t.hi() > 0.0 && t.hi() < f64::MIN_POSITIVE);
        assert_eq!(medium.transmittance_enclosure(cx, a, a, &[550.0], GridEnclosureBudget::default()).unwrap(), [Interval::point(1.0)]);
        let a = Point3::new(2.0, 2.0, -1.0);
        let b = Point3::new(2.0, 2.0, 2.0);
        assert_eq!(medium.transmittance_enclosure(cx, a, b, &[550.0], GridEnclosureBudget::default()).unwrap(), [Interval::point(1.0)]);
    });
}

#[test]
fn scattering_extinction_uses_original_absorption_and_rayleigh_power_law() {
    with_cx(|cx, _| {
        let data = [2.0];
        let medium = ScatteringGrid::try_new(cx,
            grid(cx, [1; 3], &data, MassAbsorptionSpectrum::constant(0.25).unwrap()),
            MassScatteringSpectrum::rayleigh(400.0, 0.5).unwrap(), PhaseFunction::Isotropic).unwrap();
        let k = medium.extinction_coefficients_enclosure(&[400.0, 600.0]).unwrap();
        tight(k[0], 0.75);
        tight(k[1], 113.0 / 324.0); // 1/4 + (1/2)*(2/3)^4
        let a = Point3::new(0.25, 0.25, 0.0);
        let b = Point3::new(0.25, 0.25, 1.0);
        let t = medium.transmittance_enclosure(cx, a, b, &[400.0, 600.0], GridEnclosureBudget::default()).unwrap();
        tight(t[0], 0.223_130_160_148_429_82); // exp(-1.5)
        let no_scattering = ScatteringGrid::try_new(cx,
            grid(cx, [1; 3], &data, MassAbsorptionSpectrum::constant(0.25).unwrap()),
            MassScatteringSpectrum::constant(0.0).unwrap(), PhaseFunction::Isotropic).unwrap();
        assert_eq!(no_scattering.transmittance_enclosure(cx,a,b,&[550.0],GridEnclosureBudget::default()).unwrap(),
            no_scattering.absorption_grid().transmittance_enclosure(cx,a,b,&[550.0],GridEnclosureBudget::default()).unwrap());
    });
}

#[test]
fn source_domain_budgets_and_cancellation_refuse_instead_of_returning_a_prefix() {
    with_cx(|cx, gate| {
        let data = [1.0; 8];
        let medium = grid(cx, [2; 3], &data, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let a = Point3::new(0.25, 0.25, -1.0);
        let b = Point3::new(0.25, 0.25, 3.0);
        assert_eq!(medium.extinction_coefficients_enclosure(&[379.0]), Err(VolumeTransportError::Wavelength));
        assert_eq!(medium.extinction_coefficients_enclosure(&[]), Err(VolumeTransportError::Sample));
        assert_eq!(medium.column_density_enclosure(cx,a,b,GridEnclosureBudget {max_nodes_per_segment:1}), Err(VolumeTransportError::Budget));
        let tight_budget = AbsorbingGrid::try_new(cx, VolumeGrid::new([2;3], &data, [0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(), None,
            AbsorbingGridBudget {max_grid_cells:8,max_steps_per_segment:1}).unwrap();
        assert_eq!(tight_budget.column_density_enclosure(cx,a,b,GridEnclosureBudget::default()), Err(VolumeTransportError::Budget));
        gate.request();
        assert_eq!(medium.column_density_enclosure(cx,a,b,GridEnclosureBudget::default()), Err(VolumeTransportError::Cancelled));
    });
}

#[test]
fn hierarchy_prunes_without_a_dense_scan_or_density_copy() {
    with_cx(|cx, _| {
        let data = [1.0; 512];
        let medium = grid(cx, [8; 3], &data, MassAbsorptionSpectrum::constant(0.5).unwrap());
        let column = medium.column_density_enclosure(cx, Point3::new(0.25,0.25,-1.0),
            Point3::new(0.25,0.25,9.0), GridEnclosureBudget::default()).unwrap();
        tight(column.kg_per_m2, 8.0);
        assert_eq!(column.possible_cells, 8);
        assert!(column.visited_nodes < 128);
    });
}

#[test]
fn independent_optical_palettes_integrate_each_cells_actual_spectrum() {
    use crate::volumes::transport::{CellPalette, CellPaletteBudget};
    use crate::volumes::transport::scattering::ScatteringOptics;
    use std::sync::Arc;
    with_cx(|cx, _| {
        let data = [1.0, 2.0, 4.0];
        let palette_budget = CellPaletteBudget { max_models: 2, max_cells: 3 };
        let absorption = CellPalette::try_new(cx, vec![
            MassAbsorptionSpectrum::constant(0.25).unwrap(),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(),
        ], Arc::from([0, 1, 0]), palette_budget).unwrap();
        let grid = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([1, 1, 3], &data, [0.0; 3], [1.0; 3]), absorption, None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 32 }).unwrap();
        let scattering = CellPalette::try_new(cx, vec![
            ScatteringOptics {
                scattering: MassScatteringSpectrum::rayleigh(400.0, 0.5).unwrap(),
                phase: PhaseFunction::Rayleigh,
            },
            ScatteringOptics {
                scattering: MassScatteringSpectrum::constant(1.0).unwrap(),
                phase: PhaseFunction::Isotropic,
            },
        ], Arc::from([1, 0, 1]), palette_budget).unwrap();
        let medium = ScatteringGrid::try_new_indexed(cx, grid, scattering).unwrap();
        let wavelengths = [400.0, 600.0];
        assert_eq!(medium.extinction_coefficients_enclosure(&wavelengths), Err(VolumeTransportError::Palette));
        // Two independent two-model tables at two wavelengths need 8 values,
        // not a four-pair table and not a coefficient packet in every cell.
        assert!(matches!(medium.prepare_extinction_enclosure(cx, &wavelengths, 7), Err(VolumeTransportError::Budget)));
        let prepared = medium.prepare_extinction_enclosure(cx, &wavelengths, 8).unwrap();
        let a = Point3::new(0.25, 0.25, -1.0);
        let b = Point3::new(0.25, 0.25, 4.0);
        let value = prepared.transmittance(cx, a, b, GridEnclosureBudget::default()).unwrap();
        // Independent per-cell formula: (1+4)*(1/4+1) +
        // 2*((lambda-380)/400 + (1/2)*(400/lambda)^4).
        for (lane, lambda) in wavelengths.into_iter().enumerate() {
            let ratio = 400.0 / lambda;
            let tau = 6.25 + 2.0 * ((lambda - 380.0) / 400.0 + 0.5 * ratio * ratio * ratio * ratio);
            tight(value[lane], (-tau).exp());
        }
        assert_eq!(value, medium.transmittance_enclosure(cx, a, b, &wavelengths, GridEnclosureBudget::default()).unwrap());
        let reverse = prepared.transmittance(cx, b, a, GridEnclosureBudget::default()).unwrap();
        for lane in 0..2 { assert!(value[lane].intersect(reverse[lane]).is_some()); }
    });
}

#[test]
fn preparation_refuses_limits_and_cancellation_without_mutating_the_medium() {
    with_cx(|cx, gate| {
        let data = [1.0];
        let medium = grid(cx, [1; 3], &data, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let identity = medium.content_identity();
        assert!(matches!(medium.prepare_extinction_enclosure(cx, &[400.0, 600.0], 1), Err(VolumeTransportError::Budget)));
        assert!(matches!(medium.prepare_extinction_enclosure(cx, &[], 2), Err(VolumeTransportError::Sample)));
        gate.request();
        assert!(matches!(medium.prepare_extinction_enclosure(cx, &[400.0], 1), Err(VolumeTransportError::Cancelled)));
        assert_eq!(medium.content_identity(), identity);
    });
}
