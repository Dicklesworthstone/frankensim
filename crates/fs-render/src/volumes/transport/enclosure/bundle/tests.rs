use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, CellPalette, CellPaletteBudget};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid, ScatteringOptics};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 173, kernel_id: 0x424f5854, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx, &gate)
    })
}
fn boxed(p: Point3) -> [Interval; 3] { [p.x, p.y, p.z].map(Interval::point) }
fn grid<'a>(cx: &Cx<'_>, data: &'a [f64], dims: [usize; 3], max_steps: u32) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(cx, VolumeGrid::new(dims, data, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(0.25).unwrap(), None,
        AbsorbingGridBudget { max_grid_cells: data.len(), max_steps_per_segment: max_steps }).unwrap()
}

#[test]
fn singleton_boxes_preserve_original_column_and_spectral_packet_bits() {
    with_cx(|cx, _| {
        let data = [1.0, 2.0, 4.0, 8.0];
        let medium = grid(cx, &data, [1, 1, 4], 32);
        let packet = medium.prepare_extinction_enclosure(cx, &[400.0, 600.0], 2).unwrap();
        for (a, b) in [
            (Point3::new(0.25, 0.25, -1.0), Point3::new(0.75, 0.5, 5.0)),
            (Point3::new(0.5, 0.5, 5.0), Point3::new(0.5, 0.5, -1.0)),
            (Point3::new(1.0, 0.0, -1.0), Point3::new(1.0, 0.0, 5.0)),
            (Point3::new(0.25, 0.25, 0.25), Point3::new(0.25, 0.25, 0.25)),
        ] {
            let budget = GridEnclosureBudget::default();
            assert_eq!(medium.column_density_box_enclosure(cx, boxed(a), boxed(b), budget).unwrap(),
                medium.column_density_enclosure(cx, a, b, budget).unwrap());
            assert_eq!(packet.transmittance_box(cx, boxed(a), boxed(b), budget).unwrap(),
                packet.transmittance(cx, a, b, budget).unwrap());
        }
    });
}

#[test]
fn layered_bundle_contains_independent_oblique_column_integrals() {
    with_cx(|cx, _| {
        let data = [2.0, 4.0];
        let medium = grid(cx, &data, [1, 1, 2], 16);
        let start = [Interval::new(0.25, 0.75), Interval::new(0.3, 0.7), Interval::point(-1.0)];
        let end = [Interval::point(0.5), Interval::point(0.5), Interval::point(3.0)];
        let result = medium.column_density_box_enclosure(cx, start, end, GridEnclosureBudget::default()).unwrap();
        // Two full 1 m layers, crossed obliquely: integral rho ds = 6*r/4.
        for x in [0.25_f64, 0.375, 0.5, 0.625, 0.75] {
            for y in [0.3_f64, 0.5, 0.7] {
                let r = (16.0 + (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)).sqrt();
                assert!(result.kg_per_m2.contains(1.5*r));
            }
        }
        assert!(result.kg_per_m2.lo() > 5.99);
        assert!(result.kg_per_m2.hi() < 6.1);
        assert_eq!(result.possible_cells, 2);
    });
}

#[test]
fn zero_containing_directions_keep_both_cell_ownerships_and_zero_length() {
    with_cx(|cx, _| {
        let data = [1.0, 8.0];
        let medium = grid(cx, &data, [2, 1, 1], 16);
        let start = [Interval::new(0.9, 1.1), Interval::point(0.5), Interval::point(-1.0)];
        let end = [Interval::new(0.9, 1.1), Interval::point(0.5), Interval::point(2.0)];
        let result = medium.column_density_box_enclosure(cx, start, end, GridEnclosureBudget::default()).unwrap();
        assert!(result.kg_per_m2.contains(1.0));
        assert!(result.kg_per_m2.contains(8.0));
        assert!(result.uncertain_parallel_cells > 0);
        let same = [Interval::new(0.25, 0.75); 3];
        let result = medium.column_density_box_enclosure(cx, same, same, GridEnclosureBudget::default()).unwrap();
        assert_eq!(result.kg_per_m2.lo(), 0.0);
        assert!(result.kg_per_m2.contains(0.75_f64.sqrt()));
    });
}

#[test]
fn independent_indexed_optics_enclose_every_wavelength_before_exponentiation() {
    with_cx(|cx, _| {
        let data = [1.0, 2.0, 4.0];
        let budget = CellPaletteBudget { max_models: 2, max_cells: 3 };
        let absorption = CellPalette::try_new(cx, vec![
            MassAbsorptionSpectrum::constant(0.25).unwrap(),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap(),
        ], Arc::from([0, 1, 0]), budget).unwrap();
        let grid = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([1, 1, 3], &data, [0.0; 3], [1.0; 3]), absorption, None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 32 }).unwrap();
        let scattering = CellPalette::try_new(cx, vec![
            ScatteringOptics { scattering: MassScatteringSpectrum::rayleigh(400.0, 0.5).unwrap(), phase: PhaseFunction::Rayleigh },
            ScatteringOptics { scattering: MassScatteringSpectrum::constant(1.0).unwrap(), phase: PhaseFunction::Isotropic },
        ], Arc::from([1, 0, 1]), budget).unwrap();
        let medium = ScatteringGrid::try_new_indexed(cx, grid, scattering).unwrap();
        let wavelengths = [400.0, 600.0];
        let packet = medium.prepare_extinction_enclosure(cx, &wavelengths, 8).unwrap();
        let a = [Interval::new(0.25, 0.75), Interval::point(0.5), Interval::point(-1.0)];
        let b = boxed(Point3::new(0.5, 0.5, 4.0));
        let result = packet.transmittance_box(cx, a, b, GridEnclosureBudget::default()).unwrap();
        for x in [0.25_f64, 0.5, 0.75] {
            let oblique = (25.0 + (x-0.5)*(x-0.5)).sqrt() / 5.0;
            for (lane, lambda) in wavelengths.into_iter().enumerate() {
                let ratio = 400.0/lambda;
                let tau = oblique * (6.25 + 2.0*((lambda-380.0)/400.0 + 0.5*ratio*ratio*ratio*ratio));
                assert!(result[lane].contains((-tau).exp()));
            }
        }
        assert!(result[0].lo() > result[1].hi());
    });
}

#[test]
fn box_queries_preserve_cell_node_and_cancellation_refusals() {
    with_cx(|cx, gate| {
        let data = [1.0; 8];
        let medium = grid(cx, &data, [1, 1, 8], 32);
        let a = [Interval::new(0.2, 0.8), Interval::point(0.5), Interval::point(-1.0)];
        let b = boxed(Point3::new(0.5, 0.5, 9.0));
        let identity = medium.content_identity();
        assert!(matches!(medium.column_density_box_enclosure(cx, a, b,
            GridEnclosureBudget { max_nodes_per_segment: 1 }), Err(VolumeTransportError::Budget)));
        let limited = grid(cx, &data, [1, 1, 8], 1);
        assert!(matches!(limited.column_density_box_enclosure(cx, a, b,
            GridEnclosureBudget::default()), Err(VolumeTransportError::Budget)));
        let mut invalid = a; invalid[0] = Interval::WHOLE;
        assert!(matches!(medium.column_density_box_enclosure(cx, invalid, b,
            GridEnclosureBudget::default()), Err(VolumeTransportError::Ray)));
        gate.request();
        assert!(matches!(medium.column_density_box_enclosure(cx, a, b,
            GridEnclosureBudget::default()), Err(VolumeTransportError::Cancelled)));
        assert_eq!(medium.content_identity(), identity);
    });
}
