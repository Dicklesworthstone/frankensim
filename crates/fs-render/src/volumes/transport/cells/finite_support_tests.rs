//! Original cell-face traversal against independent density-column cases.
use super::*;

fn column(
    grid: &VolumeGrid<'_>,
    origin: Point3,
    direction: Vec3,
    coefficient: f64,
    max_cells: u32,
) -> Result<(f64, u32), VolumeTransportError> {
    grid_column_with_checkpoint(grid, origin, direction, coefficient, max_cells, &mut || Ok(()))
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-14 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn g1_layered_original_cells_keep_forward_reverse_and_partial_columns() {
    let grid = VolumeGrid::new([3, 1, 1], &[2.0, 0.0, 5.0], [0.0; 3], [0.5, 2.0, 2.0]);
    let a = column(&grid, Point3::new(-2.0, 1.0, 1.0), Vec3::new(1.0, 0.0, 0.0), 0.4, 3).unwrap();
    let b = column(&grid, Point3::new(4.0, 1.0, 1.0), Vec3::new(-1.0, 0.0, 0.0), 0.4, 3).unwrap();
    close(a.0, 1.4);
    close(b.0, 1.4);
    assert_eq!((a.1, b.1), (3, 3));
    let part = column(&grid, Point3::new(0.25, 1.0, 1.0), Vec3::new(1.0, 0.0, 0.0), 0.4, 3).unwrap();
    close(part.0, 1.2);
    assert_eq!(grid_forward_interval(&grid, Point3::new(-2.0, 1.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(), Some((2.0, 1.5)));
    let distant = column(&grid, Point3::new(-1e12, 1.0, 1.0), Vec3::new(1.0, 0.0, 0.0), 0.4, 3).unwrap();
    close(distant.0, a.0);
}

#[test]
fn g1_oblique_and_exact_face_rays_use_the_native_original_cell_owner() {
    let grid = VolumeGrid::new([2, 2, 1], &[1.0, 9.0, 13.0, 3.0], [0.0; 3], [1.0; 3]);
    let d = 0.5_f64.sqrt();
    let diagonal = column(&grid, Point3::new(0.0, 0.0, 0.5), Vec3::new(d, d, 0.0), 0.5, 2).unwrap();
    close(diagonal.0, 2.0 / d);
    assert_eq!(diagonal.1, 2);
    let forward = column(&grid, Point3::new(1.0, 0.5, 0.5), Vec3::new(1.0, 0.0, 0.0), 1.0, 1).unwrap();
    let reverse = column(&grid, Point3::new(1.0, 0.5, 0.5), Vec3::new(-1.0, 0.0, 0.0), 1.0, 1).unwrap();
    assert_eq!(forward, (9.0, 1));
    assert_eq!(reverse, (1.0, 1));
    let decimal = VolumeGrid::new([3, 1, 1], &[2.0, 4.0, 8.0], [0.1, 0.0, 0.0], [0.1, 1.0, 1.0]);
    let face = 2.0_f64.mul_add(0.1, 0.1);
    close(column(&decimal, Point3::new(face, 0.5, 0.5), Vec3::new(-1.0, 0.0, 0.0), 1.0, 2).unwrap().0, 0.6);
}

#[test]
fn g0_vacuum_misses_and_zero_rate_require_no_positive_cell_budget() {
    let grid = VolumeGrid::new([1, 1, 1], &[4.0], [0.0; 3], [1.0; 3]);
    let direction = Vec3::new(1.0, 0.0, 0.0);
    let miss = Point3::new(-1.0, 2.0, 0.5);
    assert_eq!(column(&grid, miss, direction, 1.0, 0).unwrap(), (0.0, 0));
    assert_eq!(grid_forward_interval(&grid, miss, direction).unwrap(), None);
    assert_eq!(column(&grid, Point3::new(-1.0, 0.5, 0.5), direction, 0.0, 0).unwrap(), (0.0, 0));
    assert_eq!(column(&grid, Point3::new(-1.0, 0.5, 0.5), direction, 1.0, 0), Err(VolumeTransportError::Budget));
}

#[test]
fn g4_exact_crossing_budget_cancellation_and_unrepresentable_depth_refuse() {
    let grid = VolumeGrid::new([2, 1, 1], &[1.0, 2.0], [0.0; 3], [1.0; 3]);
    let origin = Point3::new(-1.0, 0.5, 0.5);
    let direction = Vec3::new(1.0, 0.0, 0.0);
    assert_eq!(column(&grid, origin, direction, 1.0, 1), Err(VolumeTransportError::Budget));
    let mut polls = 0;
    let cancelled = grid_column_with_checkpoint(&grid, origin, direction, 1.0, 2, &mut || {
        polls += 1;
        if polls == 3 { Err(VolumeTransportError::Cancelled) } else { Ok(()) }
    });
    assert_eq!(cancelled, Err(VolumeTransportError::Cancelled));
    assert_eq!(column(&grid, origin, direction, -1.0, 2), Err(VolumeTransportError::ScatteringSpectrum));
    let tiny = VolumeGrid::new([1, 1, 1], &[f64::MIN_POSITIVE], [0.0; 3], [1.0; 3]);
    assert_eq!(column(&tiny, origin, direction, f64::MIN_POSITIVE, 1), Err(VolumeTransportError::Unrepresentable));
    let thick = VolumeGrid::new([1, 1, 1], &[1e308], [0.0; 3], [10.0, 1.0, 1.0]);
    assert_eq!(column(&thick, origin, direction, 1.0, 1), Err(VolumeTransportError::Unrepresentable));
}
