//! Regression for the raw-volume smoothing stage of planar wetting.
//!
//! The older angle tests manufactured an already smoothed field. These tests
//! initialize tracked mass and populations, so the production smoothing step
//! must preserve an affine plane before its normal and pressure are evaluated.

use super::*;

const WALLS: [PlanarWall3; 6] = [
    PlanarWall3::XMin,
    PlanarWall3::XMax,
    PlanarWall3::YMin,
    PlanarWall3::YMax,
    PlanarWall3::ZMin,
    PlanarWall3::ZMax,
];

fn plane_fill(position: [i64; 3], wall: PlanarWall3, angle: f64, sign: f64) -> f64 {
    let axis = wall.axis();
    let inward = if wall.upper() {
        7 - position[axis]
    } else {
        position[axis]
    };
    // Two nonzero tangential components exercise the three-dimensional
    // gradient norm. Their magnitude is exactly 0.02 analytically.
    0.5 - sign * 0.012 * (position[(axis + 1) % 3] as f64 - 3.0)
        - 0.016 * (position[(axis + 2) % 3] as f64 - 3.0)
        - 0.02 * (angle.cos() / angle.sin()) * inward as f64
}

fn raw_plane(wall: PlanarWall3, angle: f64, sign: f64) -> FreeSurface3 {
    let mut surface = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.001,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    for (&key, tile) in &mut surface.tiles {
        let slot = surface.grid.slot_of(key).unwrap();
        for lane in 0..TILE_CELLS {
            let (x, y, z) = FreeSurface3::global_of(key, lane);
            let fill = plane_fill([x, y, z], wall, angle, sign);
            assert!(fill > 0.0 && fill < 1.0);
            let density: f64 = surface.grid.populations(slot, lane).iter().sum();
            tile.cells[lane] = Cell3::Interface;
            tile.mass[lane] = density * fill;
        }
    }
    surface
}

#[test]
fn g1_raw_planar_wetting_preserves_fill_and_zero_capillary_pressure_on_all_faces() {
    let gate = CancelGate::new_clock_free();
    for wall in WALLS {
        for angle in [
            core::f64::consts::PI / 3.0,
            core::f64::consts::FRAC_PI_2,
            2.0 * core::f64::consts::PI / 3.0,
        ] {
            for sign in [-1.0, 1.0] {
                let mut surface = raw_plane(wall, angle, sign);
                surface.set_planar_contact_angle(wall, angle).unwrap();
                let mass = surface.ledger_mass();
                let raw = surface.refresh_fill(&gate).unwrap();
                assert_eq!(mass, surface.ledger_mass());
                // Both wall-adjacent layers belong to the fallback curvature
                // stencil; neither may acquire an artificial normal bend.
                for inward in [0, 1] {
                    let mut p = [3_i64; 3];
                    p[wall.axis()] = if wall.upper() { 7 - inward } else { inward };
                    let (slot, lane) = surface.grid.resolve_source(p[0], p[1], p[2]).unwrap();
                    let key = surface.grid.active_keys()[slot];
                    let expected = plane_fill(p, wall, angle, sign);
                    assert!((raw[&key][lane] - expected).abs() < 2e-15);
                    assert!(
                        (surface.fill_smooth[&key][lane] - expected).abs() < 2e-14,
                        "{wall:?}, angle={angle}, depth={inward}: smoothing bent the plane"
                    );
                    let curvature = surface.curvature(key, lane);
                    assert!(
                        curvature.abs() < 2e-12,
                        "{wall:?}, angle={angle}, depth={inward}: curvature={curvature}"
                    );
                    assert!((surface.reference_density(key, lane) - 1.0).abs() < 2e-14);
                }
            }
        }
    }
}

#[test]
fn g3_per_cell_contact_angles_share_the_uniform_raw_smoothing_operator() {
    let gate = CancelGate::new_clock_free();
    for wall in WALLS {
        for angle in [
            core::f64::consts::PI / 3.0,
            2.0 * core::f64::consts::PI / 3.0,
        ] {
            let mut uniform = raw_plane(wall, angle, -1.0);
            let mut spatial = uniform.clone();
            uniform.set_planar_contact_angle(wall, angle).unwrap();
            let axis = wall.axis();
            let mut angles = BTreeMap::new();
            for first in 0..8 {
                for second in 0..8 {
                    let mut p = [0; 3];
                    p[axis] = if wall.upper() { 7 } else { 0 };
                    p[(axis + 1) % 3] = first;
                    p[(axis + 2) % 3] = second;
                    angles.insert(p, angle);
                }
            }
            spatial.set_planar_contact_field(wall, angles, &gate).unwrap();
            let uniform_raw = uniform.refresh_fill(&gate).unwrap();
            let spatial_raw = spatial.refresh_fill(&gate).unwrap();
            assert_eq!(uniform_raw, spatial_raw);
            assert_eq!(uniform.fill_smooth, spatial.fill_smooth);
            let mut p = [3; 3];
            p[axis] = if wall.upper() { 7 } else { 0 };
            let (slot, lane) = spatial.grid.resolve_source(p[0], p[1], p[2]).unwrap();
            let key = spatial.grid.active_keys()[slot];
            assert!(spatial.curvature(key, lane).abs() < 2e-12);
            assert!((spatial.reference_density(key, lane) - 1.0).abs() < 2e-14);
        }
    }
}

#[test]
fn g1_neutral_transverse_corner_preserves_raw_plane_and_laplace_pressure() {
    let gate = CancelGate::new_clock_free();
    for wall in WALLS {
        let axis = wall.axis();
        for graph_axis in (0..3).filter(|&d| d != axis) {
            let neutral = 3 - axis - graph_axis;
            for outward in [-1.0, 1.0] {
                for angle in [1.2_f64, core::f64::consts::FRAC_PI_2, 1.9] {
                    let cotangent = if angle == core::f64::consts::FRAC_PI_2 {
                        0.0
                    } else {
                        angle.cos() / angle.sin()
                    };
                    // Independent affine volume field with n.n_wall=cos(theta),
                    // n.n_transverse=0 and the prescribed graph orientation.
                    let plane = |p: [i64; 3]| {
                        let inward = if wall.upper() { 7 - p[axis] } else { p[axis] };
                        0.5 - 0.02 * outward * (p[graph_axis] as f64 - 3.0)
                            - 0.02 * cotangent * inward as f64
                    };
                    let mut surface = FreeSurface3::new(
                        8, 8, 8, 0.8, [0.0; 3], 0.001,
                        ContactModel3::Neutral, |_, _, _| true,
                    ).unwrap();
                    let mut angles = BTreeMap::new();
                    for (&key, tile) in &mut surface.tiles {
                        let slot = surface.grid.slot_of(key).unwrap();
                        for lane in 0..TILE_CELLS {
                            let (x, y, z) = FreeSurface3::global_of(key, lane);
                            let p = [x, y, z];
                            let density: f64 = surface.grid.populations(slot, lane).iter().sum();
                            tile.cells[lane] = Cell3::Interface;
                            tile.mass[lane] = density * plane(p);
                            if p[axis] == if wall.upper() { 7 } else { 0 } {
                                angles.insert(p.map(|v| v as usize), angle);
                            }
                        }
                    }
                    surface.set_planar_contact_field(wall, angles, &gate).unwrap();
                    surface.enable_planar_neutral_corners(graph_axis, &gate).unwrap();
                    let masses = surface.mass_snapshot();
                    let populations = surface.grid().state_bits();
                    let raw = surface.refresh_fill(&gate).unwrap();
                    for edge in [0, 7] {
                        let mut ghost = [3_i64; 3];
                        ghost[axis] = if wall.upper() { 8 } else { -1 };
                        ghost[neutral] = if edge == 0 { -1 } else { 8 };
                        let value = surface.planar_angle_ghost_in(ghost, &raw).unwrap();
                        assert!((value - plane(ghost)).abs() < 2e-14);
                        for inward in [0, 1] {
                            let mut p = [3_i64; 3];
                            p[axis] = if wall.upper() { 7 - inward } else { inward };
                            p[neutral] = edge;
                            let (slot, lane) = surface.grid.resolve_source(p[0], p[1], p[2]).unwrap();
                            let key = surface.grid.active_keys()[slot];
                            assert!((surface.fill_smooth[&key][lane] - plane(p)).abs() < 2e-14,
                                "{wall:?}, graph={graph_axis}, edge={edge}, depth={inward}");
                            assert!(surface.curvature(key, lane).abs() < 2e-12,
                                "corner smoothing must not manufacture Laplace pressure");
                            assert!((surface.reference_density(key, lane) - 1.0).abs() < 2e-14);
                        }
                    }
                    assert_eq!(surface.mass_snapshot(), masses);
                    assert_eq!(surface.grid().state_bits(), populations);
                }
            }
        }
    }
}
