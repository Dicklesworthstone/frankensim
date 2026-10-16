use super::super::historical_ring_piston;
use super::*;
use core::f64::consts::TAU;

fn medium() -> Medium { Medium { density: 1.2, sound_speed: 343.0 } }

fn disc(radius: f64, rings: usize) -> Vec<AreaSample> {
    let mut cells = Vec::new();
    for m in 0..rings {
        let r0 = radius * m as f64 / rings as f64;
        let r1 = radius * (m + 1) as f64 / rings as f64;
        let rc = f64::midpoint(r0, r1);
        let sectors = 6 * (m + 1);
        let area = PI * (r1 * r1 - r0 * r0) / sectors as f64;
        for sector in 0..sectors {
            let angle = TAU * (sector as f64 + 0.5) / sectors as f64;
            cells.push(AreaSample { position_m: [rc * det::cos(angle), rc * det::sin(angle)],
                area_m2: area, mode_shape: vec![1.0] });
        }
    }
    cells
}

#[test]
fn modal_piston_matches_existing_owner_and_independent_small_ka_terms() {
    let radius = 0.4;
    let cells = disc(radius, 24);
    let area = PI * radius * radius;
    let prepared = PreparedSurface::new(&cells).unwrap();
    let air = medium();
    for ka in [0.02, 0.05, 0.1] {
        let k = ka / radius;
        let z = prepared.impedance(k * air.sound_speed, air).unwrap()[0].scale(1.0 / area);
        let existing = historical_ring_piston(radius, k, air, 24).unwrap();
        assert!((z - existing).abs() / existing.abs() < 2e-12, "piston identity at ka={ka}");
        let normalized = z.scale(1.0 / (air.density * air.sound_speed));
        let resistance = 0.5 * ka * ka;
        let reactance = -8.0 * ka / (3.0 * PI);
        assert!((normalized.re / resistance - 1.0).abs() < 0.003, "small-ka resistance at {ka}: {normalized:?}");
        assert!((normalized.im / reactance - 1.0).abs() < 0.02, "small-ka added mass at {ka}: {normalized:?}");
    }
}

#[test]
fn signed_modal_cross_load_and_coordinate_changes_are_reciprocal() {
    let radius = 0.2;
    let mut cells = disc(radius, 5);
    for cell in &mut cells {
        let [x, y] = cell.position_m.map(|v| v / radius);
        cell.mode_shape = vec![1.0 + 0.3 * x, -0.2 + 0.5 * x + 0.4 * y];
    }
    let omega = TAU * 200.0;
    let z = PreparedSurface::new(&cells).unwrap().impedance(omega, medium()).unwrap();
    assert_eq!(z[1], z[2]);
    assert!(z[1].re < 0.0 && z[1].abs() > 0.01 * z[0].abs(), "signed cross-port radiation vanished");
    for v in [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, -2.0], [-0.2, 1.0]] {
        let power = (0..2).map(|i| (0..2).map(|j| v[i] * z[2 * i + j].re * v[j]).sum::<f64>()).sum::<f64>();
        assert!(power > 0.0);
    }
    let t = [[1.0, 0.7], [-0.2, 1.3]];
    let mut changed = cells.clone();
    for cell in &mut changed {
        let old = cell.mode_shape.clone();
        cell.mode_shape = (0..2).map(|j| (0..2).map(|i| old[i] * t[i][j]).sum()).collect();
    }
    let transformed = PreparedSurface::new(&changed).unwrap().impedance(omega, medium()).unwrap();
    for row in 0..2 { for column in 0..2 {
        let mut expected = C64::ZERO;
        for i in 0..2 { for j in 0..2 { expected = expected + z[2 * i + j].scale(t[i][row] * t[j][column]); } }
        assert!((transformed[2 * row + column] - expected).abs() < 2e-12 * z[0].abs());
    } }
}

#[test]
fn modal_resistance_agrees_with_independent_hemisphere_pressure_power() {
    let radius = 0.2;
    let mut cells = disc(radius, 6);
    for cell in &mut cells { cell.mode_shape = vec![1.0, cell.position_m[0] / radius]; }
    let air = medium();
    let omega = TAU * 150.0;
    let k = omega / air.sound_speed;
    let z = PreparedSurface::new(&cells).unwrap().impedance(omega, air).unwrap();
    let velocity = [C64::new(0.3, -0.2), C64::new(-0.1, 0.4)];
    let mut work = C64::ZERO;
    for i in 0..2 { for j in 0..2 { work = work + velocity[i].conj() * z[2 * i + j] * velocity[j]; } }
    let surface_power = 0.5 * work.re;
    // Independent outgoing pressure on a distant hemisphere, integrated in
    // mu=cos(theta) and azimuth. This does not call the loading kernel.
    let (polar, azimuth) = (32, 64);
    let mut radiated_power = 0.0;
    for p in 0..polar {
        let mu = (p as f64 + 0.5) / polar as f64;
        let transverse = det::sqrt(1.0 - mu * mu);
        for a in 0..azimuth {
            let angle = TAU * (a as f64 + 0.5) / azimuth as f64;
            let direction = [transverse * det::cos(angle), transverse * det::sin(angle)];
            let mut volume = C64::ZERO;
            for cell in &cells {
                let local = velocity[0].scale(cell.mode_shape[0]) + velocity[1].scale(cell.mode_shape[1]);
                let phase = -k * (direction[0] * cell.position_m[0] + direction[1] * cell.position_m[1]);
                volume = volume + local * C64::new(det::cos(phase), det::sin(phase)).scale(cell.area_m2);
            }
            let amplitude = volume * C64::new(0.0, -omega * air.density / TAU);
            radiated_power += amplitude.abs().powi(2) / (2.0 * air.density * air.sound_speed)
                * TAU / (polar * azimuth) as f64;
        }
    }
    assert!(surface_power > 0.0);
    assert!((surface_power / radiated_power - 1.0).abs() < 0.001,
        "surface={surface_power:e}, independent outgoing={radiated_power:e}");
}

#[test]
fn geometrically_similar_radiators_preserve_physical_area_scaling() {
    let cells = disc(0.1, 4);
    let omega = TAU * 400.0;
    let base = PreparedSurface::new(&cells).unwrap().impedance(omega, medium()).unwrap()[0];
    let mut larger = cells;
    for cell in &mut larger {
        cell.position_m = cell.position_m.map(|v| 3.0 * v);
        cell.area_m2 *= 9.0;
    }
    let result = PreparedSurface::new(&larger).unwrap().impedance(omega / 3.0, medium()).unwrap()[0];
    assert!((result - base.scale(9.0)).abs() / result.abs() < 2e-13);
}

fn complete_bending_surface() -> Vec<AreaSample> {
    let (nx, ny) = (20, 12);
    let (width, height) = (0.7, 0.4);
    let mut cells = Vec::new();
    for row in 0..ny {
        for column in 0..nx {
            let x = (column as f64 + 0.5) / nx as f64;
            let y = (row as f64 + 0.5) / ny as f64;
            cells.push(AreaSample {
                position_m: [width * x, height * y],
                area_m2: width * height / (nx * ny) as f64,
                mode_shape: (0..128).map(|mode| {
                    let sign = if mode % 3 == 0 { -1.0 } else { 1.0 };
                    sign * det::sin(PI * (mode % 16 + 1) as f64 * x)
                        * det::sin(PI * (mode / 16 + 1) as f64 * y)
                }).collect(),
            });
        }
    }
    cells
}

#[test]
fn full_structural_basis_retains_every_high_index_cross_load() {
    let cells = complete_bending_surface();
    let omega = TAU * 360.0;
    let full = PreparedSurface::new(&cells).unwrap().impedance(omega, medium()).unwrap();
    assert_eq!(full.len(), 128 * 128);
    for row in 0..128 {
        assert!(full[row * 128 + row].abs() > 0.0, "missing field {row}");
        for column in 0..128 {
            assert_eq!(full[row * 128 + column], full[column * 128 + row]);
        }
    }
    // Independent integrals of widely separated fields expose a retained
    // 32-field prefix or diagonal-only expansion of the physical operator.
    let selected = [0, 31, 32, 63, 95, 127];
    let mut selected_cells = cells;
    for cell in &mut selected_cells {
        cell.mode_shape = selected.iter().map(|&mode| cell.mode_shape[mode]).collect();
    }
    let direct = PreparedSurface::new(&selected_cells).unwrap().impedance(omega, medium()).unwrap();
    let scale = full.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    for (row, &i) in selected.iter().enumerate() {
        for (column, &j) in selected.iter().enumerate() {
            assert!((full[i * 128 + j] - direct[row * selected.len() + column]).abs() < 2e-12 * scale);
        }
    }
    assert!(full[31 * 128 + 127].abs() > 1e-4 * scale, "high-field cross load vanished");
}

#[test]
fn acoustic_port_projection_matches_independent_surface_field_projection() {
    let cells = complete_bending_surface();
    let prepared = PreparedSurface::new(&cells).unwrap();
    let ports = 4;
    // Orthonormal cosine columns mix ALL 128 signed structural fields. The
    // result is compared to newly integrated projected physical surface fields,
    // not to another algebraic implementation of the same matrix projection.
    let columns: Vec<Vec<f64>> = (0..ports).map(|port| {
        (0..128).map(|mode| {
            det::sqrt(2.0 / 128.0)
                * det::cos(PI * (mode as f64 + 0.5) * (port + 1) as f64 / 128.0)
        }).collect()
    }).collect();
    let mut projected_cells = cells;
    for cell in &mut projected_cells {
        cell.mode_shape = columns.iter().map(|column| {
            cell.mode_shape.iter().zip(column).map(|(shape, weight)| shape * weight).sum()
        }).collect();
    }
    let projected = PreparedSurface::new(&projected_cells).unwrap();
    for hz in [40.0, 900.0] {
        let full = prepared.impedance(TAU * hz, medium()).unwrap();
        let direct = projected.impedance(TAU * hz, medium()).unwrap();
        let scale = direct.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
        for row in 0..ports {
            for column in 0..ports {
                let mut congruence = C64::ZERO;
                for i in 0..128 {
                    for j in 0..128 {
                        congruence = congruence + full[i * 128 + j].scale(columns[row][i] * columns[column][j]);
                    }
                }
                assert!((direct[row * ports + column] - congruence).abs() < 2e-11 * scale,
                    "projected field work differs at {hz} Hz, ({row}, {column})");
            }
        }
    }
}

#[test]
fn invalid_or_unresolved_area_fields_refuse() {
    let good = disc(0.1, 2);
    assert!(PreparedSurface::new(&[]).is_err());
    let mut duplicate = good.clone(); duplicate[1].position_m = duplicate[0].position_m;
    assert!(PreparedSurface::new(&duplicate).is_err());
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut bad = good.clone(); bad[0].area_m2 = area;
        assert!(PreparedSurface::new(&bad).is_err());
    }
    let mut missing = good.clone(); missing[0].mode_shape.clear();
    assert!(PreparedSurface::new(&missing).is_err());
    let prepared = PreparedSurface::new(&good).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, TAU * 30_000.0] {
        assert!(prepared.impedance(omega, medium()).is_err());
    }
    assert!(prepared.impedance(100.0, Medium { density: f64::INFINITY, sound_speed: 343.0 }).is_err());
}

#[test]
fn expanded_fields_refuse_excess_work_before_geometry_allocation() {
    let cell = AreaSample { position_m: [0.0, 0.0], area_m2: 1e-4, mode_shape: Vec::new() };
    let mut oversized = vec![cell; MAX_DENSE_PANELS];
    oversized[0].mode_shape = vec![1.0; MAX_MODAL_RAYLEIGH_FIELDS];
    assert!(matches!(PreparedSurface::new(&oversized), Err(HelmholtzError::BadParameter {
        what: "Rayleigh loading exceeds the per-frequency modal work limit"
    })));
    let mut excess_fields = disc(0.1, 1);
    for cell in &mut excess_fields { cell.mode_shape = vec![1.0; MAX_MODAL_RAYLEIGH_FIELDS + 1]; }
    assert!(PreparedSurface::new(&excess_fields).is_err());
}
