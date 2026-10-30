use super::super::Medium;
use super::*;
use crate::fixtures::{kuhn_cube, single_tet};
use fs_rep_mesh::TetComplex;

mod molecular;

const AIR: Medium = Medium {
    density_kg_m3: 1.2,
    sound_speed_m_s: 343.0,
};
const TRANSPORT: BoundaryLayerTransport = BoundaryLayerTransport {
    dynamic_viscosity_pa_s: 1.84e-5,
    thermal_conductivity_w_m_k: 0.0257,
    cp_j_kg_k: 1005.0,
    heat_capacity_ratio: 1.4,
};
const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
    minimum_curvature_radius_m: f64::INFINITY,
    minimum_half_gap_m: 0.05,
    maximum_layer_ratio: 0.05,
};

// Complete coordinates invoke the original full FE system, with no modal
// truncation or artificial operator in the numerical verification.
fn complete(fe: &HelmholtzFe<'_>) -> HelmholtzModes {
    let n = fe.ndof();
    HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes: (0..n)
            .map(|i| {
                let mut x = vec![0.0; n];
                x[i] = 1.0;
                x
            })
            .collect(),
        below_low: 0,
        below_high: n,
    }
}

fn box_mesh(n: usize, dimensions: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, points) = kuhn_cube(n);
    (
        mesh,
        points
            .iter()
            .map(|p| core::array::from_fn(|i| p[i] * dimensions[i]))
            .collect(),
    )
}

fn end(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-14)
        })
        .collect()
}

fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn spatial_trace_is_work_conjugate_and_uniform_port_recovers_original_fe() {
    let (mesh, points) = single_tet();
    let points: Vec<_> = points.iter().map(|p| p.map(|v| 0.2 * v)).collect();
    let fe = HelmholtzFe::new(&mesh, &points, 1, AIR).unwrap();
    let faces = fe.boundary_faces();
    let basis = complete(&fe);
    let face = faces[0];
    let shape = [0.2, -0.7, 1.4];
    let ports = [
        BoundaryMode::uniform(&[face]),
        BoundaryMode {
            traces: vec![(face, shape)],
        },
    ];
    let model = InteriorMultiport::new(&fe, &basis, &ports, false, &|| false).unwrap();
    let frequency = 330.0;
    let prepared = model.prepare(frequency, &|| false).unwrap();
    let q = [C64::new(1e-6, -0.3e-6), C64::new(-0.2e-6, 0.7e-6)];
    let response = prepared.response(&q, &|| false).unwrap();
    let p: Vec<_> = face
        .iter()
        .map(|&v| {
            let value = response.amplitudes[v as usize];
            C64::new(value[0], value[1])
        })
        .collect();
    // Independent exact triangle moment: integral(lambda_i lambda_j)/A is
    // 1/6 on the diagonal and 1/12 off it. No product of face means is used.
    let sum_p = p.iter().fold(C64::ZERO, |s, &v| s + v);
    let expected = (sum_p.scale(shape.iter().sum())
        + p.iter()
            .zip(shape)
            .fold(C64::ZERO, |s, (p, v)| s + p.scale(v)))
    .scale(1.0 / 12.0);
    assert!(relative(response.pressure_pa[1], expected) < 3e-14);
    let averaged = sum_p.scale(shape.iter().sum::<f64>() / 9.0);
    assert!((expected - averaged).abs() > 0.05 * expected.abs());
    assert!(
        response.input_power_w().abs()
            < 1e-12
                * q.iter().map(|v| v.abs()).sum::<f64>()
                * response.pressure_pa.iter().map(|v| v.abs()).sum::<f64>()
    );

    let uniform = InteriorMultiport::new(&fe, &basis, &ports[..1], false, &|| false).unwrap();
    let old = ReducedHelmholtz::new(&fe, &basis, &[&[face]], false).unwrap();
    let old_response = old.response(frequency, 0, 0.01, &[[0.0, 0.0]]).unwrap();
    let new_response = uniform
        .prepare(frequency, &|| false)
        .unwrap()
        .response(&[C64::from_re(old_response.inflow)], &|| false)
        .unwrap();
    for (old, new) in old_response.amplitudes.iter().zip(&new_response.amplitudes) {
        assert!(relative(C64::new(new[0], new[1]), C64::new(old[0], old[1])) < 2e-12);
    }
}

#[test]
fn multiple_modes_retain_reciprocity_passivity_and_actual_wall_work() {
    let (mesh, points) = box_mesh(1, [0.3, 0.15, 0.12]);
    let fe = HelmholtzFe::new(&mesh, &points, 2, AIR).unwrap();
    let left = end(&fe, 0.0);
    let right = end(&fe, 0.3);
    let mut ports = vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)];
    ports.push(BoundaryMode {
        traces: left
            .iter()
            .map(|&face| {
                (
                    face,
                    face.map(|v| 2.0 * fe.positions()[v as usize][1] / 0.15 - 1.0),
                )
            })
            .collect(),
    });
    let walls: Vec<_> = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| !left.contains(f) && !right.contains(f))
        .collect();
    let basis = complete(&fe);
    let model = InteriorMultiport::new_with_boundary_layers(
        &fe,
        &basis,
        &ports,
        false,
        &walls,
        TRANSPORT,
        GEOMETRY,
        &|| false,
    )
    .unwrap();
    let frequency = 490.0;
    let z = model.impedance(frequency, &|| false).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (z[i * 3 + j] - z[j * 3 + i]).abs()
                    < 2e-12 * z.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
            );
        }
    }
    let prepared = model.prepare(frequency, &|| false).unwrap();
    for phase in 0..7 {
        let q = [
            C64::new(1e-6, 0.2e-6),
            C64::new(-0.3e-6, 0.4e-6),
            C64::new((phase as f64 - 3.0) * 0.1e-6, 0.7e-6),
        ];
        let response = prepared.response(&q, &|| false).unwrap();
        let loss = model
            .boundary_layer_power(frequency, &response.amplitudes)
            .unwrap();
        let total = loss.viscous_w + loss.thermal_w;
        assert!(total > 0.0);
        assert!((response.input_power_w() - total).abs() < 2e-9 * total);
        for row in 0..3 {
            let from_z = (0..3).fold(C64::ZERO, |s, col| s + z[row * 3 + col] * q[col]);
            assert!(relative(response.pressure_pa[row], from_z) < 2e-12);
        }
    }
    let r = [AIR.density_kg_m3 * AIR.sound_speed_m_s / (0.15 * 0.12); 3];
    let s = model.scattering(frequency, &r, &|| false).unwrap();
    // Contractivity tested on independent complex excitations. The exact
    // energy identity above is the positive dissipation mechanism.
    for v in [
        [C64::ONE, C64::ZERO, C64::ZERO],
        [C64::new(0.2, -0.7), C64::new(0.6, 0.1), C64::new(-0.4, 0.3)],
    ] {
        let incoming: f64 = v.iter().map(|z| z.norm_sq()).sum();
        let outgoing: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .fold(C64::ZERO, |sum, j| sum + s[i * 3 + j] * v[j])
                    .norm_sq()
            })
            .sum();
        assert!(outgoing < incoming);
    }
    for i in 0..3 {
        for j in 0..3 {
            assert!((s[i * 3 + j] - s[j * 3 + i]).abs() < 2e-12);
        }
    }
    assert_eq!(s, model.scattering(frequency, &r, &|| false).unwrap());
}

fn normalized_cosine_mode(fe: &HelmholtzFe<'_>, faces: &[[u32; 3]], width: f64) -> BoundaryMode {
    let mut traces: Vec<_> = faces
        .iter()
        .map(|&face| {
            (
                face,
                face.map(|v| {
                    det::sin_cos(core::f64::consts::PI * fe.positions()[v as usize][1] / width).1
                }),
            )
        })
        .collect();
    let mut integral = 0.0;
    let mut area = 0.0;
    for (face, v) in &traces {
        let (_, _, a) = fe.boundary_geometry(*face).unwrap();
        integral +=
            a / 12.0 * (v.iter().sum::<f64>().powi(2) + v.iter().map(|x| x * x).sum::<f64>());
        area += a;
    }
    let norm = det::sqrt(area / integral);
    for (_, v) in &mut traces {
        *v = v.map(|x| norm * x);
    }
    BoundaryMode { traces }
}

#[test]
fn straight_duct_plane_and_evanescent_mode_follow_independent_mesh_ladders() {
    let dimensions = [0.1, 0.2, 0.15];
    let frequency = 500.0;
    let omega = core::f64::consts::TAU * frequency;
    let k = omega / AIR.sound_speed_m_s;
    let area = dimensions[1] * dimensions[2];
    let (sin, cos) = det::sin_cos(k * dimensions[0]);
    let z0 = omega * AIR.density_kg_m3 / (k * area);
    let plane = [C64::new(0.0, -z0 * cos / sin), C64::new(0.0, -z0 / sin)];
    let kappa = det::sqrt((core::f64::consts::PI / dimensions[1]).powi(2) - k * k);
    let ep = det::exp(kappa * dimensions[0]);
    let em = 1.0 / ep;
    let z1 = omega * AIR.density_kg_m3 / (kappa * area);
    let evanescent = [
        C64::new(0.0, z1 * (ep + em) / (ep - em)),
        C64::new(0.0, 2.0 * z1 / (ep - em)),
    ];
    let mut errors = Vec::new();
    for n in [1, 2, 3] {
        let (mesh, points) = box_mesh(n, dimensions);
        let fe = HelmholtzFe::new(&mesh, &points, 2, AIR).unwrap();
        let left = end(&fe, 0.0);
        let right = end(&fe, dimensions[0]);
        let ports = [
            BoundaryMode::uniform(&left),
            BoundaryMode::uniform(&right),
            normalized_cosine_mode(&fe, &left, dimensions[1]),
            normalized_cosine_mode(&fe, &right, dimensions[1]),
        ];
        let model = InteriorMultiport::new(&fe, &complete(&fe), &ports, false, &|| false).unwrap();
        let z = model.impedance(frequency, &|| false).unwrap();
        let pe = relative(z[0], plane[0]).max(relative(z[1], plane[1]));
        let ee = relative(z[2 * 4 + 2], evanescent[0]).max(relative(z[2 * 4 + 3], evanescent[1]));
        assert!(z[2 * 4 + 2].im > 0.0); // retained reactive evanescence
        errors.push([pe, ee]);
        println!(
            "multiport n={n} dofs={} plane_error={pe:.6e} evanescent_error={ee:.6e}",
            model.dim()
        );
    }
    assert!(errors[1][0] < errors[0][0] && errors[2][0] < errors[1][0]);
    assert!(errors[1][1] < errors[0][1] && errors[2][1] < errors[1][1]);
    assert!(errors[2][0] < 3e-5, "plane: {errors:?}");
    assert!(errors[2][1] < 0.01, "evanescent: {errors:?}");
}

#[test]
fn invalid_modes_references_and_late_cancellation_publish_no_partial_matrix() {
    use core::cell::Cell;
    let (mesh, points) = single_tet();
    let fe = HelmholtzFe::new(&mesh, &points, 1, AIR).unwrap();
    let faces = fe.boundary_faces();
    let basis = complete(&fe);
    let mode = BoundaryMode::uniform(&faces[..1]);
    for bad in [
        BoundaryMode { traces: vec![] },
        BoundaryMode {
            traces: vec![(faces[0], [0.0; 3])],
        },
        BoundaryMode {
            traces: vec![(faces[0], [f64::NAN; 3])],
        },
        BoundaryMode {
            traces: vec![(faces[0], [1.0; 3]), (faces[0], [1.0; 3])],
        },
    ] {
        assert!(InteriorMultiport::new(&fe, &basis, &[bad], false, &|| false).is_err());
    }
    assert!(
        InteriorMultiport::new(
            &fe,
            &basis,
            &vec![mode.clone(); MAX_PORT_MODES + 1],
            false,
            &|| false
        )
        .is_err()
    );
    let model = InteriorMultiport::new(&fe, &basis, &[mode], false, &|| false).unwrap();
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(model.scattering(150.0, &[bad], &|| false).is_err());
        assert!(model.prepare(bad, &|| false).is_err());
    }
    let calls = Cell::new(0);
    let original = model
        .scattering(150.0, &[800.0], &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let count = calls.get();
    calls.set(0);
    assert!(
        model
            .scattering(150.0, &[800.0], &|| {
                calls.set(calls.get() + 1);
                calls.get() == count
            })
            .is_err()
    );
    assert_eq!(
        original,
        model.scattering(150.0, &[800.0], &|| false).unwrap()
    );
}
