use super::*;
use crate::fixtures::kuhn_cube;
use fs_modal::SliceOptions;
use fs_rep_mesh::TetComplex;

pub(super) const AIR: Medium = Medium {
    density_kg_m3: 1.2,
    sound_speed_m_s: 343.0,
};
pub(super) const TRANSPORT: BoundaryLayerTransport = BoundaryLayerTransport {
    dynamic_viscosity_pa_s: 1.84e-5,
    thermal_conductivity_w_m_k: 0.0257,
    cp_j_kg_k: 1005.0,
    heat_capacity_ratio: 1.4,
};
pub(super) const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
    minimum_curvature_radius_m: f64::INFINITY,
    minimum_half_gap_m: 0.02,
    maximum_layer_ratio: 0.05,
};

pub(super) fn box_mesh(n: usize, dimensions: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, points) = kuhn_cube(n);
    (
        complex,
        points
            .iter()
            .map(|p| std::array::from_fn(|i| p[i] * dimensions[i]))
            .collect(),
    )
}

pub(super) fn end_faces(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions[v as usize][0] - x).abs() < 1e-14)
        })
        .collect()
}

pub(super) fn side_faces(fe: &HelmholtzFe<'_>, length: f64) -> Vec<[u32; 3]> {
    let ends: BTreeSet<_> = end_faces(fe, 0.0)
        .into_iter()
        .chain(end_faces(fe, length))
        .collect();
    fe.boundary_faces()
        .into_iter()
        .filter(|face| !ends.contains(face))
        .collect()
}

pub(super) fn quadratic(matrix: &Csr, x: &[f64]) -> f64 {
    let mut y = vec![0.0; x.len()];
    matrix.spmv(x, &mut y);
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

// A complete coordinate basis makes the existing Galerkin solver a full FE
// solve. No eigenbasis truncation or synthetic numerical operator is involved.
pub(super) fn complete_basis(fe: &HelmholtzFe<'_>) -> HelmholtzModes {
    let n = fe.ndof();
    HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes: (0..n)
            .map(|i| {
                let mut v = vec![0.0; n];
                v[i] = 1.0;
                v
            })
            .collect(),
        below_low: 0,
        below_high: n,
    }
}

#[test]
fn original_trace_surface_gradient_matches_triangle_energy_all_orders() {
    let (complex, positions) = box_mesh(1, [0.7, 0.3, 0.2]);
    let gradient = [1.3, -0.8, 0.4];
    for order in 1..=4 {
        let fe = HelmholtzFe::new(&complex, &positions, order, AIR).unwrap();
        for face in fe.boundary_faces() {
            let (vertices, normal, area) = fe.boundary_geometry(face).unwrap();
            let stiffness = fe.boundary_tangential_stiffness(&[face]).unwrap();
            let mut coefficients = vec![0.0; fe.ndof()];
            for (i, point) in positions.iter().enumerate() {
                coefficients[i] = 0.37 + dot3(gradient, *point);
            }
            let expected = area * (dot3(gradient, gradient) - dot3(gradient, normal).powi(2));
            let actual = quadratic(&stiffness, &coefficients);
            assert!(
                (actual - expected).abs() < 2e-13 * expected.max(1.0),
                "order={order} {actual} {expected}"
            );
            let reversed = fe
                .boundary_tangential_stiffness(&[[face[2], face[1], face[0]]])
                .unwrap();
            assert!((quadratic(&reversed, &coefficients) - actual).abs() < 1e-14);
            coefficients.fill(0.0);
            for (i, point) in positions.iter().enumerate() {
                coefficients[i] = dot3(normal, sub3(*point, vertices[0]));
            }
            assert!(
                quadratic(&stiffness, &coefficients).abs() < 1e-13,
                "normal gradient leaked"
            );
        }
        let face = fe.boundary_faces()[0];
        assert!(fe.boundary_tangential_stiffness(&[face, face]).is_err());
        assert!(
            fe.boundary_tangential_stiffness_controlled(&[face], &|| true)
                .is_err()
        );
    }
}

#[test]
fn penetration_depths_and_computed_band_refuse_unresolved_layers() {
    let (low, high) = TRANSPORT.frequency_band_hz(AIR, GEOMETRY).unwrap();
    assert!(low > 0.0 && high > low);
    let f = (low * high).sqrt();
    let a = TRANSPORT.scales(AIR, GEOMETRY, f).unwrap();
    let b = TRANSPORT.scales(AIR, GEOMETRY, 4.0 * f).unwrap();
    assert!((a.viscous_depth_m / b.viscous_depth_m - 2.0).abs() < 1e-14);
    assert!((a.thermal_depth_m / b.thermal_depth_m - 2.0).abs() < 1e-14);
    assert!(TRANSPORT.scales(AIR, GEOMETRY, low * 0.99).is_err());
    assert!(TRANSPORT.scales(AIR, GEOMETRY, high * 1.01).is_err());
    let thin_gap = BoundaryLayerGeometry {
        minimum_half_gap_m: 1e-6,
        ..GEOMETRY
    };
    assert!(TRANSPORT.frequency_band_hz(AIR, thin_gap).is_err());
    let bad = BoundaryLayerTransport {
        cp_j_kg_k: f64::NAN,
        ..TRANSPORT
    };
    assert!(bad.scales(AIR, GEOMETRY, f).is_err());
    assert!(
        TRANSPORT
            .scales(
                AIR,
                BoundaryLayerGeometry {
                    maximum_layer_ratio: 0.1,
                    ..GEOMETRY
                },
                f
            )
            .is_err()
    );
}

#[test]
fn full_fe_boundary_layer_input_work_closes_and_zero_loss_is_bitwise() {
    let (complex, positions) = box_mesh(2, [0.5, 0.04, 0.04]);
    let fe = HelmholtzFe::new(&complex, &positions, 2, AIR).unwrap();
    let source = end_faces(&fe, 0.0);
    let outlet = end_faces(&fe, 0.5);
    let walls = side_faces(&fe, 0.5);
    let modes = complete_basis(&fe);
    let ports = [source.as_slice(), outlet.as_slice()];
    let loss = ReducedHelmholtz::new_with_boundary_layers(
        &fe, &modes, &ports, false, &walls, TRANSPORT, GEOMETRY,
    )
    .unwrap();
    let baseline = ReducedHelmholtz::new(&fe, &modes, &ports, false).unwrap();
    let zero = ReducedHelmholtz::new_with_boundary_layers(
        &fe,
        &modes,
        &ports,
        false,
        &walls,
        BoundaryLayerTransport {
            dynamic_viscosity_pa_s: 0.0,
            thermal_conductivity_w_m_k: 0.0,
            ..TRANSPORT
        },
        GEOMETRY,
    )
    .unwrap();
    assert!(
        ReducedHelmholtz::new_with_boundary_layers(
            &fe, &modes, &ports, false, &source, TRANSPORT, GEOMETRY
        )
        .is_err()
    );
    for frequency in [200.0, 600.0, 1200.0] {
        let admittances = [[0.0, 0.0], [0.7, 0.0]];
        let response = loss.response(frequency, 0, 0.01, &admittances).unwrap();
        let power = loss
            .boundary_layer_power(frequency, &response.amplitudes)
            .unwrap();
        assert!(power.viscous_w > 0.0 && power.thermal_w > 0.0);
        let input = 0.5 * response.mean_pressure[0][0] * response.inflow;
        // Real pressure admittance dissipates pointwise, not merely at the
        // product of port-averaged pressure and flow. Integrate the original C.
        let mut pressure_real = vec![0.0; fe.ndof()];
        let mut pressure_imag = vec![0.0; fe.ndof()];
        for (i, p) in response.amplitudes.iter().enumerate() {
            pressure_real[i] = p[0];
            pressure_imag[i] = p[1];
        }
        let port_mass = fe.boundary_mass(&outlet).unwrap();
        let outlet_power = 0.5 * 0.7 / (AIR.density_kg_m3 * AIR.sound_speed_m_s)
            * (quadratic(&port_mass, &pressure_real) + quadratic(&port_mass, &pressure_imag));
        let defect = input - outlet_power - power.viscous_w - power.thermal_w;
        assert!(
            defect.abs() < 5e-10 * input.abs(),
            "f={frequency}, balance={defect:e}, input={input:e}"
        );
        let reference = baseline.response(frequency, 0, 0.01, &admittances).unwrap();
        let no_loss = zero.response(frequency, 0, 0.01, &admittances).unwrap();
        assert_eq!(
            reference
                .amplitudes
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            no_loss
                .amplitudes
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            zero.boundary_layer_power(frequency, &no_loss.amplitudes)
                .unwrap(),
            BoundaryLayerPower::default()
        );
        assert_ne!(response.amplitudes, reference.amplitudes);
    }
}

#[test]
fn axial_box_mode_wall_damping_converges_to_analytic_surface_integrals() {
    let dimensions = [0.5, 0.04, 0.04];
    let perimeter_over_area =
        2.0 * (dimensions[1] + dimensions[2]) / (dimensions[1] * dimensions[2]);
    let exact_lambda = (core::f64::consts::PI / dimensions[0]).powi(2);
    let mut errors = Vec::new();
    for mesh in [1usize, 2, 3] {
        let (complex, positions) = box_mesh(mesh, dimensions);
        let fe = HelmholtzFe::new(&complex, &positions, 2, AIR).unwrap();
        let modes = fe.modes((-1.0, 500.0), &SliceOptions::default()).unwrap();
        assert_eq!(modes.lambda.len(), 2);
        let p = &modes.shapes[1];
        let mass = quadratic(fe.mass(), p);
        let walls = side_faces(&fe, dimensions[0]);
        let surface_mass = quadratic(&fe.boundary_mass(&walls).unwrap(), p);
        let surface_stiffness = quadratic(&fe.boundary_tangential_stiffness(&walls).unwrap(), p);
        let frequency = AIR.sound_speed_m_s / (2.0 * dimensions[0]);
        let scales = TRANSPORT.scales(AIR, GEOMETRY, frequency).unwrap();
        let observed_viscous = 0.5 * scales.viscous_depth_m * surface_stiffness / mass;
        let exact_viscous = 0.5 * scales.viscous_depth_m * perimeter_over_area * exact_lambda;
        let observed_thermal =
            0.5 * (TRANSPORT.heat_capacity_ratio - 1.0) * scales.thermal_depth_m * surface_mass
                / mass;
        let exact_thermal = 0.5
            * (TRANSPORT.heat_capacity_ratio - 1.0)
            * scales.thermal_depth_m
            * perimeter_over_area;
        let error = (observed_viscous / exact_viscous - 1.0)
            .abs()
            .max((observed_thermal / exact_thermal - 1.0).abs());
        println!("boundary-layer axial mode mesh={mesh} error={error:e}");
        errors.push(error);
    }
    assert!(errors.windows(2).all(|w| w[1] < w[0]), "{errors:?}");
    assert!(errors[2] < 0.003, "{errors:?}");
}
