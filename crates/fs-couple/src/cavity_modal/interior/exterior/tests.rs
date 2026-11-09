use super::*;
use fs_bem::reciprocal_radiation::Budget;
use fs_feec::acoustics::helmholtz::Medium as FeMedium;
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::GasSpec;

mod sources;

const LENGTH: f64 = 0.06;
const WIDTH: f64 = 0.04;
const FREQUENCY: f64 = 480.0;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn medium(gas: &GasState) -> FeMedium {
    FeMedium {
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    }
}
fn budget() -> RadiationBudget {
    RadiationBudget {
        maximum_pressure_dofs: 128,
        maximum_panels: 128,
        maximum_trace_entries: 1024,
        maximum_projection_terms: 32768,
    }
}
fn route() -> ApertureRoute {
    ApertureRoute::Galerkin {
        budget: Budget {
            maximum_unknowns: 256,
            maximum_cache_bytes: 32 * 1024 * 1024,
            maximum_work: 200_000_000,
            angular_refinement: 0,
        },
    }
}
fn patch(fe: &HelmholtzFe<'_>, axis: usize, coordinate: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][axis] - coordinate).abs() < 1e-14)
        })
        .collect()
}
fn signed(fe: &HelmholtzFe<'_>, faces: &[[u32; 3]]) -> BoundaryMode {
    BoundaryMode {
        traces: faces
            .iter()
            .map(|&face| {
                (
                    face,
                    face.map(|v| 2.0 * fe.positions()[v as usize][1] / WIDTH - 1.0),
                )
            })
            .collect(),
    }
}
fn ports(fe: &HelmholtzFe<'_>) -> Vec<BoundaryMode> {
    let left = patch(fe, 0, 0.0);
    let right = patch(fe, 0, LENGTH);
    vec![
        BoundaryMode::uniform(&left),
        signed(fe, &left),
        BoundaryMode::uniform(&right),
        signed(fe, &right),
        BoundaryMode::uniform(&patch(fe, 2, WIDTH)),
    ]
}
fn surface(fe: &HelmholtzFe<'_>) -> SpherePanels {
    let triangles = fe
        .boundary_faces()
        .into_iter()
        .map(|face| {
            let mut triangle = face.map(|v| fe.positions()[v as usize]);
            let (_, normal, _) = fe.boundary_geometry(face).unwrap();
            let a: [f64; 3] = std::array::from_fn(|i| triangle[1][i] - triangle[0][i]);
            let b: [f64; 3] = std::array::from_fn(|i| triangle[2][i] - triangle[0][i]);
            let cross = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            if cross.iter().zip(normal).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
                triangle.swap(1, 2);
            }
            // Winding and vertex numbering need not equal the FE's ordering.
            triangle.rotate_left(1);
            triangle
        })
        .collect();
    SpherePanels::from_triangles(triangles).unwrap()
}
fn source_loads(fe: &HelmholtzFe<'_>, ports: &[BoundaryMode]) -> Vec<Vec<f64>> {
    ports
        .iter()
        .map(|port| {
            let faces: Vec<_> = port.traces.iter().map(|(face, _)| *face).collect();
            let (_, area) = fe.boundary_load(&faces).unwrap();
            fe.boundary_linear_load(&port.traces)
                .unwrap()
                .into_iter()
                .map(|b| b / area)
                .collect()
        })
        .collect()
}
fn source_apertures(
    fe: &HelmholtzFe<'_>,
    ports: &[BoundaryMode],
    exterior: &SpherePanels,
    selected: &[usize],
) -> Vec<Aperture> {
    selected
        .iter()
        .map(|&mode| {
            let faces: Vec<_> = ports[mode].traces.iter().map(|(f, _)| *f).collect();
            let mut shape = vec![0.0; exterior.areas().len()];
            for &(face, values) in &ports[mode].traces {
                let points = face.map(|v| fe.positions()[v as usize]);
                let panel = exterior
                    .triangles()
                    .unwrap()
                    .iter()
                    .position(|triangle| points.iter().all(|p| triangle.contains(p)))
                    .unwrap();
                shape[panel] = values.iter().sum::<f64>() / 3.0;
            }
            let (_, normal, _) = fe.boundary_geometry(faces[0]).unwrap();
            Aperture {
                shape,
                normal,
                reference_point: fe.positions()[faces[0][0] as usize],
            }
        })
        .collect()
}
fn near(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * scale.max(a.abs()).max(b.abs()),
        "{a:?} != {b:?}"
    );
}

/// Independent full mixed block, using original CSR K/M and original FE trace
/// integration directly. It does not use the adapter's scattering, pressure
/// response, reduced matrices or port-closure assembly.
fn direct_block(
    fe: &HelmholtzFe<'_>,
    ports: &[BoundaryMode],
    radiating: &[usize],
    driven: &[usize],
    source: &[C64],
    z_plus: &[C64],
) -> Vec<C64> {
    direct_block_with_load(fe, ports, radiating, driven, source, z_plus, None)
}

fn direct_block_with_load(
    fe: &HelmholtzFe<'_>,
    ports: &[BoundaryMode],
    radiating: &[usize],
    driven: &[usize],
    source: &[C64],
    z_plus: &[C64],
    distributed: Option<&[C64]>,
) -> Vec<C64> {
    let pressure_dofs = fe.ndof();
    assert!(distributed.is_none_or(|load| load.len() == pressure_dofs));
    let nr = radiating.len();
    let n = pressure_dofs + nr;
    let omega = core::f64::consts::TAU * FREQUENCY;
    let k = omega / fe.medium().sound_speed_m_s;
    let iwrho = C64::new(0.0, omega * fe.medium().density_kg_m3);
    let loads = source_loads(fe, ports);
    let mut matrix = vec![C64::ZERO; n * n];
    let mut rhs = vec![C64::ZERO; n];
    for i in 0..pressure_dofs {
        for (a, scale) in [(fe.stiffness(), 1.0), (fe.mass(), -k * k)] {
            let (columns, values) = a.row(i);
            for (&j, &v) in columns.iter().zip(values) {
                matrix[i * n + j] = matrix[i * n + j] + C64::from_re(scale * v);
            }
        }
        for (&port, &q) in driven.iter().zip(source) {
            rhs[i] = rhs[i] + iwrho * q.scale(loads[port][i]);
        }
        if let Some(load) = distributed {
            rhs[i] = rhs[i] + iwrho * load[i];
        }
        for (r, &port) in radiating.iter().enumerate() {
            matrix[i * n + pressure_dofs + r] = iwrho.scale(-loads[port][i]);
            matrix[(pressure_dofs + r) * n + i] = C64::from_re(loads[port][i]);
        }
    }
    for r in 0..nr {
        for s in 0..nr {
            matrix[(pressure_dofs + r) * n + pressure_dofs + s] = z_plus[r * nr + s];
        }
    }
    lu_complex_with_cancel(&matrix, n, || false)
        .unwrap()
        .solve(&mut rhs);
    rhs
}

#[test]
fn complete_pressure_field_matches_original_fe_bem_block_and_reaction_power() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    // Order two includes edge pressure DOFs: retaining only vertices fails.
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    assert!(fe.ndof() > points.len());
    let ports = ports(&fe);
    let surface = surface(&fe);
    let selected = [2, 3, 4];
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &selected,
        &surface,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(model.interior().dim(), fe.ndof());
    assert_eq!(model.driven_modes(), &[0, 1]);
    let reference = gas.density * gas.sound_speed / (WIDTH * WIDTH);
    let mut prepared = model
        .prepare(FREQUENCY, &[reference; 5], &route(), &gate)
        .unwrap();
    let feeds = [C64::new(1e-6, -0.3e-6), C64::new(0.4e-6, 0.7e-6)];
    let actual = prepared.response(&feeds, &gate).unwrap();
    for (mode, &prescribed) in feeds.iter().enumerate() {
        near(actual.interior.inward_flow_m3_s[mode], prescribed, 1e-6);
    }
    let native = aperture_impedance(
        &surface,
        &source_apertures(&fe, &ports, &surface, &selected),
        &[core::f64::consts::TAU * FREQUENCY],
        Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        &route(),
        &gate,
    )
    .unwrap();
    let areas = model.interior().port_areas_m2();
    let z: Vec<_> = native.rows[0]
        .impedance
        .iter()
        .enumerate()
        .map(|(entry, z)| {
            z.conj()
                .scale(1.0 / (areas[selected[entry / 3]] * areas[selected[entry % 3]]))
        })
        .collect();
    for (&a, &b) in prepared.reaction_impedance().iter().zip(&z) {
        near(a, b, reference);
    }
    // Off-diagonal coupling between distinct physical openings is material.
    assert!(
        z[2].abs() > 1e-4 * z[0].abs(),
        "missing right/top mutual reaction"
    );
    let expected = direct_block(&fe, &ports, &selected, &[0, 1], &feeds, &z);
    assert_eq!(actual.interior.amplitudes.len(), fe.ndof());
    for (&p, &expected) in actual.interior.amplitudes.iter().zip(&expected) {
        near(C64::new(p[0], p[1]), expected, 1.0);
    }
    let loads = source_loads(&fe, &ports);
    for (mode, load) in loads.iter().enumerate() {
        let pressure = load
            .iter()
            .zip(&expected)
            .fold(C64::ZERO, |sum, (&b, &p)| sum + p.scale(b));
        near(actual.interior.pressure_pa[mode], pressure, 1.0);
    }
    for (r, &mode) in selected.iter().enumerate() {
        near(
            actual.interior.inward_flow_m3_s[mode],
            expected[fe.ndof() + r],
            1e-6,
        );
        near(
            actual.exterior_pressure_pa[r],
            actual.interior.pressure_pa[mode],
            1.0,
        );
    }
    assert!(
        actual.outward_flow_m3_s[1].abs() > 1e-12,
        "zero-net mode was discarded"
    );
    let mut native_power = 0.0;
    for r in 0..3 {
        for s in 0..3 {
            native_power += 0.5
                * (actual.outward_flow_m3_s[r].conj() * z[r * 3 + s] * actual.outward_flow_m3_s[s])
                    .re;
        }
    }
    assert!(native_power > 0.0);
    assert!((actual.radiation_power_w - native_power).abs() < 1e-10 * native_power);
    assert!((actual.source_power_w - native_power).abs() < 2e-8 * native_power);
    assert!(actual.interior.input_power_w().abs() < 2e-8 * native_power);
    assert!(actual.port_relative_residual < 1e-10);
    assert!(actual.power_balance_residual_w.abs() < 2e-8 * native_power);

    // A finite reactive opening differs from an ideal pressure-release cap.
    let release = direct_block(&fe, &ports, &selected, &[0, 1], &feeds, &[C64::ZERO; 9]);
    let field_change = expected[..fe.ndof()]
        .iter()
        .zip(&release)
        .map(|(&p, &q)| (p - q).abs())
        .fold(0.0_f64, f64::max);
    assert!(field_change > 1e-5);

    // Reference coordinates can change strongly without changing the field.
    let other = model
        .prepare(
            FREQUENCY,
            &[
                0.25 * reference,
                3.0 * reference,
                2.0 * reference,
                0.5 * reference,
                4.0 * reference,
            ],
            &route(),
            &gate,
        )
        .unwrap();
    let other = other.response(&feeds, &gate).unwrap();
    for (&a, &b) in actual
        .interior
        .amplitudes
        .iter()
        .zip(&other.interior.amplitudes)
    {
        near(C64::new(a[0], a[1]), C64::new(b[0], b[1]), 1.0);
    }
    // A wrong closure can still return a finite native FE field. Admission
    // must check the ORIGINAL prescribed flows and interface pressures.
    let mut identity = vec![C64::ZERO; 25];
    for i in 0..5 {
        identity[i * 5 + i] = C64::ONE;
    }
    prepared.closure = lu_complex_with_cancel(&identity, 5, || false).unwrap();
    assert!(matches!(
        prepared.response(&feeds, &gate),
        Err(RadiationError::Invalid(_))
    ));
}

#[test]
fn signed_original_traces_survive_zero_net_flow_and_exact_geometry_binding() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &[2, 3, 4],
        &exterior,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    let expected = source_apertures(&fe, &ports, &exterior, &[2, 3, 4]);
    for (a, b) in model.apertures().iter().zip(&expected) {
        for (&x, &y) in a.shape.iter().zip(&b.shape) {
            assert!((x - y).abs() < 1e-12 * y.abs().max(1.0));
        }
    }
    let net: Vec<f64> = model
        .apertures()
        .iter()
        .map(|a| {
            a.shape
                .iter()
                .zip(exterior.areas())
                .map(|(v, area)| v * area)
                .sum()
        })
        .collect();
    assert!((net[0] - WIDTH * WIDTH).abs() < 1e-16);
    assert!(net[1].abs() < 1e-14);
    assert!(model.apertures()[1].shape.iter().any(|v| *v < 0.0));
    assert!(model.apertures()[1].shape.iter().any(|v| *v > 0.0));

    let moved = SpherePanels::from_triangles(
        exterior
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| triangle.map(|p| [p[0] + 0.001, p[1], p[2]]))
            .collect(),
    )
    .unwrap();
    assert!(matches!(
        RadiatingInterior::new(&fe, &ports, &[2, 3, 4], &moved, &gas, None, budget(), &gate),
        Err(RadiationError::Invalid(_))
    ));
    // A prescribed interior feed on a radiating face would omit its exterior
    // self/mutual work, so selecting only one of the shared-face modes refuses.
    assert!(matches!(
        RadiatingInterior::new(&fe, &ports, &[2, 4], &exterior, &gas, None, budget(), &gate),
        Err(RadiationError::Invalid(_))
    ));
    let unresolved = [
        ports[0].clone(),
        BoundaryMode {
            traces: vec![(ports[2].traces[0].0, [1.0, -1.0, 0.0])],
        },
    ];
    assert!(matches!(
        RadiatingInterior::new(
            &fe,
            &unresolved,
            &[1],
            &exterior,
            &gas,
            None,
            budget(),
            &gate
        ),
        Err(RadiationError::Invalid(_))
    ));
}

#[test]
fn original_full_field_wall_loss_balances_source_and_radiated_work() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let port_faces: BTreeSet<_> = ports
        .iter()
        .flat_map(|p| p.traces.iter().map(|(f, _)| face_key(*f)))
        .collect();
    let wall_faces: Vec<_> = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| !port_faces.contains(&face_key(*f)))
        .collect();
    let geometry = BoundaryLayerGeometry {
        minimum_curvature_radius_m: f64::INFINITY,
        minimum_half_gap_m: WIDTH / 2.0,
        maximum_layer_ratio: 0.05,
    };
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &[2, 3, 4],
        &exterior,
        &gas,
        Some((&wall_faces, geometry)),
        budget(),
        &gate,
    )
    .unwrap();
    let reference = gas.density * gas.sound_speed / (WIDTH * WIDTH);
    let prepared = model
        .prepare(FREQUENCY, &[reference; 5], &route(), &gate)
        .unwrap();
    let actual = prepared
        .response(&[C64::new(1e-6, 0.2e-6), C64::new(0.3e-6, -0.4e-6)], &gate)
        .unwrap();
    let original = model
        .interior()
        .boundary_layer_power(FREQUENCY, &actual.interior.amplitudes)
        .unwrap();
    assert_eq!(actual.boundary_layer_power, original);
    assert!(original.viscous_w > 0.0 && original.thermal_w > 0.0);
    let loss = original.viscous_w + original.thermal_w;
    let error = actual.source_power_w - actual.radiation_power_w - loss;
    assert!(
        error.abs() < 2e-8 * actual.source_power_w,
        "power defect {error:e}"
    );
    assert!((actual.interior.input_power_w() - loss).abs() < 2e-8 * actual.source_power_w);
    assert!(actual.port_relative_residual < 1e-10);
    assert!(actual.power_balance_residual_w.abs() < 2e-8 * actual.source_power_w);
    assert!(
        prepared
            .response(&[C64::new(f64::NAN, 0.0), C64::ZERO], &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        prepared.response(&[C64::ZERO; 2], &gate),
        Err(RadiationError::Cancelled)
    ));
}

#[test]
fn full_pressure_cap_and_source_gas_are_enforced_before_binding() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let gate = CancelGate::new();
    let mut limit = budget();
    limit.maximum_pressure_dofs = fe.ndof() - 1;
    assert!(matches!(
        RadiatingInterior::new(&fe, &ports, &[2, 3, 4], &exterior, &gas, None, limit, &gate),
        Err(RadiationError::Interior(_))
    ));
    let warmer = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.15, 101325.0).unwrap();
    assert!(matches!(
        RadiatingInterior::new(
            &fe,
            &ports,
            &[2, 3, 4],
            &exterior,
            &warmer,
            None,
            budget(),
            &gate
        ),
        Err(RadiationError::Invalid(_))
    ));
}
