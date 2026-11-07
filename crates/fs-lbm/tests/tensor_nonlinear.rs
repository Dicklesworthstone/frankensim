//! Nonlinear full-tensor conduction through the public native owner.
//!
//! These cases exercise original interval laws, actual surface temperatures
//! and local/global heat accounting. They make no discrete maximum-principle
//! or general nonlinear convergence claim for rotated anisotropy.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CellSink, ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance,
    RotatedSolidConductivity, SolidMaterial, TensorConductivityLaw, TensorNonlinearConfig,
    ThermalFace, ThermalSetup, Voxel, VoxelDomain, solve_nonlinear_tensor_conduction,
    solve_tensor_conduction,
};

const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn center(domain: &VoxelDomain, cell: usize) -> [f64; 3] {
    let [x, y, z] = domain.coords(cell);
    domain.center(x, y, z)
}

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        ..TensorNonlinearConfig::default()
    }
}

#[test]
fn nonlinear_transformed_affine_patch_keeps_transverse_heat_on_graded_cells() {
    let material = [SolidMaterial::new("original linear-in-T law", 1.0)
        .with_orthotropic([6.0, 2.0, 3.0])
        .with_conductivity_table(&[(200.0, 2.0 / 3.0), (600.0, 2.0)])];
    let law = RotatedSolidConductivity::new(&material[0], ROTATION).unwrap();
    let base = law.conductivity_at_temperature(300.0).unwrap().as_array();
    let gradient_u = [0.0, 30.0, -10.0];
    let flux = base.map(|row| -dot(row, gradient_u));
    // U(T)=T^2/600. An affine U gives an exact continuum heat flux,
    // including nonzero x heat despite zero x temperature gradient.
    let temperature = |point| (350.0_f64.powi(2) + 600.0 * dot(point, gradient_u)).sqrt();
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.2, 0.35, 0.6], vec![0.4, 0.7], vec![0.3, 0.8]],
        |_| Voxel::Solid(0),
    )
    .unwrap();
    let faces = std::array::from_fn(|f| {
        ThermalFace::HeatFlux(if f & 1 == 0 {
            flux[f / 2]
        } else {
            -flux[f / 2]
        })
    });
    let mut setup = ThermalSetup::new(faces);
    setup
        .fixed_temperature
        .push((0, temperature(center(&domain, 0))));
    let laws: Vec<&dyn TensorConductivityLaw> = vec![&law; domain.cell_count()];
    let config = controls();
    let result = solve_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &laws,
        &vec![350.0; domain.cell_count()],
        &config,
        &CancelGate::new(),
    )
    .unwrap();
    let solved = &result.solution;
    for c in 0..domain.cell_count() {
        assert!((solved.temperature[c] - temperature(center(&domain, c))).abs() < 3e-7);
        for f in 0..6 {
            let a = f / 2;
            let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
            let mut point = center(&domain, c);
            point[a] += 0.5 * sign * domain.widths(c)[a];
            let expected = sign * domain.face_area(c, a) * flux[a];
            assert!((solved.face_temperature_k[c][f] - temperature(point)).abs() < 3e-7);
            assert!((solved.outward_heat_w[c][f] - expected).abs() < 3e-7);
            if let Some(n) = domain.neighbor(c, f) {
                assert_eq!(
                    solved.face_temperature_k[c][f],
                    solved.face_temperature_k[n][f ^ 1]
                );
                assert!(
                    (solved.outward_heat_w[c][f] + solved.outward_heat_w[n][f ^ 1]).abs() < 3e-7
                );
            }
        }
    }
    assert!((solved.temperature[0] - solved.temperature[1]).abs() < 3e-7);
    assert!(solved.outward_heat_w[0][1].abs() > 1.0);
    assert!(result.newton_iterations > 0);
    assert!(result.newton_iterations <= config.max_newton_iterations);
    assert!(solved.report.max_heat_residual <= config.linear.tolerance);
    assert!(solved.report.balance.relative_residual < 2e-9);
}

#[test]
fn independent_rotated_laws_resolve_finite_contact_surfaces_and_close_heat() {
    let material = [
        SolidMaterial::new("anisotropic A", 1.0).with_axis_conductivity_tables([
            &[(280.0, 2.0), (350.0, 8.0), (460.0, 5.0)],
            &[(260.0, 4.0), (390.0, 2.0), (500.0, 9.0)],
            &[(300.0, 3.0), (450.0, 4.0)],
        ]),
        SolidMaterial::new("anisotropic B", 1.0).with_axis_conductivity_tables([
            &[(250.0, 1.0), (340.0, 3.0), (500.0, 10.0)],
            &[(280.0, 7.0), (420.0, 3.0)],
            &[(280.0, 4.0), (460.0, 8.0)],
        ]),
    ];
    let domain = VoxelDomain::graded_from_fn([vec![0.3, 0.7], vec![0.2, 0.5], vec![0.4]], |x| {
        Voxel::Solid(if x[0] < 0.3 { 0 } else { 1 })
    })
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(430.0);
    faces[1] = ThermalFace::Temperature(290.0);
    let mut setup = ThermalSetup::new(faces);
    let resistance = 0.025;
    for c in 0..domain.cell_count() {
        if domain.coords(c)[0] == 0 {
            setup.face_contacts.push(FaceContactResistance {
                cells: [c, c + 1],
                resistance_m2_k_w: resistance,
            });
        }
    }
    let config = controls();
    let solve = |rotation| {
        let adapted = [
            RotatedSolidConductivity::new(&material[0], rotation).unwrap(),
            RotatedSolidConductivity::new(&material[1], rotation).unwrap(),
        ];
        let laws: Vec<&dyn TensorConductivityLaw> = (0..domain.cell_count())
            .map(|c| &adapted[domain.coords(c)[0]] as &dyn TensorConductivityLaw)
            .collect();
        solve_nonlinear_tensor_conduction(
            &domain,
            &material,
            &setup,
            &laws,
            &vec![360.0; domain.cell_count()],
            &config,
            &CancelGate::new(),
        )
        .unwrap()
    };
    let result = solve(ROTATION);
    let aligned = solve(IDENTITY);
    let solved = result.solution;
    let mut total_contact_heat = 0.0;
    for c in 0..domain.cell_count() {
        assert!(solved.outward_heat_w[c].iter().sum::<f64>().abs() < 3e-7);
        if domain.coords(c)[0] == 0 {
            let q = solved.outward_heat_w[c][1];
            let other = solved.outward_heat_w[c + 1][0];
            let drop = solved.face_temperature_k[c][1] - solved.face_temperature_k[c + 1][0];
            let expected = domain.face_area(c, 0) / resistance * drop;
            assert!((q + other).abs() < 3e-7);
            assert!((q - expected).abs() < 3e-7);
            total_contact_heat += q;
        }
    }
    assert!(total_contact_heat > 0.0);
    let maximum_orientation_change = solved
        .temperature
        .iter()
        .zip(&aligned.solution.temperature)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(maximum_orientation_change > 1e-3);
    let law = RotatedSolidConductivity::new(&material[0], ROTATION).unwrap();
    let cold = law.conductivity_at_temperature(300.0).unwrap().as_array();
    let hot = law.conductivity_at_temperature(420.0).unwrap().as_array();
    assert!((cold[0][1] - hot[0][1]).abs() > 0.01);
    assert!(solved.report.max_heat_residual <= config.linear.tolerance);
    assert!(solved.report.balance.residual_w.abs() < 3e-7);
    assert!(solved.report.balance.relative_residual < 2e-9);
}

#[test]
fn constant_laws_match_frozen_solver_with_contacts_sources_sink_and_film() {
    let domain = VoxelDomain::graded_from_fn([vec![0.1, 0.3, 0.6], vec![0.2], vec![0.3]], |_| {
        Voxel::Solid(0)
    })
    .unwrap();
    let material = [SolidMaterial::new("constant owner", 1.0)];
    let tensor =
        ConductivityTensor::new([[5.0, 1.0, 0.3], [1.0, 3.0, 0.2], [0.3, 0.2, 2.0]]).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(400.0);
    faces[1] = ThermalFace::Convective {
        h: 8.0,
        ambient: 330.0,
    };
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![1.0, 0.5, 2.0];
    setup.cell_sinks.push(CellSink {
        cell: 1,
        conductance_w_k: 0.4,
        temperature_k: 310.0,
    });
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.05,
    });
    let config = controls();
    let frozen = solve_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[tensor; 3],
        &config.linear,
        &CancelGate::new(),
    )
    .unwrap();
    let laws: [&dyn TensorConductivityLaw; 3] = [&tensor; 3];
    let nonlinear = solve_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &laws,
        &[350.0; 3],
        &config,
        &CancelGate::new(),
    )
    .unwrap()
    .solution;
    for c in 0..3 {
        assert!((nonlinear.temperature[c] - frozen.temperature[c]).abs() < 3e-7);
        for f in 0..6 {
            assert!(
                (nonlinear.face_temperature_k[c][f] - frozen.face_temperature_k[c][f]).abs() < 3e-7
            );
            assert!((nonlinear.outward_heat_w[c][f] - frozen.outward_heat_w[c][f]).abs() < 3e-7);
        }
    }
    assert!(
        (nonlinear.report.balance.sink_outflow_w - frozen.report.balance.sink_outflow_w).abs()
            < 3e-7
    );
    let film_heat = 8.0 * domain.face_area(2, 0) * (nonlinear.face_temperature_k[2][1] - 330.0);
    assert!((nonlinear.outward_heat_w[2][1] - film_heat).abs() < 3e-7);
    assert!(nonlinear.report.max_heat_residual <= config.linear.tolerance);
}

struct GappedSource {
    tensor: ConductivityTensor,
}

impl TensorConductivityLaw for GappedSource {
    fn conductivity_at_temperature(&self, t: f64) -> Result<ConductivityTensor, ChtError> {
        if t.is_finite() && ((250.0..=340.0).contains(&t) || (360.0..=500.0).contains(&t)) {
            Ok(self.tensor)
        } else {
            Err(ChtError::InvalidInput {
                field: "test.original_tensor_source",
                reason: "the original source has no support here".into(),
            })
        }
    }

    fn mean_conductivity(&self, a: f64, b: f64) -> Result<ConductivityTensor, ChtError> {
        self.conductivity_at_temperature(a)?;
        self.conductivity_at_temperature(b)?;
        if a.min(b) < 360.0 && a.max(b) > 340.0 {
            return Err(ChtError::InvalidInput {
                field: "test.original_tensor_source",
                reason: "the closed temperature interval crosses a source-support gap".into(),
            });
        }
        Ok(self.tensor)
    }
}

#[test]
fn refuses_unfunded_work_invalid_states_source_gaps_and_cancellation() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [SolidMaterial::new("nonlinear owner", 1.0)
        .with_conductivity_table(&[(250.0, 1.0), (450.0, 3.0)])];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let laws: [&dyn TensorConductivityLaw; 2] = [&law; 2];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.faces[0] = ThermalFace::Temperature(300.0);
    setup.faces[1] = ThermalFace::Temperature(400.0);
    let config = controls();
    let solve = |laws: &[&dyn TensorConductivityLaw],
                 initial: &[f64],
                 config: &TensorNonlinearConfig,
                 gate: &CancelGate| {
        solve_nonlinear_tensor_conduction(&domain, &material, &setup, laws, initial, config, gate)
    };
    let no_newton = TensorNonlinearConfig {
        max_newton_iterations: 0,
        ..controls()
    };
    assert!(solve(&laws, &[320.0; 2], &no_newton, &CancelGate::new()).is_err());
    let no_linear = TensorNonlinearConfig {
        linear: EnergyConfig {
            max_iterations: 0,
            ..config.linear
        },
        ..controls()
    };
    assert!(solve(&laws, &[320.0; 2], &no_linear, &CancelGate::new()).is_err());
    assert!(solve(&laws[..1], &[320.0; 2], &config, &CancelGate::new()).is_err());
    assert!(solve(&laws, &[320.0], &config, &CancelGate::new()).is_err());
    assert!(solve(&laws, &[0.0, 320.0], &config, &CancelGate::new()).is_err());
    let gapped = GappedSource {
        tensor: ConductivityTensor::new([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
            .unwrap(),
    };
    // Both boundary endpoints are supported; the intervening original law is not.
    assert!(gapped.conductivity_at_temperature(300.0).is_ok());
    assert!(gapped.conductivity_at_temperature(400.0).is_ok());
    assert!(matches!(
        solve(
            &[&gapped, &gapped],
            &[320.0; 2],
            &config,
            &CancelGate::new()
        ),
        Err(ChtError::InvalidInput { .. }),
    ));
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        solve(&laws, &[320.0; 2], &config, &gate).unwrap_err(),
        ChtError::Cancelled
    );
}
