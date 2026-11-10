//! Original nonlinear source/current, rigid midpoint and immutable mass oracles.
use super::*;
use fs_couple::circuit_enthalpy::conduction::spatial::mechanical::response::{
    SpatialConductionResponseKickControls, SpatialConductionResponseKickError,
    kick_response_candidate, kick_response_with_materials_candidate,
};
use fs_em::spatial::conduction::response::{
    ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits,
    SpatialConductionResponse, SpatialConductionResponseControls,
};
use fs_material::conductor::response::{ConductorMagnitude, ScalarConductorResponse};
use fs_qty::{QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;

#[derive(Clone, Copy)]
enum DensityLaw {
    Constant,
    Temperature,
    ElectricalValue,
    ElectricalApplicability,
}

fn context(temperature: f64, coordinate: ConductorMagnitude) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap()
}

fn nonlinear_fixture(
    coordinate: ConductorMagnitude,
    density_law: DensityLaw,
) -> (SpatialConductionResponse, ResolvedBody) {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut claims = ClaimSet::new();
    let temperature_domain = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        800.0,
    );
    let (density_value, density_domain, density_interpolation) = match density_law {
        DensityLaw::Constant => (
            PropertyValue::Scalar {
                value: 3.0,
                dims: Density::DIMS,
            },
            ValidityDomain::unconstrained(),
            InterpolationPolicy::ConstantWithinValidity,
        ),
        DensityLaw::Temperature => (
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 2.9), (800.0, 3.5)],
                dims: Density::DIMS,
            },
            temperature_domain.clone(),
            InterpolationPolicy::LinearInside,
        ),
        DensityLaw::ElectricalValue => (
            PropertyValue::Curve {
                abscissa: "operating".into(),
                abscissa_dims: coordinate.quantity().dims(),
                knots: vec![(0.0, 3.0), (20.0, 3.4)],
                dims: Density::DIMS,
            },
            ValidityDomain::unconstrained().with_quantity(
                "operating",
                coordinate.quantity(),
                0.0,
                20.0,
            ),
            InterpolationPolicy::LinearInside,
        ),
        DensityLaw::ElectricalApplicability => (
            PropertyValue::Scalar {
                value: 3.0,
                dims: Density::DIMS,
            },
            ValidityDomain::unconstrained().with_quantity(
                "operating",
                coordinate.quantity(),
                0.0,
                20.0,
            ),
            InterpolationPolicy::ConstantWithinValidity,
        ),
    };
    for (key, value, validity, interpolation) in [
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            PropertyValue::Surface2D {
                axes: ["T".into(), "operating".into()],
                axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                samples: [[200.0, 0.0], [800.0, 0.0], [200.0, 20.0], [800.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, 1.0 + 0.005 * (p[0] - 300.0) + 0.5 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
            temperature_domain.with_quantity("operating", coordinate.quantity(), 0.0, 20.0),
            InterpolationPolicy::LinearInside,
        ),
        (
            PropertyKey::new("density", Density::DIMS),
            density_value,
            density_domain,
            density_interpolation,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured nonlinear original-body momentum oracle".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured nonlinear conductor".into(),
            phase: "solid".into(),
            process: "fixed body momentum oracle".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let point = context(300.0, coordinate);
    let original = resolve_conductor_state_with_basis(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        basis,
    )
    .unwrap();
    let source =
        ScalarConductorResponse::new(&card, original, &[("operating", coordinate)]).unwrap();
    let (geometry, _) = fixture(1.0, 3.0, &[-0.5, 0.5]);
    let positions = geometry.positions_m().to_vec();
    let tets = geometry.mesh().tets.clone();
    let problem =
        SpatialConductionResponse::new(positions, tets, vec![source; 6], assembly(), &mut || false)
            .unwrap()
            .with_relaxation(0.5)
            .unwrap();
    let keys: Vec<_> = (0..6).map(|i| format!("tet-{i:04}")).collect();
    let regions: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .map(|(i, tet)| RegionBinding {
            key: &keys[i],
            material: &card,
            point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::tetrahedron(tet.map(|v| {
                let p = problem.positions_m()[v as usize];
                Vec3::new(p[0], p[1], p[2])
            }))
            .unwrap(),
        })
        .collect();
    (problem, resolve_body(&regions, 6, || false).unwrap())
}

fn mass<'a>(problem: &SpatialConductionResponse, body: &'a ResolvedBody) -> TetrahedralMass<'a> {
    let keys: Vec<_> = (0..6).map(|i| format!("tet-{i:04}")).collect();
    let names: Vec<_> = keys.iter().map(String::as_str).collect();
    let tets: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .map(|t| t.map(|v| v as usize))
        .collect();
    TetrahedralMass::bind(body, problem.positions_m(), &tets, &names, 100, || false).unwrap()
}

fn faces(problem: &SpatialConductionResponse) -> [Vec<usize>; 2] {
    [-0.5, 0.5].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect()
    })
}

fn boundary(problem: &SpatialConductionResponse, voltage: f64) -> SpatialConductionBoundary {
    SpatialConductionBoundary {
        voltage: faces(problem)
            .iter()
            .zip([0.0, voltage])
            .flat_map(|(faces, value)| {
                faces.iter().map(move |&face| ConductionVoltageFace {
                    face,
                    voltage_v: [value; 3],
                })
            })
            .collect(),
        current: Vec::new(),
    }
}

fn nonlinear_controls() -> SpatialConductionResponseKickControls {
    SpatialConductionResponseKickControls {
        physical: SpatialConductionKickControls {
            maximum_field_solves: 20_000,
            ..controls()
        },
        response: SpatialConductionResponseControls {
            maximum_iterations: 100,
            maximum_field_solves: 100,
            absolute_conductivity_tolerance_s_m: 1e-12,
            relative_conductivity_tolerance: 1e-12,
            current_density_tolerance_a_m2: 1e-11,
        },
        maximum_newton_iterations: 20,
        angular_velocity_scale_rad_s: 1.0,
        linear_velocity_scale_m_s: 1.0,
        difference_step: 1e-4,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            linear_restart: 6,
            ..NewtonKrylovConfig::default()
        },
    }
}

// Independent constitutive algebra: rho=r+|E'|/2, or rho=r+|J|/2.
// For the latter, use the stable root of (r+|J|/2)J=E'.
fn current(coordinate: ConductorMagnitude, r: f64, field: f64) -> f64 {
    match coordinate {
        ConductorMagnitude::ElectricField => field / (r + 0.5 * field.abs()),
        ConductorMagnitude::CurrentDensity => {
            2.0 * field / (r + (r * r + 2.0 * field.abs()).sqrt())
        }
    }
}

fn bar_midpoint(
    coordinate: ConductorMagnitude,
    r: f64,
    old: f64,
    voltage: f64,
    external: f64,
    dt: f64,
) -> (f64, f64) {
    // Monotone scalar momentum balance, independent of numerical Jacobians,
    // FE assembly and the production native impulse implementation.
    let residual = |mid| {
        mid - old - dt / 6.0 * (-1.5 * current(coordinate, r, 1.5 * mid - voltage) + external)
    };
    let mut low = -20.0;
    let mut high = 20.0;
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let mid = 0.5 * (low + high);
    (mid, current(coordinate, r, 1.5 * mid - voltage))
}

#[test]
fn nonlinear_braking_generation_and_motoring_close_original_current_and_native_storage() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let (problem, body) = nonlinear_fixture(coordinate, DensityLaw::Constant);
        let mass = mass(&problem, &body);
        for (voltage, external) in [(0.0, 0.0), (0.4, 0.0), (2.4, 0.2)] {
            let boundary = boundary(&problem, voltage);
            let dt = 0.3;
            let candidate = kick_response_candidate(
                &problem,
                &mass,
                &boundary,
                &[[0.0, 0.0, 1.5]; 6],
                dt,
                translating(1.0),
                Wrench::new(Vec3::default(), Vec3::new(0.0, external, 0.0)),
                assembly(),
                nonlinear_controls(),
                &mut || false,
            )
            .unwrap();
            let report = candidate.report();
            let (mid, j) = bar_midpoint(coordinate, 1.0, 1.0, voltage, external, dt);
            let after = 2.0 * mid - 1.0;
            near(report.midpoint_twist_body.linear.y, mid);
            near(report.mechanical.base_after.twist_body.linear.y, after);
            near(report.lorentz_wrench_body.force.y, -1.5 * j);
            near(report.joule_heat_j, dt * (1.5 * mid - voltage) * j);
            near(report.terminal_work_j, -dt * voltage * j);
            near(report.mechanical_work_j, -dt * 1.5 * mid * j);
            near(report.external_work_j, dt * external * mid);
            near(
                report.mechanical.kinetic_energy_change,
                1.5 * (after - 1.0) * (after + 1.0),
            );
            near(report.mechanical_work_residual_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            for value in report.field.field.current_density_a_m2() {
                near(value[0], j);
                near(value[1], 0.0);
                near(value[2], 0.0);
            }
            assert!(report.field_solves > report.field.field_solves);
            assert!(report.source_evaluations > report.field.source_evaluations);
            assert!(report.field_solves <= nonlinear_controls().physical.maximum_field_solves);
            assert_eq!(
                report.mechanical.base_after.world_from_base,
                Se3::identity()
            );
            if voltage == 0.4 {
                assert!(report.terminal_work_j < 0.0);
            }
            if voltage == 2.4 {
                assert!(report.mechanical_work_j > 0.0);
            }
        }
        assert_eq!(
            problem.conductors()[0].material().query().axes()["operating"],
            0.0
        );
    }
}

#[test]
fn original_total_current_terminal_drive_continues_with_changed_electrical_coordinates() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let (mut problem, body) = nonlinear_fixture(coordinate, DensityLaw::Constant);
    let mass = mass(&problem, &body);
    let faces = faces(&problem);
    let terminal_faces = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let currents = [-1.25, 1.25];
    let drive = ConductionTerminalDrive {
        terminals: &terminal_faces,
        control: ConductionTerminalControl::Currents(&currents),
        limits: ConductionTerminalLimits {
            maximum_terminals: 2,
            maximum_workspace_values: 20_000,
        },
    };
    let mut state = translating(-0.4);
    for _ in 0..2 {
        let old = state.twist_body.linear.y;
        let candidate = kick_response_candidate(
            &problem,
            &mass,
            drive,
            &[[0.0, 0.0, 1.5]; 6],
            0.2,
            state,
            Wrench::default(),
            assembly(),
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
        let (next, report) = candidate.into_parts();
        let after = old + 0.2 * 1.5 * 1.25 / 3.0;
        let mid = 0.5 * (old + after);
        let voltage = (1.0 + 0.5 * 1.25) * 1.25 + 1.5 * mid;
        near(report.mechanical.base_after.twist_body.linear.y, after);
        near(report.terminal_work_j, 0.2 * 1.25 * voltage);
        near(report.joule_heat_j, 0.2 * (1.0 + 0.5 * 1.25) * 1.25 * 1.25);
        let terminal = report.field.terminals.as_ref().unwrap();
        for (&a, &b) in terminal.current_a.iter().zip(&currents) {
            near(a, b);
        }
        near(report.terminal_work_j, 0.2 * terminal.supplied_power_w);
        near(report.total_energy_residual_j, 0.0);
        assert_eq!(
            next.conductors()[0].material().card_identity(),
            problem.conductors()[0].material().card_identity()
        );
        near(
            next.conductors()[0].material().query().axes()["operating"],
            1.25,
        );
        state = report.mechanical.base_after;
        problem = next;
    }
    // Same original density query remains at its initial electrical coordinate.
    assert_eq!(
        mass.reference().regions()[0].material().query().axes()["operating"],
        0.0
    );
}

#[test]
fn original_density_pins_admit_external_temperature_refresh_and_refuse_field_dependence() {
    let coordinate = ConductorMagnitude::ElectricField;
    for law in [
        DensityLaw::Constant,
        DensityLaw::Temperature,
        DensityLaw::ElectricalValue,
        DensityLaw::ElectricalApplicability,
    ] {
        let (problem, body) = nonlinear_fixture(coordinate, law);
        let mass = mass(&problem, &body);
        let boundary = boundary(&problem, 0.0);
        let seeds = problem
            .sources()
            .iter()
            .map(|source| source.at_query_point(&context(600.0, coordinate)).unwrap())
            .collect();
        let attempted = kick_response_with_materials_candidate(
            &problem,
            seeds,
            &mass,
            &boundary,
            &[[0.0, 0.0, 1.5]; 6],
            0.2,
            translating(1.0),
            Wrench::default(),
            assembly(),
            nonlinear_controls(),
            &mut || false,
        );
        match law {
            DensityLaw::Constant => {
                let candidate = attempted.unwrap();
                let (mid, _) = bar_midpoint(coordinate, 2.5, 1.0, 0.0, 0.0, 0.2);
                near(candidate.report().midpoint_twist_body.linear.y, mid);
                assert_eq!(
                    candidate.problem().conductors()[0]
                        .material()
                        .query()
                        .axes()["T"],
                    600.0
                );
            }
            _ => assert!(matches!(
                attempted,
                Err(SpatialConductionResponseKickError::Kick(
                    SpatialConductionKickError::Invalid(_)
                ))
            )),
        }
        assert_eq!(
            problem.conductors()[0].material().query().axes()["T"],
            300.0
        );
        assert_eq!(
            mass.reference().regions()[0].material().query().axes()["T"],
            300.0
        );
    }
}

#[test]
fn full_rigid_nonlinear_feedback_matches_independent_original_cell_moments() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (problem, body) = nonlinear_fixture(coordinate, DensityLaw::Constant);
    let mass = mass(&problem, &body);
    let boundary = SpatialConductionBoundary {
        voltage: problem
            .exterior_faces()
            .map(|face| ConductionVoltageFace {
                face,
                voltage_v: [0.0; 3],
            })
            .collect(),
        current: Vec::new(),
    };
    let old = [0.7, -0.3, 0.2, 0.1, 1.0, -0.2];
    let initial = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(
            Vec3::new(old[0], old[1], old[2]),
            Vec3::new(old[3], old[4], old[5]),
        ),
    );
    let flux = Vec3::new(0.2, -0.4, 1.3);
    let dt = 0.2;
    let centroids: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .map(|tet| {
            let mut c = Vec3::default();
            for &v in tet {
                let p = problem.positions_m()[v as usize];
                c = c + Vec3::new(p[0] * 0.25, p[1] * 0.25, p[2] * 0.25);
            }
            c
        })
        .collect();
    // At zero potential on every original vertex E=0. Analytic rho(|E'|)
    // and direct cell moments give the six-axis nonlinear wrench. Cube mass=3,
    // principal inertia=1/2; this fixed-point oracle does not use FE or impulses.
    let force = |mid: [f64; 6]| {
        let angular = Vec3::new(mid[0], mid[1], mid[2]);
        let linear = Vec3::new(mid[3], mid[4], mid[5]);
        let mut torque = Vec3::default();
        let mut force = Vec3::default();
        for &c in &centroids {
            let e = (linear + angular.cross(c)).cross(flux);
            let norm = (e.x * e.x + e.y * e.y + e.z * e.z).sqrt();
            let j = Vec3::new(e.x, e.y, e.z).scale(1.0 / (1.0 + 0.5 * norm));
            let cell_force = j.cross(flux).scale(1.0 / 6.0);
            torque = torque + c.cross(cell_force);
            force = force + cell_force;
        }
        [torque.x, torque.y, torque.z, force.x, force.y, force.z]
    };
    let mut mid = old;
    for _ in 0..100 {
        let wrench = force(mid);
        mid =
            core::array::from_fn(|i| old[i] + 0.5 * dt * wrench[i] / if i < 3 { 0.5 } else { 3.0 });
    }
    let candidate = kick_response_candidate(
        &problem,
        &mass,
        &boundary,
        &[[flux.x, flux.y, flux.z]; 6],
        dt,
        initial,
        Wrench::default(),
        assembly(),
        nonlinear_controls(),
        &mut || false,
    )
    .unwrap();
    let report = candidate.report();
    for (actual, expected) in report.midpoint_twist_body.to_array().into_iter().zip(mid) {
        near(actual, expected);
    }
    for (actual, expected) in report
        .lorentz_wrench_body
        .to_array()
        .into_iter()
        .zip(force(mid))
    {
        near(actual, expected);
    }
    near(report.terminal_work_j, 0.0);
    near(report.mechanical_work_j, -report.joule_heat_j);
    near(
        report.mechanical.kinetic_energy_change,
        -report.joule_heat_j,
    );
    near(report.total_energy_residual_j, 0.0);
}

#[test]
fn budget_and_cancellation_refuse_atomically_and_do_not_poison_retry() {
    let (problem, body) =
        nonlinear_fixture(ConductorMagnitude::ElectricField, DensityLaw::Constant);
    let mass = mass(&problem, &body);
    let boundary = boundary(&problem, 0.0);
    let run = |controls, cancel: &mut dyn FnMut() -> bool| {
        kick_response_candidate(
            &problem,
            &mass,
            &boundary,
            &[[0.0, 0.0, 1.5]; 6],
            0.2,
            translating(1.0),
            Wrench::default(),
            assembly(),
            controls,
            &mut || cancel(),
        )
    };
    let mut bounded = nonlinear_controls();
    bounded.physical.maximum_field_solves = 1;
    assert!(matches!(
        run(bounded, &mut || false),
        Err(SpatialConductionResponseKickError::FieldSolveBudget { limit: 1, .. })
    ));
    assert!(matches!(
        run(nonlinear_controls(), &mut || true),
        Err(SpatialConductionResponseKickError::Cancelled)
    ));
    let mut calls = 0;
    let reference = run(nonlinear_controls(), &mut || {
        calls += 1;
        false
    })
    .unwrap();
    let cutoff = calls / 2;
    let mut interrupted = 0;
    assert!(
        run(nonlinear_controls(), &mut || {
            interrupted += 1;
            interrupted >= cutoff
        })
        .is_err()
    );
    let retry = run(nonlinear_controls(), &mut || false).unwrap();
    assert_eq!(retry.report().mechanical, reference.report().mechanical);
    assert_eq!(
        retry.report().field.field.current_density_a_m2(),
        reference.report().field.field.current_density_a_m2()
    );
    assert_eq!(
        problem.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}

#[test]
fn owned_nonlinear_load_handoff_reuses_accepted_current_and_checks_local_velocity() {
    use fs_em::spatial::conduction::motion::{SpatialConductionMotion, SpatialLorentzLoadError};
    let (problem, _) = nonlinear_fixture(ConductorMagnitude::ElectricField, DensityLaw::Constant);
    let boundary = boundary(&problem, 0.0);
    let candidate = problem
        .solve_with_motion_candidate(
            &boundary,
            SpatialConductionMotion {
                velocity_m_s: &[[0.0, 1.0, 0.0]; 6],
                flux_density_t: &[[0.0, 0.0, 1.5]; 6],
            },
            assembly(),
            controls().field,
            nonlinear_controls().response,
            &mut || false,
        )
        .unwrap();
    let solves = candidate.report().field_solves;
    let selected = [5, 3, 1, 4, 2, 0];
    let mut nodal = vec![[0.0, 1.0, 0.0]; problem.positions_m().len()];
    let loads = candidate
        .lorentz_loads_for_velocity(&selected, &nodal, 0.0, controls().loads, &mut || false)
        .unwrap();
    let wrench = loads.wrench(&Se3::identity(), [0.0; 3]).unwrap();
    near(
        wrench.pairing(translating(1.0).twist_body),
        candidate.report().field.mechanical_power_w(),
    );
    near(
        wrench.force.y,
        -1.5 * current(ConductorMagnitude::ElectricField, 1.0, 1.5),
    );
    nodal[0][1] += 0.1;
    assert!(matches!(
        candidate
            .lorentz_loads_for_velocity(&selected, &nodal, 0.0, controls().loads, &mut || false,),
        Err(SpatialLorentzLoadError::Invalid(_))
    ));
    assert!(
        candidate
            .lorentz_loads(&[0, 0], controls().loads, &mut || false,)
            .is_err()
    );
    assert_eq!(candidate.report().field_solves, solves);
}
