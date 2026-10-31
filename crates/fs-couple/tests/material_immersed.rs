//! G1/G3/G4: complete source substitutions drive real native fluid/body steps.
//! Synthetic cards exercise the implementation, not measured-material accuracy.
use fs_couple::material_immersed::{
    AdditionalMaterialSphere, FilmGeometry, ImmersedDomain, MaterialImmersedError,
    RetainedFluidSpheres, SphereStateSi,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::{E3, PartialSaturationError3, Q3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};
use fs_time::So3;

const RADIUS: f64 = 0.006;
const DX: f64 = 0.01;
const DT: f64 = 0.1;

fn card(label: &str, properties: &[(&str, Dims, f64)], phase: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: if name == "thermal_conductivity" {
                    PropertyKey::with_quantity(
                        name,
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity,
                            ValueForm::Static,
                        )),
                    )
                } else {
                    PropertyKey::new(name, dims)
                },
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: format!("synthetic {label}"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: label.into(),
            phase: phase.into(),
            process: "synthetic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn body(density: f64) -> ResolvedBody {
    let card = card(
        "additional sphere",
        &[("density", Density::DIMS, density)],
        "solid",
    );
    resolve_body(
        &[RegionBinding {
            key: "sphere",
            geometry: VolumeMoments::sphere(RADIUS).unwrap(),
            material: &card,
            point: &point(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn fluid(density: f64, viscosity: f64) -> ResolvedLiquidState {
    let card = card(
        "fluid",
        &[
            ("density", Density::DIMS, density),
            ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, viscosity),
            (
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                1000.0,
            ),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        ],
        "liquid",
    );
    resolve_liquid_state(&card, &point(), MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn domain() -> ImmersedDomain {
    ImmersedDomain {
        cells: [8; 3],
        periodic: [true; 3],
        cell_size_m: DX,
        timestep_s: DT,
        maximum_speed_m_s: 0.001,
        maximum_cells: 512,
        samples_per_axis: 4,
    }
}
fn film() -> FilmGeometry {
    FilmGeometry {
        cutoff_gap_m: 0.003,
        roughness_gap_m: 1e-5,
        tangential_gaps_m: None,
        maximum_pairs: 1,
        maximum_resistance_modes: 8,
    }
}
fn initial() -> [SphereStateSi; 2] {
    [0.03, 0.043].map(|x| SphereStateSi {
        center_m: [x, 0.04, 0.04],
        orientation: So3::identity(),
        velocity_m_s: [if x < 0.04 { 0.0001 } else { -0.0001 }, 0.0, 0.0],
        angular_momentum_body_kg_m2_s: [0.0; 3],
    })
}
fn assembly<'a>(
    fluid: &'a ResolvedLiquidState,
    body: &'a ResolvedBody,
) -> RetainedFluidSpheres<'a, ResolvedLiquidState, 2> {
    RetainedFluidSpheres::new(
        fluid,
        domain(),
        [AdditionalMaterialSphere::new(body, RADIUS).unwrap(); 2],
        initial(),
        film(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; delta={:.4e}",
        a - b
    );
}
fn snapshot(model: &RetainedFluidSpheres<'_, ResolvedLiquidState, 2>) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                result.push(model.grid().populations(x, y, z));
            }
        }
    }
    result
}
fn momentum(model: &RetainedFluidSpheres<'_, ResolvedLiquidState, 2>) -> [f64; 3] {
    // Independently sum ACTUAL populations, including retained interior fluid,
    // and the source-derived additional carrier momentum. Do not subtract it.
    let mut total = [0.0; 3];
    for populations in snapshot(model) {
        for q in 0..Q3 {
            let e = [E3[q].0, E3[q].1, E3[q].2];
            for axis in 0..3 {
                total[axis] += populations[q] * f64::from(e[axis]);
            }
        }
    }
    for body in model.bodies_lattice() {
        for axis in 0..3 {
            total[axis] += body.mass * body.velocity[axis];
        }
    }
    total
}

#[test]
fn live_material_refresh_preserves_history_and_changes_bulk_and_film_together() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let carrier = body(3000.0);
    let gate = CancelGate::new();
    let mut model = assembly(&first, &carrier);
    let mut held = assembly(&first, &carrier);
    for _ in 0..2 {
        let a = model
            .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
            .unwrap();
        let b = held.step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate).unwrap();
        assert_eq!(a, b);
    }
    let populations = snapshot(&model);
    let bodies = *model.bodies_lattice();
    let before = momentum(&model);
    let mass_unit = model.mass_kg(1.0).unwrap();
    let energy_unit = model.energy_j(1.0).unwrap();
    model.rebind_fluid(&second, &gate).unwrap();
    assert!(core::ptr::eq(model.scaling().material(), &second));
    assert_eq!(model.accepted_steps(), 2);
    assert_eq!(snapshot(&model), populations);
    assert_eq!(model.bodies_lattice(), &bodies);
    assert_eq!(momentum(&model), before);
    assert_eq!(model.mass_kg(1.0).unwrap(), mass_unit);
    assert_eq!(model.energy_j(1.0).unwrap(), energy_unit);
    close(
        model.grid().viscosity(),
        0.22 / 1000.0 * DT / DX.powi(2),
        1e-16,
    );

    let report = model
        .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
        .unwrap();
    held.step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate).unwrap();
    assert_ne!(snapshot(&model), snapshot(&held));
    assert_ne!(model.bodies_lattice(), held.bodies_lattice());
    assert_eq!(report.coupling.lubrication.active_pairs, 1);
    let gap = (report.bodies[1].midpoint_center[0] - report.bodies[0].midpoint_center[0]) * DX
        - 2.0 * RADIUS;
    let relative = (report.coupling.bodies[0].velocity_midpoint[0]
        - report.coupling.bodies[1].velocity_midpoint[0])
        * DX
        / DT;
    // SI asymptote independently uses the NEW viscosity and ACTUAL midpoint
    // gap/speed. Updating only the bulk or only a source receipt fails here.
    let resistance = 6.0
        * core::f64::consts::PI
        * 0.22
        * (RADIUS / 2.0).powi(2)
        * (1.0 / gap.max(film().roughness_gap_m) - 1.0 / film().cutoff_gap_m);
    close(
        model
            .impulse_ns(report.coupling.lubrication.body_impulses[0][0])
            .unwrap(),
        -resistance * relative * DT,
        2e-20,
    );
    close(
        model
            .energy_j(report.coupling.lubrication.dissipated_energy)
            .unwrap(),
        resistance * relative.powi(2) * DT,
        2e-23,
    );
    let after = momentum(&model);
    for axis in 0..3 {
        close(after[axis], before[axis], 3e-13);
    }
    assert_eq!(model.accepted_steps(), 3);
}

#[test]
fn rejected_refresh_and_same_source_retry_keep_the_original_trajectory() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let foreign_density = fluid(2000.0, 0.22);
    let unresolved_viscosity = fluid(1000.0, 1e-20);
    let carrier = body(3000.0);
    let gate = CancelGate::new();
    let stopped = CancelGate::new();
    stopped.request();
    let mut model = assembly(&first, &carrier);
    let mut reference = assembly(&first, &carrier);
    let loads = [[2e-6, 0.0, 0.0, 0.0, 0.0, 1e-10], [0.0; 6]];
    model.step_forces_si(loads, [0.0; 3], &gate).unwrap();
    reference.step_forces_si(loads, [0.0; 3], &gate).unwrap();
    let populations = snapshot(&model);
    let bodies = *model.bodies_lattice();
    let collision = model.grid().collision_model();
    assert!(matches!(
        model.rebind_fluid(&foreign_density, &gate),
        Err(MaterialImmersedError::Invalid(
            "live fluid density changes require conservative mass and momentum transfer"
        ))
    ));
    assert!(matches!(
        model.rebind_fluid(&second, &stopped),
        Err(MaterialImmersedError::Native(
            PartialSaturationError3::Cancelled
        ))
    ));
    assert!(matches!(
        model.rebind_fluid(&unresolved_viscosity, &gate),
        Err(MaterialImmersedError::Scaling(_))
    ));
    assert!(core::ptr::eq(model.scaling().material(), &first));
    assert_eq!(snapshot(&model), populations);
    assert_eq!(model.bodies_lattice(), &bodies);
    assert_eq!(model.grid().collision_model(), collision);
    assert_eq!(model.accepted_steps(), 1);
    model.rebind_fluid(&first, &gate).unwrap();
    assert_eq!(
        model.step_forces_si(loads, [0.0; 3], &gate).unwrap(),
        reference.step_forces_si(loads, [0.0; 3], &gate).unwrap()
    );
    assert_eq!(snapshot(&model), snapshot(&reference));
    assert_eq!(model.bodies_lattice(), reference.bodies_lattice());
    // A later valid intervention is independent of all refused trials.
    model.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert_eq!(
        model.step_forces_si(loads, [0.0; 3], &gate).unwrap(),
        reference.step_forces_si(loads, [0.0; 3], &gate).unwrap()
    );
}

#[test]
fn native_refresh_preserves_collision_family_higher_rates_and_populations() {
    use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, FaceBoundary3};
    let models = [
        (
            CollisionModel3::Bgk { tau: 0.8 },
            CollisionModel3::Bgk { tau: 1.25 },
        ),
        (
            CollisionModel3::CentralMoment {
                second_order_rate: 1.25,
                higher_order_rate: 0.9,
            },
            CollisionModel3::CentralMoment {
                second_order_rate: 0.8,
                higher_order_rate: 0.9,
            },
        ),
        (
            CollisionModel3::ReducedCumulant {
                second_order_rate: 1.25,
                third_order_rate: 0.9,
                fourth_order_rate: 1.1,
            },
            CollisionModel3::ReducedCumulant {
                second_order_rate: 0.8,
                third_order_rate: 0.9,
                fourth_order_rate: 1.1,
            },
        ),
    ];
    for (initial, updated) in models {
        let boundaries = BoundarySpec3::new([FaceBoundary3::Periodic; 6]);
        let mut grid =
            BoundaryGrid3::with_collision_model(4, 4, 4, initial, [1e-6, -2e-6, 3e-6], boundaries);
        let mut populations = grid.populations(1, 1, 1);
        populations[1] += 1e-4;
        populations[2] -= 1e-4;
        grid.set_populations(1, 1, 1, &populations);
        grid.rebind_shear_relaxation(1.25).unwrap();
        assert_eq!(grid.collision_model(), updated);
        assert_eq!(grid.populations(1, 1, 1), populations);
        assert_eq!(grid.boundaries(), boundaries);
        for tau in [0.5, 0.0, f64::INFINITY, f64::NAN] {
            assert!(grid.rebind_shear_relaxation(tau).is_err());
            assert_eq!(grid.collision_model(), updated);
            assert_eq!(grid.populations(1, 1, 1), populations);
        }
    }
}

#[test]
fn material_substitutions_drive_bulk_film_inertia_and_the_real_motion_together() {
    let mut endpoints = Vec::new();
    for (fluid_density, viscosity, carrier_density) in [
        (1000.0, 0.11, 3000.0),
        (1000.0, 0.22, 3000.0),
        (1000.0, 0.11, 6000.0),
        (2000.0, 0.11, 3000.0),
    ] {
        let source = fluid(fluid_density, viscosity);
        let carrier = body(carrier_density);
        let mut model = assembly(&source, &carrier);
        assert!(core::ptr::eq(model.scaling().material(), &source));
        assert!(core::ptr::eq(model.sources()[0].material(), &carrier));
        close(
            model.grid().viscosity(),
            viscosity / fluid_density * DT / (DX * DX),
            1e-16,
        );
        let mass = carrier_density * (4.0 * core::f64::consts::PI / 3.0) * RADIUS.powi(3);
        close(
            model.mass_kg(model.bodies_lattice()[0].mass).unwrap(),
            mass,
            1e-17,
        );
        close(
            model
                .inertia_kg_m2(model.bodies_lattice()[0].inertia.matrix().m[0])
                .unwrap(),
            0.4 * mass * RADIUS.powi(2),
            1e-22,
        );
        let before = momentum(&model);
        let report = model
            .step_forces_si([[0.0; 6]; 2], [0.0; 3], &CancelGate::new())
            .unwrap();
        assert_eq!(report.coupling.lubrication.active_pairs, 1);
        let gap = (report.bodies[1].midpoint_center[0] - report.bodies[0].midpoint_center[0]) * DX
            - 2.0 * RADIUS;
        let relative = (report.coupling.bodies[0].velocity_midpoint[0]
            - report.coupling.bodies[1].velocity_midpoint[0])
            * DX
            / DT;
        // Independent SI Stokes asymptote with the resolved-cutoff resistance
        // subtracted. This checks rho/dx/dt powers, force sign and actual gap.
        let resistance = 6.0
            * core::f64::consts::PI
            * viscosity
            * (RADIUS / 2.0).powi(2)
            * (1.0 / gap.max(film().roughness_gap_m) - 1.0 / film().cutoff_gap_m);
        let expected_impulse = -resistance * relative * DT;
        close(
            model
                .impulse_ns(report.coupling.lubrication.body_impulses[0][0])
                .unwrap(),
            expected_impulse,
            2e-20,
        );
        close(
            model
                .energy_j(report.coupling.lubrication.dissipated_energy)
                .unwrap(),
            resistance * relative * relative * DT,
            2e-23,
        );
        close(
            model
                .energy_j(report.coupling.lubrication.work_residual)
                .unwrap(),
            0.0,
            1e-23,
        );
        let after = momentum(&model);
        for axis in 0..3 {
            close(after[axis], before[axis], 3e-13);
        }
        for step in report.bodies {
            close(model.energy_j(step.work_residual).unwrap(), 0.0, 2e-22);
            for axis in 0..3 {
                close(
                    model.impulse_ns(step.momentum_residual[axis]).unwrap(),
                    0.0,
                    1e-19,
                );
            }
        }
        for _ in 0..2 {
            model
                .step_forces_si([[0.0; 6]; 2], [0.0; 3], &CancelGate::new())
                .unwrap();
        }
        assert_eq!(model.accepted_steps(), 3);
        endpoints.push(model.state_si(0).unwrap().velocity_m_s[0]);
    }
    for changed in &endpoints[1..] {
        assert!((changed - endpoints[0]).abs() > 1e-9);
    }
}

#[test]
fn held_physical_loads_enter_once_and_failed_steps_preserve_the_complete_trajectory() {
    let source = fluid(1000.0, 0.11);
    let carrier = body(3000.0);
    let mut model = assembly(&source, &carrier);
    let mut reference = assembly(&source, &carrier);
    let input = [[2e-6, 0.0, 0.0, 0.0, 0.0, 1e-10], [0.0; 6]];
    let before = momentum(&model);
    let accepted = model
        .step_forces_si(input, [0.001, -0.002, 0.003], &CancelGate::new())
        .unwrap();
    let after = momentum(&model);
    close(
        model.impulse_ns(after[0] - before[0]).unwrap(),
        input[0][0] * DT,
        1e-16,
    );
    assert_eq!(
        accepted,
        reference
            .step_forces_si(input, [0.001, -0.002, 0.003], &CancelGate::new())
            .unwrap()
    );
    let populations = snapshot(&model);
    let bodies = *model.bodies_lattice();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        model.step_forces_si(input, [0.0; 3], &gate),
        Err(MaterialImmersedError::Native(
            PartialSaturationError3::Cancelled
        ))
    ));
    // Finite force passes unit conversion, then the actual native candidate
    // leaves its low-Mach/geometry domain and must refuse atomically.
    assert!(
        model
            .step_forces_si([[1.0; 6]; 2], [0.0; 3], &CancelGate::new())
            .is_err()
    );
    assert_eq!(model.accepted_steps(), 1);
    assert_eq!(model.bodies_lattice(), &bodies);
    assert_eq!(snapshot(&model), populations);
    assert_eq!(
        model
            .step_forces_si(input, [0.0; 3], &CancelGate::new())
            .unwrap(),
        reference
            .step_forces_si(input, [0.0; 3], &CancelGate::new())
            .unwrap()
    );
    assert_eq!(snapshot(&model), snapshot(&reference));
    assert_eq!(model.bodies_lattice(), reference.bodies_lattice());
}

#[test]
fn geometry_and_complete_source_queries_cannot_be_silently_replaced() {
    let source = fluid(1000.0, 0.11);
    let carrier = body(3000.0);
    assert!(AdditionalMaterialSphere::new(&carrier, 2.0 * RADIUS).is_err());
    let bindings = [AdditionalMaterialSphere::new(&carrier, RADIUS).unwrap(); 2];
    let mut invalid = domain();
    invalid.maximum_cells = 511;
    assert!(
        RetainedFluidSpheres::new(
            &source,
            invalid,
            bindings,
            initial(),
            film(),
            &CancelGate::new()
        )
        .is_err()
    );
    invalid = domain();
    invalid.cells[0] = 7;
    assert!(
        RetainedFluidSpheres::new(
            &source,
            invalid,
            bindings,
            initial(),
            film(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut invalid_film = film();
    invalid_film.roughness_gap_m = 1.0;
    assert!(
        RetainedFluidSpheres::new(
            &source,
            domain(),
            bindings,
            initial(),
            invalid_film,
            &CancelGate::new()
        )
        .is_err()
    );
    let source_card = card(
        "unsupported fluid",
        &[
            ("density", Density::DIMS, 1000.0),
            ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, 0.11),
            (
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                1000.0,
            ),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        ],
        "liquid",
    );
    assert!(
        resolve_liquid_state(
            &source_card,
            &QueryPoint::new().with("T", 400.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}
