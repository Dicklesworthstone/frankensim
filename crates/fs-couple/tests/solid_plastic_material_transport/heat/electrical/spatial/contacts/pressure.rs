//! Solved compression closes material pressure/resistance/heat/stress feedback.
//! Symmetric manufactured caloric receivers keep both actual tractions equal;
//! the surface law explicitly owns mean temperature and static pressure.
use super::*;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
};

const ALPHA: f64 = 2e-5;
const COMPRESSION: f64 = 1e-4;
const THREE_BULK: f64 = 2e9;
const NINE_BULK: f64 = 6e9;
const R0: f64 = 0.001;
const T_SLOPE: f64 = 0.01;
const PRESSURE: [SpatialElectricalContactPressure<'static>; 1] =
    [SpatialElectricalContactPressure {
        contact: 0,
        pressure_axis: "p",
        agreement_tolerance_pa: 1.0,
    }];

fn controls() -> SpatialElectricalSolidSolve {
    let mut result = contact_controls();
    // A 1e-9 J derivative probe changes local traction by <0.05 Pa. The
    // explicit 1 Pa physical agreement gate remains tight during these probes.
    result.coupled.coupled.difference_step_j = 1e-9;
    result
}

fn fixture(beta: f64, accepted_pressure: f64, wrong_kind: bool) -> ContactFixture {
    let material = fixture_with_properties(
        EMISSIVITY,
        &[
            scalar(
                "electrical_resistivity",
                ELECTRICAL_RESISTIVITY_DIMS,
                BULK_RHO,
            ),
            scalar(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                ALPHA,
            ),
        ],
    );
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pressure_quantity = if wrong_kind {
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static))
    } else {
        pressure_quantity()
    };
    let point = QueryPoint::new()
        .with_quantity("Ta", temperature, 300.0)
        .unwrap()
        .with_quantity("Tb", temperature, 300.0)
        .unwrap()
        .with_quantity("Tm", temperature, 300.0)
        .unwrap()
        .with_quantity("p", pressure_quantity, accepted_pressure)
        .unwrap();
    let original = interface(&material.cards[1], R0);
    let mut claim = scalar(ER, ER_DIMS, R0);
    if !wrong_kind {
        claim.value = PropertyValue::Surface2D {
            axes: ["Tm".into(), "p".into()],
            axis_dims: [Temperature::DIMS, fs_qty::Pressure::DIMS],
            samples: [
                [280.0, 0.0],
                [320.0, 0.0],
                [280.0, 600000.0],
                [320.0, 600000.0],
            ]
            .map(|p| (p, R0 * (1.0 + T_SLOPE * (p[0] - 300.0) - beta * p[1])))
            .to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: ER_DIMS,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("Tm", temperature, 280.0, 320.0)
            .with_quantity("p", pressure_quantity, 0.0, 600000.0);
    }
    // For the role-refusal case a constant law admits the original field, but
    // the new mechanical owner still must not interpret a Stress axis as Pressure.
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    ContactFixture {
        material,
        interface,
        point,
        fraction: 0.5,
        mean_temperature_axis: Some("Tm"),
        expansion: true,
        thermal_contact_area_m2: None,
    }
}

fn imposed(compression: f64, translation: f64) -> Vec<(usize, f64)> {
    specimen_mesh()
        .0
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |axis| {
                (
                    3 * node + axis,
                    -compression * x[axis] + if axis == 0 { translation } else { 0.0 },
                )
            })
        })
        .collect()
}

fn endpoint(
    old: [f64; 6],
    old_compression: f64,
    compression: f64,
    beta: f64,
    current: f64,
    dt: f64,
) -> [f64; 6] {
    let stretch = 1.0 - compression;
    let bulk_heat = dt * current * current * BULK_RHO * L / (3.0 * AREA * stretch);
    let contact_heat = dt * current * current * R0 / (AREA * stretch * stretch);
    let slope = T_SLOPE - beta * THREE_BULK * ALPHA;
    core::array::from_fn(|cell| {
        // Original reversible heat: -Told*alpha*delta(tr(sigma)). Both
        // reference mass and stress-free caloric Cp remain unchanged.
        let capacity = MASS * 5.0 - VOLUME * NINE_BULK * old[cell] * ALPHA * ALPHA;
        let compression_heat =
            VOLUME * NINE_BULK * old[cell] * ALPHA * (compression - old_compression);
        let surface = if cell == 2 || cell == 3 {
            0.5 * contact_heat
        } else {
            0.0
        };
        let base = 1.0 - beta * THREE_BULK * compression + slope * (old[cell] - 300.0);
        old[cell] + (bulk_heat + compression_heat + surface * base) / (capacity - surface * slope)
    })
}

fn check(
    thermal: &DepositionSimulation<'_>,
    solid: &ThermoplasticSolid<'_>,
    problem: &SpatialConduction,
    report: &ElectricalSolidReport<HeatNetworkSample, SpatialConductionThermalStep>,
    old: [f64; 6],
    old_compression: f64,
    compression: f64,
    beta: f64,
    current: f64,
    dt: f64,
) -> [f64; 6] {
    let actual = IDS.map(|cell| thermal.temperature_k(cell).unwrap());
    let expected = endpoint(old, old_compression, compression, beta, current, dt);
    for i in 0..6 {
        close(actual[i], expected[i], 3e-6);
    }
    close(actual[2], actual[3], 1e-6);
    let pressure = -solid.solution().element_stress[2][0];
    close(pressure, -solid.solution().element_stress[3][0], 1.0);
    close(
        pressure,
        THREE_BULK * (compression + ALPHA * (actual[2] - 300.0)),
        0.2,
    );
    assert!(pressure > 1e5);
    let source = &problem.contacts()[0].source;
    // Source pressure must be the actual first traction, not authored preload,
    // a temperature-only estimate, or a second independently evaluated model.
    assert_eq!(source.resolved().state_coordinate("p").unwrap(), pressure);
    close(
        source.resolved().state_coordinate("Ta").unwrap(),
        actual[2],
        1e-10,
    );
    close(
        source.resolved().state_coordinate("Tb").unwrap(),
        actual[3],
        1e-10,
    );
    let mean = 0.5 * (actual[2] + actual[3]);
    close(
        source.resolved().state_coordinate("Tm").unwrap(),
        mean,
        1e-10,
    );
    let specific_r = R0 * (1.0 + T_SLOPE * (mean - 300.0) - beta * pressure);
    close(source.resistance_ohm_m2(), specific_r, 1e-12);
    let stretch = 1.0 - compression;
    let bulk_r = BULK_RHO * L / (AREA * stretch);
    let contact_r = specific_r / (AREA * stretch * stretch);
    let step = &report.circuit_thermal.electrical;
    close(
        step.field.voltage_v[9],
        current * (2.0 * bulk_r + contact_r),
        1e-9,
    );
    close(
        step.contacts[0].response.power_w,
        current * current * contact_r,
        1e-10,
    );
    let joule = dt * current * current * (2.0 * bulk_r + contact_r);
    let reversible: f64 = (0..6)
        .map(|i| {
            VOLUME
                * NINE_BULK
                * old[i]
                * ALPHA
                * (compression - old_compression + ALPHA * (actual[i] - old[i]))
        })
        .sum();
    close(step.supplied_work_j, joule, 1e-9);
    close(report.circuit_thermal.deposited_joule_j, joule, 1e-9);
    close(report.deposited_mechanical_heat_j, reversible, 1e-8);
    close(report.physical_plastic_heat_j, 0.0, 1e-12);
    close(
        report.circuit_thermal.thermal.stored_energy_change_j,
        joule + reversible,
        1e-8,
    );
    close(report.combined_residual_j, 0.0, 1e-8);
    actual
}

#[test]
fn native_thermal_compression_changes_interface_resistance_and_continues_with_inertia() {
    let mut observations = Vec::new();
    for beta in [0.0, 1e-6] {
        with_contact_fixture(
            fixture(beta, 0.0, false),
            |cx, thermal, solid, problem, prepared, bindings, electrodes| {
                let first = solid
                    .advance_spatial_electrical_pressure_contact_candidate(
                        thermal,
                        problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        1.0,
                        &imposed(COMPRESSION, 0.0),
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                let temperatures = check(
                    first.thermal(),
                    first.solid(),
                    first.problem(),
                    first.report(),
                    [300.0; 6],
                    0.0,
                    COMPRESSION,
                    beta,
                    CURRENT,
                    1.0,
                );
                observations.push(temperatures[2]);
                let (thermal, solid, problem, _) = first.into_parts();
                let second = solid
                    .advance_spatial_electrical_pressure_contact_candidate(
                        &thermal,
                        &problem,
                        prepared,
                        current_drive(electrodes, -CURRENT),
                        bindings,
                        &PRESSURE,
                        1.25,
                        &imposed(COMPRESSION, 0.0),
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                let temperatures = check(
                    second.thermal(),
                    second.solid(),
                    second.problem(),
                    second.report(),
                    temperatures,
                    COMPRESSION,
                    COMPRESSION,
                    beta,
                    -CURRENT,
                    0.25,
                );
                let (thermal, solid, problem, _) = second.into_parts();
                let velocity = 0.1;
                let dt = 0.1;
                let moving =
                    ThermoplasticDynamics::try_new(solid, vec![[velocity, 0.0, 0.0]; 12], cx)
                        .unwrap();
                let next = moving
                    .advance_spatial_electrical_pressure_contact_candidate(
                        &thermal,
                        &problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        1.25 + dt,
                        &imposed(COMPRESSION, velocity * dt),
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                check(
                    next.thermal(),
                    next.dynamics().solid(),
                    next.problem(),
                    &next.report().material,
                    temperatures,
                    COMPRESSION,
                    COMPRESSION,
                    beta,
                    CURRENT,
                    dt,
                );
                for v in next.dynamics().velocity_mps() {
                    close(v[0], velocity, 1e-12);
                }
                close(next.report().inertia.kinetic_change_j, 0.0, 1e-12);
                close(next.report().inertia.numerical_loss_j, 0.0, 1e-12);
                assert_eq!(moving.time_s(), 1.25);
            },
        );
    }
    // Only pressure sensitivity in the interface DATA changes. Compression
    // lowers R'', reducing heat and therefore the actual constrained stress.
    assert!(observations[0] > observations[1] + 1.0);
}

#[test]
fn mismatched_tensile_unprescribed_and_stale_pressure_are_atomic_refusals() {
    with_contact_fixture(
        fixture(1e-6, 0.0, false),
        |cx, thermal, solid, problem, prepared, bindings, electrodes| {
            let cells = thermal.cells().clone();
            let solution = solid.solution().clone();
            let sources = problem.contacts().to_vec();
            let geometry = problem.positions_m().to_vec();
            let mut mismatch = imposed(COMPRESSION, 0.0);
            // Extra axial compression of B leaves its interface coordinates fixed,
            // but produces an independently different native normal traction.
            for (dof, u) in &mut mismatch {
                if *dof >= 27 && *dof % 3 == 0 {
                    *u -= COMPRESSION * L;
                }
            }
            let mut incomplete = imposed(COMPRESSION, 0.0);
            incomplete.retain(|(dof, _)| *dof != 9);
            for (prescribed, text) in [
                (mismatch, "two-sided native normal tractions do not agree"),
                (imposed(-COMPRESSION, 0.0), "tensile normal traction"),
                (incomplete, "every paired-face displacement prescribed"),
            ] {
                let error = solid
                    .advance_spatial_electrical_pressure_contact_candidate(
                        thermal,
                        problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        1.0,
                        &prescribed,
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .err()
                    .expect("pressure must refuse");
                assert!(format!("{error:?}").contains(text));
            }
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(solid.solution(), &solution);
            assert_eq!(problem.contacts(), &sources);
            assert_eq!(problem.positions_m(), &geometry);
            assert_eq!(solid.time_s(), 0.0);
        },
    );
    for (accepted, wrong_kind, message) in [
        (10.0, false, "accepted electrical contacts must match"),
        (0.0, true, "original static pressure coordinate"),
    ] {
        with_contact_fixture(
            fixture(1e-6, accepted, wrong_kind),
            |cx, thermal, solid, problem, prepared, bindings, electrodes| {
                let cells = thermal.cells().clone();
                let solution = solid.solution().clone();
                let contacts = problem.contacts().to_vec();
                let error = solid
                    .advance_spatial_electrical_pressure_contact_candidate(
                        thermal,
                        problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        1.0,
                        &imposed(COMPRESSION, 0.0),
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .err()
                    .expect("foreign pressure source must refuse");
                assert!(format!("{error:?}").contains(message));
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(solid.solution(), &solution);
                assert_eq!(problem.contacts(), &contacts);
            },
        );
    }
}
