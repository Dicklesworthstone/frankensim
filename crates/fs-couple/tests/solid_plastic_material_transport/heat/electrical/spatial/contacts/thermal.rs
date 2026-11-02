//! Same native traction, two source properties, one finite heat/current solve.
//! Nonzero interface heat transfer uses original reference-area Fourier paths;
//! electrical resistance uses the actual prescribed X+u interface metric.
use super::*;

const R0: f64 = 0.001;
const T_SLOPE: f64 = 0.01;
const P_SLOPE: f64 = 1e-6;
const COMPRESSION: f64 = 1e-4;
const THREE_BULK: f64 = 2e9;
const PRESSURE: [SpatialElectricalContactPressure<'static>; 1] =
    [SpatialElectricalContactPressure {
        contact: 0,
        pressure_axis: "p",
        agreement_tolerance_pa: 1.0,
    }];
const THERMAL: [SpatialElectricalThermalContact; 1] = [SpatialElectricalThermalContact {
    contact: 0,
    thermal_link: CONTACT_LINK,
}];

fn fixture(thermal_r0: f64, area: f64) -> ContactFixture {
    let material = fixture_with_property(
        EMISSIVITY,
        Some(scalar(
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            BULK_RHO,
        )),
    );
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let point = QueryPoint::new()
        .with_quantity("Ta", temperature, 300.0)
        .unwrap()
        .with_quantity("Tb", temperature, 300.0)
        .unwrap()
        .with_quantity("p", pressure_quantity(), 0.0)
        .unwrap();
    let original = interface(&material.cards[1], R0);
    // This declared A-side film law is not a guessed mean interface state.
    // Unequal Joule heat splits produce a real thermal flux, and changing only
    // thermal R'' then changes Ta and therefore the electrical voltage.
    let mut electrical = scalar(ER, ER_DIMS, R0);
    electrical.value = PropertyValue::Surface2D {
        axes: ["Ta".into(), "p".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: [
            [280.0, 0.0],
            [320.0, 0.0],
            [280.0, 600000.0],
            [320.0, 600000.0],
        ]
        .map(|q| (q, R0 * (1.0 + T_SLOPE * (q[0] - 300.0) - P_SLOPE * q[1])))
        .to_vec(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: ER_DIMS,
    };
    electrical.interpolation = InterpolationPolicy::LinearInside;
    electrical.validity = ValidityDomain::unconstrained()
        .with_quantity("Ta", temperature, 280.0, 320.0)
        .with_quantity("p", pressure_quantity(), 0.0, 600000.0);
    let mut heat = scalar(CONTACT_R, CONTACT_DIMS, thermal_r0);
    heat.value = PropertyValue::Curve {
        abscissa: "p".into(),
        abscissa_dims: Pressure::DIMS,
        knots: vec![(0.0, thermal_r0), (600000.0, 0.4 * thermal_r0)],
        dims: CONTACT_DIMS,
    };
    heat.interpolation = InterpolationPolicy::LinearInside;
    heat.validity =
        ValidityDomain::unconstrained().with_quantity("p", pressure_quantity(), 0.0, 600000.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(electrical).unwrap();
    claims.insert_claim(heat).unwrap();
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
        fraction: FRACTION,
        mean_temperature_axis: None,
        expansion: false,
        thermal_contact_area_m2: Some(area),
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

fn conductance(thermal_r0: f64, pressure: f64) -> f64 {
    let bulk = AREA / ((L / 4.0) / 0.01 + (L / 4.0) / 0.01);
    1.0 / (1.0 / bulk + thermal_r0 * (1.0 - P_SLOPE * pressure) / AREA)
}

fn endpoint(old: [f64; 6], compression: f64, thermal_r0: f64, dt: f64) -> [f64; 6] {
    // Independent 2x2 backward-Euler balances of the two finite contact cells.
    // E is constant and alpha=0, so unequal temperatures do not fabricate a
    // traction mismatch or reversible heat. Pressure is 3K*compression.
    let pressure = THREE_BULK * compression;
    let stretch = 1.0 - compression;
    let capacity = MASS * 5.0;
    let bulk = dt * CURRENT * CURRENT * BULK_RHO * L / (3.0 * AREA * stretch);
    let contact = dt * CURRENT * CURRENT * R0 / (AREA * stretch * stretch);
    let g = dt * conductance(thermal_r0, pressure);
    let base = 1.0 - P_SLOPE * pressure + T_SLOPE * (old[2] - 300.0);
    let a = capacity + g - FRACTION * contact * T_SLOPE;
    let b = -g;
    let c = -g - (1.0 - FRACTION) * contact * T_SLOPE;
    let d = capacity + g;
    let rhs_a = bulk + FRACTION * contact * base + g * (old[3] - old[2]);
    let rhs_b = bulk + (1.0 - FRACTION) * contact * base - g * (old[3] - old[2]);
    let mut next = old.map(|t| t + bulk / capacity);
    next[2] += (rhs_a * d - b * rhs_b) / (a * d - b * c) - bulk / capacity;
    next[3] += (a * rhs_b - c * rhs_a) / (a * d - b * c) - bulk / capacity;
    next
}

fn check(
    thermal: &DepositionSimulation<'_>,
    solid: &ThermoplasticSolid<'_>,
    problem: &SpatialConduction,
    report: &ElectricalSolidReport<HeatNetworkSample, SpatialConductionThermalStep>,
    binding: &ConductionContactThermal<'_>,
    old: [f64; 6],
    compression: f64,
    thermal_r0: f64,
    current: f64,
    dt: f64,
) -> [f64; 6] {
    let actual = IDS.map(|cell| thermal.temperature_k(cell).unwrap());
    let expected = endpoint(old, compression, thermal_r0, dt);
    for i in 0..6 {
        close(actual[i], expected[i], 3e-6);
    }
    assert!(actual[3] > actual[2] + 1.0);
    let pressure = -solid.solution().element_stress[2][0];
    close(pressure, -solid.solution().element_stress[3][0], 1e-6);
    close(pressure, THREE_BULK * compression, 1e-6);
    let electrical = &problem.contacts()[0].source;
    assert_eq!(
        electrical.resolved().state_coordinate("p").unwrap(),
        pressure
    );
    let source = &report.circuit_thermal.source.material;
    assert_eq!(source.thermal_contacts.len(), 1);
    let contact = &source.thermal_contacts[0];
    assert_eq!(contact.cells, [IDS[2], IDS[3]]);
    assert_eq!(
        contact.interface.card_identity(),
        binding.point.interface.content_hash()
    );
    let receipt = contact.interface.receipt();
    binding
        .point
        .interface
        .claims()
        .verify_receipt(receipt)
        .unwrap();
    let coordinate = |axis: &str| {
        receipt
            .query_point
            .iter()
            .find(|(name, _)| name == axis)
            .unwrap()
            .1
    };
    assert_eq!(coordinate("p"), pressure);
    for (side, axis) in ["Ta", "Tb"].into_iter().enumerate() {
        close(coordinate(axis), actual[side + 2], 1e-10);
        close(
            contact.phases[side].temperature_k(),
            actual[side + 2],
            1e-10,
        );
    }
    close(
        contact.interface.value_m2_k_per_w(),
        thermal_r0 * (1.0 - P_SLOPE * pressure),
        1e-12,
    );
    close(contact.area_m2, AREA, 0.0);
    close(thermal.links()[&CONTACT_LINK].area_m2, AREA, 0.0);
    close(
        contact.bulk_conductance_w_k,
        AREA / (L / (2.0 * 0.01)),
        1e-14,
    );
    close(
        contact.conductance_w_k,
        conductance(thermal_r0, pressure),
        1e-14,
    );
    let stretch = 1.0 - compression;
    let bulk_r = BULK_RHO * L / (AREA * stretch);
    let specific_r = R0 * (1.0 - P_SLOPE * pressure + T_SLOPE * (actual[2] - 300.0));
    let contact_r = specific_r / (AREA * stretch * stretch);
    close(electrical.resistance_ohm_m2(), specific_r, 1e-12);
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
    let exchange = dt * contact.conductance_w_k * (actual[3] - actual[2]);
    assert!(exchange > 1e-6);
    let bulk_heat = dt * current * current * bulk_r / 3.0;
    let contact_heat = dt * current * current * contact_r;
    close(
        MASS * 5.0 * (actual[2] - old[2]),
        bulk_heat + FRACTION * contact_heat + exchange,
        1e-9,
    );
    close(
        MASS * 5.0 * (actual[3] - old[3]),
        bulk_heat + (1.0 - FRACTION) * contact_heat - exchange,
        1e-9,
    );
    let joule = 6.0 * bulk_heat + contact_heat;
    close(step.supplied_work_j, joule, 1e-9);
    close(report.circuit_thermal.deposited_joule_j, joule, 1e-9);
    close(
        report.circuit_thermal.thermal.stored_energy_change_j,
        joule,
        1e-9,
    );
    close(report.deposited_mechanical_heat_j, 0.0, 1e-12);
    close(report.combined_residual_j, 0.0, 1e-8);
    actual
}

#[test]
fn one_solved_pressure_drives_heat_and_current_with_material_swap_and_inertial_continuation() {
    let mut observations = Vec::new();
    for thermal_r0 in [0.5, 2.0] {
        with_contact_fixture(
            fixture(thermal_r0, AREA),
            |cx, thermal, solid, problem, prepared, bindings, electrodes| {
                let first = solid
                    .advance_spatial_electrothermal_contact_candidate(
                        thermal,
                        problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        &THERMAL,
                        1.0,
                        &imposed(COMPRESSION, 0.0),
                        &[],
                        IncrementSettings::default(),
                        contact_controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                let old = check(
                    first.thermal(),
                    first.solid(),
                    first.problem(),
                    first.report(),
                    &bindings[0],
                    [300.0; 6],
                    COMPRESSION,
                    thermal_r0,
                    CURRENT,
                    1.0,
                );
                observations.push((
                    old[2],
                    first.report().circuit_thermal.electrical.field.voltage_v[9],
                ));
                let pin = first
                    .report()
                    .circuit_thermal
                    .source
                    .material
                    .thermal_contacts[0]
                    .interface
                    .receipt()
                    .selected;
                let (thermal, solid, problem, _) = first.into_parts();
                let second = solid
                    .advance_spatial_electrothermal_contact_candidate(
                        &thermal,
                        &problem,
                        prepared,
                        current_drive(electrodes, -CURRENT),
                        bindings,
                        &PRESSURE,
                        &THERMAL,
                        1.25,
                        &imposed(2.0 * COMPRESSION, 0.0),
                        &[],
                        IncrementSettings::default(),
                        contact_controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                let old = check(
                    second.thermal(),
                    second.solid(),
                    second.problem(),
                    second.report(),
                    &bindings[0],
                    old,
                    2.0 * COMPRESSION,
                    thermal_r0,
                    -CURRENT,
                    0.25,
                );
                assert_eq!(
                    second
                        .report()
                        .circuit_thermal
                        .source
                        .material
                        .thermal_contacts[0]
                        .interface
                        .receipt()
                        .selected,
                    pin
                );
                let (thermal, solid, problem, _) = second.into_parts();
                let velocity = 0.1;
                let dt = 0.1;
                let moving =
                    ThermoplasticDynamics::try_new(solid, vec![[velocity, 0.0, 0.0]; 12], cx)
                        .unwrap();
                let next = moving
                    .advance_spatial_electrothermal_contact_candidate(
                        &thermal,
                        &problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        &PRESSURE,
                        &THERMAL,
                        1.25 + dt,
                        &imposed(2.0 * COMPRESSION, velocity * dt),
                        &[],
                        IncrementSettings::default(),
                        contact_controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                check(
                    next.thermal(),
                    next.dynamics().solid(),
                    next.problem(),
                    &next.report().material,
                    &bindings[0],
                    old,
                    2.0 * COMPRESSION,
                    thermal_r0,
                    CURRENT,
                    dt,
                );
                close(next.report().inertia.kinetic_change_j, 0.0, 1e-12);
                close(next.report().inertia.numerical_loss_j, 0.0, 1e-12);
                assert_eq!(moving.time_s(), 1.25);
            },
        );
    }
    // Only the thermal resistance DATA changed, with electrical law unchanged.
    // Larger thermal resistance reduces transfer to A, its temperature and V.
    assert!(observations[0].0 > observations[1].0 + 0.05);
    assert!(observations[0].1 > observations[1].1 + 1e-4);
}

#[test]
fn foreign_joint_contact_owners_and_reference_metric_refuse_atomically() {
    for area in [AREA, 2.0 * AREA] {
        with_contact_fixture(
            fixture(0.5, area),
            |cx, thermal, solid, problem, prepared, bindings, electrodes| {
                let cells = thermal.cells().clone();
                let links = thermal.links().clone();
                let solution = solid.solution().clone();
                let contacts = problem.contacts().to_vec();
                let geometry = problem.positions_m().to_vec();
                let foreign_card = fixture(2.0, area).interface;
                let mut foreign = bindings.to_vec();
                foreign[0].point.interface = &foreign_card;
                for (bindings, mappings, pressures, message) in [
                    (
                        bindings,
                        THERMAL.to_vec(),
                        PRESSURE.to_vec(),
                        "reference face area",
                    ),
                    (
                        bindings,
                        vec![SpatialElectricalThermalContact {
                            contact: 0,
                            thermal_link: 999,
                        }],
                        PRESSURE.to_vec(),
                        "original pinned thermal source binding",
                    ),
                    (
                        foreign.as_slice(),
                        THERMAL.to_vec(),
                        PRESSURE.to_vec(),
                        "same ordered interface",
                    ),
                    (
                        bindings,
                        THERMAL.to_vec(),
                        vec![],
                        "original solved-pressure and thermal owners",
                    ),
                ] {
                    if area == AREA && message == "reference face area" {
                        continue;
                    }
                    let error = solid
                        .advance_spatial_electrothermal_contact_candidate(
                            thermal,
                            problem,
                            prepared,
                            current_drive(electrodes, CURRENT),
                            bindings,
                            &pressures,
                            &mappings,
                            1.0,
                            &imposed(COMPRESSION, 0.0),
                            &[],
                            IncrementSettings::default(),
                            contact_controls(),
                            &[0.0; 6],
                            cx,
                        )
                        .err()
                        .expect("foreign pressure/heat correspondence must refuse");
                    assert!(format!("{error:?}").contains(message), "{error:?}");
                }
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(thermal.links(), &links);
                assert_eq!(solid.solution(), &solution);
                assert_eq!(problem.contacts(), &contacts);
                assert_eq!(problem.positions_m(), &geometry);
                assert_eq!(solid.time_s(), 0.0);
            },
        );
    }
}
