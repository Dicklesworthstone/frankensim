//! Unequal native normal/transverse stretch changes BOTH thermal half-paths
//! and area. This is a manufactured scalar planar-path test, not a material fit.
use super::thermal::{P_SLOPE, PRESSURE, R0, T_SLOPE, THERMAL, fixture};
use super::*;
use fs_couple::solid_enthalpy::plastic::electrical::spatial::SpatialThermalGeometrySolve;

const HEAT_R0: f64 = 0.5;
const FIRST: [f64; 3] = [-0.003, 0.006, 0.0054];
const SECOND: [f64; 3] = [-0.004, 0.008, 0.0072];
const GEOMETRY: SpatialThermalGeometrySolve = SpatialThermalGeometrySolve {
    maximum_displacement_gradient: 0.03,
    orthogonality_tolerance: 1e-10,
};

fn imposed(strain: [f64; 3], translation: f64, shear: f64) -> Vec<(usize, f64)> {
    specimen_mesh()
        .0
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |axis| {
                (
                    3 * node + axis,
                    strain[axis] * x[axis]
                        + if axis == 0 {
                            translation + shear * x[1]
                        } else {
                            0.0
                        },
                )
            })
        })
        .collect()
}

// Independent small-strain Hooke stress for E=1.2 GPa, nu=.2; alpha=0.
fn pressure(strain: [f64; 3]) -> f64 {
    -(1e9 * strain[0] + (1e9 / 3.0) * strain.into_iter().sum::<f64>())
}

fn geometry(strain: [f64; 3]) -> (f64, f64, f64) {
    let area = AREA * (1.0 + strain[1]) * (1.0 + strain[2]);
    let distance = (L / 4.0) * (1.0 + strain[0]);
    let bulk = area / (2.0 * distance / 0.01);
    let heat = 1.0 / (1.0 / bulk + HEAT_R0 * (1.0 - P_SLOPE * pressure(strain)) / area);
    (area, bulk, heat)
}

fn endpoint(old: [f64; 6], strain: [f64; 3], dt: f64) -> [f64; 6] {
    let (area, _, conductance) = geometry(strain);
    let capacity = MASS * 5.0;
    let bulk = dt * CURRENT * CURRENT * BULK_RHO * L * (1.0 + strain[0]) / (3.0 * area);
    let film = dt * CURRENT * CURRENT * R0 / area;
    let g = dt * conductance;
    let base = 1.0 - P_SLOPE * pressure(strain) + T_SLOPE * (old[2] - 300.0);
    // The original two finite caloric balances form this exact 2x2 BE system.
    let a = capacity + g - FRACTION * film * T_SLOPE;
    let b = -g;
    let c = -g - (1.0 - FRACTION) * film * T_SLOPE;
    let d = capacity + g;
    let rhs_a = bulk + FRACTION * film * base + g * (old[3] - old[2]);
    let rhs_b = bulk + (1.0 - FRACTION) * film * base - g * (old[3] - old[2]);
    let mut next = old.map(|t| t + bulk / capacity);
    next[2] = old[2] + (rhs_a * d - b * rhs_b) / (a * d - b * c);
    next[3] = old[3] + (a * rhs_b - c * rhs_a) / (a * d - b * c);
    next
}

#[allow(clippy::too_many_arguments)]
fn check(
    thermal: &DepositionSimulation<'_>,
    solid: &ThermoplasticSolid<'_>,
    field: &SpatialConduction,
    report: &ElectricalSolidReport<HeatNetworkSample, SpatialConductionThermalStep>,
    binding: &ConductionContactThermal<'_>,
    old: [f64; 6],
    strain: [f64; 3],
    dt: f64,
) -> [f64; 6] {
    let actual = IDS.map(|cell| thermal.temperature_k(cell).unwrap());
    let expected = endpoint(old, strain, dt);
    for i in 0..6 {
        close(actual[i], expected[i], 3e-6);
        close(thermal.cells()[&IDS[i]].mass_kg, MASS, 1e-18);
        close(thermal.cells()[&IDS[i]].geometry.volume_m3, VOLUME, 1e-20);
    }
    let (area, bulk_conductance, conductance) = geometry(strain);
    let contact = &report.circuit_thermal.source.material.thermal_contacts[0];
    assert_eq!(contact.cells, [IDS[2], IDS[3]]);
    assert_eq!(
        contact.interface.card_identity(),
        binding.point.interface.content_hash()
    );
    binding
        .point
        .interface
        .claims()
        .verify_receipt(contact.interface.receipt())
        .unwrap();
    close(contact.area_m2, area, 1e-16);
    close(contact.bulk_conductance_w_k, bulk_conductance, 1e-14);
    close(contact.conductance_w_k, conductance, 1e-14);
    // Reference geometry is inventory, not an already stretched next-step input.
    close(thermal.links()[&CONTACT_LINK].area_m2, AREA, 0.0);
    let p = -solid.solution().element_stress[2][0];
    close(p, pressure(strain), 2e-6);
    assert_eq!(
        field.contacts()[0].source.resolved().state_coordinate("p"),
        Some(p)
    );
    for (axis, expected) in [("p", p), ("Ta", actual[2]), ("Tb", actual[3])] {
        let q = contact
            .interface
            .receipt()
            .query_point
            .iter()
            .find(|(name, _)| name == axis)
            .unwrap()
            .1;
        close(q, expected, 1e-10);
    }
    let bulk_r = BULK_RHO * L * (1.0 + strain[0]) / area;
    let film_r = R0 * (1.0 - P_SLOPE * p + T_SLOPE * (actual[2] - 300.0)) / area;
    close(
        contact.interface.value_m2_k_per_w(),
        HEAT_R0 * (1.0 - P_SLOPE * p),
        1e-12,
    );
    let step = &report.circuit_thermal.electrical;
    close(
        step.field.voltage_v[9],
        CURRENT * (2.0 * bulk_r + film_r),
        1e-9,
    );
    let film_heat = dt * CURRENT * CURRENT * film_r;
    let bulk_heat = dt * CURRENT * CURRENT * bulk_r / 3.0;
    let exchange = dt * conductance * (actual[3] - actual[2]);
    assert!(exchange > 1e-6);
    close(
        MASS * 5.0 * (actual[2] - old[2]),
        bulk_heat + FRACTION * film_heat + exchange,
        1e-9,
    );
    close(
        MASS * 5.0 * (actual[3] - old[3]),
        bulk_heat + (1.0 - FRACTION) * film_heat - exchange,
        1e-9,
    );
    let work = 6.0 * bulk_heat + film_heat;
    close(step.supplied_work_j, work, 1e-9);
    close(report.circuit_thermal.deposited_joule_j, work, 1e-9);
    close(
        report.circuit_thermal.thermal.stored_energy_change_j,
        work,
        1e-9,
    );
    close(report.deposited_mechanical_heat_j, 0.0, 1e-12);
    close(report.combined_residual_j, 0.0, 1e-8);
    actual
}

#[test]
fn actual_area_and_normal_half_paths_close_heat_current_and_inertial_continuation() {
    with_contact_fixture(
        fixture(HEAT_R0, AREA),
        |cx, thermal, solid, field, prepared, bindings, electrodes| {
            let first = solid
                .advance_spatial_electrothermal_geometry_candidate(
                    thermal,
                    field,
                    prepared,
                    current_drive(electrodes, CURRENT),
                    bindings,
                    &PRESSURE,
                    &THERMAL,
                    GEOMETRY,
                    1.0,
                    &imposed(FIRST, 0.0, 0.0),
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
                FIRST,
                1.0,
            );
            let pin = first
                .report()
                .circuit_thermal
                .source
                .material
                .thermal_contacts[0]
                .interface
                .receipt()
                .selected;
            let (thermal, solid, field, _) = first.into_parts();
            let second = solid
                .advance_spatial_electrothermal_geometry_candidate(
                    &thermal,
                    &field,
                    prepared,
                    current_drive(electrodes, CURRENT),
                    bindings,
                    &PRESSURE,
                    &THERMAL,
                    GEOMETRY,
                    1.25,
                    &imposed(SECOND, 0.0, 0.0),
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
                SECOND,
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
            let (thermal, solid, field, _) = second.into_parts();
            let velocity = 0.01;
            let dt = 0.1;
            let moving =
                ThermoplasticDynamics::try_new(solid, vec![[velocity, 0.0, 0.0]; 12], cx).unwrap();
            let third = moving
                .advance_spatial_electrothermal_geometry_candidate(
                    &thermal,
                    &field,
                    prepared,
                    current_drive(electrodes, CURRENT),
                    bindings,
                    &PRESSURE,
                    &THERMAL,
                    GEOMETRY,
                    1.25 + dt,
                    &imposed(SECOND, velocity * dt, 0.0),
                    &[],
                    IncrementSettings::default(),
                    contact_controls(),
                    &[0.0; 6],
                    cx,
                )
                .unwrap();
            check(
                third.thermal(),
                third.dynamics().solid(),
                third.problem(),
                &third.report().material,
                &bindings[0],
                old,
                SECOND,
                dt,
            );
            close(third.report().inertia.kinetic_change_j, 0.0, 1e-12);
            close(third.report().inertia.numerical_loss_j, 0.0, 1e-12);
            assert_eq!(moving.time_s(), 1.25);
        },
    );
}

#[test]
fn nonprincipal_shear_and_foreign_reference_half_paths_refuse_without_publication() {
    for (paths, shear, message) in [
        ([L / 4.0; 2], 0.003, "principal reference normal"),
        ([L / 3.0, L / 4.0], 0.0, "centroid-to-face"),
    ] {
        with_contact_fixture_paths(
            fixture(HEAT_R0, AREA),
            paths,
            |cx, thermal, solid, field, prepared, bindings, electrodes| {
                let cells = thermal.cells().clone();
                let solution = solid.solution().clone();
                let contacts = field.contacts().to_vec();
                let positions = field.positions_m().to_vec();
                let result = solid.advance_spatial_electrothermal_geometry_candidate(
                    thermal,
                    field,
                    prepared,
                    current_drive(electrodes, CURRENT),
                    bindings,
                    &PRESSURE,
                    &THERMAL,
                    GEOMETRY,
                    1.0,
                    &imposed(FIRST, 0.0, shear),
                    &[],
                    IncrementSettings::default(),
                    contact_controls(),
                    &[0.0; 6],
                    cx,
                );
                let error = match result {
                    Ok(_) => panic!("unsupported thermal path accepted"),
                    Err(error) => error,
                };
                assert!(format!("{error:?}").contains(message), "{error:?}");
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(solid.solution(), &solution);
                assert_eq!(field.contacts(), contacts);
                assert_eq!(field.positions_m(), positions);
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(solid.time_s(), 0.0);
            },
        );
    }
}
