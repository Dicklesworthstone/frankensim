//! G1/G3/G4: original radial volume, spatial Joule heat and material feedback.
//! Manufactured laws verify this coupled model, not measured conductor data.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "axisymmetric_conduction_enthalpy/contacts.rs"]
mod contacts;
#[path = "axisymmetric_conduction_enthalpy/response.rs"]
mod response;
#[path = "axisymmetric_conduction_enthalpy/terminals.rs"]
mod terminals;

use common::{boundary, config, material, material_with_electrical_basis, near, point};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    conduction::{ConductionDrive, ConductionRegion, ConductionThermalNetwork},
    material::{ConductionLink, MaterialCell},
};
use fs_em::conduction::{
    ConductionControls, ConductionCurrentTrace, ConductionPowerSource,
    axisymmetric::AxisymmetricConduction,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::{
    conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis},
    phase::ResolvedHeatCapacityCurve,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use fs_thermal::enthalpy::*;

const REGIONS: [ConductionRegion; 4] = [
    ConductionRegion {
        triangle: 3,
        cell: CellId(10),
    },
    ConductionRegion {
        triangle: 0,
        cell: CellId(30),
    },
    ConductionRegion {
        triangle: 2,
        cell: CellId(10),
    },
    ConductionRegion {
        triangle: 1,
        cell: CellId(30),
    },
];
const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 16,
        max_dofs: 32,
        max_triplets: 144,
    }
}
fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 100,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two axial thermal layers of one cylinder",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
        ],
        vec![[0, 1, 3], [0, 3, 2], [2, 3, 5], [2, 5, 4]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}
fn problem(
    cards: &[MaterialCard; 2],
    bases: [ElectricalPropertyBasis; 2],
    temperature: f64,
    angle: f64,
) -> AxisymmetricConduction {
    let sources = [0, 0, 1, 1]
        .map(|cell| {
            resolve_conductor_state_with_basis(
                &cards[cell],
                &point(temperature),
                MaterialPropertySelection::SingleClaimOnly,
                bases[cell],
            )
            .unwrap()
        })
        .to_vec();
    AxisymmetricConduction::new(mesh(angle), sources, limits(), &mut || false).unwrap()
}
fn initial(curves: &[ResolvedHeatCapacityCurve; 2], angle: f64) -> DepositionSimulation<'_> {
    let area = angle / 2.0; // radius=1; each axial layer has length=1.
    let mut events = Vec::new();
    for (index, id) in [CellId(30), CellId(10)].into_iter().enumerate() {
        let knots = curves[index].curve().knots();
        let cp = (knots.last().unwrap().specific_enthalpy_j_kg - knots[0].specific_enthalpy_j_kg)
            / 1800.0;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: index as u64,
                body_id: index as u64,
                material_id: MaterialId(index as u64),
                mass_kg: area, // authored density=1 kg/m3, actual swept volume.
                specific_h_j_kg: 100.0 * cp,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, index as f64],
                    end_m: [0.0, 0.0, index as f64 + 1.0],
                    volume_m3: area,
                    surface_area_m2: 2.0 * area + angle,
                },
                boundary: boundary(),
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: area,
            conductance_w_k: 999.0, // replaced by the ORIGINAL material link.
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> Vec<MaterialCell<'a>> {
    [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect()
}
fn network<'a>(cells: &'a [MaterialCell<'a>]) -> ConductionThermalNetwork<'a> {
    ConductionThermalNetwork {
        cells,
        links: &LINKS,
        regions: &REGIONS,
        volume_tolerance_m3: 1e-12,
    }
}

#[test]
fn local_current_heating_and_fourier_transfer_keep_the_actual_sector_volume() {
    let alpha = 0.002;
    let (a, ca) = material(1.0, 1.0, alpha);
    let (b, cb) = material_with_electrical_basis(
        ElectricalPropertyBasis::Conductivity,
        1.0 / 3.0,
        2.0,
        0.0,
        1.0,
    );
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    for angle in [0.7, 2.0, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let thermal = initial(&curves, angle);
        let original = problem(
            &cards,
            [
                ElectricalPropertyBasis::Resistivity,
                ElectricalPropertyBasis::Conductivity,
            ],
            300.0,
            angle,
        );
        let currents = [ConductionCurrentTrace {
            vertices: [4, 5],
            inward_current_density_a_m2: 2.0,
        }];
        for drive in [
            ConductionDrive::Current {
                positive: &[4, 5],
                negative: &[0, 1],
                current_a: 2.0 * area,
            },
            ConductionDrive::CurrentDensity {
                voltage: &[(0, 0.0), (1, 0.0)],
                currents: &currents,
            },
        ] {
            let next = network(&bindings)
                .step_axisymmetric_candidate(
                    &thermal,
                    0.5,
                    &original,
                    drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            // Divide the independent two-cell BE balance by swept area:
            // (1.5-2*alpha)dT1-.5dT2=2; -.5dT1+2.5dT2=6.
            let a11 = 1.5 - 2.0 * alpha;
            let determinant = a11 * 2.5 - 0.25;
            let d1 = 8.0 / determinant;
            let d2 = (6.0 * a11 + 1.0) / determinant;
            let rho = 1.0 + alpha * d1;
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + d1,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + d2,
            );
            near(next.report().source.links[0].conductance_w_k, area);
            let field = &next.report().electrical;
            near(field.voltage_v[4], 2.0 * (rho + 3.0));
            near(field.voltage_v[5], field.voltage_v[4]);
            near(field.heat_j, 2.0 * area * (rho + 3.0));
            near(field.supplied_work_j, field.heat_j);
            near(next.report().deposited_joule_j, field.heat_j);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                area * (d1 + 2.0 * d2),
            );
            for triangle in 0..4 {
                let resistivity = if triangle < 2 { rho } else { 3.0 };
                near(field.electric_field_v_m[triangle][0], 0.0);
                near(field.current_density_a_m2[triangle][0], 0.0);
                near(field.current_density_a_m2[triangle][1], -2.0);
                near(
                    field.cell_power_w[triangle],
                    4.0 * resistivity * original.mesh().face_measure(triangle).unwrap(),
                );
                assert_eq!(
                    next.problem().conductors()[triangle]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                    original.conductors()[triangle]
                        .source_property()
                        .answer()
                        .receipt
                        .selected
                );
            }
            // Equal meridional triangle areas have DIFFERENT swept volumes.
            near(field.cell_power_w[0], 2.0 * field.cell_power_w[1]);
            if matches!(drive, ConductionDrive::CurrentDensity { .. }) {
                near(field.nodal_applied_current_a[4].unwrap(), 2.0 * area / 3.0);
                near(field.nodal_applied_current_a[5].unwrap(), 4.0 * area / 3.0);
            }
            near(next.problem().conductors()[0].resistivity_ohm_m(), rho);
            assert_eq!(
                next.problem().conductors()[2].source_basis(),
                ElectricalPropertyBasis::Conductivity
            );
            assert!(
                next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
            );
        }
    }
}

#[test]
fn voltage_and_regulated_power_follow_temperature_and_continue_the_same_history() {
    let alpha = 0.005;
    let (a, ca) = material(1.0, 2.0, alpha);
    let (b, cb) = material(1.0, 2.0, alpha);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let angle = 1.3;
    let area = angle / 2.0;
    let mut thermal = initial(&curves, angle);
    let mut original = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        angle,
    );
    let voltage = [(0, 0.0), (1, 0.0), (4, 4.0), (5, 4.0)];
    let fixed_voltage = network(&bindings)
        .step_axisymmetric_candidate(
            &thermal,
            0.5,
            &original,
            ConductionDrive::Voltage(&voltage),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // delta*(1+alpha*delta)=dt*V^2/(4*cp*rho0)=1.
    let exact = 2.0 / (1.0 + (1.0 + 4.0 * alpha).sqrt());
    for id in [CellId(30), CellId(10)] {
        near(
            fixed_voltage.thermal().temperature_k(id).unwrap(),
            300.0 + exact,
        );
    }
    near(fixed_voltage.report().electrical.heat_j, 4.0 * area * exact);

    let source = ConductionPowerSource {
        power_w: 8.0 * area,
        maximum_voltage_v: 100.0,
        maximum_current_a: 100.0,
    };
    for step in 1..=3 {
        let next = network(&bindings)
            .step_axisymmetric_candidate(
                &thermal,
                step as f64 * 0.5,
                &original,
                ConductionDrive::Power {
                    positive: &[4, 5],
                    negative: &[0, 1],
                    source,
                },
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let rise = step as f64;
        for id in [CellId(30), CellId(10)] {
            near(next.thermal().temperature_k(id).unwrap(), 300.0 + rise);
        }
        near(next.report().electrical.heat_j, 4.0 * area);
        near(next.report().electrical.supplied_work_j, 4.0 * area);
        near(
            next.report().electrical.voltage_v[4],
            4.0 * (1.0 + alpha * rise).sqrt(),
        );
        near(next.thermal().energy_j() - thermal.energy_j(), 4.0 * area);
        (thermal, original, _) = next.into_parts();
    }
}

#[test]
fn invalid_volume_source_compliance_and_late_cancellation_preserve_both_owners() {
    let (a, ca) = material(1.0, 2.0, 0.005);
    let (b, cb) = material(1.0, 2.0, 0.005);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 2.0);
    let original = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        2.0,
    );
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let saved = thermal.cells().clone();
    let sources = original.conductors().to_vec();
    let drive = ConductionDrive::Current {
        positive: &[4, 5],
        negative: &[0, 1],
        current_a: 2.0,
    };
    let mut polls = 0;
    let accepted = layout
        .step_axisymmetric_candidate(
            &thermal,
            0.5,
            &original,
            drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        layout.step_axisymmetric_candidate(
            &thermal,
            0.5,
            &original,
            drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));

    for bad in [
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            3.0,
        ),
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            301.0,
            2.0,
        ),
    ] {
        assert!(
            layout
                .step_axisymmetric_candidate(
                    &thermal,
                    0.5,
                    &bad,
                    drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .is_err()
        );
    }
    // Initial P=8 W needs 4 V; the hot endpoint needs sqrt(16.08)>4.005 V.
    assert!(
        layout
            .step_axisymmetric_candidate(
                &thermal,
                0.5,
                &original,
                ConductionDrive::Power {
                    positive: &[4, 5],
                    negative: &[0, 1],
                    source: ConductionPowerSource {
                        power_w: 8.0,
                        maximum_voltage_v: 4.005,
                        maximum_current_a: 100.0,
                    },
                },
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    let mut controls = config();
    controls.maximum_circuit_trials = 1;
    assert!(
        layout
            .step_axisymmetric_candidate(
                &thermal,
                0.5,
                &original,
                drive,
                limits(),
                field_controls(),
                controls,
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.conductors(), &sources);
    let retry = layout
        .step_axisymmetric_candidate(
            &thermal,
            0.5,
            &original,
            drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
    assert_eq!(
        retry.report().electrical.voltage_v,
        accepted.report().electrical.voltage_v
    );
}
