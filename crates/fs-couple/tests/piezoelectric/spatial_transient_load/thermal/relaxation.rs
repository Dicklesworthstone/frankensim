//! Source-bound ceramic heat, separate conductor destinations and atomic thermal steps.
#[path = "relaxation/support_motion.rs"]
mod support_motion;
#[path = "relaxation/acoustic.rs"]
mod acoustic;

use super::*;
use crate::harmonic::relaxation::fixture as memory_fixture;
use fs_couple::piezoelectric::spatial::transient::load::thermal::ThermalCeramicCell;
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;
use fs_qty::Time;

const CERAMIC_DENSITY_KG_M3: f64 = 7500.0;
const CERAMIC_CP_J_KG_K: f64 = 1.0;
const POLE_INCREMENT_F_M: f64 = 0.8 * EPSILON[2];
const POLE_TIME_S: f64 = 1e-6;
const MEMORY_BUDGET: usize = 256;

fn ceramic(
    temperature_dependent_time: bool,
) -> (
    MaterialCard,
    ResolvedHeatCapacityCurve,
    TensorPiezoelectric,
    RelaxingTensorPiezoelectric,
) {
    ceramic_with_mechanical_memory(temperature_dependent_time, None)
}

fn ceramic_with_mechanical_memory(
    temperature_dependent_time: bool,
    mechanical: Option<(f64, f64)>,
) -> (
    MaterialCard,
    ResolvedHeatCapacityCurve,
    TensorPiezoelectric,
    RelaxingTensorPiezoelectric,
) {
    let (reference, memory) =
        memory_fixture(mechanical, Some((POLE_INCREMENT_F_M, POLE_TIME_S)), FRAME);
    let mut claims = ClaimSet::new();
    // Declare one complete synthetic card; all coupled owners then select
    // their original roles from that SAME card and typed physical context.
    for property in reference
        .material()
        .properties()
        .iter()
        .chain(memory.material().properties())
    {
        let name = property.requirement().name();
        let dims = property.requirement().dims();
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar {
                value: property.value_si(),
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored ceramic memory and caloric fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        if temperature_dependent_time && name == "memory_tau_p" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.5 * POLE_TIME_S), (400.0, 3.5 * POLE_TIME_S)],
                dims: Time::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            value: PropertyValue::Scalar {
                value: CERAMIC_CP_J_KG_K,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored constant ceramic caloric fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal ceramic with caloric storage".into(),
            phase: "solid".into(),
            process: "explicit frozen permittivity, polarization memory and finite heat capacity"
                .into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    let base = TensorPiezoelectric::resolve(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    let loss = RelaxingTensorPiezoelectric::resolve(
        base.clone(),
        &card,
        memory.mechanical_properties().to_vec(),
        memory.dielectric_properties().to_vec(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, curve, base, loss)
}

fn isolated_cells(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let masses = [CERAMIC_DENSITY_KG_M3 * VOLUME, VOLUME];
    let specific_enthalpies = [20.0 * CERAMIC_CP_J_KG_K, 40.0];
    let events = CELLS
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: masses[i],
                specific_h_j_kg: specific_enthalpies[i],
                geometry: CellGeometry {
                    start_m: [i as f64 * LENGTH, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * LENGTH, 0.0, 0.0],
                    volume_m3: VOLUME,
                    surface_area_m2: 6.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-7,
            energy_tolerance_j: 1e-16,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn caloric_cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
) -> [MaterialCell<'a>; 2] {
    core::array::from_fn(|i| MaterialCell {
        cell: CELLS[i],
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    })
}

fn ceramic_routes(tetrahedra: usize) -> Vec<ThermalCeramicCell> {
    (0..tetrahedra)
        .map(|piezoelectric_cell| ThermalCeramicCell {
            piezoelectric_cell,
            thermal_cell: CELLS[0],
        })
        .collect()
}

fn heat_journey(with_resistor: bool) {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (ceramic_card, ceramic_curve, base, loss) = ceramic(false);
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            field_controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![loss], MEMORY_BUDGET, field_controls(), cx)
        .unwrap();
        let (resistor_card, resistor_curve) = thermal_material(2.0, 4.0, 0.0, 400.0);
        let cards = [ceramic_card, resistor_card];
        let curves = [ceramic_curve, resistor_curve];
        let context = point(300.0);
        let cells = caloric_cells(&cards, &curves, &context);
        let resistors = if with_resistor {
            vec![CellResistor {
                destination: ResistorCell {
                    branch: LOAD.resistor_branch,
                    cell: CELLS[1],
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
                length_m: LENGTH,
                area_m2: AREA,
            }]
        } else {
            vec![]
        };
        let bindings = if with_resistor { vec![LOAD] } else { vec![] };
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut thermal = isolated_cells(&curves);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let routes = ceramic_routes(cube.tets.len());
        let initial_energy = thermal.energy_j() + body.state().total_energy_j();
        let cp = POLE_INCREMENT_F_M * LENGTH;
        let dt = 0.4 * POLE_TIME_S;
        let h = 0.5 * dt;
        let conductance = if with_resistor { 1.0 / RESISTANCE } else { 0.0 };
        let (mut voltage, mut polarization) = (0.0, 0.0);
        let (mut ceramic_heat, mut resistor_heat, mut supplied) = (0.0, 0.0, 0.0);
        for current in [3e-3, 1e-3, 0.0] {
            // Independently eliminate C*Vdot + Pdot + V/R = I and
            // tau*Pdot + P = Cp*V. Thermal cells integrate these actual losses.
            let vm = (CAPACITANCE * voltage / h + current + polarization / (POLE_TIME_S + h))
                / (CAPACITANCE / h + conductance + cp / (POLE_TIME_S + h));
            let pm = (POLE_TIME_S * polarization + h * cp * vm) / (POLE_TIME_S + h);
            let next_voltage = 2.0 * vm - voltage;
            let next_polarization = 2.0 * pm - polarization;
            let heat = POLE_TIME_S * (next_polarization - polarization).powi(2) / (cp * dt);
            let joule = dt * conductance * vm * vm;
            let stored = 0.5 * CAPACITANCE * next_voltage * next_voltage
                + 0.5 * next_polarization * next_polarization / cp;
            let target = thermal.time_s() + dt;
            let before_t = CELLS.map(|cell| thermal.temperature_k(cell).unwrap());
            let before_time = body.state().time_s();
            let before_memory = body.state().material_memory_energy_j();
            let candidate = body
                .step_with_material_thermal_loads_candidate(
                    &thermal,
                    target,
                    &prepared,
                    &bindings,
                    &routes,
                    "T",
                    &[Drive::Voltage(0.0), Drive::Current(current)],
                    &zero,
                    &charge,
                    field_controls(),
                    thermal_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let report = candidate.report();
            near(report.electrical.native.material_heat_j, heat, 2e-16, 3e-8);
            near(
                report.electrical.native.dielectric_material_heat_j,
                heat,
                2e-16,
                3e-8,
            );
            assert_eq!(report.electrical.native.mechanical_material_heat_j, 0.0);
            near(report.electrical.loads.joule_heat_j, joule, 2e-16, 3e-8);
            near(report.deposited_joule_j, heat + joule, 2e-16, 3e-8);
            near(report.non_electrical_heat_j, 0.0, 2e-16, 0.0);
            near(report.total_energy_residual_j, 0.0, 2e-16, 0.0);
            assert_eq!(report.electrical.ceramic_destinations, routes);
            assert_eq!(
                report.electrical.destinations,
                resistors.iter().map(|r| r.destination).collect::<Vec<_>>()
            );
            let sample = report
                .source
                .cells
                .iter()
                .find(|sample| sample.cell == CELLS[0])
                .unwrap();
            assert_eq!(
                sample.phase.material_card_identity(),
                cards[0].content_hash()
            );
            near(
                candidate.body().state().electrodes()[1].voltage_v,
                next_voltage,
                2e-8,
                3e-8,
            );
            near(
                candidate.body().state().electrodes()[1].charge_c,
                CAPACITANCE * next_voltage + next_polarization,
                2e-18,
                3e-8,
            );
            near(
                candidate.body().state().stored_energy_j(),
                stored,
                2e-16,
                3e-8,
            );
            near(
                report.electrical.loads.electrical_source_work_j,
                dt * current * vm,
                2e-16,
                3e-8,
            );
            ceramic_heat += heat;
            resistor_heat += joule;
            supplied += dt * current * vm;
            for (cell, expected_temperature) in [
                (
                    CELLS[0],
                    300.0 + ceramic_heat / (CERAMIC_DENSITY_KG_M3 * VOLUME * CERAMIC_CP_J_KG_K),
                ),
                (CELLS[1], 300.0 + resistor_heat / (2.0 * VOLUME)),
            ] {
                near(
                    candidate.thermal().temperature_k(cell).unwrap(),
                    expected_temperature,
                    2e-6,
                    0.0,
                );
            }
            near(
                candidate.thermal().energy_j() + candidate.body().state().total_energy_j(),
                initial_energy + supplied,
                4e-16,
                0.0,
            );
            // Discarded Newton proposals never heat or advance the accepted histories.
            assert_eq!(body.state().time_s(), before_time);
            assert_eq!(body.state().material_memory_energy_j(), before_memory);
            assert_eq!(
                CELLS.map(|cell| thermal.temperature_k(cell).unwrap()),
                before_t
            );
            assert_eq!(candidate.thermal().time_s(), target);
            assert_eq!(candidate.body().state().time_s(), target);
            (thermal, body, _) = candidate.into_parts();
            voltage = next_voltage;
            polarization = next_polarization;
        }
        assert!(ceramic_heat > 0.0 && thermal.temperature_k(CELLS[0]).unwrap() > 300.0);
        if with_resistor {
            assert!(resistor_heat > 0.0 && thermal.temperature_k(CELLS[1]).unwrap() > 300.0);
        } else {
            assert_eq!(thermal.temperature_k(CELLS[1]).unwrap(), 300.0);
        }
    });
}

#[test]
fn g1_current_driven_debye_heat_reaches_the_named_finite_ceramic_cell() {
    heat_journey(false);
}

#[test]
fn g1_ceramic_and_resistor_heat_reach_distinct_cells_and_close_source_work() {
    heat_journey(true);
}

#[test]
fn g4_bad_routes_and_temperature_dependent_ceramic_sources_leave_both_histories_retryable() {
    with_cx(|cx, _| {
        for changing_source in [false, true] {
            let cube = cube(1);
            let n = cube.nodes.len();
            let (ceramic_card, ceramic_curve, base, loss) = ceramic(changing_source);
            assert!(loss.frequency_coefficients(100.0).is_ok());
            let materials = [base];
            let model =
                Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
            let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
            let zero = vec![[0.0; 3]; n];
            let charge = vec![0.0; n];
            let body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &zero,
                &zero,
                &charge,
                field_controls(),
                cx,
            )
            .unwrap()
            .with_relaxation(vec![loss], MEMORY_BUDGET, field_controls(), cx)
            .unwrap();
            let (resistor_card, resistor_curve) = thermal_material(2.0, 4.0, 0.0, 400.0);
            let cards = [ceramic_card, resistor_card];
            let curves = [ceramic_curve, resistor_curve];
            let context = point(300.0);
            let cells = caloric_cells(&cards, &curves, &context);
            let network = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let thermal = isolated_cells(&curves);
            let prepared = network.prepare(&thermal, || false).unwrap();
            let routes = ceramic_routes(cube.tets.len());
            let advance = |routes: &[ThermalCeramicCell]| {
                body.step_with_material_thermal_loads_candidate(
                    &thermal,
                    0.4 * POLE_TIME_S,
                    &prepared,
                    &[],
                    routes,
                    "T",
                    &[Drive::Voltage(0.0), Drive::Current(3e-3)],
                    &zero,
                    &charge,
                    field_controls(),
                    thermal_controls(),
                    &[0.0; 2],
                    cx,
                )
            };
            if changing_source {
                // tau(T) has the original value at 300 K and is a valid fixed-T
                // material. Heating cannot silently freeze that varying source.
                assert!(advance(&routes).is_err());
            } else {
                let mut duplicate = routes.clone();
                duplicate[1] = duplicate[0];
                let mut outside_mesh = routes.clone();
                outside_mesh[0].piezoelectric_cell = cube.tets.len();
                let mut inactive = routes.clone();
                inactive[0].thermal_cell = CellId(999);
                let wrong_material = routes
                    .iter()
                    .map(|route| ThermalCeramicCell {
                        thermal_cell: CELLS[1],
                        ..*route
                    })
                    .collect::<Vec<_>>();
                for bad in [
                    routes[..routes.len() - 1].to_vec(),
                    duplicate,
                    outside_mesh,
                    inactive,
                    wrong_material,
                ] {
                    assert!(advance(&bad).is_err());
                }
                let first = advance(&routes).unwrap();
                let retry = advance(&routes).unwrap();
                assert!(retry.thermal().temperature_k(CELLS[0]).unwrap() > 300.0);
                assert_eq!(
                    retry.body().state().material_memory_energy_j(),
                    first.body().state().material_memory_energy_j()
                );
                assert_eq!(retry.thermal().energy_j(), first.thermal().energy_j());
            }
            assert_eq!(body.state().time_s(), 0.0);
            assert_eq!(body.state().stored_energy_j(), 0.0);
            assert_eq!(body.state().material_memory_energy_j(), 0.0);
            assert_eq!(body.state().electrodes()[1].charge_c, 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(
                CELLS.map(|cell| thermal.temperature_k(cell).unwrap()),
                [300.0; 2]
            );
        }
    });
}
