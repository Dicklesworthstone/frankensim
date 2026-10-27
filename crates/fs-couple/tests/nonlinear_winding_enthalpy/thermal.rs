//! G1/G3/G4: voltage-driven mu(T), reversible field heat and original enthalpy.
//! The two square-star modes reduce independently to B=mu(T)*I. Synthetic
//! constitutive data exercise coupling; they are not measured coil validation.

use super::*;
#[path = "../../../fs-em/tests/winding_ohmic/common.rs"]
mod linear;

use fs_couple::circuit_enthalpy::eddy::thermal::MagneticThermalLaw;
use fs_em::{
    planar::nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    winding::{
        nonlinear::{NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindings},
        ohmic::WindingConductor,
    },
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::thermal::magnetic::ThermalMagneticMaterial,
};
use fs_phs::circuit::{Branch, CircuitGraph};

const VOLTAGE: f64 = 20.0;
const END_H: f64 = 0.2;
const HEAT_CAPACITY: f64 = 4.0 * VOLUME;
const RHO_SLOPE: f64 = 0.005;

fn card(slope: f64) -> MaterialCard {
    if slope == 0.0 {
        linear::card(
            1.0 / 18.0,
            RHO_SLOPE,
            1.0,
            ElectricalPropertyBasis::Resistivity,
        )
    } else {
        linear::affine_magnetic_card(
            1.0 / 18.0,
            RHO_SLOPE,
            1.0,
            slope,
            ElectricalPropertyBasis::Resistivity,
        )
    }
}

fn problem(card: &MaterialCard) -> PlanarNonlinearMagnetostatics {
    // Reuse the existing manufactured two-square geometry and winding profile.
    let geometry = field::problem(&field::source(false, 18.0), DEPTH);
    let material = linear::bundle_at(
        card,
        &linear::absolute_point(300.0),
        ElectricalPropertyBasis::Resistivity,
    );
    PlanarNonlinearMagnetostatics::new(
        geometry.mesh().clone(),
        vec![IsotropicMagneticMaterial::Linear(material); 8],
        linear::limits(),
        &mut || false,
    )
    .unwrap()
}

fn circuit(problem: &PlanarNonlinearMagnetostatics, initial: f64) -> NonlinearWindingCircuit<'_> {
    let profiles = field::profiles(1);
    let conductors = problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, law)| WindingConductor {
            triangle,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                law.material().clone(),
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap(),
        })
        .collect();
    NonlinearWindingCircuit::new(
        field::initial_field(problem, &profiles, &[initial]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: SERIES_OHM }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![END_H],
        },
        conductors,
        NonlinearWindingInitial {
            currents_a: vec![initial],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![VOLTAGE],
        },
        field::config(),
        &mut || false,
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &linear::absolute_point(200.0),
            upper: &linear::absolute_point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}

fn thermal(curve: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MATERIAL,
                mass_kg: VOLUME,
                specific_h_j_kg: 4.0 * (300.0 - 200.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [2.0, 0.0, 0.0],
                    volume_m3: VOLUME,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn laws(card: &MaterialCard, problem: &PlanarNonlinearMagnetostatics) -> Vec<MagneticThermalLaw> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, material)| MagneticThermalLaw {
            triangle,
            law: ThermalMagneticMaterial::new(card, material.material(), "T", [290.0, 360.0])
                .unwrap(),
            maximum_flux_density_t: 15.0,
        })
        .collect()
}

fn mu(t: f64, slope: f64) -> f64 {
    1.0 + slope * (t - 300.0)
}
fn field_u(t: f64, b: f64, slope: f64) -> f64 {
    VOLUME * b * b / (2.0 * mu(t, slope)) * (1.0 + t * slope / mu(t, slope))
}
fn resistance(t: f64) -> f64 {
    72.0 * DEPTH / 18.0 * (1.0 + RHO_SLOPE * (t - 300.0)) + SERIES_OHM
}

struct Reference {
    temperature: f64,
    current: f64,
    b: f64,
    stored: f64,
    heat: f64,
    numerical: f64,
    parameter: f64,
    reversible: f64,
}

fn reference(i0: f64, t0: f64, slope: f64) -> Reference {
    let b0 = mu(t0, slope) * i0;
    let current = |t| {
        (2.0 * DEPTH * b0 + END_H * i0 + DT * VOLTAGE)
            / (2.0 * DEPTH * mu(t, slope) + END_H + DT * resistance(t))
    };
    let numerical = |t, i| {
        VOLUME * (mu(t, slope) * i - b0).powi(2) / (2.0 * mu(t, slope))
            + 0.5 * END_H * (i - i0).powi(2)
    };
    let residual = |t| {
        let i = current(t);
        HEAT_CAPACITY * (t - t0) + field_u(t, mu(t, slope) * i, slope) - field_u(t0, b0, slope)
            + 0.5 * END_H * (i * i - i0 * i0)
            + numerical(t, i)
            - DT * VOLTAGE * i
    };
    let (mut lo, mut hi) = (290.0, 360.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..90 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let t = f64::midpoint(lo, hi);
    let i = current(t);
    let b = mu(t, slope) * i;
    let parameter = VOLUME * b0 * b0 * (1.0 / mu(t, slope) - 1.0 / mu(t0, slope)) / 2.0;
    let delta_ts = VOLUME * slope / 2.0
        * (t * b * b / mu(t, slope).powi(2) - t0 * b0 * b0 / mu(t0, slope).powi(2));
    Reference {
        temperature: t,
        current: i,
        b,
        stored: VOLUME * b * b / (2.0 * mu(t, slope)) + 0.5 * END_H * i * i,
        heat: DT * resistance(t) * i * i,
        numerical: numerical(t, i),
        parameter,
        reversible: -parameter - delta_ts,
    }
}

#[test]
fn permeability_data_changes_voltage_current_and_thermal_history_with_physical_energy_closure() {
    let mut outcomes = Vec::new();
    for slope in [0.0, 0.004] {
        let card = card(slope);
        let curve = caloric(&card);
        let mut thermal = thermal(&curve);
        let problem = problem(&card);
        let mut circuit = circuit(&problem, 0.0);
        let bindings = laws(&card, &problem);
        let context = linear::absolute_point(300.0);
        let cells = cells(&card, &curve, &context);
        let regions = routes();
        let layout = network(&cells, &regions);
        let mut supplied = 0.0;
        let mut numerical = 0.0;
        for n in 1..=2 {
            let expected = reference(
                circuit.currents_a()[0],
                thermal.temperature_k(CELL).unwrap(),
                slope,
            );
            let candidate = layout
                .step_nonlinear_thermal_candidate(
                    &thermal,
                    n as f64 * DT,
                    &circuit,
                    &[VOLTAGE],
                    &bindings,
                    None,
                    controls(),
                    field::controls(),
                    linear::limits(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let report = candidate.report();
            let magnetic = report.source.magnetic_caloric.as_ref().unwrap();
            let balance = candidate.caloric_balance().unwrap();
            near(
                candidate.thermal().temperature_k(CELL).unwrap(),
                expected.temperature,
            );
            near(candidate.circuit().currents_a()[0], expected.current);
            near(report.electrical.stored_after_j, expected.stored);
            near(report.electrical.heat_j, expected.heat);
            near(
                report.electrical.numerical_dissipation_j,
                expected.numerical,
            );
            near(report.electrical.parameter_work_j, expected.parameter);
            near(magnetic.reversible_heat_j, expected.reversible);
            near(report.deposited_joule_j, expected.heat);
            near(report.non_electrical_heat_j, expected.reversible);
            near(balance.external_heat_j, 0.0);
            near(balance.total_energy_residual_j, 0.0);
            near(
                report.thermal.stored_energy_change_j,
                expected.heat + expected.reversible,
            );
            for sample in &magnetic.fields {
                assert_eq!(sample.cell, CELL);
                near(sample.flux_density_t[1], expected.b);
                near(sample.temperature_k[1], expected.temperature);
            }
            near(
                magnetic.field_internal_energy_j[1],
                field_u(expected.temperature, expected.b, slope),
            );
            supplied += DT * VOLTAGE * expected.current;
            numerical += expected.numerical;
            near(
                field_u(expected.temperature, expected.b, slope)
                    + 0.5 * END_H * expected.current.powi(2)
                    + HEAT_CAPACITY * (expected.temperature - 300.0)
                    + numerical,
                supplied,
            );
            // Continue the exact accepted field, original circuit and caloric
            // history. The next interval must not restart the magnetic state.
            let (next_thermal, next_circuit, _) = candidate.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
        }
        outcomes.push((
            circuit.currents_a()[0],
            thermal.temperature_k(CELL).unwrap(),
        ));
    }
    assert!((outcomes[0].0 - outcomes[1].0).abs() > 1e-3);
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 0.1);
}

#[test]
fn constant_permeability_matches_existing_solver_and_replays_the_same_endpoint() {
    let card = card(0.0);
    let curve = caloric(&card);
    let thermal = thermal(&curve);
    let problem = problem(&card);
    let circuit = circuit(&problem, 1.0);
    let bindings = laws(&card, &problem);
    let context = linear::absolute_point(300.0);
    let cells = cells(&card, &curve, &context);
    let regions = routes();
    let layout = network(&cells, &regions);
    let run = || {
        layout
            .step_nonlinear_thermal_candidate(
                &thermal,
                DT,
                &circuit,
                &[VOLTAGE],
                &bindings,
                None,
                controls(),
                field::controls(),
                linear::limits(),
                &[0.0],
                || false,
            )
            .unwrap()
    };
    let first = run();
    let replay = run();
    let fixed = layout
        .step_nonlinear_candidate(
            &thermal,
            DT,
            &circuit,
            &[VOLTAGE],
            controls(),
            field::controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(first.circuit().currents_a(), replay.circuit().currents_a());
    assert_eq!(
        first.circuit().field().potential_wb_m(),
        replay.circuit().field().potential_wb_m()
    );
    assert_eq!(first.thermal().cells(), replay.thermal().cells());
    near(
        first.circuit().currents_a()[0],
        fixed.circuit().currents_a()[0],
    );
    near(
        first.thermal().temperature_k(CELL).unwrap(),
        fixed.thermal().temperature_k(CELL).unwrap(),
    );
    near(
        first
            .report()
            .source
            .magnetic_caloric
            .as_ref()
            .unwrap()
            .reversible_heat_j,
        0.0,
    );
    near(first.report().electrical.parameter_work_j, 0.0);
    assert_eq!(circuit.time_s(), 0.0);
    near(thermal.temperature_k(CELL).unwrap(), 300.0);
}

#[test]
fn missing_source_capacity_envelope_and_final_cancellation_preserve_all_histories() {
    let card = card(0.004);
    let curve = caloric(&card);
    let thermal = thermal(&curve);
    let problem = problem(&card);
    let circuit = circuit(&problem, 1.0);
    let before_field = circuit.field().potential_wb_m().to_vec();
    let before_cells = thermal.cells().clone();
    let bindings = laws(&card, &problem);
    let context = linear::absolute_point(300.0);
    let cells = cells(&card, &curve, &context);
    let regions = routes();
    let layout = network(&cells, &regions);
    let run = |bindings: &[MagneticThermalLaw]| {
        layout.step_nonlinear_thermal_candidate(
            &thermal,
            DT,
            &circuit,
            &[VOLTAGE],
            bindings,
            None,
            controls(),
            field::controls(),
            linear::limits(),
            &[0.0],
            || false,
        )
    };
    assert!(run(&bindings[..7]).is_err());
    let mut changed = bindings.clone();
    changed[0].maximum_flux_density_t = 0.9;
    assert!(run(&changed).is_err());
    let mut changed = bindings.clone();
    for binding in &mut changed {
        binding.maximum_flux_density_t = 100.0;
    }
    assert!(run(&changed).is_err());
    let foreign = self::card(0.003);
    assert!(run(&laws(&foreign, &self::problem(&foreign))).is_err());
    let mut checkpoints = 0;
    let success = layout
        .step_nonlinear_thermal_candidate(
            &thermal,
            DT,
            &circuit,
            &[VOLTAGE],
            &bindings,
            None,
            controls(),
            field::controls(),
            linear::limits(),
            &[0.0],
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut visits = 0;
    assert!(matches!(
        layout.step_nonlinear_thermal_candidate(
            &thermal,
            DT,
            &circuit,
            &[VOLTAGE],
            &bindings,
            None,
            controls(),
            field::controls(),
            linear::limits(),
            &[0.0],
            || {
                visits += 1;
                visits == checkpoints
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(success.circuit().time_s(), DT);
    assert_eq!(circuit.time_s(), 0.0);
    assert_eq!(circuit.currents_a(), &[1.0]);
    assert_eq!(circuit.field().potential_wb_m(), before_field);
    assert_eq!(thermal.cells(), &before_cells);
    assert_eq!(thermal.time_s(), 0.0);
}
