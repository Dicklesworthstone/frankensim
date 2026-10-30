use super::*;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::{
    ResistorCell,
    material::{CellResistor, MaterialCell},
};
use crate::material_circuit::MaterialCircuitSolve;
use crate::speaker::{PanelMode, TsCard, VoiceCoilMagnetics};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const COIL: CellId = CellId(41);
const MASS: f64 = 0.001;

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}

// Authored, deliberately scaled analytic fixtures; not qualified copper data.
fn material(
    alpha: f64,
    basis: ElectricalPropertyBasis,
    upper: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    material_with_properties(alpha, basis, upper, 1.0, 1.0)
}

fn material_with_properties(
    alpha: f64,
    basis: ElectricalPropertyBasis,
    upper: f64,
    density_kg_m3: f64,
    resistance_scale: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new(basis.property_name(), basis.dimensions());
    // length/area = 10 m^-1, so either source gives R(300 K)=0.8 ohm.
    let base = if basis == ElectricalPropertyBasis::Resistivity {
        0.08 * resistance_scale
    } else {
        12.5 / resistance_scale
    };
    let mut claims = ClaimSet::new();
    for (key, a, b, limit) in [
        (heat.clone(), 1.0, 1.0, 2000.0),
        (density.clone(), density_kg_m3, density_kg_m3, 2000.0),
        (
            electrical,
            base * (1.0 - 100.0 * alpha),
            base * (1.0 + 1700.0 * alpha),
            upper,
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, a), (2000.0, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    limit,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic coil thermal fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}
fn initial(
    curve: &ResolvedHeatCapacityCurve,
    temperature: f64,
    htc: f64,
) -> DepositionSimulation<'_> {
    let density = curve
        .curve()
        .state_at_specific_enthalpy(temperature - 200.0)
        .unwrap()
        .bulk_density_kg_m3();
    initial_with_mass(curve, temperature, htc, MASS * density)
}

fn initial_with_mass(
    curve: &ResolvedHeatCapacityCurve,
    temperature: f64,
    htc: f64,
    mass_kg: f64,
) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: COIL,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg,
                specific_h_j_kg: temperature - 200.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: MASS,
                    surface_area_m2: 0.06,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: htc,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bindings<'a>(
    material: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> ([MaterialCell<'a>; 1], [CellResistor; 1]) {
    (
        [MaterialCell {
            cell: COIL,
            material_id: MaterialId(0),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }],
        [CellResistor {
            destination: ResistorCell {
                branch: 0,
                cell: COIL,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 0.1,
            area_m2: 0.01,
        }],
    )
}
fn driver(resistance: f64, nonlinear: bool, cabinet: bool) -> MagneticTsDriver {
    let card = TsCard {
        re_ohm: resistance,
        le_h: 0.006,
        bl: 1.2,
        mms_kg: 0.7,
        cms_m_per_n: 0.9,
        rms_n_s_m: 0.1,
        sd_m2: 0.001,
        x_max_m: 0.1,
    };
    let profile = VoiceCoilMagnetics::new(
        [-0.05, 0.05],
        0.01,
        if nonlinear {
            [1.2, -0.4, -0.2]
        } else {
            [1.2, 0.0, 0.0]
        },
        if nonlinear {
            [0.006, 0.0015, 0.0008]
        } else {
            [0.006, 0.0, 0.0]
        },
    )
    .unwrap();
    MagneticTsDriver::new(
        card,
        cabinet.then_some(0.02),
        cabinet.then_some(PanelMode {
            frequency_hz: 2.0,
            mass_kg: 0.05,
            area_m2: 0.001,
            damping_ratio: 0.02,
        }),
        profile,
    )
    .unwrap()
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 200,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-6,
        cell_energy_tolerance_j: 1e-9,
        thermal_energy_tolerance_j: 1e-9,
        total_energy_tolerance_j: 2e-9,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-10,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}

#[test]
fn self_heating_refreshes_actual_resistance_and_compresses_magnetic_current() {
    let (card, curve) = material(0.005, ElectricalPropertyBasis::Resistivity, 2000.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = initial(&curve, 300.0, 0.0);
    let mut coupled = ThermalMagneticTsDriver::bind(
        driver(0.8, false, false),
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let mut cold = driver(0.8, false, false);
    let start_energy = thermal.energy_j();
    let (mut supplied, mut external_loss, mut heat) = (0.0, 0.0, 0.0);
    for tick in 1..=20 {
        cold.step(0.4, 0.005).unwrap();
        let proposal = coupled
            .step_candidate(
                &thermal,
                f64::from(tick) * 0.005,
                0.4,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let temperature = proposal.thermal().temperature_k(COIL).unwrap();
        let report = proposal.report();
        let step = &report.electrical;
        near(
            step.coil.resistance_ohm(),
            0.8 * (1.0 + 0.005 * (temperature - 300.0)),
            2e-14,
        );
        near(
            proposal.driver().driver().coil_resistance_ohm(),
            step.coil.resistance_ohm(),
            0.0,
        );
        near(
            step.coil_heat_j,
            0.005 * step.coil.resistance_ohm() * step.native.y[0].powi(2),
            2e-15,
        );
        near(report.deposited_joule_j, step.coil_heat_j, 1e-15);
        near(report.total_energy_residual_j, 0.0, 2e-9);
        supplied += step.native.supplied;
        external_loss += step.mechanical_acoustic_loss_j;
        heat += step.coil_heat_j;
        let (t, d, _) = proposal.into_parts();
        thermal = t;
        coupled = d;
    }
    assert!(thermal.temperature_k(COIL).unwrap() > 310.0);
    assert!(coupled.driver().coil_resistance_ohm() > 0.84);
    assert!(coupled.driver().coil_current_a() < cold.coil_current_a());
    assert!(coupled.driver().magnetic_force_n() < cold.magnetic_force_n());
    near(thermal.energy_j() - start_energy, heat, 2e-9);
    near(
        coupled.driver().stored_energy_j() + heat + external_loss,
        supplied,
        3e-9,
    );
    assert_eq!(coupled.driver().time_s(), thermal.time_s());
}

#[test]
fn constant_source_reduces_to_original_nonlinear_speaker_without_heating_mechanical_loss() {
    // Exercise the explicit conductivity source route and all five cabinet states.
    let (card, curve) = material(0.0, ElectricalPropertyBasis::Conductivity, 2000.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = initial(&curve, 300.0, 0.0);
    let mut coupled = ThermalMagneticTsDriver::bind(
        driver(0.8, true, true),
        &thermal,
        &network,
        ElectricalPropertyBasis::Conductivity,
        || false,
    )
    .unwrap();
    let mut fixed = driver(0.8, true, true);
    let start = thermal.energy_j();
    let (mut coil_heat, mut other_loss) = (0.0, 0.0);
    for tick in 1..=12 {
        let (reference, pressure) = fixed.step(0.4, 0.005).unwrap();
        let proposal = coupled
            .step_candidate(
                &thermal,
                f64::from(tick) * 0.005,
                0.4,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let step = &proposal.report().electrical;
        for (&a, &b) in step.native.x.iter().zip(&reference.x) {
            near(a, b, 3e-13);
        }
        near(step.pressure_pa, pressure, 1e-12);
        near(
            step.native.dissipated,
            step.coil_heat_j + step.mechanical_acoustic_loss_j,
            2e-15,
        );
        assert_eq!(
            step.coil.source_basis(),
            ElectricalPropertyBasis::Conductivity
        );
        coil_heat += step.coil_heat_j;
        other_loss += step.mechanical_acoustic_loss_j;
        let (t, d, _) = proposal.into_parts();
        thermal = t;
        coupled = d;
    }
    assert!(other_loss > 1e-6);
    near(thermal.energy_j() - start, coil_heat, 2e-9);
    assert!((thermal.energy_j() - start - coil_heat - other_loss).abs() > 1e-6);
}

#[test]
fn authored_cooling_uses_original_finite_thermal_store_and_source_temperature() {
    let (card, curve) = material(0.005, ElectricalPropertyBasis::Resistivity, 2000.0);
    let context = point(400.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 400.0, 1.0);
    let coupled = ThermalMagneticTsDriver::bind(
        driver(1.2, false, false),
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let proposal = coupled
        .step_candidate(&thermal, 0.01, 0.0, controls(), &[0.0], || false)
        .unwrap();
    // Independent backward-Euler cooling solution: C dT/dt=-h A(T-Tamb).
    let expected = (MASS * 400.0 + 0.01 * 0.06 * 300.0) / (MASS + 0.01 * 0.06);
    near(
        proposal.thermal().temperature_k(COIL).unwrap(),
        expected,
        2e-8,
    );
    near(
        proposal.driver().coil().resistance_ohm(),
        0.8 * (1.0 + 0.005 * (expected - 300.0)),
        1e-10,
    );
    near(proposal.report().electrical.coil_heat_j, 0.0, 0.0);
    near(proposal.driver().driver().stored_energy_j(), 0.0, 0.0);
    assert!(proposal.report().non_electrical_heat_j < 0.0);
    near(proposal.report().total_energy_residual_j, 0.0, 2e-9);
}

#[test]
fn source_failure_cancellation_and_work_exhaustion_publish_neither_history() {
    let (card, curve) = material(0.005, ElectricalPropertyBasis::Resistivity, 301.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 300.0, 0.0);
    let coupled = ThermalMagneticTsDriver::bind(
        driver(0.8, false, false),
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let saved = thermal.cells().clone();
    let state = coupled.driver().x.clone();
    // Heating requests a temperature beyond the pinned electrical source;
    // no extrapolation, stale resistance or cold-coil fallback is accepted.
    assert!(matches!(
        coupled.step_candidate(&thermal, 0.01, 2.0, controls(), &[0.0], || false),
        Err(CircuitEnthalpyError::Source(_))
    ));
    assert!(matches!(
        coupled.step_candidate(&thermal, 0.01, 0.1, controls(), &[0.0], || true),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut budget = controls();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        coupled.step_candidate(&thermal, 0.01, 0.1, budget, &[0.0], || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(coupled.driver().x, state);
    assert_eq!(coupled.driver().time_s(), 0.0);
    near(coupled.coil().resistance_ohm(), 0.8, 2e-15);
}

fn source_mechanics() -> SpeakerMechanics {
    SpeakerMechanics {
        non_coil_moving_mass_kg: 0.002,
        cms_m_per_n: 0.9,
        rms_n_s_m: 0.1,
        sd_m2: 0.001,
        x_max_m: 0.1,
    }
}

fn source_profile() -> VoiceCoilMagnetics {
    VoiceCoilMagnetics::new([-0.05, 0.05], 0.01, [1.2, 0.0, 0.0], [0.006, 0.0, 0.0]).unwrap()
}

#[test]
fn material_coil_heating_retains_distinct_cavity_and_exterior_media() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (card, curve) = material(0.005, basis, 2000.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = initial(&curve, 300.0, 0.0);
    let exterior = crate::speaker::declared_air();
    let cavity = GasState::try_new(
        &fs_material::gas_species::GasSpecies::Helium.spec(),
        310.0,
        exterior.pressure,
    )
    .unwrap();
    let mut coupled = ThermalMagneticTsDriver::from_material_with_media(
        source_mechanics(),
        Some(0.02),
        None,
        exterior,
        Some(cavity),
        source_profile(),
        &thermal,
        &network,
        basis,
        || false,
    )
    .unwrap();
    let original_mechanics = coupled.driver().storage.mechanical.clone();
    let initial_resistance = coupled.driver().coil_resistance_ohm();
    let mut emitted = false;
    for tick in 1..=5 {
        let proposal = coupled
            .step_candidate(
                &thermal,
                f64::from(tick) * 0.001,
                0.2,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let accepted = proposal.driver().driver();
        assert_eq!(accepted.exterior_gas(), &exterior);
        assert_eq!(accepted.cavity_gas(), Some(&cavity));
        assert_eq!(accepted.storage.mechanical, original_mechanics);
        assert_eq!(accepted.time_s(), proposal.thermal().time_s());
        assert_eq!(
            proposal.report().electrical.coil.material().card_identity(),
            card.content_hash()
        );
        near(proposal.report().total_energy_residual_j, 0.0, 2e-9);
        emitted |= proposal.report().electrical.pressure_pa.abs() > 1e-9;
        let (next_thermal, next_driver, _) = proposal.into_parts();
        thermal = next_thermal;
        coupled = next_driver;
    }
    assert!(emitted);
    assert!(thermal.temperature_k(COIL).unwrap() > 300.0);
    assert!(coupled.driver().coil_resistance_ohm() > initial_resistance);
}

#[test]
fn source_coil_mass_and_resistance_swaps_reach_native_motion_and_pressure() {
    let mut first_steps = Vec::new();
    for (density, resistance_scale, basis) in [
        (1.0, 1.0, ElectricalPropertyBasis::Resistivity),
        (4.0, 1.0, ElectricalPropertyBasis::Resistivity),
        (1.0, 2.0, ElectricalPropertyBasis::Resistivity),
        (1.0, 1.0, ElectricalPropertyBasis::Conductivity),
    ] {
        let (card, curve) = material_with_properties(0.0, basis, 2000.0, density, resistance_scale);
        let context = point(300.0);
        let (cells, resistors) = bindings(&card, &curve, &context);
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut thermal = initial(&curve, 300.0, 0.0);
        let mechanics = source_mechanics();
        let mut coupled = ThermalMagneticTsDriver::from_material(
            mechanics,
            None,
            None,
            crate::speaker::declared_air(),
            source_profile(),
            &thermal,
            &network,
            basis,
            || false,
        )
        .unwrap();
        // Independent geometry/source inventory: same length*area in every run.
        let coil_mass = density * 0.1 * 0.01;
        let total_mass = mechanics.non_coil_moving_mass_kg + coil_mass;
        let resistance = 0.08 * resistance_scale * 0.1 / 0.01;
        near(coupled.coil_mass_kg(), coil_mass, 1e-17);
        near(coupled.moving_assembly_mass_kg(), total_mass, 1e-17);
        near(coupled.coil().resistance_ohm(), resistance, 2e-15);
        near(coupled.driver().card.bl, 1.2, 0.0);
        near(coupled.driver().card.le_h, 0.006, 0.0);
        assert_eq!(
            coupled.coil().material().card_identity(),
            card.content_hash()
        );

        let mut native = MagneticTsDriver::new_with_gas(
            TsCard {
                re_ohm: resistance,
                le_h: 0.006,
                bl: 1.2,
                mms_kg: total_mass,
                cms_m_per_n: 0.9,
                rms_n_s_m: 0.1,
                sd_m2: 0.001,
                x_max_m: 0.1,
            },
            None,
            None,
            crate::speaker::declared_air(),
            source_profile(),
        )
        .unwrap();
        for tick in 1..=5 {
            let (reference, pressure) = native.step(0.05, 0.001).unwrap();
            let proposal = coupled
                .step_candidate(
                    &thermal,
                    f64::from(tick) * 0.001,
                    0.05,
                    controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            for (&actual, &expected) in proposal
                .report()
                .electrical
                .native
                .x
                .iter()
                .zip(&reference.x)
            {
                near(actual, expected, 3e-13);
            }
            near(proposal.report().electrical.pressure_pa, pressure, 1e-12);
            near(proposal.report().total_energy_residual_j, 0.0, 2e-9);
            if tick == 1 {
                first_steps.push([
                    proposal.driver().driver().coil_current_a(),
                    proposal.driver().driver().excursion_m(),
                    proposal.report().electrical.pressure_pa,
                ]);
            }
            let (next_thermal, next_driver, _) = proposal.into_parts();
            thermal = next_thermal;
            coupled = next_driver;
            near(coupled.coil_mass_kg(), coil_mass, 1e-17);
            if tick == 2 {
                // Rebinding accepted nonzero motion must preserve linkage,
                // momentum, displacement and clock, not recreate a cold driver.
                let state = coupled.driver().x.clone();
                coupled = ThermalMagneticTsDriver::bind(
                    coupled.driver,
                    &thermal,
                    &network,
                    basis,
                    || false,
                )
                .unwrap();
                assert_eq!(coupled.driver().x, state);
                assert_eq!(coupled.driver().time_s(), thermal.time_s());
            }
        }
    }
    let [base, dense, resistive, conductivity] = first_steps.as_slice() else {
        panic!("all four source substitutions must run");
    };
    assert!((dense[0] - base[0]).abs() > 1e-8); // Reciprocal back EMF also changes.
    assert!(dense[1] > 0.0 && dense[1] < 0.75 * base[1]);
    assert!(dense[2] > 0.0 && dense[2] < 0.75 * base[2]);
    for axis in 0..3 {
        assert!(resistive[axis] > 0.0 && resistive[axis] < base[axis]);
        near(conductivity[axis], base[axis], 1e-12);
    }
}

#[test]
fn material_construction_preserves_thermal_clock_and_refuses_stale_coil_mass() {
    let (card, curve) = material(0.0, ElectricalPropertyBasis::Resistivity, 2000.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&card, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 300.0, 0.0);
    let make = |mechanics, thermal: &DepositionSimulation<'_>, cancel| {
        ThermalMagneticTsDriver::from_material(
            mechanics,
            None,
            None,
            crate::speaker::declared_air(),
            source_profile(),
            thermal,
            &network,
            ElectricalPropertyBasis::Resistivity,
            move || cancel,
        )
    };
    let coupled = make(source_mechanics(), &thermal, false).unwrap();
    let saved = thermal.cells().clone();
    for mass in [-0.1, f64::NAN, f64::INFINITY] {
        let mut mechanics = source_mechanics();
        mechanics.non_coil_moving_mass_kg = mass;
        assert!(make(mechanics, &thermal, false).is_err());
    }
    assert!(matches!(
        make(source_mechanics(), &thermal, true),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let changed_mass = initial_with_mass(&curve, 300.0, 0.0, 2.0 * MASS);
    assert!(make(source_mechanics(), &changed_mass, false).is_err());
    assert!(matches!(
        coupled.step_candidate(&changed_mass, 0.001, 0.05, controls(), &[0.0], || false,),
        Err(CircuitEnthalpyError::Invalid(_))
    ));
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(coupled.driver().stored_energy_j(), 0.0);

    let advanced = thermal
        .advance_to(
            0.01,
            StepConfig {
                max_dt_s: 0.01,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let fresh = make(source_mechanics(), &advanced, false).unwrap();
    assert_eq!(fresh.driver().time_s(), 0.01);
    assert_eq!(fresh.driver().stored_energy_j(), 0.0);
    assert_eq!(fresh.driver().coil_current_a(), 0.0);
    let next = fresh
        .step_candidate(&advanced, 0.011, 0.05, controls(), &[0.0], || false)
        .unwrap();
    assert_eq!(next.driver().driver().time_s(), 0.011);
    assert!(next.driver().driver().excursion_m() > 0.0);
}
