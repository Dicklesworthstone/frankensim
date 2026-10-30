//! Independent biased-coil energy, thermal-feedback and trajectory fixtures.
//! All material values below are authored analytic tests, not measured grades.
use super::*;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    material::{CellResistor, MaterialCell, MaterialNetwork},
};
use crate::magnetostrictive::{MagnetostrictiveLimits, MagnetostrictiveRod};
use crate::material_circuit::MaterialCircuitSolve;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::MAGNETIC_PERMEABILITY_DIMS,
    magnetostrictive::*,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

mod heat;

const COIL: CellId = CellId(41);
const MASS: f64 = 0.001;
const CAPACITY: f64 = 0.01;
const RESISTANCE: f64 = 0.8;
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn card(claims: ClaimSet, chemistry: &str) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "authored analytic regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn rod(d: f64, bias: f64) -> MagnetostrictiveRod {
    let roles = [
        (COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS, d),
        (PERMEABILITY_33_T_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 2.0),
        ("density", Density::DIMS, 3.0),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in roles {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic biased rod".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let material = card(claims, "synthetic rod");
    let context = point(300.0)
        .with_quantity(
            "bias_H",
            QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
            bias,
        )
        .unwrap();
    let law = LongitudinalMagnetostrictive::resolve(
        &material,
        &context,
        "bias_H",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MagnetostrictiveRod::new(
        law,
        0.2,
        0.005,
        3,
        MagnetostrictiveLimits {
            maximum_abs_strain: 0.05,
            maximum_abs_field_increment_a_m: 3.0,
        },
    )
    .unwrap()
}
fn winding(
    alpha: f64,
    basis: ElectricalPropertyBasis,
    upper: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new(basis.property_name(), basis.dimensions());
    let base = if basis == ElectricalPropertyBasis::Resistivity {
        0.08
    } else {
        12.5
    };
    let mut claims = ClaimSet::new();
    for (key, low, high, limit) in [
        (heat.clone(), 10.0, 10.0, 600.0),
        (density.clone(), 1.0, 1.0, 600.0),
        (
            electrical,
            base * (1.0 - 10.0 * alpha),
            base * (1.0 + 300.0 * alpha),
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
                    knots: vec![(290.0, low), (600.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    290.0,
                    limit,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic thermal winding".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let material = card(claims, "synthetic winding");
    let curve = ResolvedHeatCapacityCurve::try_new(
        &material,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(600.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (material, curve)
}
fn initial(curve: &ResolvedHeatCapacityCurve, conductance: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: COIL,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: MASS,
                specific_h_j_kg: 100.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: 0.001,
                    surface_area_m2: 0.06,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: conductance / 0.06,
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
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 200,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-6,
        cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 2e-10,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-10,
            power_tolerance_w: 2e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            ..NewtonKrylovConfig::default()
        },
    }
}

#[test]
fn dc_bias_alone_heats_the_finite_winding_even_when_incremental_loss_is_zero() {
    let (material, curve) = winding(0.0, ElectricalPropertyBasis::Resistivity, 600.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&material, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 0.0);
    let rod = rod(0.2, 10.0);
    let bias = rod.bias_current_a();
    let voltage = RESISTANCE * bias;
    let driver = ThermalMagnetostrictiveCircuit::new(
        rod,
        0.0,
        0.0,
        0.0,
        voltage,
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let dt = 0.02;
    let candidate = driver
        .step_candidate(&thermal, dt, voltage, 0.0, controls(), &[0.0], || false)
        .unwrap();
    let s = &candidate.report().electrical;
    let heat = dt * RESISTANCE * bias * bias;
    near(s.incremental.incremental_resistor_dissipation_j, 0.0, 1e-24);
    near(s.coil_heat_j, heat, 1e-12);
    near(s.electrical_supply_j, heat, 1e-12);
    near(s.bias_storage_change_j, 0.0, 1e-14);
    near(
        candidate.thermal().temperature_k(COIL).unwrap(),
        300.0 + heat / CAPACITY,
        1e-8,
    );
    near(candidate.driver().current_a().unwrap(), bias, 1e-12);
    near(candidate.report().total_energy_residual_j, 0.0, 2e-10);
    near(thermal.temperature_k(COIL).unwrap(), 300.0, 1e-12);
    near(driver.time_s(), 0.0, 0.0);
}

#[test]
fn total_terminal_work_closes_with_bias_linkage_storage_and_reciprocal_motion() {
    let (material, curve) = winding(0.0, ElectricalPropertyBasis::Resistivity, 600.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&material, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 0.0);
    let rod = rod(0.2, 10.0);
    let (m, k, l, theta, bias) = (
        rod.modal_mass_kg(),
        rod.fixed_current_stiffness_n_m(),
        rod.clamped_inductance_h(),
        rod.coupling_wb_m(),
        rod.bias_current_a(),
    );
    let (x0, v0, i0, dt, force) = (0.0001, -0.003, 0.01, 0.01, 0.02);
    let lambda0 = l * i0 + theta * x0;
    let voltage = -0.05;
    let a = 2.0 * l / dt + RESISTANCE;
    let electrical = voltage - RESISTANCE * bias + 2.0 * l * i0 / dt;
    let vm = (force - k * x0 + 2.0 * m * v0 / dt + theta * electrical / a)
        / (2.0 * m / dt + k * dt / 2.0 + theta * theta / a);
    let im = (electrical - theta * vm) / a;
    let (x1, v1, i1) = (x0 + dt * vm, 2.0 * vm - v0, 2.0 * im - i0);
    let lambda1 = l * i1 + theta * x1;
    let expected_storage = 0.5 * m * (v1 * v1 - v0 * v0)
        + 0.5 * k * (x1 * x1 - x0 * x0)
        + 0.5 * l * (i1 * i1 - i0 * i0)
        + bias * (lambda1 - lambda0);
    let driver = ThermalMagnetostrictiveCircuit::new(
        rod,
        lambda0,
        x0,
        v0,
        voltage,
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let candidate = driver
        .step_candidate(&thermal, dt, voltage, force, controls(), &[0.0], || false)
        .unwrap();
    let s = &candidate.report().electrical;
    near(candidate.driver().circuit().extension_m(), x1, 2e-12);
    near(candidate.driver().circuit().velocity_m_s(), v1, 2e-11);
    near(candidate.driver().current_a().unwrap(), bias + i1, 2e-11);
    near(s.midpoint_current_a, bias + im, 2e-11);
    near(s.stored_change_j, expected_storage, 1e-11);
    near(s.bias_storage_change_j, bias * (lambda1 - lambda0), 1e-11);
    near(s.electrical_supply_j, dt * voltage * (bias + im), 1e-11);
    near(s.coil_heat_j, dt * RESISTANCE * (bias + im).powi(2), 1e-11);
    near(s.mechanical_supply_j, force * (x1 - x0), 1e-12);
    assert!(s.coil_heat_j > 100.0 * s.incremental.incremental_resistor_dissipation_j);
    assert!(s.bias_storage_change_j.abs() > 1e-5);
    assert!(s.electrical_supply_j < 0.0);
    assert!(s.stored_change_j < 0.0);
    near(s.energy_residual_j, 0.0, 1e-10);
    near(candidate.report().total_energy_residual_j, 0.0, 2e-10);
}

#[test]
fn temperature_dependent_sources_close_the_same_step_without_coefficient_lag() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let alpha = 0.02;
        let (material, curve) = winding(alpha, basis, 600.0);
        let context = point(300.0);
        let (cells, resistors) = bindings(&material, &curve, &context);
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let thermal = initial(&curve, 0.012);
        let rod = rod(0.0, 10.0);
        let (l, bias) = (rod.clamped_inductance_h(), rod.bias_current_a());
        let voltage = RESISTANCE * bias;
        let dt = 0.2;
        let r = |t: f64| {
            if basis == ElectricalPropertyBasis::Resistivity {
                RESISTANCE * (1.0 + alpha * (t - 300.0))
            } else {
                RESISTANCE / (1.0 + alpha * (t - 300.0))
            }
        };
        let im = |t: f64| (2.0 * l * bias + dt * voltage) / (2.0 * l + dt * r(t));
        let residual =
            |t: f64| CAPACITY * (t - 300.0) + dt * 0.012 * (t - 300.0) - dt * r(t) * im(t).powi(2);
        let (mut lo, mut hi) = (300.0, 400.0);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..60 {
            let mid = f64::midpoint(lo, hi);
            if residual(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let expected = f64::midpoint(lo, hi);
        let driver = ThermalMagnetostrictiveCircuit::new(
            rod,
            0.0,
            0.0,
            0.0,
            voltage,
            &thermal,
            &network,
            basis,
            || false,
        )
        .unwrap();
        let candidate = driver
            .step_candidate(&thermal, dt, voltage, 0.0, controls(), &[0.0], || false)
            .unwrap();
        near(
            candidate.thermal().temperature_k(COIL).unwrap(),
            expected,
            2e-7,
        );
        near(
            candidate.driver().coil().resistance_ohm(),
            r(expected),
            2e-9,
        );
        near(
            candidate.driver().current_a().unwrap(),
            2.0 * im(expected) - bias,
            2e-9,
        );
        if basis == ElectricalPropertyBasis::Resistivity {
            assert!(candidate.driver().current_a().unwrap() < bias);
        } else {
            assert!(candidate.driver().current_a().unwrap() > bias);
        }
        near(candidate.report().total_energy_residual_j, 0.0, 2e-10);
    }
}

fn reference(
    rod: &MagnetostrictiveRod,
    alpha: f64,
    voltage: f64,
    conductance: f64,
    end: f64,
) -> [f64; 4] {
    let (m, k, l, theta, bias) = (
        rod.modal_mass_kg(),
        rod.fixed_current_stiffness_n_m(),
        rod.clamped_inductance_h(),
        rod.coupling_wb_m(),
        rod.bias_current_a(),
    );
    let derivative = |s: [f64; 4]| {
        let r = RESISTANCE * (1.0 + alpha * (s[3] - 300.0));
        let current = bias + s[2];
        [
            s[1],
            (-k * s[0] + theta * s[2]) / m,
            (voltage - r * current - theta * s[1]) / l,
            (r * current * current - conductance * (s[3] - 300.0)) / CAPACITY,
        ]
    };
    // Independent continuous ODE RK4 reference, confined to this test.
    let n = 16_384;
    let dt = end / n as f64;
    let mut s = [0.0, 0.0, 0.0, 300.0];
    for _ in 0..n {
        let a = derivative(s);
        let b = derivative(std::array::from_fn(|j| s[j] + 0.5 * dt * a[j]));
        let c = derivative(std::array::from_fn(|j| s[j] + 0.5 * dt * b[j]));
        let d = derivative(std::array::from_fn(|j| s[j] + dt * c[j]));
        for j in 0..4 {
            s[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    s
}
#[test]
fn coupled_motion_and_thermal_feedback_converge_to_the_independent_continuous_equations() {
    let alpha = 0.02;
    let (material, curve) = winding(alpha, ElectricalPropertyBasis::Resistivity, 600.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&material, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let prototype = rod(0.2, 10.0);
    let voltage = RESISTANCE * prototype.bias_current_a() + 0.015;
    let expected = reference(&prototype, alpha, voltage, 0.012, 0.2);
    let mut errors = Vec::new();
    for n in [40, 80, 160] {
        let mut thermal = initial(&curve, 0.012);
        let mut driver = ThermalMagnetostrictiveCircuit::new(
            prototype.clone(),
            0.0,
            0.0,
            0.0,
            voltage,
            &thermal,
            &network,
            ElectricalPropertyBasis::Resistivity,
            || false,
        )
        .unwrap();
        for k in 1..=n {
            let candidate = driver
                .step_candidate(
                    &thermal,
                    0.2 * k as f64 / n as f64,
                    voltage,
                    0.0,
                    controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let (next_thermal, next_driver, report) = candidate.into_parts();
            near(report.total_energy_residual_j, 0.0, 2e-10);
            thermal = next_thermal;
            driver = next_driver;
        }
        let observed = [
            driver.circuit().extension_m(),
            driver.circuit().velocity_m_s(),
            driver.current_a().unwrap() - prototype.bias_current_a(),
            thermal.temperature_k(COIL).unwrap(),
        ];
        let error = (observed[3] - expected[3]).abs()
            + 10.0 * (observed[2] - expected[2]).abs()
            + 100.0 * (observed[0] - expected[0]).abs()
            + (observed[1] - expected[1]).abs();
        errors.push(error);
    }
    assert!(expected[0].abs() > 1e-6);
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (0.8..1.2).contains(&order),
            "errors {errors:?}, order {order}"
        );
    }
}

#[test]
fn simultaneous_bias_and_winding_orientation_reversal_preserves_motion_and_heat() {
    let (material, curve) = winding(0.02, ElectricalPropertyBasis::Resistivity, 600.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&material, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 0.0);
    let mut results = Vec::new();
    for sign in [1.0, -1.0] {
        let rod = rod(sign * 0.2, sign * 10.0);
        let voltage = RESISTANCE * rod.bias_current_a() + sign * 0.05;
        let driver = ThermalMagnetostrictiveCircuit::new(
            rod,
            0.0,
            0.0,
            0.0,
            voltage,
            &thermal,
            &network,
            ElectricalPropertyBasis::Resistivity,
            || false,
        )
        .unwrap();
        let c = driver
            .step_candidate(&thermal, 0.04, voltage, 0.0, controls(), &[0.0], || false)
            .unwrap();
        results.push([
            c.driver().circuit().extension_m(),
            c.driver().circuit().velocity_m_s(),
            sign * c.driver().current_a().unwrap(),
            c.thermal().temperature_k(COIL).unwrap(),
            c.report().electrical.coil_heat_j,
            c.report().electrical.electrical_supply_j,
        ]);
    }
    for j in 0..6 {
        near(results[0][j], results[1][j], 2e-10);
    }
}

#[test]
fn source_limits_real_budgets_and_late_cancellation_preserve_both_accepted_histories() {
    let (material, curve) = winding(0.02, ElectricalPropertyBasis::Resistivity, 600.0);
    let context = point(300.0);
    let (cells, resistors) = bindings(&material, &curve, &context);
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curve, 0.0);
    let rod = rod(0.2, 10.0);
    let voltage = RESISTANCE * rod.bias_current_a() + 0.02;
    let driver = ThermalMagnetostrictiveCircuit::new(
        rod,
        0.0,
        0.0,
        0.0,
        voltage,
        &thermal,
        &network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    let mut calls = 0;
    let accepted = driver
        .step_candidate(&thermal, 0.02, voltage, 0.0, controls(), &[0.0], || {
            calls += 1;
            false
        })
        .unwrap();
    let mut checked = 0;
    assert!(matches!(
        driver.step_candidate(&thermal, 0.02, voltage, 0.0, controls(), &[0.0], || {
            checked += 1;
            checked >= calls
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut bounded = controls();
    bounded.maximum_circuit_trials = 1;
    assert!(matches!(
        driver.step_candidate(&thermal, 0.02, voltage, 0.0, bounded, &[0.0], || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    near(driver.time_s(), 0.0, 0.0);
    near(driver.circuit().flux_linkage_increment_wb(), 0.0, 0.0);
    near(thermal.temperature_k(COIL).unwrap(), 300.0, 1e-12);
    let repeated = driver
        .step_candidate(&thermal, 0.02, voltage, 0.0, controls(), &[0.0], || false)
        .unwrap();
    assert_eq!(accepted.driver().circuit(), repeated.driver().circuit());
    assert_eq!(accepted.thermal().cells(), repeated.thermal().cells());

    let (limited_material, limited_curve) =
        winding(0.02, ElectricalPropertyBasis::Resistivity, 300.01);
    let (limited_cells, limited_resistors) = bindings(&limited_material, &limited_curve, &context);
    let limited_network = MaterialNetwork {
        cells: &limited_cells,
        resistors: &limited_resistors,
        links: &[],
    };
    let limited_thermal = initial(&limited_curve, 0.0);
    let limited_driver = ThermalMagnetostrictiveCircuit::new(
        self::rod(0.2, 10.0),
        0.0,
        0.0,
        0.0,
        voltage,
        &limited_thermal,
        &limited_network,
        ElectricalPropertyBasis::Resistivity,
        || false,
    )
    .unwrap();
    assert!(matches!(
        limited_driver.step_candidate(
            &limited_thermal,
            0.1,
            voltage,
            0.0,
            controls(),
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::Source(
            MagnetostrictiveThermalError::Material(_)
        ))
    ));
    near(limited_thermal.temperature_k(COIL).unwrap(), 300.0, 1e-12);
}
