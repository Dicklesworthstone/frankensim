//! Source-backed gas composition: native calorics, both port tangents,
//! analytic thermal paths and a piston assembled without pairing code.
#![cfg(feature = "component")]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_gas::component::{ClosedGasParcel, GasCavity, GasEnergyError, GasEnergySolve};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortTangent, Proposal};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_system::elements::{Force, Mass, Motion, Restriction, SourceAccount, Waveform};
use fs_system::{AdmissionError, NodeId, SolveConfig, System, SystemBuilder, SystemError};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

const T0: f64 = 400.0;
const P0: f64 = 100_000.0;
const V0: f64 = 0.001;
const MOLAR_MASS: f64 = 0.028;

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1.0e-7,
        maximum_iterations: 128,
    }
}

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey { seed: 113, kernel_id: 7, tile: 0, iteration: 0 };
        f(&Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn provenance() -> Provenance {
    Provenance {
        source: "authored gas Component oracle; no measured-data claim".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

/// One authored species: Cp/R=3.5+slope*T. Formation datum is independently
/// variable; it contributes no work, heat capacity or pressure change.
fn parcel(slope: f64, formation_j_mol: f64) -> ClosedGasParcel {
    let id = SpeciesId::new("fixture-gas").unwrap();
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[3] = slope;
    coefficients[7] = formation_j_mol / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: P0, dims: Pressure::DIMS });
    let species = Nasa9StandardStateModelV1::new(
        id.clone(),
        MolarMass::new(MOLAR_MASS),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("gas-component-fixture").unwrap(),
        ).unwrap(),
        vec![Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            sources: vec![],
            provenance: provenance(),
        }).unwrap()],
    ).unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2.0e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), 200.0, 2000.0)
                .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1.0e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: provenance(),
        }).unwrap();
    }
    let transport = SpeciesTransportMaterial {
        species: id,
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "fixture-gas".into(), phase: "gas".into(),
                process: "authored".into(), revision: 0,
            },
            claims,
            vec![],
        ).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        ).unwrap(),
        vec![transport],
        None,
    ).unwrap();
    let point = QueryPoint::new()
        .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), T0).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), P0).unwrap()
        .with("specimen-lot", 42.0).unwrap();
    ClosedGasParcel::new(model, point, V0).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.16e} vs {expected:.16e}, allowance {tolerance:.3e}",
    );
}

fn step(
    cavity: &GasCavity,
    old: &ClosedGasParcel,
    dt: f64,
    inputs: [f64; 2],
    cx: &Cx<'_>,
) -> (ClosedGasParcel, [f64; 2]) {
    let mut next = old.clone();
    let mut output = [f64::NAN; 2];
    assert!(matches!(
        cavity.propose(old, dt, &inputs, &mut next, &mut output, cx).unwrap(),
        Proposal::Step(_)
    ));
    (next, output)
}

#[test]
fn native_source_steps_conform_and_both_inputs_have_analytic_tangents() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for slope in [0.0, 0.002] {
            let old = parcel(slope, 0.0);
            let stale = parcel(0.0, 12.0);
            let cavity = GasCavity::new(old.clone(), controls()).unwrap();
            assert_eq!(cavity.class(), Class::Smooth);
            for inputs in [[0.0, 0.0], [0.05, 0.0], [-0.03, 2.0], [0.0, 10.0]] {
                let case = Case { state: &old, stale: &stale, dt_s: 0.001, inputs: &inputs };
                check_component(&cavity, &case, cx).unwrap();
                check_port_tangent(&cavity, &case, 2.0e-5, 1.0e-5, cx).unwrap();
                let (next, outputs) = step(&cavity, &old, case.dt_s, inputs, cx);
                let t = next.state().state().temperature;
                // Independent integrals for the explicitly authored Cv/R law.
                let scale = old.mass_kg() * R / MOLAR_MASS;
                let du = scale * (t - T0) * (2.5 + 0.5 * slope * (t + T0));
                let ds = scale * (2.5 * (t / T0).ln() + slope * (t - T0)
                    + (next.volume_m3() / V0).ln());
                close(cavity.stored_energy_j(&next).unwrap(), du, 5.0e-10);
                close(ds, case.dt_s * inputs[1] / outputs[1], 5.0e-12);
            }
            // Check each column separately against independently perturbed
            // pressure and temperature, so a small heat derivative cannot
            // hide beside the much larger compression-flow derivative.
            let inputs = [0.02, 5.0];
            for (column, perturbation) in [(0, 1.0e-5), (1, 0.1)] {
                let mut plus = inputs;
                let mut minus = inputs;
                plus[column] += perturbation;
                minus[column] -= perturbation;
                let (_, a) = step(&cavity, &old, 0.002, plus, cx);
                let (_, b) = step(&cavity, &old, 0.002, minus, cx);
                let mut direction = [0.0; 2];
                direction[column] = 1.0;
                let mut tangent = [0.0; 2];
                cavity.jvp(&old, 0.002, &inputs, &direction, &mut tangent, cx).unwrap();
                for row in 0..2 {
                    let fd = (a[row] - b[row]) / (2.0 * perturbation);
                    close(tangent[row], fd, 2.0e-5 * fd.abs() + 1.0e-5);
                }
            }
        }
        // Analytic coincidence limit, including both off-diagonal entries.
        let old = parcel(0.0, 0.0);
        let cavity = GasCavity::new(old.clone(), controls()).unwrap();
        let dt = 0.001;
        let cv = 2.5 * R / MOLAR_MASS;
        let expected = [
            [dt * 1.4 * P0 / (2.0 * V0), dt * R / MOLAR_MASS / (2.0 * V0 * cv)],
            [dt * 0.4 * T0 / (2.0 * V0), dt / (2.0 * old.mass_kg() * cv)],
        ];
        for column in 0..2 {
            let mut direction = [0.0; 2];
            direction[column] = 1.0;
            let mut output = [0.0; 2];
            cavity.jvp(&old, dt, &[0.0; 2], &direction, &mut output, cx).unwrap();
            for row in 0..2 {
                close(output[row], expected[row][column], 1.0e-12 * expected[row][column].abs());
            }
        }
    });
}

#[test]
fn constant_capacity_adiabatic_and_fixed_volume_heat_follow_closed_forms() {
    let old = parcel(0.0, 0.0);
    let original = old.clone();
    let cavity = GasCavity::new(old.clone(), controls()).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for ratio in [0.8, 1.0, 1.25] {
            let dt = 0.1;
            let volume = V0 * ratio;
            let (next, _) = step(&cavity, &old, dt, [(V0 - volume) / dt, 0.0], cx);
            close(next.state().state().temperature, T0 * ratio.powf(-0.4), 2.0e-8);
            close(next.state().state().pressure, P0 * ratio.powf(-1.4), 1.0e-5);
            assert_eq!(next.mass_kg().to_bits(), old.mass_kg().to_bits());
            assert_eq!(next.point().axes().get("specimen-lot"), Some(&42.0));
        }
        let capacity = old.mass_kg() * 2.5 * R / MOLAR_MASS;
        for heat in [-20.0, 0.0, 30.0] {
            let (next, output) = step(&cavity, &old, 0.1, [0.0, heat / 0.1], cx);
            let temperature = T0 + heat / capacity;
            let theta = if heat == 0.0 { T0 } else { (temperature - T0) / (temperature / T0).ln() };
            close(next.state().state().temperature, temperature, 2.0e-8);
            close(output[1], theta, 2.0e-8);
            close(output[0], P0 * theta / T0, 1.0e-5);
            close(cavity.stored_energy_j(&next).unwrap(), heat, 5.0e-10);
            close(heat / output[1], capacity * (temperature / T0).ln(), 5.0e-12);
            assert_eq!(next.volume_m3().to_bits(), V0.to_bits());
        }
    });
    assert_eq!(old, original);
}

#[test]
fn energy_measurement_and_dynamics_ignore_large_formation_reference_shifts() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut reference: Option<[f64; 4]> = None;
        for formation in [0.0, 1.0e12, -1.0e12] {
            let old = parcel(0.002, formation);
            let cavity = GasCavity::new(old.clone(), controls()).unwrap();
            assert_eq!(cavity.stored_energy_j(&old).unwrap().to_bits(), 0.0_f64.to_bits());
            let (next, output) = step(&cavity, &old, 0.01, [0.005, 100.0], cx);
            let values = [
                next.state().state().temperature, output[0], output[1],
                cavity.stored_energy_j(&next).unwrap(),
            ];
            if let Some(expected) = reference {
                for (a, b) in values.into_iter().zip(expected) {
                    close(a, b, 1.0e-10 * b.abs().max(1.0));
                }
            } else {
                reference = Some(values);
            }
        }
    });
}

/// One graph: gas, inertia, area transformer and an insulated heat boundary.
/// No custom gas/mechanics integration is used by the system under test.
fn piston(dt: f64) -> [f64; 2] {
    let gas = parcel(0.0, 0.0);
    let mass_kg = 1.0;
    let area = 0.01;
    let mut builder = SystemBuilder::new();
    let body = builder.add_with_tangent(
        "piston", Mass::new(mass_kg).unwrap(),
        Motion { position_m: V0 / area, velocity_m_s: 0.0 },
    ).unwrap();
    let cavity = builder.add_with_tangent(
        "gas", GasCavity::new(gas.clone(), controls()).unwrap(), gas.clone(),
    ).unwrap();
    builder.transformer(
        builder.port(body, "p").unwrap(), builder.port(cavity, "p").unwrap(), area,
    ).unwrap();
    builder.boundary(
        builder.port(cavity, "heat").unwrap(), Box::new(|_, _, values| values[0] = 0.0),
    ).unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig { relative_tolerance: 1.0e-11, ..SolveConfig::default() });
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..(0.008 / dt).round() as usize {
            system.step(dt, cx).unwrap();
        }
    });
    let motion = system.state::<Motion>(body).unwrap();
    let final_gas = system.state::<ClosedGasParcel>(cavity).unwrap();
    let volume = final_gas.volume_m3();
    let temperature = final_gas.state().state().temperature;
    close(volume, area * motion.position_m, 5.0e-15);
    close(temperature, T0 * (V0 / volume).powf(0.4), 3.0e-7);
    let du = gas.mass_kg() * 2.5 * R / MOLAR_MASS * (temperature - T0);
    close(du + 0.5 * mass_kg * motion.velocity_m_s.powi(2), 0.0, 1.0e-7);
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-7);
    close(system.time_s(), 0.008, 1.0e-15);
    [volume, motion.velocity_m_s]
}

#[test]
fn generic_vacuum_piston_refines_to_the_independent_adiabatic_ode() {
    // dV/dt=A*v, dv/dt=A*p0*(V0/V)^gamma/M. RK4 uses only this
    // independent ideal-gas equation, never a material proposal or tangent.
    let mut reference = [V0, 0.0];
    let derivative = |state: [f64; 2]| [
        0.01 * state[1],
        0.01 * P0 * (V0 / state[0]).powf(1.4),
    ];
    let h = 1.0e-6;
    for _ in 0..8000 {
        let a = derivative(reference);
        let b = derivative([reference[0] + 0.5*h*a[0], reference[1] + 0.5*h*a[1]]);
        let c = derivative([reference[0] + 0.5*h*b[0], reference[1] + 0.5*h*b[1]]);
        let d = derivative([reference[0] + h*c[0], reference[1] + h*c[1]]);
        for i in 0..2 {
            reference[i] += h * (a[i] + 2.0*b[i] + 2.0*c[i] + d[i]) / 6.0;
        }
    }
    let coarse = piston(0.0001);
    let fine = piston(0.00005);
    let norm = |value: [f64; 2]| {
        ((value[0] - reference[0]) / V0).abs()
            .max(((value[1] - reference[1]) / reference[1]).abs())
    };
    assert!(norm(fine) < 0.3 * norm(coarse), "coarse {}, fine {}", norm(coarse), norm(fine));
    assert!(norm(fine) < 2.0e-5, "fine error {}", norm(fine));
}

const PISTON_KG: f64 = 1.0;
const PISTON_M2: f64 = 0.01;

/// A single-acting piston about its rest volume, from a graph: the gas on
/// one face of a mass, a constant force `P0 * area` on the other, no heat.
/// It starts at rest volume moving outward at `rate_m_s`. With `insulated`
/// the gas is the one-port cavity, a scalar monotone relation; without, it
/// is the two-port cavity with its heat port held at zero. With `tangent`
/// the gas is added with the tangent it offers; without, the kernel is left
/// to difference it. Returns the system, the mass and the gas.
fn spring_piston(insulated: bool, tangent: bool, rate_m_s: f64) -> (System, NodeId, NodeId) {
    let gas = parcel(0.0, 0.0);
    let mut builder = SystemBuilder::new();
    let body = builder.add_with_tangent(
        "piston", Mass::new(PISTON_KG).unwrap(),
        Motion { position_m: V0 / PISTON_M2, velocity_m_s: rate_m_s },
    ).unwrap();
    let closed = GasCavity::new(gas.clone(), controls()).unwrap();
    let cavity = match (insulated, tangent) {
        (true, true) => builder.add_with_tangent("gas", closed.insulated().unwrap(), gas),
        (true, false) => builder.add("gas", closed.insulated().unwrap(), gas),
        (false, true) => builder.add_with_tangent("gas", closed, gas),
        (false, false) => builder.add("gas", closed, gas),
    }
    .unwrap();
    if !insulated {
        builder.boundary(
            builder.port(cavity, "heat").unwrap(), Box::new(|_, _, values| values[0] = 0.0),
        ).unwrap();
    }
    let outside = builder.add_with_tangent(
        "outside", Force::new(Waveform::constant(-P0 * PISTON_M2).unwrap()), SourceAccount::default(),
    ).unwrap();
    // The mass sets the velocity; the gas sees it through the area, the
    // outside force directly.
    builder.common_flow_ratios(
        builder.port(body, "p").unwrap(),
        &[(builder.port(cavity, "p").unwrap(), PISTON_M2), (builder.port(outside, "p").unwrap(), 1.0)],
    ).unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig { relative_tolerance: 1.0e-11, ..SolveConfig::default() });
    (system, body, cavity)
}

/// The angular frequency [rad/s] of [`spring_piston`] on the insulated
/// cavity over `periods` periods, from the times at which the mass's
/// velocity passes through zero going down, each found by a straight line
/// between the two steps around it.
fn piston_frequency(dt: f64, rate_m_s: f64, periods: usize) -> f64 {
    let gas_kg = parcel(0.0, 0.0).mass_kg();
    let (mut system, body, cavity) = spring_piston(true, true, rate_m_s);
    let mut crossings = Vec::new();
    let mut before = rate_m_s;
    let mut steps = 0_u32;
    with_cx(&CancelGate::new_clock_free(), |cx| {
        while crossings.len() <= periods {
            system.step(dt, cx).unwrap();
            let after = system.state::<Motion>(body).unwrap().velocity_m_s;
            if before > 0.0 && after <= 0.0 {
                crossings.push(dt * (f64::from(steps) + before / (before - after)));
            }
            before = after;
            steps += 1;
            assert!(steps < 100_000, "no oscillation");
        }
    });
    let omega = core::f64::consts::TAU * periods as f64 / (crossings[periods] - crossings[0]);
    // Two accounts. What the gas owner's energy solve leaves in a step the
    // cavity books as numerical loss: at most its tolerance per kilogram
    // and the rounding of the two p V the step's volume is formed from.
    // That booked loss is what a cavity returning a wrong pressure would
    // have to hide in, so it is held to that many steps' worth. The
    // ledger's residual is what nobody booked; it cancels a booked loss
    // identically and so says nothing about it, and is held to the same
    // sum. The motion's whole energy is m v^2 / 2 = 5e-5 J.
    let ledger = system.ledger();
    let booked_j = ledger.nodes[cavity.index()].numerical_loss_j.value();
    let ledger_j = ledger.energy_residual_j();
    let a_step_j =
        gas_kg * controls().energy_tolerance_j_kg + 64.0 * f64::EPSILON * 2.0 * P0 * V0;
    let budget_j = f64::from(steps) * a_step_j;
    eprintln!(
        "piston, dt {dt:e} s: {steps} steps, omega {omega:.12} rad/s, gas numerical loss booked \
         {booked_j:.3e} J, ledger residual {ledger_j:.3e} J (budget for each {budget_j:.3e} J)"
    );
    close(booked_j, 0.0, budget_j);
    close(ledger_j, 0.0, budget_j);
    omega
}

#[test]
fn an_insulated_cavity_is_the_closed_cavity_without_heat_and_its_pressure_rises_with_the_flow() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for slope in [0.0, 0.002] {
            let old = parcel(slope, 0.0);
            let stale = parcel(0.0, 12.0);
            let closed = GasCavity::new(old.clone(), controls()).unwrap();
            let cavity = closed.clone().insulated().unwrap();
            let ports: Vec<&str> = cavity.ports().iter().map(|port| port.id().as_str()).collect();
            assert_eq!(ports, ["p"]);
            // One input and one output, so a scalar class can be true of
            // it; the two-port cavity keeps the class of a two-input map.
            assert_eq!((cavity.class(), closed.class()), (Class::Monotone1d, Class::Smooth));
            assert_eq!(cavity.stored_energy_j(&old).unwrap().to_bits(), 0.0_f64.to_bits());
            let dt = 0.001;
            let mut below: Option<f64> = None;
            // Flows that take up to a twentieth of the litre in the step,
            // either way, in rising order.
            for flow in [-0.05, -0.02, -1.0e-6, 0.0, 1.0e-6, 0.02, 0.05] {
                let case = Case { state: &old, stale: &stale, dt_s: dt, inputs: &[flow] };
                check_component(&cavity, &case, cx).unwrap();
                check_port_tangent(&cavity, &case, 2.0e-5, 1.0e-5, cx).unwrap();
                // The step is the closed cavity's with no heat, to the bit.
                let (expected, outputs) = step(&closed, &old, dt, [flow, 0.0], cx);
                let mut next = stale.clone();
                let mut pressure = [f64::NAN];
                let proposal = cavity.propose(&old, dt, &[flow], &mut next, &mut pressure, cx).unwrap();
                assert!(matches!(proposal, Proposal::Step(_)));
                assert_eq!(next, expected);
                assert_eq!(pressure[0].to_bits(), outputs[0].to_bits());
                // Its slope, in either direction, is the closed cavity's
                // entry for pressure against flow, and it is positive.
                let (mut forward, mut back, mut pair) = ([0.0], [0.0], [0.0; 2]);
                cavity.jvp(&old, dt, &[flow], &[1.0], &mut forward, cx).unwrap();
                cavity.vjp(&old, dt, &[flow], &[1.0], &mut back, cx).unwrap();
                closed.jvp(&old, dt, &[flow, 0.0], &[1.0, 0.0], &mut pair, cx).unwrap();
                assert_eq!(forward[0].to_bits(), pair[0].to_bits());
                assert_eq!(back[0].to_bits(), pair[0].to_bits());
                assert!(pair[0] > 0.0, "{flow}: {}", pair[0]);
                // And the pressure itself rises from each flow to the next.
                assert!(below.is_none_or(|lower| pressure[0] > lower), "{flow}: {}", pressure[0]);
                below = Some(pressure[0]);
            }
        }

        // What it refuses, leaving the candidate and the outputs alone.
        let old = parcel(0.0, 0.0);
        let cavity = GasCavity::new(old.clone(), controls()).unwrap().insulated().unwrap();
        let stale = parcel(0.002, 7.0);
        let mut next = stale.clone();
        let (mut one, mut two) = ([7.0], [7.0, 8.0]);
        for inputs in [vec![], vec![0.0, 0.0], vec![f64::NAN]] {
            assert!(cavity.propose(&old, 0.1, &inputs, &mut next, &mut one, cx).is_err());
            assert!(cavity.jvp(&old, 0.1, &inputs, &[1.0], &mut one, cx).is_err());
        }
        // Two outputs, two seeds, two slots.
        assert!(cavity.propose(&old, 0.1, &[0.0], &mut next, &mut two, cx).is_err());
        assert!(cavity.jvp(&old, 0.1, &[0.0], &[1.0, 0.0], &mut one, cx).is_err());
        assert!(cavity.vjp(&old, 0.1, &[0.0], &[1.0], &mut two, cx).is_err());
        assert_eq!(next, stale);
        assert_eq!(one[0].to_bits(), 7.0_f64.to_bits());
        assert_eq!(two.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
        // A cavity has an insulating wall or a stream port, not both.
        assert!(cavity.clone().streaming(1.0e-15).is_err());
        let streaming = GasCavity::new(old.clone(), controls()).unwrap().streaming(1.0e-15).unwrap();
        assert!(streaming.insulated().is_err());
    });
}

#[test]
fn a_piston_on_the_insulated_cavity_is_the_one_on_the_closed_cavity_without_heat_to_the_bit() {
    // The same graph twice: the gas as the one-port cavity, and as the
    // two-port cavity with its heat port held at zero. The kernel sees the
    // same unknowns, the same slopes and the same proposals, so every state
    // is equal to the bit, step for step; nothing is allowed for.
    let (dt, rate) = (2.5e-4, 0.01);
    let (mut one, body_one, gas_one) = spring_piston(true, true, rate);
    let (mut two, body_two, gas_two) = spring_piston(false, true, rate);
    assert_eq!(one.plan_loops(), two.plan_loops());
    let mut iterations = 0_u32;
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for step in 0..400 {
            let a = one.step(dt, cx).unwrap();
            let b = two.step(dt, cx).unwrap();
            assert_eq!(a.iterations, b.iterations, "step {step}");
            iterations += a.iterations;
            let wall_one = one.state::<Motion>(body_one).unwrap();
            let wall_two = two.state::<Motion>(body_two).unwrap();
            assert_eq!(wall_one.position_m.to_bits(), wall_two.position_m.to_bits(), "step {step}");
            assert_eq!(wall_one.velocity_m_s.to_bits(), wall_two.velocity_m_s.to_bits(), "step {step}");
            assert_eq!(
                one.state::<ClosedGasParcel>(gas_one).unwrap(),
                two.state::<ClosedGasParcel>(gas_two).unwrap(),
                "step {step}"
            );
        }
    });
    // Nearly two periods (1.88): the wall is inside its rest position and
    // slower than it began, and the solve iterated.
    let wall = one.state::<Motion>(body_one).unwrap();
    assert!(wall.position_m < V0 / PISTON_M2 - 1.0e-5, "{}", wall.position_m);
    assert!((wall.velocity_m_s - rate).abs() > 1.0e-3 && iterations >= 400, "{} {iterations}", wall.velocity_m_s);
    assert_eq!(
        one.ledger().energy_residual_j().to_bits(),
        two.ledger().energy_residual_j().to_bits()
    );
    eprintln!(
        "insulated against closed with no heat: 400 steps equal to the bit, {iterations} iterations, \
         wall at {:.9e} m/s, ledger residual {:.3e} J",
        wall.velocity_m_s,
        one.ledger().energy_residual_j()
    );
}

#[test]
fn an_insulated_cavity_in_a_loop_needs_the_tangent_it_offers() {
    // What its class costs today. The kernel differences a component that
    // offers no tangent only if its class is smooth or exactly linear. The
    // two-port cavity is smooth and steps without its tangent; the
    // insulated one is a scalar monotone relation, and without its tangent
    // its first step is refused, naming it and what it lacks.
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut smooth, _, _) = spring_piston(false, false, 0.01);
        assert!(smooth.step(2.5e-4, cx).unwrap().iterations >= 1);
        let (mut monotone, _, _) = spring_piston(true, false, 0.01);
        let error = monotone.step(2.5e-4, cx).unwrap_err();
        assert!(
            matches!(&error, SystemError::Unsupported { node, .. } if node == "gas"),
            "{error}"
        );
        assert_eq!(monotone.steps(), 0);
    });
}

#[test]
fn a_small_oscillation_about_rest_has_the_frequency_of_the_adiabatic_gas_spring() {
    // omega0^2 = gamma p A^2 / (m V), with gamma and p from the parcel at
    // its state point.
    let state = *parcel(0.0, 0.0).state().state();
    let gamma = state.gamma;
    close(gamma, 1.4, 1.0e-12);
    let omega0 = (gamma * state.pressure * PISTON_M2 * PISTON_M2 / (PISTON_KG * V0)).sqrt();
    let (dt, rate, periods) = (2.5e-4, 0.01, 4_usize);

    // What the assertion has to allow, each from its expansion.
    //
    // The amplitude. With x the displacement and l = V / A, the force is
    // p A ((1 + x / l)^-gamma - 1), so
    //   x'' + omega0^2 x = -alpha x^2 - beta x^3,
    //   alpha = -omega0^2 (gamma + 1) / (2 l),
    //   beta = omega0^2 (gamma + 1) (gamma + 2) / (6 l^2),
    // and to second order in the amplitude a the frequency moves by
    //   (3 beta / (8 omega0) - 5 alpha^2 / (12 omega0^3)) a^2
    //     = -omega0 (gamma + 1) (2 gamma - 1) / 48 * (a / l)^2
    // (Landau and Lifshitz, Mechanics, section 28). Here a = rate / omega0
    // and a / l = 8.5e-4, so that is 6.4e-8 of omega0; twice it is allowed,
    // for the orders above.
    let swing = rate / omega0 * PISTON_M2 / V0;
    let amplitude = 2.0 * (gamma + 1.0) * (2.0 * gamma - 1.0) / 48.0 * swing * swing;
    // The crossing times. A straight line through two samples of a sine a
    // fraction s and 1 - s of a step either side of its zero misses the
    // zero by dt (omega dt)^2 s (1 - s) (1 - 2 s) / 6, at most
    // 0.0161 dt (omega dt)^2. The velocity is odd about each of its zeros
    // (they are the turning points; the step is a reversible map, so the
    // discrete solution keeps that symmetry), so the second harmonic the
    // gas adds,
    // h = (gamma + 1) / 6 * (a / l) of the first in the velocity, adds no
    // quadratic term there: it changes the slope by at most 2 h and the
    // cubic term by at most 8 h of theirs. Two crossings over `periods`
    // periods.
    let harmonic = (gamma + 1.0) / 6.0 * swing;
    let turn = |dt: f64| {
        let angle = omega0 * dt;
        2.0 * 0.0161 * angle.powi(3) * (1.0 + 8.0 * harmonic) / (1.0 - 2.0 * harmonic)
            / (core::f64::consts::TAU * periods as f64)
    };
    // The step. Linearised about rest, the cavity's discrete-gradient
    // pressure is the pressure at the mean volume and the mass's rule is
    // the midpoint rule, so the graph is the implicit midpoint rule on the
    // oscillator: second order, with the frequency
    //   (2 / dt) atan(omega0 dt / 2) = omega0 (1 - (omega0 dt)^2 / 12 + ...),
    // a lag of 7.3e-5 at this step.
    let midpoint = |dt: f64| 2.0 / dt * (0.5 * omega0 * dt).atan();
    let lag = (omega0 * dt).powi(2) / 12.0;

    // Together 1.6e-7 (1.3e-7 and 0.3e-7), a four-hundredth of the lag.
    let allowed = amplitude + turn(dt);
    assert!(allowed < 2.0e-7 && lag > 400.0 * allowed, "{allowed:e} {lag:e}");
    // With errors e1 and e2 of at most `allowed` at the two steps, the
    // ratio of the lags is (lag + e1) / (lag / 4 + e2), within this of 4.
    let ratio_allowed = 5.0 * allowed / (0.25 * lag - allowed);

    let measured = piston_frequency(dt, rate, periods);
    let halved = piston_frequency(0.5 * dt, rate, periods);
    let ratio = (omega0 - measured) / (omega0 - halved);
    eprintln!(
        "gas spring: omega0 {omega0:.12} rad/s; at dt {measured:.12}, {:.3e} of omega0 from it \
         (lag {lag:.3e}) and {:.3e} from the midpoint rule's {:.12} (allowed {allowed:.3e}); \
         at dt / 2 {halved:.12}, {:.3e} from the midpoint rule's (allowed {:.3e}); \
         ratio of the lags {ratio:.5} (allowed 4 +- {ratio_allowed:.4})",
        measured / omega0 - 1.0,
        (measured - midpoint(dt)) / omega0,
        midpoint(dt),
        (halved - midpoint(0.5 * dt)) / omega0,
        amplitude + turn(0.5 * dt),
    );
    // Against the closed form, within the step's own lag and the rest.
    assert!((measured / omega0 - 1.0).abs() <= lag + allowed, "{measured} vs {omega0}");
    // The lag is the midpoint rule's, and nothing else is left.
    assert!((measured - midpoint(dt)).abs() <= allowed * omega0, "{measured} vs {}", midpoint(dt));
    assert!(
        (halved - midpoint(0.5 * dt)).abs() <= (amplitude + turn(0.5 * dt)) * omega0,
        "{halved} vs {}",
        midpoint(0.5 * dt)
    );
    // Second order: at half the step the lag is a quarter.
    assert!((ratio - 4.0).abs() <= ratio_allowed, "{ratio}");
}

#[test]
fn source_refusals_bad_inputs_and_cancellation_leave_candidate_outputs_untouched() {
    let old = parcel(0.0, 0.0);
    let cavity = GasCavity::new(old.clone(), controls()).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let original_error = old.discrete_gradient_heat_volume(V0, 2000.0, controls(), || false)
            .unwrap_err();
        let mut next = parcel(0.002, 7.0);
        let stale = next.clone();
        let mut outputs = [7.0, 8.0];
        let error = cavity.propose(&old, 1.0, &[0.0, 2000.0], &mut next, &mut outputs, cx).unwrap_err();
        assert_eq!(error, original_error);
        assert_eq!(next, stale);
        assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
        for (dt, inputs) in [(0.0, vec![0.0, 0.0]), (0.1, vec![0.1, 0.0]), (0.1, vec![f64::NAN, 0.0]), (0.1, vec![0.0])] {
            assert!(cavity.propose(&old, dt, &inputs, &mut next, &mut outputs, cx).is_err());
            assert_eq!(next, stale);
            assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
        }
        assert!(cavity.stored_energy_j(&stale).is_err(), "different source identity");
        assert!(cavity.propose(&old, 0.1, &[0.0; 2], &mut next, &mut [1.0], cx).is_err());
        // A subsequent valid call still produces the original native endpoint.
        let (recovered, _) = step(&cavity, &old, 0.01, [0.0, 100.0], cx);
        close(recovered.state().state().temperature, T0 + 1.0 / 0.625, 2.0e-8);
        cx.cancel_gate().request();
        assert!(matches!(
            cavity.propose(&old, 0.01, &[0.0; 2], &mut next, &mut outputs, cx),
            Err(GasEnergyError::Cancelled)
        ));
        assert_eq!(next, stale);
        assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
        let mut tangent = [7.0, 8.0];
        assert!(matches!(
            cavity.jvp(&old, 0.01, &[0.0; 2], &[1.0; 2], &mut tangent, cx),
            Err(GasEnergyError::Cancelled)
        ));
        assert_eq!(tangent.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
    });
}

/// One authored species for the mixtures: Cp/R = cp_r + slope T, with a
/// formation energy, so that a stream's enthalpy is mostly datum.
fn authored_species(
    name: &str, molar_mass: f64, cp_r: f64, slope: f64, formation_j_mol: f64,
) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    coefficients[7] = formation_j_mol / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: P0, dims: Pressure::DIMS });
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(molar_mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("gas-component-fixture").unwrap(),
        ).unwrap(),
        vec![Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            sources: vec![],
            provenance: provenance(),
        }).unwrap()],
    ).unwrap()
}

fn authored_transport(name: &str) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (property, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2.0e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(property, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), 200.0, 2000.0)
                .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1.0e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: provenance(),
        }).unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(), phase: "gas".into(),
                process: "authored".into(), revision: 0,
            },
            claims,
            vec![],
        ).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

const HEAVY_KG_MOL: f64 = 0.028;
const LIGHT_KG_MOL: f64 = 0.004;

/// A litre at 400 K and `pressure` of two species: mole fraction `heavy` of
/// a gas whose heat capacity rises with temperature, the rest a light
/// monatomic one.
fn mixture(heavy: f64, pressure: f64) -> ClosedGasParcel {
    let (mut species, mut fractions) = (Vec::new(), Vec::new());
    // Formation energies of opposite sign, each far above the thermal part
    // (about 4e5 J/kg at 400 K): 7.1e6 and -1.25e7 J/kg.
    if heavy > 0.0 {
        species.push(authored_species("heavy", HEAVY_KG_MOL, 3.5, 0.001, 2.0e5));
        fractions.push(heavy);
    }
    species.push(authored_species("light", LIGHT_KG_MOL, 2.5, 0.0, -5.0e4));
    fractions.push(1.0 - heavy);
    catalog_parcel(species, fractions, pressure)
}

fn catalog_parcel(
    species: Vec<Nasa9StandardStateModelV1>,
    fractions: Vec<f64>,
    pressure: f64,
) -> ClosedGasParcel {
    let transport = species
        .iter()
        .map(|model| authored_transport(model.species().as_str()))
        .collect();
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            T0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
        .unwrap();
    ClosedGasParcel::new(model, point, V0).unwrap()
}

/// `|a - b|`, which must be a number: a NaN would pass a `max`.
fn gap(a: f64, b: f64) -> f64 {
    let difference = (a - b).abs();
    assert!(difference.is_finite(), "{a} against {b}");
    difference
}

/// Mass [kg] of the heavy species in a mixture.
fn heavy_kg(gas: &ClosedGasParcel) -> f64 {
    let moles = gas.model().thermodynamics().mole_composition().fractions()[0];
    let heavy = moles * HEAVY_KG_MOL;
    gas.mass_kg() * heavy / (heavy + (1.0 - moles) * LIGHT_KG_MOL)
}

/// What a run of [`exchange`] found.
struct Exchange {
    /// Largest differences between graph and reference: mass and heavy mass
    /// (relative to the mass), temperature [K], pressure (relative), and,
    /// for the second cavity, the temperature and pressure its heat and
    /// wall ports returned [K, relative].
    worst: [f64; 6],
    /// Steps in which gas went from the first cavity to the second, and
    /// the other way.
    steps: [u32; 2],
    start: [ClosedGasParcel; 2],
    reference: [ClosedGasParcel; 2],
    ledger_j: f64,
}

/// Two litres of mixtures joined by a linear restriction, the second heated
/// at `heat_w` and its wall taking `compression_m3_s` of volume from it.
///
/// In the graph each is the gas owner's cavity with a stream port, and the
/// restriction is the kernel's ideal element, which carries the enthalpy
/// and the species of the side the gas leaves. The reference does the same
/// with the owner's own calls on two parcels, a step at a time: the pair
/// transfer of `dt * conductance * (p_a - p_b)` from the side at the higher
/// pressure, then the heat and volume step on the second; with `heat_first`
/// the reference takes them in the other order, which the cavity does not.
fn exchange(
    start: [ClosedGasParcel; 2],
    heat_w: f64,
    compression_m3_s: f64,
    steps: u32,
    heat_first: bool,
) -> Exchange {
    let (conductance, dt) = (2.0e-9, 1.0e-3);
    let mass_tolerance_kg = 1.0e-15;

    let mut builder = SystemBuilder::new();
    let catalog = GasCavity::new(start[0].clone(), controls()).unwrap()
        .streaming(mass_tolerance_kg).unwrap();
    let pipe = builder
        .add_with_tangent("pipe", Restriction::for_stream(conductance, &catalog.ports()[2]).unwrap(), ())
        .unwrap();
    let mut nodes = Vec::new();
    for (index, (name, gas)) in ["a", "b"].into_iter().zip(&start).enumerate() {
        let cavity = GasCavity::new(gas.clone(), controls()).unwrap()
            .streaming(mass_tolerance_kg).unwrap();
        let ports: Vec<&str> = cavity.ports().iter().map(|port| port.id().as_str()).collect();
        assert_eq!(ports, ["p", "heat", "stream"]);
        assert_eq!(cavity.ports()[2].width(), 4);
        // A streaming cavity offers no tangent, so it is added without one.
        let node = builder.add(name, cavity, gas.clone()).unwrap();
        builder.bond(
            builder.port(node, "stream").unwrap(), builder.port(pipe, name).unwrap(),
        ).unwrap();
        let held = if index == 1 { [compression_m3_s, heat_w] } else { [0.0, 0.0] };
        for (port, value) in ["p", "heat"].into_iter().zip(held) {
            builder.boundary(
                builder.port(node, port).unwrap(), Box::new(move |_, _, values| values[0] = value),
            ).unwrap();
        }
        nodes.push(node);
    }
    let mut system = builder.seal().unwrap();
    let returned = |system: &fs_system::System, port: &str| {
        system.last_output(system.port(nodes[1], port).unwrap())[0]
    };

    let solve = GasMassTransferSolve {
        donor: controls(),
        receiver: controls(),
        mass_tolerance_kg,
        energy_tolerance_j: 1.0e-9,
    };
    let pressure = |gas: &ClosedGasParcel| gas.state().state().pressure;
    let kelvin = |gas: &ClosedGasParcel| gas.state().state().temperature;
    let mut reference = start.clone();
    let mut worst = [0.0_f64; 6];
    let mut directions = [0_u32; 2];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for step in 0..steps {
            let rate = conductance * (pressure(&reference[0]) - pressure(&reference[1]));
            let (from, to) = if rate >= 0.0 { (0, 1) } else { (1, 0) };
            directions[from] += u32::from(rate != 0.0);
            let cross = |reference: &mut [ClosedGasParcel; 2]| {
                let moved = reference[from]
                    .transfer_mass_to(&reference[to], dt * rate.abs(), solve, || false)
                    .unwrap();
                let (donor, receiver, _) = moved.into_parts();
                reference[from] = donor;
                reference[to] = receiver;
            };
            // Heat and wall work on what the second cavity holds.
            let heat_and_wall = |reference: &mut [ClosedGasParcel; 2]| {
                let trial = reference[1]
                    .discrete_gradient_heat_volume(
                        reference[1].volume_m3() - dt * compression_m3_s,
                        dt * heat_w,
                        controls(),
                        || false,
                    )
                    .unwrap();
                let returned = (
                    trial.thermodynamic_mean().state().temperature,
                    trial.step().work_pressure_pa,
                );
                reference[1] = trial.into_parts().0;
                returned
            };
            let (mean_k, work_pa) = if heat_first {
                let returned = heat_and_wall(&mut reference);
                cross(&mut reference);
                returned
            } else {
                cross(&mut reference);
                heat_and_wall(&mut reference)
            };

            system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
            for (node, expected) in nodes.iter().zip(&reference) {
                let held = system.state::<ClosedGasParcel>(*node).unwrap();
                let mass = expected.mass_kg();
                worst[0] = worst[0].max(gap(held.mass_kg(), mass) / mass);
                worst[1] = worst[1].max(gap(heavy_kg(held), heavy_kg(expected)) / mass);
                worst[2] = worst[2].max(gap(kelvin(held), kelvin(expected)));
                worst[3] = worst[3].max(gap(pressure(held), pressure(expected)) / pressure(expected));
            }
            worst[4] = worst[4].max(gap(returned(&system, "heat"), mean_k));
            worst[5] = worst[5].max(gap(returned(&system, "p"), work_pa) / work_pa);
        }
    });
    Exchange {
        worst,
        steps: directions,
        ledger_j: system.ledger().energy_residual_j(),
        start,
        reference,
    }
}

#[test]
fn two_streaming_cavities_exchange_gas_as_the_owners_pair_transfer_does() {
    // 2 bar and 1 bar, nothing else: no heat, no wall work.
    let run = exchange([mixture(0.75, 2.0e5), mixture(0.25, 1.0e5)], 0.0, 0.0, 40, false);
    let Exchange { worst, start, reference, .. } = &run;
    assert_eq!(run.steps, [40, 0]);
    // Both sides solve the owner's energy equation to its 1e-7 J/kg, which
    // is about 1e-10 K a solve, through different entry points: the cavity
    // from the stream's enthalpy in the species' datum, where the
    // formation energies dominate, the pair transfer relative to the
    // sources. Measured over the 40 steps: masses within 8e-16, 1.7e-11 K
    // and 4.3e-14 of the pressure. The bounds are what 40 solves may leave.
    assert!(worst[0] <= 1.0e-13 && worst[1] <= 1.0e-13, "{worst:?}");
    assert!(worst[2] <= 1.0e-8 && worst[3] <= 1.0e-10, "{worst:?}");

    // What happened, from the reference alone: the donor cooled as it
    // expanded and kept its composition; the receiver warmed and gained
    // heavy gas; mass and each species are conserved.
    let kelvin = |gas: &ClosedGasParcel| gas.state().state().temperature;
    assert!(kelvin(&reference[0]) < T0 - 0.9 && kelvin(&reference[1]) > T0 + 2.5);
    let fraction = |gas: &ClosedGasParcel| heavy_kg(gas) / gas.mass_kg();
    close(fraction(&reference[0]), fraction(&start[0]), 1.0e-13);
    assert!(fraction(&reference[1]) > fraction(&start[1]) + 1.0e-3);
    let total = |gases: &[ClosedGasParcel; 2], of: fn(&ClosedGasParcel) -> f64| of(&gases[0]) + of(&gases[1]);
    close(total(reference, ClosedGasParcel::mass_kg), total(start, ClosedGasParcel::mass_kg), 1.0e-18);
    close(total(reference, heavy_kg), total(start, heavy_kg), 1.0e-17);
    // The graph's own accounts: what one cavity lost in enthalpy the other
    // gained. Measured 2.4e-15 J.
    close(run.ledger_j, 0.0, 1.0e-11);
}

#[test]
fn a_heated_and_compressed_cavity_turns_the_flow_round_and_keeps_the_owners_order() {
    // The second cavity starts 5 kPa below the first, is heated at 400 W
    // and loses a fifth of a litre a second to its wall, so the gas goes
    // into it and then, as its pressure passes the first's, back out. Each
    // step it first takes what crossed and then the heat and the wall
    // work: the reference does them in that order.
    let start = || [mixture(0.75, 1.05e5), mixture(0.25, 1.0e5)];
    let run = exchange(start(), 400.0, 2.0e-4, 80, false);
    let Exchange { worst, reference, .. } = &run;
    // Twenty steps in, sixty back out; the second cavity ends at 472 K.
    assert!(run.steps[0] > 10 && run.steps[1] > 10, "{:?}", run.steps);
    assert!(reference[1].state().state().temperature > T0 + 60.0);
    // Measured: masses within 5.5e-14, 9.0e-10 K and 2.2e-12 of the
    // pressure, after 160 of the owner's solves at about 1e-10 K each.
    assert!(worst[0] <= 1.0e-12 && worst[1] <= 1.0e-12, "{worst:?}");
    assert!(worst[2] <= 5.0e-8 && worst[3] <= 2.0e-10, "{worst:?}");
    // The heat port returned the heat-and-volume step's mean temperature
    // and the wall port its work pressure: measured 4.6e-10 K and 1.1e-12.
    assert!(worst[4] <= 5.0e-8 && worst[5] <= 2.0e-10, "{worst:?}");
    // These bounds tell the order: a reference that heats and compresses
    // first and transfers after is another trajectory. Measured against
    // it: 3.7e-3 K in the state and 2.0e-2 K in the returned temperature.
    let other = exchange(start(), 400.0, 2.0e-4, 80, true);
    assert!(other.worst[2] > 1.0e-4 && other.worst[4] > 1.0e-4, "{:?}", other.worst);
    // After the turn the first cavity holds some of the second's lighter
    // mixture.
    let fraction = |gas: &ClosedGasParcel| heavy_kg(gas) / gas.mass_kg();
    assert!(fraction(&reference[0]) < fraction(&run.start[0]) - 1.0e-6);
    // The ledger closes on the heat supplied and the wall's work. Measured
    // 9.1e-13 J of the 32 J of heat.
    close(run.ledger_j, 0.0, 1.0e-10);
}

#[test]
fn a_streaming_cavity_conforms_with_gas_entering_and_leaving_and_refuses_what_it_cannot_do() {
    let gas = mixture(0.6, 1.5e5);
    let stale = mixture(0.3, 1.0e5);
    assert!(GasCavity::new(gas.clone(), controls()).unwrap().streaming(0.0).is_err());
    assert!(GasCavity::new(gas.clone(), controls()).unwrap().streaming(f64::NAN).is_err());
    let cavity = GasCavity::new(gas.clone(), controls()).unwrap().streaming(1.0e-15).unwrap();
    // Asking twice adds the port once.
    assert_eq!(cavity.clone().streaming(1.0e-14).unwrap().ports().len(), 3);
    let enthalpy = gas.state().specific_enthalpy_j_kg();
    let fraction = heavy_kg(&gas) / gas.mass_kg();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        // Compression, heat, and a stream in each direction: its own gas
        // leaving, and heavy gas alone arriving with twice its enthalpy.
        for (rate, heavy, carried) in [(-2.0e-4, fraction, enthalpy), (3.0e-4, 1.0, 2.0 * enthalpy)] {
            let inputs = [
                2.0e-4, 15.0, rate, rate * carried, rate * heavy, rate * (1.0 - heavy),
            ];
            let case = Case { state: &gas, stale: &stale, dt_s: 1.0e-3, inputs: &inputs };
            let Ok(Observed::Stepped { outputs, .. }) = check_component(&cavity, &case, cx) else {
                panic!("no conforming step at {inputs:?}");
            };
            // The stream side is what the cavity held when the step began.
            close(outputs[2], gas.state().state().pressure, 0.0);
            close(outputs[3], enthalpy, 0.0);
            close(outputs[4], fraction, 1.0e-15);
            close(outputs[4] + outputs[5], 1.0, 1.0e-15);
            // The mass moved by what the species rates carried.
            let mut next = stale.clone();
            let mut slot = [0.0; 6];
            cavity.propose(&gas, 1.0e-3, &inputs, &mut next, &mut slot, cx).unwrap();
            close(next.mass_kg(), gas.mass_kg() + 1.0e-3 * rate, 1.0e-18);
            close(heavy_kg(&next), heavy_kg(&gas) + 1.0e-3 * rate * heavy, 1.0e-18);
            close(next.volume_m3(), V0 - 2.0e-7, 1.0e-18);
        }
        // With nothing crossing it is the closed cavity: the same endpoint
        // and the same pressure and temperature, to the bit.
        let closed = GasCavity::new(gas.clone(), controls()).unwrap();
        let (mut alone, mut shut) = (stale.clone(), stale.clone());
        let (mut six, mut two) = ([0.0; 6], [0.0; 2]);
        cavity.propose(&gas, 1.0e-3, &[2.0e-4, 15.0, 0.0, 0.0, 0.0, 0.0], &mut alone, &mut six, cx).unwrap();
        closed.propose(&gas, 1.0e-3, &[2.0e-4, 15.0], &mut shut, &mut two, cx).unwrap();
        assert_eq!(alone, shut);
        assert_eq!([six[0].to_bits(), six[1].to_bits()], two.map(f64::to_bits));
        // An increment too small to change what is held is none: the mass
        // is the same to the bit and the step is taken. (The owner refuses
        // such an increment; a flow dying away would end there.)
        let faint = 1.0e-22;
        let mut next = stale.clone();
        let step = cavity.propose(
            &gas, 1.0e-3,
            &[0.0, 0.0, faint, faint * enthalpy, faint * fraction, faint * (1.0 - fraction)],
            &mut next, &mut six, cx,
        ).unwrap();
        assert!(matches!(step, Proposal::Step(_)));
        assert_eq!(next.mass_kg().to_bits(), gas.mass_kg().to_bits());

        // No tangent is offered, in either direction, and it says so.
        let mut outputs = [7.0; 6];
        let no_tangent = |result: Result<(), GasEnergyError>| {
            matches!(result, Err(GasEnergyError::Invalid(why)) if why.contains("no port tangent"))
        };
        assert!(no_tangent(cavity.jvp(&gas, 1.0e-3, &[0.0; 6], &[1.0; 6], &mut outputs, cx)));
        assert!(no_tangent(cavity.vjp(&gas, 1.0e-3, &[0.0; 6], &[1.0; 6], &mut outputs, cx)));
        // A step that would draw more light gas than the cavity holds, or
        // all of its gas, is shortened; the candidate is left alone.
        let mut next = stale.clone();
        for inputs in [
            [0.0, 0.0, -10.0, -10.0 * enthalpy, 0.0, -10.0],
            [0.0, 0.0, -10.0, -10.0 * enthalpy, -10.0 * fraction, -10.0 * (1.0 - fraction)],
        ] {
            let asked = cavity.propose(&gas, 1.0e-3, &inputs, &mut next, &mut outputs, cx).unwrap();
            let Proposal::Shorten { dt_s } = asked else { panic!("not shortened: {inputs:?}") };
            assert_eq!(dt_s.to_bits(), 0.5e-3_f64.to_bits(), "{inputs:?}");
            assert_eq!(next, stale);
            assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(); 6]);
        }
        // Refusals, before the exchange and after it (a wall that would
        // take more than the whole volume), leave the candidate alone too.
        for inputs in [
            vec![0.0; 5],
            vec![0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0],
            vec![2.0, 0.0, 1.0e-4, 1.0e-4 * enthalpy, 1.0e-4 * fraction, 1.0e-4 * (1.0 - fraction)],
        ] {
            assert!(cavity.propose(&gas, 1.0e-3, &inputs, &mut next, &mut outputs, cx).is_err());
            assert_eq!(next, stale);
            assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(); 6]);
        }

        // A cavity that holds none of a species of its reference returns a
        // zero fraction for it and takes it in: the reference is the
        // catalog. Pure light gas, with heavy gas arriving.
        let light = mixture(0.0, 1.0e5);
        assert_eq!(light.model().thermodynamics().components().len(), 1);
        let arriving = [0.0, 0.0, 2.0e-4, 2.0e-4 * enthalpy, 2.0e-4, 0.0];
        let case = Case { state: &light, stale: &stale, dt_s: 1.0e-3, inputs: &arriving };
        let Ok(Observed::Stepped { outputs: held, .. }) = check_component(&cavity, &case, cx) else {
            panic!("no conforming step from light gas alone");
        };
        assert_eq!([held[4].to_bits(), held[5].to_bits()], [0.0_f64.to_bits(), 1.0_f64.to_bits()]);
        let mut mixed = stale.clone();
        cavity.propose(&light, 1.0e-3, &arriving, &mut mixed, &mut six, cx).unwrap();
        assert_eq!(mixed.model().thermodynamics().components().len(), 2);
        close(heavy_kg(&mixed), 2.0e-7, 1.0e-20);
        close(mixed.mass_kg(), light.mass_kg() + 2.0e-7, 1.0e-18);

        // Cancellation, last: the gate stays requested.
        cx.cancel_gate().request();
        assert!(matches!(
            cavity.propose(&gas, 1.0e-3, &arriving, &mut next, &mut outputs, cx),
            Err(GasEnergyError::Cancelled)
        ));
        assert_eq!(next, stale);
        assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(); 6]);
    });
}

fn connect_streams(a: GasCavity, b: GasCavity) -> Result<(), AdmissionError> {
    let mut builder = SystemBuilder::new();
    let pipe = Restriction::for_stream(2.0e-9, &a.ports()[2]).unwrap();
    let pipe = builder.add_with_tangent("pipe", pipe, ()).unwrap();
    for (name, cavity) in [("a", a), ("b", b)] {
        let old = cavity.reference().clone();
        let node = builder.add(name, cavity, old).unwrap();
        builder.bond(
            builder.port(node, "stream").unwrap(),
            builder.port(pipe, name).unwrap(),
        )?;
    }
    Ok(())
}

fn streaming(gas: ClosedGasParcel) -> GasCavity {
    GasCavity::new(gas, controls())
        .unwrap()
        .streaming(1.0e-15)
        .unwrap()
}

#[test]
fn stream_catalog_identity_ignores_live_mixture_and_state() {
    let a = streaming(mixture(0.75, 1.0e5));
    let b = mixture(0.25, 1.5e5)
        .discrete_gradient_heat_volume(1.1 * V0, 2.0, controls(), || false)
        .unwrap()
        .into_parts()
        .0;
    assert_ne!(b.state().state().temperature, T0);
    let b = streaming(b);
    assert_eq!(a.ports()[2].stream_basis(), b.ports()[2].stream_basis());
    connect_streams(a, b).unwrap();
}

#[test]
fn equal_width_streams_refuse_different_species_sources_and_molar_masses() {
    let model = authored_species("fixture", MOLAR_MASS, 3.5, 0.001, 100.0);
    let original = || streaming(catalog_parcel(vec![model.clone()], vec![1.0], P0));
    let different_mass = Nasa9StandardStateModelV1::new(
        model.species().clone(),
        MolarMass::new(2.0 * MOLAR_MASS),
        model.convention().clone(),
        model.regions().to_vec(),
    )
    .unwrap();
    let different_reference = Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("different-energy-reference").unwrap(),
        )
        .unwrap(),
        model.regions().to_vec(),
    )
    .unwrap();
    for other in [
        authored_species("other", MOLAR_MASS, 3.5, 0.001, 100.0),
        authored_species("fixture", MOLAR_MASS, 3.5, 0.001, 300.0),
        different_mass,
        different_reference,
    ] {
        let changed_state = catalog_parcel(vec![other], vec![1.0], P0);
        let cavity = original();
        assert!(cavity.stored_energy_j(&changed_state).is_err());
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let mut next = cavity.reference().clone();
            let stale = next.clone();
            let mut outputs = [17.0; 5];
            assert!(
                cavity
                    .propose(
                        &changed_state,
                        1.0e-3,
                        &[0.0; 5],
                        &mut next,
                        &mut outputs,
                        cx,
                    )
                    .is_err()
            );
            assert_eq!(next, stale);
            assert_eq!(outputs, [17.0; 5]);
        });
        let other = streaming(changed_state);
        assert_eq!(original().ports()[2].width(), other.ports()[2].width());
        assert!(matches!(
            connect_streams(original(), other),
            Err(AdmissionError::StreamBasis { .. })
        ));
    }
}

#[test]
fn the_same_source_cards_with_another_caloric_datum_do_not_connect() {
    let model = authored_species("fixture", MOLAR_MASS, 3.5, 0.0, 100.0);
    let mut lower = model.regions()[0].card().clone();
    let mut upper = lower.clone();
    lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 600.0);
    upper.validity = ValidityDomain::unconstrained().with("T", 600.0, 2000.0);
    upper.parameters.get_mut("a7").unwrap().value += 75.0 / R;
    let raw = Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        model.convention().clone(),
        vec![
            Nasa9RegionV1::from_card(lower).unwrap(),
            Nasa9RegionV1::from_card(upper).unwrap(),
        ],
    )
    .unwrap();
    let corrected = raw
        .with_caloric_continuity(
            Temperature::new(800.0),
            fs_thermochem::MolarEnergyQuantityV1::new(100.0),
            fs_thermochem::MolarThermalQuantityV1::new(1.0),
        )
        .unwrap();
    assert!(
        corrected.caloric_reference().unwrap().regions()[0]
            .enthalpy_offset()
            .value()
            > 74.0
    );
    let raw = streaming(catalog_parcel(vec![raw], vec![1.0], P0));
    let corrected = streaming(catalog_parcel(vec![corrected], vec![1.0], P0));
    let a = raw.ports()[2].stream_basis().unwrap();
    let b = corrected.ports()[2].stream_basis().unwrap();
    assert_eq!(a.medium(), b.medium());
    assert_ne!(a.enthalpy_reference(), b.enthalpy_reference());
    assert!(matches!(
        connect_streams(raw, corrected),
        Err(AdmissionError::StreamBasis { .. })
    ));
}
