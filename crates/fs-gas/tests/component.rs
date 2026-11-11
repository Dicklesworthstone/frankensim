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
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec};
use fs_system::elements::{Mass, Motion, Restriction};
use fs_system::{SolveConfig, SystemBuilder};
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

/// One authored species for the mixtures: Cp/R = cp_r + slope T, no
/// formation energy.
fn authored_species(name: &str, molar_mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
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
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                authored_species("heavy", HEAVY_KG_MOL, 3.5, 0.001),
                authored_species("light", LIGHT_KG_MOL, 2.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![heavy, 1.0 - heavy]).unwrap(),
        ).unwrap(),
        vec![authored_transport("heavy"), authored_transport("light")],
        None,
    ).unwrap();
    let point = QueryPoint::new()
        .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), T0).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure).unwrap();
    ClosedGasParcel::new(model, point, V0).unwrap()
}

/// Mass [kg] of the heavy species in a mixture.
fn heavy_kg(gas: &ClosedGasParcel) -> f64 {
    let moles = gas.model().thermodynamics().mole_composition().fractions()[0];
    let heavy = moles * HEAVY_KG_MOL;
    gas.mass_kg() * heavy / (heavy + (1.0 - moles) * LIGHT_KG_MOL)
}

#[test]
fn two_streaming_cavities_exchange_gas_as_the_owners_pair_transfer_does() {
    // Two litres of different mixtures at 2 bar and 1 bar, joined by a
    // linear restriction. In the graph each is the gas owner's cavity with a
    // stream port, and the restriction is the kernel's ideal element, which
    // carries the enthalpy and the species of the side the gas leaves. The
    // reference moves the same mass between the two parcels with the
    // owner's own pair transfer, a step at a time.
    let (conductance, dt, steps) = (2.0e-9, 1.0e-3, 40_u32);
    let start = [mixture(0.75, 2.0e5), mixture(0.25, 1.0e5)];
    let mass_tolerance_kg = 1.0e-15;

    let mut builder = SystemBuilder::new();
    let pipe = builder
        .add_with_tangent("pipe", Restriction::new(conductance, 2).unwrap(), ())
        .unwrap();
    let mut nodes = Vec::new();
    for (name, gas) in ["a", "b"].into_iter().zip(&start) {
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
        for held in ["p", "heat"] {
            builder.boundary(
                builder.port(node, held).unwrap(), Box::new(|_, _, values| values[0] = 0.0),
            ).unwrap();
        }
        nodes.push(node);
    }
    let mut system = builder.seal().unwrap();

    let solve = GasMassTransferSolve {
        donor: controls(),
        receiver: controls(),
        mass_tolerance_kg,
        energy_tolerance_j: 1.0e-9,
    };
    let pressure = |gas: &ClosedGasParcel| gas.state().state().pressure;
    let kelvin = |gas: &ClosedGasParcel| gas.state().state().temperature;
    let mut reference = start.clone();
    // Largest differences: mass, heavy mass (both relative to the mass),
    // temperature [K], pressure (relative).
    let mut worst = [0.0_f64; 4];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for step in 0..steps {
            let rate = conductance * (pressure(&reference[0]) - pressure(&reference[1]));
            assert!(rate > 0.0, "step {step}");
            let moved = reference[0]
                .transfer_mass_to(&reference[1], dt * rate, solve, || false)
                .unwrap();
            let (donor, receiver, _) = moved.into_parts();
            reference = [donor, receiver];
            system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
            for (node, expected) in nodes.iter().zip(&reference) {
                let held = system.state::<ClosedGasParcel>(*node).unwrap();
                let mass = expected.mass_kg();
                worst[0] = worst[0].max((held.mass_kg() - mass).abs() / mass);
                worst[1] = worst[1].max((heavy_kg(held) - heavy_kg(expected)).abs() / mass);
                worst[2] = worst[2].max((kelvin(held) - kelvin(expected)).abs());
                worst[3] = worst[3]
                    .max((pressure(held) - pressure(expected)).abs() / pressure(expected));
            }
        }
    });
    // Measured: the masses equal to the bit, 5.7e-13 K and 1.3e-15 of the
    // pressure. Both sides solve the owner's energy equation to its
    // 1e-7 J/kg, through different entry points.
    assert!(worst[0] <= 1.0e-15 && worst[1] <= 1.0e-15, "{worst:?}");
    assert!(worst[2] <= 1.0e-10 && worst[3] <= 1.0e-13, "{worst:?}");

    // What happened, from the reference alone: the donor cooled as it
    // expanded (to 399.07 K) and kept its composition; the receiver warmed
    // (to 402.52 K) and gained heavy gas; mass and each species are
    // conserved.
    assert!(kelvin(&reference[0]) < T0 - 0.9 && kelvin(&reference[1]) > T0 + 2.5);
    let fraction = |gas: &ClosedGasParcel| heavy_kg(gas) / gas.mass_kg();
    close(fraction(&reference[0]), fraction(&start[0]), 1.0e-13);
    assert!(fraction(&reference[1]) > fraction(&start[1]) + 1.0e-3);
    let total = |gases: &[ClosedGasParcel; 2], of: fn(&ClosedGasParcel) -> f64| of(&gases[0]) + of(&gases[1]);
    close(total(&reference, ClosedGasParcel::mass_kg), total(&start, ClosedGasParcel::mass_kg), 1.0e-18);
    close(total(&reference, heavy_kg), total(&start, heavy_kg), 1.0e-17);
    // The graph's own accounts: no heat, no wall work, what one cavity lost
    // in enthalpy the other gained. Measured 2.5e-16 J.
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-12);
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
        // No tangent is offered, and the refusals leave the candidate alone.
        let mut next = stale.clone();
        let mut outputs = [7.0; 6];
        assert!(cavity.jvp(&gas, 1.0e-3, &[0.0; 6], &[1.0; 6], &mut outputs, cx).is_err());
        for inputs in [
            vec![0.0; 5],
            vec![0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0],
            // More light gas drawn than the cavity holds.
            vec![0.0, 0.0, -10.0, -10.0 * enthalpy, 0.0, -10.0],
        ] {
            assert!(cavity.propose(&gas, 1.0e-3, &inputs, &mut next, &mut outputs, cx).is_err());
            assert_eq!(next, stale);
            assert_eq!(outputs.map(f64::to_bits), [7.0_f64.to_bits(); 6]);
        }
    });
}
