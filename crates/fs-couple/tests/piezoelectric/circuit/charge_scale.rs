//! Authored ceramic-scale constants, not a measured PZT material card.
//! The charge support is much smaller than the native generic Newton probe.
use super::*;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::*;
use fs_qty::Density;

const LENGTH: f64 = 1e-3;
const AREA: f64 = 1e-12;
const DENSITY: f64 = 7500.0;
const COMPLIANCE: f64 = 20e-12;
const PERMITTIVITY: f64 = 20e-9;

fn micro_rod(d: f64) -> PiezoelectricRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, COMPLIANCE),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, d),
        (PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, PERMITTIVITY),
        ("density", Density::DIMS, DENSITY),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored microscopic reciprocal ceramic fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic ceramic".into(), phase: "solid".into(),
        process: "analytic small-electrode regression".into(), revision: 0,
    }, claims, vec![]).unwrap();
    PiezoelectricRod::new(LongitudinalPiezoelectric::resolve(&card, &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly).unwrap(), LENGTH, AREA,
        SmallSignalLimits { max_abs_strain: 1e-3, max_abs_electric_field_v_m: 1e4 }).unwrap()
}

fn relative(actual: f64, expected: f64, scale: f64) {
    assert!((actual-expected).abs() <= 2e-10*scale,
        "actual {actual:e}, expected {expected:e}, physical scale {scale:e}");
}

// Independently eliminate the physical midpoint equations, with no native
// matrix, solver, material-response callback or state-coordinate ordering.
fn midpoint(x: f64, v: f64, q: f64, d: f64, resistance: f64, force: f64, h: f64)
    -> [f64; 4] {
    let mass = DENSITY*AREA*LENGTH/3.0;
    let stiffness = AREA/(COMPLIANCE*LENGTH);
    let theta = d*AREA/(COMPLIANCE*LENGTH);
    let capacitance = (PERMITTIVITY-d*d/COMPLIANCE)*AREA/LENGTH;
    let denominator = 1.0+h*h*stiffness/(4.0*mass);
    let a = (x+0.5*h*v+h*h*force/(4.0*mass))/denominator;
    let b = h*h*theta/(4.0*mass*denominator);
    let voltage = (q-theta*a)/(capacitance+h/(2.0*resistance)+theta*b);
    let midpoint_x = a+b*voltage;
    [2.0*midpoint_x-x, 4.0*(midpoint_x-x)/h-v,
        q-h*voltage/resistance, voltage]
}

fn model(d: f64, resistance: f64) -> PiezoelectricCircuit {
    PiezoelectricCircuit::new(CircuitGraph {
        node_count: 2,
        branches: vec![
            (1,0,Branch::Capacitor { farads: f64::NAN }),
            (1,0,Branch::Resistor { ohms: f64::NAN }),
        ], transformers: vec![],
    }, vec![source(1,resistance)], vec![], vec![6e-18], 0,
        micro_rod(d), 1e-10, 0.0, &[], || false).unwrap()
}

#[test]
fn g1_subfemtocoulomb_circuit_matches_independent_loaded_ceramic_equations() {
    let controls = MaterialCircuitSolve { energy_tolerance_j: 1e-26, power_tolerance_w: 1e-18 };
    let h = 2e-8;
    let force = 1e-8;
    for d in [-400e-12, 0.0, 400e-12] {
        for resistance in [1e9, 4e9] {
            let mut body = model(d,resistance);
            let (mut x,mut v,mut q) = (1e-10,0.0,6e-18);
            for _ in 0..8 {
                let expected = midpoint(x,v,q,d,resistance,force,h);
                let original = body.clone();
                let trial = body.step_candidate(vec![source(1,resistance)], &[], force,
                    h, controls, || false).unwrap();
                relative(trial.model().extension_m(), expected[0], 1e-10);
                relative(trial.model().velocity_m_s(), expected[1], 1e-3);
                relative(trial.model().charge_c(), expected[2], 1e-17);
                relative(trial.step().heat_j, h*expected[3]*expected[3]/resistance, 1e-17);
                assert!(trial.step().energy_residual_j.abs() <= controls.energy_tolerance_j);
                assert_eq!(body, original);
                (x,v,q) = (expected[0],expected[1],expected[2]);
                body = trial.into_parts().0;
            }
        }
    }
}

#[test]
fn g3_exact_circuit_tangent_preserves_material_envelope_and_retry() {
    let body = model(400e-12,1e9);
    let original = body.clone();
    let controls = MaterialCircuitSolve { energy_tolerance_j: 1e-26, power_tolerance_w: 1e-18 };
    assert!(body.step_candidate(vec![source(1,1e9)], &[], 1.0,
        2e-8, controls, || false).is_err());
    assert!(body.step_candidate(vec![source(1,1e9)], &[], 1e-8,
        2e-8, controls, || true).is_err());
    assert_eq!(body,original);
    let accepted = body.step_candidate(vec![source(1,1e9)], &[], 1e-8,
        2e-8, controls, || false).unwrap();
    let retry = body.step_candidate(vec![source(1,1e9)], &[], 1e-8,
        2e-8, controls, || false).unwrap();
    assert_eq!(accepted.model(),retry.model());
}

// This oracle integrates the physical ceramic equations directly, without the
// descriptor, source-response callbacks, transformed coordinates or cached flow.
fn exact_reference(mut y: [f64; 4], d: f64, resistance: f64, force: f64, h: f64) -> [f64; 4] {
    let mass = DENSITY * AREA * LENGTH / 3.0;
    let stiffness = AREA / (COMPLIANCE * LENGTH);
    let theta = d * stiffness;
    let capacitance = (PERMITTIVITY - d * d / COMPLIANCE) * AREA / LENGTH;
    let derivative = |state: [f64; 4]| {
        let voltage = (state[2] - theta * state[0]) / capacitance;
        [
            state[1],
            (force - stiffness * state[0] + theta * voltage) / mass,
            -voltage / resistance,
            voltage * voltage / resistance,
        ]
    };
    y[3] = 0.0;
    let dt = h / 512.0;
    for _ in 0..512 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let e = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..4 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + e[i]) / 6.0;
        }
    }
    y
}

#[test]
fn g1_prepared_ceramic_circuit_preserves_subfemtocoulomb_state_and_physical_loss() {
    use fs_couple::material_circuit::piezoelectric::prepared::PreparedPiezoelectricCircuit;
    use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};
    let controls = MaterialCircuitSolve {
        energy_tolerance_j: 1e-26,
        power_tolerance_w: 1e-18,
    };
    for d in [-400e-12, 0.0, 400e-12] {
        let mut prepared = PreparedPiezoelectricCircuit::new(
            model(d, 1e9),
            &[],
            2e-8,
            LinearSplitBudget {
                max_modes: 8,
                max_omega_dt: 3.0,
            },
            DenseExponentialBudget {
                max_dimension: 16,
                max_dense_bytes: 8 * 1024 * 1024,
                max_multiply_terms: 100_000_000,
            },
            || false,
        )
        .unwrap();
        let mut expected = [1e-10, 0.0, 6e-18, 0.0];
        for _ in 0..8 {
            expected = exact_reference(expected, d, 1e9, 1e-8, 2e-8);
            let report = prepared.step(&[], 1e-8, controls, || false).unwrap();
            relative(prepared.model().extension_m(), expected[0], 1e-10);
            relative(prepared.model().velocity_m_s(), expected[1], 1e-3);
            relative(prepared.model().charge_c(), expected[2], 1e-17);
            relative(report.heat_j, expected[3], 1e-18);
            assert_eq!(report.electrical_supply_j, 0.0);
            assert!(report.energy_residual_j.abs() <= controls.energy_tolerance_j);
            assert_eq!(prepared.resistor_heat().unwrap()[0].heat_j, report.heat_j);
        }
    }
}
