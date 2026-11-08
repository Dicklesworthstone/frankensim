//! Physical hold/dwell paths through the typed spatial material-point protocol.
//! Coefficients are synthetic; these are independent discrete-law references.

use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::creep::recovery::{
    FiniteRecoveringCreepState, FiniteRecoveringNortonCreep,
};
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{
    FiniteInelasticMaterial, FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Material,
    FiniteJ2Parameters, FiniteJ2Response, FiniteJ2State,
};

const MU: f64 = 80e6;
const DT: f64 = 0.25;
const ETA: f64 = 40e6;
const A: f64 = 0.025;

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn kinetics() -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: ETA,
        stress_scale_pa: 2e6,
        exponent: 1.0,
    }
}

fn norton() -> FiniteNortonCreep {
    FiniteNortonCreep::new(MU, 160e6, limits(), kinetics(), DT).unwrap()
}

fn plastic() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: 160e6,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 10e6,
        },
        limits(),
    )
    .unwrap()
}

fn stretch(axial: f64) -> [f64; 9] {
    [
        axial.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial).exp(),
    ]
}

fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + 2e-10 * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}",
    );
}

/// This generic consumer cannot inspect or reconstruct a concrete state.
fn candidate<M: FiniteInelasticMaterial>(
    material: &M,
    f: &[f64; 9],
    old: &M::State,
) -> FiniteJ2Response<M::State> {
    // A line-search call and repeated candidate evaluation must consume no
    // extra physical time and must retain the same accepted old memory.
    let (piola, potential) = material.stress_and_potential(f, old).unwrap();
    let response = material.response(f, old).unwrap();
    assert_eq!(piola, response.first_piola_pa);
    assert_eq!(potential, response.incremental_potential_j_m3);
    material
        .validate_committed(f, &response.state, 1e-4)
        .unwrap();
    near(
        material.stored_energy(f, &response.state).unwrap(),
        response.stored_energy_j_m3,
        1e-7,
    );
    response
}

#[test]
fn generic_norton_hold_retains_clock_and_relaxes_below_any_yield_threshold() {
    let law = norton();
    assert_eq!(FiniteInelasticMaterial::duration_s(&law), Some(DT));
    let f = stretch(0.05);
    let mut state = FiniteInelasticMaterial::initial_state(&law);
    let ratio = 1.0 / (1.0 + 3.0 * MU * DT / ETA);
    let mut q = 3.0 * MU * 0.05;
    for _ in 0..9 {
        let old = state.clone();
        let response = candidate(&law, &f, &old);
        assert_eq!(old, state);
        q *= ratio;
        near(response.equivalent_kirchhoff_pa, q, 1e-5);
        near(response.state.alpha, 0.05 - q / (3.0 * MU), 1e-13);
        near(response.plastic_increment, DT * q / ETA, 1e-13);
        let heat = q * response.plastic_increment;
        near(response.inelastic_dissipation_j_m3(), heat, 1e-7);
        near(response.flow_potential_j_m3, 0.5 * heat, 1e-7);
        assert!(response.plastic_increment > 0.0);
        state = response.state;
    }
    assert!(q < 1e4);
}

#[test]
fn generic_dwell_and_reload_keep_recovering_resistance_separate_from_permanent_path() {
    let primary = norton()
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap();
    let law = FiniteRecoveringNortonCreep::new(primary, 2.0).unwrap();
    assert_eq!(FiniteInelasticMaterial::duration_s(&law), Some(DT));
    let free = [2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];
    let mut state = FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [0.25, 4.0, 1.0, 0.0, 0.0, 0.0],
            alpha: 1.2,
            creep_alpha: 0.0,
        },
        resistance_strain: 0.3,
    };
    let permanent = state.native.clone();
    let g = 1.0 + 2.0 * DT;
    let mut resistance = 0.3;
    for _ in 0..5 {
        let response = candidate(&law, &free, &state);
        resistance /= g;
        near(response.state.resistance_strain, resistance, 1e-15);
        assert_eq!(response.state.native, permanent);
        assert_eq!(response.inelastic_dissipation_j_m3(), 0.0);
        assert_eq!(response.flow_potential_j_m3, 0.0);
        state = response.state;
    }

    let increment = stretch(0.04);
    let reload = [
        2.0 * increment[0],
        0.0,
        0.0,
        0.0,
        0.5 * increment[4],
        0.0,
        0.0,
        0.0,
        increment[8],
    ];
    let response = candidate(&law, &reload, &state);
    // For n=m=1 and zero floor the BE kinetic equation is this quadratic
    // in creep increment, with r_new=(r_old+dc)/g. No accumulated alpha
    // enters it; the permanent path still advances by dc.
    let c = DT * A * g / ETA;
    let b = A * g + resistance + 3.0 * MU * c;
    let forcing = c * 3.0 * MU * 0.04;
    let dc = 2.0 * forcing / (b + (b * b + 4.0 * forcing).sqrt());
    near(response.plastic_increment, dc, 1e-13);
    near(response.state.native.alpha, permanent.alpha + dc, 1e-13);
    near(
        response.state.resistance_strain,
        (resistance + dc) / g,
        1e-13,
    );
    assert!(response.state.resistance_strain < response.state.native.alpha);
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            FiniteInelasticMaterial::validate_committed(&law, &reload, &state, tolerance).is_err()
        );
    }
}

struct AdmittedNative(FiniteJ2);
impl FiniteJ2Material for AdmittedNative {
    fn finite_j2(&self) -> &FiniteJ2 {
        &self.0
    }
}

#[test]
fn generic_mechanisms_retain_their_equilibrium_and_overstress_validation() {
    let f = stretch(0.05);
    let rate_independent = AdmittedNative(plastic());
    assert_eq!(rate_independent.duration_s(), None);
    let virgin = rate_independent.initial_state();
    assert!(matches!(
        rate_independent.validate_committed(&f, &virgin, 1e-4),
        Err(FiniteJ2Error::YieldResidual { .. })
    ));
    let returned = candidate(&rate_independent, &f, &virgin);
    near(
        returned.flow_potential_j_m3,
        returned.inelastic_dissipation_j_m3(),
        1e-7,
    );

    let perzyna = FiniteOverstressJ2::new(plastic(), kinetics(), DT).unwrap();
    assert_eq!(FiniteInelasticMaterial::duration_s(&perzyna), Some(DT));
    let returned = candidate(&perzyna, &f, &perzyna.initial_state());
    assert!(returned.equivalent_kirchhoff_pa > 2e6 + 10e6 * returned.state.alpha);
    assert!(returned.flow_potential_j_m3 < returned.inelastic_dissipation_j_m3());

    let combined = FiniteNortonJ2::new(plastic(), kinetics(), DT, 3.0).unwrap();
    assert_eq!(FiniteInelasticMaterial::duration_s(&combined), Some(DT));
    let returned = candidate(&combined, &f, &combined.initial_state());
    assert!(returned.plastic_increment > 0.0 && returned.creep_increment > 0.0);
    near(
        returned.equivalent_kirchhoff_pa,
        2e6 + 10e6 * returned.state.alpha,
        1e-5,
    );
    assert!(returned.plastic_dissipation_j_m3 > 0.0);
    assert!(returned.creep_dissipation_j_m3 > 0.0);
    assert!(returned.flow_potential_j_m3 < returned.inelastic_dissipation_j_m3());
}
