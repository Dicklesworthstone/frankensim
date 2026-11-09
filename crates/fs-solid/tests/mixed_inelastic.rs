//! Native mixed kinetic endpoints against independent uniaxial equations.
//! Synthetic coefficients verify the implementation, not material calibration
//! or general mesh stability. Lateral motion is solved, not prescribed.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::StrainHardeningCreep;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::creep::recovery::FiniteRecoveringNortonCreep;
use fs_material::plastic::finite::mixed::MixedFiniteInelasticMaterial;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State};
use fs_solid::hyper3::mixed::plastic::{
    MixedPlasticLawField, MixedPlasticTetError, MixedPlasticTetSolid,
};
use fs_solid::hyper3::mixed::{MixedUpError, MixedUpSettings};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

const MU: f64 = 80.0;
const YIELD: f64 = 2.0;
const HARDENING: f64 = 8.0;
const ETA: f64 = 80.0;
const DT: f64 = 0.2;

fn with_gate(gate: &CancelGate, f: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 64,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    with_gate(&CancelGate::new(), f);
}

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.25,
        maximum_elastic_log_strain: 0.7,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn backbone(bulk: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: bulk,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: HARDENING,
        },
        limits(),
    )
    .unwrap()
}

fn kinetics() -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: ETA,
        stress_scale_pa: YIELD,
        exponent: 1.0,
    }
}

fn norton(bulk: f64, dt: f64) -> FiniteNortonCreep {
    FiniteNortonCreep::new(MU, bulk, limits(), kinetics(), dt).unwrap()
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 40,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: MU,
        maximum_assembly_bytes: 32 * 1024 * 1024,
        maximum_factor_bytes: 16 * 1024 * 1024,
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 2e-8,
        maximum_endpoint_work_defect_j: 10.0,
        yield_tolerance_pa: 2e-8,
    }
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn solid<'a, L: MixedFiniteInelasticMaterial>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a L,
    exact: bool,
    cx: &Cx<'_>,
) -> MixedPlasticTetSolid<'a, L> {
    MixedPlasticTetSolid::try_new(
        nodes,
        elements,
        MixedPlasticLawField::Uniform(material),
        if exact {
            VolumetricTreatment::ExactIncompressible
        } else {
            VolumetricTreatment::OriginalMaterial
        },
        TetAssemblyBudget::standard(),
        4,
        64,
        cx,
    )
    .unwrap()
}

fn supports(nodes: &[[f64; 3]], axial_log_stretch: f64) -> Vec<(usize, f64)> {
    let mut result = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        result.push((3 * node, axial_log_stretch.exp_m1() * x[0]));
        for axis in 1..3 {
            if x[axis] == 0.0 {
                result.push((3 * node + axis, 0.0));
            }
        }
    }
    result
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.14e}, expected={expected:.14e}, tolerance={tolerance:e}"
    );
}

fn elastic_constants(bulk: f64, exact: bool) -> (f64, f64) {
    if exact {
        (3.0 * MU, 0.5)
    } else {
        (
            9.0 * MU / (3.0 + MU / bulk),
            (3.0 - 2.0 * MU / bulk) / (6.0 + 2.0 * MU / bulk),
        )
    }
}

fn check_metric(state: &FiniteJ2State, net_flow: f64) {
    let metric = [
        (-2.0 * net_flow).exp(),
        net_flow.exp(),
        net_flow.exp(),
        0.,
        0.,
        0.,
    ];
    for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip(metric) {
        near(*actual, expected, 3e-8);
    }
}

// Zero lateral stress implies tau_x=E*(a-net), b=-nu*(a-net)-net/2.
// The scalar expectations use no constitutive response or native tangent.
fn check_uniaxial<L: MixedFiniteInelasticMaterial>(
    solid: &MixedPlasticTetSolid<'_, L>,
    a: f64,
    net: f64,
    plastic: f64,
    bulk: f64,
    exact: bool,
) {
    let (young, poisson) = elastic_constants(bulk, exact);
    let b = -poisson * (a - net) - 0.5 * net;
    let expected_f = [a.exp(), 0., 0., 0., b.exp(), 0., 0., 0., b.exp()];
    let tau = young * (a - net);
    let pressure = tau / (3.0 * (a + 2.0 * b).exp());
    for f in solid.deformation_gradients() {
        for (actual, expected) in f.iter().zip(expected_f) {
            near(*actual, expected, 3e-8);
        }
    }
    for p in solid.pressure_pa() {
        near(*p, pressure, 5e-7);
    }
    for p in solid.first_piola_pa() {
        near(p[0], tau / a.exp(), 5e-7);
        for value in &p[1..] {
            near(*value, 0.0, 5e-7);
        }
    }
    near(
        solid.stored_energy_j(),
        solid.reference_volume_m3()
            * (0.5 * young * (a - net).powi(2) + 0.5 * HARDENING * plastic.powi(2)),
        3e-8,
    );
    if exact {
        assert_eq!(solid.bulk_energy_j(), 0.0);
    }
}

#[test]
fn linear_perzyna_relaxes_without_projecting_accepted_overstress_to_yield() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        for (bulk, exact) in [(160.0, false), (8e8, false), (1e300, true)] {
            let material = FiniteOverstressJ2::new(backbone(bulk), kinetics(), DT).unwrap();
            let mut accepted = solid(&nodes, &elements, &material, exact, cx);
            assert_eq!(accepted.duration_s(), Some(DT));
            let (young, _) = elastic_constants(bulk, exact);
            let a = 0.06;
            let prescribed = supports(&nodes, a);
            let mut alpha = 0.0;
            let mut previous_q = f64::INFINITY;
            for interval in 0..4 {
                let old = alpha;
                let increment = (young * a - (young + HARDENING) * old - YIELD)
                    / (young + HARDENING + ETA / DT);
                alpha += increment;
                let q = young * (a - alpha);
                assert!(increment > 0.0 && q > YIELD + HARDENING * alpha);
                assert!(q < previous_q);
                let step = accepted
                    .advance(&prescribed, &[], [0.0; 3], None, settings(), controls(), cx)
                    .unwrap();
                check_uniaxial(&step.solid, a, alpha, alpha, bulk, exact);
                for state in step.solid.quadrature_states() {
                    near(state.alpha, alpha, 3e-8);
                    assert_eq!(state.creep_alpha, 0.0);
                    check_metric(state, alpha);
                }
                let v = step.solid.reference_volume_m3();
                let viscous = ETA * increment / DT;
                near(
                    step.energy.plastic_dissipation_j,
                    v * (YIELD + viscous) * increment,
                    3e-8,
                );
                assert_eq!(step.energy.creep_dissipation_j, 0.0);
                near(
                    step.energy.flow_potential_j,
                    v * (YIELD + 0.5 * viscous) * increment,
                    3e-8,
                );
                assert!(step.energy.inelastic_dissipation_j() > step.energy.flow_potential_j);
                if interval > 0 {
                    assert!(step.energy.stored_after_j < step.energy.stored_before_j);
                }
                assert!(step.energy.equilibrium_work_residual_j.abs() < 2e-8);
                for old_state in accepted.quadrature_states() {
                    near(old_state.alpha, old, 3e-8);
                }
                previous_q = q;
                accepted = step.solid;
            }
        }
    });
}

#[test]
fn linear_norton_sustained_stress_creep_then_relaxation_matches_scalar_flow() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        for (bulk, exact) in [(160.0, false), (8e8, false), (1e300, true)] {
            let material = norton(bulk, DT);
            let mut accepted = solid(&nodes, &elements, &material, exact, cx);
            let (young, _) = elastic_constants(bulk, exact);
            let mut creep = 0.0;
            let mut a = 0.0;
            for interval in 0..4 {
                let old = creep;
                // Prescribe the independent constant-Kirchhoff-stress creep
                // solution, then hold its final extension for relaxation.
                if interval < 3 {
                    a = creep + DT * 4.0 / ETA + 4.0 / young;
                }
                let increment = young * (a - creep) / (young + ETA / DT);
                creep += increment;
                let q = young * (a - creep);
                if interval < 3 {
                    near(q, 4.0, 2e-14);
                } else {
                    assert!(q > 0.0 && q < 4.0);
                }
                let step = accepted
                    .advance(
                        &supports(&nodes, a),
                        &[],
                        [0.0; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_uniaxial(&step.solid, a, creep, 0.0, bulk, exact);
                for state in step.solid.quadrature_states() {
                    near(state.alpha, creep, 3e-8);
                    assert!(state.alpha > old);
                    assert_eq!(state.creep_alpha, 0.0);
                    check_metric(state, creep);
                }
                let heat = step.solid.reference_volume_m3() * q * increment;
                // Sole Norton keeps its native plastic-named primary channel.
                near(step.energy.plastic_dissipation_j, heat, 3e-8);
                assert_eq!(step.energy.creep_dissipation_j, 0.0);
                near(step.energy.flow_potential_j, 0.5 * heat, 3e-8);
                assert!(step.energy.inelastic_dissipation_j() > 0.0);
                accepted = step.solid;
            }
        }
    });
}

#[test]
fn combined_creep_and_yield_retain_separate_histories_and_physical_heat() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        for (bulk, exact) in [(160.0, false), (8e8, false), (1e300, true)] {
            let material = FiniteNortonJ2::new(backbone(bulk), kinetics(), DT, 3.0).unwrap();
            let mut accepted = solid(&nodes, &elements, &material, exact, cx);
            let (young, _) = elastic_constants(bulk, exact);
            let (mut plastic, mut creep) = (0.0, 0.0);
            for (interval, a) in [0.06, 0.10, 0.10, 0.10].into_iter().enumerate() {
                let trial = young * (a - plastic - creep);
                let yield_old = YIELD + HARDENING * plastic;
                let relaxation = 1.0 + young * DT / ETA;
                let dp =
                    ((trial - relaxation * yield_old) / (young + relaxation * HARDENING)).max(0.0);
                let q = (trial - young * dp) / relaxation;
                let dc = DT * q / ETA;
                plastic += dp;
                creep += dc;
                let step = accepted
                    .advance(
                        &supports(&nodes, a),
                        &[],
                        [0.0; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_uniaxial(&step.solid, a, plastic + creep, plastic, bulk, exact);
                for state in step.solid.quadrature_states() {
                    near(state.alpha, plastic, 3e-8);
                    near(state.creep_alpha, creep, 3e-8);
                    check_metric(state, plastic + creep);
                }
                let v = step.solid.reference_volume_m3();
                near(step.energy.plastic_dissipation_j, v * YIELD * dp, 3e-8);
                near(step.energy.creep_dissipation_j, v * q * dc, 3e-8);
                near(
                    step.energy.inelastic_dissipation_j(),
                    v * (YIELD * dp + q * dc),
                    3e-8,
                );
                near(
                    step.energy.flow_potential_j,
                    v * (YIELD * dp + 0.5 * q * dc),
                    3e-8,
                );
                if interval < 2 {
                    assert!(dp > 0.0 && step.energy.plastic_dissipation_j > 0.0);
                } else {
                    assert_eq!(dp, 0.0);
                    near(step.energy.plastic_dissipation_j, 0.0, 3e-10);
                }
                assert!(dc > 0.0 && step.energy.creep_dissipation_j > 0.0);
                assert!(step.energy.equilibrium_work_residual_j.abs() < 2e-8);
                accepted = step.solid;
            }
        }
    });
}

#[test]
fn recovering_creep_keeps_kinetic_memory_distinct_from_permanent_shape() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        for (bulk, exact) in [(8e8, false), (1e300, true)] {
            let material = FiniteRecoveringNortonCreep::new(
                norton(bulk, DT)
                    .with_strain_hardening(StrainHardeningCreep {
                        reference_strain: 0.025,
                        exponent: 1.0,
                        residual_rate_fraction: 0.0,
                    })
                    .unwrap(),
                1.5,
            )
            .unwrap();
            let mut accepted = solid(&nodes, &elements, &material, exact, cx);
            let (young, _) = elastic_constants(bulk, exact);
            let (mut creep, mut resistance) = (0.0, 0.0);
            for interval in 0..3 {
                let a = if interval < 2 {
                    0.06
                } else {
                    creep + 0.02 / young
                };
                let old_resistance = resistance;
                let trial = young * (a - creep);
                let divisor = 1.0 + 1.5 * DT;
                let scale = 0.025 * divisor;
                let modulus = ETA / DT;
                // q=M*dc*(1+(r_old+dc)/scale), q+E*dc=q_trial.
                let quadratic = modulus / scale;
                let linear = young + modulus * (1.0 + resistance / scale);
                let dc =
                    2.0 * trial / (linear + (linear * linear + 4.0 * quadratic * trial).sqrt());
                let q = trial - young * dc;
                creep += dc;
                resistance = (resistance + dc) / divisor;
                let step = accepted
                    .advance(
                        &supports(&nodes, a),
                        &[],
                        [0.0; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_uniaxial(&step.solid, a, creep, 0.0, bulk, exact);
                for state in step.solid.quadrature_states() {
                    near(state.native.alpha, creep, 3e-8);
                    near(state.resistance_strain, resistance, 3e-8);
                    assert!(state.resistance_strain < state.native.alpha);
                    assert_eq!(state.native.creep_alpha, 0.0);
                    check_metric(&state.native, creep);
                }
                let v = step.solid.reference_volume_m3();
                near(step.energy.inelastic_dissipation_j(), v * q * dc, 3e-8);
                near(
                    step.energy.flow_potential_j,
                    v * modulus
                        * (0.5 * (1.0 + old_resistance / scale) * dc * dc
                            + dc.powi(3) / (3.0 * scale)),
                    3e-8,
                );
                if interval == 2 {
                    assert!(resistance < old_resistance && dc > 0.0);
                }
                accepted = step.solid;
            }
        }
    });
}

#[test]
fn clocks_and_failed_proposals_do_not_consume_a_material_interval() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let mismatched: Vec<_> = (0..64)
            .map(|q| norton(160.0, if q == 63 { 2.0 * DT } else { DT }))
            .collect();
        assert!(matches!(
            MixedPlasticTetSolid::try_new(
                &nodes,
                &elements,
                MixedPlasticLawField::PerQuadrature(&mismatched),
                VolumetricTreatment::OriginalMaterial,
                TetAssemblyBudget::standard(),
                4,
                64,
                cx,
            ),
            Err(MixedPlasticTetError::Native(
                MixedUpError::InvalidInput { .. }
            ))
        ));
        let material = FiniteOverstressJ2::new(backbone(160.0), kinetics(), DT).unwrap();
        let virgin = solid(&nodes, &elements, &material, false, cx);
        let accepted = virgin
            .advance(
                &supports(&nodes, 0.06),
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let history = accepted.quadrature_states().to_vec();
        let displacement = accepted.displacement_m().to_vec();
        let pressure = accepted.pressure_pa().to_vec();
        let deformation = accepted.deformation_gradients().to_vec();
        let stress = accepted.first_piola_pa().to_vec();
        let storage = accepted.stored_energy_j();
        let prescribed = supports(&nodes, 0.06);
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 0.0;
        assert!(matches!(
            accepted.advance(&prescribed, &[], [0.0; 3], None, settings(), strict, cx),
            Err(MixedPlasticTetError::EnergyGate { .. })
        ));
        let cancelled = CancelGate::new();
        cancelled.request();
        with_gate(&cancelled, |cx| {
            assert!(matches!(
                accepted.advance(&prescribed, &[], [0.0; 3], None, settings(), controls(), cx),
                Err(MixedPlasticTetError::Native(MixedUpError::Cancelled))
            ));
        });
        assert_eq!(accepted.quadrature_states(), history);
        assert_eq!(accepted.displacement_m(), displacement);
        assert_eq!(accepted.pressure_pa(), pressure);
        assert_eq!(accepted.deformation_gradients(), deformation);
        assert_eq!(accepted.first_piola_pa(), stress);
        assert_eq!(accepted.stored_energy_j(), storage);
        assert_eq!(accepted.duration_s(), Some(DT));
        let retry = accepted
            .advance(&prescribed, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        let repeat = accepted
            .advance(&prescribed, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        assert_eq!(
            retry.solid.quadrature_states(),
            repeat.solid.quadrature_states()
        );
        assert_eq!(retry.solid.displacement_m(), repeat.solid.displacement_m());
        assert_eq!(retry.solid.pressure_pa(), repeat.solid.pressure_pa());
        assert_eq!(retry.energy, repeat.energy);
        let (young, _) = elastic_constants(160.0, false);
        let old = history[0].alpha;
        let increment =
            (young * 0.06 - (young + HARDENING) * old - YIELD) / (young + HARDENING + ETA / DT);
        for state in retry.solid.quadrature_states() {
            near(state.alpha, old + increment, 3e-8);
        }
    });
}
