//! Mixed finite plasticity against independent coaxial Hencky/J2 equations.
//! These manufactured coefficients verify numerics, not material calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
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

#[path = "mixed_plastic/bending.rs"]
mod bending;

fn with_gate(gate: &CancelGate, f: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 63,
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

fn law(bulk: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: bulk,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: HARDENING,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
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

fn supports(nodes: &[[f64; 3]], axial_log_stretch: Option<f64>) -> Vec<(usize, f64)> {
    let mut result = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        for axis in 0..3 {
            if axis == 0
                && let Some(a) = axial_log_stretch
            {
                result.push((3 * node, a.exp_m1() * x[0]));
            } else if x[axis] == 0.0 {
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

// Coaxial Hencky strains e_x=a-p, e_t=b+p/2 and zero lateral
// Kirchhoff stress give tau_x=E*(a-p), b=-nu*(a-p)-p/2.
// On a monotone plastic branch |tau_x|=yield+H*alpha. No native
// material update is used to compute these expected results.
fn check_uniaxial(
    solid: &MixedPlasticTetSolid<'_>,
    a: f64,
    plastic: f64,
    alpha: f64,
    bulk: f64,
    exact: bool,
) {
    let (young, poisson) = if exact {
        (3.0 * MU, 0.5)
    } else {
        (
            9.0 * MU / (3.0 + MU / bulk),
            (3.0 - 2.0 * MU / bulk) / (6.0 + 2.0 * MU / bulk),
        )
    };
    let b = -poisson * (a - plastic) - 0.5 * plastic;
    let f = [a.exp(), 0., 0., 0., b.exp(), 0., 0., 0., b.exp()];
    let tau = young * (a - plastic);
    let j = (a + 2.0 * b).exp();
    let pressure = tau / (3.0 * j);
    for gradient in solid.deformation_gradients() {
        for c in 0..9 {
            near(gradient[c], f[c], 3e-8);
        }
    }
    for p in solid.pressure_pa() {
        near(*p, pressure, 5e-7);
    }
    for stress in solid.first_piola_pa() {
        near(stress[0], tau / a.exp(), 5e-7);
        for value in &stress[1..] {
            near(*value, 0.0, 5e-7);
        }
    }
    let metric = [
        (-2.0 * plastic).exp(),
        plastic.exp(),
        plastic.exp(),
        0.,
        0.,
        0.,
    ];
    for state in solid.quadrature_states() {
        near(state.alpha, alpha, 3e-8);
        near(state.creep_alpha, 0.0, 0.0);
        for c in 0..6 {
            near(state.inverse_plastic_cauchy_green[c], metric[c], 3e-8);
        }
    }
    let v = solid.reference_volume_m3();
    let elastic = 0.5 * young * (a - plastic).powi(2);
    near(
        solid.stored_energy_j(),
        v * (elastic + 0.5 * HARDENING * alpha * alpha),
        2e-8,
    );
    if exact {
        assert_eq!(solid.bulk_energy_j(), 0.0);
    }
}

#[test]
fn original_bulk_and_exact_constraint_retain_plastic_unload_and_reverse_history() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        for (bulk, exact) in [(160.0, false), (8e8, false), (1e300, true)] {
            let material = law(bulk);
            let treatment = if exact {
                VolumetricTreatment::ExactIncompressible
            } else {
                VolumetricTreatment::OriginalMaterial
            };
            let mut solid = MixedPlasticTetSolid::try_new(
                &nodes,
                &elements,
                MixedPlasticLawField::Uniform(&material),
                treatment,
                TetAssemblyBudget::standard(),
                4,
                64,
                cx,
            )
            .unwrap();
            let young = if exact {
                3.0 * MU
            } else {
                9.0 * MU / (3.0 + MU / bulk)
            };
            let mut plastic = 0.0;
            for a in [0.002, 0.04, 0.08] {
                let next_plastic = ((young * a - YIELD) / (young + HARDENING)).max(0.0);
                let step = solid
                    .advance(
                        &supports(&nodes, Some(a)),
                        &[],
                        [0.0; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_uniaxial(&step.solid, a, next_plastic, next_plastic, bulk, exact);
                near(
                    step.energy.plastic_dissipation_j,
                    solid.reference_volume_m3() * YIELD * (next_plastic - plastic),
                    2e-8,
                );
                near(
                    step.energy.flow_potential_j,
                    step.energy.plastic_dissipation_j,
                    2e-8,
                );
                near(step.energy.equilibrium_work_residual_j, 0.0, 2e-8);
                solid = step.solid;
                plastic = next_plastic;
            }
            let released = solid
                .advance(
                    &supports(&nodes, None),
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_uniaxial(&released.solid, plastic, plastic, plastic, bulk, exact);
            near(released.energy.plastic_dissipation_j, 0.0, 2e-10);
            // The accepted loading endpoint remains available to fork/retry.
            check_uniaxial(&solid, 0.08, plastic, plastic, bulk, exact);
            let a = -0.03;
            let reverse_increment =
                (-young * (a - plastic) - YIELD - HARDENING * plastic) / (young + HARDENING);
            assert!(reverse_increment > 0.0);
            let reverse = released
                .solid
                .advance(
                    &supports(&nodes, Some(a)),
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_uniaxial(
                &reverse.solid,
                a,
                plastic - reverse_increment,
                plastic + reverse_increment,
                bulk,
                exact,
            );
            near(
                reverse.energy.plastic_dissipation_j,
                solid.reference_volume_m3() * YIELD * reverse_increment,
                3e-8,
            );
        }
    });
}

#[test]
fn quadrature_bindings_preserve_distinct_hardening_histories() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let base = law(160.0);
        let stronger = FiniteJ2::new(
            FiniteJ2Parameters {
                yield_stress_pa: 2.0 * YIELD,
                ..*base.parameters()
            },
            *base.limits(),
        )
        .unwrap();
        let laws: Vec<_> = (0..64)
            .map(|q| {
                if q % 2 == 0 {
                    base.clone()
                } else {
                    stronger.clone()
                }
            })
            .collect();
        let solid = MixedPlasticTetSolid::try_new(
            &nodes,
            &elements,
            MixedPlasticLawField::PerQuadrature(&laws),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let a: f64 = 0.04;
        let stretch = [a.exp(), (-0.5 * a).exp(), (-0.5 * a).exp()];
        let prescribed: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, x)| (0..3).map(move |c| (3 * node + c, (stretch[c] - 1.0) * x[c])))
            .collect();
        let step = solid
            .advance(&prescribed, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        assert_eq!(step.solid.quadrature_states().len(), 64);
        for (q, state) in step.solid.quadrature_states().iter().enumerate() {
            let yield_radius = if q % 2 == 0 { YIELD } else { 2.0 * YIELD };
            let expected = (3.0 * MU * a - yield_radius) / (3.0 * MU + HARDENING);
            near(state.alpha, expected, 2e-12);
        }
        assert!(solid.quadrature_states().iter().all(|s| s.alpha == 0.0));
        assert!(
            MixedPlasticTetSolid::try_new(
                &nodes,
                &elements,
                MixedPlasticLawField::PerQuadrature(&laws[..63]),
                VolumetricTreatment::OriginalMaterial,
                TetAssemblyBudget::standard(),
                4,
                64,
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn refusal_and_cancellation_preserve_the_entire_accepted_endpoint() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(160.0);
        let virgin = MixedPlasticTetSolid::try_new(
            &nodes,
            &elements,
            MixedPlasticLawField::Uniform(&material),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let accepted = virgin
            .advance(
                &supports(&nodes, Some(0.04)),
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let histories = accepted.quadrature_states().to_vec();
        let displacements = accepted.displacement_m().to_vec();
        let pressures = accepted.pressure_pa().to_vec();
        let energy = accepted.stored_energy_j();
        let mut invalid_settings = settings();
        invalid_settings.load_steps = 2;
        assert!(
            accepted
                .advance(
                    &supports(&nodes, Some(0.08)),
                    &[],
                    [0.0; 3],
                    None,
                    invalid_settings,
                    controls(),
                    cx
                )
                .is_err()
        );
        let mut no_work_error = controls();
        no_work_error.maximum_endpoint_work_defect_j = 0.0;
        assert!(matches!(
            accepted.advance(
                &supports(&nodes, Some(0.08)),
                &[],
                [0.0; 3],
                None,
                settings(),
                no_work_error,
                cx
            ),
            Err(MixedPlasticTetError::EnergyGate { .. })
        ));
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(matches!(
                accepted.advance(
                    &supports(&nodes, Some(0.08)),
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls(),
                    cancelled
                ),
                Err(MixedPlasticTetError::Native(MixedUpError::Cancelled))
            ));
        });
        assert_eq!(accepted.quadrature_states(), histories);
        assert_eq!(accepted.displacement_m(), displacements);
        assert_eq!(accepted.pressure_pa(), pressures);
        assert_eq!(accepted.stored_energy_j(), energy);
        let retry = accepted
            .advance(
                &supports(&nodes, Some(0.08)),
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let repeat = accepted
            .advance(
                &supports(&nodes, Some(0.08)),
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            retry.solid.quadrature_states(),
            repeat.solid.quadrature_states()
        );
        assert_eq!(retry.solid.displacement_m(), repeat.solid.displacement_m());
        assert_eq!(retry.solid.pressure_pa(), repeat.solid.pressure_pa());
        assert_eq!(retry.energy, repeat.energy);
    });
}
