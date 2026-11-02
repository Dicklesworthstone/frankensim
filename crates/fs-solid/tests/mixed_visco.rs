//! Retained mixed histories against native material-point and traction oracles.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{
    BranchEnergy, FiniteStrainBranch, ViscoHyperelastic, ViscoRemainderState,
};
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::visco::{MixedViscoLawField, MixedViscoTetError, MixedViscoTetSolid};
use fs_solid::hyper3::visco::FiniteViscoStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn law(bulk: f64, tau: f64) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 1.0,
                c01: 0.0,
                kappa: bulk,
            },
            5.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 4.0 },
            tau_s: tau,
        }],
    )
    .unwrap()
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 6.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn controls() -> FiniteViscoStepControls {
    FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 2e-9,
        maximum_endpoint_work_defect_j: 0.5,
    }
}

fn mesh(count: usize, cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    let mut nodes = Vec::new();
    let mut elements = Vec::new();
    for e in 0..count {
        let start = nodes.len();
        for p in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]] {
            nodes.push([p[0] + 3.0 * e as f64, p[1], p[2]]);
        }
        elements.push([start, start + 1, start + 2, start + 3]);
    }
    promote_p1_to_p2(&nodes, &elements, None, cx).unwrap()
}

fn supports(nodes: &[[f64; 3]], held_stretch: Option<f64>) -> Vec<(usize, f64)> {
    let mut values = Vec::new();
    for (node, p) in nodes.iter().enumerate() {
        let local_x = p[0] - 3.0 * (p[0] / 3.0).floor();
        for axis in 0..3 {
            if axis == 0
                && let Some(stretch) = held_stretch
            {
                values.push((3 * node, (stretch - 1.0) * local_x));
            } else if (if axis == 0 { local_x } else { p[axis] }).abs() < 1e-12 {
                values.push((3 * node + axis, 0.0));
            }
        }
    }
    values
}

// The inclined face has A=sqrt(3)/2 and nominal traction P/sqrt(3).
// Its exact quadratic force is P/6 at each midside; corner weights vanish.
fn traction(element: &[usize; 10], nominal: f64) -> Vec<(usize, f64)> {
    TET10_EDGES
        .iter()
        .enumerate()
        .filter(|(_, edge)| !edge.contains(&0))
        .map(|(edge, _)| (3 * element[4 + edge], nominal / 6.0))
        .collect()
}

fn deformation(stretch: f64) -> [f64; 9] {
    [
        stretch,
        0.,
        0.,
        0.,
        stretch.sqrt().recip(),
        0.,
        0.,
        0.,
        stretch.sqrt().recip(),
    ]
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn match_point(solid: &MixedViscoTetSolid<'_>, expected: &ViscoRemainderState, pressure: f64) {
    for state in solid.quadrature_states() {
        for (a, b) in state
            .deformation_gradient
            .iter()
            .zip(expected.deformation_gradient)
        {
            near(*a, b, 3e-9);
        }
        for (a, b) in state
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(expected.inverse_viscous_cauchy_green.iter().flatten())
        {
            near(*a, *b, 4e-9);
        }
        near(state.stored, expected.stored, 3e-9);
        near(state.dissipated, expected.dissipated, 3e-9);
    }
    for p in solid.pressure_pa() {
        near(*p, pressure, 3e-8);
    }
}

#[test]
fn held_stretch_relaxes_native_histories_independently_of_unused_bulk_penalty() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(1, cx);
        let stretch = 1.2;
        let f = deformation(stretch);
        let held = supports(&nodes, Some(stretch));
        let mut baseline = None;
        for bulk in [1e3, 1e300] {
            let material = law(bulk, 0.2);
            let mut solid = MixedViscoTetSolid::try_new(
                &nodes,
                &elements,
                MixedViscoLawField::Uniform(&material),
                VolumetricTreatment::ExactIncompressible,
                TetAssemblyBudget::standard(),
                4,
                64,
                cx,
            )
            .unwrap();
            let mut point = material.initial_remainder_state();
            for dt in [0.0, 0.1, 0.1] {
                let response = material.step_remainder_response(&point, &f, dt).unwrap();
                let expected_p = -response.first_piola_pa[4] / stretch.sqrt();
                let step = solid
                    .advance(dt, &held, &[], [0.; 3], None, settings(), controls(), cx)
                    .unwrap();
                match_point(&step.solid, &response.state, expected_p);
                near(
                    step.energy.remainder_after_j,
                    response.state.stored / 6.0,
                    8e-10,
                );
                near(
                    step.energy.dissipated_j,
                    response.dissipated_pa / 6.0,
                    8e-10,
                );
                assert_eq!(step.energy.bulk_after_j, 0.0);
                if dt > 0.0 {
                    assert!(step.energy.dissipated_j > 0.0);
                    assert!(step.energy.stored_after_j < step.energy.stored_before_j);
                    near(step.energy.endpoint_work_defect_j, 0.0, 2e-9);
                }
                point = response.state;
                solid = step.solid;
            }
            let result = (
                solid.displacement_m().to_vec(),
                solid.pressure_pa().to_vec(),
                solid.quadrature_states().to_vec(),
                solid.stored_energy_j(),
            );
            if let Some(ref prior) = baseline {
                assert_eq!(&result, prior);
            }
            baseline = Some(result);
        }
    });
}

fn traction_point(
    material: &ViscoHyperelastic,
    old: &ViscoRemainderState,
    dt: f64,
    nominal: f64,
) -> ([f64; 9], f64, ViscoRemainderState) {
    let frozen = material.freeze_remainder(old, dt).unwrap();
    let mut bounds = [1.0, 2.0];
    for _ in 0..50 {
        let stretch = 0.5 * (bounds[0] + bounds[1]);
        let f = deformation(stretch);
        let (p, _) = frozen.remainder_stress_and_potential(&f).unwrap();
        let pressure = -p[4] / stretch.sqrt();
        let axial = p[0] + pressure / stretch;
        bounds[usize::from(axial > nominal)] = stretch;
    }
    let stretch = 0.5 * (bounds[0] + bounds[1]);
    let f = deformation(stretch);
    let result = material.step_remainder_response(old, &f, dt).unwrap();
    (f, -result.first_piola_pa[4] / stretch.sqrt(), result.state)
}

#[test]
fn maintained_traction_creeps_with_free_lateral_motion_and_independent_pressure() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(1, cx);
        let material = law(1e300, 0.2);
        let rollers = supports(&nodes, None);
        let force = traction(&elements[0], 0.6);
        let mut solid = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::Uniform(&material),
            VolumetricTreatment::ExactIncompressible,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let mut point = material.initial_remainder_state();
        let mut previous_stretch = 1.0;
        for dt in [0.0, 0.05, 0.05, 0.1] {
            let (f, pressure, expected) = traction_point(&material, &point, dt, 0.6);
            let step = solid
                .advance(
                    dt,
                    &rollers,
                    &force,
                    [0.; 3],
                    None,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            match_point(&step.solid, &expected, pressure);
            for stress in step.solid.first_piola_pa() {
                near(stress[0], 0.6, 2e-8);
                near(stress[4], 0.0, 2e-8);
                near(stress[8], 0.0, 2e-8);
            }
            assert!(f[0] > previous_stretch);
            near(step.equilibrium.min_det_f, 1.0, 3e-9);
            let reaction: f64 = step
                .equilibrium
                .reactions_n
                .iter()
                .filter(|(d, _)| d % 3 == 0)
                .map(|(_, force)| force)
                .sum();
            near(reaction, -0.3, 2e-9);
            if dt > 0.0 {
                assert!(step.energy.dissipated_j > 0.0);
            }
            previous_stretch = f[0];
            point = expected;
            solid = step.solid;
        }
        near(solid.time_s(), 0.2, 1e-15);
    });
}

#[test]
fn original_bulk_dual_has_physical_storage_and_elastic_limit_does_not_age() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(1, cx);
        let material = ViscoHyperelastic::new(law(20.0, 0.2).equilibrium, vec![]).unwrap();
        let rollers = supports(&nodes, None);
        let force = traction(&elements[0], 0.6);
        let initial = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::PerElement(core::slice::from_ref(&material)),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let loaded = initial
            .advance(
                0.0,
                &rollers,
                &force,
                [0.; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let p = loaded.solid.pressure_pa()[0];
        for pressure in loaded.solid.pressure_pa() {
            near(*pressure, p, 3e-9);
        }
        near(
            loaded.energy.bulk_after_j,
            p * p / (2.0 * 20.0 * 6.0),
            2e-10,
        );
        assert!(loaded.energy.bulk_after_j > 0.0);
        near(
            loaded.energy.stored_after_j,
            loaded.solid.quadrature_states()[0].stored / 6.0 + loaded.energy.bulk_after_j,
            2e-10,
        );
        near(
            loaded.energy.mixed_potential_j,
            loaded.energy.stored_after_j + loaded.energy.pressure_constraint_work_j,
            2e-10,
        );
        let held = loaded
            .solid
            .advance(
                0.7,
                &rollers,
                &force,
                [0.; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held.energy.dissipated_j, 0.0);
        near(
            held.energy.stored_after_j,
            loaded.energy.stored_after_j,
            2e-10,
        );
        for (a, b) in held
            .solid
            .displacement_m()
            .iter()
            .flatten()
            .zip(loaded.solid.displacement_m().iter().flatten())
        {
            near(*a, *b, 2e-9);
        }
    });
}

#[test]
fn heterogeneous_quadrature_spectra_and_rejected_endpoint_remain_transactional() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(2, cx);
        let materials: Vec<_> = [0.05, 0.5]
            .into_iter()
            .flat_map(|tau| core::iter::repeat_n(law(1e30, tau), 64))
            .collect();
        let held = supports(&nodes, Some(1.15));
        let initial = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::PerQuadrature(&materials),
            VolumetricTreatment::ExactIncompressible,
            TetAssemblyBudget::standard(),
            4,
            128,
            cx,
        )
        .unwrap();
        let jump = initial
            .advance(0.0, &held, &[], [0.; 3], None, settings(), controls(), cx)
            .unwrap();
        let relaxed = jump
            .solid
            .advance(0.1, &held, &[], [0.; 3], None, settings(), controls(), cx)
            .unwrap();
        for e in 0..2 {
            let expected = materials[64 * e]
                .step_remainder_response(
                    &jump.solid.quadrature_states()[64 * e],
                    &deformation(1.15),
                    0.1,
                )
                .unwrap();
            for state in &relaxed.solid.quadrature_states()[64 * e..64 * (e + 1)] {
                near(state.stored, expected.state.stored, 3e-9);
                near(state.dissipated, expected.state.dissipated, 3e-9);
            }
        }
        assert_ne!(
            relaxed.solid.quadrature_states()[0].inverse_viscous_cauchy_green,
            relaxed.solid.quadrature_states()[64].inverse_viscous_cauchy_green
        );
        let before = (
            relaxed.solid.displacement_m().to_vec(),
            relaxed.solid.pressure_pa().to_vec(),
            relaxed.solid.quadrature_states().to_vec(),
            relaxed.solid.time_s(),
        );
        let changed = supports(&nodes, Some(1.25));
        let strict = FiniteViscoStepControls {
            maximum_endpoint_work_defect_j: 0.0,
            ..controls()
        };
        let refusal =
            relaxed
                .solid
                .advance(0.02, &changed, &[], [0.; 3], None, settings(), strict, cx);
        assert!(
            matches!(
                refusal,
                Err(MixedViscoTetError::EnergyGate {
                    quantity: "signed endpoint integration defect",
                    ..
                })
            ),
            "{refusal:?}"
        );
        assert_eq!(
            before,
            (
                relaxed.solid.displacement_m().to_vec(),
                relaxed.solid.pressure_pa().to_vec(),
                relaxed.solid.quadrature_states().to_vec(),
                relaxed.solid.time_s()
            )
        );
        let retry = relaxed
            .solid
            .advance(
                0.02,
                &changed,
                &[],
                [0.; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        near(retry.solid.time_s(), 0.12, 1e-15);
        assert!(retry.energy.endpoint_work_defect_j.abs() > 1e-5);
    });
}
