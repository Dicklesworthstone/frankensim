//! Native inertia/plasticity and bounded fatigue on actual accepted samples.
//! Synthetic coefficients are analytic implementation references, not test data.
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::{
    CriticalPlaneError, CriticalPlaneLimits, CycleClosureTolerance,
};
use fs_material::fatigue::strain_life::source::{
    CmbPowerBasis, CmbProperty, CmbSourceRequest, resolve_cmb,
};
use fs_material::fatigue::strain_life::{
    CmbParameters, PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainStressSample,
};
use fs_material::fatigue::{FatigueError, FatigueScope, RainflowLimits};
use fs_material::plastic::source::{HardeningPropertyBasis, ResolvedJ2Plasticity, j2_requirements};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::plastic3::dynamics::fatigue::{
    J2DynamicFatigueConfig, J2DynamicFatigueRun, J2DynamicFatigueTarget,
};
use fs_solid::plastic3::fatigue::{
    J2CriticalPlaneGauge, J2FatigueConfig, J2FatigueError, J2FatigueGauge, J2FatiguePath,
};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FRAME: ContentHash = ContentHash([19; 32]);

fn material(plastic: bool, density_dims: Option<Dims>) -> ResolvedJ2Plasticity {
    material_with_source(plastic, density_dims, None).0
}
fn material_with_source(
    plastic: bool,
    density_dims: Option<Dims>,
    calibration: Option<&CmbSourceRequest>,
) -> (ResolvedJ2Plasticity, Option<Arc<StrainLifeCurve>>) {
    let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let (temperature_domain, point) = if let Some(request) = calibration {
        (
            ValidityDomain::unconstrained().with_quantity(
                &request.temperature_axis,
                absolute_temperature,
                300.0,
                600.0,
            ),
            QueryPoint::new()
                .with_quantity(&request.temperature_axis, absolute_temperature, 300.0)
                .unwrap(),
        )
    } else {
        (
            ValidityDomain::unconstrained().with("T", 300.0, 600.0),
            QueryPoint::new().with("T", 300.0).unwrap(),
        )
    };
    let mut claims = ClaimSet::new();
    let mut values = vec![
        (
            "young_modulus",
            Pressure::DIMS,
            if plastic { 1000.0 } else { 60.0 },
        ),
        ("poisson_ratio", Dims::NONE, if plastic { 0.3 } else { 0.0 }),
        (
            "yield_stress",
            Pressure::DIMS,
            if plastic { 1.0 } else { 10.0 },
        ),
        ("isotropic_hardening_modulus", Pressure::DIMS, 20.0),
        ("kinematic_hardening_modulus", Pressure::DIMS, 80.0),
    ];
    if let Some(dims) = density_dims {
        values.push(("density", dims, 6.0));
    }
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: temperature_domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native dynamic fatigue fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    if let Some(request) = calibration {
        assert_eq!(request.power_basis, CmbPowerBasis::Reversals);
        for (property, value) in [
            (
                CmbProperty::FatigueStrengthCoefficient,
                if plastic { 10.0 } else { 0.6 },
            ),
            (CmbProperty::FatigueStrengthExponent, -0.5),
            (CmbProperty::FatigueDuctilityCoefficient, 0.04),
            (CmbProperty::FatigueDuctilityExponent, -0.5),
        ] {
            let key = property.property_key();
            claims
                .insert_claim(PropertyClaim {
                    value: PropertyValue::Scalar {
                        value,
                        dims: key.dims(),
                    },
                    key,
                    validity: temperature_domain
                        .clone()
                        .with_quantity(
                            &request.stress_ratio_axis,
                            QuantitySpec::dimensional(Dims::NONE),
                            -1.0,
                            -1.0,
                        )
                        .with_quantity(
                            request.power_basis.life_axis(),
                            QuantitySpec::dimensional(Dims::NONE),
                            1.0,
                            2e12,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic common CMB fit".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dynamic metal".into(),
            phase: "solid".into(),
            process: "synthetic dynamic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    let mut requirements = j2_requirements(basis).unwrap();
    if let Some(dims) = density_dims {
        requirements.push(
            ScalarPropertyRequirement::try_new(
                "density",
                dims,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        );
    }
    let resolved = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let curve = calibration.map(|request| {
        let fit_point = point
            .clone()
            .with_quantity(
                &request.stress_ratio_axis,
                QuantitySpec::dimensional(Dims::NONE),
                -1.0,
            )
            .unwrap()
            .with_quantity(
                request.power_basis.life_axis(),
                QuantitySpec::dimensional(Dims::NONE),
                100.0,
            )
            .unwrap();
        Arc::new(
            resolve_cmb(
                &card,
                &fit_point,
                request,
                [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 5],
            )
            .unwrap(),
        )
    });
    (
        ResolvedJ2Plasticity::try_from_material(resolved, basis, 0.05).unwrap(),
        curve,
    )
}

fn curve(law: &ResolvedJ2Plasticity, minimum: f64, maximum: f64) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic equal-exponent CMB reference; not measured".into(),
                material: "synthetic dynamic metal".into(),
                process: "synthetic dynamic coupon".into(),
                surface_and_geometry: "smooth fixed-plane coupon".into(),
                environment: "synthetic fixture".into(),
                load_basis:
                    "signed total mechanical strain amplitude; full cycles; named mean policy"
                        .into(),
                temperature_k: 300.0,
            },
            CmbParameters {
                youngs_modulus_pa: law.law().elastic.youngs,
                fatigue_strength_coefficient_pa: 0.01 * law.law().elastic.youngs,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.04,
                fatigue_ductility_exponent: -0.5,
            },
            minimum,
            maximum,
        )
        .unwrap(),
    )
}
fn gauge(
    law: &ResolvedJ2Plasticity,
    curve: Arc<StrainLifeCurve>,
    mean_stress: StrainLifeMeanStress,
) -> J2FatigueGauge {
    J2FatigueGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: law.material().identity(),
        curve,
        mean_stress,
    }
}
fn config(steps: usize) -> J2DynamicFatigueConfig {
    J2DynamicFatigueConfig {
        fatigue: J2FatigueConfig {
            history: RainflowLimits {
                maximum_samples: (steps + 1) as u64,
                maximum_residue: 64,
            },
            maximum_gauges: 4,
            maximum_total_residue: 256,
        },
        maximum_steps: steps,
    }
}
fn fixed(free_axial: bool) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&d| !free_axial || d != 3)
        .map(|d| (d, if d == 3 { 1.0 } else { 0.0 }))
        .collect()
}
fn problem<'a>(
    law: &'a ResolvedJ2Plasticity,
    fixed: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: fixed,
        nodal_forces_n: &[],
        load_path: path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 9,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}"
    );
}
fn same_state(a: &J2DynamicFatigueRun<'_, '_>, b: &J2DynamicFatigueRun<'_, '_>) {
    assert_eq!(a.solution(), b.solution());
    assert_eq!(a.velocity_mps(), b.velocity_mps());
    assert_eq!(a.time_s(), b.time_s());
    assert_eq!(a.accepted_steps(), b.accepted_steps());
    assert_eq!(a.last_step(), b.last_step());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(a.summaries(), b.summaries());
    assert_eq!(a.critical_plane_results(), b.critical_plane_results());
}

fn critical_gauge(
    law: &ResolvedJ2Plasticity,
    element: usize,
    normals: Vec<[f64; 3]>,
) -> J2CriticalPlaneGauge {
    J2CriticalPlaneGauge {
        element,
        material_state_identity: law.material().identity(),
        curve: curve(law, 0.5, 1e12),
        normals_in_mesh_frame: normals,
        limits: CriticalPlaneLimits {
            maximum_samples: 32,
            maximum_planes: 8,
            maximum_projection_work: 256,
        },
        closure: CycleClosureTolerance {
            absolute_strain: 1e-11,
            absolute_stress_pa: 1e-8,
            relative: 1e-10,
        },
    }
}

#[test]
fn g1_native_closed_tensor_cycle_selects_rotated_critical_plane_and_preserves_fixed_history() {
    let law = material(false, Some(Density::DIMS));
    let amplitude = 0.003;
    let loads = [amplitude, -amplitude, 0.0];
    // Actual affine native displacement gives eps=e*n*n^T with n=(1,1,0)/sqrt(2).
    // nu=0 gives sigma=E*eps. No tensor samples are supplied to the observer.
    let fixed: Vec<_> = NODES
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |c| (3 * node + c, if c < 2 { 0.5 * (x[0] + x[1]) } else { 0.0 }))
        })
        .collect();
    let p = problem(&law, &fixed, &loads);
    let dynamic_p = problem(&law, &fixed, &[1.0]);
    let coarse = critical_gauge(&law, 0, vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    let fine = critical_gauge(
        &law,
        0,
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
    );
    with_cx(&CancelGate::new(), |cx| {
        let fixed_gauge = gauge(
            &law,
            curve(&law, 0.5, 1e12),
            StrainLifeMeanStress::Uncorrected,
        );
        let mut native =
            J2FatiguePath::new(&p, FRAME, vec![fixed_gauge.clone()], config(3).fatigue, cx)
                .unwrap();
        let mut observed = J2FatiguePath::new_with_critical_planes(
            &p,
            FRAME,
            vec![fixed_gauge],
            vec![coarse.clone(), fine.clone()],
            config(3).fatigue,
            cx,
        )
        .unwrap();
        assert_eq!(observed.critical_plane_results(), vec![None, None]);
        observed.advance_cycle(3, cx).unwrap();
        native.advance(3, cx).unwrap();
        assert_eq!(observed.solution(), native.solution());
        assert_eq!(observed.summaries(), native.summaries());
        let results = observed.critical_plane_results();
        let (coarse_result, result) = (results[0].unwrap(), results[1].unwrap());
        assert_eq!(coarse_result.selected_plane, 0);
        assert_eq!(result.selected_plane, 2);
        assert_eq!(result.samples, 4);
        close(result.strain_amplitude, amplitude, 2e-12);
        close(result.maximum_normal_stress_pa, 60.0 * amplitude, 2e-10);
        close(result.minimum_normal_stress_pa, -60.0 * amplitude, 2e-10);
        // b=c=-1/2: N = sigma_f*(sigma_f/E+epsilon_f)/(2*sigma_max*epsilon_a).
        let life = 0.6 * 0.05 / (2.0 * (60.0 * amplitude) * amplitude);
        close(result.cycles_to_failure.unwrap(), life, 2e-9);
        close(result.damage_per_cycle, 1.0 / life, 2e-12);
        close(
            coarse_result.damage_per_cycle,
            result.damage_per_cycle / 4.0,
            2e-12,
        );
        assert_eq!(observed.finish(cx).unwrap(), native.finish(cx).unwrap());
        assert_eq!(
            observed.restart_critical_plane_cycle(cx),
            Err(J2FatigueError::Finished)
        );

        // Critical-only observation uses the same source binding and actual
        // inertial owner, with no artificial scalar/rainflow gauge required.
        let mut dynamic = J2DynamicFatigueRun::new_with_critical_planes(
            &dynamic_p,
            FRAME,
            vec![],
            vec![coarse, fine],
            config(3),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let mut whole = dynamic.clone();
        let targets: Vec<_> = loads
            .iter()
            .enumerate()
            .map(|(i, &load_factor)| J2DynamicFatigueTarget {
                time_s: 0.1 * (i + 1) as f64,
                load_factor,
            })
            .collect();
        assert!(matches!(
            dynamic.advance_cycle(&targets[..2], cx),
            Err(J2FatigueError::CriticalPlane {
                error: CriticalPlaneError::NotClosed { .. },
                ..
            })
        ));
        same_state(&dynamic, &whole);
        whole.advance_cycle(&targets, cx).unwrap();
        dynamic.advance(&targets[..2], cx).unwrap();
        assert_eq!(dynamic.critical_plane_results(), vec![None, None]);
        dynamic.advance_cycle(&targets[2..], cx).unwrap();
        same_state(&dynamic, &whole);
        assert_eq!(dynamic.critical_plane_results(), results);
        assert_eq!(
            dynamic.critical_plane_curve(1).unwrap().parameters(),
            curve(&law, 0.5, 1e12).parameters()
        );
        assert!(dynamic.finish(cx).unwrap().is_empty());
        assert_eq!(
            dynamic.restart_critical_plane_cycle(cx),
            Err(J2FatigueError::Finished)
        );
    });
}

#[test]
fn g3_actual_plastic_nonclosure_rolls_back_all_gauges_and_warmup_restart_preserves_history() {
    let laws = [
        material(false, Some(Density::DIMS)),
        material(true, Some(Density::DIMS)),
    ];
    let nodes: Vec<_> = NODES
        .into_iter()
        .chain(NODES.map(|mut x| {
            x[0] += 2.0;
            x
        }))
        .collect();
    let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
    // Opposite prescribed extension keeps the elastic bank in tension during
    // the later small elastic cycle about the plastic specimen's warm state.
    let fixed: Vec<_> = (0..24)
        .map(|d| {
            (
                d,
                if d == 3 {
                    -1.0
                } else if d == 15 {
                    1.0
                } else {
                    0.0
                },
            )
        })
        .collect();
    let base = -0.0005;
    let delta = 0.00002;
    let loads = [0.006, -0.006, 0.0, base, base + delta, base - delta, base];
    let p = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::PerElement(&laws),
        prescribed_m: &fixed,
        nodal_forces_n: &[],
        load_path: &loads,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let dynamic_p = SmallStrainTetProblem {
        load_path: &[1.0],
        ..p.clone()
    };
    let mut fixed_gauge = gauge(
        &laws[1],
        curve(&laws[1], 0.5, 1e12),
        StrainLifeMeanStress::Uncorrected,
    );
    fixed_gauge.element = 1;
    let banks = vec![
        critical_gauge(&laws[0], 0, vec![[1.0, 0.0, 0.0]]),
        critical_gauge(&laws[1], 1, vec![[1.0, 0.0, 0.0]]),
    ];
    with_cx(&CancelGate::new(), |cx| {
        let mut run = J2DynamicFatigueRun::new_with_critical_planes(
            &dynamic_p,
            FRAME,
            vec![fixed_gauge.clone()],
            banks.clone(),
            config(7),
            0.0,
            vec![[0.0; 3]; 8],
            cx,
        )
        .unwrap();
        let targets: Vec<_> = loads
            .iter()
            .enumerate()
            .map(|(i, &load_factor)| J2DynamicFatigueTarget {
                time_s: 0.1 * (i + 1) as f64,
                load_factor,
            })
            .collect();
        let initial = run.clone();
        // Gauge zero closes; the actual residual stress at plastic gauge one
        // prevents a representative cycle. No mechanics or fixed damage leaks.
        assert!(matches!(
            run.advance_cycle(&targets[..3], cx),
            Err(J2FatigueError::CriticalPlane {
                gauge: 1,
                error: CriticalPlaneError::NotClosed { .. },
            })
        ));
        same_state(&run, &initial);
        run.advance(&targets[..4], cx).unwrap();
        assert!(run.solution().element_states[1].alpha > 0.0);
        let warm = run.clone();
        assert!(matches!(
            run.finish(cx),
            Err(J2FatigueError::CriticalPlane { .. })
        ));
        same_state(&run, &warm);
        run.restart_critical_plane_cycle(cx).unwrap();
        same_state(&run, &warm);
        run.advance_cycle(&targets[4..], cx).unwrap();
        assert_eq!(
            run.solution().element_states,
            warm.solution().element_states
        );
        assert_eq!(run.summaries()[0].samples, 8);
        assert!(
            run.critical_plane_results()
                .iter()
                .all(|r| r.is_some_and(|r| r.samples == 4 && r.damage_per_cycle > 0.0))
        );

        let mut quasistatic = J2FatiguePath::new_with_critical_planes(
            &p,
            FRAME,
            vec![fixed_gauge],
            banks,
            config(7).fatigue,
            cx,
        )
        .unwrap();
        let initial = quasistatic.clone();
        assert!(matches!(
            quasistatic.advance_cycle(3, cx),
            Err(J2FatigueError::CriticalPlane { gauge: 1, .. })
        ));
        assert_eq!(quasistatic.solution(), initial.solution());
        assert_eq!(quasistatic.summaries(), initial.summaries());
        assert_eq!(
            quasistatic.critical_plane_results(),
            initial.critical_plane_results()
        );
        assert_eq!(quasistatic.accepted_increments(), 0);
        quasistatic.advance(4, cx).unwrap();
        let fixed_history = quasistatic.summaries();
        quasistatic.restart_critical_plane_cycle(cx).unwrap();
        assert_eq!(quasistatic.summaries(), fixed_history);
        quasistatic.advance_cycle(3, cx).unwrap();
        assert_eq!(quasistatic.samples(), run.samples());
        assert_eq!(
            quasistatic.critical_plane_results(),
            run.critical_plane_results()
        );
        assert_eq!(quasistatic.finish(cx).unwrap(), run.finish(cx).unwrap());
    });
}

#[test]
fn g1_unforced_vibration_counts_actual_fatigue_and_resolves_backward_euler_error() {
    let law = material(false, Some(Density::DIMS));
    let fixed = fixed(true);
    let p = problem(&law, &fixed, &[1.0]);
    with_cx(&CancelGate::new(), |cx| {
        let calibration = curve(&law, 0.5, 1e12);
        let mut previous_error = f64::INFINITY;
        for steps in [65, 130, 260] {
            let dt = 1.3 / steps as f64;
            let mut initial_v = vec![[0.0; 3]; 4];
            initial_v[1][0] = 0.01;
            let mut run = J2DynamicFatigueRun::new(
                &p,
                FRAME,
                vec![gauge(
                    &law,
                    calibration.clone(),
                    StrainLifeMeanStress::Uncorrected,
                )],
                config(steps),
                0.0,
                initial_v,
                cx,
            )
            .unwrap();
            let mut independent = StrainLifeMiner::new(
                calibration.clone(),
                StrainLifeMeanStress::Uncorrected,
                config(steps).fatigue.history,
            )
            .unwrap();
            independent.push(0.0).unwrap();
            // Unit tet, nu=0: K=E*V=10 N/m; M_11=rho*V/10=0.1 kg.
            // Thus omega=10/s. These scalar recurrences are independent of FE.
            let (mut u, mut v, mut time) = (0.0_f64, 0.01_f64, 0.0_f64);
            for index in 1..=steps {
                let next_time = index as f64 * dt;
                let interval = next_time - time;
                let next_u = (u + interval * v) / (1.0 + 100.0 * interval * interval);
                let next_v = (v - 100.0 * interval * u) / (1.0 + 100.0 * interval * interval);
                let report = run.advance_to(next_time, 0.0, cx).unwrap();
                close(run.solution().displacement_m[1][0], next_u, 2e-12);
                close(run.velocity_mps()[1][0], next_v, 2e-11);
                close(run.samples()[0].mechanical_strain, next_u, 2e-12);
                close(run.samples()[0].normal_stress_pa, 60.0 * next_u, 2e-10);
                close(report.energy.kinetic_before_j, 0.05 * v * v, 2e-14);
                close(report.energy.kinetic_after_j, 0.05 * next_v * next_v, 2e-14);
                close(
                    report.energy.numerical_loss_j,
                    0.05 * (next_v - v).powi(2),
                    2e-14,
                );
                // Both BE kinetic loss and endpoint elastic work cause the
                // numerical energy decay; neither is Miner heat or plasticity.
                let old_energy = 0.05 * v * v + 5.0 * u * u;
                let next_energy = 0.05 * next_v * next_v + 5.0 * next_u * next_u;
                close(
                    old_energy - next_energy,
                    report.energy.numerical_loss_j + 5.0 * (next_u - u).powi(2),
                    2e-14,
                );
                independent.push(next_u).unwrap();
                (u, v, time) = (next_u, next_v, next_time);
            }
            assert_eq!(run.solution().element_states[0].alpha, 0.0);
            let expected = independent.finish().unwrap();
            let result = run.finish(cx).unwrap()[0];
            assert_eq!(result.samples, (steps + 1) as u64);
            assert_eq!(result.full_cycles, expected.full_cycles);
            assert_eq!(result.half_cycles, expected.half_cycles);
            assert!(result.miner_damage > 0.0);
            close(result.miner_damage, expected.miner_damage, 2e-11);
            let time_error =
                (u - 0.001 * (10.0 * time).sin()).hypot((v - 0.01 * (10.0 * time).cos()) / 10.0);
            assert!(time_error < previous_error);
            previous_error = time_error;
            assert_eq!(run.finish(cx).unwrap()[0], result);
            assert_eq!(
                run.advance_to(time + dt, 0.0, cx),
                Err(J2FatigueError::Finished)
            );
        }
    });
}

#[test]
fn g1_cyclic_native_plastic_history_stress_and_inertial_support_work_are_preserved() {
    let law = material(true, Some(Density::DIMS));
    let fixed = fixed(false);
    let p = problem(&law, &fixed, &[1.0]);
    with_cx(&CancelGate::new(), |cx| {
        let calibration = curve(&law, 0.5, 1e12);
        let loads = [0.006, -0.006, 0.006, -0.006, 0.006];
        let corrected = [
            StrainLifeMeanStress::Morrow,
            StrainLifeMeanStress::SmithWatsonTopper,
        ];
        let observations = [
            StrainLifeMeanStress::Uncorrected,
            corrected[0],
            corrected[1],
        ]
        .into_iter()
        .map(|policy| gauge(&law, calibration.clone(), policy))
        .collect();
        let mut run = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            observations,
            config(loads.len()),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let mut paired: Vec<_> = corrected
            .into_iter()
            .map(|policy| {
                PairedStrainLifeMiner::new(
                    calibration.clone(),
                    policy,
                    config(loads.len()).fatigue.history,
                )
                .unwrap()
            })
            .collect();
        for counter in &mut paired {
            counter
                .push(StrainStressSample {
                    mechanical_strain: 0.0,
                    normal_stress_pa: 0.0,
                })
                .unwrap();
        }
        let mut independent = StrainLifeMiner::new(
            calibration,
            StrainLifeMeanStress::Uncorrected,
            config(loads.len()).fatigue.history,
        )
        .unwrap();
        independent.push(0.0).unwrap();
        let (lambda, mu) = law.law().elastic.lame();
        let (mut plastic, mut accumulated, mut old_e, mut old_v, mut old_time) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for (index, epsilon) in loads.into_iter().enumerate() {
            // Independent axisymmetric strain-controlled J2/Prager return:
            // eps=(e,0,0), eps_p=(p,-p/2,-p/2), q_trial=|2mu e-(3mu+Hkin)p|.
            let shifted_trial = 2.0 * mu * epsilon - (3.0 * mu + 80.0) * plastic;
            let excess = shifted_trial.abs() - (1.0 + 20.0 * accumulated);
            if excess > 0.0 {
                let delta = excess / (3.0 * mu + 100.0);
                plastic += shifted_trial.signum() * delta;
                accumulated += delta;
            }
            let stress = (lambda + 2.0 * mu) * epsilon - 2.0 * mu * plastic;
            let time = 0.1 * (index + 1) as f64;
            let interval = time - old_time;
            let velocity = (epsilon - old_e) / interval;
            let acceleration = (velocity - old_v) / interval;
            let report = run.advance_to(time, epsilon, cx).unwrap();
            close(run.samples()[0].mechanical_strain, epsilon, 1e-13);
            close(run.samples()[0].normal_stress_pa, stress, 2e-10);
            close(
                run.solution().element_states[0].plastic_strain[0],
                plastic,
                2e-12,
            );
            close(
                run.solution().element_states[0].plastic_strain[1],
                -0.5 * plastic,
                2e-12,
            );
            close(run.solution().element_states[0].alpha, accumulated, 2e-12);
            close(
                law.law().backstress_pa(&run.solution().element_states[0])[0],
                (2.0 / 3.0) * 80.0 * plastic,
                2e-11,
            );
            close(
                run.solution().increments[0].reactions_n[3],
                stress / 6.0 + 0.1 * acceleration,
                2e-10,
            );
            close(
                report.energy.kinetic_after_j,
                0.05 * velocity * velocity,
                2e-13,
            );
            close(
                report.energy.inertial_work_j,
                report.energy.kinetic_change_j + report.energy.numerical_loss_j,
                2e-13,
            );
            independent.push(epsilon).unwrap();
            for counter in &mut paired {
                counter
                    .push(StrainStressSample {
                        mechanical_strain: epsilon,
                        normal_stress_pa: stress,
                    })
                    .unwrap();
            }
            (old_e, old_v, old_time) = (epsilon, velocity, time);
        }
        assert!(accumulated > 0.01);
        assert_eq!(run.material(0).unwrap(), &law);
        assert_eq!(run.reference_density_kg_m3(), &[6.0]);
        let results = run.finish(cx).unwrap();
        assert_eq!(results[0], independent.finish().unwrap());
        for (actual, expected) in results[1..].iter().zip(&mut paired) {
            let expected = expected.finish().unwrap();
            assert_eq!(actual.samples, expected.samples);
            assert_eq!(actual.full_cycles, expected.full_cycles);
            assert_eq!(actual.half_cycles, expected.half_cycles);
            close(actual.miner_damage, expected.miner_damage, 2e-11);
        }
    });
}

#[test]
fn g1_named_mean_stress_models_use_same_actual_pairs_in_static_and_dynamic_owners() {
    let law = material(false, Some(Density::DIMS));
    let fixed = fixed(false);
    let loads = [0.003, 0.001, 0.003, 0.001, 0.003];
    let p = problem(&law, &fixed, &[1.0]);
    let static_p = problem(&law, &fixed, &loads);
    with_cx(&CancelGate::new(), |cx| {
        for policy in [
            StrainLifeMeanStress::Morrow,
            StrainLifeMeanStress::SmithWatsonTopper,
        ] {
            let calibration = curve(&law, 0.5, 1e12);
            let observation = gauge(&law, calibration.clone(), policy);
            let mut static_run = J2FatiguePath::new(
                &static_p,
                FRAME,
                vec![observation.clone()],
                config(loads.len()).fatigue,
                cx,
            )
            .unwrap();
            let mut dynamic = J2DynamicFatigueRun::new(
                &p,
                FRAME,
                vec![observation],
                config(loads.len()),
                0.0,
                vec![[0.0; 3]; 4],
                cx,
            )
            .unwrap();
            let mut independent = PairedStrainLifeMiner::new(
                calibration,
                policy,
                config(loads.len()).fatigue.history,
            )
            .unwrap();
            independent
                .push(StrainStressSample {
                    mechanical_strain: 0.0,
                    normal_stress_pa: 0.0,
                })
                .unwrap();
            let targets: Vec<_> = loads
                .iter()
                .enumerate()
                .map(|(i, &e)| J2DynamicFatigueTarget {
                    time_s: 0.1 * (i + 1) as f64,
                    load_factor: e,
                })
                .collect();
            dynamic.advance(&targets, cx).unwrap();
            static_run.advance(loads.len(), cx).unwrap();
            for e in loads {
                independent
                    .push(StrainStressSample {
                        mechanical_strain: e,
                        normal_stress_pa: 60.0 * e,
                    })
                    .unwrap();
            }
            let expected = independent.finish().unwrap();
            let dynamic_summary = dynamic.finish(cx).unwrap()[0];
            let static_summary = static_run.finish(cx).unwrap()[0];
            assert_eq!(dynamic_summary, static_summary);
            close(dynamic_summary.miner_damage, expected.miner_damage, 2e-12);
            assert_eq!(dynamic.samples(), static_run.samples());
            assert!(dynamic_summary.miner_damage > 0.0);
        }
    });
}

#[test]
fn g3_chunks_clones_and_late_gauge_refusal_preserve_whole_native_motion() {
    let law = material(true, Some(Density::DIMS));
    let fixed = fixed(false);
    let p = problem(&law, &fixed, &[1.0]);
    let loads = [0.006, -0.006, 0.006, -0.006, 0.006];
    let targets: Vec<_> = loads
        .iter()
        .enumerate()
        .map(|(i, &load_factor)| J2DynamicFatigueTarget {
            time_s: 0.1 * (i + 1) as f64,
            load_factor,
        })
        .collect();
    with_cx(&CancelGate::new(), |cx| {
        let mut whole = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            vec![gauge(
                &law,
                curve(&law, 0.5, 1e12),
                StrainLifeMeanStress::Uncorrected,
            )],
            config(8),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let mut chunks = whole.clone();
        whole.advance(&targets, cx).unwrap();
        chunks.advance(&targets[..2], cx).unwrap();
        let mut restored = chunks.clone();
        chunks.advance(&targets[2..], cx).unwrap();
        for target in &targets[2..] {
            restored
                .advance_to(target.time_s, target.load_factor, cx)
                .unwrap();
        }
        same_state(&whole, &chunks);
        same_state(&whole, &restored);
        assert_eq!(whole.finish(cx).unwrap(), chunks.finish(cx).unwrap());

        let mut guarded = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            vec![
                gauge(
                    &law,
                    curve(&law, 0.5, 1e12),
                    StrainLifeMeanStress::Uncorrected,
                ),
                gauge(
                    &law,
                    curve(&law, 200.0, 1e12),
                    StrainLifeMeanStress::Uncorrected,
                ),
            ],
            config(8),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        guarded.advance_to(0.05, 0.001, cx).unwrap();
        let before = guarded.clone();
        assert!(matches!(
            guarded.advance(&targets, cx),
            Err(J2FatigueError::Gauge { gauge: 1, .. })
        ));
        same_state(&guarded, &before);
        // Earlier accepted state is still usable after rejected plastic motion.
        guarded.advance_to(0.1, 0.0015, cx).unwrap();
        assert_eq!(guarded.accepted_steps(), 2);
    });
}

#[test]
fn g0_time_budget_cancellation_strain_and_compressive_swt_refusals_are_atomic() {
    let law = material(false, Some(Density::DIMS));
    let fixed = fixed(false);
    let p = problem(&law, &fixed, &[1.0]);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut run = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            vec![gauge(
                &law,
                curve(&law, 0.5, 1e12),
                StrainLifeMeanStress::SmithWatsonTopper,
            )],
            config(5),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let initial = run.clone();
        for time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                run.advance_to(time, 0.001, cx),
                Err(J2FatigueError::Invalid(_))
            ));
            same_state(&run, &initial);
        }
        assert!(matches!(
            run.advance_to(0.1, 0.051, cx),
            Err(J2FatigueError::Invalid(_))
        ));
        same_state(&run, &initial);
        run.advance_to(0.1, -0.001, cx).unwrap();
        let before_finish = run.clone();
        assert!(matches!(
            run.finish(cx),
            Err(J2FatigueError::Gauge {
                error: FatigueError::NoTensileStress { .. },
                ..
            })
        ));
        same_state(&run, &before_finish);
        assert!(matches!(
            run.advance(
                &[
                    J2DynamicFatigueTarget {
                        time_s: 0.2,
                        load_factor: -0.003
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.3,
                        load_factor: -0.001
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.4,
                        load_factor: -0.003
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.5,
                        load_factor: -0.001
                    },
                ],
                cx
            ),
            Err(J2FatigueError::Gauge {
                error: FatigueError::NoTensileStress { .. },
                ..
            })
        ));
        same_state(&run, &before_finish);
        assert_eq!(
            run.advance(
                &[J2DynamicFatigueTarget {
                    time_s: 1.0,
                    load_factor: 0.0
                }; 5],
                cx
            ),
            Err(J2FatigueError::Budget)
        );
        same_state(&run, &before_finish);
        gate.request();
        assert_eq!(
            run.advance_to(0.2, 0.001, cx),
            Err(J2FatigueError::Cancelled)
        );
        same_state(&run, &before_finish);
    });
}

#[test]
fn g0_original_density_identity_initial_support_velocity_and_path_are_required() {
    with_cx(&CancelGate::new(), |cx| {
        let fixed = fixed(false);
        for dims in [None, Some(Dims::NONE)] {
            let law = material(false, dims);
            let p = problem(&law, &fixed, &[1.0]);
            assert!(matches!(
                J2DynamicFatigueRun::new(
                    &p,
                    FRAME,
                    vec![gauge(
                        &law,
                        curve(&law, 0.5, 1e12),
                        StrainLifeMeanStress::Uncorrected
                    )],
                    config(4),
                    0.0,
                    vec![[0.0; 3]; 4],
                    cx
                ),
                Err(J2FatigueError::Invalid(_))
            ));
        }
        let law = material(false, Some(Density::DIMS));
        let p = problem(&law, &fixed, &[1.0]);
        let observation = gauge(
            &law,
            curve(&law, 0.5, 1e12),
            StrainLifeMeanStress::Uncorrected,
        );
        let mut bad_velocity = vec![[0.0; 3]; 4];
        bad_velocity[1][0] = 0.01;
        assert!(matches!(
            J2DynamicFatigueRun::new(
                &p,
                FRAME,
                vec![observation.clone()],
                config(4),
                0.0,
                bad_velocity,
                cx
            ),
            Err(J2FatigueError::Invalid(_))
        ));
        let mut wrong = observation.clone();
        wrong.material_state_identity = ContentHash([7; 32]);
        assert!(matches!(
            J2DynamicFatigueRun::new(
                &p,
                FRAME,
                vec![wrong],
                config(4),
                0.0,
                vec![[0.0; 3]; 4],
                cx
            ),
            Err(J2FatigueError::SourceMismatch { .. })
        ));
        let mut insufficient = config(4);
        insufficient.fatigue.history.maximum_samples = 4;
        assert!(matches!(
            J2DynamicFatigueRun::new(
                &p,
                FRAME,
                vec![observation.clone()],
                insufficient,
                0.0,
                vec![[0.0; 3]; 4],
                cx
            ),
            Err(J2FatigueError::Budget)
        ));
        let wrong_path = problem(&law, &fixed, &[0.0, 1.0]);
        assert!(matches!(
            J2DynamicFatigueRun::new(
                &wrong_path,
                FRAME,
                vec![observation],
                config(4),
                0.0,
                vec![[0.0; 3]; 4],
                cx
            ),
            Err(J2FatigueError::Invalid(_))
        ));

        // A source-declared typed absolute-temperature axis need not be named T.
        // The same original modulus claim must bind the CMB and actual J2 law.
        let request = CmbSourceRequest {
            power_basis: CmbPowerBasis::Reversals,
            temperature_axis: "specimen-temperature".into(),
            stress_ratio_axis: "calibration-stress-ratio".into(),
            surface_and_geometry: "smooth fixed-plane coupon".into(),
            environment: "synthetic fixture".into(),
            load_basis: "signed total mechanical strain amplitude; fully reversed calibration"
                .into(),
        };
        let (source_law, source_curve) =
            material_with_source(false, Some(Density::DIMS), Some(&request));
        assert_eq!(source_law.material().state_coordinate("T"), None);
        let source_curve = source_curve.unwrap();
        assert_eq!(
            source_curve.source().unwrap().request().temperature_axis,
            "specimen-temperature"
        );
        let source_gauge = gauge(
            &source_law,
            source_curve.clone(),
            StrainLifeMeanStress::Morrow,
        );
        let source_problem = problem(&source_law, &fixed, &[1.0]);
        let static_problem = problem(&source_law, &fixed, &[0.001]);
        let mut dynamic = J2DynamicFatigueRun::new(
            &source_problem,
            FRAME,
            vec![source_gauge.clone()],
            config(4),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let mut quasistatic = J2FatiguePath::new(
            &static_problem,
            FRAME,
            vec![source_gauge],
            config(4).fatigue,
            cx,
        )
        .unwrap();
        dynamic.advance_to(0.1, 0.001, cx).unwrap();
        quasistatic.advance(1, cx).unwrap();
        assert_eq!(dynamic.samples(), quasistatic.samples());
        let result = dynamic.finish(cx).unwrap();
        assert_eq!(result, quasistatic.finish(cx).unwrap());
        assert!(result[0].miner_damage > 0.0);

        // Critical-only banks bind the same immutable fit and typed temperature
        // receipts; they cannot substitute a curve from a different native state.
        let bank = J2CriticalPlaneGauge {
            curve: source_curve,
            ..critical_gauge(&source_law, 0, vec![[1.0, 0.0, 0.0]])
        };
        let mut wrong_bank = bank.clone();
        wrong_bank.material_state_identity = ContentHash([7; 32]);
        assert!(matches!(
            J2DynamicFatigueRun::new_with_critical_planes(
                &source_problem,
                FRAME,
                vec![],
                vec![wrong_bank],
                config(3),
                0.0,
                vec![[0.0; 3]; 4],
                cx
            ),
            Err(J2FatigueError::SourceMismatch { gauge: 0, .. })
        ));
        let mut critical = J2DynamicFatigueRun::new_with_critical_planes(
            &source_problem,
            FRAME,
            vec![],
            vec![bank],
            config(3),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        critical
            .advance_cycle(
                &[
                    J2DynamicFatigueTarget {
                        time_s: 0.1,
                        load_factor: 0.001,
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.2,
                        load_factor: -0.001,
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.3,
                        load_factor: 0.0,
                    },
                ],
                cx,
            )
            .unwrap();
        assert!(
            critical.critical_plane_results()[0]
                .unwrap()
                .damage_per_cycle
                > 0.0
        );
        assert_eq!(
            critical
                .critical_plane_curve(0)
                .unwrap()
                .source()
                .unwrap()
                .request()
                .temperature_axis,
            "specimen-temperature"
        );
        assert!(critical.finish(cx).unwrap().is_empty());
    });
}
