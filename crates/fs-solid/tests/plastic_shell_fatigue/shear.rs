//! Shear life from original shell coordinates and accepted material histories.
//! Synthetic coefficients verify implementation, not measured component life.

use super::*;
use fs_material::fatigue::critical_plane::fatemi_socie::source::{
    FatemiSociePowerBasis, FatemiSocieProperty, FatemiSocieSourceRequest,
    resolve_isotropic_fatemi_socie,
};
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieCycle, FatemiSocieError, FatemiSocieLimits,
    FatemiSocieNormalStress,
};
use fs_material::fatigue::critical_plane::{CycleClosureTolerance, TensorFatigueSample};
use fs_material::plastic::chaboche::source::{BackstressProperties, resolve_chaboche_branches};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY, resolve_voce_j2,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solid::plastic_shell::SmallStrainShellSolution;
use fs_solid::plastic_shell::fatigue::{ShellFatemiSocieConfig, ShellFatemiSocieGauge};

const G: f64 = E / (2.0 * (1.0 + NU));
const GAMMA: f64 = 0.001;

pub(super) struct ShearFixture {
    pub(super) law: ResolvedJ2Plasticity,
    pub(super) calibration: Arc<FatemiSocieCalibration>,
    pub(super) card: MaterialCard,
}

fn source_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

pub(super) fn fixture(revision: u32) -> ShearFixture {
    let base = card_with_properties(
        revision,
        source_temperature(),
        &[
            (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.8),
            (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 120.0),
            ("fast_backstress_modulus", Pressure::DIMS, 80.0),
            ("fast_backstress_recovery", Dims::NONE, 40.0),
            ("slow_backstress_modulus", Pressure::DIMS, 40.0),
            ("slow_backstress_recovery", Dims::NONE, 10.0),
        ],
    );
    let request = FatemiSocieSourceRequest {
        power_basis: FatemiSociePowerBasis::Reversals,
        normal_stress: FatemiSocieNormalStress::TensionOnly,
        temperature_axis: "T".into(),
        shear_stress_ratio_axis: "torsion-R".into(),
        surface_and_geometry: "synthetic smooth shell torsion coupon".into(),
        environment: "synthetic fixed environment".into(),
        load_basis: "engineering shear strain amplitude".into(),
    };
    let scalar = QuantitySpec::dimensional(Dims::NONE);
    let domain = ValidityDomain::unconstrained()
        .with_quantity("T", source_temperature(), 300.0, 300.0)
        .with_quantity(&request.shear_stress_ratio_axis, scalar, -1.0, -1.0)
        .with_quantity(request.power_basis.life_axis(), scalar, 2.0, 2e8);
    let mut claims = base.claims().clone();
    for (property, value) in [
        (
            FatemiSocieProperty::ShearFatigueStrengthCoefficient,
            0.1 * G,
        ),
        (FatemiSocieProperty::ShearFatigueStrengthExponent, -0.5),
        (FatemiSocieProperty::ShearFatigueDuctilityCoefficient, 0.1),
        (FatemiSocieProperty::ShearFatigueDuctilityExponent, -0.5),
        (FatemiSocieProperty::NormalStressSensitivity, 0.2),
    ] {
        let key = property.property_key(request.normal_stress);
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity: domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "one synthetic torsional fatigue fit".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::supersede(&base, claims, vec![]).unwrap();
    let law = material(&card, source_temperature());
    let query = point(source_temperature())
        .with_quantity(&request.shear_stress_ratio_axis, scalar, -1.0)
        .unwrap()
        .with_quantity(request.power_basis.life_axis(), scalar, 100.0)
        .unwrap();
    let calibration = Arc::new(
        resolve_isotropic_fatemi_socie(
            &card,
            &query,
            &request,
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 8],
        )
        .unwrap(),
    );
    ShearFixture {
        law,
        calibration,
        card,
    }
}

pub(super) fn fs_gauge<L: SourceBoundFatigueLaw>(
    law: &L,
    calibration: &Arc<FatemiSocieCalibration>,
) -> ShellFatemiSocieGauge {
    ShellFatemiSocieGauge {
        element: 0,
        surface_point: 0,
        thickness_point: 1,
        material_state_identity: law.fatigue_source().identity(),
        calibration: calibration.clone(),
        normals_in_facet_frame: vec![[1.0, 0.0, 0.0]],
        limits: FatemiSocieLimits {
            maximum_samples: 16,
            maximum_planes: 4,
            maximum_projection_work: 64,
            maximum_circle_work: 4096,
        },
        closure: CycleClosureTolerance {
            absolute_strain: 1e-12,
            absolute_stress_pa: 1e-10,
            relative: 1e-10,
        },
    }
}

pub(super) fn shear_config() -> ShellFatemiSocieConfig {
    ShellFatemiSocieConfig {
        maximum_total_samples: 128,
        maximum_total_projection_work: 512,
        maximum_total_circle_work: 32768,
    }
}

/// A symmetric displacement gradient gives eps_xy=gamma/2 with no rigid spin.
pub(super) fn shear_field(mesh: &ShellMesh, gamma: f64) -> Vec<(usize, f64)> {
    let frame = mesh.facet(0).unwrap().frame;
    let origin = mesh.nodes[0];
    let mut prescribed = Vec::new();
    for (node, point) in mesh.nodes.iter().enumerate() {
        let x: f64 = (0..3).map(|c| frame[0][c] * (point[c] - origin[c])).sum();
        let y: f64 = (0..3).map(|c| frame[1][c] * (point[c] - origin[c])).sum();
        for c in 0..6 {
            let value = if c < 3 {
                0.5 * gamma * (frame[0][c] * y + frame[1][c] * x)
            } else {
                0.0
            };
            prescribed.push((6 * node + c, value));
        }
    }
    prescribed
}

fn normal_gauge<L: SourceBoundFatigueLaw>(law: &L) -> ShellFatigueGauge {
    let mut g = gauge(law, 1, StrainLifeMeanStress::Uncorrected);
    g.normals_in_facet_frame = vec![[1.0, 1.0, 0.0]];
    g
}

fn native_sample<L: SourceBoundFatigueLaw>(
    law: &L,
    solution: &SmallStrainShellSolution<L::State>,
) -> TensorFatigueSample {
    let point = &solution.section_states[0].points[1];
    law.accepted_fatigue_sample(&point.strain, &point.history)
        .unwrap()
}

fn same<L: SourceBoundFatigueLaw>(
    a: &ShellFatiguePath<'_, '_, L>,
    b: &ShellFatiguePath<'_, '_, L>,
) {
    assert_eq!(a.solution(), b.solution());
    assert_eq!(a.velocity(), b.velocity());
    assert_eq!(a.elapsed_time_s(), b.elapsed_time_s());
    assert_eq!(a.next_step(), b.next_step());
    assert_eq!(a.sample_count(), b.sample_count());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(a.plane_summaries(0), b.plane_summaries(0));
    assert_eq!(a.results(), b.results());
    assert_eq!(a.fatemi_socie_samples(), b.fatemi_socie_samples());
    assert_eq!(a.fatemi_socie_results(), b.fatemi_socie_results());
    assert_eq!(
        a.fatemi_socie_miner_summaries(),
        b.fatemi_socie_miner_summaries()
    );
    assert_eq!(a.is_finished(), b.is_finished());
}

#[test]
fn g1_shell_shear_miner_keeps_cycle_totals_and_rolls_back_exhausted_cycle_exposure() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = shear_field(&mesh, GAMMA);
        let factors = [1.0, -1.0, 0.0, 0.5, -0.5, 0.0, 1.0, -1.0, 0.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut g = fs_gauge(&f.law, &f.calibration);
        g.normals_in_facet_frame = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let mut run = ShellFatiguePath::new_with_fatemi_socie(
            &p,
            vec![normal_gauge(&f.law)],
            vec![g],
            config(),
            shear_config(),
            cx,
        )
        .unwrap();
        assert_eq!(run.fatemi_socie_miner_summaries(), vec![None]);
        assert!(run.enable_fatemi_socie_miner(0, cx).is_err());
        assert_eq!(run.fatemi_socie_miner_summaries(), vec![None]);
        run.enable_fatemi_socie_miner(2, cx).unwrap();
        run.advance_next(cx).unwrap();
        run.advance_next(cx).unwrap();
        run.advance_next_cycle(cx).unwrap();
        let first = run.clone();
        let result = run.fatemi_socie_results()[0].unwrap();
        assert_eq!(run.close_fatemi_socie_cycle(cx).unwrap(), vec![result]);
        same(&run, &first);
        close(
            run.fatemi_socie_miner_summaries()[0].unwrap().miner_damage,
            1.0 / 20000.0,
            1e-12,
        );

        let mut resumed = run.clone();
        for branch in [&mut run, &mut resumed] {
            branch.restart_fatemi_socie_cycle(cx).unwrap();
            assert_eq!(
                branch.fatemi_socie_miner_summaries(),
                first.fatemi_socie_miner_summaries()
            );
            assert!(branch.enable_fatemi_socie_miner(100, cx).is_err());
            branch.advance_next(cx).unwrap();
            branch.advance_next(cx).unwrap();
            branch.advance_next_cycle(cx).unwrap();
        }
        same(&run, &resumed);
        let summary = run.fatemi_socie_miner_summaries()[0].unwrap();
        assert_eq!(summary.completed_cycles, 2);
        // Halving engineering shear amplitude quadruples this synthetic life.
        close(summary.miner_damage, 1.0 / 20000.0 + 1.0 / 80000.0, 1e-12);
        let planes = run.fatemi_socie_plane_damage(0).unwrap();
        close(planes[0].miner_damage, summary.miner_damage, 1e-14);
        assert_eq!(planes[1].miner_damage, 0.0);

        run.restart_fatemi_socie_cycle(cx).unwrap();
        run.advance_next(cx).unwrap();
        run.advance_next(cx).unwrap();
        let prefix = run.clone();
        assert!(matches!(
            run.advance_next_cycle(cx),
            Err(ShellFatigueError::FatemiSocie {
                error: FatemiSocieError::Budget,
                ..
            })
        ));
        same(&run, &prefix);
        gate.request();
        assert_eq!(
            run.advance_next_cycle(cx),
            Err(ShellFatigueError::Cancelled)
        );
        assert_eq!(
            run.restart_fatemi_socie_cycle(cx),
            Err(ShellFatigueError::Cancelled)
        );
        same(&run, &prefix);
    });
}

#[test]
fn elastic_shell_shear_matches_engineering_amplitude_and_analytic_life_with_normal_bank() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = shear_field(&mesh, GAMMA);
        let factors = [1.0, -1.0, 0.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut g = fs_gauge(&f.law, &f.calibration);
        g.normals_in_facet_frame = vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut run = ShellFatiguePath::new_with_fatemi_socie(
            &p,
            vec![normal_gauge(&f.law)],
            vec![g],
            config(),
            shear_config(),
            cx,
        )
        .unwrap();
        let mut direct = StrainLifeMiner::new(
            curve(0.5),
            StrainLifeMeanStress::Uncorrected,
            limits().rainflow,
        )
        .unwrap();
        direct.push(0.0).unwrap();
        for (index, factor) in factors.into_iter().enumerate() {
            if index == 2 {
                run.advance_next_cycle(cx).unwrap();
            } else {
                run.advance_next(cx).unwrap();
            }
            let sample = run.fatemi_socie_samples()[0];
            close(
                sample.total_mechanical_strain[3],
                0.5 * GAMMA * factor,
                1e-14,
            );
            close(sample.stress_pa[3], G * GAMMA * factor, 1e-10);
            assert_eq!(sample, run.samples()[0]);
            assert_eq!(
                run.solution().section_states[0].points[1].history.alpha,
                0.0
            );
            direct.push(0.5 * GAMMA * factor).unwrap();
        }
        let result = run.fatemi_socie_results()[0].unwrap();
        assert_eq!(result.samples, 4);
        assert_eq!(result.selected_plane, 0);
        close(result.engineering_shear_strain_amplitude, GAMMA, 1e-14);
        close(result.normal_stress_correction, 1.0, 1e-14);
        // gamma_a=.2/sqrt(2N), so tensor shear mistaken for engineering shear
        // would overestimate life by four: 80000 instead of 20000 cycles.
        close(
            result.cycles_to_failure.unwrap(),
            0.5 * (0.2 / GAMMA).powi(2),
            1e-6,
        );
        assert_eq!(run.fatemi_socie_material_source(0), Some(f.law.material()));
        assert_eq!(
            run.fatemi_socie_calibration(0),
            Some(f.calibration.as_ref())
        );
        let closed = run.clone();
        assert_eq!(run.close_fatemi_socie_cycle(cx).unwrap(), vec![result]);
        same(&run, &closed);
        close(
            run.finish(cx).unwrap()[0].summary.miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-13,
        );
        assert_eq!(run.fatemi_socie_results(), vec![Some(result)]);
    });
}

fn conditioned_cycle<L: SourceBoundFatigueLaw>(
    law: &L,
    calibration: &Arc<FatemiSocieCalibration>,
    cx: &Cx<'_>,
) -> L::State {
    let section = section(law.clone());
    let mesh = mesh();
    let fixed = shear_field(&mesh, 0.02);
    let factors = [1.0, 0.9, 1.0, 0.9, 1.0];
    let p = problem(&mesh, &section, &fixed, &[], &factors);
    let g = fs_gauge(law, calibration);
    let mut run = ShellFatiguePath::new_with_fatemi_socie(
        &p,
        vec![normal_gauge(law)],
        vec![g.clone()],
        config(),
        shear_config(),
        cx,
    )
    .unwrap();
    run.advance_next(cx).unwrap();
    let mut native = problem(&mesh, &section, &fixed, &[], &[1.0])
        .solve(cx)
        .unwrap();
    assert_eq!(run.solution(), &native);
    assert!(law.accepted_plastic_strain(&native.section_states[0].points[1].history)[3] > 0.001);
    // The real plastic preload is retained. These small subsequent loops are
    // elastic; a transient hardening loop is not labelled cyclically stabilized.
    for _ in 0..2 {
        let before = run.clone();
        run.restart_fatemi_socie_cycle(cx).unwrap();
        assert_eq!(run.solution(), before.solution());
        assert_eq!(run.plane_summaries(0), before.plane_summaries(0));
        assert_eq!(run.sample_count(), before.sample_count());
        let mut direct = FatemiSocieCycle::new(
            g.calibration.clone(),
            g.normals_in_facet_frame.clone(),
            g.limits,
            g.closure,
        )
        .unwrap();
        direct.push(native_sample(law, &native)).unwrap();
        for factor in [0.9, 1.0] {
            native = problem(&mesh, &section, &fixed, &[], &[factor])
                .solve_from(&native.coordinates, &native.section_states, cx)
                .unwrap();
            direct.push(native_sample(law, &native)).unwrap();
            if factor == 1.0 {
                run.advance_next_cycle(cx).unwrap();
            } else {
                run.advance_next(cx).unwrap();
            }
            assert_eq!(run.solution(), &native);
            assert_eq!(run.fatemi_socie_samples()[0], native_sample(law, &native));
        }
        assert_eq!(
            run.fatemi_socie_results(),
            vec![Some(direct.finish().unwrap())]
        );
    }
    assert_eq!(run.sample_count(), 6);
    run.finish(cx).unwrap();
    native.section_states[0].points[1].history.clone()
}

#[test]
fn conditioned_voce_and_chaboche_windows_keep_complete_native_history_and_restart() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let query = point(source_temperature());
        let basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
        let voce = resolve_voce_j2(
            &f.card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            basis,
            0.05,
        )
        .unwrap();
        assert!(conditioned_cycle(&voce, &f.calibration, cx).alpha > 0.0);
        let chaboche = resolve_chaboche_branches(
            &f.card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            basis,
            0.05,
            [
                BackstressProperties::new("fast_backstress_modulus", "fast_backstress_recovery"),
                BackstressProperties::new("slow_backstress_modulus", "slow_backstress_recovery"),
            ],
        )
        .unwrap();
        let history = conditioned_cycle(&chaboche, &f.calibration, cx);
        assert!(history.plastic.alpha > 0.0);
        assert!(history.backstress_pa.iter().all(|branch| branch[3] > 0.01));
        assert_ne!(history.backstress_pa[0], history.backstress_pa[1]);
    });
}

#[test]
fn shear_cycle_backward_euler_keeps_native_motion_energy_and_clone_continuation() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = shear_field(&mesh, GAMMA);
        let factors = [1.0, -1.0, 0.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1.0), cx)
            .unwrap();
        let mut run = ShellFatiguePath::new_with_fatemi_socie(
            &p,
            vec![normal_gauge(&f.law)],
            vec![fs_gauge(&f.law, &f.calibration)],
            config(),
            shear_config(),
            cx,
        )
        .unwrap();
        run.enable_fatemi_socie_miner(2, cx).unwrap();
        let mut native = problem(&mesh, &section, &fixed, &[], &[0.0])
            .solve(cx)
            .unwrap();
        let mut velocity = vec![[0.0; 6]; 3];
        let dt = 0.01;
        let mut resumed = None;
        for (index, factor) in factors.into_iter().enumerate() {
            let prescribed: Vec<_> = fixed.iter().map(|&(d, v)| (d, factor * v)).collect();
            let endpoint = problem(&mesh, &section, &prescribed, &[], &[1.0]);
            let mut step = endpoint
                .advance_backward_euler(Some(&native), &velocity, &reference, dt, cx)
                .unwrap();
            let energy = if index == 2 {
                run.advance_next_backward_euler_cycle(&reference, dt, cx)
            } else {
                run.advance_next_backward_euler(&reference, dt, cx)
            }
            .unwrap()
            .unwrap();
            assert_eq!(energy, step.inertia);
            step.solution.increments[0].load_factor = factor;
            assert_eq!(run.solution(), &step.solution);
            assert_eq!(run.velocity(), step.velocity);
            assert_eq!(
                run.fatemi_socie_samples()[0],
                native_sample(&f.law, &step.solution)
            );
            native = step.solution;
            velocity = step.velocity;
            if index == 0 {
                resumed = Some(run.clone());
            } else if let Some(restarted) = &mut resumed {
                let copied = if index == 2 {
                    restarted.advance_next_backward_euler_cycle(&reference, dt, cx)
                } else {
                    restarted.advance_next_backward_euler(&reference, dt, cx)
                }
                .unwrap();
                assert_eq!(copied, Some(energy));
                same(&run, restarted);
            }
        }
        close(run.elapsed_time_s(), 3.0 * dt, 1e-15);
        let damage = run.fatemi_socie_miner_summaries()[0].unwrap();
        assert_eq!(damage.completed_cycles, 1);
        close(damage.miner_damage, 1.0 / 20000.0, 1e-12);
        close(
            run.fatemi_socie_results()[0]
                .unwrap()
                .cycles_to_failure
                .unwrap(),
            20000.0,
            1e-6,
        );
        assert_eq!(
            run.finish(cx).unwrap(),
            resumed.unwrap().finish(cx).unwrap()
        );
    });
}

#[test]
fn nonclosure_and_late_second_shear_refusal_restore_open_prefix_before_smaller_retry() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = shear_field(&mesh, GAMMA);
        let factors = [1.0, -1.0, 0.0, -1.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let wide = fs_gauge(&f.law, &f.calibration);
        let mut narrow = wide.clone();
        narrow.calibration = Arc::new(
            FatemiSocieCalibration::new(
                f.calibration.scope().clone(),
                f.calibration.parameters(),
                30000.0,
                1e8,
            )
            .unwrap(),
        );
        let mut run = ShellFatiguePath::new_with_fatemi_socie(
            &p,
            vec![normal_gauge(&f.law)],
            vec![wide, narrow],
            config(),
            shear_config(),
            cx,
        )
        .unwrap();
        run.enable_fatemi_socie_miner(2, cx).unwrap();
        run.advance_next(cx).unwrap();
        let prefix = run.clone();
        assert!(matches!(
            run.advance_next_cycle(cx),
            Err(ShellFatigueError::FatemiSocie {
                error: FatemiSocieError::NotClosed { .. },
                ..
            })
        ));
        same(&run, &prefix);
        run.advance_next(cx).unwrap();
        let prefix = run.clone();
        assert!(matches!(
            run.advance_next_cycle(cx),
            Err(ShellFatigueError::FatemiSocie { gauge: 1, .. })
        ));
        same(&run, &prefix);
        assert_eq!(run.fatemi_socie_results(), vec![None, None]);
        run.restart_fatemi_socie_cycle(cx).unwrap();
        assert_eq!(run.solution(), prefix.solution());
        assert_eq!(run.plane_summaries(0), prefix.plane_summaries(0));
        run.advance_next(cx).unwrap();
        run.advance_next_cycle(cx).unwrap();
        for result in run.fatemi_socie_results().into_iter().map(Option::unwrap) {
            assert_eq!(result.samples, 3);
            close(result.cycles_to_failure.unwrap(), 80000.0, 1e-5);
        }
        assert_eq!(run.sample_count(), 5);
        run.finish(cx).unwrap();
    });
}

#[test]
fn shear_source_binding_aggregate_work_and_window_storage_are_bounded() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = shear_field(&mesh, GAMMA);
        let factors = [1.0, -1.0, 0.0, -1.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let g = fs_gauge(&f.law, &f.calibration);
        let mut wrong_yield = g.clone();
        let mut parameters = f.calibration.parameters();
        parameters.initial_yield_stress_pa *= 2.0;
        wrong_yield.calibration = Arc::new(
            FatemiSocieCalibration::new(f.calibration.scope().clone(), parameters, 1.0, 1e8)
                .unwrap(),
        );
        let mut wrong_modulus = g.clone();
        parameters = f.calibration.parameters();
        parameters.shear_modulus_pa *= 2.0;
        wrong_modulus.calibration = Arc::new(
            FatemiSocieCalibration::new(f.calibration.scope().clone(), parameters, 1.0, 1e8)
                .unwrap(),
        );
        let mut wrong_card = g.clone();
        wrong_card.calibration = fixture(1).calibration;
        for bad in [wrong_yield, wrong_modulus, wrong_card] {
            assert!(matches!(
                ShellFatiguePath::new_with_fatemi_socie(
                    &p,
                    vec![],
                    vec![bad],
                    config(),
                    shear_config(),
                    cx
                ),
                Err(ShellFatigueError::SourceMismatch { .. })
            ));
        }
        for controls in [
            ShellFatemiSocieConfig {
                maximum_total_samples: 31,
                ..shear_config()
            },
            ShellFatemiSocieConfig {
                maximum_total_projection_work: 127,
                ..shear_config()
            },
            ShellFatemiSocieConfig {
                maximum_total_circle_work: 8191,
                ..shear_config()
            },
        ] {
            assert!(matches!(
                ShellFatiguePath::new_with_fatemi_socie(
                    &p,
                    vec![],
                    vec![g.clone(), g.clone()],
                    config(),
                    controls,
                    cx
                ),
                Err(ShellFatigueError::Budget)
            ));
        }
        let mut limited = g;
        limited.limits.maximum_samples = 3;
        let mut run = ShellFatiguePath::new_with_fatemi_socie(
            &p,
            vec![],
            vec![limited],
            config(),
            shear_config(),
            cx,
        )
        .unwrap();
        run.advance_next(cx).unwrap();
        run.advance_next(cx).unwrap();
        let prefix = run.clone();
        assert!(matches!(
            run.advance_next(cx),
            Err(ShellFatigueError::FatemiSocie {
                error: FatemiSocieError::Budget,
                ..
            })
        ));
        same(&run, &prefix);
        // Storage is bounded per selected representative window. A restart
        // keeps mechanics but releases the old window before sampling further.
        run.restart_fatemi_socie_cycle(cx).unwrap();
        run.advance_next(cx).unwrap();
        run.advance_next_cycle(cx).unwrap();
        assert_eq!(run.fatemi_socie_results()[0].unwrap().samples, 3);
        assert!(run.finish(cx).unwrap().is_empty());
    });
}
