//! Original frame equilibria and physical clocks drive local ductile initiation.
//! Synthetic coefficients provide implementation oracles, not alloy validation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::ductile::{
    DuctileDomain, DuctileState, JohnsonCookInitiation, JohnsonCookParameters,
    TabulatedDuctileLocus,
};
use fs_material::fiber::MenegottoPintoSteel;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::reduction::{UniaxialJ2, UniaxialTabulated, UniaxialVoce};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening, TabulatedJ2};
use fs_material::plastic::voce::VoceJ2;
use fs_solid::ductile::DuctileMaterialBinding;
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::ductile::{
    FrameDuctileConfig, FrameDuctileError, FrameDuctileGauge, FrameDuctilePath,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const E: f64 = 200e9;
const SY: f64 = 200e6;
const HI: f64 = 2e9;
const HK: f64 = 18e9;
const AREA: f64 = 1e-3;
const ARM: f64 = 0.1;
const PATTERN: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 91,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    });
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn backbone() -> J2Plasticity {
    J2Plasticity::new_mixed(IsotropicElastic::new(E, 0.3, 0.04).unwrap(), SY, HI, HK).unwrap()
}
fn metal() -> FiberLaw {
    FiberLaw::J2(UniaxialJ2::new(backbone()).unwrap())
}
fn frame(law: &FiberLaw, bending: bool) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    let left = frame.add_node(0.0, 0.0);
    let right = frame.add_node(1.0, 0.0);
    frame.fix(left, [true; 3]);
    frame.fix(
        right,
        if bending {
            [true, false, false]
        } else {
            [false, true, true]
        },
    );
    frame.add_mass(right, 100.0);
    frame
        .add_element(left, right, 3, 0.0, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: law.clone(),
                    })
                    .to_vec(),
            ))
        })
        .unwrap();
    frame
}
fn config() -> FrameDuctileConfig {
    FrameDuctileConfig {
        maximum_gauges: 6,
        maximum_accepted_steps: 1000,
        maximum_absolute_axial_strain: 0.04,
    }
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-0.5, 0.5],
        strain_rate_s_inv: [1e-12, 1e8],
        temperature_k: [300.0, 300.0],
    }
}
fn constant(capacity: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-0.5, capacity), (0.5, capacity)], domain()).unwrap()
}
fn gauges<L>(law: &L) -> Vec<FrameDuctileGauge<'_, L>> {
    (0..3)
        .flat_map(|section| {
            (0..2).map(move |fiber| FrameDuctileGauge {
                element: 0,
                section,
                fiber,
                law,
                temperature_k: 300.0,
                material_source: None,
            })
        })
        .collect()
}
fn path<'a, L: DuctileMaterialBinding>(
    frame: Frame2d,
    law: &'a L,
    cx: &Cx<'_>,
) -> FrameDuctilePath<'a, L> {
    FrameDuctilePath::new(frame, gauges(law), 10.0, config(), cx).unwrap()
}
// Independent associative Prager return in the axial stress subspace.
fn prager(strain: f64, plastic: &mut f64, alpha: &mut f64) -> f64 {
    let trial = E * (strain - *plastic) - HK * *plastic;
    let excess = trial.abs() - SY - HI * *alpha;
    if excess > 0.0 {
        let increment = excess / (E + HI + HK);
        *plastic += trial.signum() * increment;
        *alpha += increment;
    }
    E * (strain - *plastic)
}

#[test]
fn reversed_plastic_flow_uses_local_physical_stress_and_retained_path_length() {
    with_cx(|_, cx| {
        // epsilon_f = 0.03 - 0.03 eta: 0.02 in tension, 0.04 in compression.
        let locus =
            TabulatedDuctileLocus::new(vec![(-0.5, 0.045), (0.5, 0.015)], domain()).unwrap();
        let mut run = path(frame(&metal(), false), &locus, cx);
        let mut batched = run.clone();
        let targets = [0.003, 0.0025, -0.003, 0.002];
        let (mut plastic, mut alpha, mut omega) = (0.0, 0.0, 0.0);
        for target in targets {
            let old_alpha = alpha;
            let stress = prager(target, &mut plastic, &mut alpha);
            omega += (alpha - old_alpha) / if stress > 0.0 { 0.02 } else { 0.04 };
            run.displacement_control(&PATTERN, 3, &[target], 0.125, cx)
                .unwrap();
            for (sample, state) in run.samples().iter().zip(run.points()) {
                close(sample.axial_strain, target, 1e-11);
                close(sample.stress_pa, stress, 2e-9);
                close(sample.equivalent_plastic_strain, alpha, 1e-11);
                close(state.accumulated_plastic_strain(), alpha, 1e-11);
                close(state.utilization(), omega, 1e-9);
            }
        }
        assert!(
            alpha > plastic.abs(),
            "reverse plastic flow must not heal consumed ductility"
        );
        batched
            .displacement_control(&PATTERN, 3, &targets, 0.125, cx)
            .unwrap();
        close(batched.points()[0].utilization(), omega, 1e-9);
        assert_eq!(run.time_s(), batched.time_s());
        assert_eq!(run.accepted_steps(), 4);
        let mut resumed = run.clone();
        run.displacement_control(&PATTERN, 3, &[0.0015], 0.125, cx)
            .unwrap();
        resumed
            .displacement_control(&PATTERN, 3, &[0.0015], 0.125, cx)
            .unwrap();
        assert_eq!(run.points(), resumed.points());
        assert_eq!(run.samples(), resumed.samples());
        assert_eq!(run.frame().displacements(), resumed.frame().displacements());
    });
}

#[test]
fn opposing_plastic_bending_fibers_do_not_disappear_in_a_zero_axial_resultant() {
    with_cx(|_, cx| {
        let locus =
            TabulatedDuctileLocus::new(vec![(-0.5, 0.045), (0.5, 0.015)], domain()).unwrap();
        let mut run = path(frame(&metal(), true), &locus, cx);
        let pattern = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        run.displacement_control(&pattern, 5, &[0.03], 0.25, cx)
            .unwrap();
        close(run.frame().element_forces(0)[0], 0.0, 1e-5);
        for section in 0..3 {
            let positive = run.samples()[2 * section];
            let negative = run.samples()[2 * section + 1];
            assert!(positive.stress_pa > 0.0 && negative.stress_pa < 0.0);
            close(positive.stress_pa, -negative.stress_pa, 1e-9);
            assert!(positive.equivalent_plastic_strain > 0.0);
            close(
                run.points()[2 * section].utilization(),
                2.0 * run.points()[2 * section + 1].utilization(),
                1e-9,
            );
        }
    });
}

#[test]
fn actual_duration_drives_rate_and_late_initiation_rolls_back_the_complete_block() {
    with_cx(|_, cx| {
        let rate = JohnsonCookInitiation::new(
            JohnsonCookParameters {
                d1: 0.04,
                d2: 0.0,
                d3: 0.0,
                d4: 0.02,
                d5: 0.0,
                reference_strain_rate_s_inv: 0.01,
                transition_temperature_k: 300.0,
                melting_temperature_k: 1000.0,
            },
            domain(),
        )
        .unwrap();
        let mut fast = path(frame(&metal(), false), &rate, cx);
        let mut slow = fast.clone();
        fast.displacement_control(&PATTERN, 3, &[0.003], 0.1, cx)
            .unwrap();
        slow.displacement_control(&PATTERN, 3, &[0.003], 10.0, cx)
            .unwrap();
        assert_eq!(fast.frame().displacements(), slow.frame().displacements());
        assert_eq!(fast.samples(), slow.samples());
        let alpha = (E * 0.003 - SY) / (E + HI + HK);
        for (run, dt) in [(&fast, 0.1), (&slow, 10.0)] {
            let capacity = 0.04 * (1.0 + 0.02 * (alpha / dt / 0.01).ln());
            close(run.points()[0].utilization(), alpha / capacity, 1e-9);
            assert_eq!(run.time_s(), 10.0 + dt);
        }
        assert!(slow.points()[0].utilization() > fast.points()[0].utilization());

        let narrow = constant(0.003);
        let broad = constant(0.03);
        let mut selected = gauges(&broad);
        // The later gauge is the only point reaching initiation.
        selected[5].law = &narrow;
        let mut run =
            FrameDuctilePath::new(frame(&metal(), false), selected, 10.0, config(), cx).unwrap();
        let mut unchanged = run.clone();
        match run
            .displacement_control(&PATTERN, 3, &[0.0015, 0.006], 0.125, cx)
            .unwrap_err()
        {
            FrameDuctileError::Initiated {
                gauge,
                step,
                interval_s,
                candidate,
            } => {
                assert_eq!((gauge, step), (5, 1));
                assert_eq!(interval_s, [10.125, 10.25]);
                assert!(candidate.state.initiated());
                assert!(candidate.initiation_fraction.unwrap() > 0.0);
            }
            other => panic!("expected local initiation, got {other:?}"),
        }
        assert_eq!(run.points(), &[DuctileState::virgin(); 6]);
        assert_eq!(run.samples(), unchanged.samples());
        assert_eq!(
            run.frame().displacements(),
            unchanged.frame().displacements()
        );
        assert_eq!(run.time_s(), 10.0);
        assert_eq!(run.accepted_steps(), 0);
        for checkpoint in [&mut run, &mut unchanged] {
            checkpoint
                .displacement_control(&PATTERN, 3, &[0.0015, 0.003], 0.125, cx)
                .unwrap();
        }
        assert_eq!(run.samples(), unchanged.samples());
        assert_eq!(run.points(), unchanged.points());
    });
}

#[test]
fn native_load_ramp_and_newmark_keep_original_motion_work_and_all_subincrements() {
    with_cx(|_, cx| {
        let locus = constant(0.5);
        let mut run = path(frame(&metal(), false), &locus, cx);
        let mut native = run.frame().clone();
        let load = [0.0, 0.0, 0.0, 600_000.0, 0.0, 0.0];
        let expected = native.static_load(&load, 4).unwrap();
        let actual = run.static_load(&load, 4, 0.125, cx).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(run.accepted_steps(), 4);
        assert_eq!(run.time_s(), 10.5);
        let (mut plastic, mut alpha) = (0.0, 0.0);
        for step in &expected {
            prager(step.u[3], &mut plastic, &mut alpha);
        }
        close(run.points()[0].utilization(), alpha / 0.5, 1e-9);

        let mut run = path(frame(&metal(), false), &locus, cx);
        let mut native = run.frame().clone();
        let mut split = run.clone();
        let mut ground = vec![0.0; 101];
        ground[1..51].fill(-20_000.0);
        let dt = 1e-5;
        let damping = Rayleigh { a0: 1.0, a1: 1e-6 };
        let expected = native
            .newmark_with_load(&ground, dt, (1.0, 0.0), damping, &[0.0; 6])
            .unwrap();
        let actual = run
            .newmark_with_load(&ground, dt, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(actual, expected);
        let (mut plastic, mut alpha) = (0.0, 0.0);
        for u in &expected.u {
            prager(u[3], &mut plastic, &mut alpha);
        }
        assert!(
            alpha > 0.0,
            "native dynamic trajectory must exercise yielding"
        );
        close(run.points()[0].utilization(), alpha / 0.5, 1e-9);
        assert_eq!(run.accepted_steps(), 100);
        split
            .newmark_with_load(&ground[..51], dt, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        split
            .newmark_with_load(&ground[50..], dt, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(split.time_s(), run.time_s());
        assert_eq!(split.accepted_steps(), run.accepted_steps());
        close(
            split.points()[0].utilization(),
            run.points()[0].utilization(),
            1e-9,
        );
        close(
            split.frame().displacements()[3],
            run.frame().displacements()[3],
            1e-10,
        );

        let narrow = constant(1e-6);
        let mut refused = path(frame(&metal(), false), &narrow, cx);
        let unchanged = refused.clone();
        assert!(matches!(
            refused.newmark_with_load(&ground, dt, (1.0, 0.0), damping, &[0.0; 6], cx),
            Err(FrameDuctileError::Initiated { .. })
        ));
        assert_eq!(
            refused.frame().displacements(),
            unchanged.frame().displacements()
        );
        assert_eq!(refused.frame().velocities(), unchanged.frame().velocities());
        assert_eq!(
            refused.frame().accelerations(),
            unchanged.frame().accelerations()
        );
        assert_eq!(refused.points(), unchanged.points());
        assert_eq!(refused.time_s(), unchanged.time_s());
    });
}

fn branches<const N: usize>() -> FiberLaw {
    FiberLaw::chaboche_branches(
        UniaxialChaboche::new(
            ChabocheJ2::<N>::new(
                VoceJ2::new(backbone(), 80e6, 500.0).unwrap(),
                std::array::from_fn(|i| {
                    BackstressBranch::new(6e9 + 2e9 * i as f64, 50.0 * (i + 1).pow(2) as f64)
                        .unwrap()
                }),
            )
            .unwrap(),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn native_nonlinear_hardening_families_keep_their_original_committed_stress() {
    with_cx(|_, cx| {
        let table = TabulatedHardening::new(vec![
            HardeningPoint {
                plastic_strain: 0.0,
                yield_stress_pa: SY,
            },
            HardeningPoint {
                plastic_strain: 0.05,
                yield_stress_pa: SY + HI * 0.05,
            },
        ])
        .unwrap();
        let tabulated = TabulatedJ2::new(
            J2Plasticity::new_mixed(IsotropicElastic::new(E, 0.3, 0.04).unwrap(), SY, 0.0, HK)
                .unwrap(),
            table,
        )
        .unwrap();
        let laws = [
            FiberLaw::Voce(
                UniaxialVoce::new(VoceJ2::new(backbone(), 80e6, 500.0).unwrap()).unwrap(),
            ),
            FiberLaw::Tabulated(UniaxialTabulated::new(tabulated).unwrap()),
            branches::<1>(),
            branches::<2>(),
            branches::<3>(),
            branches::<4>(),
            branches::<5>(),
            branches::<6>(),
            branches::<7>(),
            branches::<8>(),
        ];
        let locus = constant(0.5);
        for law in laws {
            let mut run = path(frame(&law, false), &locus, cx);
            let mut native = run.frame().clone();
            let targets = [0.003, -0.0025, 0.002];
            let expected = native.displacement_control(&PATTERN, 3, &targets).unwrap();
            assert_eq!(
                run.displacement_control(&PATTERN, 3, &targets, 0.125, cx)
                    .unwrap(),
                expected
            );
            let physical_stress = native.element_forces(0)[0] / (2.0 * AREA);
            for (sample, point) in run.samples().iter().zip(run.points()) {
                close(sample.stress_pa, physical_stress, 1e-8);
                assert!(sample.equivalent_plastic_strain > 0.0);
                close(
                    point.utilization(),
                    sample.equivalent_plastic_strain / 0.5,
                    1e-10,
                );
            }
        }
    });
}

#[test]
fn unsupported_histories_preload_budgets_and_cancellation_do_not_publish() {
    with_cx(|gate, cx| {
        let locus = constant(0.5);
        let mut prior = frame(&metal(), false);
        prior.displacement_control(&PATTERN, 3, &[0.003]).unwrap();
        assert!(matches!(
            FrameDuctilePath::new(prior, gauges(&locus), 0.0, config(), cx),
            Err(FrameDuctileError::Invalid(_))
        ));
        let empirical = FiberLaw::Steel(MenegottoPintoSteel::new(E, SY, 0.1).unwrap());
        assert!(matches!(
            FrameDuctilePath::new(frame(&empirical, false), gauges(&locus), 0.0, config(), cx),
            Err(FrameDuctileError::Invalid(_))
        ));
        let mut duplicates = gauges(&locus);
        duplicates[1] = duplicates[0];
        assert!(matches!(
            FrameDuctilePath::new(frame(&metal(), false), duplicates, 0.0, config(), cx),
            Err(FrameDuctileError::Invalid(_))
        ));
        let mut controls = config();
        controls.maximum_accepted_steps = 1;
        let mut run =
            FrameDuctilePath::new(frame(&metal(), false), gauges(&locus), 0.0, controls, cx)
                .unwrap();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.0005, 0.001], 0.1, cx),
            Err(FrameDuctileError::Budget)
        ));
        assert_eq!(run.time_s(), 0.0);
        run.displacement_control(&PATTERN, 3, &[0.0005], 0.1, cx)
            .unwrap();
        let unchanged = run.clone();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.001], f64::NAN, cx),
            Err(FrameDuctileError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.001], 0.1, cx),
            Err(FrameDuctileError::Cancelled)
        ));
        assert_eq!(
            run.frame().displacements(),
            unchanged.frame().displacements()
        );
        assert_eq!(run.points(), unchanged.points());
        assert_eq!(run.time_s(), unchanged.time_s());
    });
}

mod source {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_material::ductile::source::{
        DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
        STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
    };
    use fs_material::plastic::source::HardeningPropertyBasis;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Density, Dims, Pressure, QuantitySpec};

    fn temperature() -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ))
    }
    fn point() -> QueryPoint {
        QueryPoint::new()
            .with_quantity("T", temperature(), 300.0)
            .unwrap()
    }
    fn card(capacity: f64) -> MaterialCard {
        let fixed = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 300.0);
        let provenance = Provenance {
            source: "synthetic native frame ductility".into(),
            license: "test-only".into(),
            artifact: None,
        };
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, 7800.0),
            ("young_modulus", Pressure::DIMS, E),
            ("poisson_ratio", Dims::NONE, 0.3),
            ("yield_stress", Pressure::DIMS, SY),
            ("isotropic_hardening_modulus", Pressure::DIMS, HI),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: fixed.clone(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: STRESS_TRIAXIALITY_AXIS.into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(-1.0, capacity), (1.0, capacity)],
                    dims: Dims::NONE,
                },
                validity: fixed
                    .with_quantity(
                        STRESS_TRIAXIALITY_AXIS,
                        QuantitySpec::dimensional(Dims::NONE),
                        -1.0,
                        1.0,
                    )
                    .with_quantity(
                        PLASTIC_STRAIN_RATE_AXIS,
                        QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                        1e-8,
                        1.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::TabulatedOnly,
                observations: vec![],
                provenance,
            })
            .unwrap();
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic frame alloy".into(),
                phase: "solid".into(),
                process: "synthetic fixed condition".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn original_card_state_temperature_and_native_law_cannot_be_substituted() {
        with_cx(|_, cx| {
            let high_card = card(0.02);
            let low_card = card(0.001);
            let resolve = |card| {
                FiberMaterial::resolve_j2(
                    card,
                    &point(),
                    "density",
                    HardeningPropertyBasis::PlasticStrainModulus,
                    0.04,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap()
            };
            let high = resolve(&high_card);
            let low = resolve(&low_card);
            let calibration_point = point()
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -1.0,
                )
                .unwrap()
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    0.01,
                )
                .unwrap();
            let request = DuctileSourceRequest {
                temperature_axis: "T".into(),
                maximum_points: 2,
            };
            let locus = resolve_ductile_locus(
                &high_card,
                &calibration_point,
                &request,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap();
            assert!(matches!(
                FrameDuctilePath::new(frame(high.law(), false), gauges(&locus), 0.0, config(), cx),
                Err(FrameDuctileError::SourceMismatch { .. })
            ));
            let mut bound = gauges(&locus);
            for gauge in &mut bound {
                gauge.material_source = Some(&high);
            }
            let mut accepted =
                FrameDuctilePath::new(frame(high.law(), false), bound.clone(), 0.0, config(), cx)
                    .unwrap();
            accepted
                .displacement_control(&PATTERN, 3, &[0.003], 1.0, cx)
                .unwrap();
            assert!(accepted.points()[0].utilization() > 0.0);
            bound[0].material_source = Some(&low);
            assert!(matches!(
                FrameDuctilePath::new(frame(high.law(), false), bound.clone(), 0.0, config(), cx),
                Err(FrameDuctileError::SourceMismatch { gauge: 0, .. })
            ));
            bound[0].material_source = Some(&high);
            bound[0].temperature_k = 301.0;
            assert!(matches!(
                FrameDuctilePath::new(frame(high.law(), false), bound.clone(), 0.0, config(), cx),
                Err(FrameDuctileError::SourceMismatch { gauge: 0, .. })
            ));
            bound[0].temperature_k = 300.0;
            assert!(matches!(
                FrameDuctilePath::new(frame(&metal(), false), bound, 0.0, config(), cx),
                Err(FrameDuctileError::SourceMismatch { .. })
            ));
        });
    }
}
