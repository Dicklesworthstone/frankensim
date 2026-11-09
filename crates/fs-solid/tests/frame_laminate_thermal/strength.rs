//! Native accepted stress, hot strength and immutable selected-claim checks.

use super::*;
use fs_material::laminate::strength::{
    LaminaFace, LaminaFailureMode, LaminaMaximumStress, LaminaStrengthProperties,
    LaminaStrengthSelection,
};
use fs_solid::SectionModel;

const NAMES: [&str; 5] = [
    "tension_1",
    "compression_1",
    "tension_2",
    "compression_2",
    "shear_12",
];
const LIMITS: [f64; 5] = [1e6, 2e6, 2e6, 3e6, 0.7e6];

fn properties() -> LaminaStrengthProperties {
    LaminaStrengthProperties {
        tensile: [NAMES[0].into(), NAMES[2].into()],
        compressive: [NAMES[1].into(), NAMES[3].into()],
        shear: NAMES[4].into(),
    }
}
fn selections(layers: usize) -> Vec<LaminaStrengthSelection> {
    vec![
        LaminaStrengthSelection {
            properties: properties(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        };
        layers
    ]
}
fn weakening(t: f64) -> f64 {
    1.0 - 0.005 * (t - 300.0)
}
fn claim(name: &str, value: f64, lo: f64, hi: f64, provenance: &str) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, Pressure::DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(lo, value * weakening(lo)), (hi, value * weakening(hi))],
            dims: Pressure::DIMS,
        },
        interpolation: InterpolationPolicy::LinearInside,
        validity: ValidityDomain::unconstrained().with("T", lo, hi),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: provenance.into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn strength_source(plies: &[Ply], competing_claim: bool) -> Arc<LaminateThermalSource> {
    let original = source(plies, ReferenceSurface::Midsurface);
    let inputs = original
        .layers()
        .iter()
        .enumerate()
        .map(|(layer, original)| {
            let mut input = original.input().clone();
            let mut claims = input.card.claims().clone();
            for (mode, name) in NAMES.into_iter().enumerate() {
                let high = if competing_claim && mode == 0 {
                    340.0
                } else {
                    400.0
                };
                claims
                    .insert_claim(claim(
                        name,
                        LIMITS[mode] * (layer + 1) as f64,
                        260.0,
                        high,
                        "original strength source",
                    ))
                    .unwrap();
            }
            if competing_claim {
                claims
                    .insert_claim(claim(
                        NAMES[0],
                        10e6,
                        340.0,
                        400.0,
                        "different high-temperature strength source",
                    ))
                    .unwrap();
            }
            input.card = Arc::new(
                MaterialCard::assemble(
                    MaterialStateId {
                        chemistry: format!("thermal strength ply {layer}"),
                        phase: "solid".into(),
                        process: "numerical strength oracle".into(),
                        revision: 0,
                    },
                    claims,
                    vec![],
                )
                .unwrap(),
            );
            input
        })
        .collect();
    Arc::new(
        LaminateThermalSource::try_new(inputs, ReferenceSurface::Midsurface, || false).unwrap(),
    )
}

#[test]
fn free_bimetal_strength_finds_failure_hidden_by_zero_resultants() {
    with_cx(|cx| {
        let plies = bimetal();
        let temperatures = [340.0, 360.0];
        let width = 0.0005;
        let mut original = member(
            strength_source(&plies, false),
            width,
            0.1,
            0.0,
            &vec![temperatures.to_vec(); 3],
            cx,
        );
        original.static_load(&[0.0; 6], 1).unwrap();
        let force = original.element_forces(0);
        near(force[0], 0.0, 1.0);
        near(force[1], 0.0, 0.01);
        near(force[2], 0.0, 0.01);
        let bound = original
            .with_laminate_strength(0, &selections(2), cx)
            .unwrap();
        assert_eq!(bound.displacements(), original.displacements());
        assert_eq!(bound.element_forces(0), original.element_forces(0));
        assert_eq!(
            bound.inertial_forces(&[1.0; 6]).unwrap(),
            original.inertial_forces(&[1.0; 6]).unwrap()
        );
        let critical = bound.assess_laminate_element_strength(0, cx).unwrap();

        // Independent axial-only free-sided beam integrals and signed strengths.
        let (k, thermal, faces) = oracle(&plies, &temperatures, width);
        let [eps, curvature] = solve(k, thermal);
        let mut worst = (
            0.0,
            0,
            LaminaFace::Bottom,
            LaminaFailureMode::LongitudinalTension,
        );
        for (layer, ((&p, &t), positions)) in
            plies.iter().zip(&temperatures).zip(&faces).enumerate()
        {
            for (face, &y) in [LaminaFace::Bottom, LaminaFace::Top]
                .into_iter()
                .zip(positions)
            {
                let stress = young(p, t) * (eps - y * curvature - free(p, t));
                let mode = if stress >= 0.0 { 0 } else { 1 };
                let utilization = stress.abs() / (LIMITS[mode] * (layer + 1) as f64 * weakening(t));
                if utilization > worst.0 {
                    worst = (
                        utilization,
                        layer,
                        face,
                        if mode == 0 {
                            LaminaFailureMode::LongitudinalTension
                        } else {
                            LaminaFailureMode::LongitudinalCompression
                        },
                    );
                }
            }
        }
        near(critical.utilization(), worst.0, worst.0);
        assert_eq!(critical.element, 0);
        assert!(critical.station < 3);
        assert_eq!(critical.layer, worst.1);
        assert_eq!(critical.face, worst.2);
        assert_eq!(critical.critical_mode(), Some(worst.3));
        assert!(critical.reaches_limit());
        let station = bound.laminate_thermal_section(0, critical.station).unwrap();
        assert!(std::ptr::eq(
            station.strength().unwrap().thermal().prepared().as_ref(),
            station.prepared().source()
        ));
        assert!(
            original
                .laminate_thermal_section(0, 0)
                .unwrap()
                .strength()
                .is_none()
        );
    });
}

#[test]
fn hot_strength_changes_at_identical_mechanics_and_uses_exact_shared_snapshot() {
    with_cx(|cx| {
        let ply = Ply {
            h: 0.005,
            rho: 4300.0,
            young: 2e9,
            poisson: 0.31,
            alpha: 0.0,
            softening: 0.0,
        };
        let width = 0.0005;
        let mut frame = member(
            strength_source(&[ply], false),
            width,
            0.1,
            0.0,
            &vec![vec![300.0]; 3],
            cx,
        );
        let mut load = [0.0; 6];
        load[3] = 0.8e6 * width * ply.h;
        frame.static_load(&load, 1).unwrap();
        let cold = frame.with_laminate_strength(0, &selections(1), cx).unwrap();
        let before = cold.assess_laminate_element_strength(0, cx).unwrap();
        let hot = cold
            .with_laminate_temperatures(0, &vec![vec![350.0]; 3], cx)
            .unwrap();
        let after = hot.assess_laminate_element_strength(0, cx).unwrap();
        near(before.utilization(), 0.8, 1.0);
        near(after.utilization(), 0.8 / weakening(350.0), 1.0);
        assert!(!before.reaches_limit());
        assert!(after.reaches_limit());
        assert_eq!(hot.displacements(), cold.displacements());
        assert_eq!(hot.velocities(), cold.velocities());
        assert_eq!(
            hot.inertial_forces(&[1.0; 6]).unwrap(),
            cold.inertial_forces(&[1.0; 6]).unwrap()
        );
        let old = cold.laminate_thermal_section(0, 0).unwrap();
        let current = hot.laminate_thermal_section(0, 0).unwrap();
        assert_eq!(current.prepared().tangent(), old.prepared().tangent());
        let (eps, kappa) = hot.section_deformations(0)[0];
        assert!(
            old.strength()
                .unwrap()
                .assess_beam(current.prepared(), eps, kappa, || false)
                .is_err()
        );
        assert!(std::ptr::eq(
            current.strength().unwrap().thermal().prepared().as_ref(),
            current.prepared().source()
        ));
        near(
            current
                .strength()
                .unwrap()
                .layer_strength(0)
                .unwrap()
                .limits_pa()[0],
            LIMITS[0] * weakening(350.0),
            LIMITS[0],
        );
    });
}

#[test]
fn original_strength_claims_cannot_switch_or_rebind_and_failed_refresh_is_atomic() {
    with_cx(|cx| {
        let ply = Ply {
            h: 0.005,
            rho: 4300.0,
            young: 2e9,
            poisson: 0.31,
            alpha: 0.0,
            softening: 0.0,
        };
        let source = strength_source(&[ply], true);
        let original = member(
            Arc::clone(&source),
            0.0005,
            0.1,
            0.0,
            &vec![vec![300.0]; 3],
            cx,
        );
        let bound = original
            .with_laminate_strength(0, &selections(1), cx)
            .unwrap();
        // A different claim really is usable hot; the originally selected one
        // ends at 340 K. Thermal mechanics alone remains admitted at 350 K.
        let hot_source = source.at_temperatures(&[350.0], || false).unwrap();
        let hot_strength = LaminaMaximumStress::resolve(
            &source.layers()[0].input().card,
            &hot_source.layers()[0].elasticity().query(),
            properties(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        near(hot_strength.limits_pa()[0], 10e6 * weakening(350.0), 10e6);
        assert!(
            original
                .with_laminate_temperatures(0, &vec![vec![350.0]; 3], cx)
                .is_ok()
        );
        assert!(
            bound
                .with_laminate_temperatures(0, &[vec![310.0], vec![320.0], vec![350.0]], cx)
                .is_err()
        );
        assert!(bound.with_laminate_strength(0, &selections(1), cx).is_err());
        assert_eq!(bound.displacements(), original.displacements());
        assert_eq!(bound.element_forces(0), original.element_forces(0));
        for station in 0..3 {
            let current = bound.laminate_thermal_section(0, station).unwrap();
            assert_eq!(
                current.source().layers()[0]
                    .elasticity()
                    .state_coordinate("T"),
                Some(300.0)
            );
            assert_eq!(
                current
                    .strength()
                    .unwrap()
                    .layer_strength(0)
                    .unwrap()
                    .resolved()
                    .state_coordinate("T"),
                Some(300.0)
            );
        }
        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 21,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                original.with_laminate_strength(0, &selections(1), &cancelled),
                Err(SolidError::Cancelled)
            ));
            assert!(matches!(
                bound.assess_laminate_frame_strength(&cancelled),
                Err(SolidError::Cancelled)
            ));
            assert!(matches!(
                bound.with_laminate_temperatures(0, &vec![vec![320.0]; 3], &cancelled),
                Err(SolidError::Cancelled)
            ));
        });
    });
}

#[test]
fn whole_frame_reports_real_member_and_refuses_any_unchecked_member() {
    with_cx(|cx| {
        let plies = bimetal();
        let src = strength_source(&plies, false);
        let mut frame = member(
            Arc::clone(&src),
            0.0005,
            0.1,
            0.0,
            &vec![vec![300.0; 2]; 3],
            cx,
        );
        let i = frame.add_node(0.0, 0.2);
        let j = frame.add_node(0.1, 0.2);
        frame.fix(i, [true; 3]);
        frame
            .add_laminate_thermal_element(
                i,
                j,
                src,
                0.0005,
                &[vec![310.0, 330.0], vec![330.0, 350.0], vec![350.0, 370.0]],
                0.02,
                cx,
            )
            .unwrap();
        frame.static_load(&[0.0; 12], 1).unwrap();
        let partially_bound = frame.with_laminate_strength(1, &selections(2), cx).unwrap();
        let one = partially_bound
            .assess_laminate_element_strength(1, cx)
            .unwrap();
        assert!(partially_bound.assess_laminate_frame_strength(cx).is_err());
        let full = partially_bound
            .with_laminate_strength(0, &selections(2), cx)
            .unwrap();
        let all = full.assess_laminate_frame_strength(cx).unwrap();
        assert_eq!(all, one);
        assert_eq!(all.element, 1);
        assert!(all.utilization() > 0.0);
        let empty = Frame2d::new(Geometry::Linear);
        assert!(empty.assess_laminate_frame_strength(cx).is_err());
        let mut ordinary = Frame2d::new(Geometry::Linear);
        ordinary.add_node(0.0, 0.0);
        ordinary.add_node(1.0, 0.0);
        ordinary
            .add_element(0, 1, 3, 1.0, &|| SectionModel::Elastic {
                ea: 100.0,
                ei: 1.0,
            })
            .unwrap();
        assert!(ordinary.assess_laminate_frame_strength(cx).is_err());
    });
}
