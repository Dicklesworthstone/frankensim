//! Source loss -> finite enthalpy -> refreshed interface -> continuing motion.
#![cfg(feature = "dissipative-coupling")]
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, So3, Twist, Vec3};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedLinearInterface, ScalarAdmissibility,
    ScalarPropertyRequirement, interface_stiffness_per_area_dims,
    interface_viscosity_per_area_dims, resolve_interface_state_point,
};
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::ResolvedBody;
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::motion::damping::{
    DampedRigidPair, DamperAttachment, DampingError, DampingSolve, MaterialDamper,
};
use fs_mbd_material::thermal::heating::motion::{HeatedMotionControls, HeatedRigidBody};
use fs_qty::Temperature;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::{CellId, DepositionError};

fn interface(solid: &MaterialCard, slope: f64, upper: f64, stiffness: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("k", interface_stiffness_per_area_dims(), [stiffness; 2]),
        (
            "c",
            interface_viscosity_per_area_dims(),
            [1., 1. + slope * (upper - 300.)],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., values[0]), (upper, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300., upper)
                    .with("pressure", 1e5, 2e5)
                    .with("omega", 0., 1e4),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic dashpot fixture, not measured data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = SurfaceSpec {
        material: solid.id().clone(),
        texture_frame: "attachment".into(),
    };
    InterfaceSystemCard::assemble(
        surface.clone(),
        surface,
        SystemContext {
            medium: "fixture".into(),
            third_body: None,
            environment: "fixed pressure".into(),
            history: "unchanged synthetic layer".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn pair<'a>(
    reference: &'a ResolvedBody,
    solid: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    interface: &'a InterfaceSystemCard,
    offset: f64,
    spin: f64,
    watts: [f64; 2],
) -> DampedRigidPair<'a> {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            "k",
            interface_stiffness_per_area_dims(),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "c",
            interface_viscosity_per_area_dims(),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    let source = resolve_interface_state_point(
        interface,
        &point(300., 1e5).with("omega", 0.).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let bodies = [0, 1].map(|i| {
        let sign = if i == 0 { 1. } else { -1. };
        HeatedRigidBody::new(
            heated(reference, solid, chart, watts[i], vec![]),
            Se3::from_parts(So3::identity(), Vec3::new(-sign, 0., 0.)).unwrap(),
            Twist::new(Vec3::new(0., 0., spin), Vec3::new(sign, 0., 0.)),
            RigidBodyIntegrator::new(Gravity::ZERO),
            HeatedMotionControls {
                maximum_step_s: 0.125,
                relative_tensor_tolerance: 1e-12,
                transfer: ConfigurationTransferControls {
                    relative_momentum_tolerance: 1e-12,
                    maximum_parameter_work_j: 1.,
                },
                maximum_numerical_energy_change_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap()
    });
    DampedRigidPair::new(
        bodies,
        MaterialDamper {
            card: interface,
            source: ResolvedLinearInterface::new(source, "k", "c", 1.).unwrap(),
            attachments: [DamperAttachment {
                cell: CellId(0),
                reference_point_m: Vec3::new(0., offset, 0.),
            }; 2],
            temperature_side: 0,
            first_heat_fraction: 0.5,
        },
    )
    .unwrap()
}

fn solve() -> DampingSolve {
    DampingSolve {
        maximum_cells: 2,
        maximum_trials: 512,
        maximum_newton_iterations: 32,
        difference_step_j: 1e-5,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 1e-8,
        maximum_relaxation_step: 2.,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-11,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn unchanged(old: &DampedRigidPair<'_>, saved: &DampedRigidPair<'_>) {
    for (a, b) in old.bodies().iter().zip(saved.bodies()) {
        assert_eq!(a.state(), b.state());
        assert_eq!(a.thermal().cells(), b.thermal().cells());
        assert_eq!(a.thermal().time_s(), b.thermal().time_s());
        assert_eq!(a.configuration(), b.configuration());
    }
    assert_eq!(old.damper().source, saved.damper().source);
}

#[test]
fn central_impulse_matches_independent_momentum_heat_and_free_flight() {
    let solid = card(10., 0., 500.);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 0., 500., 0.);
    let (bodies, mut damper) =
        pair(&reference, &solid, &chart, &interface, 0., 0., [2., 3.]).into_parts();
    damper.first_heat_fraction = 0.25;
    let old = DampedRigidPair::new(bodies, damper).unwrap();
    let saved = old.clone();
    let (next, audit) = old.step_to(0.1, solve(), &[0.3, 0.6], || false).unwrap();
    // Two 2 kg cubes, c=1 Ns/m: momentum balance gives J and midpoint work.
    let impulse = 0.1 * 2. / 1.05;
    let speed = 1. - impulse / 2.;
    let heat = 2. * (1. - speed * speed);
    for (i, body) in next.bodies().iter().enumerate() {
        let sign = if i == 0 { 1. } else { -1. };
        close(body.state().linear_momentum_world().x, sign * 2. * speed);
        close(
            body.state().pose().position_world().x,
            -sign + sign * speed * 0.1,
        );
        let q = [0.2, 0.3][i] + [0.25, 0.75][i] * heat;
        close(body.thermal().energy_j(), q);
        close(
            body.thermal().temperature_k(CellId(0)).unwrap(),
            300. + q / 20.,
        );
        close(audit.heat_by_body_j[i], [0.25, 0.75][i] * heat);
    }
    close(audit.dissipated_heat_j, heat);
    close(audit.impulse.kinetic_energy_change_j, -heat);
    close(audit.external_heat_j, 0.5);
    close(audit.combined_residual_j, 0.);
    assert_eq!(audit.parameter_work_j, [0.; 2]);
    let (_, second) = next.step_to(0.2, solve(), &[0.3, 0.6], || false).unwrap();
    close(second.external_heat_j, 0.5); // Generated heat did not become a persistent source.
    unchanged(&old, &saved);
}

#[test]
fn heat_capacity_only_swap_changes_temperature_damping_and_the_continuing_trajectory() {
    let mut endpoints = Vec::new();
    for cp in [0.5, 5.] {
        let solid = card(cp, 0., 500.);
        let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
        let reference = reference(&solid, 1);
        let interface = interface(&solid, 1., 500., 0.);
        let initial = pair(&reference, &solid, &chart, &interface, 0., 0., [0.; 2]);
        let mut actual = initial.clone();
        let (mut speed, mut rise): (f64, f64) = (1., 0.);
        for step in 1..=4 {
            // Independent scalar bisection: total heat D, C=m Cp per body,
            // c=1+(T-300), equal partition, midpoint relative velocity.
            let (mut lo, mut hi) = (0., 2. * speed * speed);
            for _ in 0..80 {
                let q = 0.5 * (lo + hi);
                let c = 1. + rise + q / (4. * cp);
                let predicted = c * 0.1 * (2. * speed / (1. + 0.05 * c)).powi(2);
                if q < predicted {
                    lo = q;
                } else {
                    hi = q;
                }
            }
            let q = 0.5 * (lo + hi);
            rise += q / (4. * cp);
            let c = 1. + rise;
            speed *= (1. - 0.05 * c) / (1. + 0.05 * c);
            let (next, audit) = actual
                .step_to(step as f64 * 0.1, solve(), &[0.1; 2], || false)
                .unwrap();
            close(audit.dissipated_heat_j, q);
            close(audit.source.damping_n_s_m(), c);
            close(
                next.bodies()[0].thermal().temperature_k(CellId(0)).unwrap(),
                300. + rise,
            );
            close(
                next.bodies()[0].state().linear_momentum_world().x,
                2. * speed,
            );
            close(audit.combined_residual_j, 0.);
            for (original, refreshed) in initial
                .damper()
                .source
                .resolved()
                .properties()
                .iter()
                .zip(audit.source.resolved().properties())
            {
                assert_eq!(
                    original.answer().receipt.selected,
                    refreshed.answer().receipt.selected
                );
            }
            actual = next;
        }
        let mut replay = initial;
        for step in 1..=4 {
            replay = replay
                .step_to(step as f64 * 0.1, solve(), &[0.1; 2], || false)
                .unwrap()
                .0;
        }
        unchanged(&actual, &replay);
        endpoints.push(speed);
    }
    assert!(endpoints[0] < endpoints[1]);
}

#[test]
fn off_centre_damping_conserves_world_momenta_and_separates_expansion_work() {
    let solid = card(1., 1e-2, 500.);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 0., 500., 0.);
    let old = pair(&reference, &solid, &chart, &interface, 0.25, 1., [0.; 2]);
    let (next, audit) = old.step_to(0.05, solve(), &[0.1; 2], || false).unwrap();
    // I=1/3 kg m²; g=2(1/m+r²/I)=1.375 kg^-1.
    let impulse = 0.05 * 2. / (1. + 0.5 * 0.05 * 1.375);
    let heat = 0.05 * (2.0_f64 / (1. + 0.5 * 0.05 * 1.375)).powi(2);
    let stretch = 1. + 1e-2 * heat / 4.;
    close(audit.dissipated_heat_j, heat);
    let before = old.bodies().each_ref().map(|b| b.diagnostics().unwrap());
    let after = next.bodies().each_ref().map(|b| b.diagnostics().unwrap());
    let p0 = before[0]
        .linear_momentum_world
        .add(before[1].linear_momentum_world);
    let p1 = after[0]
        .linear_momentum_world
        .add(after[1].linear_momentum_world);
    let l0 = before[0]
        .angular_momentum_world
        .add(before[1].angular_momentum_world);
    let l1 = after[0]
        .angular_momentum_world
        .add(after[1].angular_momentum_world);
    for (a, b) in [
        (p0.x, p1.x),
        (p0.y, p1.y),
        (p0.z, p1.z),
        (l0.x, l1.x),
        (l0.y, l1.y),
        (l0.z, l1.z),
    ] {
        close(a, b);
    }
    for i in 0..2 {
        let spin_momentum = 1. / 3. + [1., -1.][i] * 0.25 * impulse;
        let kinetic_after_impulse = 1.5 * spin_momentum * spin_momentum;
        close(
            audit.parameter_work_j[i],
            kinetic_after_impulse * (1. / (stretch * stretch) - 1.),
        );
        close(next.bodies()[i].thermal().energy_j(), 0.5 * heat);
        assert!(audit.parameter_work_j[i] < 0.);
    }
    close(audit.combined_residual_j, 0.);
}

#[test]
fn late_cancellation_and_material_or_solver_refusals_leave_both_histories_intact() {
    let solid = card(1., 0., 500.);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 1., 500., 0.);
    let old = pair(&reference, &solid, &chart, &interface, 0., 0., [0.; 2]);
    let saved = old.clone();
    let mut checkpoints = 0;
    old.step_to(0.1, solve(), &[0.1; 2], || {
        checkpoints += 1;
        false
    })
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        old.step_to(0.1, solve(), &[0.1; 2], || {
            calls += 1;
            calls == checkpoints
        }),
        Err(DampingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(calls, checkpoints);
    let mut limited = solve();
    limited.maximum_trials = 1;
    assert!(matches!(
        old.step_to(0.1, limited, &[0.1; 2], || false),
        Err(DampingError::Budget(_))
    ));
    limited = solve();
    limited.maximum_relaxation_step = 1e-4;
    assert!(matches!(
        old.step_to(0.1, limited, &[0.1; 2], || false),
        Err(DampingError::Budget(_))
    ));
    unchanged(&old, &saved);
    let narrow = self::interface(&solid, 1., 300.001, 0.);
    let narrow_pair = pair(&reference, &solid, &chart, &narrow, 0., 0., [0.; 2]);
    assert!(matches!(
        narrow_pair.step_to(0.1, solve(), &[0.1; 2], || false),
        Err(DampingError::Material(_))
    ));
    close(narrow_pair.bodies()[0].thermal().time_s(), 0.);
    let spring = self::interface(&solid, 0., 500., 1.);
    let (bodies, mut damper) = old.into_parts();
    damper.card = &spring;
    let requirements = [
        ScalarPropertyRequirement::try_new(
            "k",
            interface_stiffness_per_area_dims(),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "c",
            interface_viscosity_per_area_dims(),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    damper.source = ResolvedLinearInterface::new(
        resolve_interface_state_point(
            &spring,
            &point(300., 1e5).with("omega", 0.).unwrap(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        "k",
        "c",
        1.,
    )
    .unwrap();
    assert!(matches!(
        DampedRigidPair::new(bodies, damper),
        Err(DampingError::Invalid(_))
    ));
}
