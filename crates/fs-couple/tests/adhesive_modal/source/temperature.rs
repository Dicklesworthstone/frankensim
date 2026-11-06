//! G0/G3: original thermal normal/junction sources in actual native friction.
use super::*;
use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveTangential;
use fs_contact::normal_patch::adhesive::{
    SourceJkrSphereContact,
    temperature::{AdhesiveTemperatures, AdhesiveThermalControl},
};
use fs_contact::normal_patch::maugis::{MaugisEnvelope, SourceMaugisSphereContact};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    ModalFrictionError, adhesive::temperature::AdhesiveModalSourceRefresh,
};
use fs_phs::RelaxationBranch;

#[path = "../../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod thermal_fixture;

fn temperatures(t: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0 + 0.5 * (t - 300.0),
        surface_b_k: 300.0 - 0.25 * (t - 300.0),
        interface_k: t,
    }
}
fn bound_native(source: SourceAdhesiveTangential, memory: bool) -> FrictionModalSystem {
    let native = mechanics_with_binding(1.0, AdhesiveModalBinding::from_source(source).unwrap());
    if memory {
        native
            .with_material_relaxation(
                vec![RelaxationBranch {
                    projection: vec![1.0 / MASS[0].sqrt(), 0.0, -1.0 / MASS[1].sqrt(), 0.0],
                    stiffness: 10.0,
                    relaxation_time_s: 0.002,
                }],
                Some(vec![0.001]),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    } else {
        native
    }
}
fn source_binding(s: &FrictionModalSystem) -> &SourceAdhesiveTangential {
    s.adhesive_binding().unwrap().source().unwrap()
}
fn source_drive(s: &FrictionModalSystem) -> FrictionDrive {
    FrictionDrive {
        normal_force_n: s
            .adhesive_binding()
            .unwrap()
            .tangential()
            .contact()
            .repulsive_normal_force_n(),
        ..drive()
    }
}
fn native_bits(s: &FrictionModalSystem) -> Vec<u64> {
    states(s)
        .iter()
        .flat_map(|x| {
            [
                x.displacement_m_sqrt_kg.to_bits(),
                x.velocity_m_sqrt_kg_per_s.to_bits(),
            ]
        })
        .chain(s.material_memory().iter().map(|x| x.to_bits()))
        .chain(s.radiation_state().iter().map(|x| x.to_bits()))
        .chain([
            s.samples_rendered(),
            s.sample_period_s().to_bits(),
            s.total_energy_j().unwrap().to_bits(),
        ])
        .collect()
}
fn shear(normal: &BoundNormalInterfaceState) -> SourceAdhesiveShear {
    SourceAdhesiveShear::resolve(normal, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn jkr_original_normal_and_junction_temperatures_change_next_native_force_without_fake_heat() {
    let fixture = thermal_fixture::Fixture::with_options(10.0, 0.05, true, 320.0);
    let normal = fixture.jkr();
    let original_shear = shear(normal.normal_state());
    let force =
        -0.25 * PI * normal.effective_radius_m() * fixture.normal.adhesion_energy_j_per_m2();
    for held_load in [false, true] {
        for t in [300.0, 285.0, 315.0] {
            let refresh = if held_load {
                normal
                    .refresh_at_held_load(force, fixture.cards(), temperatures(t))
                    .unwrap()
            } else {
                normal
                    .refresh_at_held_indentation(
                        normal.evaluate_attached(0.0).unwrap().state(),
                        fixture.cards(),
                        temperatures(t),
                    )
                    .unwrap()
            };
            let source = original_shear
                .bind_jkr(normal.normal_state(), refresh.before())
                .unwrap();
            let mut s = bound_native(source.clone(), true);
            let mut unchanged = bound_native(source, true);
            for _ in 0..2 {
                s.step(&[0.0; 2], source_drive(&s)).unwrap();
                unchanged.step(&[0.0; 2], source_drive(&unchanged)).unwrap();
            }
            let bits = native_bits(&s);
            let last = s.last_frame().copied();
            let network_frame = s.coupled_frame().cloned();
            let port = s.port().clone();
            let old_binding = s.adhesive_binding().unwrap().clone();
            let report = s
                .apply_adhesive_source_refresh(
                    AdhesiveModalSourceRefresh::Jkr {
                        original: &normal,
                        refreshed: &refresh,
                    },
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            assert_eq!(native_bits(&s), bits);
            assert_eq!(s.last_frame().copied(), last);
            assert_eq!(s.coupled_frame(), network_frame.as_ref());
            assert_eq!(s.port(), &port);
            assert_eq!(report.sample, 2);
            assert_eq!(report.normal_control, refresh.control());
            assert_eq!(
                report.original_normal_state_identity,
                normal.normal_state().identity()
            );
            assert_eq!(
                report.normal_state_identity,
                refresh.source().normal_state().identity()
            );
            assert_eq!(
                report.contact_before,
                refresh.before().friction_contact_state().unwrap()
            );
            assert_eq!(
                report.contact_after,
                refresh.response().friction_contact_state().unwrap()
            );
            assert_eq!(report.temperature_before_k, Some(300.0));
            assert_eq!(report.temperature_after_k, Some(t));
            let accepted = source_binding(&s);
            assert_eq!(accepted.normal_state(), refresh.source().normal_state());
            assert_eq!(
                accepted.tangential(),
                s.adhesive_binding().unwrap().tangential()
            );
            assert_eq!(
                accepted.shear_source().resolved().state_coordinate("p"),
                Some(2e5)
            );
            for (before, after) in original_shear
                .resolved()
                .properties()
                .iter()
                .zip(accepted.shear_source().resolved().properties())
            {
                assert_eq!(
                    before.answer().receipt.selected,
                    after.answer().receipt.selected
                );
            }
            // Independent manufactured curve: tau_kinetic(T)=2000-20(T-300).
            let expected = -(2000.0 - 20.0 * (t - 300.0)) * report.contact_after.area_m2();
            let next = *s.step(&[0.0; 2], source_drive(&s)).unwrap();
            let old = *unchanged.step(&[0.0; 2], source_drive(&unchanged)).unwrap();
            assert_eq!(next.regime, ModalFrictionRegime::Sliding);
            close(next.traction_n, expected, 1e-11);
            close(
                next.friction_dissipation_j,
                next.traction_n * next.slip_distance_m,
                1e-18,
            );
            assert!(next.energy_residual_j.abs() <= next.energy_tolerance_j);
            if t == 300.0 {
                assert_eq!(s.adhesive_binding().unwrap(), &old_binding);
                assert_eq!(next, old);
                assert_eq!(native_bits(&s), native_bits(&unchanged));
            } else {
                assert_ne!(refresh.parameter_work_j(), 0.0);
                if held_load {
                    assert_ne!(refresh.mechanical_work_j(), 0.0);
                } else {
                    assert_eq!(refresh.mechanical_work_j(), 0.0);
                }
                assert_ne!(next.traction_n.to_bits(), old.traction_n.to_bits());
                assert_ne!(native_bits(&s), native_bits(&unchanged));
            }
        }
    }
}

fn dmt(fixture: &thermal_fixture::Fixture) -> SourceDmtSphereContact {
    SourceDmtSphereContact::new(
        fixture.normal.clone(),
        0.01,
        1e-6,
        DmtEnvelope {
            maximum_indentation_m: 1e-6,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}
fn maugis(fixture: &thermal_fixture::Fixture) -> SourceMaugisSphereContact {
    let radius = 0.01;
    let a0 = (PI * fixture.normal.adhesion_energy_j_per_m2() * radius * radius
        / (4.0 * fixture.normal.elastic().reduced_modulus_pa() / 3.0))
        .cbrt();
    let d0 = a0 * a0 / radius;
    SourceMaugisSphereContact::new(
        fixture.normal.clone(),
        radius,
        2.0 * d0 / (PI * 0.2),
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}
fn check_finite_refresh(
    source: SourceAdhesiveTangential,
    refresh: AdhesiveModalSourceRefresh<'_>,
    successor: &BoundNormalInterfaceState,
    contact: AdhesiveContactState,
    control: AdhesiveThermalControl,
) {
    let mut s = bound_native(source, false);
    s.step(&[0.0; 2], source_drive(&s)).unwrap();
    let bits = native_bits(&s);
    let last = s.last_frame().copied();
    let network_frame = s.coupled_frame().cloned();
    let report = s
        .apply_adhesive_source_refresh(refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(native_bits(&s), bits);
    assert_eq!(s.last_frame().copied(), last);
    assert_eq!(s.coupled_frame(), network_frame.as_ref());
    assert_eq!(report.normal_control, control);
    assert_eq!(report.contact_after, contact);
    assert_eq!(source_binding(&s).normal_state(), successor);
    assert_eq!(source_binding(&s).temperature_kelvin(), Some(315.0));
    assert_eq!(
        source_binding(&s).tangential(),
        s.adhesive_binding().unwrap().tangential()
    );
    let p = momentum(&s);
    let next = *s.step(&[0.0; 2], source_drive(&s)).unwrap();
    close(next.traction_n, -1700.0 * contact.area_m2(), 1e-12);
    close(momentum(&s), p, 1e-14);
    assert!(next.energy_residual_j.abs() <= next.energy_tolerance_j);
    if contact.area_m2() == 0.0 {
        assert_eq!(next.regime, ModalFrictionRegime::Released);
        assert_eq!(next.friction_dissipation_j, 0.0);
    } else {
        assert_eq!(next.regime, ModalFrictionRegime::Sliding);
        assert!(next.friction_dissipation_j > 0.0);
    }
}

#[test]
fn finite_range_thermal_sources_preserve_true_solid_area_and_external_normal_control() {
    let fixture = thermal_fixture::Fixture::with_options(1e5, 0.05, false, 320.0);
    let normal = dmt(&fixture);
    let source = shear(normal.normal_state());
    for delta in [0.2e-6, -0.5e-6, -2e-6] {
        let refresh = normal
            .refresh_at_held_indentation(delta, fixture.cards(), temperatures(315.0))
            .unwrap();
        check_finite_refresh(
            source
                .bind_dmt(normal.normal_state(), refresh.before())
                .unwrap(),
            AdhesiveModalSourceRefresh::Dmt {
                original: &normal,
                refreshed: &refresh,
            },
            refresh.source().normal_state(),
            refresh.response().friction_contact_state().unwrap(),
            refresh.control(),
        );
    }
    let force = normal.evaluate(0.2e-6).unwrap().normal_force_n;
    let refresh = normal
        .refresh_at_held_load(force, fixture.cards(), temperatures(315.0))
        .unwrap();
    assert_ne!(refresh.mechanical_work_j(), 0.0);
    check_finite_refresh(
        source
            .bind_dmt(normal.normal_state(), refresh.before())
            .unwrap(),
        AdhesiveModalSourceRefresh::Dmt {
            original: &normal,
            refreshed: &refresh,
        },
        refresh.source().normal_state(),
        refresh.response().friction_contact_state().unwrap(),
        refresh.control(),
    );

    let fixture = thermal_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let normal = maugis(&fixture);
    let source = shear(normal.normal_state());
    let bounds = normal.noncontact_limits().unwrap().unwrap();
    for before in [
        normal.evaluate_contact_energy(0.0).unwrap(),
        normal
            .evaluate_noncontact(
                0.25 * bounds.minimum_indentation_m + 0.75 * bounds.first_contact_indentation_m,
            )
            .unwrap(),
        normal
            .evaluate_detached(-2.0 * normal.interaction_range_m())
            .unwrap(),
    ] {
        let refresh = normal
            .refresh_at_held_indentation(before, fixture.cards(), temperatures(315.0))
            .unwrap();
        assert_eq!(refresh.response().phase(), before.phase());
        check_finite_refresh(
            source
                .bind_maugis_equilibrium(&normal, refresh.before())
                .unwrap(),
            AdhesiveModalSourceRefresh::Maugis {
                original: &normal,
                refreshed: &refresh,
            },
            refresh.source().normal_state(),
            refresh.response().friction_contact_state().unwrap(),
            refresh.control(),
        );
    }
    let force = 0.1 * PI * normal.effective_radius_m() * fixture.normal.adhesion_energy_j_per_m2();
    let refresh = normal
        .refresh_at_held_load(force, fixture.cards(), temperatures(315.0))
        .unwrap();
    assert_ne!(refresh.mechanical_work_j(), 0.0);
    check_finite_refresh(
        source
            .bind_maugis_equilibrium(&normal, refresh.before())
            .unwrap(),
        AdhesiveModalSourceRefresh::Maugis {
            original: &normal,
            refreshed: &refresh,
        },
        refresh.source().normal_state(),
        refresh.response().friction_contact_state().unwrap(),
        refresh.control(),
    );
}

#[test]
fn changed_temperature_keeps_the_occupied_sliding_branch_inside_the_new_static_cone() {
    let fixture = thermal_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let normal = fixture.jkr();
    let before = normal.evaluate_attached(0.0).unwrap();
    let source = shear(normal.normal_state())
        .bind_jkr(normal.normal_state(), &before)
        .unwrap();
    let mut s = bound_native(source, false);
    assert_eq!(
        s.step(&[0.0; 2], source_drive(&s)).unwrap().regime,
        ModalFrictionRegime::Sliding
    );
    let refresh = normal
        .refresh_at_held_indentation(before.state(), fixture.cards(), temperatures(315.0))
        .unwrap();
    s.apply_adhesive_source_refresh(
        AdhesiveModalSourceRefresh::Jkr {
            original: &normal,
            refreshed: &refresh,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let area = refresh.response().contact_area_m2;
    let static_capacity = 2550.0 * area;
    let kinetic_capacity = 1700.0 * area;
    // Exact free-body held-force mobility; the target lies between the new
    // kinetic and static capacities. A reset would incorrectly capture it.
    let mobility = 0.5 * DT * (1.0 / MASS[0] + 1.0 / MASS[1]);
    let drive = FrictionDrive {
        speed_m_s: relative_velocity(&s) - mobility * 0.5 * (static_capacity + kinetic_capacity),
        ..source_drive(&s)
    };
    let next = *s.step(&[0.0; 2], drive).unwrap();
    assert_eq!(next.regime, ModalFrictionRegime::Sliding);
    assert_eq!(
        (next.roots.unwrap().found, next.roots.unwrap().selected),
        (2, 1)
    );
    close(next.traction_n, -kinetic_capacity, 1e-12);
    assert!(next.friction_dissipation_j > 0.0);
    assert!(next.energy_residual_j.abs() <= next.energy_tolerance_j);
}

#[test]
fn stale_geometry_foreign_source_and_receiver_refusals_leave_every_native_owner_replayable() {
    let fixture = thermal_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let normal = fixture.jkr();
    let before = normal.evaluate_attached(0.0).unwrap();
    let source = shear(normal.normal_state())
        .bind_jkr(normal.normal_state(), &before)
        .unwrap();
    let mut s = bound_native(source.clone(), true);
    let mut replay = bound_native(source, true);
    s.step(&[0.0; 2], source_drive(&s)).unwrap();
    replay.step(&[0.0; 2], source_drive(&replay)).unwrap();
    let refresh = normal
        .refresh_at_held_indentation(before.state(), fixture.cards(), temperatures(315.0))
        .unwrap();
    let stale = normal
        .refresh_at_held_indentation(
            normal
                .evaluate_attached(0.2 * normal.indentation_scale_m())
                .unwrap()
                .state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let foreign = SourceJkrSphereContact::new(
        normal.normal_state().clone(),
        1.1 * normal.effective_radius_m(),
        normal.interaction_range_m(),
        normal.envelope(),
    )
    .unwrap();
    let bits = native_bits(&s);
    let last = s.last_frame().copied();
    let native = s.coupled_frame().cloned();
    let binding = s.adhesive_binding().unwrap().clone();
    for rejected in [
        AdhesiveModalSourceRefresh::Jkr {
            original: &normal,
            refreshed: &stale,
        },
        AdhesiveModalSourceRefresh::Jkr {
            original: &foreign,
            refreshed: &refresh,
        },
    ] {
        assert!(
            s.apply_adhesive_source_refresh(rejected, &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(native_bits(&s), bits);
        assert_eq!(s.adhesive_binding().unwrap(), &binding);
    }
    let request = AdhesiveModalSourceRefresh::Jkr {
        original: &normal,
        refreshed: &refresh,
    };
    let refused: Result<(), ModalFrictionError> = s.apply_adhesive_source_refresh_with_admission(
        request,
        &CancelGate::new_clock_free(),
        |report, candidate| {
            assert_eq!(
                candidate.source().unwrap().normal_state(),
                refresh.source().normal_state()
            );
            assert_eq!(
                report.contact_after,
                refresh.response().friction_contact_state().unwrap()
            );
            Err(ModalFrictionError::Law(
                fs_tribo::TriboError::InvalidInput {
                    field: "thermal receiver refused",
                },
            ))
        },
    );
    assert!(refused.is_err());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), ModalFrictionError> =
        s.apply_adhesive_source_refresh_with_admission(request, &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(cancelled.is_err());
    assert!(s.apply_adhesive_source_refresh(request, &gate).is_err());
    assert_eq!(native_bits(&s), bits);
    assert_eq!(s.adhesive_binding().unwrap(), &binding);
    assert_eq!(s.last_frame().copied(), last);
    assert_eq!(s.coupled_frame(), native.as_ref());
    let mut declared = mechanics(1.0, before.contact_area_m2);
    let declared_bits = native_bits(&declared);
    assert!(
        declared
            .apply_adhesive_source_refresh(request, &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(native_bits(&declared), declared_bits);
    s.apply_adhesive_source_refresh(request, &CancelGate::new_clock_free())
        .unwrap();
    replay
        .apply_adhesive_source_refresh(request, &CancelGate::new_clock_free())
        .unwrap();
    let accepted = s.adhesive_binding().unwrap().clone();
    assert!(
        s.apply_adhesive_source_refresh(request, &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(native_bits(&s), bits);
    assert_eq!(s.adhesive_binding().unwrap(), &accepted);
    assert_eq!(
        *s.step(&[0.0; 2], source_drive(&s)).unwrap(),
        *replay.step(&[0.0; 2], source_drive(&replay)).unwrap()
    );
    assert_eq!(native_bits(&s), native_bits(&replay));
}
