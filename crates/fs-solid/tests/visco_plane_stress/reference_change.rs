//! G0/G1: reference events preserve original thin-solid histories and volume actuation.
use super::*;
use fs_material::hyper::reference::ReferenceDeformation;

fn reference(scale: f64, shear: f64) -> ReferenceDeformation {
    let change = scale - 1.0;
    ReferenceDeformation::try_new([
        scale,
        shear,
        0.5 * shear,
        0.0,
        1.0 - 0.2 * change,
        -0.25 * shear,
        0.0,
        0.0,
        1.0 + 0.7 * change,
    ])
    .unwrap()
}

fn fields() -> ([ViscoHyperelastic; 2], [ViscoHyperelastic; 2]) {
    let first = law();
    let mut second = law().with_reference_deformation(reference(1.01, -0.01));
    second.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 6.0 };
    second.branches[0].tau_s = 0.35;
    let old = [first, second];
    let mut target = [
        old[0]
            .clone()
            .with_reference_deformation(reference(1.04, 0.03)),
        old[1]
            .clone()
            .with_reference_deformation(reference(0.98, -0.02)),
    ];
    target[0].branches[0].tau_s *= 0.5;
    target[1].branches[0].tau_s *= 1.4;
    (old, target)
}

fn shape(x: f64, y: f64) -> [f64; 2] {
    [0.11 * x + 0.06 * x * y, -0.025 * y + 0.04 * x * y]
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

#[test]
fn g0_heterogeneous_reference_events_preserve_q1_memory_and_account_for_constraint_actuation() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.4, 0.7, 2, 1);
        let (old, target) = fields();
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| shape(*x, *y)).collect();
        let thickness = 0.03;
        for exact in [false, true] {
            let body = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::PerElement(&old),
                mode(exact),
                thickness,
                &initial,
                local(),
                cx,
            )
            .unwrap()
            .advance(0.04, &boundary(&shape), &[], newton(), controls(), cx)
            .unwrap()
            .solid;
            let before = snapshot(&body);
            let event = body
                .change_reference(ViscoelasticLawField::PerElement(&target), controls(), cx)
                .unwrap();
            assert_eq!(snapshot(&body), before);
            assert_eq!(event.solid.time_s(), body.time_s());
            assert_eq!(event.solid.displacement_m(), body.displacement_m());
            assert_eq!(event.solid.material_points(), body.material_points());
            assert_eq!(event.solid.reference_thickness_m(), thickness);
            assert_eq!(event.solid.material_points().len(), 8);
            let mut stored = 0.0;
            let mut fixed_reference = 0.0;
            let mut native_loading = 0.0;
            let mut constraint_reference = 0.0;
            let mut transverse_work = 0.0;
            for (index, point) in body.material_points().iter().enumerate() {
                let volume = point.reference_area_m2 * thickness;
                let (old_f, new_f, point_stored, reference_increment, loading_increment) = if exact
                {
                    let previous = &body.point_remainder_states().unwrap()[index];
                    let next = &event.solid.point_remainder_states().unwrap()[index];
                    assert_eq!(
                        next.inverse_viscous_cauchy_green,
                        previous.inverse_viscous_cauchy_green
                    );
                    assert_eq!(next.dissipated, previous.dissipated);
                    let old_volume = old[point.element].reference_deformation().volume_ratio();
                    let new_volume = target[point.element].reference_deformation().volume_ratio();
                    close(determinant(&next.deformation_gradient), new_volume, 2e-13);
                    constraint_reference += volume
                        * event.solid.point_pressure_pa().unwrap()[index]
                        * (new_volume - old_volume);
                    (
                        previous.deformation_gradient,
                        next.deformation_gradient,
                        next.stored,
                        next.reference_work - previous.reference_work,
                        next.work - previous.work,
                    )
                } else {
                    let previous = &body.point_states().unwrap()[index];
                    let next = &event.solid.point_states().unwrap()[index];
                    assert_eq!(
                        next.inverse_viscous_cauchy_green,
                        previous.inverse_viscous_cauchy_green
                    );
                    assert_eq!(next.dissipated, previous.dissipated);
                    (
                        previous.deformation_gradient,
                        next.deformation_gradient,
                        next.stored,
                        next.reference_work - previous.reference_work,
                        next.work - previous.work,
                    )
                };
                for component in [0, 1, 3, 4] {
                    assert_eq!(new_f[component], old_f[component]);
                }
                for component in [2, 5, 8] {
                    let stress = event.solid.point_first_piola_pa()[index][component];
                    close(stress, 0.0, 1e-8);
                    transverse_work += volume * stress * (new_f[component] - old_f[component]);
                }
                stored += volume * point_stored;
                fixed_reference += volume * reference_increment;
                native_loading += volume * loading_increment;
            }
            close(event.solid.stored_energy_j(), stored, 2e-12);
            close(event.energy.stored_before_j, body.stored_energy_j(), 2e-12);
            close(event.energy.stored_after_j, stored, 2e-12);
            close(event.energy.reference_work_j, fixed_reference, 2e-12);
            close(
                event.energy.native_transverse_loading_work_j,
                native_loading,
                2e-12,
            );
            close(
                event.energy.constraint_reference_work_j,
                constraint_reference,
                2e-12,
            );
            close(
                event.energy.transverse_endpoint_work_j,
                transverse_work,
                2e-12,
            );
            close(event.energy.material_balance_error_j, 0.0, 2e-12);
            close(
                event.energy.stored_after_j - event.energy.stored_before_j,
                event.energy.reference_work_j + constraint_reference + transverse_work
                    - event.energy.endpoint_work_defect_j,
                2e-12,
            );
            if exact {
                assert!(constraint_reference.abs() > 1e-8);
                let histories = event.solid.point_remainder_states().unwrap();
                assert_ne!(
                    histories[0].inverse_viscous_cauchy_green,
                    histories[1].inverse_viscous_cauchy_green
                );
            } else {
                assert_eq!(event.energy.constraint_reference_work_j, 0.0);
                let histories = event.solid.point_states().unwrap();
                assert_ne!(
                    histories[0].inverse_viscous_cauchy_green,
                    histories[1].inverse_viscous_cauchy_green
                );
            }
        }
    });
}

#[test]
fn g0_late_reference_gate_and_second_element_refusal_preserve_the_complete_accepted_body() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.4, 0.7, 2, 1);
        let (old, target) = fields();
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| shape(*x, *y)).collect();
        for exact in [false, true] {
            let body = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::PerElement(&old),
                mode(exact),
                0.03,
                &initial,
                local(),
                cx,
            )
            .unwrap();
            let before = snapshot(&body);
            let proposal = body
                .change_reference(ViscoelasticLawField::PerElement(&target), controls(), cx)
                .unwrap();
            let limit = 0.5 * proposal.energy.endpoint_work_defect_j.abs();
            assert!(limit > 1e-10);
            let strict = FiniteViscoStepControls {
                maximum_endpoint_work_defect_j: limit,
                ..controls()
            };
            match body.change_reference(ViscoelasticLawField::PerElement(&target), strict, cx) {
                Err(FiniteViscoPlaneError::EnergyGate {
                    value_j,
                    maximum_abs_j,
                    ..
                }) => {
                    close(value_j, proposal.energy.endpoint_work_defect_j, 2e-12);
                    assert_eq!(maximum_abs_j, limit);
                }
                other => panic!("expected late aggregate reference-work gate, got {other:?}"),
            }
            assert_eq!(snapshot(&body), before);
            let mut incompatible = target.clone();
            incompatible[1].branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 7.0 };
            assert!(matches!(
                body.change_reference(
                    ViscoelasticLawField::PerElement(&incompatible),
                    controls(),
                    cx
                ),
                Err(FiniteViscoPlaneError::Material {
                    element: 1,
                    local_point: 0,
                    ..
                })
            ));
            assert_eq!(snapshot(&body), before);
            assert!(matches!(
                body.change_reference(
                    ViscoelasticLawField::PerElement(&target[..1]),
                    controls(),
                    cx
                ),
                Err(FiniteViscoPlaneError::Native(
                    SolidError::InvalidInput { .. }
                ))
            ));
            assert_eq!(snapshot(&body), before);
        }
    });
}

#[test]
fn g1_positive_time_after_reference_event_uses_new_laws_without_resetting_memory_or_clock() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.4, 0.7, 2, 1);
        let (old, target) = fields();
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| shape(*x, *y)).collect();
        for exact in [false, true] {
            let body = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::PerElement(&old),
                mode(exact),
                0.03,
                &initial,
                local(),
                cx,
            )
            .unwrap()
            .advance(0.04, &boundary(&shape), &[], newton(), controls(), cx)
            .unwrap()
            .solid;
            let changed = body
                .change_reference(ViscoelasticLawField::PerElement(&target), controls(), cx)
                .unwrap()
                .solid;
            let before = snapshot(&changed);
            let dt = 0.06;
            let resumed = changed
                .advance(dt, &boundary(&shape), &[], newton(), controls(), cx)
                .unwrap();
            assert_eq!(snapshot(&changed), before);
            close(resumed.solid.time_s(), 0.04 + dt, 1e-15);
            assert!(resumed.energy.dissipated_j > 0.0);
            check_energy(&resumed.energy);
            for (index, point) in resumed.solid.material_points().iter().enumerate() {
                let native = &target[point.element];
                if exact {
                    let previous = &changed.point_remainder_states().unwrap()[index];
                    let next = &resumed.solid.point_remainder_states().unwrap()[index];
                    let f = previous.deformation_gradient;
                    let expected = native
                        .step_incompressible_plane_stress_response(
                            previous,
                            &point.in_plane,
                            dt,
                            Some([f[2], f[5]]),
                            local(),
                        )
                        .unwrap();
                    assert_eq!(*next, expected.material.state);
                    assert_eq!(next.reference_work, previous.reference_work);
                    assert_ne!(
                        next.inverse_viscous_cauchy_green,
                        previous.inverse_viscous_cauchy_green
                    );
                    close(
                        resumed.solid.point_pressure_pa().unwrap()[index],
                        expected.plane_stress.pressure_pa,
                        2e-12,
                    );
                } else {
                    let previous = &changed.point_states().unwrap()[index];
                    let next = &resumed.solid.point_states().unwrap()[index];
                    let f = previous.deformation_gradient;
                    let expected = native
                        .step_plane_stress_response(
                            previous,
                            &point.in_plane,
                            dt,
                            Some([f[2], f[5], f[8]]),
                            local(),
                        )
                        .unwrap();
                    assert_eq!(*next, expected.material.state);
                    assert_eq!(next.reference_work, previous.reference_work);
                    assert_ne!(
                        next.inverse_viscous_cauchy_green,
                        previous.inverse_viscous_cauchy_green
                    );
                }
            }
        }
    });
}
