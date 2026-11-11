//! G1/G3/G4: sampled strength over the complete native affine load ball.
use super::*;
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::AffineLoadSpaceConstraints3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::stress::AffineStressLimit3;
use fs_topopt::sdf3::stress::StressOptions3;
use std::cell::Cell;

fn strength(limit: f64) -> AffineStressLimit3 {
    AffineStressLimit3 {
        limit,
        options: StressOptions3 {
            // Both direct density coefficients matter at these densities.
            relaxation_power: 2.0,
            aggregation_power: 2.0,
            ..Default::default()
        },
    }
}

fn constraints(stress: AffineStressLimit3) -> AffineLoadSpaceConstraints3<'static> {
    AffineLoadSpaceConstraints3 {
        stress: Some(stress),
        ..Default::default()
    }
}

fn evaluate(
    study: &mut CutDensityStudy3,
    rho: &[f64],
    dead: &[f64],
    basis: &[&[f64]],
    stress: AffineStressLimit3,
    control: &mut SolveControl<'_>,
) -> AffineLoadSpaceEvaluation3 {
    study
        .evaluate_robust_affine_constraints(
            rho,
            dead,
            basis,
            RobustOptions3::default(),
            &policy(),
            constraints(stress),
            control,
        )
        .unwrap()
}

fn von_mises(s: &[f64; 6]) -> f64 {
    (0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
        + 3.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2)))
    .sqrt()
}

#[test]
fn affine_stress_envelope_uses_every_native_point_and_the_complete_material_gradient() {
    let (mut study, dead, first, second) = fixture();
    let rho = raw(&study);
    let basis = [&first[..], &second[..]];
    let params = study.params();
    let scales = study.operator().scales().to_vec();
    let limit = strength(0.5);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let ordinary = study
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let ordinary_solves = control.work().linear_solves;
    let full = evaluate(&mut study, &rho, &dead, &basis, limit, &mut control);
    let stress = full.stress.as_ref().unwrap();
    assert_eq!(stress.adjoints.len(), basis.len() + 1);
    assert_eq!(
        control.work().linear_solves - ordinary_solves,
        ordinary_solves + basis.len() + 2
    );
    close_vector(&full.offset_displacement, &ordinary.offset_displacement);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(stress.point_count, stress.witnesses.len());
    assert_eq!(stress.point_fields, stress.point_count * (basis.len() + 1));
    close(
        stress.constraint,
        stress.aggregate / limit.limit - 1.0,
        1e-12,
    );

    let mut native = operator();
    native.set_scales(&full.physical.eroded_scales).unwrap();
    let fields: Vec<_> = std::iter::once(&full.offset_displacement)
        .chain(&full.basis_displacements)
        .map(|u| {
            native
                .bulk_stress(u, limit.options.max_points, poll)
                .unwrap()
        })
        .collect();
    assert!(fields.iter().all(|field| field.len() == stress.point_count));
    let mut relaxed_max = 0.0_f64;
    let mut physical_max = 0.0_f64;
    let mut sum = 0.0;
    let mut relaxed_direct = 0.0;
    let mut physical_direct = 0.0;
    let mut direction_separation = 0.0_f64;
    for (index, witness) in stress.witnesses.iter().enumerate() {
        let point = &fields[0][index];
        assert_eq!(witness.cell, point.cell);
        assert_eq!(witness.position, point.position);
        assert_eq!(witness.weight, point.weight);
        let coefficients = &witness.oracle.coefficients;
        close(dot(coefficients, coefficients), 1.0, 1e-9);
        let tensor = |a: &[f64]| {
            std::array::from_fn(|j| {
                fields[0][index].reference_stress[j]
                    + a.iter()
                        .enumerate()
                        .map(|(k, a)| a * fields[k + 1][index].reference_stress[j])
                        .sum::<f64>()
            })
        };
        let reference = tensor(coefficients);
        close_vector(&witness.reference_stress, &reference);
        let vm = von_mises(&reference);
        close(witness.reference_von_mises, vm, 1e-9);
        close(witness.worst_reference_von_mises, vm, 1e-7);
        let r = full.physical.eroded_rho[point.cell];
        let k = full.physical.eroded_scales[point.cell];
        let relaxed = r.powf(limit.options.relaxation_power) * vm;
        let physical = k * vm;
        close(witness.relaxed_von_mises, relaxed, 1e-9);
        close(witness.physical_von_mises, physical, 1e-9);
        close(witness.worst_relaxed_von_mises, relaxed, 1e-7);
        close(witness.worst_physical_von_mises, physical, 1e-7);
        relaxed_max = relaxed_max.max(relaxed);
        physical_max = physical_max.max(physical);
        sum += relaxed.powi(2) + physical.powi(2);
        if r > 0.0 && r < 1.0 {
            relaxed_direct += relaxed / stress.aggregate * vm * 2.0 * r;
            physical_direct += physical / stress.aggregate
                * vm
                * (1.0 - params.e_min)
                * params.penal
                * r.powf(params.penal - 1.0);
        }
        direction_separation =
            direction_separation.max((1.0 - dot(coefficients, &full.worst_coefficients)).abs());
        // This is a dense numerical diagnostic. The production adversary comes
        // from its quadratic ball oracle, independently at EVERY retained point.
        let mut sampled = 0.0_f64;
        for j in 0..720 {
            let angle = std::f64::consts::TAU * j as f64 / 720.0;
            sampled = sampled.max(von_mises(&tensor(&[angle.cos(), angle.sin()])));
        }
        assert!(sampled <= witness.worst_reference_von_mises * (1.0 + 1e-8));
        close(sampled, witness.worst_reference_von_mises, 2e-4);
    }
    close(stress.aggregate, sum.sqrt(), 1e-7);
    close(stress.sampled_relaxed_max, relaxed_max, 1e-9);
    close(stress.sampled_physical_max, physical_max, 1e-9);
    assert!(stress.aggregate >= relaxed_max && stress.aggregate >= physical_max);
    assert!(relaxed_direct > 1e-4 * stress.aggregate);
    assert!(physical_direct > 1e-4 * stress.aggregate);
    assert!(direction_separation > 1e-4);

    let gradient = stress.gradient.as_ref().unwrap();
    // Holding the present native body force fixed preserves every current field
    // and stress, but removes the baseline adjoint's material-load pullback.
    let mut frozen = CutDensityStudy3::new(operator(), 0.15, params)
        .with_physical_regions(study.physical_regions().unwrap().to_vec())
        .unwrap();
    let frozen = evaluate(
        &mut frozen,
        &rho,
        &full.offset_force,
        &basis,
        limit,
        &mut control,
    );
    let frozen = frozen.stress.as_ref().unwrap();
    close(stress.aggregate, frozen.aggregate, 1e-9);
    let difference = gradient
        .iter()
        .zip(frozen.gradient.as_ref().unwrap())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let gradient_scale = gradient.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
    assert!(difference > 1e-4 * gradient_scale);
    for direction in [
        vec![1.0; rho.len()],
        (0..rho.len()).map(|i| (i % 3) as f64 * 0.5 - 0.5).collect(),
    ] {
        let h = 2e-5;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = evaluate(&mut study, &plus, &dead, &basis, limit, &mut control);
        let b = evaluate(&mut study, &minus, &dead, &basis, limit, &mut control);
        let fd =
            (a.stress.unwrap().aggregate - b.stress.unwrap().aggregate) / (2.0 * h * limit.limit);
        close(fd, dot(gradient, &direction), 3e-4);
    }
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn affine_mass_archive_and_fresh_recheck_enforce_strength_response_and_compliance_together() {
    let (mut study, dead, first, second) = fixture();
    let rho = vec![0.7; study.cells()];
    let basis = [&first[..], &second[..]];
    let scales = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = evaluate(&mut study, &rho, &dead, &basis, strength(1.0), &mut control);
    let baseline_stress = baseline.stress.as_ref().unwrap();
    let cap = AffineLoadSpace3::new(basis.len(), 2.0 * baseline.worst_compliance).unwrap();
    let q: Vec<_> = first.iter().map(|f| -f).collect();
    let center = dot(&q, &baseline.offset_displacement);
    let columns: Vec<_> = baseline
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let lower = center - dot(&columns, &columns).sqrt();
    assert!(lower < 0.0);
    let limits = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&q),
        sense: ResponseLimitSense3::AtLeast,
        bound: 2.0 * lower,
        scale: lower.abs(),
    }];
    let responses = AffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    let mut requirements = AffineLoadSpaceConstraints3 {
        responses: Some(responses),
        stress: Some(strength(0.5 * baseline_stress.sampled_physical_max)),
        ..Default::default()
    };
    let options = AffineLoadSpaceDesignOptions3 {
        density_floor: 0.6,
        ..Default::default()
    };
    {
        let mut rejected = AffineLoadSpaceDesign3::new_with_constraints(
            &mut study,
            &dead,
            &basis,
            &rho,
            cap.clone(),
            mass_goal(),
            options,
            requirements,
            &mut control,
        )
        .unwrap();
        let e = &rejected.accepted().evaluation;
        assert!(e.maximum <= 1.0);
        assert!(e.responses.as_ref().unwrap().physically_feasible());
        assert!(!e.stress.as_ref().unwrap().physically_feasible());
        assert!(rejected.best_feasible().is_none());
        assert!(rejected.recheck_best().unwrap().is_none());
    }
    let admitted = strength(2.0 * baseline_stress.aggregate);
    requirements.stress = Some(admitted);
    let mut design = AffineLoadSpaceDesign3::new_with_constraints(
        &mut study,
        &dead,
        &basis,
        &rho,
        cap.clone(),
        mass_goal(),
        options,
        requirements,
        &mut control,
    )
    .unwrap();
    assert!(design.gradient_check().passed());
    let start_mass = design
        .accepted()
        .evaluation
        .physical
        .dilated_mass_kg
        .unwrap();
    design.run(4).unwrap();
    let best = design.best_feasible().unwrap();
    assert!(best.evaluation.physical.dilated_mass_kg.unwrap() < start_mass);
    let best_rho = best.rho.clone();
    let best_solves = best.evaluation.work.linear_solves;
    let history = design.history().len();
    let optimizer_work = design.optimizer_work();
    let checked = design.recheck_best().unwrap().unwrap();
    assert_eq!(checked.rho, best_rho);
    assert_eq!(design.history().len(), history);
    assert_eq!(design.optimizer_work(), optimizer_work);
    assert!(checked.evaluation.work.linear_solves > best_solves);
    let stress = checked.evaluation.stress.as_ref().unwrap();
    assert_eq!(stress.limit.limit, admitted.limit);
    assert!(stress.physically_feasible());
    assert!(stress.aggregate <= admitted.limit);
    assert!(stress.sampled_relaxed_max <= admitted.limit);
    assert!(stress.sampled_physical_max <= admitted.limit);
    assert!(
        checked
            .evaluation
            .responses
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(checked.evaluation.worst_compliance <= cap.reference());
    audit(&checked.evaluation, &dead, &basis);
    let inventory: f64 = design
        .study()
        .operator()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&checked.evaluation.physical.dilated_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(
        checked.evaluation.physical.dilated_mass_kg.unwrap(),
        inventory,
        1e-12,
    );
    drop(design);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn affine_stress_counts_all_point_fields_and_restores_material_state_on_stops() {
    let mut native = operator();
    let scales: Vec<_> = (0..native.scales().len())
        .map(|i| 0.6 + 0.02 * i as f64)
        .collect();
    native.set_scales(&scales).unwrap();
    let (mut study, dead, first, second) = fixture_on_operator(native);
    let rho = raw(&study);
    let basis = [&first[..], &second[..]];
    let fingerprint = study.material_body_load_fingerprint();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let limit = strength(0.5);
    let full = evaluate(&mut study, &rho, &dead, &basis, limit, &mut control);
    let stress = full.stress.as_ref().unwrap();
    assert_eq!(stress.point_fields, stress.point_count * (basis.len() + 1));
    let mut limited = limit;
    limited.options.max_points = stress.point_fields - 1;
    assert!(matches!(
        study.evaluate_robust_affine_constraints(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            constraints(limited),
            &mut control,
        ),
        Err(EvaluationStop::Breakdown {
            stage: "sdf3-affine-stress-point-budget"
        })
    ));
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.material_body_load_fingerprint(), fingerprint);
    limited.options.max_points = stress.point_fields;
    let exact = evaluate(&mut study, &rho, &dead, &basis, limited, &mut control);
    let exact = exact.stress.as_ref().unwrap();
    close(exact.aggregate, stress.aggregate, 1e-12);
    close_vector(
        exact.gradient.as_ref().unwrap(),
        stress.gradient.as_ref().unwrap(),
    );
    for stage in ["sdf3-affine-stress-adjoint", "sdf3-affine-stress-publish"] {
        let reached = Cell::new(false);
        let mut stop = |progress: SolveProgress| {
            if progress.stage == stage {
                reached.set(true);
                ControlFlow::Break(())
            } else {
                poll()
            }
        };
        let mut stopped = SolveControl::new(SolveBudget::default(), &mut stop);
        assert!(matches!(
            study.evaluate_robust_affine_constraints(
                &rho,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                constraints(limit),
                &mut stopped,
            ),
            Err(EvaluationStop::Cancelled)
        ));
        assert!(reached.get());
        assert!(stopped.work().linear_iterations > 0);
        assert_eq!(study.operator().scales(), scales);
        assert_eq!(study.material_body_load_fingerprint(), fingerprint);
    }
    let replay = evaluate(&mut study, &rho, &dead, &basis, limit, &mut control);
    let replay = replay.stress.as_ref().unwrap();
    close(replay.aggregate, stress.aggregate, 1e-12);
    close_vector(
        replay.gradient.as_ref().unwrap(),
        stress.gradient.as_ref().unwrap(),
    );
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.material_body_load_fingerprint(), fingerprint);
}
