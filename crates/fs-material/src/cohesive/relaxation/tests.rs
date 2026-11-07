use super::*;

fn law() -> IntactViscoelasticBond {
    IntactViscoelasticBond::new(
        GeneralizedMaxwell::new(1e6, vec![(2e6, 0.1), (3e5, 0.8)]).unwrap(),
        GeneralizedMaxwell::new(4e5, vec![(8e5, 0.25)]).unwrap(),
        0.001,
        ViscoelasticBondLimits {
            maximum_jump_m: 0.0001,
            maximum_rate_m_s: 1.0,
            maximum_elapsed_s: 10.0,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

#[test]
fn held_separation_relaxes_with_independent_normal_and_shear_memory() {
    let law = law();
    let ramp = law.prepare(0.02).unwrap();
    let hold = law.prepare(0.3).unwrap();
    let gap = [2e-5, -1e-5, 0.5e-5];
    let first = law.bonded_point().trial(&ramp, gap, 1e-10).unwrap();
    let next = first.point.trial(&hold, gap, 1e-10).unwrap();
    assert_eq!(next.work_j_m2, 0.0);
    assert!(next.viscous_dissipated_j_m2 > 0.0);
    close(
        next.stored_change_j_m2,
        -next.viscous_dissipated_j_m2,
        2e-13,
    );
    for i in 0..3 {
        let s = &law.spectra[usize::from(i > 0)];
        let epsilon = gap[i] / law.thickness;
        let expected = s.e_inf * epsilon
            + s.terms
                .iter()
                .map(|&(e, tau)| {
                    e * epsilon * (tau / 0.02) * (1.0 - (-0.02 / tau).exp()) * (-0.3 / tau).exp()
                })
                .sum::<f64>();
        close(next.traction_pa[i], expected, 3e-13);
        assert!(next.traction_pa[i].abs() < first.traction_pa[i].abs());
    }
}

#[test]
fn complete_memory_work_and_heat_match_the_original_scalar_owner_on_reversals() {
    let law = law();
    let ramp = law.prepare(0.07).unwrap();
    let mut point = law.bonded_point();
    let mut originals =
        core::array::from_fn::<_, 3, _>(|i| law.spectra[usize::from(i > 0)].state());
    for gap in [
        [2e-5, 1e-5, -1e-5],
        [-1e-5, -2e-5, 1e-5],
        [0.0, 0.0, 0.0],
        [1e-5, -1e-5, -2e-5],
    ] {
        let step = point.trial(&ramp, gap, 1e-10).unwrap();
        let mut work = 0.0;
        let mut heat = 0.0;
        for i in 0..3 {
            let before = originals[i].clone();
            let stress = law.spectra[usize::from(i > 0)].step(
                &mut originals[i],
                gap[i] / law.thickness,
                0.07,
            );
            close(stress, step.traction_pa[i], 2e-13);
            work += (originals[i].work - before.work) * law.thickness;
            heat += (originals[i].dissipated - before.dissipated) * law.thickness;
        }
        close(work, step.work_j_m2, 2e-13);
        close(heat, step.viscous_dissipated_j_m2, 2e-13);
        close(
            step.point.stored_energy_j_m2() - point.stored_energy_j_m2(),
            step.stored_change_j_m2,
            2e-13,
        );
        point = step.point;
    }
}

#[test]
fn exact_ramp_subdivision_and_shear_rotation_preserve_the_same_physics() {
    let law = law();
    let start = law
        .bonded_point()
        .trial(&law.prepare(0.04).unwrap(), [2e-5, 1e-5, -1e-5], 1e-10)
        .unwrap()
        .point;
    let end = [-1e-5, -2e-5, 1.5e-5];
    let duration = 0.6;
    let whole = start
        .trial(&law.prepare(duration).unwrap(), end, 1e-10)
        .unwrap();
    let mut split = start.clone();
    let mut work = 0.0;
    let mut heat = 0.0;
    let small = law.prepare(duration / 40.0).unwrap();
    let old = start.jump_m();
    for j in 1..=40 {
        let gap = core::array::from_fn(|i| old[i] + (end[i] - old[i]) * f64::from(j) / 40.0);
        let step = split.trial(&small, gap, 1e-10).unwrap();
        work += step.work_j_m2;
        heat += step.viscous_dissipated_j_m2;
        split = step.point;
    }
    close(work, whole.work_j_m2, 3e-12);
    close(heat, whole.viscous_dissipated_j_m2, 3e-12);
    for i in 0..3 {
        close(split.traction_pa()[i], whole.traction_pa[i], 3e-12);
    }
    let rotate = |g: [f64; 3]| [g[0], 0.6 * g[1] - 0.8 * g[2], 0.8 * g[1] + 0.6 * g[2]];
    let rotated = law
        .bonded_point()
        .trial(&law.prepare(0.04).unwrap(), rotate(old), 1e-10)
        .unwrap()
        .point
        .trial(&law.prepare(duration).unwrap(), rotate(end), 1e-10)
        .unwrap();
    close(rotated.work_j_m2, whole.work_j_m2, 3e-13);
    close(
        rotated.viscous_dissipated_j_m2,
        whole.viscous_dissipated_j_m2,
        3e-13,
    );
    for i in 0..3 {
        close(
            rotated.mean_traction_pa[i],
            rotate(whole.mean_traction_pa)[i],
            3e-13,
        );
    }
}

#[test]
fn distinct_endpoint_and_mean_tangents_match_independent_perturbations() {
    let law = law();
    let ramp = law.prepare(0.12).unwrap();
    let point = law
        .bonded_point()
        .trial(&ramp, [1e-5, 2e-5, -1e-5], 1e-10)
        .unwrap()
        .point;
    let end = [-2e-5, 1e-5, 3e-5];
    assert_ne!(ramp.endpoint_stiffness, ramp.mean_stiffness);
    let held = point.held_mean_traction(&ramp).unwrap();
    let base = point.trial(&ramp, end, 1e-10).unwrap();
    for j in 0..3 {
        close(
            base.mean_traction_pa[j],
            held[j] + ramp.mean_stiffness[j] * (end[j] - point.jump_m()[j]),
            3e-13,
        );
        let mut plus = end;
        let mut minus = end;
        plus[j] += 1e-9;
        minus[j] -= 1e-9;
        let a = point.trial(&ramp, plus, 1e-10).unwrap();
        let b = point.trial(&ramp, minus, 1e-10).unwrap();
        for i in 0..3 {
            let expected = if i == j {
                ramp.endpoint_stiffness[i]
            } else {
                0.0
            };
            close(
                (a.traction_pa[i] - b.traction_pa[i]) / 2e-9,
                expected,
                3e-11,
            );
            let expected = if i == j { ramp.mean_stiffness[i] } else { 0.0 };
            close(
                (a.mean_traction_pa[i] - b.mean_traction_pa[i]) / 2e-9,
                expected,
                3e-11,
            );
        }
    }
}

#[test]
fn source_law_rate_range_and_time_refusals_leave_accepted_memory_retryable() {
    let law = law();
    let copy = law.clone();
    let point = law.bonded_point();
    assert!(
        point
            .trial(&copy.prepare(0.1).unwrap(), [0.0; 3], 1e-10)
            .is_err()
    );
    assert!(
        point
            .trial(&law.prepare(0.1).unwrap(), [1.01e-4, 0.0, 0.0], 1e-10)
            .is_err()
    );
    assert!(
        point
            .trial(&law.prepare(1e-8).unwrap(), [1e-5, 0.0, 0.0], 1e-10)
            .is_err()
    );
    let late = point
        .trial(&law.prepare(9.0).unwrap(), [1e-5, 0.0, 0.0], 1e-10)
        .unwrap()
        .point;
    assert!(
        late.trial(&law.prepare(2.0).unwrap(), [0.0; 3], 1e-10)
            .is_err()
    );
    assert_eq!(point.jump_m(), [0.0; 3]);
    assert_eq!(point.elapsed_s(), 0.0);
    let prepared = law.prepare(0.1).unwrap();
    let a = point.trial(&prepared, [1e-5, 0.0, 0.0], 1e-10).unwrap();
    let b = point.trial(&prepared, [1e-5, 0.0, 0.0], 1e-10).unwrap();
    assert_eq!(a.traction_pa, b.traction_pa);
    assert_eq!(a.work_j_m2, b.work_j_m2);
    assert_eq!(a.viscous_dissipated_j_m2, b.viscous_dissipated_j_m2);
}
