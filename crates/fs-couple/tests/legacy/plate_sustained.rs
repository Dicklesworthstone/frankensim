//! G0/G1/G3/G4: original midpoint trajectories under sustained normal contact.

use fs_contact::multicontact::MultiImpactConfig;
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::SolveControls;
use fs_couple::plate_impact::PlateRigidImpactBudget;
use fs_couple::plate_impact::multiple::{PlateRigidImpactSetBudget, PlateRigidImpactSetControls};
use fs_couple::plate_sustained::{
    PlateSustainedBudget, PlateSustainedControls, PlateSustainedError, PlateSustainedLoads,
    PlateSustainedPair, PlateSustainedReport, PlateSustainedSite, PlateSustainedSourceError,
    PlateSustainedSurface,
};
use fs_la::factor::lu;
use fs_plate::loading::PlateLoadBudget;
use fs_plate::{PlateMesh, PlateModel};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateLoad,
};
use fs_sparse::{Coo, Csr};
use fs_time::galpha::impulse::ImpulseControls;
use fs_time::galpha::initialization::InitialSolveConfig;
use fs_time::galpha::{ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState};
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
struct Fixture {
    mesh: PlateMesh,
    model: PlateModel,
    mass: Csr,
    stiffness: Csr,
}

// Independent exact P1 area integrals on the unit right triangle: M_ii=1/2,
// M_ij=1/4 and total transverse mass3. The rank-one elastic strain annihilates
// rigid translation and tilt, so stiffness cannot invent external momentum.
fn fixture(supported: bool, kappa: f64) -> Fixture {
    let mesh = PlateMesh {
        nodes: vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        tris: vec![[0, 1, 2]],
    };
    let mut mass = Coo::new(9, 9);
    for i in 0..3 {
        for j in 0..3 {
            for c in 0..3 {
                mass.push(
                    3 * i + c,
                    3 * j + c,
                    (if c == 0 { 1.0 } else { 1.0 / 12.0 }) * if i == j { 0.5 } else { 0.25 },
                );
            }
        }
    }
    let mass = mass.assemble();
    let strain = [-1.0, -0.5, 0.0, 1.0, -0.5, 0.0, 0.0, 0.0, 0.0];
    let mut stiffness = Coo::new(9, 9);
    for i in 0..9 {
        for j in 0..9 {
            if strain[i] != 0.0 && strain[j] != 0.0 && kappa != 0.0 {
                stiffness.push(i, j, kappa * strain[i] * strain[j]);
            }
        }
    }
    let stiffness = stiffness.assemble();
    let mut free = 0;
    let map: Vec<_> = (0..9)
        .map(|i| {
            if supported && i < 3 {
                None
            } else {
                let r = free;
                free += 1;
                Some(r)
            }
        })
        .collect();
    let reduce = |matrix: &Csr| {
        let mut reduced = Coo::new(free, free);
        for i in 0..9 {
            if let Some(r) = map[i] {
                let (cols, values) = matrix.row(i);
                for (&j, &v) in cols.iter().zip(values) {
                    if let Some(c) = map[j] {
                        reduced.push(r, c, v);
                    }
                }
            }
        }
        reduced.assemble()
    };
    let model = PlateModel {
        m: reduce(&mass),
        k: reduce(&stiffness),
        dof_map: map,
        free,
    };
    Fixture {
        mesh,
        model,
        mass,
        stiffness,
    }
}

struct HeldLoad<'a> {
    full: &'a [f64],
    map: &'a [Option<usize>],
    refuse: Cell<bool>,
}
impl PlateLoad for HeldLoad<'_> {
    fn parameter_count(&self) -> usize {
        0
    }
    fn forcing(&self, _: f64, output: &mut [f64]) -> Result<(), String> {
        if self.refuse.get() {
            return Err("original source refused".into());
        }
        output.fill(0.0);
        for (i, &r) in self.map.iter().enumerate() {
            if let Some(r) = r {
                output[r] = self.full[i];
            }
        }
        Ok(())
    }
    fn forcing_vjp(&self, _: f64, _: &[f64], output: &mut [f64]) -> Result<(), String> {
        output.fill(0.0);
        Ok(())
    }
}
fn surface<'a, L: PlateLoad + ?Sized>(
    f: &'a Fixture,
    load: &'a L,
    mass: f64,
    alpha: f64,
    beta: f64,
) -> PlateSustainedSurface<'a, L> {
    PlateSustainedSurface {
        dynamics: PlateDynamics::new(
            &f.model,
            PlateDynamicsParameters {
                stiffness_scale: 1.0,
                mass_scale: mass,
                mass_damping_per_s: alpha,
                stiffness_damping_s: beta,
            },
            load,
            PlateDynamicsBudget {
                max_dofs: 9,
                max_nonzeros: 100,
                max_full_dofs: 9,
                max_load_parameters: 0,
            },
            &mut || false,
        )
        .unwrap(),
        mesh: &f.mesh,
        full_mass: &f.mass,
        full_stiffness: &f.stiffness,
        reference_height_m: 0.0,
    }
}
fn controls() -> PlateSustainedControls {
    let mut step = ImplicitSolveConfig::default();
    step.newton.absolute_tolerance = 1e-13;
    step.newton.relative_tolerance = 0.0;
    step.newton.linear_restart = 9;
    step.newton.max_linear_cycles = 4;
    let mut contact = SolveControls::default();
    contact.newton.absolute_tolerance = 1e-11;
    contact.newton.relative_tolerance = 0.0;
    PlateSustainedControls {
        mechanics: PlateRigidImpactSetControls {
            impulse: ImpulseControls {
                mass_solve: InitialSolveConfig {
                    restart: 9,
                    max_cycles: 4,
                    tolerance: 1e-13,
                },
                absolute_impulse_tolerance: 1e-11,
                relative_impulse_tolerance: 1e-11,
                work_tolerance: 1e-10,
            },
            contact: MultiImpactConfig::default(),
            point_tolerance_m: 1e-12,
            speed_tolerance_m_s: 1e-10,
            linear_impulse_tolerance_ns: 1e-10,
            angular_impulse_tolerance_nms: 1e-10,
            energy_tolerance_j: 1e-9,
        },
        step,
        contact,
        gap_tolerance_m: 1e-10,
        kinematic_tolerance: 1e-11,
    }
}
fn budget() -> PlateSustainedBudget {
    let owner = PlateRigidImpactSetBudget {
        plate: PlateRigidImpactBudget {
            location: PlateLoadBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            max_full_mass_nonzeros: 100,
            max_workspace_components: 100000,
        },
        max_contacts: 3,
        max_mobility_products: 1000,
    };
    PlateSustainedBudget {
        positive: owner,
        negative: owner,
        max_stiffness_nonzeros: 200,
        max_workspace_components: 1000000,
        max_history_entries: 64,
    }
}
fn sites(points: &[(f64, f64)]) -> Vec<PlateSustainedSite> {
    points
        .iter()
        .map(|&(x, y)| PlateSustainedSite {
            point_xy: [x, y],
            scales: NaturalMapScales {
                normal_n_per_m: 1e4,
                tangent_n_per_m: 1e4,
            },
        })
        .collect()
}
fn state(n: usize) -> SecondOrderState {
    SecondOrderState::new(0.0, &vec![0.0; n], &vec![0.0; n], &vec![0.0; n])
}
fn audit(report: &PlateSustainedReport) {
    for &v in &report.momentum_balance_residual {
        near(v, 0.0, 1e-10);
    }
    for e in [
        report.positive.energy_residual_j,
        report.negative.energy_residual_j,
        report.energy_residual_j,
    ] {
        near(e, 0.0, 1e-9);
    }
    near(
        report.positive.contact_work_j + report.negative.contact_work_j,
        0.0,
        1e-9,
    );
    assert!(
        report.positive.damping_midpoint_j >= -1e-12
            && report.negative.damping_midpoint_j >= -1e-12
    );
    for b in &report.gap_bounds_m {
        assert!(b[0] >= -1e-10 && b[0] <= b[1]);
    }
}

#[test]
fn loaded_closed_contact_repeats_on_original_native_timelines_and_restarts_stale_acceleration() {
    let f = fixture(false, 0.0);
    let full = [-2.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0];
    let load = HeldLoad {
        full: &full,
        map: &f.model.dof_map,
        refuse: Cell::new(false),
    };
    let pair = PlateSustainedPair::prepare(
        surface(&f, &load, 1.0, 0.0, 0.0),
        surface(&f, &NoPlateLoad, 2.0, 0.0, 0.0),
        &sites(&f.mesh.nodes),
        controls(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut positive = state(9);
    let mut negative = state(9);
    positive.a.fill(17.0);
    negative.a.fill(-11.0);
    let h = 0.05;
    let mut forces = vec![0.0; 3];
    for step in 1..=5 {
        let before = positive.clone();
        let report = pair
            .step(
                &mut positive,
                &mut negative,
                h,
                PlateSustainedLoads::held(&full),
                PlateSustainedLoads::held(&[0.0; 9]),
                &forces,
                &mut || false,
            )
            .unwrap();
        audit(&report);
        let t = step as f64 * h;
        for r in 0..9 {
            for s in [&positive, &negative] {
                near(s.q[r], if r % 3 == 0 { -t * t / 3.0 } else { 0.0 }, 3e-10);
                near(s.v[r], if r % 3 == 0 { -2.0 * t / 3.0 } else { 0.0 }, 3e-10);
                near(s.a[r], if r % 3 == 0 { -2.0 / 3.0 } else { 0.0 }, 3e-9);
            }
        }
        for (i, force) in report.contact.forces_n.iter().enumerate() {
            near(force[0], 4.0 / 3.0, 3e-9);
            near(force[1], 0.0, 1e-12);
            near(force[2], 0.0, 1e-12);
            forces[i] = force[0];
            near(report.relative_normal_speeds_after[i], 0.0, 1e-10);
        }
        near(report.positive.external_impulse[0], -6.0 * h, 1e-12);
        near(
            report.positive.external_work_midpoint_j,
            -6.0 * (positive.q[0] - before.q[0]),
            1e-11,
        );
        for s in [&positive, &negative] {
            near(s.t, t, 1e-15);
            assert_eq!(s.steps, step);
            assert_eq!(s.history.len(), step);
            assert_eq!(s.history[step - 1].step, step - 1);
        }
    }
}

#[test]
fn common_motion_rayleigh_loss_matches_midpoint_and_second_order_refinement() {
    let f = fixture(false, 0.0);
    let mut errors = Vec::new();
    for steps in [4usize, 8, 16] {
        let h = 0.4 / steps as f64;
        let pair = PlateSustainedPair::prepare(
            surface(&f, &NoPlateLoad, 1.0, 0.9, 0.0),
            surface(&f, &NoPlateLoad, 2.0, 0.3, 0.0),
            &sites(&f.mesh.nodes),
            controls(),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut positive = state(9);
        let mut negative = state(9);
        for r in [0, 3, 6] {
            positive.v[r] = 1.0;
            negative.v[r] = 1.0;
        }
        let factor = (1.0 - 0.25 * h) / (1.0 + 0.25 * h);
        let mut exact = 1.0;
        let mut forces = [0.4; 3];
        for _ in 0..steps {
            let before = exact;
            exact *= factor;
            let middle = 0.5 * (before + exact);
            let report = pair
                .step(
                    &mut positive,
                    &mut negative,
                    h,
                    PlateSustainedLoads::held(&[0.0; 9]),
                    PlateSustainedLoads::held(&[0.0; 9]),
                    &forces,
                    &mut || false,
                )
                .unwrap();
            audit(&report);
            for r in [0, 3, 6] {
                near(positive.v[r], exact, 3e-10);
                near(negative.v[r], exact, 3e-10);
            }
            for (i, force) in report.contact.forces_n.iter().enumerate() {
                near(force[0], 0.4 * middle, 3e-9);
                forces[i] = force[0];
            }
            // Mass damping transfers momentum to the modeled environment.
            // Its nonzero impulse must not be hidden inside support/contact.
            near(report.positive.damping_impulse[0], 2.7 * h * middle, 1e-10);
            near(report.negative.damping_impulse[0], 1.8 * h * middle, 1e-10);
            near(
                report.positive.momentum_change[0] + report.negative.momentum_change[0],
                9.0 * (exact - before),
                1e-10,
            );
            near(
                report.positive.damping_midpoint_j + report.negative.damping_midpoint_j,
                4.5 * h * middle * middle,
                1e-10,
            );
            assert!(
                report.positive.support_impulses.is_empty()
                    && report.negative.support_impulses.is_empty()
            );
            near(positive.a[0], forces[0] - 0.9 * positive.v[0], 3e-9);
            near(negative.a[0], -0.5 * forces[0] - 0.3 * negative.v[0], 3e-9);
        }
        errors.push((positive.v[0] - (-0.2_f64).exp()).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            (3.95..4.05).contains(&(pair[0] / pair[1])),
            "refinement errors {errors:?}"
        );
    }
}

#[test]
fn supported_elastic_plates_match_independent_dense_midpoint_saddle_and_complete_reactions() {
    let positive_fixture = fixture(true, 10.0);
    let negative_fixture = fixture(true, 2.0);
    let originals = [&positive_fixture, &negative_fixture];
    let parameters = [(1.0, 0.07, 0.03), (1.7, 0.11, 0.05)];
    let pair = PlateSustainedPair::prepare(
        surface(originals[0], &NoPlateLoad, 1.0, 0.07, 0.03),
        surface(originals[1], &NoPlateLoad, 1.7, 0.11, 0.05),
        &sites(&[(1.0, 0.0)]),
        controls(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut positive = state(6);
    let mut negative = state(6);
    positive.q[0] = 0.1;
    negative.q[0] = 0.1;
    let before = [positive.clone(), negative.clone()];
    let h = 0.05;
    // H v1 -/+ h B^T lambda = -h K q0, B(v1+ - v1-)=0.
    // This independent dense13x13 saddle solve uses complete original matrices,
    // without the production natural map, contact response or cached mobility.
    let mut saddle = vec![0.0; 13 * 13];
    let mut rhs = vec![0.0; 13];
    for owner in 0..2 {
        let f = originals[owner];
        let (scale, alpha, beta) = parameters[owner];
        let mass = f.model.m.to_dense();
        let stiffness = f.model.k.to_dense();
        for r in 0..6 {
            for c in 0..6 {
                let m = scale * mass[6 * r + c];
                let k = stiffness[6 * r + c];
                saddle[13 * (6 * owner + r) + 6 * owner + c] =
                    m + 0.5 * h * (alpha * m + beta * k) + 0.25 * h * h * k;
                rhs[6 * owner + r] -= h * k * before[owner].q[c];
            }
        }
        let sign = if owner == 0 { 1.0 } else { -1.0 };
        saddle[13 * (6 * owner) + 12] = -h * sign;
        saddle[13 * 12 + 6 * owner] = sign;
    }
    lu(&saddle, 13).unwrap().solve(&mut rhs);
    assert!(rhs[12] > 0.0);
    let report = pair
        .step(
            &mut positive,
            &mut negative,
            h,
            PlateSustainedLoads::held(&[0.0; 9]),
            PlateSustainedLoads::held(&[0.0; 9]),
            &[0.0],
            &mut || false,
        )
        .unwrap();
    audit(&report);
    near(report.contact.forces_n[0][0], rhs[12], 3e-9);
    for (owner, (actual, receipt)) in [(&positive, &report.positive), (&negative, &report.negative)]
        .into_iter()
        .enumerate()
    {
        let f = originals[owner];
        let (scale, alpha, beta) = parameters[owner];
        let sign = if owner == 0 { 1.0 } else { -1.0 };
        for r in 0..6 {
            near(actual.v[r], rhs[6 * owner + r], 3e-10);
            near(
                actual.q[r],
                before[owner].q[r] + 0.5 * h * rhs[6 * owner + r],
                3e-10,
            );
        }
        let mut velocity = vec![0.0; 9];
        let mut displacement = vec![0.0; 9];
        for (i, &r) in f.model.dof_map.iter().enumerate() {
            if let Some(r) = r {
                velocity[i] = actual.v[r];
                displacement[i] = 0.5 * (before[owner].q[r] + actual.q[r]);
            }
        }
        let mass = f.mass.to_dense();
        let stiffness = f.stiffness.to_dense();
        let mut expected_damping = 0.0;
        for &(r, impulse) in &receipt.support_impulses {
            let expected: f64 = (0..9)
                .map(|c| {
                    let m = scale * mass[9 * r + c];
                    let k = stiffness[9 * r + c];
                    m * velocity[c]
                        + 0.5 * h * (alpha * m + beta * k) * velocity[c]
                        + h * k * displacement[c]
                })
                .sum();
            near(impulse, expected, 3e-10);
        }
        assert!(receipt.support_impulses.iter().any(|(_, j)| j.abs() > 1e-4));
        let mut endpoint_rhs = vec![0.0; 6];
        let reduced_mass: Vec<_> = f.model.m.to_dense().iter().map(|v| v * scale).collect();
        let reduced_stiffness = f.model.k.to_dense();
        for r in 0..6 {
            for c in 0..6 {
                endpoint_rhs[r] -= (alpha * reduced_mass[6 * r + c]
                    + beta * reduced_stiffness[6 * r + c])
                    * actual.v[c]
                    + reduced_stiffness[6 * r + c] * actual.q[c];
                expected_damping += 0.25
                    * h
                    * actual.v[r]
                    * (alpha * reduced_mass[6 * r + c] + beta * reduced_stiffness[6 * r + c])
                    * actual.v[c];
            }
        }
        endpoint_rhs[0] += sign * rhs[12];
        lu(&reduced_mass, 6).unwrap().solve(&mut endpoint_rhs);
        for (&a, &b) in actual.a.iter().zip(&endpoint_rhs) {
            near(a, b, 3e-9);
        }
        near(receipt.damping_midpoint_j, expected_damping, 1e-10);
        near(
            receipt.contact_work_j,
            sign * rhs[12] * (actual.q[0] - before[owner].q[0]),
            1e-10,
        );
    }
}

#[test]
fn mode_source_budget_cancellation_and_late_restart_refusals_preserve_both_trajectories() {
    let f = fixture(false, 0.0);
    let full = [-2.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0];
    let load = HeldLoad {
        full: &full,
        map: &f.model.dof_map,
        refuse: Cell::new(false),
    };
    let make = |limits| {
        PlateSustainedPair::prepare(
            surface(&f, &load, 1.0, 0.0, 0.0),
            surface(&f, &NoPlateLoad, 2.0, 0.0, 0.0),
            &sites(&f.mesh.nodes),
            controls(),
            limits,
            &mut || false,
        )
    };
    let mut small = budget();
    small.max_workspace_components = 1;
    assert!(make(small).is_err());
    let mut overflow = budget();
    overflow.positive.plate.max_workspace_components = usize::MAX;
    assert!(make(overflow).is_err());
    let pair = make(budget()).unwrap();
    let initial = (state(9), state(9));
    let h = 0.05;
    for kind in 0..5 {
        let mut candidate = initial.clone();
        match kind {
            0 => candidate.0.q[0] = 0.01,
            1 => candidate.0.q[0] = -0.01,
            2 => candidate.0.v[0] = -0.1,
            3 => candidate.0.v[0] = 0.1,
            _ => candidate.1.t = 0.01,
        }
        let before = candidate.clone();
        assert!(
            pair.step(
                &mut candidate.0,
                &mut candidate.1,
                h,
                PlateSustainedLoads::held(&full),
                PlateSustainedLoads::held(&[0.0; 9]),
                &[0.0; 3],
                &mut || false
            )
            .is_err()
        );
        assert_eq!(candidate, before);
    }
    let mut history_budget = budget();
    history_budget.max_history_entries = 0;
    let bounded = make(history_budget).unwrap();
    let mut candidate = initial.clone();
    assert!(
        bounded
            .step(
                &mut candidate.0,
                &mut candidate.1,
                h,
                PlateSustainedLoads::held(&full),
                PlateSustainedLoads::held(&[0.0; 9]),
                &[0.0; 3],
                &mut || false
            )
            .is_err()
    );
    assert_eq!(candidate, initial);
    let mut wrong = full;
    wrong[3] = 0.0;
    assert!(
        pair.step(
            &mut candidate.0,
            &mut candidate.1,
            h,
            PlateSustainedLoads {
                start: &full,
                midpoint: &wrong,
                end: &full
            },
            PlateSustainedLoads::held(&[0.0; 9]),
            &[0.0; 3],
            &mut || false
        )
        .is_err()
    );
    assert_eq!(candidate, initial);
    load.refuse.set(true);
    assert!(matches!(
        pair.step(
            &mut candidate.0,
            &mut candidate.1,
            h,
            PlateSustainedLoads::held(&full),
            PlateSustainedLoads::held(&[0.0; 9]),
            &[0.0; 3],
            &mut || false
        ),
        Err(PlateSustainedError::Source(
            PlateSustainedSourceError::Load(_)
        ))
    ));
    assert_eq!(candidate, initial);
    load.refuse.set(false);
    let mut polls = 0;
    pair.step(
        &mut candidate.0,
        &mut candidate.1,
        h,
        PlateSustainedLoads::held(&full),
        PlateSustainedLoads::held(&[0.0; 9]),
        &[0.0; 3],
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for limit in [0, polls / 2, polls - 1] {
        let mut trial = initial.clone();
        let mut count = 0;
        assert!(
            pair.step(
                &mut trial.0,
                &mut trial.1,
                h,
                PlateSustainedLoads::held(&full),
                PlateSustainedLoads::held(&[0.0; 9]),
                &[0.0; 3],
                &mut || {
                    let stop = count == limit;
                    count += 1;
                    stop
                }
            )
            .is_err()
        );
        assert_eq!(trial, initial);
    }
    let mut retry = initial.clone();
    pair.step(
        &mut retry.0,
        &mut retry.1,
        h,
        PlateSustainedLoads::held(&full),
        PlateSustainedLoads::held(&[0.0; 9]),
        &[0.0; 3],
        &mut || false,
    )
    .unwrap();
    assert_eq!(retry, candidate);

    // A separating load releases initially closed contact without a tensile
    // reaction. The next open state belongs to the original free time owner.
    let separating = full.map(|v| -v);
    let outward = HeldLoad {
        full: &separating,
        map: &f.model.dof_map,
        refuse: Cell::new(false),
    };
    let outward_dynamics = surface(&f, &outward, 1.0, 0.0, 0.0).dynamics;
    let detached = PlateSustainedPair::prepare(
        surface(&f, &outward, 1.0, 0.0, 0.0),
        surface(&f, &NoPlateLoad, 2.0, 0.0, 0.0),
        &sites(&f.mesh.nodes),
        controls(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut states = initial.clone();
    let report = detached
        .step(
            &mut states.0,
            &mut states.1,
            h,
            PlateSustainedLoads::held(&separating),
            PlateSustainedLoads::held(&[0.0; 9]),
            &[0.0; 3],
            &mut || false,
        )
        .unwrap();
    audit(&report);
    for i in 0..3 {
        near(report.contact.forces_n[i][0], 0.0, 1e-12);
        near(report.gap_bounds_m[i][1], h * h, 1e-12);
    }
    let before = states.clone();
    assert!(
        detached
            .step(
                &mut states.0,
                &mut states.1,
                h,
                PlateSustainedLoads::held(&separating),
                PlateSustainedLoads::held(&[0.0; 9]),
                &[0.0; 3],
                &mut || false
            )
            .is_err()
    );
    assert_eq!(states, before);
    OperatorGeneralizedAlpha::new(9, h, 1.0, controls().step)
        .step(&mut states.0, &outward_dynamics, &separating)
        .unwrap();
    near(states.0.q[0], 4.0 * h * h, 1e-11);

    // All admitted source samples are finite and authentic. Only the final
    // native endpoint acceleration is unrepresentable, after contact/work
    // acceptance; neither candidate may publish on this late restart failure.
    struct EndSpike(f64);
    impl PlateLoad for EndSpike {
        fn parameter_count(&self) -> usize {
            0
        }
        fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
            output.fill(0.0);
            if time == self.0 {
                output[0] = f64::MAX;
            }
            Ok(())
        }
        fn forcing_vjp(&self, _: f64, _: &[f64], output: &mut [f64]) -> Result<(), String> {
            output.fill(0.0);
            Ok(())
        }
    }
    let spike = EndSpike(h);
    let late = PlateSustainedPair::prepare(
        surface(&f, &spike, 1.0, 0.0, 0.0),
        surface(&f, &NoPlateLoad, 2.0, 0.0, 0.0),
        &sites(&f.mesh.nodes),
        controls(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut end = [0.0; 9];
    end[0] = f64::MAX;
    let mut candidate = initial.clone();
    assert!(matches!(
        late.step(
            &mut candidate.0,
            &mut candidate.1,
            h,
            PlateSustainedLoads {
                start: &[0.0; 9],
                midpoint: &[0.0; 9],
                end: &end
            },
            PlateSustainedLoads::held(&[0.0; 9]),
            &[0.0; 3],
            &mut || false
        ),
        Err(PlateSustainedError::Source(
            PlateSustainedSourceError::Initialization(_)
        ))
    ));
    assert_eq!(candidate, initial);
}
