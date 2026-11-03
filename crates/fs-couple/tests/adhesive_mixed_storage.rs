//! Native assemblies retain each original adhesive and unilateral potential.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::dmt::{DmtEnvelope, SourceDmtSphereContact};
use fs_contact::normal_patch::maugis::continuation::MaugisContinuationBudget;
use fs_contact::normal_patch::maugis::mechanical::MaugisAttachedPotential;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, MaugisWorkBudget, SourceMaugisSphereContact,
};
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_contact::normal_patch::{
    FiniteGapChartSamplingAuthority, FiniteGapGrid, FiniteGapResponseCurve,
    FiniteGapResponseCurveAuthority, FiniteGapResponseCurveRequest, NormalLawSource,
    NormalLawSourceError, build_finite_gap_response_curve,
};
use fs_couple::unilateral_contact::binding::{BoundNormalContact, BoundNormalStorage};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage};

fn adhesive() -> (MaugisAttachedPotential, f64, f64) {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let radius = 0.01;
    let a0 =
        (PI * 1e-7 * radius * radius / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    let f0 = PI * 1e-7 * radius;
    let source = SourceMaugisSphereContact::new(
        normal,
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
    .unwrap();
    let potential = MaugisAttachedPotential::new(
        source,
        MaugisContinuationBudget {
            work: MaugisWorkBudget {
                absolute_tolerance_j: 1e-12 * f0 * d0,
                relative_tolerance: 1e-11,
                maximum_evaluations: 2049,
            },
            energy_absolute_tolerance_j: 1e-10 * f0 * d0,
            energy_relative_tolerance: 1e-10,
        },
    )
    .unwrap()
    .with_noncontact_continuation()
    .unwrap();
    (potential, d0, f0)
}

fn point(d0: f64, chi: Option<f64>) -> SourcePointContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let mut envelope = common::envelope();
    envelope.maximum_indentation_m = 20.0 * d0;
    SourcePointContact::new(
        common::bind(&a, &b, chi, 0.0),
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.01,
        },
        envelope,
    )
    .unwrap()
}

fn dmt(owner: &MaugisAttachedPotential, d0: f64) -> SourceDmtSphereContact {
    SourceDmtSphereContact::new(
        owner.source().normal_state().clone(),
        0.01,
        10.0 * d0,
        DmtEnvelope {
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
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}

fn table(owner: &MaugisAttachedPotential, d0: f64, f0: f64) -> FiniteGapResponseCurve {
    let radius = owner.source().effective_radius_m();
    let width = (radius * d0).sqrt();
    let grid = FiniteGapGrid {
        cells_x: 11,
        cells_y: 11,
        cell_width_m: width,
        cell_depth_m: width,
    };
    let mut gap = Vec::new();
    for index in 0..121 {
        let x = ((index % 11) as f64 - 5.0) * width;
        let y = ((index / 11) as f64 - 5.0) * width;
        gap.push((x * x + y * y) / (2.0 * radius));
    }
    let request = FiniteGapResponseCurveRequest {
        gap_source_identity: fs_blake3::hash_domain("test/mixed-adhesion", b"paraboloid"),
        gap_source_authority: FiniteGapChartSamplingAuthority::Estimate,
        grid,
        undeformed_gap_m: gap,
        reduced_modulus_pa: owner.source().normal_state().elastic().reduced_modulus_pa(),
        approach_nodes_m: (0..=16).map(|i| 0.25 * i as f64 * d0).collect(),
        maximum_active_set_iterations: 256,
        complementarity_tolerance_m: 1e-9 * d0,
        boundary_clearance_cells: 2,
        absolute_energy_tolerance_j: 1e-10 * f0 * d0,
        relative_energy_tolerance: 0.05,
    };
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        build_finite_gap_response_curve(&request, &cx).unwrap()
    })
}

fn storage(masses: &[f64], contacts: Vec<BoundNormalContact>) -> BoundNormalStorage {
    let n = 2 * masses.len();
    let mut q = vec![0.0; n * n];
    for (i, mass) in masses.iter().enumerate() {
        q[(2 * i + 1) * n + 2 * i + 1] = 1.0 / mass;
    }
    BoundNormalStorage::new(
        Box::new(QuadraticStorage::new(q, n).unwrap()),
        masses.len(),
        contacts,
    )
    .unwrap()
}

fn system(storage: BoundNormalStorage, bodies: usize) -> PortHamiltonian {
    let n = 2 * bodies;
    let mut j = vec![0.0; n * n];
    for i in 0..bodies {
        j[2 * i * n + 2 * i + 1] = 1.0;
        j[(2 * i + 1) * n + 2 * i] = -1.0;
    }
    PortHamiltonian::new(n, 0, j, vec![0.0; n * n], vec![], Box::new(storage)).unwrap()
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= scale * tolerance,
        "{actual:e} != {expected:e}; scaled defect {:e}",
        (actual - expected).abs() / scale
    );
}

#[test]
fn original_hertz_and_hunt_crossley_storage_work_keep_the_inactive_gap_fraction() {
    let (_, d0, f0) = adhesive();
    for chi in [None, Some(0.4)] {
        let point = point(d0, chi);
        let source = NormalLawSource::from(point.clone());
        let k = point.stiffness_n_per_m_3_2();
        let energy = |d: f64| 0.4 * k * d.max(0.0).powf(2.5);
        for (a, b) in [(-2.0 * d0, d0), (d0, 3.0 * d0), (3.0 * d0, -d0)] {
            let force = source.mean_force_between(a, b).unwrap().unwrap();
            close(force * (b - a), energy(b) - energy(a), f0 * d0, 3e-14);
            close(
                force,
                source.mean_force_between(b, a).unwrap().unwrap(),
                f0,
                3e-14,
            );
        }
        assert_eq!(
            source.mean_force_between(-d0, -2.0 * d0).unwrap(),
            Some(0.0)
        );
        assert_eq!(source.mean_force_between(0.0, 0.0).unwrap(), Some(0.0));
        let tiny = point.conservative_path_average(d0, d0.next_up()).unwrap();
        close(tiny.mean_force_n, k * d0.powf(1.5), f0, 2e-14);
        assert!(tiny.work_j > 0.0);
        if chi.is_none() {
            assert_eq!(
                point.conservative_path_average(d0, 2.0 * d0).unwrap(),
                point.elastic_path_average(d0, 2.0 * d0, 1e-3).unwrap()
            );
        } else {
            assert!(point.elastic_path_average(d0, 2.0 * d0, 1e-3).is_err());
        }
        assert!(point.conservative_path_average(-d0, d0).is_err());
        assert!(source.mean_force_between(d0, 21.0 * d0).is_err());
        assert!(source.mean_force_between(f64::NAN, d0).is_err());
    }
}

#[test]
fn retained_bem_curve_supplies_exact_local_work_without_promoting_its_authority() {
    let (owner, d0, f0) = adhesive();
    let table = table(&owner, d0, f0);
    let identity = table.identity;
    let source = NormalLawSource::from(table.clone());
    assert_eq!(source.identity(), identity);
    assert_eq!(
        table.response_authority,
        FiniteGapResponseCurveAuthority::Estimate
    );
    for (a, b) in [(0.1 * d0, 3.6 * d0), (2.9 * d0, 0.6 * d0), (-d0, 2.1 * d0)] {
        let ua = table.evaluate(a.max(0.0)).unwrap().reversible_energy_j;
        let ub = table.evaluate(b).unwrap().reversible_energy_j;
        let mean = source.mean_force_between(a, b).unwrap().unwrap();
        close(mean * (b - a), ub - ua, f0 * d0, 3e-14);
        assert_eq!(mean, source.mean_force_between(b, a).unwrap().unwrap());
        let mid = a.midpoint(b);
        let split = source.mean_force_between(a, mid).unwrap().unwrap() * (mid - a)
            + source.mean_force_between(mid, b).unwrap().unwrap() * (b - mid);
        close(split, mean * (b - a), f0 * d0, 3e-14);
    }
    for at in [0.6 * d0, d0, 2.5 * d0] {
        close(
            table.mean_force_between(at, at.next_up()).unwrap(),
            table.evaluate(at).unwrap().normal_force_n,
            f0,
            2e-14,
        );
    }
    assert!(table.mean_force_between(0.0, (4.0 * d0).next_up()).is_err());
    assert!(table.mean_force_between(f64::NAN, d0).is_err());
    let mut corrupt = table;
    corrupt.nodes[3].curve_reversible_energy_j += f0 * d0;
    assert!(corrupt.mean_force_between(0.0, d0).is_err());
}

#[test]
fn all_four_original_normal_owners_share_reciprocal_three_body_native_work() {
    let (owner, d0, f0) = adhesive();
    let point = NormalLawSource::from(point(d0, None));
    let table = NormalLawSource::from(table(&owner, d0, f0));
    let dmt = dmt(&owner, d0);
    let contacts = vec![
        BoundNormalContact::new(owner.clone().into(), vec![1.0, -1.0, 0.0], 0.0).unwrap(),
        BoundNormalContact::new(point.clone(), vec![0.0, 1.0, -1.0], 0.2 * d0).unwrap(),
        BoundNormalContact::new(table.clone(), vec![0.0, 1.0, -1.0], 0.4 * d0).unwrap(),
        BoundNormalContact::new(dmt.clone().into(), vec![1.0, 0.0, -1.0], 0.0).unwrap(),
    ];
    let bound = storage(&[1.0, 3.0, 2.0], contacts);
    assert!(bound.discrete_gradient_rule().is_some());
    assert!(bound.coordinate_scale(4).unwrap() > 0.0);
    let p0 = (f0 * d0).sqrt();
    let a = [d0, 0.2 * p0, 0.0, -0.3 * p0, d0, 0.4 * p0];
    let b = [1.8 * d0, 0.3 * p0, 0.3 * d0, -0.2 * p0, -0.5 * d0, 0.1 * p0];
    let effort = fs_phs::discrete_gradient(&bound, &a, &b);
    let fm = owner.path(d0, 1.5 * d0).unwrap().mean_force_n;
    // Independent endpoint energies on finite paths expose missing inactive
    // fractions: both unilateral sites start OPEN, with different clearances.
    let fh = point.energy_and_force(0.6 * d0).unwrap().0 / (1.8 * d0);
    let ft = table.energy_and_force(0.4 * d0).unwrap().0 / (1.8 * d0);
    let fd = dmt.reversible_work_between(0.0, 2.3 * d0).unwrap() / (2.3 * d0);
    close(effort[0], fm + fd, f0, 3e-10);
    close(effort[2], -fm + fh + ft, f0, 3e-10);
    close(effort[4], -fh - ft - fd, f0, 3e-10);
    close(effort[0] + effort[2] + effort[4], 0.0, f0, 3e-14);
    for (i, mass) in [1.0, 3.0, 2.0].iter().enumerate() {
        close(
            effort[2 * i + 1],
            (0.5 * a[2 * i + 1] + 0.5 * b[2 * i + 1]) / mass,
            p0,
            2e-14,
        );
    }
    let work: f64 = effort
        .iter()
        .zip(a.iter().zip(b))
        .map(|(e, (a, b))| e * (b - a))
        .sum();
    close(
        work,
        bound.hamiltonian(&b) - bound.hamiltonian(&a),
        f0 * d0,
        5e-10,
    );
}

#[test]
fn native_mixed_patches_make_and_lose_unilateral_contact_with_conserved_momentum_and_energy() {
    let (owner, d0, f0) = adhesive();
    for partner in [
        NormalLawSource::from(point(d0, None)),
        NormalLawSource::from(table(&owner, d0, f0)),
    ] {
        let bound = storage(
            &[1.0, 3.0],
            vec![
                BoundNormalContact::new(owner.clone().into(), vec![1.0, -1.0], 0.0).unwrap(),
                BoundNormalContact::new(partner, vec![1.0, -1.0], 0.5 * d0).unwrap(),
            ],
        );
        let system = system(bound, 2);
        let mut workspace = StepWorkspace::new(&system).unwrap();
        let limits = owner.noncontact_continuation().unwrap();
        let start = limits
            .minimum_indentation_m
            .midpoint(limits.first_contact_indentation_m);
        let mut x = [0.75 * start, 0.0, -0.25 * start, 0.0];
        let initial = system.hamiltonian(&x);
        let dt = 0.03 * (0.75 * d0 / f0).sqrt();
        let p0 = (0.75 * f0 * d0).sqrt();
        let mut closed = false;
        let mut reopened = false;
        for _ in 0..400 {
            let mut next = [f64::NAN; 4];
            let report = workspace
                .step_into(&system, &x, &[], dt, &mut next, &mut [])
                .unwrap();
            x = next;
            if x[0] - x[2] > 0.5 * d0 {
                closed = true;
            } else {
                reopened |= closed;
            }
            close(x[1] + x[3], 0.0, p0, 3e-9);
            close(x[0] + 3.0 * x[2], 0.0, d0, 3e-8);
            close(system.hamiltonian(&x), initial, f0 * d0, 5e-8);
            assert_eq!(report.dissipated, 0.0);
            assert_eq!(report.supplied, 0.0);
        }
        assert!(closed && reopened);
    }
}

#[test]
fn original_rate_loss_and_source_refusals_survive_mixed_native_binding() {
    let (owner, d0, f0) = adhesive();
    let point = point(d0, Some(0.4));
    let bound = storage(
        &[1.0, 3.0],
        vec![
            BoundNormalContact::new(owner.clone().into(), vec![1.0, -1.0], 0.0).unwrap(),
            BoundNormalContact::new(point.clone().into(), vec![1.0, -1.0], 0.5 * d0).unwrap(),
        ],
    );
    let x = [d0, 0.1, 0.0, 0.0];
    let effort = fs_phs::discrete_gradient(&bound, &x, &x);
    let mut flow = [f64::NAN; 4];
    bound
        .dissipative_flow_into(&x, &effort, 1e-3, &mut flow)
        .unwrap();
    let loss_force = 0.4 * 0.1 * point.stiffness_n_per_m_3_2() * (0.5 * d0).powf(1.5);
    close(flow[1], loss_force, f0, 2e-14);
    assert_eq!(flow[0], 0.0);
    assert_eq!(flow[2], 0.0);
    assert_eq!(flow[1], -flow[3]);
    let mut too_fast = effort.clone();
    too_fast[1] = point.envelope().maximum_absolute_rate_m_per_s.next_up();
    assert!(
        bound
            .dissipative_flow_into(&x, &too_fast, 1e-3, &mut flow)
            .is_err()
    );
    assert!(
        bound
            .dissipative_flow_into(&x, &effort, (1e-3_f64).next_up(), &mut flow)
            .is_err()
    );

    let table = NormalLawSource::from(table(&owner, d0, f0));
    assert!(matches!(
        table.mean_force_between(d0, 5.0 * d0),
        Err(NormalLawSourceError::FiniteGap(_))
    ));
    let bound = storage(
        &[1.0, 3.0],
        vec![
            BoundNormalContact::new(owner.into(), vec![1.0, -1.0], 0.0).unwrap(),
            BoundNormalContact::new(table, vec![1.0, -1.0], 0.5 * d0).unwrap(),
        ],
    );
    let system = system(bound, 2);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let mut output = [17.0; 4];
    let dt = 0.01 * (0.75 * d0 / f0).sqrt();
    assert!(
        workspace
            .step_into(
                &system,
                &[5.0 * d0, 0.0, 0.0, 0.0],
                &[],
                dt,
                &mut output,
                &mut []
            )
            .is_err()
    );
    assert_eq!(output, [17.0; 4]);
    let accepted = [d0, 0.0, 0.0, 0.0];
    let reference = fs_phs::step(&system, &accepted, &[], dt).unwrap();
    workspace
        .step_into(&system, &accepted, &[], dt, &mut output, &mut [])
        .unwrap();
    for i in 0..4 {
        close(
            output[i],
            reference.x[i],
            if i % 2 == 0 { d0 } else { (f0 * d0).sqrt() },
            3e-8,
        );
    }
}
