//! Actual field/conductor coupling against independent scalar P1 balances.
//! Manufactured coefficients test physics/accounting, not catalog validity.

#![cfg(feature = "magnetic-hysteresis")]

#[path = "hysteretic_eddy/winding.rs"]
mod winding;

use fs_em::planar::{
    eddy::{
        EddyConductor,
        hysteretic::{HystereticEddyControls, HystereticEddyError, PlanarHystereticEddy},
    },
    nonlinear::{
        NonlinearMagneticControls,
        hysteretic::{
            PlanarHystereticControls, PlanarHystereticMagnetostatics, PlanarHystereticState,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        hysteresis::{
            source::{
                HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis,
                ResolvedVectorHysteresis,
            },
            vector::VectorHysteresisControls,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{MaterialPropertySelection, resolve_material_state_point},
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use std::sync::Arc;

const MU: f64 = 0.01;
const CHI: f64 = 0.2;
const PINNING: f64 = 2.0;
const SATURATION: f64 = 1.0;
const DIRECTIONS: [[f64; 2]; 4] = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];

fn sources(
    mu: f64,
    pinning: bool,
    sigma: f64,
) -> (ResolvedVectorHysteresis, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut roles = HysteresisSources {
        background_permeability: PropertyKey::new("background_mu", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![],
    };
    let mut values = vec![
        (roles.background_permeability.clone(), mu),
        (roles.field_limit.clone(), 30.0),
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            sigma,
        ),
    ];
    if pinning {
        let branch = HysteresisBranchSources {
            susceptibility: PropertyKey::new("branch_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("branch_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("branch_saturation", FLUX_DENSITY_DIMS),
        };
        values.extend([
            (branch.susceptibility.clone(), CHI),
            (branch.pinning_field.clone(), PINNING),
            (branch.saturation_polarization.clone(), SATURATION),
        ]);
        roles.branches.push(branch);
    }
    let mut claims = ClaimSet::new();
    for (key, value) in values {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured conducting magnetic play fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic conducting magnetic law".into(),
                phase: "solid".into(),
                process: "manufactured transient field".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let point = QueryPoint::new();
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        resolve_material_state_point(
            &card,
            &point,
            &[basis.requirement().unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        basis,
    )
    .unwrap();
    let scalar = ResolvedScalarHysteresis::resolve(
        card,
        &point,
        roles,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        ResolvedVectorHysteresis::from_isotropic_source(scalar).unwrap(),
        conductor,
    )
}

fn square(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "hysteretic-conductor-square",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

fn paired_squares(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "balanced-hysteretic-conductor",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
        ],
        vec![
            [0, 1, 6],
            [1, 2, 6],
            [2, 3, 6],
            [3, 0, 6],
            [1, 4, 7],
            [4, 5, 7],
            [5, 2, 7],
            [2, 1, 7],
        ],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

fn problem(
    mesh: TriComplex2,
    mu: f64,
    pinning: bool,
    sigma: f64,
) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let (magnetic, conductor) = sources(mu, pinning, sigma);
    let count = mesh.faces().len();
    let owner = PlanarHystereticMagnetostatics::new(
        mesh,
        vec![magnetic; count],
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 100,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap();
    let conductors = (0..count)
        .map(|triangle| EddyConductor {
            triangle,
            region: 7,
            conductor: conductor.clone(),
        })
        .collect();
    (owner, conductors)
}

fn controls() -> HystereticEddyControls {
    HystereticEddyControls {
        magnetic: PlanarHystereticControls {
            nonlinear: NonlinearMagneticControls {
                absolute_residual_tolerance_a_m: 1e-9,
                relative_residual_tolerance: 1e-10,
                max_newton_iterations: 50,
                max_backtracks: 48,
                linear: CheckedCgConfig {
                    tolerance: 1e-13,
                    max_iterations: 100,
                    max_corrections: 2,
                },
            },
            material: VectorHysteresisControls {
                relative_tolerance: 1e-12,
                absolute_field_tolerance_a_m: 1e-12,
                absolute_flux_tolerance_t: 1e-13,
                work_tolerance_j_m3: 2e-10,
                maximum_branch_iterations: 100,
                maximum_inverse_iterations: 80,
                maximum_backtracks: 48,
            },
            energy_tolerance_j: 2e-8,
        },
        current_balance_tolerance_a: 1e-8,
    }
}

fn boundary(count: usize, gauge: f64) -> Vec<(usize, f64)> {
    (0..count).map(|v| (v, gauge)).collect()
}

fn magnetized_square(
    owner: &PlanarHystereticMagnetostatics,
    gauge: f64,
    h: f64,
) -> PlanarHystereticState<'_> {
    let zero = owner.demagnetized(gauge, &mut || false).unwrap();
    owner
        .candidate(
            &zero,
            &[6.0 * h; 4],
            &boundary(4, gauge),
            controls().magnetic,
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-8 * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

fn play(old_z: f64, h: f64) -> f64 {
    old_z
        .clamp(CHI * (h - PINNING), CHI * (h + PINNING))
        .clamp(-SATURATION, SATURATION)
}

fn storage(h: f64, z: f64) -> f64 {
    0.5 * MU * h * h + z * z / (2.0 * CHI)
}

// This explicit scalar play formula and its monotone balance are independent
// of the production vector inverse, tangent, and field Newton implementation.
fn relaxation(old_z: f64, old_b: f64, damping: f64) -> (f64, f64, f64) {
    let mut lo = -30.0;
    let mut hi = 30.0;
    for _ in 0..90 {
        let h = 0.5 * (lo + hi);
        let residual = h + damping * (MU * h + play(old_z, h) - old_b);
        if residual < 0.0 {
            lo = h;
        } else {
            hi = h;
        }
    }
    let h = 0.5 * (lo + hi);
    let z = play(old_z, h);
    (h, z, MU * h + z)
}

#[test]
fn conducting_magnet_relaxes_with_separate_pinning_and_joule_heat_and_retains_remanence() {
    let depth = 0.3;
    for sigma in [18.0, 36.0] {
        for (gauge, sign) in [(0.0, 1.0), (100.0, -1.0)] {
            let (owner, conductors) = problem(square(depth), MU, true, sigma);
            let initial = magnetized_square(&owner, gauge, sign * 12.0);
            let mut state = PlanarHystereticEddy::new(
                &initial,
                conductors,
                &[],
                &boundary(4, gauge),
                0.0,
                controls(),
                &mut || false,
            )
            .unwrap();
            let (mut old_h, mut old_z, mut old_b) = (sign * 12.0, sign, sign * 1.12);
            let mut pinning_heat = 0.0;
            let mut joule_heat = 0.0;
            // B=2*(A_center-gauge), magnetic force=2*depth*H, and
            // C44=sigma*depth/18 after the floating-region mean projection.
            for dt in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0] {
                let (h, z, b) = relaxation(old_z, old_b, sigma / (72.0 * dt));
                let delta_a = 0.5 * (b - old_b);
                let before = depth * storage(old_h, old_z);
                let after = depth * storage(h, z);
                let q_pinning = depth * PINNING * (z - old_z).abs();
                let q_joule = sigma * depth / 18.0 * delta_a * delta_a / dt;
                let field_work = depth * h * (b - old_b);
                let numerical = field_work - (after - before) - q_pinning;
                let next = state
                    .step_to_candidate(
                        state.time_s() + dt,
                        state.conductors().to_vec(),
                        &[0.0; 4],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                let step = next.step();
                near(next.state().potential_wb_m()[4] - gauge, 0.5 * b);
                near(next.state().energy_j(), after);
                near(step.stored_before_j, before);
                near(step.stored_after_j, after);
                near(step.hysteresis_heat_j, q_pinning);
                near(step.joule_heat_j, q_joule);
                near(step.heat_j, q_pinning + q_joule);
                near(step.numerical_dissipation_j, numerical);
                near(step.supplied_work_j, 0.0);
                near(
                    after - before + step.heat_j + step.numerical_dissipation_j,
                    0.0,
                );
                near(step.energy_residual_j, 0.0);
                near(step.magnetic.field_work_j, field_work);
                near(step.magnetic.irreversible_heat_j, q_pinning);
                assert!(step.numerical_dissipation_j >= 0.0);
                assert_eq!(step.region_currents_a.len(), 1);
                assert_eq!(step.region_currents_a[0].0, 7);
                near(step.region_currents_a[0].1, 0.0);
                for (cell, direction) in DIRECTIONS.iter().enumerate() {
                    let magnetic = &next.state().magnetic_state().cell_states()[cell];
                    for axis in 0..2 {
                        near(magnetic.field_a_m()[axis], h * direction[axis]);
                        near(magnetic.polarizations_t()[0][axis], z * direction[axis]);
                        near(
                            next.state().magnetic_state().magnetic_flux_density_t()[cell][axis],
                            b * direction[axis],
                        );
                    }
                    let heat = &step.cells[cell];
                    near(heat.heat_j, q_joule / 4.0);
                    near(heat.volume_m3, depth / 4.0);
                    near(
                        step.magnetic.cells[cell].irreversible_heat_j,
                        q_pinning / 4.0,
                    );
                    for (local, &vertex) in owner.mesh().faces()[cell].iter().enumerate() {
                        let electric =
                            delta_a / dt * if vertex == 4 { -2.0 / 3.0 } else { 1.0 / 3.0 };
                        near(heat.electric_field_v_m[local], electric);
                        near(heat.current_density_a_m2[local], sigma * electric);
                    }
                }
                pinning_heat += step.hysteresis_heat_j;
                joule_heat += step.joule_heat_j;
                (state, _) = next.into_parts();
                (old_h, old_z, old_b) = (h, z, b);
            }
            assert!(pinning_heat > 0.3 && joule_heat > 0.0);
            assert!(sign * old_b > 0.39 && old_h.abs() < 1e-4);
            assert!(
                (state.magnetic_state().cell_states()[0].polarizations_t()[0][0] - sign * 0.4)
                    .abs()
                    < 1e-6
            );
        }
    }
}

#[test]
fn branchless_limit_recovers_exact_backward_euler_conductor_decay() {
    let depth = 0.3;
    let gauge = 17.0;
    for (mu, sigma) in [(1.0, 18.0), (0.4, 36.0)] {
        let (owner, conductors) = problem(square(depth), mu, false, sigma);
        let initial = magnetized_square(&owner, gauge, -2.0 / mu);
        let mut state = PlanarHystereticEddy::new(
            &initial,
            conductors,
            &[],
            &boundary(4, gauge),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let c = sigma * depth / 18.0;
        let k = 4.0 * depth / mu;
        let mut old_a = -1.0;
        for dt in [0.01, 0.1, 0.3] {
            let a = c * old_a / (c + dt * k);
            let delta = a - old_a;
            let next = state
                .step_to_candidate(
                    state.time_s() + dt,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(next.state().potential_wb_m()[4] - gauge, a);
            near(next.step().stored_before_j, 0.5 * k * old_a * old_a);
            near(next.step().stored_after_j, 0.5 * k * a * a);
            near(next.step().joule_heat_j, c * delta * delta / dt);
            near(next.step().numerical_dissipation_j, 0.5 * k * delta * delta);
            assert_eq!(next.step().hysteresis_heat_j, 0.0);
            near(next.step().heat_j, next.step().joule_heat_j);
            near(next.step().supplied_work_j, 0.0);
            near(next.step().energy_residual_j, 0.0);
            assert!(
                next.state()
                    .magnetic_state()
                    .cell_states()
                    .iter()
                    .all(|cell| cell.polarizations_t().is_empty())
            );
            (state, _) = next.into_parts();
            old_a = a;
        }
    }
}

fn driven_minor_loops(dt: f64) -> (f64, f64, f64) {
    let depth = 0.3;
    let sigma = 18.0;
    let gauge = -11.0;
    let (owner, conductors) = problem(paired_squares(depth), MU, true, sigma);
    let initial = owner.demagnetized(gauge, &mut || false).unwrap();
    let mut state = PlanarHystereticEddy::new(
        &initial,
        conductors,
        &[],
        &boundary(6, gauge),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let (mut old_h, mut old_z, mut old_b) = (0.0, 0.0, 0.0);
    let mut total = (0.0, 0.0, 0.0);
    let mut return_state = None;
    for (index, h) in [12.0, 0.0, 6.0, 1.0, 5.5, 1.0, 6.0, 0.0, -12.0, 0.0]
        .into_iter()
        .enumerate()
    {
        let z = play(old_z, h);
        let b = MU * h + z;
        let delta_b = b - old_b;
        // Opposite center potentials make the connected conductor's mean
        // dA/dt zero. Each center then has mass sigma*depth/6. The impressed
        // current includes exactly the amount needed to traverse this H path.
        let drive_h = h + sigma / (24.0 * dt) * delta_b;
        let mut currents = [6.0 * drive_h; 8];
        currents[4..].fill(-6.0 * drive_h);
        let next = state
            .step_to_candidate(
                state.time_s() + dt,
                state.conductors().to_vec(),
                &currents,
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = next.step();
        let stored_before = 2.0 * depth * storage(old_h, old_z);
        let stored_after = 2.0 * depth * storage(h, z);
        let pinning = 2.0 * depth * PINNING * (z - old_z).abs();
        let joule = sigma * depth / (12.0 * dt) * delta_b * delta_b;
        let field_work = 2.0 * depth * h * delta_b;
        let numerical = field_work - (stored_after - stored_before) - pinning;
        near(next.state().potential_wb_m()[6] - gauge, 0.5 * b);
        near(next.state().potential_wb_m()[7] - gauge, -0.5 * b);
        near(step.stored_before_j, stored_before);
        near(step.stored_after_j, stored_after);
        near(step.hysteresis_heat_j, pinning);
        near(step.joule_heat_j, joule);
        near(step.heat_j, pinning + joule);
        near(step.numerical_dissipation_j, numerical);
        near(step.supplied_work_j, 2.0 * depth * drive_h * delta_b);
        near(step.magnetic.field_work_j, field_work);
        near(step.energy_residual_j, 0.0);
        near(
            stored_after - stored_before + step.heat_j + step.numerical_dissipation_j,
            step.supplied_work_j,
        );
        assert_eq!(step.region_currents_a.len(), 1);
        near(step.region_currents_a[0].1, 0.0);
        for cell in 0..8 {
            let sign = if cell < 4 { 1.0 } else { -1.0 };
            let memory = &next.state().magnetic_state().cell_states()[cell];
            for axis in 0..2 {
                near(
                    memory.field_a_m()[axis],
                    sign * h * DIRECTIONS[cell % 4][axis],
                );
                near(
                    memory.polarizations_t()[0][axis],
                    sign * z * DIRECTIONS[cell % 4][axis],
                );
            }
            near(step.cells[cell].heat_j, joule / 8.0);
            near(step.magnetic.cells[cell].irreversible_heat_j, pinning / 8.0);
            for (local, &vertex) in owner.mesh().faces()[cell].iter().enumerate() {
                let electric = if vertex >= 6 {
                    -sign * 0.5 * delta_b / dt
                } else {
                    0.0
                };
                near(step.cells[cell].electric_field_v_m[local], electric);
                near(
                    step.cells[cell].current_density_a_m2[local],
                    sigma * electric,
                );
            }
        }
        if index == 1 {
            return_state = Some(next.state().magnetic_state().cell_states().to_vec());
        } else if index == 7 {
            for (actual, returned) in next
                .state()
                .magnetic_state()
                .cell_states()
                .iter()
                .zip(return_state.as_ref().unwrap())
            {
                for axis in 0..3 {
                    near(
                        actual.flux_density_t()[axis],
                        returned.flux_density_t()[axis],
                    );
                    near(
                        actual.polarizations_t()[0][axis],
                        returned.polarizations_t()[0][axis],
                    );
                }
            }
        }
        total.0 += step.hysteresis_heat_j;
        total.1 += step.joule_heat_j;
        total.2 += step.numerical_dissipation_j;
        (state, _) = next.into_parts();
        (old_h, old_z, old_b) = (h, z, b);
    }
    total
}

#[test]
fn balanced_current_minor_loops_keep_rate_independent_memory_and_resolved_eddy_loss() {
    let fast = driven_minor_loops(0.2);
    let slow = driven_minor_loops(2.0);
    assert!(fast.0 > 0.0 && fast.1 > 0.0 && fast.2 > 0.0);
    near(fast.0, slow.0);
    near(fast.1, 10.0 * slow.1);
    // The same finite B/H path has the same numerical endpoint remainder;
    // neither pinning nor that remainder is an empirical frequency loss map.
    near(fast.2, slow.2);
}

#[test]
fn failed_cancelled_and_foreign_source_trials_leave_the_accepted_magnetic_history_intact() {
    let (owner, conductors) = problem(square(0.3), MU, true, 18.0);
    let initial = magnetized_square(&owner, 0.0, 12.0);
    let state = PlanarHystereticEddy::new(
        &initial,
        conductors,
        &[],
        &boundary(4, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let snapshot = state.magnetic_state().clone();
    let mut limited = controls();
    limited.magnetic.nonlinear.max_newton_iterations = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[0.0; 4],
                limited,
                &mut || false,
            )
            .is_err()
    );
    let mut polls = 0;
    let cancelled = state.step_to_candidate(
        0.1,
        state.conductors().to_vec(),
        &[0.0; 4],
        controls(),
        &mut || {
            polls += 1;
            polls >= 12
        },
    );
    assert!(matches!(cancelled, Err(HystereticEddyError::Cancelled)));
    let (_, foreign) = sources(MU, true, 19.0);
    let foreign = state
        .conductors()
        .iter()
        .map(|source| EddyConductor {
            triangle: source.triangle,
            region: source.region,
            conductor: foreign.clone(),
        })
        .collect();
    assert!(
        state
            .step_to_candidate(0.1, foreign, &[0.0; 4], controls(), &mut || false,)
            .is_err()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb_m(), snapshot.potential_wb_m());
    assert_eq!(state.magnetic_state().cell_states(), snapshot.cell_states());
    assert_eq!(
        state.magnetic_state().magnetic_flux_density_t(),
        snapshot.magnetic_flux_density_t()
    );
    assert_eq!(state.energy_j(), snapshot.stored_energy_j());
    let retry = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.state().time_s(), 0.1);
    assert!(retry.state().energy_j() < state.energy_j());
    assert!(retry.step().hysteresis_heat_j > 0.0 && retry.step().joule_heat_j > 0.0);
    assert_eq!(state.magnetic_state().cell_states(), snapshot.cell_states());
}
